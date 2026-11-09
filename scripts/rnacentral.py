"""
RNAcentral lookups for Stockholm sequences named with RNAcentral identifiers.

Sequence names such as ``URS0000572555_1268238/1-81`` refer to an RNAcentral
unique sequence (URS), optionally qualified by an NCBI taxid, with 1-based
coordinates inside that RNAcentral sequence. These are validated against the
RNAcentral REST API instead of NCBI.

Two endpoints are used:
- /api/v1/rna/<URS>/ (fast, reliable): sequence and is_active. This is the
  definitive check of the sequence data.
- /api/v1/rna/<URS>_<taxid>/ (species-specific): 200 if the URS is annotated
  in that taxid, 404 if not. For URS found in many species this endpoint can
  be very slow or return HTTP 500, in which case the taxid is reported as
  unverified rather than failing an otherwise correct sequence.

API rules (https://rnacentral.org/help/rest-api): /api/v1/rna/ is limited to
30 requests/minute per IP (one request every 2 s recommended); on HTTP 429
honour Retry-After, otherwise back off with jitter; send a descriptive
User-Agent.
"""

import http.client
import json
import random
import re
import time
import urllib.error
import urllib.request

from scripts import config

API_BASE = 'https://rnacentral.org/api/v1/rna'

# URS + 10 hex digits, optional _taxid
_URS_RE = re.compile(r'^(URS[0-9A-F]{10})(?:_(\d+))?$')

# Per-run caches
_seq_cache = {}     # urs -> entry dict or None
_taxid_cache = {}   # (urs, taxid) -> 'ok' | 'wrong' | 'unknown'
_last_request = [0.0]

# Ids whose taxid could not be verified (species endpoint errored/timed out).
# Read by the caller to report a summary.
unverified_taxids = []


class RNAcentralError(RuntimeError):
    """Persistent failure talking to RNAcentral (not a 404)."""


def is_rnacentral_id(accession):
    """Return True if the accession looks like an RNAcentral id (URS or URS_taxid)."""
    return bool(_URS_RE.match(accession or ''))


def split_rnacentral_id(accession):
    """Split 'URS..._taxid' into (urs, taxid or None)."""
    m = _URS_RE.match(accession)
    if not m:
        raise ValueError(f"Not an RNAcentral id: {accession}")
    return m.group(1), m.group(2)


def _throttle():
    wait = config.RNACENTRAL_REQUEST_DELAY - (time.monotonic() - _last_request[0])
    if wait > 0:
        time.sleep(wait)
    _last_request[0] = time.monotonic()


def _backoff(attempt, retry_after=None):
    try:
        return float(retry_after)
    except (TypeError, ValueError):
        return (2 ** attempt) * config.RNACENTRAL_REQUEST_DELAY + random.uniform(0, 1)


def _get_json(url, max_retries, timeout):
    """
    GET a JSON document. Returns the parsed object, or None on HTTP 404.
    Raises RNAcentralError when the request keeps failing.
    """
    request = urllib.request.Request(url, headers={
        'Accept': 'application/json',
        'User-Agent': config.RNACENTRAL_USER_AGENT,
    })
    last_error = None
    for attempt in range(max_retries + 1):
        _throttle()
        try:
            with urllib.request.urlopen(request, timeout=timeout) as resp:
                return json.loads(resp.read().decode('utf-8'))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last_error = e
            retryable = e.code == 429 or e.code >= 500
            if retryable and (attempt < max_retries or e.code == 429):
                time.sleep(_backoff(attempt, e.headers.get('Retry-After') if e.headers else None))
                if attempt < max_retries:
                    continue
            break
        except (http.client.HTTPException, OSError, ValueError) as e:
            # OSError: URLError, timeouts, connection resets.
            # HTTPException: IncompleteRead, RemoteDisconnected (not wrapped by urllib).
            # ValueError: malformed JSON or undecodable body.
            last_error = e
            if attempt < max_retries:
                time.sleep(_backoff(attempt))
                continue
            break
    raise RNAcentralError(f"RNAcentral request failed for {url}: {last_error}")


def fetch_entry(urs):
    """
    Fetch the (species-independent) RNAcentral entry for a bare URS.

    Returns:
        dict with keys sequence (uppercase, U->T) and is_active,
        or None if RNAcentral has no such URS (404).
    Raises:
        RNAcentralError on persistent network/API failure.
    """
    if urs in _seq_cache:
        return _seq_cache[urs]
    data = _get_json(f"{API_BASE}/{urs}/?format=json",
                     max_retries=config.RNACENTRAL_MAX_RETRIES,
                     timeout=config.RNACENTRAL_TIMEOUT)
    entry = None
    if data and data.get('sequence'):
        entry = {
            'sequence': data['sequence'].upper().replace('U', 'T'),
            'is_active': data.get('is_active', True),
        }
    _seq_cache[urs] = entry
    return entry


def check_taxid(urs, taxid):
    """
    Check that a URS is annotated in a taxid.

    Returns:
        'ok' (annotated), 'wrong' (404: not annotated in this taxid) or
        'unknown' (endpoint errored or timed out).
    """
    key = (urs, taxid)
    if key in _taxid_cache:
        return _taxid_cache[key]
    try:
        data = _get_json(f"{API_BASE}/{urs}_{taxid}/?format=json",
                         max_retries=config.RNACENTRAL_TAXID_MAX_RETRIES,
                         timeout=config.RNACENTRAL_TAXID_TIMEOUT)
        result = 'ok' if data is not None else 'wrong'
    except RNAcentralError:
        result = 'unknown'
    _taxid_cache[key] = result
    return result


def validate_sequence(accession, start, end, stockholm_seq, verbose=False):
    """
    Check that stockholm_seq equals the RNAcentral sequence[start..end]
    (1-based, inclusive) and, if a taxid is given, that the URS is annotated
    in that taxid.

    Args:
        accession: 'URS...' or 'URS..._taxid'
        start, end: coordinates from the Stockholm name
        stockholm_seq: ungapped sequence (any case, U or T)

    Returns:
        tuple (ok: bool, failure: None | 'not_found' | 'mismatch')
    """
    urs, taxid = split_rnacentral_id(accession)

    try:
        entry = fetch_entry(urs)
    except RNAcentralError as e:
        if verbose:
            print(f"    RNAcentral lookup failed for {urs}: {e}")
        return False, 'not_found'

    if entry is None:
        if verbose:
            print(f"    {urs}: NOT FOUND in RNAcentral")
        return False, 'not_found'

    if not entry['is_active']:
        if verbose:
            print(f"    {urs}: RNAcentral entry is inactive (obsolete)")
        return False, 'not_found'

    ref = entry['sequence']
    query = stockholm_seq.upper().replace('U', 'T')

    if start > end:
        if verbose:
            print(f"    {accession}: reverse coordinates ({start}-{end}) are not valid for an RNAcentral sequence")
        return False, 'mismatch'
    if start < 1 or end > len(ref):
        if verbose:
            print(f"    {accession}: coordinates {start}-{end} outside RNAcentral sequence length {len(ref)}")
        return False, 'mismatch'

    sub = ref[start - 1:end]
    if sub != query:
        if verbose:
            if len(sub) != len(query):
                print(f"    Length: Stockholm={len(query)}, RNAcentral={len(sub)}")
            else:
                for i, (a, b) in enumerate(zip(query, sub)):
                    if a != b:
                        print(f"    First difference at position {i + 1}: Stockholm={a}, RNAcentral={b}")
                        break
        return False, 'mismatch'

    # Sequence matches; check taxid annotation if one is given
    if taxid:
        status = check_taxid(urs, taxid)
        if status == 'wrong':
            if verbose:
                print(f"    {accession}: {urs} is not annotated in taxid {taxid}")
            return False, 'not_found'
        if status == 'unknown':
            if accession not in unverified_taxids:
                unverified_taxids.append(accession)
            if verbose:
                print(f"    {accession}: sequence OK, but taxid {taxid} could not be verified "
                      f"(RNAcentral species endpoint unavailable) - keeping")

    return True, None


def find_coordinates(accession, stockholm_seq):
    """
    Locate an ungapped sequence inside its RNAcentral entry.

    Returns:
        'ACCESSION/START-END' on success, or None.
    Raises:
        RNAcentralError on persistent network/API failure.
    """
    urs, _ = split_rnacentral_id(accession)
    entry = fetch_entry(urs)
    if entry is None or not entry['is_active']:
        return None
    query = stockholm_seq.upper().replace('U', 'T')
    pos = entry['sequence'].find(query)
    if pos == -1:
        return None
    return f"{accession}/{pos + 1}-{pos + len(query)}"
