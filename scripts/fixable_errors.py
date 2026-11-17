"""
Fixable error handlers for Stockholm files.

These errors can be automatically corrected.
"""

import errno
import io
import os
import random
import sys
import time
import re
import subprocess
import warnings

# Suppress BioPython warnings
warnings.filterwarnings("ignore", module="Bio")
warnings.filterwarnings("ignore", message=".*FASTA.*")

# Import configuration
from scripts import config
from scripts import rnacentral


class RunAbortError(RuntimeError):
    """A failure that must stop the run: carrying on would remove sequences
    for reasons that have nothing to do with the sequences themselves."""


class LocalBlastError(RunAbortError):
    """Local blastn could not run (missing binary or database)."""


class NCBIUnavailableError(RunAbortError):
    """NCBI could not be used: it kept refusing requests (HTTP 429 rate limit
    or 5xx), could not be reached at all (no connection, DNS, timeouts), or
    rejected the API key."""


class CacheError(RunAbortError):
    """The sequence cache folder cannot be written (read-only, full or over quota)."""


# Errors from writing a file that mean the cache folder is unusable. Permission
# errors count only when they name a file, so a refused network connection is
# never mistaken for a cache problem.
_DISK_FULL_ERRNOS = {errno.ENOSPC, errno.EROFS, getattr(errno, 'EDQUOT', errno.ENOSPC)}


def _is_cache_write_error(e):
    if not isinstance(e, OSError):
        return False
    if e.errno in _DISK_FULL_ERRNOS:
        return True
    return e.errno in (errno.EACCES, errno.EPERM) and getattr(e, 'filename', None) is not None


def _cache_error(e, folder):
    return CacheError(
        f"Cannot write to the sequence cache folder {folder} ({e.strerror or e}). "
        f"Set RFAM_SEED_QC_CACHE to a writable folder with free space")


_cache_folder = []   # resolved once per run


def _cache_dir():
    """Folder where downloaded and local-database records are cached.

    $RFAM_SEED_QC_CACHE if set, else $XDG_CACHE_HOME/rfam-seed-qc/fasta
    (default ~/.cache/rfam-seed-qc/fasta). Kept outside the install folder,
    which is read-only in containers and shared installs. If it cannot be
    written the run stops: otherwise every lookup would fail and valid
    sequences would be removed as "not found".
    """
    if _cache_folder:
        return _cache_folder[0]
    import tempfile
    folder = os.environ.get('RFAM_SEED_QC_CACHE', '').strip()
    if not folder:
        base = (os.environ.get('XDG_CACHE_HOME', '').strip()
                or os.path.join(os.path.expanduser('~'), '.cache'))
        folder = os.path.join(base, 'rfam-seed-qc', 'fasta')
    folder = os.path.abspath(os.path.expanduser(folder))
    try:
        os.makedirs(folder, exist_ok=True)
        fd, probe = tempfile.mkstemp(dir=folder, suffix='.part')
        os.close(fd)
        os.remove(probe)
    except OSError as e:
        raise _cache_error(e, folder) from e
    _cache_folder.append(folder)
    return folder


def _cache_tempfile(folder):
    """Temporary file in the cache folder; a write failure stops the run."""
    import tempfile
    try:
        return tempfile.mkstemp(dir=folder, suffix='.part')
    except OSError as e:
        raise _cache_error(e, folder) from e


def _cache_store(tmp_path, fasta_file):
    """Move a finished download into place (atomic); a failure stops the run."""
    try:
        os.replace(tmp_path, fasta_file)
    except OSError as e:
        raise _cache_error(e, os.path.dirname(fasta_file)) from e


class RNAcentralUnavailableError(RunAbortError):
    """RNAcentral could not be reached or kept failing (not a 404)."""


def _ncbi_refused(e):
    """True for HTTP errors that mean 'try later' (rate limit, server trouble)."""
    import urllib.error
    return isinstance(e, urllib.error.HTTPError) and (e.code == 429 or e.code >= 500)


def _is_network_error(e):
    """True if the server could not be reached at all (no connection, DNS
    failure, timeout, dropped connection). An HTTP error status is not one:
    the server answered."""
    import http.client
    import urllib.error
    if isinstance(e, urllib.error.HTTPError):
        return False
    return isinstance(e, (OSError, http.client.HTTPException))


def parse_sequence_identifier(seq_name):
    """
    Parse sequence identifier to extract accession and coordinates.
    
    Format: ACCESSION/START-END (e.g., AF228364.1/1-74)
    
    Args:
        seq_name: Sequence identifier string
        
    Returns:
        tuple: (accession, coordinates) or (seq_name, None) if no coordinates found
    """
    if '/' in seq_name:
        parts = seq_name.split('/', 1)
        return parts[0], parts[1]
    return seq_name, None


def find_missing_coordinates(sequence_entries):
    """
    Find sequences that are missing start/stop coordinates.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples

    Returns:
        list: Sequence names missing coordinates
    """
    missing = []
    for seq_name, _ in sequence_entries:
        _, coords = parse_sequence_identifier(seq_name)
        if coords is None:
            missing.append(seq_name)
    return missing


class AccessionNotFoundError(RuntimeError):
    """NCBI efetch reported that the identifier does not exist."""


def _read_head(path, size=256):
    """Return the first bytes of a file as text, ignoring NULs (NCBI error
    replies can be UTF-16 encoded)."""
    with open(path, 'rb') as f:
        head = f.read(size)
        if head.startswith((b'\xff\xfe', b'\xfe\xff')):
            return head.decode('utf-16', errors='replace')
        return head.replace(b'\x00', b'').decode('utf-8', errors='replace')


def _is_fasta_file(path):
    """True if the file exists and looks like FASTA (first non-blank char is '>')."""
    try:
        return os.path.getsize(path) > 0 and _read_head(path).lstrip().startswith('>')
    except OSError:
        return False


# NCBI nucleotide accession formats, optionally versioned: INSDC (X03205,
# AF265547, AC277891), WGS/TSA (JAHKSV010000123) and RefSeq (NR_038927,
# NC_000913, NZ_CP045771, NZ_JAHKSV010000123). With --local, only these fall
# back to NCBI when missing from the local database; other labels (species
# names, metagenome contigs such as SRS013261_C1604114, plain numbers) go
# straight to BLAST, so offline runs do not wait on network timeouts.
_ACCESSION_RE = re.compile(r'^(?:[A-Z]{2}_[A-Z]{0,6}|[A-Z]{1,6})\d{5,12}(\.\d+)?$')

# An accession at the start of a sequence name, followed by anything that is
# not part of it: NZ_CP176678.1_2667588, CP061243.1_2605127, CP061243.1/1-84.
_ACCESSION_HINT_RE = re.compile(
    r'^((?:[A-Z]{2}_[A-Z]{0,6}|[A-Z]{1,6})\d{5,12})(?:\.(\d+))?(?![A-Za-z0-9])')


def _accession_core(accession):
    """Accession without version or RefSeq NZ_ prefix, for comparing records:
    NZ_CP176678.1 and CP176678.2 both give CP176678. RefSeq NZ_ records are
    copies of the GenBank record with the same core, and nt holds the GenBank
    one."""
    accession = re.sub(r'\.\d+$', '', accession)
    return accession[3:] if accession.startswith('NZ_') else accession


def _accession_hint(seq_name):
    """
    The accession a sequence name starts with, if any, e.g. from
    NZ_CP176678.1_2667588. Used to prefer the submitter's own record.

    Returns:
        (core, candidates, start) or None.
        core: as in _accession_core (CP176678).
        candidates: the accessions to look up, GenBank form first, then the
            name's own form, then unversioned (CP176678.1, NZ_CP176678.1,
            CP176678), in case the database holds another version. A
            different version is safe to use: the sequence must still match
            it exactly.
        start: the number after the accession when the name is just
            ACCESSION_NUMBER (2667588), else None. It is only used when an
            exact match starts exactly there (see _locate_in_record).
        None for names that do not start with an accession (species names,
        RNAcentral ids).
    """
    accession = parse_sequence_identifier(seq_name)[0]
    if rnacentral.is_rnacentral_id(accession):
        return None
    m = _ACCESSION_HINT_RE.match(accession)
    if not m:
        return None
    version = f'.{m.group(2)}' if m.group(2) else ''
    named = m.group(1) + version
    core = _accession_core(m.group(1))
    candidates = []
    for candidate in (core + version, named, core):
        if candidate not in candidates:
            candidates.append(candidate)
    start = None
    if seq_name == accession:  # no /coordinates in the name
        number = re.fullmatch(r'_(\d+)', accession[m.end():])
        if number and int(number.group(1)) > 0:
            start = int(number.group(1))
    return core, candidates, start


def _fetch_local(identifier, fasta_file, fasta_folder):
    """With --local, copy the record from the local BLAST database into the
    FASTA cache. Returns True on success, False if the database lacks it.

    Raises LocalBlastError for any other failure (unreadable volume, broken
    index, I/O error): treating those as "not found" could remove valid
    sequences, so the run stops instead.
    """
    if ',' in identifier:
        return False  # blastdbcmd -entry would read a comma as a list of entries
    fd, tmp_path = _cache_tempfile(fasta_folder)
    try:
        with os.fdopen(fd, 'wb') as out:
            proc = subprocess.run(
                ['blastdbcmd', '-db', config.BLAST_LOCAL_DB, '-entry', identifier,
                 '-target_only', '-outfmt', '%f'],
                stdout=out, stderr=subprocess.PIPE, text=True)
        if proc.returncode == 0 and _is_fasta_file(tmp_path):
            _cache_store(tmp_path, fasta_file)  # atomic: cache only valid FASTA
            return True
        if 'Entry not found' in proc.stderr:
            return False
        raise LocalBlastError(f"blastdbcmd failed for {identifier}: "
                              f"{proc.stderr.strip() or 'exit code %d' % proc.returncode}")
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def get_fasta_file(identifier, max_retries=3, allow_ncbi=True):
    """
    Download fasta file from NCBI if not already downloaded.
    Retries on transient network failures (longer waits when NCBI rate-limits).

    With --local (config.BLAST_LOCAL_DB set), the record is read from the
    local BLAST database first; NCBI is used only for real-looking accessions
    the local database lacks, and only if allow_ncbi is True.

    NCBI efetch answers an unknown identifier with HTTP 200 and an
    'Error: ...' body, so the download is checked to be real FASTA before it
    is cached. Error replies are never cached, and invalid files left in the
    cache by older versions are discarded and re-downloaded.

    Downloads use Python's urllib (no external curl/wget needed). A network
    failure is never reported as a missing accession.

    Raises:
        AccessionNotFoundError: NCBI says the identifier does not exist.
        NCBIUnavailableError: NCBI kept refusing requests, could not be
            reached, or rejected the API key (stops the run).
        CacheError: the sequence cache cannot be written (stops the run).
        RuntimeError: any other failure that persists after the retries,
            such as a reply that is not FASTA.
    """
    import http.client
    import shutil
    import urllib.error
    import urllib.request

    fasta_folder = _cache_dir()
    fasta_file = os.path.join(fasta_folder, identifier + '.fasta')

    # A versioned accession (AF228364.1) never changes sequence, so its cached
    # copy is valid for any source. An unversioned one (AF228364) means "the
    # current version", which can differ between sources and snapshots, so
    # with --local it is always resolved in the selected database (cheap).
    versioned = re.search(r'\.\d+$', identifier) is not None
    if _is_fasta_file(fasta_file):
        if versioned or not config.BLAST_LOCAL_DB:
            return fasta_file
    elif os.path.exists(fasta_file):
        os.remove(fasta_file)  # stale error reply or empty file from an earlier run

    if config.BLAST_LOCAL_DB:
        # --local: read the record from the local database; fall back to NCBI
        # only for real-looking accessions the local database lacks (nt omits
        # some divisions, e.g. most WGS contigs).
        if _fetch_local(identifier, fasta_file, fasta_folder):
            return fasta_file
        if not _ACCESSION_RE.match(identifier):
            raise AccessionNotFoundError(
                f"{identifier} is not an accession and is not in the local database")
        if not allow_ncbi:
            raise AccessionNotFoundError(f"{identifier} is not in the local database")

    from urllib.parse import urlencode
    params = {'db': 'nuccore', 'id': identifier, 'rettype': 'fasta', 'retmode': 'text',
              'tool': config.NCBI_TOOL}
    if config.NCBI_EMAIL:
        params['email'] = config.NCBI_EMAIL
    if config.NCBI_API_KEY:
        # Rate limit then applies to the key, not a shared (cluster) IP
        params['api_key'] = config.NCBI_API_KEY
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(params)
    request = urllib.request.Request(url, headers={'User-Agent': config.NCBI_TOOL})
    last_error = None
    for attempt in range(1, max_retries + 1):
        time.sleep(config.NCBI_REQUEST_DELAY_WITH_KEY if config.NCBI_API_KEY
                   else config.NCBI_REQUEST_DELAY)
        fd, tmp_path = _cache_tempfile(fasta_folder)
        os.close(fd)
        wait = 2 * attempt  # backoff before retry
        try:
            ok = False
            try:
                with urllib.request.urlopen(request, timeout=config.NCBI_TIMEOUT) as resp, \
                        open(tmp_path, 'wb') as out:
                    shutil.copyfileobj(resp, out)
                ok = True
            except urllib.error.HTTPError as e:
                # Keep the body: NCBI explains errors in it
                last_error = e
                if e.code == 429 or e.code >= 500:
                    # Rate-limited or overloaded (common from shared cluster
                    # IPs): honour Retry-After, else back off for longer.
                    try:
                        wait = min(float(e.headers.get('Retry-After')), 120)
                    except (TypeError, ValueError, AttributeError):
                        wait = 10 * attempt + random.uniform(0, 3)
                try:
                    with open(tmp_path, 'wb') as out:
                        out.write(e.read())
                except (http.client.HTTPException, OSError):
                    pass
                if e.code == 400 and config.NCBI_API_KEY and 'API key' in _read_head(tmp_path):
                    # A bad key makes every request fail: stop rather than
                    # report every sequence as missing. Never echo the key.
                    raise NCBIUnavailableError(
                        "NCBI rejected the API key in NCBI_API_KEY (invalid). "
                        "Fix or unset NCBI_API_KEY")
            except (http.client.HTTPException, OSError, ValueError) as e:
                # OSError: URLError, timeouts, connection resets.
                # HTTPException: IncompleteRead, RemoteDisconnected.
                if _is_cache_write_error(e):
                    # Disk full or over quota while saving: not a network
                    # problem, and retrying or skipping would drop sequences
                    raise _cache_error(e, fasta_folder) from e
                last_error = e

            if ok and _is_fasta_file(tmp_path):
                _cache_store(tmp_path, fasta_file)  # atomic: cache only valid FASTA
                return fasta_file
            if ok and os.path.getsize(tmp_path) > 0:
                # NCBI answers an unknown id with HTTP 200 and an 'Error: ...'
                # body: a definitive answer, retrying will not help. An error
                # status (429, 5xx) is never read this way, even if its body
                # starts with 'Error', so it gets the retries below.
                head = _read_head(tmp_path).strip()
                if head.startswith('Error'):
                    raise AccessionNotFoundError(f"NCBI has no record for {identifier}")
            if ok:
                last_error = 'unexpected reply (not FASTA)'
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        if attempt < max_retries:
            time.sleep(wait)
    if _ncbi_refused(last_error):
        # NCBI is up but refusing (rate limit / overload). Treating this as a
        # missing record would remove valid sequences: stop the run instead.
        raise NCBIUnavailableError(
            f"NCBI kept refusing requests for {identifier} ({last_error}) after {max_retries} "
            f"attempts. Try again later, set NCBI_API_KEY, or use --local")
    if _is_network_error(last_error):
        # NCBI could not be reached at all (offline node, DNS, firewall,
        # timeouts). That says nothing about the sequence: treating it as
        # "not found" would remove every sequence in the file.
        where = " (it is not in the local database)" if config.BLAST_LOCAL_DB else ""
        raise NCBIUnavailableError(
            f"Could not reach NCBI to fetch {identifier}{where} ({last_error}) after "
            f"{max_retries} attempts. Check the network connection and try again")
    raise RuntimeError(f"Failed to download FASTA for {identifier} after {max_retries} attempts"
                       f" ({last_error})")


def get_accession_version(fasta_file):
    """Get the version suffix (e.g. '.2') from the FASTA header.

    Returns '' when the header has no version, e.g. a record in a custom
    local BLAST database named Homo_sapiens_chr12: inventing '.1' would
    produce a name that exists nowhere and fails validation. NCBI headers
    are always versioned.
    """
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line.split()[0][1:]  # Remove '>'
                m = re.search(r'(\.\d+)$', name)
                return m.group(1) if m else ''
    return ''


def _hit_accession(alignment):
    """Accession of a BLAST hit, or None if it cannot be read."""
    # Format: "gi|xxx|ref|NZ_CP045771.1|" or "gi|xxx|gb|CP041614.1|"
    # (web API), or "ref|NZ_CP045771.1|" (local BLAST+ databases);
    # tpg/tpe/tpd are third-party annotation records in nt
    accession_match = re.search(r'(?:^|\|)(ref|gb|emb|dbj|tpg|tpe|tpd)\|([^|]+)\|', alignment.title)
    hit_id = alignment.hit_id
    if accession_match:
        return accession_match.group(2).rstrip('|')
    if config.BLAST_LOCAL_DB and (hit_id.startswith('lcl|') or '|' not in hit_id):
        # Custom database built with makeblastdb -parse_seqids: the hit id is
        # the record's own name (Homo_sapiens_chr12), which blastdbcmd can
        # retrieve for validation
        return hit_id[4:] if hit_id.startswith('lcl|') else hit_id
    # Try format: "NZ_CP045771.1 description"
    accession_match = re.search(r'^([A-Z]{1,2}_?[A-Z]*\d+\.\d+)', hit_id)
    if accession_match:
        return accession_match.group(1).rstrip('|')
    return None


def blast_search(sequence, verbose=False, min_identity=None, min_coverage=None, max_evalue=None,
                 exact_only=False, prefer_accession=None, prefer_start=None):
    """
    Perform BLAST search for a sequence that couldn't be found directly.

    Args:
        sequence: DNA/RNA sequence string (without gaps)
        verbose: Print progress
        min_identity: Minimum identity percentage (default from config)
        min_coverage: Minimum query coverage percentage (default from config)
        max_evalue: Maximum e-value to accept (default from config)
        exact_only: Accept only an exact full-length match (all bases
            identical, no gaps, whole query). The best inexact hit is still
            reported, but not returned.
        prefer_accession: The accession the sequence was named with, if any
            (e.g. CP176678 from NZ_CP176678.1_2667588). An exact match from
            this record (same accession, ignoring version and NZ_ prefix) is
            chosen over exact matches from other records, which are often
            other species with an identical sequence. It never wins over an
            exact match with an inexact one.
        prefer_start: The number in a name such as NZ_CP176678.1_2667588.
            Among exact matches from prefer_accession, one whose start (as
            written in the name, so the higher coordinate on the reverse
            strand) is exactly this number is chosen; otherwise the first.

    Returns:
        tuple: (accession_with_coords, is_accurate) or (None, False) if no good hit
    """
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio.Seq import Seq

    # Use config defaults if not specified
    if min_identity is None:
        min_identity = config.BLAST_MIN_IDENTITY
    if min_coverage is None:
        min_coverage = config.BLAST_MIN_COVERAGE
    if max_evalue is None:
        max_evalue = config.BLAST_MAX_EVALUE

    # Convert RNA to DNA for BLAST
    dna_seq = str(Seq(sequence).back_transcribe())

    if verbose:
        print(f"    Running BLAST search ({len(dna_seq)} bp)...")

    result_handle = None
    blast_records = None

    try:
        if config.BLAST_LOCAL_DB:
            # Local BLAST+ with the same task and scoring as the web search
            # below, so a local run and an API run differ only in database.
            gap_open, gap_extend = config.BLAST_GAP_COSTS.split()
            proc = subprocess.run(
                ['blastn', '-task', 'blastn', '-db', config.BLAST_LOCAL_DB,
                 '-word_size', str(config.BLAST_WORD_SIZE),
                 '-reward', str(config.BLAST_MATCH_REWARD),
                 '-penalty', str(config.BLAST_MISMATCH_PENALTY),
                 '-gapopen', gap_open, '-gapextend', gap_extend,
                 '-max_target_seqs', str(config.BLAST_HITLIST_SIZE),
                 '-evalue', str(max_evalue),
                 '-num_threads', str(config.BLAST_NUM_THREADS),
                 '-outfmt', '5'],
                input=f">query\n{dna_seq}\n", capture_output=True, text=True
            )
            if proc.returncode != 0:
                # A broken local setup is not "no hit": stop rather than let
                # every sequence be removed.
                raise LocalBlastError(f"local blastn failed: {proc.stderr.strip()}")
            result_handle = io.StringIO(proc.stdout)
        else:
            # Identify ourselves to NCBI (sent with every qblast request)
            NCBIWWW.tool = config.NCBI_TOOL
            NCBIWWW.email = config.NCBI_EMAIL

            # Standard blastn (not megablast) against nt, with the task's scoring
            # parameters set explicitly rather than relying on NCBI's defaults.
            # MEGABLAST is deliberately not sent. Limit results for speed.
            # Note: NCBI's server currently searches core_nt when asked for nt.
            # NCBI may throttle (429) or be overloaded (5xx): retry with a
            # long backoff, then stop the run rather than report "no hit".
            for attempt in range(1, 4):
                try:
                    result_handle = NCBIWWW.qblast(
                        "blastn", "nt", dna_seq,
                        word_size=config.BLAST_WORD_SIZE,
                        nucl_reward=config.BLAST_MATCH_REWARD,
                        nucl_penalty=config.BLAST_MISMATCH_PENALTY,
                        gapcosts=config.BLAST_GAP_COSTS,
                        hitlist_size=config.BLAST_HITLIST_SIZE,
                        expect=max_evalue,
                        format_type="XML"
                    )
                    break
                except Exception as e:
                    if not (_ncbi_refused(e) or _is_network_error(e)):
                        raise
                    if attempt == 3:
                        raise NCBIUnavailableError(
                            f"NCBI BLAST kept refusing or failing requests ({e}) after 3 attempts. "
                            f"Try again later, or use --local") from e
                    if verbose:
                        print(f"    NCBI BLAST request failed ({e}), retrying...")
                    time.sleep(config.NCBI_BLAST_RETRY_WAIT * attempt)
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)

        # Look at every hit returned (up to BLAST_HITLIST_SIZE), not just the
        # first. Choose, in order:
        #   1. an exact match from prefer_accession (the submitter's record):
        #      identical sequences in other species are common, and BLAST's
        #      order among equally good hits says nothing about which is right;
        #   2. the first exact full-length match (all bases identical, no
        #      gaps, whole query): only exact matches survive the later
        #      NCBI/MD5 check, and the exact one is not always listed first;
        #   3. unless exact_only, the first acceptable hit in BLAST's order.
        # Within that record, an exact match starting exactly at prefer_start
        # is used if there is one (the copy the submitter named), otherwise
        # the first in BLAST's order.
        prefer_core = _accession_core(prefer_accession) if prefer_accession else None
        first_ok = None
        first_exact = None
        preferred = None
        preferred_at_start = None
        n_exact = 0
        for alignment in blast_record.alignments:
            accession = _hit_accession(alignment)
            if not accession:
                continue
            for hsp in alignment.hsps:
                # Calculate identity and coverage. BLAST query coordinates are
                # 1-based and inclusive, so the span covers end - start + 1 bases.
                identity = (hsp.identities / hsp.align_length) * 100
                coverage = (hsp.query_end - hsp.query_start + 1) / len(dna_seq) * 100
                if not (identity >= min_identity and coverage >= min_coverage
                        and hsp.expect <= max_evalue):
                    continue

                # Coordinates follow the strand (start > end on the reverse strand)
                hit = (f"{accession}/{hsp.sbjct_start}-{hsp.sbjct_end}",
                       identity, coverage, hsp.expect)
                if first_ok is None:
                    first_ok = hit
                if hsp.identities == len(dna_seq) == hsp.align_length:
                    n_exact += 1
                    if first_exact is None:
                        first_exact = hit
                    if prefer_core and _accession_core(accession) == prefer_core:
                        if preferred is None:
                            preferred = hit
                        if (preferred_at_start is None and prefer_start is not None
                                and hsp.sbjct_start == prefer_start):
                            preferred_at_start = hit

        preferred = preferred_at_start or preferred
        chosen = preferred or first_exact or (None if exact_only else first_ok)
        if chosen:
            new_accession, identity, coverage, evalue = chosen
            if verbose:
                print(f"    BLAST hit: {new_accession} (identity={identity:.1f}%, coverage={coverage:.1f}%, e={evalue:.2e})")
                if preferred:
                    where = (f" starting at {prefer_start}" if preferred is preferred_at_start
                             else "")
                    print(f"    Chosen as an exact match from the named accession {prefer_core}"
                          f"{where} ({n_exact} exact match(es) in total)")
                elif prefer_core and first_exact:
                    print(f"    No exact match from the named accession {prefer_core} among the "
                          f"{n_exact} exact match(es); using the first one")
            return new_accession, True

        if verbose:
            if first_ok:
                # exact_only: report the near miss so curators can see it
                new_accession, identity, coverage, evalue = first_ok
                print(f"    No exact BLAST match (best: {new_accession}, identity={identity:.1f}%, "
                      f"coverage={coverage:.1f}%)")
            else:
                print("    No accurate BLAST hit found")
        return None, False

    except RunAbortError:
        raise
    except Exception as e:
        if config.BLAST_LOCAL_DB:
            # Locally, a real "no hit" is valid XML with no hits and never
            # gets here; any exception (blastn missing, unreadable output) is
            # a failure, and must not be mistaken for "no hit".
            raise LocalBlastError(f"local BLAST failed: {type(e).__name__}: {e}") from e
        if _is_network_error(e):
            # Web BLAST could not reach NCBI: not "no hit", so stop rather
            # than remove the sequence
            raise NCBIUnavailableError(
                f"Could not reach NCBI BLAST ({e}). Check the network connection and "
                f"try again, or use --local") from e
        if verbose:
            print(f"    BLAST search failed: {e}")
        return None, False

    finally:
        # NCBIXML.parse is a generator and we return as soon as the first record
        # gives us a hit, leaving it suspended. Its eventual garbage collection
        # re-finalizes the XML parser, and when NCBI truncates a response that
        # surfaces as an "Exception ignored ... ExpatError" traceback that no
        # caller can catch. Close it here instead, where the torn tail is ours
        # to discard: the record we needed is already parsed.
        if blast_records is not None:
            try:
                blast_records.close()
            except Exception:
                pass
        if result_handle is not None:
            try:
                result_handle.close()
            except Exception:
                pass


def _locate_in_record(identifier, ungapped, allow_ncbi=True, start=None):
    """
    Find an ungapped sequence in a record.

    If start is given (the number in a name such as NZ_CP176678.1_2667588)
    and an exact match begins exactly there, on either strand, that match is
    used: a sequence can occur more than once in a genome, and this is the
    copy the submitter named. The number is never used otherwise (no nearest
    match): if no match begins there, the first occurrence is used, forward
    strand first, then reverse complement.

    Returns:
        'ACCESSION/START-END' (START > END on the reverse strand), or None if
        the record does not contain the sequence.
    Raises:
        whatever get_fasta_file raises (record missing, NCBI unreachable, ...)
    """
    from Bio import SeqIO
    from Bio.Seq import Seq

    fasta_file = get_fasta_file(identifier, allow_ncbi=allow_ncbi)
    fasta = SeqIO.read(fasta_file, "fasta")

    # Get version if not in the identifier
    version = get_accession_version(fasta_file) if '.' not in identifier else ''

    # Convert RNA to DNA for matching
    s1 = Seq(ungapped).back_transcribe()
    s2 = s1.reverse_complement()
    n = len(s1)
    name = f'{identifier}{version}'

    if start is not None:
        # Forward strand: START-(START+n-1); reverse strand: START-(START-n+1)
        if fasta.seq[start - 1:start - 1 + n] == s1:
            return f'{name}/{start}-{start + n - 1}'
        if start >= n and fasta.seq[start - n:start] == s2:
            return f'{name}/{start}-{start - n + 1}'

    pos = fasta.seq.find(s1)
    if pos != -1:
        return f'{name}/{pos + 1}-{pos + n}'
    pos = fasta.seq.find(s2)
    if pos != -1:
        return f'{name}/{pos + n}-{pos + 1}'
    return None


def _find_in_named_accession(seq_name, hint, ungapped, verbose=False):
    """
    Look for the sequence in the accession its name starts with (hint from
    _accession_hint), trying the GenBank form first (CP176678.1, then
    NZ_CP176678.1, then unversioned CP176678). Names that are themselves the
    accession were already looked up directly and are skipped. A number
    after the accession picks the copy only as in _locate_in_record.

    With --local only the local database is used: these lookups are a
    shortcut before BLAST, and must not make offline runs depend on NCBI.

    Returns:
        'ACCESSION/START-END' or None.
    Raises:
        RunAbortError (NCBI unreachable in web mode, broken local database,
        cache not writable), as for any other lookup.
    """
    if not hint:
        return None
    _, candidates, start = hint
    for candidate in candidates:
        if candidate == seq_name:
            continue
        try:
            found = _locate_in_record(candidate, ungapped, start=start,
                                      allow_ncbi=not config.BLAST_LOCAL_DB)
        except RunAbortError:
            raise
        except Exception as e:
            if verbose:
                print(f"  {seq_name}: named accession {candidate} not available ({e})")
            continue
        if found:
            return found
        if verbose:
            print(f"  {seq_name}: sequence not found in named accession {candidate}")
    return None


def fix_missing_coordinates(filepath, verbose=False, use_blast_fallback=True):
    """
    Fix missing coordinates by looking up sequences in NCBI.
    Falls back to BLAST search if direct lookup fails.

    Args:
        filepath: Path to Stockholm file
        verbose: Print progress messages
        use_blast_fallback: Use BLAST search for sequences not found directly

    Returns:
        tuple: (id_mapping dict, sequences_to_remove set)
    """
    from types import SimpleNamespace
    from scripts import parser

    mapping = {}
    to_remove = set()

    # Read rows with the validator's own parser, not Biopython's AlignIO, which
    # rejects files the validator accepts (a duplicate row, an SS_cons of the
    # wrong length) and would crash the run instead of fixing them.
    with open(filepath) as f:
        entries = parser.parse_stockholm_file(f.readlines())['sequence_entries']
    align = [SimpleNamespace(id=name, seq=seq) for name, seq in entries]

    # Records that already carry coordinates need no lookup, so they are settled
    # up front and left out of the progress total.
    pending = []
    for record in align:
        _, coords = parse_sequence_identifier(record.id)
        if coords is not None:
            mapping[record.id] = record.id
        else:
            pending.append(record)

    total = len(pending)
    show_bar = total > 0 and not verbose and sys.stderr.isatty()

    for position, record in enumerate(pending, start=1):
        if verbose:
            print(f"  [{position}/{total} {position / total * 100:.1f}%] {record.id}")
        elif show_bar:
            sys.stderr.write('\r  Resolving coordinates ' + _render_progress(position - 1, total))
            sys.stderr.flush()

        # Remove gaps for sequence matching
        ungapped = str(record.seq).replace('.', '').replace('-', '').upper()

        # RNAcentral ids (URS..., URS..._taxid): locate within the RNAcentral
        # sequence instead of NCBI. On failure fall through to BLAST.
        if rnacentral.is_rnacentral_id(record.id):
            try:
                new_accession = rnacentral.find_coordinates(record.id, ungapped)
            except rnacentral.RNAcentralError as e:
                # RNAcentral unreachable: not a reason to drop the sequence
                raise RNAcentralUnavailableError(f"{e}. Try again later") from e
            if new_accession:
                if verbose:
                    print(f"  {record.id} -> {new_accession} (RNAcentral)")
                mapping[record.id] = new_accession
                continue
            if verbose:
                print(f"  {record.id}: sequence not found in RNAcentral entry")

        # A name that starts with an accession followed by other text
        # (NZ_CP176678.1_2667588) is not an NCBI record: the accession in it
        # is looked up instead (below). Without --local the whole name is not
        # sent to NCBI (a wasted request that always fails); with --local it
        # is still tried, as a custom database may hold records named that
        # way, but its expected failure is not reported.
        hint = _accession_hint(record.id)
        named_only = hint is not None and not _ACCESSION_RE.match(record.id)

        try:
            if rnacentral.is_rnacentral_id(record.id):
                # Not an NCBI accession; skip straight to BLAST
                raise LookupError("RNAcentral id, NCBI lookup skipped")
            if named_only and not config.BLAST_LOCAL_DB:
                raise LookupError("not an accession, looking up the accession in the name")
            # Inside the try: a record that is missing or unreadable costs
            # only this sequence. Errors that say nothing about the sequence
            # (NCBI unreachable, cache not writable) stop the run instead, via
            # the RunAbortError clause below.
            new_accession = _locate_in_record(record.id, ungapped)
            if new_accession:
                if verbose:
                    print(f"  {record.id} -> {new_accession}")
                mapping[record.id] = new_accession
                continue

        except RunAbortError:
            raise  # broken local database or NCBI refusing: stop, don't drop the sequence
        except Exception as e:
            if verbose and not named_only:
                print(f"  Warning: Could not read fasta for {record.id}: {e}")

        # The name may start with an accession followed by other text
        # (NZ_CP176678.1_2667588): look for the sequence in that record before
        # BLAST, which may list identical copies from other species first, or
        # not list this record at all among BLAST_HITLIST_SIZE hits.
        new_accession = _find_in_named_accession(record.id, hint, ungapped, verbose)
        if new_accession:
            if verbose:
                print(f"  {record.id} -> {new_accession} (from the accession in its name)")
            mapping[record.id] = new_accession
            continue

        # Direct lookup failed, try BLAST if enabled
        if use_blast_fallback:
            if verbose:
                print(f"  {record.id} not found directly, trying BLAST...")

            # Only an exact match can pass the final validation (Rfam's MD5
            # check). Accepting a near miss here would rename the sequence,
            # fail validation, and repeat this same BLAST search for nothing.
            blast_result, is_accurate = blast_search(
                ungapped, verbose=verbose, exact_only=True,
                prefer_accession=hint[0] if hint else None,
                prefer_start=hint[2] if hint else None)

            if is_accurate and blast_result:
                if verbose:
                    print(f"  {record.id} -> {blast_result} (via BLAST)")
                mapping[record.id] = blast_result
                continue

        # Neither direct lookup nor BLAST worked - mark for removal
        if verbose:
            print(f"  {record.id} will be REMOVED (no exact match found)")
        to_remove.add(record.id)

    if show_bar:
        sys.stderr.write('\r  Resolving coordinates ' + _render_progress(total, total) + '\n')
        sys.stderr.flush()

    return mapping, to_remove


def _render_progress(done, total, width=30):
    """Render a text progress bar, e.g. '[####------] 12/40 (30.0%)'."""
    fraction = done / total if total else 1.0
    filled = int(round(fraction * width))
    bar = '#' * filled + '-' * (width - filled)
    return f'[{bar}] {done}/{total} ({fraction * 100:.1f}%)'


def check_overlap(start1, end1, start2, end2):
    """
    Check if two coordinate ranges overlap by at least 1 bp.
    Handles both forward (start < end) and reverse (start > end) strands.
    """
    # Normalize coordinates so min is always first
    a1, a2 = min(start1, end1), max(start1, end1)
    b1, b2 = min(start2, end2), max(start2, end2)
    # Check for overlap
    return a1 <= b2 and b1 <= a2


def find_overlapping_sequences(sequence_entries, preferred=None):
    """
    Find sequences from the same accession that overlap.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples
        preferred: Optional set of names to keep first when they overlap
            others (e.g. those that match NCBI). Without it, or among equally
            preferred sequences, the one with the lower start is kept.

    Returns:
        set: Sequence names to remove due to overlaps
    """
    from collections import defaultdict
    preferred = preferred or set()

    # Group sequences by accession
    accession_groups = defaultdict(list)
    for seq_name, seq_data in sequence_entries:
        accession, coords = parse_sequence_identifier(seq_name)
        if coords:
            # Parse coordinates (handle both start-end formats)
            match = re.match(r'(\d+)-(\d+)', coords)
            if match:
                start, end = int(match.group(1)), int(match.group(2))
                accession_groups[accession].append((seq_name, start, end, seq_data))

    to_remove = set()

    # Check each accession group for overlaps
    for accession, sequences in accession_groups.items():
        if len(sequences) < 2:
            continue

        # Greedy: keep preferred sequences first, then by lower start (using
        # min of start/end for reverse strand); drop any that overlap a kept one
        sorted_seqs = sorted(sequences, key=lambda x: (x[0] not in preferred, min(x[1], x[2])))

        kept = []
        for seq_name, start, end, seq_data in sorted_seqs:
            # Check if this sequence overlaps with any kept sequence
            overlaps = False
            for kept_name, kept_start, kept_end, kept_data in kept:
                if check_overlap(start, end, kept_start, kept_end):
                    overlaps = True
                    break

            if not overlaps:
                kept.append((seq_name, start, end, seq_data))
            else:
                to_remove.add(seq_name)

    return to_remove


def validate_sequences_against_ncbi(sequence_entries, verbose=False, use_blast_fallback=True):
    """
    Validate that sequences match their NCBI source at the given coordinates.
    When validation fails, optionally falls back to BLAST to find the correct
    accession and coordinates.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples
        verbose: Print progress messages
        use_blast_fallback: Use BLAST search for sequences not found or mismatched

    Returns:
        tuple: (invalid set, not_found set, mismatched set, blast_fixed dict)
            - invalid: all sequences that failed validation (after BLAST attempts)
            - not_found: sequences where accession wasn't found in NCBI (and BLAST failed)
            - mismatched: sequences where data doesn't match NCBI (and BLAST failed)
            - blast_fixed: dict mapping old seq_name to new accession/coords from BLAST
    """
    from Bio import SeqIO
    from Bio.Seq import Seq

    not_found = set()
    mismatched = set()
    blast_fixed = {}

    # Counted over every entry so the total matches the "Validating N sequence(s)
    # against NCBI" line the caller prints; entries skipped below for lack of
    # usable coordinates are rare and resolve instantly.
    total = len(sequence_entries)
    show_bar = total > 0 and not verbose and sys.stderr.isatty()

    for position, (seq_name, seq_data) in enumerate(sequence_entries, start=1):
        if verbose:
            print(f"  [{position}/{total} {position / total * 100:.1f}%] {seq_name}")
        elif show_bar:
            sys.stderr.write('\r  Validating against NCBI ' + _render_progress(position - 1, total))
            sys.stderr.flush()

        accession, coords = parse_sequence_identifier(seq_name)

        if not coords:
            # No coordinates, can't validate
            continue

        # Parse coordinates
        match = re.match(r'(\d+)-(\d+)', coords)
        if not match:
            continue

        start, end = int(match.group(1)), int(match.group(2))

        # Remove gaps from Stockholm sequence (needed for both validation and BLAST)
        stockholm_seq = seq_data.replace('.', '').replace('-', '').upper()

        # 'not_found' or 'mismatch' when the reference check fails; used to
        # categorise the sequence if BLAST also fails.
        failure = None

        # Get FASTA file
        ncbi_ok = False
        if rnacentral.is_rnacentral_id(accession):
            # RNAcentral id: validate against RNAcentral instead of NCBI
            try:
                ncbi_ok, failure = rnacentral.validate_sequence(
                    accession, start, end, stockholm_seq, verbose=verbose)
            except rnacentral.RNAcentralError as e:
                raise RNAcentralUnavailableError(f"{e}. Try again later") from e
            if verbose:
                print(f"  {seq_name}: " + ("validated OK (RNAcentral)" if ncbi_ok
                                          else f"RNAcentral validation failed ({failure})"))
        else:
            try:
                fasta_file = get_fasta_file(accession)
                fasta = SeqIO.read(fasta_file, "fasta")

                # Extract subsequence from NCBI (handle forward and reverse strand)
                # Convert to 0-based indexing
                if start < end:
                    # Forward strand
                    ncbi_subseq = str(fasta.seq[start - 1:end]).upper()
                else:
                    # Reverse strand
                    ncbi_subseq = str(fasta.seq[end - 1:start].reverse_complement()).upper()

                # Convert RNA to DNA for comparison (NCBI stores DNA)
                stockholm_seq_dna = str(Seq(stockholm_seq).back_transcribe())

                # Compare
                if stockholm_seq_dna == ncbi_subseq:
                    ncbi_ok = True
                    if verbose:
                        print(f"  {seq_name}: validated OK")
                else:
                    failure = 'mismatch'
                    if verbose:
                        print(f"  {seq_name}: MISMATCH with NCBI")
                        if len(stockholm_seq_dna) != len(ncbi_subseq):
                            print(f"    Length: Stockholm={len(stockholm_seq_dna)}, NCBI={len(ncbi_subseq)}")
                        else:
                            # Show first difference
                            for i, (a, b) in enumerate(zip(stockholm_seq_dna, ncbi_subseq)):
                                if a != b:
                                    print(f"    First difference at position {i+1}: Stockholm={a}, NCBI={b}")
                                    break

            except RunAbortError:
                raise  # broken local database or NCBI refusing: stop, don't drop the sequence
            except Exception as e:
                failure = 'not_found'
                if verbose:
                    print(f"  {seq_name}: accession NOT FOUND in NCBI ({e})")

        if ncbi_ok:
            continue

        # Validation failed — try BLAST fallback
        if use_blast_fallback:
            if verbose:
                print(f"    Trying BLAST fallback for {seq_name}...")
            hint = _accession_hint(seq_name)
            blast_result, is_accurate = blast_search(
                stockholm_seq, verbose=verbose, prefer_accession=hint[0] if hint else None)
            if is_accurate and blast_result:
                # Re-validate the Stockholm sequence against the BLAST hit's accession
                # to ensure the sequence data actually matches (required for MD5 check)
                blast_accession, blast_coords = parse_sequence_identifier(blast_result)
                blast_match = re.match(r'(\d+)-(\d+)', blast_coords) if blast_coords else None
                if blast_match:
                    blast_start, blast_end = int(blast_match.group(1)), int(blast_match.group(2))
                    try:
                        blast_fasta_file = get_fasta_file(blast_accession)
                        blast_fasta = SeqIO.read(blast_fasta_file, "fasta")
                        if blast_start < blast_end:
                            blast_ncbi_seq = str(blast_fasta.seq[blast_start - 1:blast_end]).upper()
                        else:
                            blast_ncbi_seq = str(blast_fasta.seq[blast_end - 1:blast_start].reverse_complement()).upper()
                        stockholm_seq_dna = str(Seq(stockholm_seq).back_transcribe())
                        if stockholm_seq_dna == blast_ncbi_seq:
                            if verbose:
                                print(f"    {seq_name} -> {blast_result} (via BLAST, re-validated OK)")
                            blast_fixed[seq_name] = blast_result
                            continue
                        else:
                            if verbose:
                                print(f"    BLAST hit {blast_result} found but sequence still does not match NCBI — removing")
                    except RunAbortError:
                        raise  # broken local database or NCBI refusing: stop, don't drop the sequence
                    except Exception:
                        if verbose:
                            print(f"    Could not re-validate BLAST hit {blast_result} — removing")

        # Neither direct validation nor BLAST succeeded. Categorise using the
        # reason recorded above (avoids a second lookup, and keeps RNAcentral
        # ids away from NCBI).
        if failure == 'mismatch':
            mismatched.add(seq_name)
        else:
            not_found.add(seq_name)

    if show_bar:
        sys.stderr.write('\r  Validating against NCBI ' + _render_progress(total, total) + '\n')
        sys.stderr.flush()

    invalid = not_found | mismatched
    return invalid, not_found, mismatched, blast_fixed


def find_duplicates_from_entries(sequence_entries, raw_entries=None):
    """
    Find duplicate sequences from a list of sequence entries.

    When raw_entries is provided (block-aware mode), duplicates are detected as
    the same (accession, coords, seq_data) appearing multiple times within the
    SAME block. Same name across different blocks is interleaved, not a duplicate.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples (merged if interleaved)
        raw_entries: Optional list of (seq_name, seq_data, block_index) from parser

    Returns:
        tuple: (unique_sequences_dict, duplicate_indices)
    """
    if raw_entries is not None:
        return _find_duplicates_block_aware(sequence_entries, raw_entries)

    # Legacy mode: same as before
    seen_keys = {}
    unique_sequences = {}
    duplicate_indices = []

    for idx, (seq_name, seq_data) in enumerate(sequence_entries):
        accession, coords = parse_sequence_identifier(seq_name)
        key = (accession, coords, seq_data)

        if key not in seen_keys:
            seen_keys[key] = idx
            unique_sequences[seq_name] = seq_data
        else:
            duplicate_indices.append(idx)

    return unique_sequences, duplicate_indices


def _find_duplicates_block_aware(sequence_entries, raw_entries):
    """
    Detect true duplicates: same (accession, coords, data) appearing multiple
    times within the SAME block. Repeated names across different blocks are
    interleaved format, not duplicates.
    """
    from collections import defaultdict

    # For each block, track (accession, coords, seq_data) keys
    per_block_seen = defaultdict(lambda: defaultdict(int))
    duplicate_names = set()

    for seq_name, seq_data, block_idx in raw_entries:
        accession, coords = parse_sequence_identifier(seq_name)
        key = (accession, coords, seq_data)
        per_block_seen[block_idx][key] += 1
        if per_block_seen[block_idx][key] > 1:
            duplicate_names.add(seq_name)

    # Build unique_sequences from merged entries
    unique_sequences = {}
    duplicate_indices = []

    for idx, (seq_name, seq_data) in enumerate(sequence_entries):
        unique_sequences[seq_name] = seq_data
        if seq_name in duplicate_names:
            duplicate_indices.append(idx)

    return unique_sequences, duplicate_indices


def remove_duplicates(lines):
    """
    Remove duplicate sequences from Stockholm file lines.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (corrected_lines, num_duplicates_removed)
    """
    # Separate different types of lines
    header_lines = []
    annotation_lines = []
    footer_lines = []
    sequence_entries = []  # List of (seq_name, seq_content) tuples
    
    sequences_started = False
    sequences_ended = False
    
    for line in lines:
        stripped = line.strip()
        
        if stripped == '//':
            sequences_ended = True
            footer_lines.append(line)
        elif sequences_ended:
            footer_lines.append(line)
        elif stripped.startswith('# STOCKHOLM'):
            header_lines.append(line)
            sequences_started = True
        elif stripped.startswith('#=GF') or stripped.startswith('#=GC') or stripped.startswith('#=GS') or stripped.startswith('#=GR'):
            annotation_lines.append(line)
        elif stripped.startswith('#'):
            if not sequences_started:
                header_lines.append(line)
            else:
                annotation_lines.append(line)
        elif stripped:
            # This is a sequence line
            parts = stripped.split(None, 1)
            if len(parts) >= 2:
                seq_name = parts[0]
                seq_content = parts[1].replace(' ', '')
                sequence_entries.append((seq_name, seq_content))
        else:
            # Empty line
            if not sequences_started:
                header_lines.append(line)
    
    # Find duplicates
    seen_keys = {}
    unique_entries = []
    duplicates_count = 0
    
    for seq_name, seq_content in sequence_entries:
        accession, coords = parse_sequence_identifier(seq_name)
        key = (accession, coords, seq_content)
        
        if key not in seen_keys:
            seen_keys[key] = True
            unique_entries.append((seq_name, seq_content))
        else:
            duplicates_count += 1
    
    # Reconstruct file
    corrected_lines = []
    
    # Add header lines
    corrected_lines.extend(header_lines)
    
    # Calculate max sequence name length for formatting
    max_name_len = max(len(name) for name, _ in unique_entries) if unique_entries else 30
    name_width = max(max_name_len + 2, 30)
    
    # Add unique sequences
    for seq_name, seq_content in unique_entries:
        corrected_lines.append(f"{seq_name:<{name_width}} {seq_content}\n")
    
    # Add annotation lines
    corrected_lines.extend(annotation_lines)
    
    # Add footer lines
    corrected_lines.extend(footer_lines)
    
    return corrected_lines, duplicates_count
