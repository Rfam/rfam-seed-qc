"""
Offline tests for how BLAST hits are accepted.

tests/fixtures/blast/*.xml are real BLAST+ outputs (blastn -outfmt 5, same
settings as the tool) for one SNORD19 sequence against a two-record database:
exact_hit.xml for the sequence itself, near_hit.xml for a copy with one base
changed (75/76 identical). blastn is replaced by a stub that returns them, so
no BLAST+ install or network is needed.
"""

import io
import sys
import urllib.error
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts import config, fixable_errors  # noqa: E402

BLAST = Path(__file__).resolve().parent / 'fixtures' / 'blast'
EXACT_QUERY = 'GAGGTATGAGGAAAATATCCAACAGTCTGATTGTCATTGAAGACTGGTAAGATCCAACTCTGACCTCAGCAAAACT'
NEAR_QUERY = 'GAGGTATGAGGAAAATATCCAACAGTCTGATTGTCATTGACGACTGGTAAGATCCAACTCTGACCTCAGCAAAACT'
HIT = 'AAYZ01025109.1/988-1063'


@pytest.fixture
def fake_blastn(monkeypatch):
    """Use local mode with a stub blastn that returns the given XML file."""
    monkeypatch.setattr(config, 'BLAST_LOCAL_DB', 'mini')
    calls = []

    def use(xml_name):
        xml = (BLAST / xml_name).read_text()

        def run(cmd, **kwargs):
            assert cmd[0] == 'blastn', cmd
            calls.append(cmd)
            return SimpleNamespace(returncode=0, stdout=xml, stderr='')
        monkeypatch.setattr(fixable_errors.subprocess, 'run', run)
        return calls
    return use


def test_near_miss_accepted_by_default(fake_blastn):
    """Validation's BLAST fallback keeps the previous behaviour (it re-checks
    the hit for an exact match itself)."""
    fake_blastn('near_hit.xml')
    assert fixable_errors.blast_search(NEAR_QUERY) == (HIT, True)


def test_near_miss_rejected_when_exact_only(fake_blastn, capsys):
    fake_blastn('near_hit.xml')
    assert fixable_errors.blast_search(NEAR_QUERY, verbose=True, exact_only=True) == (None, False)
    out = capsys.readouterr().out
    assert f'No exact BLAST match (best: {HIT}, identity=98.7%, coverage=100.0%)' in out


def test_exact_match_accepted_when_exact_only(fake_blastn):
    fake_blastn('exact_hit.xml')
    assert fixable_errors.blast_search(EXACT_QUERY, exact_only=True) == (HIT, True)


@pytest.mark.parametrize('xml_name, query, expected', [
    ('exact_hit.xml', EXACT_QUERY, HIT),     # exact match: coordinates added
    ('near_hit.xml', NEAR_QUERY, None),      # near miss: removed straight away
])
def test_missing_coordinates_need_an_exact_match(fake_blastn, monkeypatch, tmp_path,
                                                 xml_name, query, expected):
    """A row without coordinates gets them from BLAST only for an exact match.
    A near miss is removed here, instead of being renamed, failing validation
    and being searched again."""
    calls = fake_blastn(xml_name)

    def not_found(identifier, **kwargs):
        raise fixable_errors.AccessionNotFoundError(identifier)
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', not_found)

    sto = tmp_path / 'species.sto'
    sto.write_text(f'# STOCKHOLM 1.0\n\nOchotona_princeps {query.replace("T", "U")}\n//\n')
    mapping, to_remove = fixable_errors.fix_missing_coordinates(str(sto))

    assert len(calls) == 1
    if expected:
        assert mapping == {'Ochotona_princeps': expected} and not to_remove
    else:
        assert to_remove == {'Ochotona_princeps'} and 'Ochotona_princeps' not in mapping


# Web BLAST: a reply that is not usable BLAST XML says nothing about the
# sequence. It is retried (usually temporary), and if it never works the run
# stops (like local BLAST). It must never count as "no hit", which would
# remove the sequence.

@pytest.fixture
def fake_qblast(monkeypatch):
    """qblast returning the given replies in turn (an exception is raised)."""
    from Bio.Blast import NCBIWWW
    monkeypatch.setattr(config, 'BLAST_LOCAL_DB', None)
    monkeypatch.setattr(fixable_errors.time, 'sleep', lambda seconds: None)
    calls = []

    def use(*replies):
        def qblast(*args, **kwargs):
            reply = replies[min(len(calls), len(replies) - 1)]
            calls.append(reply)
            if isinstance(reply, Exception):
                raise reply
            return io.StringIO(reply)
        monkeypatch.setattr(NCBIWWW, 'qblast', qblast)
        return calls
    return use


def _no_hits_xml():
    xml = (BLAST / 'near_hit.xml').read_text()
    return (xml.split('<Iteration_hits>')[0] + '<Iteration_hits>\n</Iteration_hits>'
            + xml.split('</Iteration_hits>')[1])


EXACT_XML = (BLAST / 'exact_hit.xml').read_text()
TEMPORARY_FAILURES = {
    'html': ('<html><body>Service temporarily unavailable</body></html>',
             "not BLAST XML; starts with '<html><body>Service temporarily unavailable"),
    'empty': ('', "not BLAST XML; starts with ''"),
    'cut-off': (EXACT_XML[:1500], 'ExpatError'),
    'ncbi-error': (ValueError('Error message from NCBI: CPU usage limit was exceeded'),
                   'CPU usage limit was exceeded'),
    'http-503': (urllib.error.HTTPError('https://blast.ncbi.nlm.nih.gov', 503, 'Unavailable', {}, None),
                 'HTTP Error 503'),
}


@pytest.mark.parametrize('kind', TEMPORARY_FAILURES)
def test_temporary_failure_is_retried_then_succeeds(fake_qblast, kind):
    calls = fake_qblast(TEMPORARY_FAILURES[kind][0], EXACT_XML)
    assert fixable_errors.blast_search(EXACT_QUERY, exact_only=True) == (HIT, True)
    assert len(calls) == 2


@pytest.mark.parametrize('kind', TEMPORARY_FAILURES)
def test_persistent_failure_stops_the_run_after_3_attempts(fake_qblast, kind):
    reply, message = TEMPORARY_FAILURES[kind]
    calls = fake_qblast(reply)
    with pytest.raises(fixable_errors.NCBIUnavailableError) as e:
        fixable_errors.blast_search(EXACT_QUERY)
    assert 'NCBI BLAST failed 3 times' in str(e.value)
    assert message in str(e.value)
    assert len(calls) == 3


def test_http_403_stops_without_retrying(fake_qblast):
    calls = fake_qblast(urllib.error.HTTPError('https://blast.ncbi.nlm.nih.gov', 403,
                                               'Forbidden', {}, None))
    with pytest.raises(fixable_errors.NCBIUnavailableError) as e:
        fixable_errors.blast_search(EXACT_QUERY)
    assert 'NCBI BLAST returned an unusable reply (HTTP Error 403' in str(e.value)
    assert len(calls) == 1


def test_web_blast_with_no_hits_is_still_no_hit(fake_qblast):
    calls = fake_qblast(_no_hits_xml())
    assert fixable_errors.blast_search(EXACT_QUERY) == (None, False)
    assert len(calls) == 1


def test_web_blast_exact_hit(fake_qblast):
    fake_qblast(EXACT_XML)
    assert fixable_errors.blast_search(EXACT_QUERY, exact_only=True) == (HIT, True)
