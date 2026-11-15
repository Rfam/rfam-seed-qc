"""
Offline tests for how BLAST hits are accepted.

tests/fixtures/blast/*.xml are real BLAST+ outputs (blastn -outfmt 5, same
settings as the tool) for one SNORD19 sequence against a two-record database:
exact_hit.xml for the sequence itself, near_hit.xml for a copy with one base
changed (75/76 identical). blastn is replaced by a stub that returns them, so
no BLAST+ install or network is needed.
"""

import sys
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

    def not_found(identifier, max_retries=3):
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
