"""
Offline tests: when a sequence name starts with an accession
(NZ_CP176678.1_2667588), that record is preferred over identical sequences
from other records, first by looking in the record directly, then among the
exact BLAST matches.

The BLAST output is tests/fixtures/blast/exact_hit.xml with exact hits from
CP176678.1 added after the others, so they are not BLAST's first choice: one
copy at 2001-2076 and a second copy, on the reverse strand, at 5076-5001.
"""

import re
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts import config, fixable_errors  # noqa: E402

BLAST = Path(__file__).resolve().parent / 'fixtures' / 'blast'
QUERY = 'GAGGTATGAGGAAAATATCCAACAGTCTGATTGTCATTGAAGACTGGTAAGATCCAACTCTGACCTCAGCAAAACT'
FIRST_EXACT = 'AAYZ01025109.1/988-1063'
SOURCE_HIT = 'CP176678.1/2001-2076'
NAME = 'NZ_CP176678.1_2667588'


def _two_exact_hits_xml():
    """exact_hit.xml plus two exact copies in CP176678.1."""
    xml = (BLAST / 'exact_hit.xml').read_text()
    first = re.search(r'<Hit>.*?</Hit>\n', xml, re.S).group(0)

    def copy(num, hit_from, hit_to):
        return (first.replace('<Hit_num>1</Hit_num>', f'<Hit_num>{num}</Hit_num>')
                .replace('gb|AAYZ01025109.1|', 'gb|CP176678.1|')
                .replace('<Hit_accession>AAYZ01025109', '<Hit_accession>CP176678')
                .replace('<Hsp_hit-from>988<', f'<Hsp_hit-from>{hit_from}<')
                .replace('<Hsp_hit-to>1063<', f'<Hsp_hit-to>{hit_to}<'))
    extra = copy(3, 2001, 2076) + copy(4, 5076, 5001)
    return xml.replace('</Iteration_hits>', extra + '</Iteration_hits>')


@pytest.fixture
def blastn(monkeypatch):
    """Local mode with a stub blastn returning the two-exact-hit XML."""
    monkeypatch.setattr(config, 'BLAST_LOCAL_DB', 'mini')
    xml = _two_exact_hits_xml()
    calls = []

    def run(cmd, **kwargs):
        assert cmd[0] == 'blastn', cmd
        calls.append(cmd)
        return SimpleNamespace(returncode=0, stdout=xml, stderr='')
    monkeypatch.setattr(fixable_errors.subprocess, 'run', run)
    return calls


@pytest.mark.parametrize('name, expected', [
    ('NZ_CP176678.1_2667588', ('CP176678', ['CP176678.1', 'NZ_CP176678.1', 'CP176678'], 2667588)),
    ('CP061243.1_2605127', ('CP061243', ['CP061243.1', 'CP061243'], 2605127)),
    ('CP061243.1/1-84', ('CP061243', ['CP061243.1', 'CP061243'], None)),
    ('CP061243.1_12_extra', ('CP061243', ['CP061243.1', 'CP061243'], None)),
    ('CP061243.1_0', ('CP061243', ['CP061243.1', 'CP061243'], None)),
    ('AC104446', ('AC104446', ['AC104446'], None)),
    ('NC_000913.3_100', ('NC_000913', ['NC_000913.3', 'NC_000913'], 100)),
    ('Ochotona_princeps', None),
    ('URS0000572555_562', None),
    ('CP0612435abc', None),
])
def test_accession_hint(name, expected):
    assert fixable_errors._accession_hint(name) == expected


def test_first_exact_hit_without_preference(blastn):
    assert fixable_errors.blast_search(QUERY, exact_only=True) == (FIRST_EXACT, True)


def test_exact_hit_from_named_accession_preferred(blastn, capsys):
    assert fixable_errors.blast_search(QUERY, verbose=True, exact_only=True,
                                       prefer_accession='NZ_CP176678.1') == (SOURCE_HIT, True)
    assert 'exact match from the named accession CP176678 (3 exact match(es)' in capsys.readouterr().out


@pytest.mark.parametrize('start, expected', [
    (5076, 'CP176678.1/5076-5001'),   # the second copy starts exactly there
    (None, SOURCE_HIT),               # no number: first copy
    (5000, SOURCE_HIT),               # no copy starts there: first, not nearest
])
def test_named_start_picks_the_copy(blastn, start, expected):
    assert fixable_errors.blast_search(QUERY, exact_only=True, prefer_accession='CP176678',
                                       prefer_start=start) == (expected, True)


def test_named_start_alone_does_not_prefer_other_records(blastn):
    """AAYZ01025109.1 starts at 988, but is not the named accession."""
    assert fixable_errors.blast_search(QUERY, exact_only=True, prefer_accession='CP176678',
                                       prefer_start=988) == (SOURCE_HIT, True)


def test_inexact_hit_from_named_accession_never_wins(blastn):
    """AC104446.2 is in the hits, but only as a partial match."""
    assert fixable_errors.blast_search(QUERY, prefer_accession='AC104446') == (FIRST_EXACT, True)


def _fasta(tmp_path, name, seq):
    path = tmp_path / f'{name}.fasta'
    path.write_text(f'>{name} test record\n{seq}\n')
    return str(path)


def _sto(tmp_path, name):
    sto = tmp_path / 'input.sto'
    sto.write_text(f'# STOCKHOLM 1.0\n\n{name} {QUERY.replace("T", "U")}\n//\n')
    return str(sto)


def test_missing_coordinates_found_in_named_accession(blastn, monkeypatch, tmp_path):
    """The named record contains the sequence: no BLAST search, and with
    --local only the local database is asked (no NCBI fallback)."""
    record = _fasta(tmp_path, 'CP176678.1', 'AAAAA' + QUERY + 'CCCCC')
    asked = []

    def get_fasta_file(identifier, allow_ncbi=True, **kwargs):
        asked.append((identifier, allow_ncbi))
        if identifier == 'CP176678.1':
            return record
        raise fixable_errors.AccessionNotFoundError(identifier)
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', get_fasta_file)

    mapping, to_remove = fixable_errors.fix_missing_coordinates(_sto(tmp_path, NAME))

    assert mapping == {NAME: 'CP176678.1/6-81'} and not to_remove
    assert asked == [(NAME, True), ('CP176678.1', False)]
    assert not blastn


def _two_copies(tmp_path, name='CP176678.1'):
    """Forward copy at 6-81, reverse-strand copy at 200-125."""
    from Bio.Seq import Seq
    seq = 'AAAAA' + QUERY + 'C' * 43 + str(Seq(QUERY).reverse_complement()) + 'CCCC'
    return _fasta(tmp_path, name, seq)


@pytest.mark.parametrize('name, expected', [
    ('NZ_CP176678.1_200', 'CP176678.1/200-125'),   # reverse copy starts at 200
    ('NZ_CP176678.1_6', 'CP176678.1/6-81'),        # forward copy starts at 6
    ('NZ_CP176678.1_150', 'CP176678.1/6-81'),      # neither: first, not nearest
    ('NZ_CP176678.1_199', 'CP176678.1/6-81'),      # off by one: not used
])
def test_number_in_name_picks_the_copy(blastn, monkeypatch, tmp_path, name, expected):
    record = _two_copies(tmp_path)

    def get_fasta_file(identifier, **kwargs):
        if identifier == 'CP176678.1':
            return record
        raise fixable_errors.AccessionNotFoundError(identifier)
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', get_fasta_file)

    mapping, _ = fixable_errors.fix_missing_coordinates(_sto(tmp_path, name))
    assert mapping == {name: expected}
    assert not blastn


def test_web_mode_does_not_look_up_the_whole_name(monkeypatch, tmp_path):
    """Without --local, NZ_CP176678.1_2667588 is never sent to NCBI."""
    monkeypatch.setattr(config, 'BLAST_LOCAL_DB', None)
    record = _fasta(tmp_path, 'CP176678.1', 'AAAAA' + QUERY)
    asked = []

    def get_fasta_file(identifier, allow_ncbi=True, **kwargs):
        asked.append((identifier, allow_ncbi))
        return record
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', get_fasta_file)

    mapping, _ = fixable_errors.fix_missing_coordinates(_sto(tmp_path, NAME))
    assert mapping == {NAME: 'CP176678.1/6-81'}
    assert asked == [('CP176678.1', True)]


def test_other_version_of_named_accession(blastn, monkeypatch, tmp_path):
    """The database holds CP176678.2, not the CP176678.1 in the name: the
    unversioned lookup finds it, and the sequence must still match exactly."""
    record = _fasta(tmp_path, 'CP176678.2', 'AAAAA' + QUERY)
    asked = []

    def get_fasta_file(identifier, **kwargs):
        asked.append(identifier)
        if identifier == 'CP176678':
            return record
        raise fixable_errors.AccessionNotFoundError(identifier)
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', get_fasta_file)

    mapping, _ = fixable_errors.fix_missing_coordinates(_sto(tmp_path, NAME))
    assert mapping == {NAME: 'CP176678.2/6-81'}
    assert asked == [NAME, 'CP176678.1', 'NZ_CP176678.1', 'CP176678']
    assert not blastn


def test_missing_coordinates_reverse_strand_in_named_accession(blastn, monkeypatch, tmp_path):
    from Bio.Seq import Seq
    record = _fasta(tmp_path, 'CP176678.1', 'AAAAA' + str(Seq(QUERY).reverse_complement()))

    def get_fasta_file(identifier, **kwargs):
        if identifier == 'CP176678.1':
            return record
        raise fixable_errors.AccessionNotFoundError(identifier)
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', get_fasta_file)

    mapping, _ = fixable_errors.fix_missing_coordinates(_sto(tmp_path, NAME))
    assert mapping == {NAME: 'CP176678.1/81-6'}
    assert not blastn


def test_missing_coordinates_blast_prefers_named_accession(blastn, monkeypatch, tmp_path):
    """The named record cannot be read: BLAST still prefers its exact hit."""
    def not_found(identifier, **kwargs):
        raise fixable_errors.AccessionNotFoundError(identifier)
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', not_found)

    mapping, to_remove = fixable_errors.fix_missing_coordinates(_sto(tmp_path, NAME))

    assert mapping == {NAME: SOURCE_HIT} and not to_remove
    assert len(blastn) == 1


def test_validation_blast_fallback_prefers_named_accession(blastn, monkeypatch, tmp_path):
    """A row whose coordinates are wrong is renamed to the same accession at
    the right place, not to another species' identical copy."""
    record = _fasta(tmp_path, 'CP176678.1', 'G' * 2000 + QUERY + 'GGGG')
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', lambda identifier, **kwargs: record)

    name = 'CP176678.1/1-76'
    invalid, _, _, blast_fixed = fixable_errors.validate_sequences_against_ncbi([(name, QUERY)])

    assert blast_fixed == {name: SOURCE_HIT} and not invalid
