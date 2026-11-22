"""
Offline tests: --fix must not crash on files the validator accepts, and must
stop (not remove sequences) when a sequence source cannot be reached.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts import fixable_errors, rnacentral  # noqa: E402

SEQ = 'CAUGGGUGAGGUACGAGG'


def test_missing_coordinates_with_duplicate_row_and_bad_ss_cons(tmp_path, monkeypatch):
    """A row without coordinates, a duplicated row and an SS_cons of the wrong
    length: all accepted by the validator (fixable / warning), and previously a
    crash here because the file was re-read with Biopython's stricter AlignIO."""
    fasta = tmp_path / 'record.fasta'
    fasta.write_text('>AC104446.2 test record\nTTTT' + SEQ.replace('U', 'T') + 'TTTT\n')
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', lambda identifier, **kwargs: str(fasta))

    sto = tmp_path / 'input.sto'
    sto.write_text('# STOCKHOLM 1.0\n\n'
                   f'AC104446.2   {SEQ}\n'
                   f'X01.1/1-18   {SEQ}\n'
                   f'X01.1/1-18   {SEQ}\n'
                   '#=GC SS_cons <<<>>>\n//\n')
    mapping, to_remove = fixable_errors.fix_missing_coordinates(str(sto))
    assert mapping['AC104446.2'] == 'AC104446.2/5-22'
    assert not to_remove


def test_rnacentral_outage_stops_validation(monkeypatch):
    """RNAcentral unreachable must stop the run, not count the sequence as
    "not found" and send it to BLAST, which could remove it."""
    def down(urs):
        raise rnacentral.RNAcentralError('RNAcentral request failed: connection refused')
    monkeypatch.setattr(rnacentral, 'fetch_entry', down)
    with pytest.raises(fixable_errors.RNAcentralUnavailableError):
        fixable_errors.validate_sequences_against_ncbi(
            [('URS0000572555_562/1-18', SEQ)], use_blast_fallback=True)


def test_rnacentral_outage_stops_coordinate_fixing(tmp_path, monkeypatch):
    def down(urs):
        raise rnacentral.RNAcentralError('RNAcentral request failed: connection refused')
    monkeypatch.setattr(rnacentral, 'fetch_entry', down)
    sto = tmp_path / 'input.sto'
    sto.write_text(f'# STOCKHOLM 1.0\n\nURS0000572555_562 {SEQ}\n//\n')
    with pytest.raises(fixable_errors.RNAcentralUnavailableError):
        fixable_errors.fix_missing_coordinates(str(sto))


@pytest.mark.parametrize('header', ['# STOCKHOLM', '# STOCKHOLM 2.0', '# STOCKHOLM 1.0'])
def test_fix_writes_stockholm_1_0_header(tmp_path, monkeypatch, header):
    """A missing or other version is only a warning; --fix still runs and
    writes a header Infernal accepts."""
    import validate_stockholm
    monkeypatch.setattr(fixable_errors, 'validate_sequences_against_ncbi',
                        lambda seqs, **kwargs: (set(), set(), set(), {}))
    sto = tmp_path / 'input.sto'
    sto.write_text(f'{header}\n#=GF ID   test\n\n'
                   f'X01.1/1-18   {SEQ}\nX02.1/1-18   {SEQ}\n'
                   '#=GC SS_cons <<<..............>>>\n//\n')
    assert validate_stockholm.validate_stockholm_file(str(sto))['is_valid']
    success, _, output_path = validate_stockholm.fix_file(str(sto))
    assert success
    lines = Path(output_path).read_text().splitlines()
    assert lines[0] == '# STOCKHOLM 1.0'
    assert sum(line.startswith('# STOCKHOLM') for line in lines) == 1
    assert '#=GF ID   test' in lines


# Repairs can give rows the same name, or overlapping coordinates, after the
# duplicate and overlap checks: keep exactly one, never both or neither.

def _fix(tmp_path, monkeypatch, text, mapping, blast_fixed):
    import validate_stockholm
    monkeypatch.setattr(fixable_errors, 'fix_missing_coordinates',
                        lambda *args, **kwargs: (dict(mapping), set()))
    monkeypatch.setattr(fixable_errors, 'validate_sequences_against_ncbi',
                        lambda seqs, **kwargs: (set(), set(), set(), dict(blast_fixed)))
    sto = tmp_path / 'input.sto'
    sto.write_text(text)
    success, _, output_path = validate_stockholm.fix_file(str(sto))
    assert success
    lines = Path(output_path).read_text().splitlines()
    rows = [line.split()[0] for line in lines if line and not line.startswith(('#', '//'))]
    return rows, lines


def test_rows_given_the_same_name_by_coordinate_repair_keep_one(tmp_path, monkeypatch, capsys):
    # same record and coordinates, different gaps: overlap removal used to drop both
    rows, lines = _fix(tmp_path, monkeypatch,
                       '# STOCKHOLM 1.0\n\n'
                       f'ACC.1        {SEQ}--\n#=GR ACC.1 PP {"9" * 20}\n'
                       f'ACC.1_copy   -{SEQ}-\n#=GR ACC.1_copy PP {"8" * 20}\n'
                       f'OTHER.1/1-18 {SEQ}--\n//\n',
                       {'ACC.1': 'ACC.1/5-22', 'ACC.1_copy': 'ACC.1/5-22'}, {})
    assert rows == ['ACC.1/5-22', 'OTHER.1/1-18']
    gr = [line for line in lines if line.startswith('#=GR')]
    # the kept row's annotation only (all-gap columns are trimmed, so 18 left)
    assert len(gr) == 1 and gr[0].endswith('9' * 18) and '8' not in gr[0]
    assert 'same name as another row when coordinates were added' in capsys.readouterr().out


def test_blast_rescues_with_the_same_name_keep_one(tmp_path, monkeypatch, capsys):
    rows, _ = _fix(tmp_path, monkeypatch,
                   f'# STOCKHOLM 1.0\n\nBAD1.1/1-18  {SEQ}\nBAD2.1/1-18  {SEQ}\nOTHER.1/1-18 {SEQ}\n//\n',
                   {}, {'BAD1.1/1-18': 'ACC.1/5-22', 'BAD2.1/1-18': 'ACC.1/5-22'})
    assert rows == ['ACC.1/5-22', 'OTHER.1/1-18']
    assert 'BLAST placed at the same location as, or overlapping, another row' in capsys.readouterr().out


def test_blast_rescue_overlapping_a_kept_row_is_removed(tmp_path, monkeypatch):
    # the row checked under its own name wins over the BLAST rescue, even with a higher start
    rows, _ = _fix(tmp_path, monkeypatch,
                   f'# STOCKHOLM 1.0\n\nACC.1/10-27  {SEQ}\nBAD1.1/1-18  {SEQ}\n//\n',
                   {}, {'BAD1.1/1-18': 'ACC.1/5-22'})
    assert rows == ['ACC.1/10-27']


def test_blast_rescues_that_do_not_collide_are_all_kept(tmp_path, monkeypatch):
    rows, _ = _fix(tmp_path, monkeypatch,
                   f'# STOCKHOLM 1.0\n\nBAD1.1/1-18  {SEQ}\nBAD2.1/1-18  {SEQ}\n//\n',
                   {}, {'BAD1.1/1-18': 'ACC.1/5-22', 'BAD2.1/1-18': 'ACC.1/105-122'})
    assert rows == ['ACC.1/5-22', 'ACC.1/105-122']


def test_gs_lines_follow_renamed_sequences(tmp_path, monkeypatch):
    """#=GS lines must use the final name (esl-reformat rejects the old one)."""
    rows, lines = _fix(tmp_path, monkeypatch,
                       '# STOCKHOLM 1.0\n\n#=GS ACC.1 DE coordinates added\n'
                       '#=GS BAD1.1/1-18 DE renamed by BLAST\n#=GS OTHER.1/1-18 DE unchanged\n\n'
                       f'ACC.1        {SEQ}\nBAD1.1/1-18  {SEQ}\nOTHER.1/1-18 {SEQ}\n//\n',
                       {'ACC.1': 'ACC.1/5-22'}, {'BAD1.1/1-18': 'NEW.1/5-22'})
    assert rows == ['ACC.1/5-22', 'NEW.1/5-22', 'OTHER.1/1-18']
    assert [' '.join(line.split()) for line in lines if line.startswith('#=GS')] == [
        '#=GS ACC.1/5-22 DE coordinates added',
        '#=GS NEW.1/5-22 DE renamed by BLAST',
        '#=GS OTHER.1/1-18 DE unchanged']
    if shutil.which('esl-reformat'):
        out = tmp_path / 'input_corrected.sto'
        assert subprocess.run(['esl-reformat', 'stockholm', str(out)],
                              capture_output=True).returncode == 0


@pytest.mark.parametrize('name, seq, ok', [
    ('REF.1/1-4', 'ACGU', True),
    ('REF.1/4-1', 'ACGU', True),      # reverse strand: reverse complement of ACGT is ACGT
    ('REF.1/1-1', 'A', True),         # a single base is read forward, not complemented
    ('REF.1/2-3', 'CG', True),
    ('REF.1/1-999', 'ACGU', False),   # past the end of the 4-base record: not cut short
    ('REF.1/999-1', 'ACGU', False),
])
def test_source_check_coordinates(tmp_path, monkeypatch, name, seq, ok):
    fasta = tmp_path / 'ref.fasta'
    fasta.write_text('>REF.1 test\nACGT\n')
    monkeypatch.setattr(fixable_errors, 'get_fasta_file', lambda identifier, **kwargs: str(fasta))
    invalid, _, mismatched, _ = fixable_errors.validate_sequences_against_ncbi(
        [(name, seq)], use_blast_fallback=False)
    assert (name not in invalid) == ok
    if not ok:
        assert name in mismatched


@pytest.mark.parametrize('coords, expected', [
    ('1-74', (1, 74)), ('74-1', (74, 1)), ('5-5', (5, 5)),
    ('not-coordinates', None), ('1-4x', None), ('0-3', None), ('3-0', None),
    ('-1-4', None), ('1-', None), ('', None), (None, None),
])
def test_parse_coordinates(coords, expected):
    assert fixable_errors.parse_coordinates(coords) == expected


@pytest.mark.parametrize('verbose', [[], ['-v']])
def test_output_mode_stdout_is_a_clean_alignment(tmp_path, monkeypatch, capsys, verbose):
    """--fix --output-mode stdout: stdout holds only the corrected alignment
    (it can be piped); status lines and warnings go to stderr."""
    import validate_stockholm
    monkeypatch.setattr(fixable_errors, 'validate_sequences_against_ncbi',
                        lambda seqs, **kwargs: (set(), set(), set(), {}))
    monkeypatch.setattr(sys, 'stdout', sys.stdout)   # restored after main() redirects it
    sto = tmp_path / 'input.sto'
    sto.write_text(f'# STOCKHOLM 1.0\n\nX01.1/1-18   {SEQ}\nX02.1/1-18   {SEQ}\n'
                   '#=GC SS_cons <<<............>>>\n//\n')    # 18 columns, as the sequences
    monkeypatch.setattr(sys, 'argv', ['rfam-seed-qc', *verbose, '--fix', '--output-mode', 'stdout', str(sto)])
    assert validate_stockholm.main() == 0
    out, err = capsys.readouterr()
    assert out.startswith('# STOCKHOLM 1.0\n') and out.endswith('//\n')
    assert out.splitlines().count('//') == 1 and 'passed validation' not in out
    assert 'All 1 file(s) passed validation' in err
    assert not list(tmp_path.glob('*_corrected*')) and not list(tmp_path.glob('*_Report*'))
    if shutil.which('esl-reformat'):
        piped = tmp_path / 'piped.sto'
        piped.write_text(out)
        assert subprocess.run(['esl-reformat', 'stockholm', str(piped)],
                              capture_output=True).returncode == 0
