"""
Offline tests: --fix must not crash on files the validator accepts, and must
stop (not remove sequences) when a sequence source cannot be reached.
"""

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
