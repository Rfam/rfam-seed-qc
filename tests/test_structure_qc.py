"""
Tests for scripts/structure_qc.py: WUSS parsing (shared with the SS_cons
format check) and base-pair support warnings.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts import structure_qc as q  # noqa: E402
from scripts.stockholm_warnings import validate_ss_cons_format  # noqa: E402


# ---------------------------------------------------------------- parse_wuss

def test_parses_all_bracket_kinds():
    parsed = q.parse_wuss("<([{..}])>")
    assert parsed['pairs'] == [(0, 9), (1, 8), (2, 7), (3, 6)]
    assert not parsed['invalid'] and not parsed['unmatched_close'] and not parsed['unmatched_open']


def test_parses_pseudoknot_letters_beyond_d():
    # E..e and Z..z were reported as invalid before (only A-D were accepted)
    parsed = q.parse_wuss("<<EE..>>..ee:Z.z")
    assert parsed['pairs'] == [(0, 7), (1, 6), (2, 11), (3, 10), (13, 15)]
    assert validate_ss_cons_format("<<EE..>>..ee:Z.z") == []


def test_brackets_of_different_kinds_may_cross():
    # each kind is matched separately (as in Easel), like pseudoknot letters
    parsed = q.parse_wuss("<.(.>.)")
    assert parsed['pairs'] == [(0, 4), (2, 6)]
    assert validate_ss_cons_format("<.(.>.)") == []


def test_unpaired_characters():
    assert q.parse_wuss(":,_-.~")['pairs'] == []
    assert validate_ss_cons_format(":,_-.~") == []


@pytest.mark.parametrize("ss, message", [
    ("<<..>?", "invalid characters: '?' (position 6)"),
    ("<..>>", "unmatched '>' at position 5"),
    ("<<..>", "1 unmatched '<' (missing '>')"),
    ("EE..e", "1 unmatched 'E' (missing 'e')"),
    ("<<..>x", "unmatched 'x' at position 6"),    # x closes an X pseudoknot
])
def test_format_problems_are_reported(ss, message):
    assert any(message in w for w in validate_ss_cons_format(ss)), validate_ss_cons_format(ss)


# -------------------------------------------------------------- pair support

SS = "<<..>>"   # pairs 1-6 and 2-5


def entries(*seqs):
    return [(f"s{i}/1-6", s) for i, s in enumerate(seqs, 1)]


def test_fully_supported_structure_gives_no_warning():
    seqs = entries("GCAAGC", "GUAAAC", "CGUUCG", "UGAAUA")  # GC/CG/GU/UG/AU/UA
    assert q.check_pair_support(SS, seqs) == []


def test_t_is_read_as_u_and_case_is_ignored():
    assert q.check_pair_support(SS, entries("gtaaac", "GCAAGC")) == []


def test_weak_pair_is_reported():
    # pair 1-6 is G-C in one sequence and A-C, C-C in two others
    warnings = q.check_pair_support(SS, entries("GCAAGC", "ACAAGC", "CCAAGC"))
    assert len(warnings) == 1
    assert "1 of 2 SS_cons base pairs is Watson-Crick" in warnings[0]
    assert "1-6 (1/3 canonical)" in warnings[0]


def test_gapped_pair_is_reported():
    warnings = q.check_pair_support(SS, entries("GCAAGC", "-CAAG-", "GCAAG-"))
    assert any("has a gap" in w and "1-6 (2/3 gapped)" in w for w in warnings), warnings


def test_misaligned_sequence_is_reported():
    ss = "<<<..>>>"
    seqs = [("good/1-8", "GCGAACGC"), ("bad/1-8", "GCGAAGCG")]  # bad: G-G, C-C, G-G
    warnings = q.check_pair_support(ss, seqs)
    assert any("bad/1-8 (3/3)" in w and "possibly misaligned" in w for w in warnings), warnings
    assert not any("good/1-8" in w for w in warnings)


def test_one_odd_pair_does_not_flag_a_sequence():
    # 1 of 2 pairs non-canonical is above 30%, but a single pair is not a pattern
    warnings = q.check_pair_support(SS, entries("GCAAGC", "GCAAGC", "GCAAGC", "ACAAGC"))
    assert not any("misaligned" in w for w in warnings)


def test_ambiguous_bases_are_ignored():
    assert q.check_pair_support(SS, entries("NCAAGC", "GCAAGN", "GCAAGC")) == []


@pytest.mark.parametrize("ss", [None, "", "<<..>", "<<..>>x", "......", "<<..>>.."])
def test_unusable_or_pairless_ss_cons_is_skipped(ss):
    # missing, unbalanced, invalid, no pairs, or wrong length: nothing to check
    assert q.check_pair_support(ss, entries("ACAAGC", "CCAAGC")) == []


def test_long_lists_are_shortened():
    ss = "<" * 12 + ">" * 12
    seqs = [("s1/1-24", "A" * 24), ("s2/1-24", "C" * 24)]
    warnings = q.check_pair_support(ss, seqs)
    weak = next(w for w in warnings if "Watson-Crick" in w)
    assert "12 of 12 SS_cons base pairs are Watson-Crick" in weak and "... (2 more)" in weak
