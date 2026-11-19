"""
Tests for scripts/build_live_cm.py model extraction.

tests/fixtures/cm/RF00005_svn.CM is the real SVN CM file for RF00005 (tRNA),
with its COM (build command) lines removed: a CM block followed by its
filter HMM block, both named SEED, without GA/TC/NC. DESC holds the GA, TC
and NC lines of RF00005's SVN DESC file. The cmpress/cmscan test runs only
if Infernal is installed.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts.build_live_cm import extract_first_model, parse_thresholds  # noqa: E402

SVN_CM = (REPO / "tests" / "fixtures" / "cm" / "RF00005_svn.CM").read_text()
DESC = "ID   tRNA\nGA   29.00\nTC   29.00\nNC   28.90\nTP   Gene; tRNA;\n"
TRNA = "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGTCCTGTGTTCGATCCACAGAATTCGCACCA"


def _headers(model):
    return [line for line in model.splitlines()
            if line.startswith(("INFERNAL1/", "HMMER3/", "NAME", "//"))]


def test_keeps_cm_and_filter_hmm_renamed():
    model = extract_first_model(SVN_CM, "RF00005")
    assert [h.split()[0] for h in _headers(model)] == [
        "INFERNAL1/a", "NAME", "//", "HMMER3/f", "NAME", "//"]
    assert [h.split()[1] for h in _headers(model) if h.startswith("NAME")] == [
        "RF00005", "RF00005"]
    assert "SEED" not in model.split("\n", 2)[1]
    assert model.endswith("//\n")


def test_ignores_blocks_after_first_model():
    model = extract_first_model(SVN_CM + SVN_CM, "RF00005")
    assert model == extract_first_model(SVN_CM, "RF00005")


@pytest.mark.parametrize("text", [
    "",
    "INFERNAL1/a [1.1.4]\nNAME     SEED\n//\n",                       # CM only
    "HMMER3/f [3.3.2]\nNAME  SEED\n//\nINFERNAL1/a [1.1.4]\nNAME  SEED\n//\n",  # wrong order
    "INFERNAL1/a [1.1.4]\nNAME     SEED\n",                           # truncated
])
def test_rejects_incomplete_models(text):
    assert extract_first_model(text, "RF00005") is None


def test_parse_thresholds():
    assert parse_thresholds(DESC) == {"GA": 29.0, "TC": 29.0, "NC": 28.9}
    assert parse_thresholds("ID   x\nTC   29.00\n") is None   # no GA
    assert parse_thresholds("GA   high\n") is None
    assert parse_thresholds(None) is None


def test_inserts_thresholds_after_null_in_cm_block_only():
    model = extract_first_model(SVN_CM, "RF00005", parse_thresholds(DESC))
    cm, hmm = model.split("HMMER3/f", 1)
    lines = cm.splitlines()
    i = next(n for n, line in enumerate(lines) if line.startswith("NULL "))
    assert lines[i + 1:i + 4] == ["GA       29.00", "TC       29.00", "NC       28.90"]
    assert not any(line.startswith(("GA ", "TC ", "NC ")) for line in hmm.splitlines())
    # existing lines are replaced, not duplicated
    again = extract_first_model(model, "RF00005", {"GA": 30.0})
    assert [l for l in again.splitlines() if l.startswith(("GA ", "TC ", "NC "))] == [
        "GA       30.00"]


@pytest.mark.skipif(not (shutil.which("cmpress") and shutil.which("cmscan")),
                    reason="Infernal not installed")
def test_output_presses_and_scans(tmp_path):
    cm = tmp_path / "live.cm"
    cm.write_text(extract_first_model(SVN_CM, "RF00005", parse_thresholds(DESC)))
    subprocess.run(["cmpress", str(cm)], check=True, capture_output=True)
    query = tmp_path / "q.fa"
    query.write_text(f">q\n{TRNA}\n")
    tbl = tmp_path / "q.tbl"
    subprocess.run(["cmscan", "--cut_ga", "--tblout", str(tbl), str(cm), str(query)],
                   check=True, capture_output=True)
    hits = [line.split()[0] for line in tbl.read_text().splitlines()
            if not line.startswith("#")]
    assert hits == ["RF00005"]


# The validator must say when the known-family check could not run, even
# without -v: a silent failure looks exactly like "no known families matched".
import validate_stockholm  # noqa: E402

ENTRIES = [("q/1-76", TRNA)]


def test_cmscan_missing_is_reported(monkeypatch, capsys):
    monkeypatch.setattr(validate_stockholm.shutil, "which", lambda name: None)
    assert validate_stockholm.filter_known_families(ENTRIES, "db.cm") == []
    assert "known Rfam families were NOT checked" in capsys.readouterr().err


@pytest.mark.skipif(not shutil.which("cmscan"), reason="Infernal not installed")
def test_cmscan_failure_is_reported(tmp_path, capsys):
    # A database without GA thresholds: cmscan --cut_ga refuses it
    cm = tmp_path / "no_ga.cm"
    cm.write_text(extract_first_model(SVN_CM, "RF00005"))
    subprocess.run(["cmpress", str(cm)], check=True, capture_output=True)
    assert validate_stockholm.filter_known_families(ENTRIES, str(cm)) == []
    err = capsys.readouterr().err
    assert "cmscan failed, known Rfam families were NOT checked" in err
    assert "GA bit threshold unavailable" in err


def test_cmscan_exception_is_reported(monkeypatch, capsys):
    monkeypatch.setattr(validate_stockholm.shutil, "which", lambda name: "/bin/cmscan")

    def boom(*args, **kwargs):
        raise OSError("cannot run")
    monkeypatch.setattr(validate_stockholm.subprocess, "Popen", boom)
    assert validate_stockholm.filter_known_families(ENTRIES, "db.cm") == []
    assert "cmscan filtering failed, known Rfam families were NOT checked" in capsys.readouterr().err
