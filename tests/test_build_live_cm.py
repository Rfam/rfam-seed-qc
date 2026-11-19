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

from scripts.build_live_cm import (  # noqa: E402
    extract_first_model, parse_family_id, parse_thresholds)

SVN_CM = (REPO / "tests" / "fixtures" / "cm" / "RF00005_svn.CM").read_text()
DESC = "ID   tRNA\nGA   29.00\nTC   29.00\nNC   28.90\nTP   Gene; tRNA;\n"
TRNA = "GCGGATTTAGCTCAGTTGGGAGAGCGCCAGACTGAAGATCTGGAGGTCCTGTGTTCGATCCACAGAATTCGCACCA"


def _headers(model):
    return [line for line in model.splitlines()
            if line.startswith(("INFERNAL1/", "HMMER3/", "NAME", "ACC", "//"))]


def test_keeps_cm_and_filter_hmm_renamed():
    model = extract_first_model(SVN_CM, "RF00005", name="tRNA")
    assert _headers(model) == [
        "INFERNAL1/a [1.1.4 | Dec 2020]", "NAME     tRNA", "ACC      RF00005", "//",
        "HMMER3/f [3.3.2 | Nov 2020]", "NAME  tRNA", "ACC   RF00005", "//"]
    assert "SEED" not in model
    assert model.endswith("//\n")


def test_name_defaults_to_accession_and_acc_is_not_duplicated():
    model = extract_first_model(SVN_CM, "RF00005")
    again = extract_first_model(model, "RF00005", name="tRNA")
    assert [h for h in _headers(again) if h.startswith(("NAME", "ACC"))] == [
        "NAME     tRNA", "ACC      RF00005", "NAME  tRNA", "ACC   RF00005"]
    assert "NAME     RF00005" in model


def test_parse_family_id():
    assert parse_family_id(DESC) == "tRNA"
    assert parse_family_id("AC   RF00005\n") is None
    assert parse_family_id(None) is None


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
    cm.write_text(extract_first_model(SVN_CM, "RF00005", parse_thresholds(DESC), "tRNA"))
    subprocess.run(["cmpress", str(cm)], check=True, capture_output=True)
    query = tmp_path / "q.fa"
    query.write_text(f">q\n{TRNA}\n")
    tbl = tmp_path / "q.tbl"
    subprocess.run(["cmscan", "--cut_ga", "--tblout", str(tbl), str(cm), str(query)],
                   check=True, capture_output=True)
    hits = [line.split()[:2] for line in tbl.read_text().splitlines()
            if not line.startswith("#")]
    assert hits == [["tRNA", "RF00005"]]


# The validator must show known-family hits, and say when the check could
# not run, even without -v: a silent failure looks exactly like "no known
# families matched". Output goes to stdout (captured in the report), or to
# stderr with --output-mode stdout.
import validate_stockholm  # noqa: E402

ENTRIES = [("q/1-76", TRNA)]


def test_cmscan_missing_is_reported(monkeypatch, capsys):
    monkeypatch.setattr(validate_stockholm.shutil, "which", lambda name: None)
    assert validate_stockholm.filter_known_families(ENTRIES, "db.cm") == []
    assert "known Rfam families were NOT checked" in capsys.readouterr().out


@pytest.mark.skipif(not shutil.which("cmscan"), reason="Infernal not installed")
def test_cmscan_failure_is_reported(tmp_path, capsys):
    # A database without GA thresholds: cmscan --cut_ga refuses it
    cm = tmp_path / "no_ga.cm"
    cm.write_text(extract_first_model(SVN_CM, "RF00005"))
    subprocess.run(["cmpress", str(cm)], check=True, capture_output=True)
    assert validate_stockholm.filter_known_families(ENTRIES, str(cm)) == []
    out = capsys.readouterr().out
    assert "cmscan failed, known Rfam families were NOT checked" in out
    assert "GA bit threshold unavailable" in out


def test_cmscan_exception_is_reported(monkeypatch, capsys):
    monkeypatch.setattr(validate_stockholm.shutil, "which", lambda name: "/bin/cmscan")

    def boom(*args, **kwargs):
        raise OSError("cannot run")
    monkeypatch.setattr(validate_stockholm.subprocess, "Popen", boom)
    assert validate_stockholm.filter_known_families(ENTRIES, "db.cm") == []
    assert "cmscan filtering failed, known Rfam families were NOT checked" in capsys.readouterr().out


def _fake_svn(monkeypatch, cm_files, desc_files):
    from scripts import build_live_cm as b
    monkeypatch.setattr(b, "fetch_family_list", lambda: sorted(cm_files))
    monkeypatch.setattr(b, "download_cm", cm_files.get)
    monkeypatch.setattr(b, "download_desc", desc_files.get)
    return b


def test_build_fails_if_a_download_fails(monkeypatch, tmp_path, capsys):
    b = _fake_svn(monkeypatch, {"RF00005": SVN_CM, "RF00006": None}, {"RF00005": DESC})
    with pytest.raises(SystemExit) as e:
        b.build_live_cm(str(tmp_path / "live.cm"), delay=0, press=False)
    assert e.value.code == 1
    assert "could not be downloaded after retries (RF00006)" in capsys.readouterr().err


def test_build_names_from_desc_and_keeps_names_unique(monkeypatch, tmp_path, capsys):
    b = _fake_svn(monkeypatch, {"RF00005": SVN_CM, "RF00006": SVN_CM},
                  {"RF00005": DESC, "RF00006": DESC})   # same ID in both
    out = tmp_path / "live.cm"
    b.build_live_cm(str(out), delay=0, press=False)
    names = [line.split()[1] for line in out.read_text().splitlines()
             if line.startswith(("NAME", "ACC"))]
    assert names == ["tRNA", "RF00005", "tRNA", "RF00005",
                     "RF00006", "RF00006", "RF00006", "RF00006"]
    assert "ID tRNA is already used by another family" in capsys.readouterr().err


@pytest.mark.skipif(not (shutil.which("cmpress") and shutil.which("cmscan")),
                    reason="Infernal not installed")
@pytest.mark.parametrize("to_stderr", [False, True])
def test_hits_are_shown_without_verbose(tmp_path, capsys, to_stderr):
    cm = tmp_path / "live.cm"
    cm.write_text(extract_first_model(SVN_CM, "RF00005", parse_thresholds(DESC), "tRNA"))
    subprocess.run(["cmpress", str(cm)], check=True, capture_output=True)
    details = validate_stockholm.filter_known_families(
        ENTRIES, str(cm), out=sys.stderr if to_stderr else None)
    assert len(details) == 1 and "matches RF00005 (tRNA)" in details[0]
    captured = capsys.readouterr()
    shown, other = (captured.err, captured.out) if to_stderr else (captured.out, captured.err)
    assert "q/1-76: matches RF00005 (tRNA)" in shown
    assert "1 sequence(s) match known Rfam families" in shown
    assert other == ""
