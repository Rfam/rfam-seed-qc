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
    monkeypatch.setattr(b, "check_wget", lambda: None)   # downloads are faked
    return b


def test_build_fails_if_a_download_fails_and_keeps_old_database(monkeypatch, tmp_path, capsys):
    b = _fake_svn(monkeypatch, {"RF00005": SVN_CM, "RF00006": None}, {"RF00005": DESC})
    out = tmp_path / "live.cm"
    out.write_text("old database\n")
    (tmp_path / "live.cm.i1i").write_text("old index\n")
    with pytest.raises(SystemExit) as e:
        b.build_live_cm(str(out), delay=0, press=False)
    assert e.value.code == 1
    assert "could not be downloaded after retries (RF00006)" in capsys.readouterr().err
    assert out.read_text() == "old database\n"
    assert (tmp_path / "live.cm.i1i").read_text() == "old index\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["live.cm", "live.cm.i1i"]


def test_interrupted_build_keeps_old_database(monkeypatch, tmp_path):
    b = _fake_svn(monkeypatch, {"RF00005": SVN_CM}, {"RF00005": DESC})

    def interrupt(family):
        raise KeyboardInterrupt
    monkeypatch.setattr(b, "download_desc", interrupt)
    out = tmp_path / "live.cm"
    out.write_text("old database\n")
    with pytest.raises(KeyboardInterrupt):
        b.build_live_cm(str(out), delay=0, press=False)
    assert out.read_text() == "old database\n"
    assert [p.name for p in tmp_path.iterdir()] == ["live.cm"]


@pytest.mark.skipif(not shutil.which("cmpress"), reason="Infernal not installed")
def test_rebuild_into_same_path_replaces_pressed_database(monkeypatch, tmp_path):
    b = _fake_svn(monkeypatch, {"RF00005": SVN_CM}, {"RF00005": DESC})
    out = tmp_path / "live.cm"
    b.build_live_cm(str(out), delay=0, press=True)
    first = {p.name for p in tmp_path.iterdir()}
    assert first == {"live.cm", "live.cm.i1f", "live.cm.i1i", "live.cm.i1m", "live.cm.i1p"}
    b.build_live_cm(str(out), delay=0, press=True)   # cmpress refuses old indexes
    assert {p.name for p in tmp_path.iterdir()} == first


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


# --cm-db without --fix: the scan is offline, so it runs on plain validation
# too (it used to be accepted and silently ignored).

def _run_cli(*args, env=None):
    result = subprocess.run([sys.executable, str(REPO / "validate_stockholm.py"), *map(str, args)],
                            capture_output=True, text=True, env=env, timeout=120)
    return result.returncode, result.stdout + result.stderr


def _trna_alignment(tmp_path):
    sto = tmp_path / "trna.sto"
    sto.write_text(f"# STOCKHOLM 1.0\n\nX01.1/1-{len(TRNA)}   {TRNA}\n"
                   f"X02.1/1-{len(TRNA)}   {TRNA}\n//\n")
    return sto


@pytest.mark.skipif(not (shutil.which("cmpress") and shutil.which("cmscan")),
                    reason="Infernal not installed")
def test_cm_db_without_fix_reports_hits(tmp_path):
    cm = tmp_path / "live.cm"
    cm.write_text(extract_first_model(SVN_CM, "RF00005", parse_thresholds(DESC), "tRNA"))
    subprocess.run(["cmpress", str(cm)], check=True, capture_output=True)
    sto = _trna_alignment(tmp_path)
    code, output = _run_cli("--cm-db", cm, sto)
    assert code == 0, output
    assert "Valid Stockholm file" in output
    assert f"X01.1/1-{len(TRNA)}: matches RF00005 (tRNA)" in output
    assert "2 sequence(s) match known Rfam families" in output
    # validation only: nothing is written next to the input
    assert sorted(p.name for p in tmp_path.iterdir() if p.suffix != ".cm" and ".cm." not in p.name) == [
        "trna.sto"]


def test_cm_db_without_fix_warns_if_cmscan_missing(tmp_path):
    sto = _trna_alignment(tmp_path)
    code, output = _run_cli("--cm-db", tmp_path / "db.cm", sto, env={"PATH": ""})
    assert code == 0, output
    assert "cmscan not found; known Rfam families were NOT checked" in output


def test_cm_db_without_fix_skips_invalid_files(tmp_path):
    sto = tmp_path / "bad.sto"
    sto.write_text("# STOCKHOLM 1.0\n\nX01.1/1-4 ACGU\nX02.1/1-3 ACG\n//\n")
    code, output = _run_cli("--cm-db", tmp_path / "db.cm", sto, env={"PATH": ""})
    assert code == 1
    assert "INVALID" in output and "Rfam families" not in output


@pytest.mark.parametrize("version_output, ok", [
    ("GNU Wget 1.25.0 built on linux-gnu.", True),
    ("GNU Wget 1.19.1 built on linux-gnu.", True),
    ("GNU Wget 1.21 built on linux-gnu.", True),
    ("GNU Wget 1.19 built on linux-gnu.", False),
    ("GNU Wget 1.14 built on linux-gnu.", False),
    ("something else", False),
])
def test_check_wget_version(monkeypatch, version_output, ok):
    from types import SimpleNamespace
    from scripts import build_live_cm as b
    monkeypatch.setattr(b.shutil, "which", lambda name: "/usr/bin/wget")
    monkeypatch.setattr(b.subprocess, "run",
                        lambda *args, **kwargs: SimpleNamespace(stdout=version_output, returncode=0))
    if ok:
        b.check_wget()
    else:
        with pytest.raises(SystemExit) as e:
            b.check_wget()
        assert "wget 1.19.1 or later is needed" in str(e.value)


def test_check_wget_missing(monkeypatch):
    from scripts import build_live_cm as b
    monkeypatch.setattr(b.shutil, "which", lambda name: None)
    with pytest.raises(SystemExit) as e:
        b.check_wget()
    assert "wget not found in PATH" in str(e.value)
