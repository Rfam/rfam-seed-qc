"""
Tests for validate_stockholm.py, run on the small files in tests/fixtures/.

- valid/, invalid/, fixable/: validation only (no --fix). Offline and fast.
- fix/: runs --fix, which looks sequences up at NCBI (marked "network").
  Kept tiny on purpose: two NCBI lookups and no BLAST search.

Run all:          uv run pytest
Offline only:     uv run pytest -m "not network"
"""

import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / 'validate_stockholm.py'
FIXTURES = Path(__file__).resolve().parent / 'fixtures'


def run_validator(*args, cwd=None, timeout=60, cache=None, env_extra=None):
    """Run the validator and return (exit code, combined output).

    cache: folder for RFAM_SEED_QC_CACHE (default: the user's normal cache).
    env_extra: extra environment variables.
    """
    env = dict(os.environ)
    if cache is not None:
        env['RFAM_SEED_QC_CACHE'] = str(cache)
    env.update(env_extra or {})
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        capture_output=True, text=True, cwd=cwd, timeout=timeout, env=env,
    )
    return result.returncode, result.stdout + result.stderr


# Each fixture: (path in tests/fixtures, expected exit code, expected status
# line, messages that must appear in the output)
VALIDATION_CASES = [
    ('valid/valid.sto', 0, 'Valid Stockholm file', []),
    ('valid/missing_ss_cons.sto', 0, 'Valid Stockholm file',
     ['Missing 2D structure consensus annotation (#=GC SS_cons)']),
    ('invalid/empty.sto', 1, 'INVALID', ['File is empty']),
    ('invalid/missing_header.sto', 1, 'INVALID',
     ["Missing required '# STOCKHOLM 1.0' header"]),
    ('invalid/bad_version.sto', 1, 'INVALID',
     ["Non-standard Stockholm version '2.0'"]),
    ('invalid/missing_terminator.sto', 1, 'INVALID',
     ["Missing required '//' terminator"]),
    ('invalid/inconsistent_lengths.sto', 1, 'INVALID',
     ['Sequences have inconsistent lengths']),
    ('invalid/no_sequences.sto', 1, 'INVALID',
     ['No sequences found in alignment']),
    ('fixable/RF00569_with_dup.so', 0, 'Valid but has fixable issues',
     ['Found 1 duplicate sequence(s)']),
    ('fixable/RF00569_with_overlap.so', 0, 'Valid but has fixable issues',
     ['Found 1 overlapping sequence(s)']),
    ('fixable/missing_coordinates.sto', 0, 'Valid but has fixable issues',
     ['Found 1 sequence(s) missing coordinates']),
]


@pytest.mark.parametrize('fixture, exit_code, status, messages', VALIDATION_CASES,
                         ids=[case[0] for case in VALIDATION_CASES])
def test_validation(fixture, exit_code, status, messages):
    code, output = run_validator('-v', FIXTURES / fixture)
    assert code == exit_code, output
    assert status in output, output
    for message in messages:
        assert message in output, output
    if status == 'Valid Stockholm file':
        assert 'fixable issues' not in output, output


def test_every_fixture_is_tested():
    """A fixture added to valid/, invalid/ or fixable/ must be listed above."""
    listed = {case[0] for case in VALIDATION_CASES}
    on_disk = {str(p.relative_to(FIXTURES)) for folder in ('valid', 'invalid', 'fixable')
               for p in (FIXTURES / folder).iterdir() if p.is_file()}
    assert on_disk == listed


def test_version():
    """--version prints the version from the VERSION file and exits, with no file needed."""
    code, output = run_validator('--version')
    assert code == 0, output
    assert output.strip() == f"rfam-seed-qc {(REPO / 'VERSION').read_text().strip()}"


def test_several_files_fail_if_one_is_invalid():
    code, output = run_validator(FIXTURES / 'valid/valid.sto',
                                 FIXTURES / 'invalid/missing_terminator.sto')
    assert code == 1, output
    assert 'Validation failed for one or more files' in output, output


@pytest.mark.network
def test_fix_removes_duplicate_and_overlap(tmp_path):
    """--fix on two real sequences plus a duplicate and an overlapping copy.

    The overlapping copy is shifted 5 nt, so it does not match NCBI. It starts
    earlier than the real one, so it would be kept if overlaps were resolved
    by position alone: the real one, which matches NCBI, must be kept instead.
    Only the two real sequences are validated, so no BLAST search runs.
    """
    src = FIXTURES / 'fix/duplicate_and_overlap.sto'
    shutil.copy(src, tmp_path / src.name)   # outputs go next to the input

    cache = tmp_path / 'cache'
    code, output = run_validator('-v', '--fix', src.name, cwd=tmp_path, timeout=300, cache=cache)
    assert code == 0, output
    assert list(cache.glob('*.fasta')), 'downloads were not cached in RFAM_SEED_QC_CACHE'

    assert 'Removed 1 duplicate sequence(s)' in output, output
    assert 'Found 1 overlapping sequence(s) to remove' in output, output
    assert 'AC104446.2/142001-142075' not in output.split('Found 1 overlapping')[1], output
    assert 'AC104446.2/142006-142080: validated OK' in output, output
    assert 'AAYZ01025109.1/988-1063: validated OK' in output, output
    assert 'Running BLAST search' not in output, output

    corrected = tmp_path / 'duplicate_and_overlap_corrected.sto'
    report = tmp_path / 'duplicate_and_overlap_Report.txt'
    assert corrected.exists() and report.exists(), output

    text = corrected.read_text()
    names = [line.split()[0] for line in text.splitlines()
             if line.strip() and not line.startswith('#') and line.strip() != '//']
    assert names == ['AC104446.2/142006-142080', 'AAYZ01025109.1/988-1063'], text
    assert 'Processed by rfam-seed-qc v' in text, text
    assert text.startswith('# STOCKHOLM 1.0') and text.rstrip().endswith('//'), text

    # The corrected file must itself validate cleanly
    code, output = run_validator('-v', corrected)
    assert code == 0, output
    assert 'Valid Stockholm file' in output and 'fixable issues' not in output, output


@pytest.mark.skipif(hasattr(os, 'geteuid') and os.geteuid() == 0,
                    reason='root can write to read-only folders')
def test_fix_stops_if_cache_is_not_writable(tmp_path):
    """A cache that cannot be written must stop the run, not make every
    lookup fail and remove valid sequences as "not found". Offline: the run
    stops before the first download."""
    src = FIXTURES / 'fix/duplicate_and_overlap.sto'
    shutil.copy(src, tmp_path / src.name)
    cache = tmp_path / 'cache'
    cache.mkdir()
    cache.chmod(stat.S_IRUSR | stat.S_IXUSR)   # read-only
    try:
        code, output = run_validator('--fix', src.name, cwd=tmp_path, cache=cache)
    finally:
        cache.chmod(stat.S_IRWXU)
    assert code != 0, output
    assert 'Set RFAM_SEED_QC_CACHE' in output, output
    assert 'Stopping so no sequence is removed by mistake' in output, output
    assert not (tmp_path / 'duplicate_and_overlap_corrected.sto').exists(), output


def test_fix_stops_if_ncbi_is_unreachable(tmp_path):
    """No connection to NCBI (offline node, DNS, firewall) must stop the run,
    not turn every lookup into "not found" and remove every sequence. Offline:
    all connections go to a proxy port where nothing listens."""
    src = FIXTURES / 'fix/duplicate_and_overlap.sto'
    shutil.copy(src, tmp_path / src.name)
    unreachable = 'http://127.0.0.1:9'
    code, output = run_validator(
        '-v', '--fix', src.name, cwd=tmp_path, cache=tmp_path / 'cache',
        env_extra={'https_proxy': unreachable, 'HTTPS_PROXY': unreachable,
                   'no_proxy': '', 'NO_PROXY': ''})
    assert code != 0, output
    assert 'Could not reach NCBI' in output, output
    assert 'Stopping so no sequence is removed by mistake' in output, output
    assert 'NOT FOUND' not in output and 'Running BLAST search' not in output, output
    assert not (tmp_path / 'duplicate_and_overlap_corrected.sto').exists(), output
