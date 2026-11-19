#!/usr/bin/env python3
"""
Build a live Rfam CM database from the Rfam SVN repository.

Downloads the latest CM files for all Rfam families from the SVN repo, adds
the GA/TC/NC thresholds from each family's DESC file, and concatenates them
into a single pressed CM database.

Usage:
    python3 scripts/build_live_cm.py
    python3 scripts/build_live_cm.py --output /path/to/Rfam_live.cm
    python3 scripts/build_live_cm.py --delay 0.2
"""

import argparse
import re
import shutil
import subprocess
import sys
import time


SVN_BASE_URL = "https://svn.rfam.org/svn/data_repos/trunk/Families"


def check_cmpress():
    """Check that cmpress is available in PATH."""
    if shutil.which("cmpress") is None:
        print("Error: cmpress not found in PATH. Install Infernal first.", file=sys.stderr)
        print("  See: http://eddylab.org/infernal/", file=sys.stderr)
        sys.exit(1)


def fetch_family_list():
    """Fetch the list of RF##### family IDs from the SVN directory listing."""
    result = subprocess.run(
        ["wget", "-q", "-O", "-", SVN_BASE_URL + "/"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("Error: Failed to fetch family listing from SVN.", file=sys.stderr)
        print(f"  URL: {SVN_BASE_URL}/", file=sys.stderr)
        if result.stderr:
            print(f"  wget: {result.stderr.strip()}", file=sys.stderr)
        sys.exit(1)

    families = sorted(set(re.findall(r'(RF\d{5})/', result.stdout)))
    if not families:
        print("Error: No RF##### families found in SVN listing.", file=sys.stderr)
        sys.exit(1)

    return families


def download_cm(family):
    """Download the CM file for a single family. Returns CM text or None."""
    return _download(family, "CM")


def download_desc(family):
    """Download the DESC file for a single family. Returns its text or None."""
    return _download(family, "DESC")


def _download(family, filename):
    url = f"{SVN_BASE_URL}/{family}/{filename}"
    result = subprocess.run(
        ["wget", "-q", "-O", "-", url],
        capture_output=True
    )
    # Some DESC files are not UTF-8 (e.g. Latin-1 author names): replace bad
    # bytes rather than crash. Model and threshold lines are plain ASCII.
    text = result.stdout.decode("utf-8", errors="replace")
    if result.returncode != 0 or not text.strip():
        return None
    return text


def _split_blocks(text):
    """Split a model file into blocks, each ending with its '//' line."""
    blocks, current = [], []
    for line in text.splitlines():
        current.append(line)
        if line.strip() == "//":
            blocks.append("\n".join(current) + "\n")
            current = []
    return blocks


def parse_thresholds(desc_text):
    """Read the GA, TC and NC bit-score thresholds from a family DESC file.

    Returns {'GA': 29.0, 'TC': ..., 'NC': ...} (TC and NC only if present),
    or None if there is no valid GA line.
    """
    thresholds = {}
    for line in (desc_text or "").splitlines():
        m = re.match(r'^(GA|TC|NC)\s+(-?\d+(?:\.\d+)?)\s*$', line)
        if m:
            thresholds.setdefault(m.group(1), float(m.group(2)))
    return thresholds if "GA" in thresholds else None


def extract_first_model(cm_text, family, thresholds=None):
    """Extract the first complete model and rename it to the family ID.

    An Infernal 1.1 model is two blocks, each ending with '//': the CM
    ('INFERNAL1/a ...') followed by its filter HMM ('HMMER3/f ...').
    cmpress and cmscan need both, so both are kept. NAME is 'SEED' in SVN
    files; it is set to the family ID (e.g. RF00001) in both blocks, so
    names are unique when all families are pressed into one database.

    SVN CM files have no GA/TC/NC lines (Rfam keeps them in DESC), and
    cmscan --cut_ga fails on a database where any model lacks GA. With
    thresholds (from parse_thresholds), they are written into the CM block
    after its NULL line, where the released Rfam.cm has them.

    Returns None if the file does not start with a CM followed by its HMM.
    """
    blocks = _split_blocks(cm_text)
    if (len(blocks) < 2 or not blocks[0].startswith("INFERNAL1/")
            or not blocks[1].startswith("HMMER3/")):
        return None

    renamed = []
    for block in blocks[:2]:
        block, n = re.subn(r'^(NAME\s+)\S+', rf'\g<1>{family}', block,
                           count=1, flags=re.MULTILINE)
        if n != 1:
            return None
        renamed.append(block)

    if thresholds:
        lines = [line for line in renamed[0].splitlines()
                 if not re.match(r'^(GA|TC|NC)\s', line)]
        null = next((i for i, line in enumerate(lines) if line.startswith("NULL ")), None)
        if null is None:
            return None
        cutoffs = [f"{tag:<9}{thresholds[tag]:.2f}" for tag in ("GA", "TC", "NC")
                   if tag in thresholds]
        lines[null + 1:null + 1] = cutoffs
        renamed[0] = "\n".join(lines) + "\n"

    return "".join(renamed)


def build_live_cm(output_path, delay=0.1, press=True):
    """Build a live Rfam CM database from SVN."""
    if press:
        check_cmpress()

    print("Fetching family list from Rfam SVN...")
    families = fetch_family_list()
    total = len(families)
    print(f"Found {total} families.")

    downloaded = 0
    skipped = 0

    with open(output_path, "w") as out:
        for i, family in enumerate(families, 1):
            cm_text = download_cm(family)
            thresholds = parse_thresholds(download_desc(family)) if cm_text else None
            model = extract_first_model(cm_text, family, thresholds) if thresholds else None
            if model:
                out.write(model)
                downloaded += 1
            else:
                if not cm_text:
                    reason = "could not download CM"
                elif not thresholds:
                    reason = "no GA threshold in DESC (needed for cmscan --cut_ga)"
                else:
                    reason = "CM file is not a CM followed by its filter HMM"
                print(f"Warning: {family}: {reason}; skipped", file=sys.stderr, flush=True)
                skipped += 1

            if i % 100 == 0:
                print(f"Downloaded {i}/{total} families...", flush=True)

            if delay > 0:
                time.sleep(delay)

    print(f"Done. Downloaded {downloaded} CMs, skipped {skipped}.")
    print(f"Output: {output_path}")

    if press:
        print("Running cmpress...")
        result = subprocess.run(["cmpress", output_path], capture_output=True, text=True)
        if result.returncode != 0:
            print("Error: cmpress failed.", file=sys.stderr)
            if result.stderr:
                print(result.stderr, file=sys.stderr)
            sys.exit(1)
        print("cmpress completed successfully.")

    print(f"CM database ready: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build a live Rfam CM database from the SVN repository."
    )
    parser.add_argument(
        "--output", "-o",
        default="Rfam_live.cm",
        help="Output path for the CM database (default: Rfam_live.cm)"
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.1,
        help="Delay between requests in seconds (default: 0.1)"
    )
    parser.add_argument(
        "--no-press",
        action="store_true",
        help="Skip running cmpress (useful when only the raw CM file is needed)"
    )
    args = parser.parse_args()

    build_live_cm(args.output, delay=args.delay, press=not args.no_press)


if __name__ == "__main__":
    main()
