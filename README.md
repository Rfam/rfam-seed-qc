# rfam-seed-qc

[![Release](https://img.shields.io/github/v/release/Rfam/rfam-seed-qc)](https://github.com/Rfam/rfam-seed-qc/releases)
[![License](https://img.shields.io/github/license/Rfam/rfam-seed-qc)](LICENSE)
![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue)

Tools for ensuring the quality of multiple sequence alignments (MSAs) in Rfam.

## Stockholm File Validation

This repository includes a modular validation script for Stockholm format alignment files (`.so`, `.sto`, `.stk`).

### Installation

**Option A: conda (everything included).** Installs Python, BioPython and the external tools (Infernal, Easel, BLAST):

```bash
git clone https://github.com/Rfam/rfam-seed-qc && cd rfam-seed-qc
conda env create -f environment.yml
conda activate rfam-seed-qc
python3 validate_stockholm.py -v file.sto
```

**Option B: [uv](https://docs.astral.sh/uv/) (Python only).** On first use, `uv run` creates a local environment with the pinned dependencies from `uv.lock`:

```bash
git clone https://github.com/Rfam/rfam-seed-qc && cd rfam-seed-qc
uv run validate_stockholm.py -v file.sto
```

uv does not install the non-Python tools. Without [Infernal](http://eddylab.org/infernal/) (`cmscan`, `cmpress`, `esl-reformat`), validation and `--fix` still work, but `--cm-db` filtering is unavailable and the all-gap column removal is skipped with a warning. Install it with e.g. `brew install infernal`, or `conda install -c conda-forge -c bioconda infernal easel` (bioconda's `infernal` package does not include `esl-reformat`; it is in `easel`). `--local` also needs [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (`blastn`, `blastdbcmd`, `update_blastdb.pl`), e.g. `conda install -c bioconda blast`.

### Usage

```bash
python3 validate_stockholm.py [-v] [--fix] [--output-mode {stdout,file}] [--cm-db Rfam.cm] [--local [--blast-db BLAST_LOCAL_DB] [--threads N]] <file1.so> [file2.so ...]
```

Options:
- `-v, --verbose`: Print detailed validation information
- `--fix`: Attempt to fix fixable errors automatically
- `--output-mode {stdout,file}`: Output mode for fixed files (default: file)
  - `file`: Create a new file with `_corrected` suffix, plus a `_Report.txt` log of the run (starting with the tool version)
  - `stdout`: Print corrected content to stdout
  - The corrected alignment records the tool and version that produced it in a `#=GF CC   Corrected by rfam-seed-qc vX.Y.Z (...)` line
- `--cm-db <path>`: Path to Rfam CM database for filtering sequences matching known Rfam families (optional)
- `--local`: Run BLAST and sequence lookups locally with BLAST+ (`blastn`, `blastdbcmd`) instead of NCBI's web services. Needs BLAST+ (included in the conda environment) and a local database; see [Local BLAST](#local-blast)
- `--blast-db <BLAST_LOCAL_DB>`: Local BLAST database to use with `--local`: a name found in `$BLASTDB` or a path (default: `nt`). A custom database must be built with `makeblastdb -parse_seqids`, so hits can be looked up by name
- `--threads <N>`: Threads for local `blastn` with `--local` (default: the CPUs available to the job, capped by `SLURM_CPUS_PER_TASK`). Set it explicitly in containers limited by CPU quota (e.g. `docker --cpus`)

#### NCBI rate limits

NCBI allows 3 requests per second per IP address. On clusters, many users share one IP, so NCBI may refuse requests (HTTP 429) even for small runs. If NCBI keeps refusing after the retries, or rejects the API key, the run stops with an error rather than removing sequences it could not check. To avoid this:

- Use `--local` (see below), which reads sequences from a local database; or
- Set an NCBI API key (free, from your [NCBI account settings](https://www.ncbi.nlm.nih.gov/account/settings/)). The limit then applies to your key (10 requests per second) instead of the shared IP. The key is read from the environment and never written to reports:

```bash
export NCBI_API_KEY=your_key_here
```

#### Local BLAST

NCBI's web BLAST searches `core_nt` even when asked for `nt`, and each search can wait minutes in NCBI's queue. With `--local`, searches run on your machine against the full `nt`, using all available CPUs. Both modes use the same blastn settings, so `--local` is not more sensitive, but it is more complete: a sequence that is only in `nt` can be rescued with `--local` but not by the web search. Download the database once. In November 2026 `nt` was 1.1 TB to download and 1.3 TB on disk; check the current size first with `update_blastdb.pl --source ncbi --showall pretty | grep -w nt`:

```bash
export BLASTDB=/path/to/blastdb && mkdir -p $BLASTDB && cd $BLASTDB
update_blastdb.pl --source ncbi --decompress --num_threads 4 nt
```

`update_blastdb.pl` downloads one volume at a time, so this can take several days: run it as a batch job, not on a login node. If it stops, run it again and it resumes.

Then run with `--local` (it finds `nt` through `$BLASTDB`), or pass the path with `--blast-db /path/to/blastdb/nt`. There are no incremental updates: re-run `update_blastdb.pl` to refresh it. With `-v`, the database date is printed in the report. In local mode, sequences are also read from the local database (with `blastdbcmd`), so there are no NCBI rate limits and jobs can run without internet access. Only accessions missing from the local database (`nt` omits some divisions, such as most WGS contigs) are fetched from NCBI, and RNAcentral IDs are still checked against RNAcentral.

### What is validated?

**Fatal Errors** (must be fixed manually):
- Missing `# STOCKHOLM 1.0` header
- Missing `//` terminator
- No sequences found in alignment
- **All sequences must have the same length**
- **Sequences must not contain whitespace characters**

**Fixable Errors** (can be auto-corrected with `--fix`):
- Duplicate sequences (same accession, coordinates, and sequence data)
- Missing coordinates (sequences without start/end positions)
- Overlapping sequences (sequences from the same accession that overlap by >=1 bp)

**Warnings** (non-critical):
- Sequences matching known Rfam families (when `--cm-db` is provided) — flagged but not removed
- Missing 2D structure consensus annotation (`#=GC SS_cons`)
- `#=GC SS_cons` length mismatch with sequence length
- `#=GC SS_cons` invalid WUSS notation characters or unbalanced brackets
- Lines exceeding 10,000 character limit

### Modular Architecture

The validation logic is split into separate modules in the `scripts/` directory:
- `fatal_errors.py`: Errors that cannot be automatically fixed
- `fixable_errors.py`: Errors that can be automatically corrected (including coordinate fixing, overlap removal, and BLAST fallback)
- `stockholm_warnings.py`: Non-critical issues
- `parser.py`: Stockholm file parsing utilities (supports both single-block and interleaved multi-block format)
- `alignment_stats.py`: Pairwise identity computation for alignment quality assessment
- `config.py`: Configurable parameters (BLAST thresholds, cmscan E-value, NCBI request delay)

### Sequence Format

Sequences should follow the format: `ACCESSION/START-END` where:
- `ACCESSION` is the sequence identifier (e.g., from GenBank like `AF228364.1`)
- `START-END` are the coordinates indicating which portion of the original sequence is included (e.g., `1-74`)

Example: `AF228364.1/1-74`

Sequence data:
- May contain any characters except whitespace
- Gaps may be indicated by `.` or `-`

### Duplicate Detection and Removal

The script can detect and remove duplicate sequences using the `--fix` flag. Duplicates are defined as sequences that have:
1. The same accession/identifier
2. The same coordinates
3. The exact same sequence data

### Coordinate Fixing

When sequences are missing coordinates (e.g., `NZ_CP038662.1` instead of `NZ_CP038662.1/4723652-4723704`), the `--fix` flag will:

1. **NCBI Direct Lookup**: Download the FASTA sequence from NCBI and find coordinates by matching the alignment sequence
2. **BLAST Fallback**: If direct lookup fails, perform a BLAST search against NCBI's nucleotide database (accepts hits with >=95% identity, >=90% coverage, e-value <=1e-10). Searches use standard **blastn**, not megablast, with its parameters set explicitly (word size 11, match/mismatch 2/−3, gap costs 5/2) so results don't change if NCBI changes its defaults. blastn is slower than megablast but more sensitive, including for short queries. Among the hits that pass, an exact full-length match is preferred, since only exact matches pass the final check
3. **Removal**: Sequences that cannot be resolved via either method are removed from the output

With `--local`, the lookup and the BLAST search use the local database instead (see [Local BLAST](#local-blast)).

**Requirements**: BioPython (`pip install biopython`); BLAST+ for `--local`

### Overlap Detection and Removal

Sequences from the same accession (species) that overlap by at least 1 bp are detected. When using `--fix`, overlapping sequences are removed using greedy interval scheduling to keep the maximum number of non-overlapping sequences.

### Sequence Validation

When using `--fix`, all sequences are validated against NCBI to ensure they match the source data at the given coordinates:

1. **NCBI Validation**: Fetches the source sequence from NCBI and compares it against the Stockholm sequence at the given coordinates
2. **BLAST Fallback**: If the accession is not found in NCBI or the sequence doesn't match, a BLAST search is performed. If a hit is found (>=95% identity, >=90% coverage, e-value <=1e-10), the hit's accession is re-fetched from NCBI and the sequence is verified for an exact match before the rescue is accepted. This ensures compatibility with Rfam's `rfsearch.pl` MD5 validation.
3. **Removal**: Sequences that fail both NCBI validation and BLAST (or whose BLAST hit does not exactly match NCBI) are removed from the output

With `--local`, sequences are fetched from the local database, and from NCBI only when an accession is missing from it.

### RNAcentral Identifiers

Sequences named with RNAcentral identifiers (e.g. `URS0000572555_1268238/1-81`, or `URS0000572555/1-81` without a taxid) are validated against the [RNAcentral REST API](https://rnacentral.org/help/rest-api) instead of NCBI:

1. **Sequence**: the RNAcentral entry `/api/v1/rna/<URS>` must exist and be active. Coordinates are 1-based positions within the RNAcentral sequence (forward strand only), and the Stockholm sequence must exactly match that slice (U and T are treated as equal).
2. **Taxid**: if the name has a taxid, `/api/v1/rna/<URS>_<taxid>` must exist, i.e. the URS is annotated in that species. For URS found in many species this endpoint can be very slow or return HTTP 500. In that case the sequence (already validated) is kept, and the taxid is reported as unverified.
3. If validation fails (URS not found, inactive, sequence mismatch, bad coordinates, or taxid not annotated), the sequence follows the normal **BLAST fallback** path. A BLAST rescue renames it to the matching NCBI accession/coordinates; otherwise it is removed.
4. URS identifiers without coordinates are located within their RNAcentral sequence when using `--fix`.

Requests follow RNAcentral's access guidelines: one request every 2 seconds (the `/api/v1/rna/` limit is 30/minute), `Retry-After` is honoured on HTTP 429, and a descriptive User-Agent is sent. Expect about 2–4 seconds per URS sequence; taxid checks for URS found in many species can take up to about a minute. See the `RNACENTRAL_*` settings in `scripts/config.py`.

### Filtering Known Rfam Families (cmscan)

When building a new Rfam family, you can check whether sequences already belong to existing Rfam families using the `--cm-db` option. This requires:

1. **Download Rfam.cm** from https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
2. **Decompress**: `gunzip Rfam.cm.gz`
3. **Press the database**: `cmpress Rfam.cm` (generates `.i1f`, `.i1i`, `.i1m`, `.i1p` index files)
4. **Install Infernal**: `cmscan` must be available in your PATH (http://eddylab.org/infernal/)

cmscan is run with `--cut_ga` to use each model's GA (gathering) threshold, matching how Rfam curates family membership. Hits that pass GA are flagged as warnings but **not removed**. An additional E-value safety filter (`CMSCAN_MAX_EVALUE` in `scripts/config.py`, default 1e-3) is applied on top.

```bash
# Check for known families during fixing
python3 validate_stockholm.py -v --fix --cm-db Rfam.cm input.sto
```

### Building a Live CM Database from SVN

The released `Rfam.cm` is static and only updated with each Rfam release. To filter against the latest curated families (including those not yet released), you can build a live CM database directly from the Rfam SVN repository:

```bash
# Build live CM database (downloads all family CMs and runs cmpress)
python3 scripts/build_live_cm.py

# Custom output path
python3 scripts/build_live_cm.py --output /path/to/Rfam_live.cm

# Then use with the validator
python3 validate_stockholm.py -v --fix --cm-db Rfam_live.cm input.sto
```

This fetches CM files from `https://svn.rfam.org/svn/data_repos/trunk/Families/` (publicly readable, no credentials needed) and concatenates them into a single pressed CM database. Requires `wget` and `cmpress` (Infernal) in your PATH.

### Pairwise Identity

When running with `--fix -v`, the script computes and displays the average pairwise identity for each sequence in the final alignment. Sequences with identity significantly below the alignment average (>20 percentage points) are flagged as potential outliers.

### Gap Minimization

After all sequence filtering and corrections are applied, the script runs `esl-reformat --mingap stockholm` on the final alignment to remove any all-gap columns. All-gap columns can be introduced when sequences are removed during earlier steps (duplicate removal, overlap resolution, NCBI validation). This step ensures the corrected alignment has no unnecessary gap-only columns.

Requires `esl-reformat` from the [Easel](http://eddylab.org/software/easel/) library, which is distributed with [Infernal](http://eddylab.org/infernal/) source releases and Homebrew; with conda it is in the separate bioconda `easel` package (included in `environment.yml`). If `esl-reformat` is not found in PATH, this step is skipped and a warning is printed.

The report includes the number of all-gap columns removed and the resulting alignment width.

### Report Output

When using `--fix` in file output mode, a report file (`<stem>_Report.txt`) is generated alongside the corrected alignment, capturing the full verbose output of the validation and fixing process.

### Configuration

Parameters can be adjusted in `scripts/config.py`:

| Parameter | Default | Description |
|---|---|---|
| `BLAST_MIN_IDENTITY` | 95 | Minimum identity % for BLAST hits |
| `BLAST_MIN_COVERAGE` | 90 | Minimum coverage % for BLAST hits |
| `BLAST_MAX_EVALUE` | 1e-10 | Maximum e-value for BLAST hits |
| `BLAST_HITLIST_SIZE` | 250 | Number of BLAST hits to retrieve. Kept large because BLAST uses it while searching, and small values can miss the best hit |
| `BLAST_WORD_SIZE` | 11 | blastn word size |
| `BLAST_MATCH_REWARD` / `BLAST_MISMATCH_PENALTY` | 2 / -3 | blastn match/mismatch scores |
| `BLAST_GAP_COSTS` | `"5 2"` | blastn gap open/extend costs |
| `BLAST_NUM_THREADS` | CPUs available | Threads for local `blastn` (`--local` only; `--threads` overrides it) |
| `CMSCAN_MAX_EVALUE` | 1e-3 | Maximum e-value for cmscan hits (known family warnings) |
| `NCBI_REQUEST_DELAY` | 0.5 | Delay between NCBI requests in seconds |
| `NCBI_TIMEOUT` | 60 | Timeout in seconds for one NCBI FASTA download |
| `NCBI_TOOL` | `rfam-seed-qc` | Tool name sent to NCBI (BLAST and E-utilities) |
| `NCBI_EMAIL` | `None` | Your email, sent to NCBI so they can contact you before blocking heavy use (recommended) |
| `NCBI_API_KEY` | from the `NCBI_API_KEY` environment variable | Optional NCBI API key (see [NCBI rate limits](#ncbi-rate-limits)). Never put it in `config.py` |
| `NCBI_REQUEST_DELAY_WITH_KEY` | 0.1 | Delay between NCBI requests in seconds when an API key is set |
| `NCBI_BLAST_RETRY_WAIT` | 30 | Seconds before retrying a web BLAST request NCBI refused (×attempt) |
| `RNACENTRAL_REQUEST_DELAY` | 2.0 | Delay between RNAcentral requests in seconds |
| `RNACENTRAL_MAX_RETRIES` | 3 | Retries on HTTP 429/5xx/network errors |
| `RNACENTRAL_TIMEOUT` | 30 | RNAcentral request timeout in seconds |
| `RNACENTRAL_TAXID_TIMEOUT` | 60 | Timeout for the species-specific (taxid) check |
| `RNACENTRAL_TAXID_MAX_RETRIES` | 1 | Retries for the taxid check before marking it unverified |
| `RNACENTRAL_USER_AGENT` | `rfam-seed-qc (...)` | User-Agent sent to RNAcentral |

### Examples

```bash
# Validate a single file
python3 validate_stockholm.py example_valid.so

# Validate multiple files with verbose output
python3 validate_stockholm.py -v file1.so file2.so file3.so

# Fix errors and create corrected file
python3 validate_stockholm.py --fix file.so

# Fix and output to stdout
python3 validate_stockholm.py --fix --output-mode stdout file.so

# Fix with known Rfam family filtering
python3 validate_stockholm.py -v --fix --cm-db Rfam.cm file.sto

# Fix with local BLAST against nt, writing file_corrected.sto and file_Report.txt
python3 validate_stockholm.py -v --fix --output-mode file --local --blast-db /path/to/blastdb/nt --threads 16 file.sto
```

### Continuous Integration

**Curating alignments** (`validate-stockholm.yml`): push a Stockholm file (`.so`, `.sto`, `.stk`, `.stockholm`) to the `alignments/` folder on any branch except `main`. The workflow runs `validate_stockholm.py -v --fix --cm-db models/15.1/Rfam.cm` on each new or changed file (web BLAST, not `--local`) and commits the `_corrected` file and `_Report.txt` back to the branch, next to the input. Files outside `alignments/` are ignored.

**Tests** (`tests.yml`): on every push and pull request, except pushes that only change `alignments/`, the validator is run on the small files in `tests/fixtures/`:
- `valid/`, `invalid/`, `fixable/`: validation only, offline. Each file checks one case, such as a missing `//` or a duplicate sequence, for the expected exit code and message.
- `fix/`: one `--fix` run on two real sequences plus a duplicate and an overlap. It looks up two sequences at NCBI (no BLAST search), so it takes seconds.

Run the tests locally with `uv run pytest`, or `uv run pytest -m "not network"` to skip the NCBI test. When adding a fixture to `valid/`, `invalid/` or `fixable/`, add its expected result to `tests/test_validator.py`; a test fails if a fixture is not listed.

### Stockholm Format Reference

Stockholm format is used for multiple sequence alignments. Basic structure:
```
# STOCKHOLM 1.0
AF228364.1/1-74    CGGCAGAUGAUGAU-UUUACUUGGAUUCCCCUUCAGAACAUUUA
AF228365.1/1-73    CGGCAGAUGAUGAU-UUUACUUGGAUUCCCCUUCAGAACAUUU
#=GC SS_cons       <<<<_______..________.__._.______.___.___.___
//
```

The `#=GC SS_cons` line is the 2D structure consensus annotation, which represents the secondary structure of the RNA alignment. While not strictly required, it is recommended for Rfam alignments.

For more information, see the [Stockholm format specification](https://en.wikipedia.org/wiki/Stockholm_format).

## Citation

If you use rfam-seed-qc, please cite it using the metadata in [`CITATION.cff`](CITATION.cff) (GitHub's **Cite this repository** button provides APA and BibTeX formats).

When releasing a new version, update `VERSION`, `version` and `date-released` in `CITATION.cff`, and `version` in `pyproject.toml` (then run `uv lock`). A CI check fails if the versions differ or `uv.lock` is out of date.

## Contact us

If you have any questions or feedback, feel free to submit a GitHub issue or [email us](https://docs.rfam.org/contact-us.html).
