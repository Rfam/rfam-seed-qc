# rfam-seed-qc
Tools for ensuring the quality of multiple sequence alignments (MSAs) in Rfam.

## Stockholm File Validation

This repository includes a modular validation script for Stockholm format alignment files (`.so`, `.sto`, `.stk`).

### Usage

```bash
python3 validate_stockholm.py [-v] [--fix] [--output-mode {stdout,file}] [--cm-db Rfam.cm] <file1.so> [file2.so ...]
```

Options:
- `-v, --verbose`: Print detailed validation information
- `--fix`: Attempt to fix fixable errors automatically
- `--output-mode {stdout,file}`: Output mode for fixed files (default: file)
  - `file`: Create a new file with `_corrected` suffix
  - `stdout`: Print corrected content to stdout
- `--cm-db <path>`: Path to Rfam CM database for filtering sequences matching known Rfam families (optional)

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
2. **BLAST Fallback**: If direct lookup fails, perform a BLAST search against NCBI's nucleotide database (accepts hits with >=95% identity, >=90% coverage, e-value <=1e-10)
3. **Removal**: Sequences that cannot be resolved via either method are removed from the output

**Requirements**: BioPython (`pip install biopython`)

### Overlap Detection and Removal

Sequences from the same accession (species) that overlap by at least 1 bp are detected. When using `--fix`, overlapping sequences are removed using greedy interval scheduling to keep the maximum number of non-overlapping sequences.

### Sequence Validation

When using `--fix`, all sequences are validated against NCBI to ensure they match the source data at the given coordinates:

1. **NCBI Validation**: Fetches the source sequence from NCBI and compares it against the Stockholm sequence at the given coordinates
2. **BLAST Fallback**: If the accession is not found in NCBI or the sequence doesn't match, a BLAST search is performed to find the correct accession and coordinates
3. **Removal**: Sequences that fail both NCBI validation and BLAST are removed from the output

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

Requires `esl-reformat` from the [Easel](http://eddylab.org/software/easel/) library, which is distributed as part of [Infernal](http://eddylab.org/infernal/). If `esl-reformat` is not found in PATH, this step is skipped and a warning is printed.

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
| `BLAST_HITLIST_SIZE` | 5 | Number of BLAST hits to retrieve |
| `CMSCAN_MAX_EVALUE` | 1e-3 | Maximum e-value for cmscan hits (known family warnings) |
| `NCBI_REQUEST_DELAY` | 0.5 | Delay between NCBI requests in seconds |

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
```

### Continuous Integration

The repository includes a GitHub Action that automatically validates Stockholm files in pull requests. The CI check will:
- Trigger when a PR modifies `.so`, `.sto`, or `.stk` files
- Run the validation script on all changed files
- Pass only if all files are valid

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

## Contact us

If you have any questions or feedback, feel free to submit a GitHub issue or [email us](https://docs.rfam.org/en/latest/).
