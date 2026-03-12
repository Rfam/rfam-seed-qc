#!/usr/bin/env python3
"""
Stockholm format validator for MSA files.

This script validates Stockholm format alignment files (.so, .sto, .stk)
by checking for required structural elements and basic formatting.
It can automatically fix certain errors like duplicate sequences.
"""

# Suppress BioPython warnings before any imports
import warnings as stdlib_warnings
stdlib_warnings.filterwarnings("ignore", module="Bio")

import sys
import io
import shutil
import subprocess
import tempfile
import argparse
from datetime import datetime
from pathlib import Path

# Import validation modules
from scripts import fatal_errors, fixable_errors, stockholm_warnings as warnings, parser, alignment_stats
from scripts.config import CMSCAN_MAX_EVALUE


def validate_stockholm_file(filepath):
    """
    Validate a Stockholm format alignment file.
    
    Args:
        filepath: Path to the Stockholm file to validate
        
    Returns:
        dict: Validation results with errors, warnings, and fixable issues
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        return {
            'is_valid': False,
            'fatal_errors': [f"Failed to read file: {e}"],
            'fixable_errors': [],
            'warnings': [],
            'can_be_fixed': False
        }
    
    # Parse the file (handles interleaved format by merging blocks)
    parsed = parser.parse_stockholm_file(lines)
    sequence_entries = parsed['sequence_entries']
    raw_entries = parsed['raw_sequence_entries']

    # Build unique sequences dict for validation (block-aware duplicate detection)
    unique_sequences, duplicates = fixable_errors.find_duplicates_from_entries(
        sequence_entries, raw_entries=raw_entries
    )
    
    # Collect errors and warnings
    fatal_error_list = []
    fixable_error_list = []
    warning_list = []
    
    # Check for fatal errors
    is_valid, error_msg = fatal_errors.check_empty_file(lines)
    if not is_valid:
        fatal_error_list.append(error_msg)
        return {
            'is_valid': False,
            'fatal_errors': fatal_error_list,
            'fixable_errors': [],
            'warnings': [],
            'can_be_fixed': False
        }
    
    is_valid, error_msg, line_num = fatal_errors.check_header(lines)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    is_valid, error_msg = fatal_errors.check_terminator(lines)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    is_valid, error_msg = fatal_errors.check_no_sequences(unique_sequences)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    # Check sequence length consistency
    is_valid, error_msg = fatal_errors.check_sequence_lengths(unique_sequences)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    # Check sequence characters
    is_valid, error_msg = fatal_errors.check_sequence_characters(unique_sequences)
    if not is_valid:
        fatal_error_list.append(error_msg)
    
    # Check for fixable errors
    if len(duplicates) > 0:
        fixable_error_list.append(f"Found {len(duplicates)} duplicate sequence(s)")

    # Check for missing coordinates
    missing_coords = fixable_errors.find_missing_coordinates(sequence_entries)
    if missing_coords:
        fixable_error_list.append(f"Found {len(missing_coords)} sequence(s) missing coordinates")

    # Check for overlapping sequences
    overlapping = fixable_errors.find_overlapping_sequences(sequence_entries)
    if overlapping:
        fixable_error_list.append(f"Found {len(overlapping)} overlapping sequence(s)")

    # Check for warnings (use merged GC annotations for interleaved support)
    has_warning, warning_msgs = warnings.check_ss_cons_from_parsed(
        parsed['gc_annotations'], sequence_entries
    )
    if has_warning:
        if isinstance(warning_msgs, list):
            warning_list.extend(warning_msgs)
        else:
            warning_list.append(warning_msgs)
    
    has_warning, warning_msgs = warnings.check_line_length(lines)
    if has_warning:
        warning_list.extend(warning_msgs)
    
    # Determine if file can be fixed
    can_be_fixed = len(fatal_error_list) == 0 and len(fixable_error_list) > 0
    
    return {
        'is_valid': len(fatal_error_list) == 0,
        'fatal_errors': fatal_error_list,
        'fixable_errors': fixable_error_list,
        'warnings': warning_list,
        'can_be_fixed': can_be_fixed,
        'lines': lines
    }


def find_cm_db(filepath):
    """Auto-detect Rfam.cm in the same directory as the input file."""
    input_dir = Path(filepath).parent
    rfam_cm = input_dir / 'Rfam.cm'
    if rfam_cm.exists():
        return str(rfam_cm)
    return None


def filter_known_families(sequence_entries, cm_db, verbose=False, evalue_threshold=CMSCAN_MAX_EVALUE):
    """
    Check sequences against existing Rfam families using cmscan.

    Sequences with significant hits (E-value below threshold) are flagged
    as warnings but NOT removed. Weak/spurious hits are noted separately.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples
        cm_db: Path to Rfam CM database file
        verbose: Print progress
        evalue_threshold: E-value cutoff for warnings (default 1e-3; lower = stricter)

    Returns:
        list: hit_details list of warning strings for reporting
    """
    if not shutil.which('cmscan'):
        if verbose:
            print("  Warning: cmscan not found. Install Infernal to filter known Rfam families.")
        return []

    hit_details = []
    tmpdir = None

    try:
        tmpdir = tempfile.mkdtemp(prefix='rfam_cmscan_')

        # Write ungapped sequences as FASTA
        fasta_path = tmpdir + '/sequences.fasta'
        with open(fasta_path, 'w') as f:
            for name, seq in sequence_entries:
                ungapped = seq.replace('.', '').replace('-', '')
                f.write(f">{name}\n{ungapped}\n")

        # Run cmscan with --cut_ga to use each model's GA (gathering) threshold,
        # which matches how Rfam curates family membership
        tblout_path = tmpdir + '/hits.tbl'
        cmd = [
            'cmscan', '--cut_ga', '--noali', '--tblout', tblout_path,
            cm_db, fasta_path
        ]
        if verbose:
            print(f"  Scanning against known Rfam families ({cm_db}), using GA thresholds...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            if verbose:
                print(f"  Warning: cmscan failed: {result.stderr.strip()}")
            return []

        # Parse tblout — keep best hit per sequence (first occurrence)
        seen = set()
        warned = []
        with open(tblout_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) < 16:
                    continue
                # tblout format: target_name target_acc query_name query_acc ... E-value score ...
                family_name = fields[0]
                family_acc = fields[1]
                seq_name = fields[2]
                evalue_str = fields[15]
                score = fields[14]

                if seq_name in seen:
                    continue
                seen.add(seq_name)

                try:
                    evalue = float(evalue_str)
                except ValueError:
                    continue

                # --cut_ga already filters by GA threshold; this E-value check
                # is kept as a safety net for edge cases
                if evalue <= evalue_threshold:
                    warned.append(seq_name)
                    detail = f"  {seq_name}: matches {family_acc} ({family_name}), E-value={evalue_str} -- WARNING"
                    hit_details.append(detail)
                    if verbose:
                        print(detail)
                else:
                    if verbose:
                        print(f"  {seq_name}: weak hit to {family_acc} ({family_name}), E-value={evalue_str} -- kept")

        if verbose:
            if warned:
                print(f"  {len(warned)} sequence(s) match known Rfam families (see warnings above).")
            else:
                print("  No sequences significantly match known Rfam families.")

        return hit_details

    except Exception as e:
        if verbose:
            print(f"  Warning: cmscan filtering failed: {e}")
        return []
    finally:
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)


def fix_file(filepath, output_mode='file', verbose=False, cm_db=None):
    """
    Fix fixable errors in a Stockholm file.

    Args:
        filepath: Path to input Stockholm file
        output_mode: 'stdout' to print to stdout, 'file' to create _corrected file
        verbose: Print progress messages
        cm_db: Path to Rfam CM database for filtering known families (auto-detected if None)

    Returns:
        tuple: (success, num_fixes_applied, output_path)
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return False, 0, None

    # Parse to check for missing coordinates
    parsed = parser.parse_stockholm_file(lines)
    sequence_entries = parsed['sequence_entries']
    missing_coords = fixable_errors.find_missing_coordinates(sequence_entries)

    num_fixes = 0
    id_mapping = {}
    to_remove = set()

    # Check sequences against known Rfam families (warning only, no removal)
    if cm_db:
        filter_known_families(sequence_entries, cm_db, verbose=verbose)

    # Fix missing coordinates if any
    if missing_coords:
        if verbose:
            print(f"  Fixing {len(missing_coords)} missing coordinates (downloading from NCBI)...")
        id_mapping, coords_to_remove = fixable_errors.fix_missing_coordinates(filepath, verbose=verbose)
        to_remove.update(coords_to_remove)
        # Count sequences that got new coordinates (where key != value)
        num_coords_fixed = len([k for k, v in id_mapping.items() if k != v])
        num_fixes += num_coords_fixed
        if verbose and to_remove:
            print(f"  {len(to_remove)} sequence(s) will be removed (no accurate match)")

    # Use parser output for all line collections (handles interleaved format)
    header_lines = parsed['header_lines']
    gf_lines = parsed['gf_lines']
    gs_lines = parsed['gs_lines']
    gr_annotations = parsed['gr_annotations']   # {(seq_name, feature): full_data}
    gc_annotations = parsed['gc_annotations']   # {feature: full_data}

    # Count within-block duplicates already removed by the parser
    seen_in_blocks = set()
    num_parser_dupes = 0
    for seq_name, seq_data, block_idx in parsed['raw_sequence_entries']:
        key = (seq_name, block_idx)
        if key in seen_in_blocks:
            num_parser_dupes += 1
        else:
            seen_in_blocks.add(key)
    if num_parser_dupes > 0 and verbose:
        print(f"  Removed {num_parser_dupes} duplicate sequence(s)")

    # Apply coordinate fixes, remove bad sequences, and remove duplicates
    processed_sequences = []  # List of (new_name, seq_content)
    seen_sequences = {}
    num_duplicates = num_parser_dupes
    num_removed = 0

    for seq_name, seq_content in sequence_entries:
        # Skip sequences marked for removal
        if seq_name in to_remove:
            num_removed += 1
            continue

        # Apply coordinate mapping if available
        new_name = id_mapping.get(seq_name, seq_name)

        # Check for duplicates
        accession, coords = fixable_errors.parse_sequence_identifier(new_name)
        key = (accession, coords, seq_content)

        if key not in seen_sequences:
            seen_sequences[key] = True
            processed_sequences.append((new_name, seq_content))
        else:
            num_duplicates += 1

    num_fixes += num_duplicates

    # Check for overlapping sequences and remove them
    overlapping = fixable_errors.find_overlapping_sequences(processed_sequences)
    if overlapping:
        if verbose:
            print(f"  Found {len(overlapping)} overlapping sequence(s) to remove")
        # Filter out overlapping sequences
        processed_sequences = [(name, seq) for name, seq in processed_sequences if name not in overlapping]
        num_removed += len(overlapping)

    # Validate sequences against NCBI (with BLAST fallback)
    if verbose:
        print(f"  Validating {len(processed_sequences)} sequence(s) against NCBI...")
    invalid_seqs, not_found_seqs, mismatched_seqs, blast_fixed = fixable_errors.validate_sequences_against_ncbi(processed_sequences, verbose=verbose)

    # Apply BLAST fixes: update sequence names with BLAST-found accession/coords
    if blast_fixed:
        if verbose:
            print(f"  {len(blast_fixed)} sequence(s) rescued via BLAST")
        processed_sequences = [
            (blast_fixed[name], seq) if name in blast_fixed else (name, seq)
            for name, seq in processed_sequences
        ]
        num_fixes += len(blast_fixed)
        # Merge BLAST renames into id_mapping so GR/GS annotations follow
        for orig_name, _ in sequence_entries:
            mapped = id_mapping.get(orig_name, orig_name)
            if mapped in blast_fixed:
                id_mapping[orig_name] = blast_fixed[mapped]

    if invalid_seqs:
        if not_found_seqs:
            print(f"  {len(not_found_seqs)} sequence(s) not found in NCBI (invalid accession?)")
        if mismatched_seqs:
            print(f"  {len(mismatched_seqs)} sequence(s) don't match NCBI data")
        processed_sequences = [(name, seq) for name, seq in processed_sequences if name not in invalid_seqs]
        num_removed += len(invalid_seqs)

    # Build set of kept sequence names (final names after all mappings)
    kept_seq_names = {name for name, _ in processed_sequences}
    # Also track which original names were kept (for GR/GS lookup)
    kept_original_names = set()
    for seq_name, _ in sequence_entries:
        mapped = id_mapping.get(seq_name, seq_name)
        if mapped in kept_seq_names:
            kept_original_names.add(seq_name)

    # Compute column width: max of sequence names, #=GR tags, and #=GC tags
    name_width = 30
    if processed_sequences:
        name_width = max(name_width, max(len(n) for n, _ in processed_sequences) + 1)
    for (gr_seq, gr_feat) in gr_annotations:
        mapped_name = id_mapping.get(gr_seq, gr_seq)
        if mapped_name in kept_seq_names:
            tag = f"#=GR {mapped_name} {gr_feat}"
            name_width = max(name_width, len(tag) + 1)
    for feature in gc_annotations:
        tag = f"#=GC {feature}"
        name_width = max(name_width, len(tag) + 1)

    # Build final corrected lines in proper Stockholm order
    corrected_lines = []

    # 1. Header (# STOCKHOLM 1.0)
    corrected_lines.extend(header_lines)

    # 2. File-level annotations (#=GF) - keep at top
    corrected_lines.extend(gf_lines)

    # 3. Sequence-level metadata (#=GS) for kept sequences
    for orig_seq_id, ann_line in gs_lines:
        if orig_seq_id is None:
            corrected_lines.append(ann_line)
        elif orig_seq_id in kept_original_names:
            corrected_lines.append(ann_line)

    # 4. Sequences with their #=GR annotations
    for new_name, seq_content in processed_sequences:
        corrected_lines.append(f"{new_name:<{name_width}} {seq_content}\n")
        # Output #=GR annotations for this sequence
        for (gr_seq, gr_feat), gr_data in gr_annotations.items():
            mapped_name = id_mapping.get(gr_seq, gr_seq)
            if mapped_name == new_name:
                tag = f"#=GR {new_name} {gr_feat}"
                corrected_lines.append(f"{tag:<{name_width}} {gr_data}\n")

    # 5. Column-level annotations (#=GC) - after sequences, aligned to same width
    for feature, data in gc_annotations.items():
        tag = f"#=GC {feature}"
        corrected_lines.append(f"{tag:<{name_width}} {data}\n")

    corrected_lines.append("//\n")

    if verbose and num_removed > 0:
        print(f"  Removed {num_removed} sequence(s) total from output")

    # Compute pairwise identity on final alignment
    if processed_sequences and verbose:
        pid_results = alignment_stats.compute_pairwise_identity(processed_sequences)
        if pid_results:
            avg_all = sum(r[1] for r in pid_results) / len(pid_results)
            print(f"\n  Pairwise identity (avg {avg_all:.1f}%):")
            for name, avg_id, _ in pid_results:
                flag = " <-- LOW" if avg_id < avg_all - 20 else ""
                print(f"    {name:<40} {avg_id:.1f}%{flag}")

    # Apply gap minimization: remove all-gap columns from final alignment
    if processed_sequences and shutil.which('esl-reformat'):
        orig_len = len(processed_sequences[0][1])
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sto', delete=False) as tmp:
                tmp_path = tmp.name
                tmp.writelines(corrected_lines)
            result = subprocess.run(
                ['esl-reformat', '--mingap', 'pfam', tmp_path],
                capture_output=True, text=True
            )
            if result.returncode == 0 and result.stdout.strip():
                mingap_lines = result.stdout.splitlines(keepends=True)
                # Determine new alignment length from first sequence line
                new_len = orig_len
                for line in mingap_lines:
                    stripped = line.strip()
                    if stripped and not stripped.startswith('#') and stripped != '//':
                        parts = stripped.split()
                        if len(parts) >= 2:
                            new_len = len(parts[1])
                            break
                cols_removed = orig_len - new_len
                corrected_lines = mingap_lines
                if verbose:
                    print(f"\n  Gap minimization (esl-reformat --mingap): removed {cols_removed} all-gap column(s), alignment now {new_len} columns wide")
            else:
                if verbose:
                    print(f"  Warning: esl-reformat --mingap failed: {result.stderr.strip()}")
        except Exception as e:
            if verbose:
                print(f"  Warning: esl-reformat --mingap failed: {e}")
        finally:
            if tmp_path:
                Path(tmp_path).unlink(missing_ok=True)
    elif processed_sequences and not shutil.which('esl-reformat'):
        if verbose:
            print("  Warning: esl-reformat not found in PATH, skipping gap minimization")

    if output_mode == 'stdout':
        for line in corrected_lines:
            print(line, end='')
        return True, num_fixes, None
    else:
        path = Path(filepath)
        output_path = path.parent / f"{path.stem}_corrected{path.suffix}"

        with open(output_path, 'w') as f:
            f.writelines(corrected_lines)

        return True, num_fixes, str(output_path)


class TeeOutput:
    """Write to both the original stdout and a buffer."""
    def __init__(self, original):
        self.original = original
        self.buffer = io.StringIO()

    def write(self, text):
        self.original.write(text)
        self.buffer.write(text)

    def flush(self):
        self.original.flush()

    def getvalue(self):
        return self.buffer.getvalue()


def main():
    """Main function to run validation from command line."""
    parser_arg = argparse.ArgumentParser(
        description='Validate Stockholm format MSA files'
    )
    parser_arg.add_argument(
        'files',
        nargs='+',
        help='Stockholm file(s) to validate'
    )
    parser_arg.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print detailed validation information'
    )
    parser_arg.add_argument(
        '--fix',
        action='store_true',
        help='Attempt to fix fixable errors'
    )
    parser_arg.add_argument(
        '--output-mode',
        choices=['stdout', 'file'],
        default='file',
        help='Output mode for fixed files: stdout or create _corrected file (default: file)'
    )
    parser_arg.add_argument(
        '--cm-db',
        default=None,
        help='Path to Rfam CM database (e.g., Rfam.cm) to filter sequences matching known families'
    )
    args = parser_arg.parse_args()
    
    all_valid = True
    
    for filepath in args.files:
        path = Path(filepath)

        if not path.exists():
            print(f"ERROR: File not found: {filepath}")
            all_valid = False
            continue

        # Capture output for report file when fixing in file mode
        tee = None
        if args.fix and args.output_mode == 'file':
            tee = TeeOutput(sys.stdout)
            sys.stdout = tee

        # Validate the file
        result = validate_stockholm_file(filepath)

        # Handle fixing if requested (always run to validate sequences against NCBI)
        if args.fix and result['is_valid']:
            success, num_fixes, output_path = fix_file(filepath, args.output_mode, verbose=args.verbose, cm_db=args.cm_db)
            if success:
                if args.output_mode == 'file':
                    print(f"✓ Fixed {num_fixes} issue(s) in {filepath}")
                    print(f"  Corrected file saved to: {output_path}")
                else:
                    # Already printed to stdout
                    pass

                # Re-validate the fixed content
                if args.output_mode == 'file':
                    result = validate_stockholm_file(output_path)
                    filepath = output_path

        # Display results
        if result['is_valid'] and not result['fixable_errors']:
            if args.verbose or result['warnings']:
                print(f"✓ {filepath}: Valid Stockholm file")

            # Display warnings
            for warning in result['warnings']:
                print(f"  ⚠ Warning: {warning}")

        elif result['is_valid'] and result['fixable_errors']:
            print(f"⚠ {filepath}: Valid but has fixable issues")
            for error in result['fixable_errors']:
                print(f"  • {error}")

            # Display warnings
            for warning in result['warnings']:
                print(f"  ⚠ Warning: {warning}")

            if not args.fix:
                print(f"  Tip: Use --fix to automatically correct these issues")

        else:
            print(f"✗ {filepath}: INVALID")
            for error in result['fatal_errors']:
                print(f"  - {error}")

            # Display fixable errors
            for error in result['fixable_errors']:
                print(f"  • Fixable: {error}")

            # Display warnings
            for warning in result['warnings']:
                print(f"  ⚠ Warning: {warning}")

            all_valid = False

        # Write report file
        if tee is not None:
            sys.stdout = tee.original
            report_path = path.parent / f"{path.stem}_Report.txt"
            version_file = Path(__file__).parent / 'VERSION'
            try:
                version = version_file.read_text().strip()
            except Exception:
                version = 'unknown'
            header = (
                f"rfam-seed-qc v{version}\n"
                f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                f"File: {path.name}\n"
                f"{'-' * 60}\n\n"
            )
            with open(report_path, 'w') as f:
                f.write(header)
                f.write(tee.getvalue())
            print(f"  Report saved to: {report_path}")
    
    if all_valid:
        print(f"\nAll {len(args.files)} file(s) passed validation.")
        return 0
    else:
        print(f"\nValidation failed for one or more files.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
