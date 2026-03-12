"""
Fixable error handlers for Stockholm files.

These errors can be automatically corrected.
"""

import os
import time
import re
import warnings

# Suppress BioPython warnings
warnings.filterwarnings("ignore", module="Bio")
warnings.filterwarnings("ignore", message=".*FASTA.*")

# Import configuration
from scripts import config


def parse_sequence_identifier(seq_name):
    """
    Parse sequence identifier to extract accession and coordinates.
    
    Format: ACCESSION/START-END (e.g., AF228364.1/1-74)
    
    Args:
        seq_name: Sequence identifier string
        
    Returns:
        tuple: (accession, coordinates) or (seq_name, None) if no coordinates found
    """
    if '/' in seq_name:
        parts = seq_name.split('/', 1)
        return parts[0], parts[1]
    return seq_name, None


def find_missing_coordinates(sequence_entries):
    """
    Find sequences that are missing start/stop coordinates.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples

    Returns:
        list: Sequence names missing coordinates
    """
    missing = []
    for seq_name, _ in sequence_entries:
        _, coords = parse_sequence_identifier(seq_name)
        if coords is None:
            missing.append(seq_name)
    return missing


def get_fasta_file(identifier, max_retries=3):
    """
    Download fasta file from NCBI if not already downloaded.
    Retries on transient network failures.
    """
    script_location = os.path.dirname(os.path.abspath(__file__))
    fasta_folder = os.path.join(script_location, 'fasta')
    os.makedirs(fasta_folder, exist_ok=True)
    fasta_file = os.path.join(fasta_folder, identifier + '.fasta')
    if not os.path.exists(fasta_file) or os.stat(fasta_file).st_size == 0:
        import subprocess
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={identifier}&rettype=fasta&retmode=text"
        for attempt in range(1, max_retries + 1):
            time.sleep(config.NCBI_REQUEST_DELAY)
            result = subprocess.run(['curl', '-s', '-o', fasta_file, url])
            if result.returncode == 0 and os.path.exists(fasta_file) and os.stat(fasta_file).st_size > 0:
                break
            if attempt < max_retries:
                time.sleep(2 * attempt)  # backoff before retry
            else:
                raise RuntimeError(f"Failed to download FASTA for {identifier} after {max_retries} attempts")
    return fasta_file


def get_accession_version(fasta_file):
    """Get version from downloaded fasta file header."""
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line.split()[0][1:]  # Remove '>'
                if '.' in name:
                    return '.' + name.split('.')[1]
    return '.1'


def blast_search(sequence, verbose=False, min_identity=None, min_coverage=None, max_evalue=None):
    """
    Perform BLAST search for a sequence that couldn't be found directly.

    Args:
        sequence: DNA/RNA sequence string (without gaps)
        verbose: Print progress
        min_identity: Minimum identity percentage (default from config)
        min_coverage: Minimum query coverage percentage (default from config)
        max_evalue: Maximum e-value to accept (default from config)

    Returns:
        tuple: (accession_with_coords, is_accurate) or (None, False) if no good hit
    """
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio.Seq import Seq

    # Use config defaults if not specified
    if min_identity is None:
        min_identity = config.BLAST_MIN_IDENTITY
    if min_coverage is None:
        min_coverage = config.BLAST_MIN_COVERAGE
    if max_evalue is None:
        max_evalue = config.BLAST_MAX_EVALUE

    # Convert RNA to DNA for BLAST
    dna_seq = str(Seq(sequence).back_transcribe())

    if verbose:
        print(f"    Running BLAST search ({len(dna_seq)} bp)...")

    try:
        # Use blastn against nt database, limit results for speed
        result_handle = NCBIWWW.qblast(
            "blastn", "nt", dna_seq,
            hitlist_size=config.BLAST_HITLIST_SIZE,
            expect=max_evalue,
            format_type="XML"
        )
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                # Calculate identity and coverage
                identity = (hsp.identities / hsp.align_length) * 100
                coverage = (hsp.query_end - hsp.query_start) / len(dna_seq) * 100

                if identity >= min_identity and coverage >= min_coverage and hsp.expect <= max_evalue:
                    # Extract accession from hit - try multiple formats
                    hit_def = alignment.title
                    accession = None

                    # Format: "gi|xxx|ref|NZ_CP045771.1|" or "gi|xxx|gb|CP041614.1|"
                    accession_match = re.search(r'\|(ref|gb|emb|dbj)\|([^|]+)\|', hit_def)
                    if accession_match:
                        accession = accession_match.group(2).rstrip('|')
                    else:
                        # Try format: "NZ_CP045771.1 description"
                        accession_match = re.search(r'^([A-Z]{1,2}_?[A-Z]*\d+\.\d+)', alignment.hit_id)
                        if accession_match:
                            accession = accession_match.group(1)

                    if accession:
                        # Clean up accession (remove trailing |)
                        accession = accession.rstrip('|')
                        # Coordinates based on strand
                        coords = f"{hsp.sbjct_start}-{hsp.sbjct_end}"

                        new_accession = f"{accession}/{coords}"

                        if verbose:
                            print(f"    BLAST hit: {new_accession} (identity={identity:.1f}%, coverage={coverage:.1f}%, e={hsp.expect:.2e})")

                        return new_accession, True

        if verbose:
            print(f"    No accurate BLAST hit found")
        return None, False

    except Exception as e:
        if verbose:
            print(f"    BLAST search failed: {e}")
        return None, False


def fix_missing_coordinates(filepath, verbose=False, use_blast_fallback=True):
    """
    Fix missing coordinates by looking up sequences in NCBI.
    Falls back to BLAST search if direct lookup fails.

    Args:
        filepath: Path to Stockholm file
        verbose: Print progress messages
        use_blast_fallback: Use BLAST search for sequences not found directly

    Returns:
        tuple: (id_mapping dict, sequences_to_remove set)
    """
    from Bio import AlignIO, SeqIO
    from Bio.Seq import Seq

    mapping = {}
    to_remove = set()

    align = AlignIO.read(filepath, 'stockholm')
    for record in align:
        # Skip if already has coordinates
        _, coords = parse_sequence_identifier(record.id)
        if coords is not None:
            mapping[record.id] = record.id
            continue

        # Remove gaps for sequence matching
        ungapped = str(record.seq).replace('.', '').replace('-', '').upper()

        fasta_file = get_fasta_file(record.id)
        found = False

        try:
            fasta = SeqIO.read(fasta_file, "fasta")

            # Get version if not in record id
            if '.' not in record.id:
                version = get_accession_version(fasta_file)
            else:
                version = ''

            # Convert RNA to DNA for matching
            s1 = Seq(ungapped).back_transcribe()

            start = fasta.seq.find(s1)
            if start == -1:  # Try reverse complement
                s2 = s1.reverse_complement()
                start = fasta.seq.find(s2)
                if start != -1:
                    new_accession = f'{record.id}{version}/{start + len(s1)}-{start + 1}'
                    found = True
            else:
                new_accession = f'{record.id}{version}/{start + 1}-{start + len(s1)}'
                found = True

            if found:
                if verbose:
                    print(f"  {record.id} -> {new_accession}")
                mapping[record.id] = new_accession
                continue

        except Exception as e:
            if verbose:
                print(f"  Warning: Could not read fasta for {record.id}: {e}")

        # Direct lookup failed, try BLAST if enabled
        if not found and use_blast_fallback:
            if verbose:
                print(f"  {record.id} not found directly, trying BLAST...")

            blast_result, is_accurate = blast_search(ungapped, verbose=verbose)

            if is_accurate and blast_result:
                if verbose:
                    print(f"  {record.id} -> {blast_result} (via BLAST)")
                mapping[record.id] = blast_result
                continue

        # Neither direct lookup nor BLAST worked - mark for removal
        if verbose:
            print(f"  {record.id} will be REMOVED (no accurate match found)")
        to_remove.add(record.id)

    return mapping, to_remove


def check_overlap(start1, end1, start2, end2):
    """
    Check if two coordinate ranges overlap by at least 1 bp.
    Handles both forward (start < end) and reverse (start > end) strands.
    """
    # Normalize coordinates so min is always first
    a1, a2 = min(start1, end1), max(start1, end1)
    b1, b2 = min(start2, end2), max(start2, end2)
    # Check for overlap
    return a1 <= b2 and b1 <= a2


def find_overlapping_sequences(sequence_entries):
    """
    Find sequences from the same accession that overlap.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples

    Returns:
        set: Sequence names to remove due to overlaps
    """
    import random
    from collections import defaultdict

    # Group sequences by accession
    accession_groups = defaultdict(list)
    for seq_name, seq_data in sequence_entries:
        accession, coords = parse_sequence_identifier(seq_name)
        if coords:
            # Parse coordinates (handle both start-end formats)
            match = re.match(r'(\d+)-(\d+)', coords)
            if match:
                start, end = int(match.group(1)), int(match.group(2))
                accession_groups[accession].append((seq_name, start, end, seq_data))

    to_remove = set()

    # Check each accession group for overlaps
    for accession, sequences in accession_groups.items():
        if len(sequences) < 2:
            continue

        # Use greedy interval scheduling to keep maximum non-overlapping sequences
        # Sort by end position (using min of start/end for reverse strand)
        sorted_seqs = sorted(sequences, key=lambda x: min(x[1], x[2]))

        kept = []
        for seq_name, start, end, seq_data in sorted_seqs:
            # Check if this sequence overlaps with any kept sequence
            overlaps = False
            for kept_name, kept_start, kept_end, kept_data in kept:
                if check_overlap(start, end, kept_start, kept_end):
                    overlaps = True
                    break

            if not overlaps:
                kept.append((seq_name, start, end, seq_data))
            else:
                to_remove.add(seq_name)

    return to_remove


def validate_sequences_against_ncbi(sequence_entries, verbose=False, use_blast_fallback=True):
    """
    Validate that sequences match their NCBI source at the given coordinates.
    When validation fails, optionally falls back to BLAST to find the correct
    accession and coordinates.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples
        verbose: Print progress messages
        use_blast_fallback: Use BLAST search for sequences not found or mismatched

    Returns:
        tuple: (invalid set, not_found set, mismatched set, blast_fixed dict)
            - invalid: all sequences that failed validation (after BLAST attempts)
            - not_found: sequences where accession wasn't found in NCBI (and BLAST failed)
            - mismatched: sequences where data doesn't match NCBI (and BLAST failed)
            - blast_fixed: dict mapping old seq_name to new accession/coords from BLAST
    """
    from Bio import SeqIO
    from Bio.Seq import Seq

    not_found = set()
    mismatched = set()
    blast_fixed = {}

    for seq_name, seq_data in sequence_entries:
        accession, coords = parse_sequence_identifier(seq_name)

        if not coords:
            # No coordinates, can't validate
            continue

        # Parse coordinates
        match = re.match(r'(\d+)-(\d+)', coords)
        if not match:
            continue

        start, end = int(match.group(1)), int(match.group(2))

        # Remove gaps from Stockholm sequence (needed for both validation and BLAST)
        stockholm_seq = seq_data.replace('.', '').replace('-', '').upper()

        # Get FASTA file
        ncbi_ok = False
        try:
            fasta_file = get_fasta_file(accession)
            fasta = SeqIO.read(fasta_file, "fasta")

            # Extract subsequence from NCBI (handle forward and reverse strand)
            # Convert to 0-based indexing
            if start < end:
                # Forward strand
                ncbi_subseq = str(fasta.seq[start - 1:end]).upper()
            else:
                # Reverse strand
                ncbi_subseq = str(fasta.seq[end - 1:start].reverse_complement()).upper()

            # Convert RNA to DNA for comparison (NCBI stores DNA)
            stockholm_seq_dna = str(Seq(stockholm_seq).back_transcribe())

            # Compare
            if stockholm_seq_dna == ncbi_subseq:
                ncbi_ok = True
                if verbose:
                    print(f"  {seq_name}: validated OK")
            else:
                if verbose:
                    print(f"  {seq_name}: MISMATCH with NCBI")
                    if len(stockholm_seq_dna) != len(ncbi_subseq):
                        print(f"    Length: Stockholm={len(stockholm_seq_dna)}, NCBI={len(ncbi_subseq)}")
                    else:
                        # Show first difference
                        for i, (a, b) in enumerate(zip(stockholm_seq_dna, ncbi_subseq)):
                            if a != b:
                                print(f"    First difference at position {i+1}: Stockholm={a}, NCBI={b}")
                                break

        except Exception as e:
            if verbose:
                print(f"  {seq_name}: accession NOT FOUND in NCBI ({e})")

        if ncbi_ok:
            continue

        # Validation failed — try BLAST fallback
        if use_blast_fallback:
            if verbose:
                print(f"    Trying BLAST fallback for {seq_name}...")
            blast_result, is_accurate = blast_search(stockholm_seq, verbose=verbose)
            if is_accurate and blast_result:
                # Re-validate the Stockholm sequence against the BLAST hit's accession
                # to ensure the sequence data actually matches (required for MD5 check)
                blast_accession, blast_coords = parse_sequence_identifier(blast_result)
                blast_match = re.match(r'(\d+)-(\d+)', blast_coords) if blast_coords else None
                if blast_match:
                    blast_start, blast_end = int(blast_match.group(1)), int(blast_match.group(2))
                    try:
                        blast_fasta_file = get_fasta_file(blast_accession)
                        blast_fasta = SeqIO.read(blast_fasta_file, "fasta")
                        if blast_start < blast_end:
                            blast_ncbi_seq = str(blast_fasta.seq[blast_start - 1:blast_end]).upper()
                        else:
                            blast_ncbi_seq = str(blast_fasta.seq[blast_end - 1:blast_start].reverse_complement()).upper()
                        stockholm_seq_dna = str(Seq(stockholm_seq).back_transcribe())
                        if stockholm_seq_dna == blast_ncbi_seq:
                            if verbose:
                                print(f"    {seq_name} -> {blast_result} (via BLAST, re-validated OK)")
                            blast_fixed[seq_name] = blast_result
                            continue
                        else:
                            if verbose:
                                print(f"    BLAST hit {blast_result} found but sequence still does not match NCBI — removing")
                    except Exception:
                        if verbose:
                            print(f"    Could not re-validate BLAST hit {blast_result} — removing")

        # Neither direct validation nor BLAST succeeded
        # Determine which category this failure belongs to
        try:
            # If we got here via the except branch, it's not_found
            fasta_file = get_fasta_file(accession)
            SeqIO.read(fasta_file, "fasta")
            # If we can read it, it was a mismatch
            mismatched.add(seq_name)
        except Exception:
            not_found.add(seq_name)

    invalid = not_found | mismatched
    return invalid, not_found, mismatched, blast_fixed


def find_duplicates_from_entries(sequence_entries, raw_entries=None):
    """
    Find duplicate sequences from a list of sequence entries.

    When raw_entries is provided (block-aware mode), duplicates are detected as
    the same (accession, coords, seq_data) appearing multiple times within the
    SAME block. Same name across different blocks is interleaved, not a duplicate.

    Args:
        sequence_entries: List of (seq_name, seq_data) tuples (merged if interleaved)
        raw_entries: Optional list of (seq_name, seq_data, block_index) from parser

    Returns:
        tuple: (unique_sequences_dict, duplicate_indices)
    """
    if raw_entries is not None:
        return _find_duplicates_block_aware(sequence_entries, raw_entries)

    # Legacy mode: same as before
    seen_keys = {}
    unique_sequences = {}
    duplicate_indices = []

    for idx, (seq_name, seq_data) in enumerate(sequence_entries):
        accession, coords = parse_sequence_identifier(seq_name)
        key = (accession, coords, seq_data)

        if key not in seen_keys:
            seen_keys[key] = idx
            unique_sequences[seq_name] = seq_data
        else:
            duplicate_indices.append(idx)

    return unique_sequences, duplicate_indices


def _find_duplicates_block_aware(sequence_entries, raw_entries):
    """
    Detect true duplicates: same (accession, coords, data) appearing multiple
    times within the SAME block. Repeated names across different blocks are
    interleaved format, not duplicates.
    """
    from collections import defaultdict

    # For each block, track (accession, coords, seq_data) keys
    per_block_seen = defaultdict(lambda: defaultdict(int))
    duplicate_names = set()

    for seq_name, seq_data, block_idx in raw_entries:
        accession, coords = parse_sequence_identifier(seq_name)
        key = (accession, coords, seq_data)
        per_block_seen[block_idx][key] += 1
        if per_block_seen[block_idx][key] > 1:
            duplicate_names.add(seq_name)

    # Build unique_sequences from merged entries
    unique_sequences = {}
    duplicate_indices = []

    for idx, (seq_name, seq_data) in enumerate(sequence_entries):
        unique_sequences[seq_name] = seq_data
        if seq_name in duplicate_names:
            duplicate_indices.append(idx)

    return unique_sequences, duplicate_indices


def remove_duplicates(lines):
    """
    Remove duplicate sequences from Stockholm file lines.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (corrected_lines, num_duplicates_removed)
    """
    # Separate different types of lines
    header_lines = []
    annotation_lines = []
    footer_lines = []
    sequence_entries = []  # List of (seq_name, seq_content) tuples
    
    sequences_started = False
    sequences_ended = False
    
    for line in lines:
        stripped = line.strip()
        
        if stripped == '//':
            sequences_ended = True
            footer_lines.append(line)
        elif sequences_ended:
            footer_lines.append(line)
        elif stripped.startswith('# STOCKHOLM'):
            header_lines.append(line)
            sequences_started = True
        elif stripped.startswith('#=GF') or stripped.startswith('#=GC') or stripped.startswith('#=GS') or stripped.startswith('#=GR'):
            annotation_lines.append(line)
        elif stripped.startswith('#'):
            if not sequences_started:
                header_lines.append(line)
            else:
                annotation_lines.append(line)
        elif stripped:
            # This is a sequence line
            parts = stripped.split(None, 1)
            if len(parts) >= 2:
                seq_name = parts[0]
                seq_content = parts[1].replace(' ', '')
                sequence_entries.append((seq_name, seq_content))
        else:
            # Empty line
            if not sequences_started:
                header_lines.append(line)
    
    # Find duplicates
    seen_keys = {}
    unique_entries = []
    duplicates_count = 0
    
    for seq_name, seq_content in sequence_entries:
        accession, coords = parse_sequence_identifier(seq_name)
        key = (accession, coords, seq_content)
        
        if key not in seen_keys:
            seen_keys[key] = True
            unique_entries.append((seq_name, seq_content))
        else:
            duplicates_count += 1
    
    # Reconstruct file
    corrected_lines = []
    
    # Add header lines
    corrected_lines.extend(header_lines)
    
    # Calculate max sequence name length for formatting
    max_name_len = max(len(name) for name, _ in unique_entries) if unique_entries else 30
    name_width = max(max_name_len + 2, 30)
    
    # Add unique sequences
    for seq_name, seq_content in unique_entries:
        corrected_lines.append(f"{seq_name:<{name_width}} {seq_content}\n")
    
    # Add annotation lines
    corrected_lines.extend(annotation_lines)
    
    # Add footer lines
    corrected_lines.extend(footer_lines)
    
    return corrected_lines, duplicates_count
