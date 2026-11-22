"""
Utility functions for parsing Stockholm files.

Handles both non-interleaved (single-block) and interleaved (multi-block)
Stockholm format. In interleaved format, sequences, #=GR, and #=GC annotations
are split across multiple blocks separated by blank lines; data is concatenated.
"""

from collections import OrderedDict


def parse_stockholm_file(lines):
    """
    Parse Stockholm file into structured components.

    Handles interleaved format by concatenating sequence data, #=GR annotations,
    and #=GC annotations across blocks. Block boundaries are blank lines.

    Args:
        lines: List of file lines

    Returns:
        dict: Parsed file structure with merged sequences, annotations, etc.
            - sequence_entries: list of (seq_name, full_seq_data) tuples (merged)
            - raw_sequence_entries: list of (seq_name, seq_data, block_index) for dedup
            - conflicting_rows: list of (seq_name, first_line, other_line): a name
              repeated within one block with different sequence data (1-based
              line numbers). Only the first row is kept in sequence_entries.
            - gr_annotations: OrderedDict {(seq_name, feature): full_data}
            - gc_annotations: OrderedDict {feature: full_data}
            - gf_lines: list of raw #=GF lines
            - gs_lines: list of (seq_name, line) tuples for #=GS lines
            - header_lines: list of header/comment lines
            - header_found, terminator_found, ss_cons_found: booleans
            - lines: original lines (kept for backward compatibility)
    """
    header_found = False
    terminator_found = False
    ss_cons_found = False

    # Merged data structures
    sequence_order = OrderedDict()   # seq_name -> full_seq_data
    gr_annotations = OrderedDict()   # (seq_name, feature) -> full_data
    gc_annotations = OrderedDict()   # feature -> full_data

    # Raw per-line entries for duplicate detection
    raw_sequence_entries = []        # [(seq_name, seq_data, block_index), ...]
    conflicting_rows = []            # [(seq_name, first_line, other_line), ...]

    # Other line collections
    gf_lines = []
    gs_lines = []
    header_lines = []

    block_index = 0
    seen_data_in_block = False       # Have we seen seq/GR/GC in this block?
    seen_in_current_block = {}       # Seq name -> (data, line number) in the current block
    seen_gr_in_current_block = set() # (seq_name, feature) GR keys seen in current block
    seen_gc_in_current_block = set() # GC feature names seen in current block

    for line_number, line in enumerate(lines, 1):
        stripped = line.strip()

        # Blank line: may indicate block boundary
        if not stripped:
            if seen_data_in_block:
                block_index += 1
                seen_data_in_block = False
                seen_in_current_block = {}
                seen_gr_in_current_block = set()
                seen_gc_in_current_block = set()
            continue

        # Stockholm header
        if stripped.startswith('# STOCKHOLM'):
            header_found = True
            header_lines.append(line)
            continue

        # Terminator
        if stripped == '//':
            terminator_found = True
            continue

        # #=GF file-level annotations
        if stripped.startswith('#=GF'):
            gf_lines.append(line)
            continue

        # #=GC column-level annotations (concatenate across blocks)
        if stripped.startswith('#=GC'):
            parts = stripped.split(None, 2)
            if len(parts) >= 3:
                feature = parts[1]
                data = parts[2].replace(' ', '')
                if feature == 'SS_cons':
                    ss_cons_found = True
                if feature in seen_gc_in_current_block:
                    # Same GC feature in same block = duplicate, skip
                    pass
                elif feature in gc_annotations:
                    gc_annotations[feature] += data
                else:
                    gc_annotations[feature] = data
                seen_gc_in_current_block.add(feature)
            seen_data_in_block = True
            continue

        # #=GS sequence-level metadata (not split across blocks)
        if stripped.startswith('#=GS'):
            parts = stripped.split(None, 2)
            if len(parts) >= 3:
                seq_name = parts[1]
                gs_lines.append((seq_name, line))
            else:
                gs_lines.append((None, line))
            continue

        # #=GR per-residue annotations (concatenate across blocks)
        if stripped.startswith('#=GR'):
            parts = stripped.split(None, 3)
            if len(parts) >= 4:
                seq_name = parts[1]
                feature = parts[2]
                data = parts[3].replace(' ', '')
                key = (seq_name, feature)
                if key in seen_gr_in_current_block:
                    # Same GR key in same block = duplicate, skip
                    pass
                elif key in gr_annotations:
                    gr_annotations[key] += data
                else:
                    gr_annotations[key] = data
                seen_gr_in_current_block.add(key)
            seen_data_in_block = True
            continue

        # Other comment lines
        if stripped.startswith('#'):
            header_lines.append(line)
            continue

        # Sequence line
        parts = stripped.split(None, 1)
        if len(parts) >= 2:
            seq_name = parts[0]
            seq_data = parts[1].replace(' ', '')

            raw_sequence_entries.append((seq_name, seq_data, block_index))

            if seq_name in seen_in_current_block:
                # Same name in same block: a duplicate row if the data is the
                # same; otherwise a conflict (which row is right is unknown)
                first_data, first_line = seen_in_current_block[seq_name]
                if seq_data != first_data:
                    conflicting_rows.append((seq_name, first_line, line_number))
            elif seq_name in sequence_order:
                # Same name in new block = interleaved, concatenate
                sequence_order[seq_name] += seq_data
            else:
                sequence_order[seq_name] = seq_data

            seen_in_current_block.setdefault(seq_name, (seq_data, line_number))
            seen_data_in_block = True

    # Build merged sequence_entries list preserving insertion order
    sequence_entries = list(sequence_order.items())

    return {
        'header_found': header_found,
        'terminator_found': terminator_found,
        'ss_cons_found': ss_cons_found,
        'sequence_entries': sequence_entries,
        'raw_sequence_entries': raw_sequence_entries,
        'conflicting_rows': conflicting_rows,
        'gr_annotations': gr_annotations,
        'gc_annotations': gc_annotations,
        'gf_lines': gf_lines,
        'gs_lines': gs_lines,
        'header_lines': header_lines,
        'lines': lines,
    }
