"""
Warning validators for Stockholm files.

These are non-critical issues that don't prevent validation.
"""

from scripts import structure_qc


def check_ss_cons_from_parsed(gc_annotations, sequence_entries):
    """
    Check for 2D structure consensus annotation using pre-parsed merged data.

    Args:
        gc_annotations: dict of {feature: merged_data} from parser
        sequence_entries: list of (seq_name, full_seq_data) tuples (merged)

    Returns:
        tuple: (has_warning, warning_message_or_list)
    """
    ss_cons_data = gc_annotations.get('SS_cons')

    if ss_cons_data is None:
        return True, "Missing 2D structure consensus annotation (#=GC SS_cons)"

    warnings = []

    # Check length matches sequences
    ss_len = len(ss_cons_data)
    if sequence_entries:
        expected_len = len(sequence_entries[0][1])
        if ss_len != expected_len:
            warnings.append(
                f"#=GC SS_cons length ({ss_len}) does not match "
                f"sequence length ({expected_len})"
            )

    # Validate WUSS notation characters and bracket balancing
    ss_warnings = validate_ss_cons_format(ss_cons_data)
    warnings.extend(ss_warnings)

    if warnings:
        return True, warnings
    return False, None


def check_ss_cons(lines):
    """
    Check for 2D structure consensus annotation: existence and length.

    Args:
        lines: List of file lines

    Returns:
        tuple: (has_warning, warning_message_or_list)
            Returns a single string or a list of strings for multiple warnings.
    """
    ss_cons_data = None
    seq_lengths = []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith('#=GC SS_cons'):
            parts = stripped.split(None, 2)
            if len(parts) >= 3:
                ss_cons_data = parts[2]
        elif not stripped.startswith('#') and stripped and stripped != '//':
            parts = stripped.split(None, 1)
            if len(parts) >= 2:
                seq_lengths.append(len(parts[1].replace(' ', '')))

    if ss_cons_data is None:
        return True, "Missing 2D structure consensus annotation (#=GC SS_cons)"

    warnings = []

    # Check length matches sequences
    ss_len = len(ss_cons_data)
    if seq_lengths:
        expected_len = seq_lengths[0]
        if ss_len != expected_len:
            warnings.append(
                f"#=GC SS_cons length ({ss_len}) does not match "
                f"sequence length ({expected_len})"
            )

    # Validate WUSS notation characters and bracket balancing
    ss_warnings = validate_ss_cons_format(ss_cons_data)
    warnings.extend(ss_warnings)

    if warnings:
        return True, warnings
    return False, None


def validate_ss_cons_format(ss_cons_data):
    """
    Validate that #=GC SS_cons uses valid WUSS notation and has balanced brackets.
    Pseudoknots may use any letter pair (Aa ... Zz), as in Easel.

    Args:
        ss_cons_data: The SS_cons annotation string

    Returns:
        list: Warning messages (empty if valid)
    """
    parsed = structure_qc.parse_wuss(ss_cons_data)
    warnings = []

    if parsed['invalid']:
        chars_str = ', '.join(
            f"'{ch}' (position {pos})" for ch, pos in parsed['invalid'].items()
        )
        warnings.append(f"#=GC SS_cons contains invalid characters: {chars_str}")

    for ch, pos in parsed['unmatched_close']:
        warnings.append(f"#=GC SS_cons has unmatched '{ch}' at position {pos}")

    for opener, count in parsed['unmatched_open'].items():
        closer = structure_qc.WUSS_OPEN_TO_CLOSE[opener]
        warnings.append(
            f"#=GC SS_cons has {count} unmatched '{opener}' (missing '{closer}')"
        )

    return warnings


def check_line_length(lines, max_length=10000):
    """
    Check for lines exceeding Rfam alignment limits.
    
    Args:
        lines: List of file lines
        max_length: Maximum allowed line length
        
    Returns:
        tuple: (has_warning, warning_messages_list)
    """
    warnings = []
    for i, line in enumerate(lines):
        if len(line) > max_length:
            warnings.append(f"Line {i+1} exceeds {max_length} character limit (length: {len(line)})")
    
    if warnings:
        return True, warnings
    return False, []
