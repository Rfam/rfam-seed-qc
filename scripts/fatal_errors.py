"""
Fatal error validators for Stockholm files.

These errors cannot be automatically fixed and must terminate execution.
"""


def check_header(lines):
    """
    Check for the required Stockholm header line. A missing or unexpected
    version is only a warning (see stockholm_warnings.check_header_version).

    Args:
        lines: List of file lines

    Returns:
        tuple: (is_valid, error_message, line_number)
    """
    for line in lines:
        if line.strip().startswith('# STOCKHOLM'):
            return True, None, None
    return False, "Missing required '# STOCKHOLM 1.0' header", None


def check_terminator(lines):
    """
    Check for required terminator.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (is_valid, error_message)
    """
    for line in lines:
        if line.strip() == '//':
            return True, None
    return False, "Missing required '//' terminator"


def check_empty_file(lines):
    """
    Check if file is empty.
    
    Args:
        lines: List of file lines
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not lines:
        return False, "File is empty"
    return True, None


def check_no_sequences(sequences):
    """
    Check if alignment has sequences.
    
    Args:
        sequences: Dictionary of sequences
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not sequences:
        return False, "No sequences found in alignment"
    return True, None


def check_sequence_lengths(sequences):
    """
    Check that all sequences have the same length.
    
    Args:
        sequences: Dictionary of {seq_name: seq_data}
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not sequences:
        return True, None
    
    lengths = {name: len(seq) for name, seq in sequences.items()}
    unique_lengths = set(lengths.values())
    
    if len(unique_lengths) > 1:
        # Show only first 5 sequences with their lengths for readability
        length_items = list(lengths.items())[:5]
        length_sample = {name: length for name, length in length_items}
        suffix = f" (showing first {len(length_sample)} of {len(lengths)})" if len(lengths) > 5 else ""
        return False, f"Sequences have inconsistent lengths: {length_sample}{suffix}"
    return True, None


def check_sequence_characters(sequences):
    """
    Check that sequences contain only valid characters (no whitespace).
    
    Args:
        sequences: Dictionary of {seq_name: seq_data}
        
    Returns:
        tuple: (is_valid, error_message)
    """
    for seq_name, seq_data in sequences.items():
        # Check for whitespace
        if any(c.isspace() for c in seq_data):
            return False, f"Sequence '{seq_name}' contains whitespace characters"
    return True, None
