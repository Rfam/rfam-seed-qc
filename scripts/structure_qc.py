"""
Structure checks for the consensus secondary structure (#=GC SS_cons).

parse_wuss() reads WUSS notation (as Infernal/Easel do) into base pairs.
pair_support() then checks whether the aligned sequences can actually form
those pairs: Watson-Crick or G-U at both positions, without gaps. It does
not test covariation (that needs a tool such as R-scape); it finds pairs no
sequence supports and sequences that break many pairs, which often point to
a misaligned sequence or a wrong pair in SS_cons.
"""

import string

from scripts import config

# WUSS: brackets of four kinds, and pseudoknots written as letter pairs
# (A pairs with a, ... Z with z), as in Easel's esl_wuss.c
WUSS_OPEN_TO_CLOSE = {'<': '>', '(': ')', '[': ']', '{': '}'}
WUSS_OPEN_TO_CLOSE.update(zip(string.ascii_uppercase, string.ascii_lowercase))
WUSS_CLOSE_TO_OPEN = {close: open_ for open_, close in WUSS_OPEN_TO_CLOSE.items()}
WUSS_UNPAIRED = set(':,_-.~')
WUSS_VALID_CHARS = set(WUSS_OPEN_TO_CLOSE) | set(WUSS_CLOSE_TO_OPEN) | WUSS_UNPAIRED

CANONICAL_PAIRS = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
GAP_CHARS = set('.-~')

# Longest list of columns or sequences shown in one warning
MAX_LISTED = 10


def parse_wuss(ss_cons):
    """
    Parse WUSS notation into base pairs.

    Returns:
        dict with
          'pairs': sorted list of (i, j) column pairs, 0-based, i < j
          'invalid': {char: first 1-based position} for characters that are
                     not WUSS
          'unmatched_close': [(char, 1-based position)] closers with no opener
          'unmatched_open': {opener: count} openers with no closer
    """
    stacks = {opener: [] for opener in WUSS_OPEN_TO_CLOSE}
    pairs, invalid, unmatched_close = [], {}, []
    for i, ch in enumerate(ss_cons):
        if ch in WUSS_OPEN_TO_CLOSE:
            stacks[ch].append(i)
        elif ch in WUSS_CLOSE_TO_OPEN:
            stack = stacks[WUSS_CLOSE_TO_OPEN[ch]]
            if stack:
                pairs.append((stack.pop(), i))
            else:
                unmatched_close.append((ch, i + 1))
        elif ch not in WUSS_UNPAIRED:
            invalid.setdefault(ch, i + 1)
    return {
        'pairs': sorted(pairs),
        'invalid': invalid,
        'unmatched_close': unmatched_close,
        'unmatched_open': {o: len(s) for o, s in stacks.items() if s},
    }


def _residue(ch):
    """'A', 'C', 'G', 'U' (T read as U), 'gap', or None for other symbols
    (N and other ambiguity codes, which say nothing about pairing)."""
    if ch in GAP_CHARS:
        return 'gap'
    ch = ch.upper().replace('T', 'U')
    return ch if ch in 'ACGU' else None


def pair_support(sequence_entries, pairs):
    """
    Count, for each base pair, how many sequences pair canonically
    (Watson-Crick or G-U), pair non-canonically, or have a gap.

    Args:
        sequence_entries: [(name, aligned sequence)], all as long as SS_cons
        pairs: [(i, j)] 0-based columns, from parse_wuss

    Returns:
        dict with
          'pairs': [{'i', 'j', 'canonical', 'noncanonical', 'gapped', 'total'}]
          'sequences': [{'name', 'canonical', 'noncanonical', 'gapped'}]
    """
    pair_stats = [{'i': i, 'j': j, 'canonical': 0, 'noncanonical': 0, 'gapped': 0,
                   'total': len(sequence_entries)} for i, j in pairs]
    seq_stats = []
    for name, seq in sequence_entries:
        counts = {'name': name, 'canonical': 0, 'noncanonical': 0, 'gapped': 0}
        for stats in pair_stats:
            a, b = _residue(seq[stats['i']]), _residue(seq[stats['j']])
            if a == 'gap' or b == 'gap':
                kind = 'gapped'
            elif a is None or b is None:
                continue  # ambiguous base: neither support nor conflict
            elif a + b in CANONICAL_PAIRS:
                kind = 'canonical'
            else:
                kind = 'noncanonical'
            stats[kind] += 1
            counts[kind] += 1
        seq_stats.append(counts)
    return {'pairs': pair_stats, 'sequences': seq_stats}


def _columns(stats):
    return f"{stats['i'] + 1}-{stats['j'] + 1}"


def _listed(items):
    shown = ', '.join(items[:MAX_LISTED])
    return shown + (f", ... ({len(items) - MAX_LISTED} more)" if len(items) > MAX_LISTED else '')


def check_pair_support(ss_cons, sequence_entries):
    """
    Warnings for base pairs in SS_cons that the sequences do not support.
    Nothing is checked if SS_cons is not usable (wrong length, invalid or
    unbalanced notation): those problems are reported by the format checks.

    Returns:
        list of warning strings (empty if every pair is well supported)
    """
    if not ss_cons or not sequence_entries:
        return []
    if any(len(seq) != len(ss_cons) for _, seq in sequence_entries):
        return []
    parsed = parse_wuss(ss_cons)
    if parsed['invalid'] or parsed['unmatched_close'] or parsed['unmatched_open']:
        return []
    if not parsed['pairs']:
        return []

    result = pair_support(sequence_entries, parsed['pairs'])
    n_pairs = len(result['pairs'])
    min_support = config.STRUCTURE_MIN_PAIR_SUPPORT
    max_gaps = config.STRUCTURE_MAX_PAIR_GAPS
    max_mispaired = config.STRUCTURE_MAX_SEQ_MISPAIRED
    warnings = []

    weak = []
    for stats in result['pairs']:
        paired = stats['canonical'] + stats['noncanonical']
        if paired and stats['canonical'] / paired < min_support:
            weak.append(f"{_columns(stats)} ({stats['canonical']}/{paired} canonical)")
    if weak:
        warnings.append(
            f"Base-pair support: {len(weak)} of {n_pairs} SS_cons base pairs "
            f"{'is' if len(weak) == 1 else 'are'} Watson-Crick "
            f"or G-U in fewer than {min_support:.0%} of the sequences that have both bases "
            f"(columns): {_listed(weak)}")

    gappy = []
    for stats in result['pairs']:
        if stats['total'] and stats['gapped'] / stats['total'] > max_gaps:
            gappy.append(f"{_columns(stats)} ({stats['gapped']}/{stats['total']} gapped)")
    if gappy:
        warnings.append(
            f"Base-pair support: {len(gappy)} of {n_pairs} SS_cons base pairs "
            f"{'has' if len(gappy) == 1 else 'have'} a gap "
            f"in more than {max_gaps:.0%} of the sequences (columns): {_listed(gappy)}")

    mispaired = []
    for counts in result['sequences']:
        paired = counts['canonical'] + counts['noncanonical']
        # At least 2: one odd pair in a short structure is not a pattern
        if counts['noncanonical'] >= 2 and counts['noncanonical'] / paired > max_mispaired:
            mispaired.append(f"{counts['name']} ({counts['noncanonical']}/{paired})")
    if mispaired:
        warnings.append(
            f"Base-pair support: {len(mispaired)} sequence(s) cannot form more than "
            f"{max_mispaired:.0%} of the SS_cons base pairs they have both bases for "
            f"(possibly misaligned): {_listed(mispaired)}")

    return warnings
