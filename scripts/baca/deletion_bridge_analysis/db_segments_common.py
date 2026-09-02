"""
db_segments_common.py — shared segment-identification logic for the
Deletion Bridge Analysis (see DELETION_BRIDGE_ANALYSIS_PLAN.md at the
project root — read that file before touching this module; it is the
canonical, locked spec, not this docstring).

A "DB pair" is one distinct (lower_bp, higher_bp) deletion-bridge anchor
pair within a single mmc5 chain. A chain can have zero, one, or many
distinct DB pairs (verified 2026-08-13: 194/366 chains have >=1, 126/366
have >1, max 19 in a single chain, all 521 DB pairs cohort-wide are
same-chromosome) — each DB pair is treated completely independently per
the locked plan decision #2 (own letter sequence restarting at A, own
Part 3 row).

For a given DB pair, the "segments" are the intact reference stretches
between every consecutive pair of breakpoints (same chromosome, position
order) from the low anchor to the high anchor, inclusive. Labeled A, B,
C, ... in position order (max observed internal breakpoints per pair is
9, i.e. at most 10 segments — well under the 26-letter alphabet, so
overflow is not handled here, only guarded against).

This module is intentionally reused (not duplicated) by Part 1's
visualization script AND, later, Part 2 — segment identity must be
IDENTICAL between Part 1 and Part 2 (plan decision #3), which is only
guaranteed if both parts call the same functions.
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'mmc5'))
from baca_chain_visualization import load_mmc5  # noqa: E402  (reuse, don't duplicate loading)


def parse_deletion_bridge(v):
    if pd.isna(v):
        return None
    s = str(v).strip()
    return int(float(s)) if s else None


def find_db_pairs_for_chain(chain_rows):
    """
    Returns a sorted list of (low_bp, high_bp) tuples — one per distinct
    deletion-bridge pair in this chain, deduplicated (a pair recorded
    from either breakpoint's own row collapses to one entry). `low_bp`/
    `high_bp` are ordered by GENOMIC POSITION, not by breakpoint-number
    value (breakpoint numbers are just mmc5 row IDs and do not track
    position order — verified this matters: e.g. P09-1042 chain 1 has a
    pair BP18/BP28 where BP18 is at the HIGHER position). Pairs
    themselves are sorted by their low anchor's position, for
    deterministic DB1/DB2/... numbering.
    """
    bp_to_pos = dict(zip(chain_rows['Breakpoint number'].astype(int), chain_rows['position']))
    pairs = set()
    for _, r in chain_rows.iterrows():
        bp = int(r['Breakpoint number'])
        delb = parse_deletion_bridge(r['Deletion bridge partner breakpoint'])
        if delb is None or delb not in bp_to_pos:
            continue
        a, b = sorted([bp, delb], key=lambda x: bp_to_pos[x])
        pairs.add((a, b))
    return sorted(pairs, key=lambda p: (bp_to_pos[p[0]], bp_to_pos[p[1]]))


def letter_for_index(i):
    """0->'A', 1->'B', ..., 25->'Z'. Raises past 26 rather than silently
    wrapping — max observed internal breakpoints per DB pair is 9 (=> at
    most 10 segments), so this should never actually fire; a real hit
    here means the data has changed and the plan doc's empirical basis
    needs re-checking, not a silent AA/AB-style extension."""
    if i >= 26:
        raise ValueError(f'segment index {i} exceeds 26-letter alphabet — '
                          f'unexpected per DELETION_BRIDGE_ANALYSIS_PLAN.md\'s '
                          f'empirical check (max 9 internal breakpoints observed); re-verify data')
    return chr(ord('A') + i)


def get_segments_for_db_pair(chain_rows, low_bp, high_bp, chromosome):
    """
    Ordered list of dicts, one per intact reference segment strictly
    spanning from low_bp to high_bp (inclusive), in position order:
      {'letter': 'A', 'low_bp': 6, 'high_bp': 18,
       'low_pos': 40938607.0, 'high_pos': 41028309.0}
    Native/reference orientation is always (low_bp)+ ---- (high_bp)-, by
    definition — independent of whichever strand that breakpoint's own
    real rearrangement (if any) actually used elsewhere.
    """
    same_chrom = chain_rows[chain_rows['chromosome'] == chromosome]
    lo_pos = same_chrom.loc[same_chrom['Breakpoint number'].astype(int) == low_bp, 'position'].values[0]
    hi_pos = same_chrom.loc[same_chrom['Breakpoint number'].astype(int) == high_bp, 'position'].values[0]
    if not lo_pos < hi_pos:
        raise ValueError(f'low_bp {low_bp} (pos {lo_pos}) must be strictly lower position than '
                          f'high_bp {high_bp} (pos {hi_pos})')

    span_rows = same_chrom[(same_chrom['position'] >= lo_pos) & (same_chrom['position'] <= hi_pos)]
    span_rows = span_rows.drop_duplicates(subset=['Breakpoint number']).sort_values('position')
    bps = span_rows['Breakpoint number'].astype(int).tolist()
    positions = span_rows['position'].tolist()

    segments = []
    for i in range(len(bps) - 1):
        segments.append({
            'letter': letter_for_index(i),
            'low_bp': bps[i], 'high_bp': bps[i + 1],
            'low_pos': positions[i], 'high_pos': positions[i + 1],
        })
    return segments


def segment_display_label(seg, db_index=None):
    """'DB1:A' style label when db_index is given (1-based; used whenever
    a chain has >1 DB pair, per the confirmed disambiguation convention),
    else plain 'A' (single-DB chains don't need the prefix)."""
    letter = seg['letter']
    return f'DB{db_index}:{letter}' if db_index is not None else letter


def segment_span_text(seg):
    """'spans 6[+]----18[-]' — native/reference orientation, always
    low_bp+ to high_bp-."""
    return f"spans {seg['low_bp']}[+]----{seg['high_bp']}[-]"


def get_all_db_pairs_with_segments(chain_rows):
    """
    Top-level entry point. Given one chain's full row-group (as produced
    by mmc5.groupby(['Individual', 'Chain number'])), returns a list of
    dicts, one per distinct DB pair, in DB1/DB2/... order:
      {'db_index': 1, 'low_bp': 6, 'high_bp': 14, 'chromosome': 8,
       'segments': [ ...from get_segments_for_db_pair... ]}
    Empty list if the chain has no deletion bridge at all (per plan
    decision #5, such chains are skipped entirely by Part 1/2).
    """
    bp_to_chrom = dict(zip(chain_rows['Breakpoint number'].astype(int), chain_rows['chromosome']))
    pairs = find_db_pairs_for_chain(chain_rows)
    result = []
    for i, (low_bp, high_bp) in enumerate(pairs, start=1):
        chrom_a, chrom_b = bp_to_chrom[low_bp], bp_to_chrom[high_bp]
        if chrom_a != chrom_b:
            # Verified empirically 2026-08-13: 0/521 DB pairs cohort-wide
            # are cross-chromosome. Guard rather than silently mishandle
            # if this ever turns up (plan doc assumes same-chromosome only).
            raise ValueError(f'cross-chromosome DB pair found: {low_bp}(chr{chrom_a}) / '
                              f'{high_bp}(chr{chrom_b}) — not expected, re-verify data')
        segments = get_segments_for_db_pair(chain_rows, low_bp, high_bp, chrom_a)
        result.append({
            'db_index': i, 'low_bp': low_bp, 'high_bp': high_bp,
            'chromosome': chrom_a, 'segments': segments,
        })
    return result


def chain_has_db(chain_rows):
    """Cheap check used to filter which of the 366 chains Part 1/2 even
    process — per plan decision #5, no-DB chains are skipped entirely."""
    return len(find_db_pairs_for_chain(chain_rows)) > 0
