"""
part1_db_segments_visualization.py — Part 1 of the Deletion Bridge
Analysis (see DELETION_BRIDGE_ANALYSIS_PLAN.md at the project root —
read that file first, it is the canonical, locked spec).

Re-draws every mmc5 chain that has >=1 deletion bridge (194/366; the
other 172 chains are skipped entirely per the plan's decision #5) in the
same visual style as scripts/baca/mmc5/baca_chain_visualization.py (all
3 real edge types, BP number/location/strand corner table), but with
each DB pair's internal segments additionally:
  - highlighted on the chain diagram as a colored underline directly
    below the relevant chromosome track (one color per DB pair, cycling
    through a fixed palette, so multiple DB pairs on the same
    chromosome stay visually distinguishable),
  - labeled with their letter (DB1:A, DB1:B, ... / DB2:A, ... when a
    chain has >1 DB pair, per the confirmed disambiguation convention),
  - listed in a new "Segment Span Information" corner table alongside
    the existing BP-location table — this table is the DEFINITIVE
    letter -> span mapping; the on-diagram highlight is a visual aid,
    not the source of truth (some very dense chains, e.g. P09-1042
    chain 1 with 14 DB pairs sharing one chromosome, will have visually
    crowded highlights — the table is always legible regardless).

Reuses (does not duplicate) baca_chain_visualization.py's data loading,
BP-number labeling, and BP-location-table helpers, and db_segments_common
for all segment identification. Output:
  results/Deletion Bridge Analysis on Baca Chains/
    Baca Chains With DB Segments/
      Baca Chains with Segments Images/{patient_id}_chain_{n}.png
"""

import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'mmc5'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))

from baca_chain_visualization import (  # noqa: E402  (reuse, don't duplicate)
    load_mmc5,
    parse_adjacent,
    parse_deletion_bridge,
    theta_notation,
    draw_bp_number_labels,
    _bp_table_dims,
    _add_bp_location_table,
    DELETION_BRIDGE_COLOR,
    ERG_PATTERN,
)
from baca_aberration_cycle_drawing import (  # noqa: E402
    build_norm_map,
    node_x,
    draw_chr_track,
    _arc,
    _assign_rads,
    COL,
    CHR_SEP,
)
from db_segments_common import (  # noqa: E402
    get_all_db_pairs_with_segments,
    segment_display_label,
    segment_span_text,
    chain_has_db,
)

RESULTS_ROOT = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
TOP_FOLDER = os.path.join(RESULTS_ROOT, "Deletion Bridge Analysis on Baca Chains")
PART1_FOLDER = os.path.join(TOP_FOLDER, "Baca Chains With DB Segments")
OUTPUT_DIR = os.path.join(PART1_FOLDER, "Baca Chains with Segments Images")

# Cycling color palette for DB pairs (tab10-style, distinct at a glance).
# Deliberately different from the existing rearrangement/adjacency/
# deletion-bridge edge colors so segment highlights don't get confused
# with real edges.
DB_PALETTE = [
    '#D81B60', '#1E88E5', '#FFC107', '#004D40', '#8E24AA',
    '#F4511E', '#00ACC1', '#43A047', '#6D4C41', '#5E35B1',
]

SEG_ROW_OFFSET = 0.85       # y-offset below a chromosome track for the segment-highlight row
SEG_LABEL_OFFSET = 1.18     # y-offset for the segment letter label text


def _db_color(db_index):
    return DB_PALETTE[(db_index - 1) % len(DB_PALETTE)]


def draw_baca_chain_with_segments(patient_id, chain_number, chain_rows, save_dir):
    db_pairs = get_all_db_pairs_with_segments(chain_rows)
    if not db_pairs:
        return None  # no-DB chains are skipped entirely (plan decision #5)

    use_db_prefix = len(db_pairs) > 1

    bp_to_node = {
        int(r['Breakpoint number']): (int(r['chromosome']), r['position'], r['node_strand'])
        for _, r in chain_rows.iterrows()
    }

    chr_to_pos = {
        int(chrom): sorted(grp['position'].unique())
        for chrom, grp in chain_rows.groupby('chromosome')
    }
    chr_to_rows = {int(chrom): grp for chrom, grp in chain_rows.groupby('chromosome')}
    sorted_chroms = sorted(chr_to_pos.keys())
    nmaps = {c: build_norm_map(chr_to_pos[c]) for c in sorted_chroms}
    chr_y = {c: (len(sorted_chroms) - 1 - i) * CHR_SEP for i, c in enumerate(sorted_chroms)}
    y_span = (len(sorted_chroms) - 1) * CHR_SEP

    # REARRANGEMENT edges
    rearr_arc_data = []
    n_unresolved_rearr = 0
    for rearr_num, grp in chain_rows.groupby('Rearrangement number'):
        bps = grp['Breakpoint number'].astype(int).tolist()
        if len(bps) != 2:
            n_unresolved_rearr += len(bps)
            continue
        node_a, node_b = bp_to_node[bps[0]], bp_to_node[bps[1]]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya, yb = chr_y[node_a[0]], chr_y[node_b[0]]
        rearr_arc_data.append((xa, ya, xb, yb, node_a[0] != node_b[0]))

    # ADJACENCY reference edges
    seen_ref_pairs = set()
    ref_arc_data = []
    for _, r in chain_rows.iterrows():
        bp_num = int(r['Breakpoint number'])
        node_a = bp_to_node[bp_num]
        for adj_num in parse_adjacent(r['Adjacent breakpoint(s)']):
            if adj_num not in bp_to_node:
                continue
            pair_key = tuple(sorted([bp_num, adj_num]))
            if pair_key in seen_ref_pairs:
                continue
            seen_ref_pairs.add(pair_key)
            node_b = bp_to_node[adj_num]
            xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
            xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
            ya = chr_y[node_a[0]]
            ref_arc_data.append((xa, ya, xb, ya))

    # DELETION BRIDGE edges
    seen_delb_pairs = set()
    delb_arc_data = []
    for _, r in chain_rows.iterrows():
        bp_num = int(r['Breakpoint number'])
        delb_num = parse_deletion_bridge(r['Deletion bridge partner breakpoint'])
        if delb_num is None or delb_num not in bp_to_node:
            continue
        pair_key = tuple(sorted([bp_num, delb_num]))
        if pair_key in seen_delb_pairs:
            continue
        seen_delb_pairs.add(pair_key)
        node_a = bp_to_node[bp_num]
        node_b = bp_to_node[delb_num]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya = chr_y[node_a[0]]
        delb_arc_data.append((xa, ya, xb, ya))

    rads = _assign_rads(rearr_arc_data, [a[4] for a in rearr_arc_data])

    # SEGMENT highlight bars, grouped by chromosome row (for margin sizing)
    segment_draw_data = []  # (chrom, x_lo, x_hi, color, label)
    segment_table_rows = []  # (label, span_text) for the new table
    for dbp in db_pairs:
        db_idx = dbp['db_index']
        color = _db_color(db_idx)
        prefix = db_idx if use_db_prefix else None
        for seg in dbp['segments']:
            label = segment_display_label(seg, prefix)
            segment_table_rows.append((label, segment_span_text(seg)))
            nmap = nmaps[dbp['chromosome']]
            x_lo = nmap[seg['low_pos']]
            x_hi = nmap[seg['high_pos']]
            segment_draw_data.append((dbp['chromosome'], x_lo, x_hi, color, label))

    b_i = {c: len(chr_to_pos[c]) for c in sorted_chroms}
    theta = theta_notation(b_i)
    n_breakpoints = len(chain_rows)
    site_col = chain_rows['Site annotation'] if 'Site annotation' in chain_rows.columns else pd.Series([], dtype=str)
    contains_erg = bool(site_col.dropna().str.contains(ERG_PATTERN).any())
    erg_tag = '  [ERG/TMPRSS2 chain]' if contains_erg else ''

    k = len(sorted_chroms)
    X_AXES_RANGE = 1.12
    y_span_layout = max(y_span, CHR_SEP)
    fig_w = 24.0

    fig_h0 = max(7.0, k * 3.8 + 3.0)
    y_range0 = y_span_layout * 1.05 + 4.0
    scale_y0 = fig_h0 / y_range0
    scale_x = fig_w / X_AXES_RANGE

    same_row_arcs = (
        [(xa, ya, xb, yb, rad) for (xa, ya, xb, yb, is_inter), rad in zip(rearr_arc_data, rads) if not is_inter]
        + [(xa, ya, xb, yb, -0.12) for xa, ya, xb, yb in ref_arc_data]
        + [(xa, ya, xb, yb, 0.18) for xa, ya, xb, yb in delb_arc_data]
    )
    top_y = max(chr_y.values())
    bottom_y = min(chr_y.values())

    def _worst_sagitta_at(row_y):
        return max(
            (abs(rad) * abs(xb - xa) * scale_x
             for xa, ya, xb, yb, rad in same_row_arcs if ya == row_y),
            default=0.0,
        )

    # Segment highlight row sits BELOW its chromosome row at a fixed data
    # offset (SEG_LABEL_OFFSET) -- convert that fixed data-space offset
    # into the same physical-inches accounting the arc sagitta uses, so
    # it's included consistently in whichever margin (top/bottom row)
    # actually needs the extra room.
    chroms_with_segments = {c for c, xlo, xhi, col, lb in segment_draw_data}
    top_chrom = next(ch for ch in sorted_chroms if chr_y[ch] == top_y)
    bottom_chrom = next(ch for ch in sorted_chroms if chr_y[ch] == bottom_y)
    has_segments_at = {top_y: top_chrom in chroms_with_segments,
                        bottom_y: bottom_chrom in chroms_with_segments}
    seg_row_extra_in_top = (SEG_LABEL_OFFSET * scale_y0 * 1.15) if has_segments_at.get(top_y) else 0.0
    seg_row_extra_in_bottom = (SEG_LABEL_OFFSET * scale_y0 * 1.15) if has_segments_at.get(bottom_y) else 0.0

    extra_in_top = max(_worst_sagitta_at(top_y) * 1.15, seg_row_extra_in_top)
    extra_in_bottom = max(_worst_sagitta_at(bottom_y) * 1.15, seg_row_extra_in_bottom)
    extra_data_top = extra_in_top / scale_y0
    extra_data_bottom = extra_in_bottom / scale_y0

    fig_h = fig_h0 + extra_in_top + extra_in_bottom

    table_w, table_h, n_col_groups, table_max_rows_per_col, table_fontsize = _bp_table_dims(n_breakpoints, fig_h)
    seg_table_w, seg_table_h, seg_n_col_groups, seg_max_rows_per_col, seg_table_fontsize = \
        _bp_table_dims(len(segment_table_rows), fig_h, max_col_groups=2)

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_facecolor('white')
    x_margin = (X_AXES_RANGE - 1.0) / 2
    ax.set_xlim(-x_margin, 1.0 + x_margin)
    ax.set_ylim(bottom_y - 2.0 - extra_data_bottom - SEG_LABEL_OFFSET, top_y + 2.0 + extra_data_top)
    ax.axis('off')

    # Full per-DB-pair listing overflows the figure width once a chain has
    # more than a handful of DB pairs (e.g. P09-1042 chain 1 has 19) —
    # fall back to a plain count for those; the Segment Span table (not
    # the title) is the definitive per-DB-pair listing regardless.
    if len(db_pairs) <= 8:
        db_summary = ', '.join(f"DB{d['db_index']}(BP{d['low_bp']}↔BP{d['high_bp']}, "
                                f"{len(d['segments'])} segs)" for d in db_pairs)
    else:
        total_segs = sum(len(d['segments']) for d in db_pairs)
        db_summary = (f"{len(db_pairs)} deletion bridge pairs, {total_segs} total segments "
                       f"— see Segment Span Information table for the full per-pair listing")
    ax.set_title(
        f'{patient_id} — Baca Chain {chain_number}{erg_tag}   [{len(db_pairs)} deletion bridge(s)]\n'
        f'{theta}   n_breakpoints={n_breakpoints}   n_rearrangement_edges={len(rearr_arc_data)}'
        f'   n_adjacency_edges={len(ref_arc_data)}   n_deletion_bridge_edges={len(delb_arc_data)}\n'
        f'{db_summary}',
        fontsize=11, fontweight='bold', pad=14, loc='center',
    )

    for c in sorted_chroms:
        draw_chr_track(ax, c, chr_y[c], chr_to_pos[c], nmaps[c])
        draw_bp_number_labels(ax, c, chr_y[c], chr_to_rows[c], nmaps[c])

    for xa, ya, xb, yb in ref_arc_data:
        _arc(ax, xa, ya, xb, yb, -0.12, COL['cycle_ref'], 1.4, '-', 0.7, '-', 3)

    for xa, ya, xb, yb in delb_arc_data:
        _arc(ax, xa, ya, xb, yb, 0.18, DELETION_BRIDGE_COLOR, 1.6, ':', 0.85, '-', 4)

    for (xa, ya, xb, yb, is_inter), rad in zip(rearr_arc_data, rads):
        col = COL['rearr_inter'] if is_inter else COL['rearr_intra']
        _arc(ax, xa, ya, xb, yb, rad, col, 2.0, '--', 0.9, '-', 5)

    # SEGMENT highlights — thin colored underline + letter label, one
    # shared row per chromosome (not stacked per DB pair, per the design
    # note in the module docstring: bounded margin cost, table is the
    # source of truth if colors visually overlap on dense chains).
    for chrom, x_lo, x_hi, color, label in segment_draw_data:
        y = chr_y[chrom] - SEG_ROW_OFFSET
        ax.plot([x_lo, x_hi], [y, y], color=color, linewidth=4.5,
                solid_capstyle='round', alpha=0.85, zorder=6)
        mid = (x_lo + x_hi) / 2
        ax.text(mid, chr_y[chrom] - SEG_LABEL_OFFSET, label, ha='center', va='top',
                fontsize=8, fontweight='bold', color=color, zorder=7,
                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec=color, lw=0.9, alpha=0.95))

    fig.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Rearrangement edge — intrachromosomal'),
        mpatches.Patch(color=COL['rearr_inter'], label='Rearrangement edge — interchromosomal'),
        mpatches.Patch(color=COL['cycle_ref'],   label='Adjacency reference edge'),
        mpatches.Patch(color=DELETION_BRIDGE_COLOR, label='Deletion bridge edge'),
        mpatches.Patch(color=DB_PALETTE[0], label='DB-span segment (color = DB pair index)'),
    ], loc='lower center', ncol=5, fontsize=9.5, framealpha=0.95, edgecolor='#CCCCCC',
       bbox_to_anchor=(0.5, 0.0))

    right_margin = max(0.96 - table_w - seg_table_w, 0.55)
    fig.subplots_adjust(left=0.04, right=right_margin, top=0.88, bottom=0.08)
    _add_bp_location_table(fig, chain_rows, table_w, table_h, n_col_groups, table_max_rows_per_col, table_fontsize)
    _add_segment_span_table(fig, segment_table_rows, table_w, seg_table_w, seg_table_h,
                             seg_n_col_groups, seg_max_rows_per_col, seg_table_fontsize)

    out_path = os.path.join(save_dir, f'{patient_id}_chain_{chain_number}.png')
    fig.savefig(out_path, dpi=150, facecolor='white')
    plt.close(fig)
    return out_path


_TABLE_Y_TOP = 0.87


def _add_segment_span_table(fig, segment_table_rows, bp_table_w, seg_table_w, seg_table_h,
                             n_col_groups, max_rows_per_col, fontsize):
    """Second corner table, placed to the LEFT of the existing BP-location
    table (which occupies the rightmost table_w of the figure)."""
    n = len(segment_table_rows)
    left_edge = 1.0 - bp_table_w - seg_table_w
    ax_table = fig.add_axes([left_edge + 0.01, max(_TABLE_Y_TOP - seg_table_h, 0.02), seg_table_w - 0.015, seg_table_h])
    ax_table.axis('off')
    ax_table.set_title(f'Segment Span Information  (n={n})', fontsize=9.5, fontweight='bold',
                        loc='left', pad=3)

    col_width = 1.0 / n_col_groups
    for gi in range(n_col_groups):
        chunk = segment_table_rows[gi * max_rows_per_col: (gi + 1) * max_rows_per_col]
        lines = [f'{label:<8} {span}' for label, span in chunk]
        ax_table.text(
            gi * col_width, 1.0, '\n'.join(lines),
            transform=ax_table.transAxes, ha='left', va='top',
            fontsize=fontsize, family='monospace', linespacing=1.6,
        )


def run_all(save_dir=OUTPUT_DIR):
    os.makedirs(save_dir, exist_ok=True)
    df = load_mmc5()

    n_drawn, n_skipped_no_db, n_failed = 0, 0, 0
    for (patient_id, chain_number), grp in df.groupby(['Individual', 'Chain number']):
        if not chain_has_db(grp):
            n_skipped_no_db += 1
            continue
        try:
            out = draw_baca_chain_with_segments(patient_id, chain_number, grp, save_dir)
            n_drawn += 1
            print(f'{patient_id} chain {chain_number} -> {out}')
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    print(f'\nDrawn: {n_drawn}  Skipped (no DB): {n_skipped_no_db}  Failed: {n_failed}  -> {save_dir}')


if __name__ == '__main__':
    run_all()
