"""
baca_chain_visualization.py — ONE BIG IMAGE PER BACA CHAIN (mmc5), showing
the full chain structure: every real edge Baca's data gives us. This is
the chain-level companion to baca_cycle_visualization.py (which highlights
only the CLOSED CYCLES within chains that actually have one).

Supersedes baca_proposed_chains_drawing.py for the "gold mine" primary
substrate (see CLAUDE.md "FINAL DECISION: primary substrate is Baca's own
mmc5 chains" / memory cycle_completion_phase.md): that earlier script only
drew 2 of mmc5's 3 real edge types (rearrangement + adjacency). This one
adds the THIRD — deletion bridge — since that's exactly the edge type that
makes mmc5-chain closure deletion-bridge-aware (the whole reason it's now
the primary substrate over the full-genome model).

Three real edge types drawn, all from mmc5.xlsx Table S5A directly — no
raw chrom_aberrations_baca.csv data, no invented connections:
  - Rearrangement edge ('Rearrangement number'): orange (inter) / blue
    (intra) dashed arcs.
  - Adjacency reference edge ('Adjacent breakpoint(s)'): green links.
  - Deletion bridge edge ('Deletion bridge partner breakpoint'): purple
    dotted links — drawn distinctly so it's visually clear which
    reference-type connections came from an actual recorded deletion
    bridge vs. plain statistical adjacency.

Output: results/baca_proposed_chains/{patient_id}_chain_{chain_number}.png
        — one image per (Individual, Chain number) in mmc5, all 366
        chains, all 52 mmc5-covered patients.
"""

import math
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
from baca_aberration_cycle_drawing import (
    build_norm_map,
    node_x,
    draw_chr_track,
    _arc,
    _assign_rads,
    COL,
    CHR_SEP,
)

BACA_DATASET_FOLDER = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/baca_dataset"
MMC5_FILE_PATH = os.path.join(BACA_DATASET_FOLDER, "mmc5.xlsx")

OUTPUT_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results/baca_proposed_chains"

STRAND_MAP = {'Forward': '+', 'Reverse': '-'}
ERG_PATTERN = re.compile(r'ERG|TMPRSS2')
DELETION_BRIDGE_COLOR = '#6A1B9A'  # purple — distinct from green adjacency edges


def load_mmc5():
    df = pd.ExcelFile(MMC5_FILE_PATH).parse('Table S5A')
    df.columns = [c.strip() for c in df.columns]
    df = df[df['Chromosome:position'].notna()].copy()
    df['chromosome'] = df['Chromosome:position'].str.split(':').str[0].str.replace('chr', '', regex=False).astype(int)
    df['position'] = df['Chromosome:position'].str.split(':').str[1].astype(float)
    df['node_strand'] = df['Strand'].map(STRAND_MAP)
    return df


def parse_adjacent(s):
    if pd.isna(s):
        return []
    return [int(x) for x in re.findall(r'\d+', str(s))]


def parse_deletion_bridge(v):
    if pd.isna(v):
        return None
    s = str(v).strip()
    return int(float(s)) if s else None


def theta_notation(b_i):
    chroms = sorted(b_i.keys())
    k = len(chroms)
    b_str = ",".join(str(b_i[c]) for c in chroms)
    return f"Theta({k},({b_str}))"


def draw_bp_number_labels(ax, chrom, y, chrom_rows, nmap, max_positions=40):
    """
    Label each breakpoint node with its mmc5 'Breakpoint number' + strand
    (e.g. 'BP14 (-)') directly on the diagram, staggered above/below in
    position order to avoid overlap. Genomic location is deliberately NOT
    drawn here -- it's cluttered the diagram in the past and is exactly
    what joins what is hard to read; full location lives in the
    corner BP#-Location table instead (see _add_bp_location_table), so the
    on-diagram label only ever needs to carry the two things that identify
    a node in the edge lists: BP number and strand.
    """
    rows_sorted = chrom_rows.sort_values('position')
    if len(rows_sorted) > max_positions:
        return
    for i, (_, r) in enumerate(rows_sorted.iterrows()):
        pos, strand, bp_num = r['position'], r['node_strand'], int(r['Breakpoint number'])
        x = node_x(pos, strand, nmap)
        above = (i % 2 == 0)
        ly = y + 0.46 if above else y - 0.46
        va = 'bottom' if above else 'top'
        ax.text(
            x, ly, f'BP{bp_num} ({strand})', ha='center', va=va,
            fontsize=8.5, fontweight='bold', color='#1A1A1A', zorder=6,
            bbox=dict(boxstyle='round,pad=0.18', fc='white', ec='#999999', lw=0.7, alpha=0.9),
        )


_TABLE_Y_TOP, _TABLE_Y_BOTTOM = 0.87, 0.05  # figure-fraction vertical budget below the title


def _bp_table_dims(n, fig_h, max_col_groups=3):
    """
    Width/height (figure-fraction) the corner table needs for n entries.
    These figures are usually very TALL (large chains stack many
    chromosome rows and can reach 30+ inches), so there's normally far
    more vertical room than horizontal — prefer growing the table
    DOWNWARD (one tall column) over sideways, and only wrap into extra
    columns when even a single column genuinely can't fit in the
    available height at a readable font size. This avoids the table
    eating a large fraction of figure width on dense chains (up to 118
    breakpoints), which would otherwise squeeze the diagram itself.
    """
    fontsize = 8.0 if n <= 40 else (7.0 if n <= 80 else (6.0 if n <= 150 else 5.5))
    line_height_in = fontsize * 1.6 / 72.0
    title_in = 0.28

    height_budget_in = (_TABLE_Y_TOP - _TABLE_Y_BOTTOM) * fig_h
    rows_that_fit = max(6, int((height_budget_in - title_in) / line_height_in))

    n_col_groups = max(1, min(max_col_groups, math.ceil(n / rows_that_fit)))
    max_rows_per_col = math.ceil(n / n_col_groups)

    table_w = min(0.115 * n_col_groups + 0.045, 0.30)
    table_h = min((max_rows_per_col * line_height_in + title_in) / fig_h, _TABLE_Y_TOP - _TABLE_Y_BOTTOM)
    return table_w, table_h, n_col_groups, max_rows_per_col, fontsize


def _add_bp_location_table(fig, chain_rows, table_w, table_h, n_col_groups, max_rows_per_col, fontsize):
    """
    Corner reference table: BP number -> full genomic location + strand,
    sorted by BP number (matches how the on-diagram labels reference them).
    Placed in FIGURE-fraction coordinates (not data coordinates), in the
    space the caller already reserved (via _bp_table_dims + shrinking the
    main axes) so it never overlaps the actual chain diagram.
    """
    rows = chain_rows[['Breakpoint number', 'chromosome', 'position', 'node_strand']].copy()
    rows['Breakpoint number'] = rows['Breakpoint number'].astype(int)
    rows = rows.sort_values('Breakpoint number')
    entries = [
        (int(r['Breakpoint number']), int(r['chromosome']), r['position'], r['node_strand'])
        for _, r in rows.iterrows()
    ]
    n = len(entries)

    # top=0.90 is reserved for the title by the main axes' subplots_adjust;
    # start just below it so the table never collides with the (often long,
    # full-width) title text.
    ax_table = fig.add_axes(
        [1.0 - table_w + 0.01, max(_TABLE_Y_TOP - table_h, 0.02), table_w - 0.015, table_h]
    )
    ax_table.axis('off')
    ax_table.set_title(f'BP # → Location  (n={n})', fontsize=9.5, fontweight='bold',
                        loc='left', pad=3)

    col_width = 1.0 / n_col_groups
    for gi in range(n_col_groups):
        chunk = entries[gi * max_rows_per_col: (gi + 1) * max_rows_per_col]
        lines = [f'BP{num:<4} chr{chrom}:{pos:,.0f} ({strand})' for num, chrom, pos, strand in chunk]
        ax_table.text(
            gi * col_width, 1.0, '\n'.join(lines),
            transform=ax_table.transAxes, ha='left', va='top',
            fontsize=fontsize, family='monospace', linespacing=1.6,
        )


def draw_baca_chain(patient_id, chain_number, chain_rows, save_dir):
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

    # ADJACENCY reference edges ('Adjacent breakpoint(s)')
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

    # DELETION BRIDGE edges ('Deletion bridge partner breakpoint') — the
    # third real edge type, not drawn by the earlier baca_proposed_chains_drawing.py.
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

    b_i = {c: len(chr_to_pos[c]) for c in sorted_chroms}
    theta = theta_notation(b_i)
    n_breakpoints = len(chain_rows)
    site_col = chain_rows['Site annotation'] if 'Site annotation' in chain_rows.columns else pd.Series([], dtype=str)
    contains_erg = bool(site_col.dropna().str.contains(ERG_PATTERN).any())
    erg_tag = '  [ERG/TMPRSS2 chain]' if contains_erg else ''

    # Sizing — same physical-inches sagitta logic as baca_proposed_chains_drawing.py
    # (matplotlib's arc3 curvature is computed in display inches, independent
    # of the y-axis data range; a wide, deeply-curved intrachromosomal arc on
    # a small-k chain needs explicit extra room below the chromosome row or
    # it clips against the figure bottom).
    k = len(sorted_chroms)
    X_AXES_RANGE = 1.12
    y_span_layout = max(y_span, CHR_SEP)
    fig_w = 24.0

    fig_h0 = max(7.0, k * 3.8 + 3.0)
    y_range0 = y_span_layout * 1.05 + 4.0
    scale_y0 = fig_h0 / y_range0
    scale_x = fig_w / X_AXES_RANGE

    # Worst-case arc sagitta across ALL same-row edge types (rearrangement
    # intra, adjacency, AND deletion bridge — the original script only
    # checked rearrangement arcs, missing that adjacency/deletion-bridge
    # arcs can bow the OPPOSITE direction (positive vs negative rad), which
    # was clipping/pushing content off-center at the top.
    #
    # Only arcs on the TOPMOST or BOTTOMMOST chromosome row can actually
    # clip past the figure's top/bottom edge — arcs on interior rows bow
    # into the existing CHR_SEP gap between neighboring rows and don't need
    # extra figure-level room. Using the single global worst-case sagitta
    # for BOTH margins (an earlier version of this fix) made dense,
    # many-breakpoint chains (e.g. 118 breakpoints, deep wide intra arcs)
    # balloon to ~56in tall with mostly blank space — caught by visually
    # inspecting PR-07-4814 chain 3 before trusting this fix either.
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

    extra_in_top = _worst_sagitta_at(top_y) * 1.15
    extra_in_bottom = _worst_sagitta_at(bottom_y) * 1.15
    extra_data_top = extra_in_top / scale_y0
    extra_data_bottom = extra_in_bottom / scale_y0

    fig_h = fig_h0 + extra_in_top + extra_in_bottom

    # Corner BP#-Location table dimensions, computed now that fig_h is
    # final, so the main diagram's right margin can be shrunk to make
    # room for it instead of the table overlapping the chain arcs.
    table_w, table_h, n_col_groups, table_max_rows_per_col, table_fontsize = _bp_table_dims(n_breakpoints, fig_h)

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_facecolor('white')
    # Symmetric margins on both axes — content (chromosome tracks span
    # roughly x in [0,1], y in [0, y_span_layout]) is centered in the frame
    # rather than pushed to one side by an asymmetric legend/arc margin.
    x_margin = (X_AXES_RANGE - 1.0) / 2
    ax.set_xlim(-x_margin, 1.0 + x_margin)
    ax.set_ylim(bottom_y - 2.0 - extra_data_bottom, top_y + 2.0 + extra_data_top)
    ax.axis('off')

    ax.set_title(
        f'{patient_id} — Baca Chain {chain_number}{erg_tag}\n'
        f'{theta}   n_breakpoints={n_breakpoints}   n_rearrangement_edges={len(rearr_arc_data)}'
        f'   n_adjacency_edges={len(ref_arc_data)}   n_deletion_bridge_edges={len(delb_arc_data)}'
        f'   rearrangement_ends_outside_chain={n_unresolved_rearr}',
        fontsize=12, fontweight='bold', pad=14, loc='center',
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

    # Legend placed OUTSIDE the data area (below the axes, centered) so it
    # never forces an asymmetric data margin the way loc='lower right'
    # (inside the axes) did — that was the other source of off-center content.
    fig.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Rearrangement edge — intrachromosomal'),
        mpatches.Patch(color=COL['rearr_inter'], label='Rearrangement edge — interchromosomal'),
        mpatches.Patch(color=COL['cycle_ref'],   label='Adjacency reference edge'),
        mpatches.Patch(color=DELETION_BRIDGE_COLOR, label='Deletion bridge edge'),
    ], loc='lower center', ncol=4, fontsize=10, framealpha=0.95, edgecolor='#CCCCCC',
       bbox_to_anchor=(0.5, 0.0))

    # Right margin shrunk by the table's reserved width so the corner
    # BP#-Location table sits in its own space, never overlapping the arcs.
    fig.subplots_adjust(left=0.04, right=max(0.96 - table_w, 0.6), top=0.90, bottom=0.08)
    _add_bp_location_table(fig, chain_rows, table_w, table_h, n_col_groups, table_max_rows_per_col, table_fontsize)

    out_path = os.path.join(save_dir, f'{patient_id}_chain_{chain_number}.png')
    fig.savefig(out_path, dpi=150, facecolor='white')
    plt.close(fig)
    return out_path


def run_all(save_dir=OUTPUT_DIR):
    os.makedirs(save_dir, exist_ok=True)
    df = load_mmc5()

    n_drawn, n_failed = 0, 0
    for (patient_id, chain_number), grp in df.groupby(['Individual', 'Chain number']):
        try:
            out = draw_baca_chain(patient_id, chain_number, grp, save_dir)
            n_drawn += 1
            print(f'{patient_id} chain {chain_number} -> {out}')
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    print(f'\nDrawn: {n_drawn}  Failed: {n_failed}  -> {save_dir}')


if __name__ == '__main__':
    run_all()
