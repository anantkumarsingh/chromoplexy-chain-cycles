"""
baca_proposed_chains_drawing.py — visualize BACA'S OWN chains from mmc5.xlsx.

This is the deferred visualization phase for the work in
baca_chain_comparison.py / BACA_CHAIN_VS_OUR_CHAIN_COMPARISON.md: one image
per (patient, Chain number) pair in mmc5's Table S5A — Baca's pre-computed,
statistically-filtered chromoplexy chain groupings, NOT our own ERG-anchored
or full-genome reconstructions.

Each image shows every chromosome touched by that one chain, schematically
spaced (equal spacing, not genomic coordinates — same rationale as
baca_aberration_cycle_drawing.py: real DSB positions often cluster within a
few hundred bp of each other and genomic spacing would collapse them onto
the same pixel).

mmc5 actually encodes TWO distinct edge types, verified directly (not
assumed) before drawing:
  - 'Rearrangement number': two breakpoint rows sharing the same number are
    the two ends of one physical rearrangement — confirmed 1549 same-chrom +
    671 cross-chrom pairs across the dataset, matching the inter/intra split
    we'd expect from real structural variants. This is the REARRANGEMENT
    edge (drawn as orange/blue dashed arcs, same convention as
    baca_aberration_cycle_drawing.py).
  - 'Adjacent breakpoint(s)': verified 100% same-chromosome (2654/2654 links,
    zero cross-chromosome) — this is NOT a rearrangement partner, it's the
    genomically-NEXT breakpoint on the same chromosome. This is exactly a
    REFERENCE edge in our own AMG terminology, given directly by Baca rather
    than inferred from position order. Drawn as thin green links.
  (An earlier version of this script wrongly treated 'Adjacent breakpoint(s)'
  as the rearrangement edge — caught by checking same-vs-cross-chromosome
  counts per patient before trusting it.)

Because mmc5 gives BOTH edge types directly, Baca's own chain data could in
principle support real AMG cycle-closure computation without inferring
reference edges from genomic position order the way our other scripts do —
that's a follow-up worth revisiting, not implemented here. This script is
purely a visualization of Baca's chain exactly as he defined it.

Output: results/baca_proposed_chains/{patient_id}_cycle_{chain_number}.png
"""

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
    draw_bp_labels,
    _arc,
    _assign_rads,
    _ax_defaults,
    COL,
    CHR_SEP,
)

BACA_DATASET_FOLDER = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/baca_dataset"
MMC5_FILE_PATH = os.path.join(BACA_DATASET_FOLDER, "mmc5.xlsx")

OUTPUT_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results/baca_proposed_chains"

STRAND_MAP = {'Forward': '+', 'Reverse': '-'}
ERG_PATTERN = re.compile(r'ERG|TMPRSS2')


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


def theta_notation(b_i):
    chroms = sorted(b_i.keys())
    k = len(chroms)
    b_str = ",".join(str(b_i[c]) for c in chroms)
    return f"Theta({k},({b_str}))"


def draw_baca_chain(patient_id, chain_number, chain_rows, save_dir):
    bp_to_node = {
        int(r['Breakpoint number']): (int(r['chromosome']), r['position'], r['node_strand'])
        for _, r in chain_rows.iterrows()
    }
    bp_to_rearr_num = dict(zip(chain_rows['Breakpoint number'].astype(int), chain_rows['Rearrangement number']))

    chr_to_pos = {
        int(chrom): sorted(grp['position'].unique())
        for chrom, grp in chain_rows.groupby('chromosome')
    }
    sorted_chroms = sorted(chr_to_pos.keys())
    nmaps = {c: build_norm_map(chr_to_pos[c]) for c in sorted_chroms}
    chr_y = {c: (len(sorted_chroms) - 1 - i) * CHR_SEP for i, c in enumerate(sorted_chroms)}
    y_span = (len(sorted_chroms) - 1) * CHR_SEP

    # REARRANGEMENT edges — two breakpoints sharing the same Rearrangement
    # number are the two ends of one real structural variant (verified
    # inter/intra split against the raw file's Class column elsewhere).
    nodes_with_rearr = set()
    rearr_arc_data = []
    n_unresolved_rearr = 0
    for rearr_num, grp in chain_rows.groupby('Rearrangement number'):
        bps = grp['Breakpoint number'].astype(int).tolist()
        if len(bps) != 2:
            n_unresolved_rearr += len(bps)  # this rearrangement's other end isn't in this chain
            continue
        node_a, node_b = bp_to_node[bps[0]], bp_to_node[bps[1]]
        nodes_with_rearr.add(node_a)
        nodes_with_rearr.add(node_b)
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya, yb = chr_y[node_a[0]], chr_y[node_b[0]]
        rearr_arc_data.append((xa, ya, xb, yb, node_a[0] != node_b[0]))

    # REFERENCE edges — 'Adjacent breakpoint(s)' is always same-chromosome
    # (verified 2654/2654 across the dataset): the next breakpoint along the
    # same chromosome, given directly by Baca rather than inferred from
    # position order.
    seen_ref_pairs = set()
    ref_arc_data = []
    for _, r in chain_rows.iterrows():
        bp_num = int(r['Breakpoint number'])
        node_a = bp_to_node[bp_num]
        for adj_num in parse_adjacent(r['Adjacent breakpoint(s)']):
            if adj_num not in bp_to_node:
                continue  # neighbor breakpoint not part of this chain in mmc5
            pair_key = tuple(sorted([bp_num, adj_num]))
            if pair_key in seen_ref_pairs:
                continue
            seen_ref_pairs.add(pair_key)
            node_b = bp_to_node[adj_num]
            xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
            xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
            ya = chr_y[node_a[0]]
            ref_arc_data.append((xa, ya, xb, ya))

    rads = _assign_rads(rearr_arc_data, [a[4] for a in rearr_arc_data])

    b_i = {c: len(chr_to_pos[c]) for c in sorted_chroms}
    theta = theta_notation(b_i)
    n_breakpoints = len(chain_rows)
    site_col = chain_rows['Site annotation'] if 'Site annotation' in chain_rows.columns else pd.Series([], dtype=str)
    contains_erg = bool(site_col.dropna().str.contains(ERG_PATTERN).any())
    erg_tag = '  [ERG/TMPRSS2 chain]' if contains_erg else ''

    # Rearrangement-edge curvature (matplotlib's arc3) computes its sagitta
    # in DISPLAY/PHYSICAL space: sagitta_inches = |rad| * dx_data * scale_x,
    # where scale_x = fig_w / xlim_range is FIXED once fig_w is chosen — i.e.
    # the sagitta in inches does not depend on fig_h or the y-axis data
    # range at all. A wide-spanning intrachromosomal arc with deep curvature
    # (rad up to -0.78) on a k=1-2 chain (small y_span, so little baseline
    # bottom padding) needs that many physical inches of room below the
    # chromosome row, or it clips.
    #
    # Fix: extend fig_h and the bottom y-limit by the exact same amount in
    # a matched inches<->data-unit ratio (scale_y0, computed from the
    # baseline sizing below) so the look of the baseline (well-tuned for
    # multi-chromosome chains) is preserved exactly, with the extra room
    # added only at the bottom where the arc actually needs it — instead of
    # inflating fig_h independently of the data range, which leaves most of
    # the figure as wasted blank space.
    k = len(sorted_chroms)
    X_AXES_RANGE = 1.12  # matches _ax_defaults' xlim(-0.02, 1.10)
    y_span_layout = max(y_span, CHR_SEP)
    fig_w = 24.0

    fig_h0 = max(7.0, k * 3.8 + 3.0)
    y_range0 = y_span_layout * 1.05 + 4.0
    scale_y0 = fig_h0 / y_range0  # inches per data-unit, baseline (unclipped) sizing
    scale_x = fig_w / X_AXES_RANGE

    worst_sagitta_in = max(
        (abs(rad) * abs(xb - xa) * scale_x
         for (xa, ya, xb, yb, is_inter), rad in zip(rearr_arc_data, rads)
         if not is_inter),
        default=0.0,
    )
    extra_in = worst_sagitta_in * 1.15  # 15% safety margin
    extra_data = extra_in / scale_y0

    fig_h = fig_h0 + extra_in
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_facecolor('white')
    ax.set_xlim(-0.02, 1.10)
    ax.set_ylim(-1.8 - y_span_layout * 0.05 - extra_data, y_span_layout + 2.2)
    ax.axis('off')

    ax.set_title(
        f'{patient_id} — Baca Chain {chain_number}{erg_tag}\n'
        f'{theta}   n_breakpoints={n_breakpoints}   n_rearrangement_edges={len(rearr_arc_data)}'
        f'   n_reference_edges={len(ref_arc_data)}   rearrangement_ends_outside_chain={n_unresolved_rearr}',
        fontsize=12, fontweight='bold', pad=14,
    )

    for c in sorted_chroms:
        draw_chr_track(ax, c, chr_y[c], chr_to_pos[c], nmaps[c])
        draw_bp_labels(ax, c, chr_y[c], chr_to_pos[c], nmaps[c], nodes_with_rearr)

    for xa, ya, xb, yb in ref_arc_data:
        _arc(ax, xa, ya, xb, yb, -0.12, COL['cycle_ref'], 1.4, '-', 0.7, '-', 3)

    for (xa, ya, xb, yb, is_inter), rad in zip(rearr_arc_data, rads):
        col = COL['rearr_inter'] if is_inter else COL['rearr_intra']
        _arc(ax, xa, ya, xb, yb, rad, col, 2.0, '--', 0.9, '-', 5)

    ax.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Rearrangement edge — intrachromosomal'),
        mpatches.Patch(color=COL['rearr_inter'], label='Rearrangement edge — interchromosomal'),
        mpatches.Patch(color=COL['cycle_ref'],   label='Reference edge (Adjacent breakpoint)'),
    ], loc='lower right', fontsize=9, framealpha=0.95, edgecolor='#CCCCCC')

    plt.tight_layout(pad=2.0)
    out_path = os.path.join(save_dir, f'{patient_id}_cycle_{chain_number}.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
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
