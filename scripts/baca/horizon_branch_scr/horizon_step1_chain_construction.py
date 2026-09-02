"""
horizon_step1_chain_construction.py — HORIZON branch, Step 1.
HORIZON branch started 2026-08-24.

Rebuilds all 366 Baca mmc5 chains (SAME 'Chain number' breakpoint
membership as the existing mmc5-primary substrate — see CLAUDE.md
"FINAL DECISION: primary substrate is Baca's own mmc5 chains") using ONLY:

  - Rearrangement edges ('Rearrangement number' pairs) — real, physical
    fusion junctions, exactly as before.
  - Our OWN genomic-position-order-inferred reference edges, via
    _build_amg_nodes_and_edges (the same function used for the ERG-chain /
    full-genome model BEFORE Phase 1.5 merged in mmc5's adjacency +
    deletion-bridge edges — i.e. "the original way").

This branch deliberately NEVER reads or uses mmc5's 'Adjacent
breakpoint(s)' or 'Deletion bridge partner breakpoint' columns, anywhere,
per explicit user instruction: those edges may exist in the chain, they
are simply never drawn or used to traverse/close cycles in this branch.

OPEN-END DEFINITION — LOCKED 2026-08-24 (confirmed directly with user,
also documented in CLAUDE.md's HORIZON section — treat that as the
canonical write-up, this is a pointer to it):
  Every breakpoint gap has two ends, a '-' (left) end and a '+' (right)
  end. The end Baca's data actually recorded always has a rearrangement
  edge, AND every end always has a reference edge too (either to a real
  neighboring breakpoint, or a "telomere edge with respect to this chain"
  if nothing further exists in this chain's own data — NOT a claim about
  the true chromosome end, since only this chain's subset of the genome
  is visible here). So the recorded end is NEVER open.
  The breakpoint's OTHER end — same position, opposite strand, never
  independently recorded — has no rearrangement edge, so it IS open,
  unless that exact opposite end also happens to be independently
  recorded as its own separate real breakpoint elsewhere in this chain.
  "Open" therefore means "has no rearrangement edge," full stop — nothing
  to do with whether a reference edge resolves to a numbered breakpoint,
  which was an earlier, incorrect definition used in this script before
  this date (see the old CONTEXT_LOG entries / conversation history for
  that superseded attempt — 'MISSING_REARR'/'dangling' terminology no
  longer applies anywhere in this branch).

Per explicit user instruction, the CHAIN DIAGRAM ITSELF IS UNCHANGED by
this definition — the gray chromosome-track bars already represent the
reference/telomere structure cleanly; no new reference-edge arrows are
drawn. Only the CSV outputs below reflect the corrected open-end concept.

Outputs:
  results/horizon_branch_res/horizon_csv_data/horizon_chain_structure_summary.csv
    — 1 row per chain (366 rows), including n_open_ends (count) and
    open_ends_bp_strand_list (exact "BP{num}({strand})" list)
  results/horizon_branch_res/horizon_csv_data/horizon_breakpoint_edges.csv
    — 1 row per recorded breakpoint (~4,440 rows), each also reporting its
    own OTHER end's strand and open/closed status
  results/horizon_branch_res/horizon_chain_cycle_visualizations/chains/
    {patient_id}_horizon_chain_{n}.png — 1 image per chain (366 images,
    UNCHANGED visually from before this definition fix)
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

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
from baca_aberration_clinical_analysis import (
    BACA_CLINICAL_PHENO_FILE_PATH,
    load_ets_status_from_mmc5,
    _build_amg_nodes_and_edges,
)
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

RESULTS_ROOT = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results/horizon_branch_res"
CSV_DIR = os.path.join(RESULTS_ROOT, "horizon_csv_data")
CHAIN_IMG_DIR = os.path.join(RESULTS_ROOT, "horizon_chain_cycle_visualizations", "chains")

CHAIN_SUMMARY_CSV = os.path.join(CSV_DIR, "horizon_chain_structure_summary.csv")
BREAKPOINT_EDGES_CSV = os.path.join(CSV_DIR, "horizon_breakpoint_edges.csv")

STRAND_MAP = {'Forward': '+', 'Reverse': '-'}
ERG_PATTERN = re.compile(r'ERG|TMPRSS2')
GENE_PATTERN = re.compile(r'([A-Z][A-Z0-9\-]+)\([+-]\)')

REF_EDGE_COLOR = COL['cycle_ref']  # green — same color the rest of this
                                    # codebase already uses for reference edges


# ── Load ─────────────────────────────────────────────────────────────────────

def load_mmc5_chains():
    df = pd.ExcelFile(MMC5_FILE_PATH).parse('Table S5A')
    df.columns = [c.strip() for c in df.columns]
    df = df[df['Chromosome:position'].notna()].copy()
    df['chromosome'] = df['Chromosome:position'].str.split(':').str[0].str.replace('chr', '', regex=False).astype(int)
    df['position'] = df['Chromosome:position'].str.split(':').str[1].astype(float)
    df['node_strand'] = df['Strand'].map(STRAND_MAP)
    return df


def extract_genes(site_strings):
    genes = set()
    for s in site_strings:
        if pd.isna(s):
            continue
        genes.update(GENE_PATTERN.findall(str(s)))
    return sorted(genes)


# ── Edge computation (rearrangement + our own inferred reference only) ───────

def build_rearr_rows(chain_rows):
    """
    Convert mmc5's one-row-per-breakpoint chain rows into the BP1/BP2
    one-row-per-rearrangement shape _build_amg_nodes_and_edges expects.
    Every 'Rearrangement number' group within a chain has exactly 2 rows
    — verified cohort-wide (0/366 chains have a rearrangement end whose
    partner falls outside its own chain) — so this is always well-defined.
    """
    out = []
    for _, g in chain_rows.groupby('Rearrangement number'):
        g = g.sort_values('Breakpoint number')
        r0, r1 = g.iloc[0], g.iloc[1]
        out.append({
            'Breakpoint 1 chromosome': r0['chromosome'],
            'Breakpoint 1 position': r0['position'],
            'Breakpoint 1 strand': r0['node_strand'],
            'Breakpoint 2 chromosome': r1['chromosome'],
            'Breakpoint 2 position': r1['position'],
            'Breakpoint 2 strand': r1['node_strand'],
        })
    return pd.DataFrame(out)


def compute_chain_edges(chain_rows):
    """
    Returns:
      bp_to_node                : Breakpoint number -> (chr,pos,strand),
                                   always 1:1 (every row has its own bp
                                   number regardless of node collisions).
      node_to_bps                : (chr,pos,strand) -> sorted list of
                                   Breakpoint numbers AT that exact node.
                                   Usually length 1. A small number of real
                                   mmc5 rows (verified: 2/366 chains, 1
                                   extra breakpoint each — e.g. PR-4240
                                   chain 2's BP40 and BP531) record TWO
                                   different breakpoint numbers, belonging
                                   to two different rearrangements, at the
                                   exact same (chromosome, position,
                                   strand) — a genuine multi-breakpoint
                                   coincidence in the data, not a parsing
                                   bug. Both breakpoints are topologically
                                   the same AMG node, so they always get
                                   the identical reference-edge outcome.
      rearr_adj, ref_adj         : as returned by _build_amg_nodes_and_edges
      bp_rearrangement_partner  : bp_number -> bp_number. Computed directly
                                   from chain_rows via 'Rearrangement
                                   number' grouping (NOT via rearr_adj/node
                                   identity) so it stays correct per-
                                   breakpoint even at a colliding node,
                                   where _build_amg_nodes_and_edges's own
                                   node-keyed rearr_adj can only keep one
                                   of two genuinely different rearrangement
                                   partners for that shared node (its own
                                   docstring already documents "last
                                   assignment wins" for this situation).
      bp_reference_partner      : bp_number -> bp_number(s, '/'-joined if
                                   the target node also collides) |
                                   'TELOMERE (this chain)' | a plain
                                   "chr:pos(strand)" string.
        LOCKED DEFINITION (confirmed with user 2026-08-24 — see CLAUDE.md
        HORIZON section): a recorded breakpoint's own end ALWAYS has a
        reference edge — the position-order neighbor computed here is
        never "missing," it just doesn't always land on another
        INDEPENDENTLY recorded breakpoint. When it doesn't, we report the
        actual neighboring coordinate/strand instead of a bp number, so
        the connection is still fully visible — but this NEVER makes the
        recorded end "open." Whether an end is open is a completely
        separate question, computed in run_all() from a different rule
        (does THIS bp's own OTHER, un-recorded side have a rearrangement
        edge or not) — see run_all()'s open-ends block for that rule.
        'TELOMERE (this chain)' means the reference edge points off the
        end of what's recorded FOR THIS CHAIN specifically — not a claim
        about the true chromosome end, since only this chain's own subset
        of breakpoints is visible here (not the full genome).
    """
    bp_to_node = {}
    node_to_bps = {}
    for _, r in chain_rows.iterrows():
        bp_num = int(r['Breakpoint number'])
        node = (int(r['chromosome']), r['position'], r['node_strand'])
        bp_to_node[bp_num] = node
        node_to_bps.setdefault(node, []).append(bp_num)
    for node in node_to_bps:
        node_to_bps[node].sort()

    rearr_rows = build_rearr_rows(chain_rows)
    all_nodes, rearr_adj, ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(rearr_rows)

    bp_rearrangement_partner = {}
    for _, g in chain_rows.groupby('Rearrangement number'):
        g = g.sort_values('Breakpoint number')
        bp_a, bp_b = int(g.iloc[0]['Breakpoint number']), int(g.iloc[1]['Breakpoint number'])
        bp_rearrangement_partner[bp_a] = bp_b
        bp_rearrangement_partner[bp_b] = bp_a

    bp_reference_partner = {}
    for bp_num, node in bp_to_node.items():
        ref_target = ref_adj.get(node)
        if ref_target == 'TELOMERE':
            bp_reference_partner[bp_num] = 'TELOMERE (this chain)'
        elif ref_target in node_to_bps:
            partners = node_to_bps[ref_target]
            bp_reference_partner[bp_num] = '/'.join(str(b) for b in partners) if len(partners) > 1 else partners[0]
        else:
            # A real genomic neighbor exists, it just isn't independently
            # recorded as its own breakpoint in this chain — show exactly
            # where it is rather than a bp number. This does NOT make the
            # recorded end (bp_num) open — see module note above.
            t_chrom, t_pos, t_strand = ref_target
            bp_reference_partner[bp_num] = f'chr{t_chrom}:{t_pos:,.0f}({t_strand}) [no separate BP# here]'

    return bp_to_node, node_to_bps, rearr_adj, ref_adj, bp_rearrangement_partner, bp_reference_partner


# ── Drawing ──────────────────────────────────────────────────────────────────

def draw_bp_number_labels(ax, chrom, y, chrom_rows, nmap, max_positions=40):
    """BP number + strand, staggered above/below in position order."""
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


_TABLE_Y_TOP, _TABLE_Y_BOTTOM = 0.87, 0.05


def _bp_table_dims(n, fig_h, max_col_groups=3):
    fontsize = 8.0 if n <= 40 else (7.0 if n <= 80 else (6.0 if n <= 150 else 5.5))
    line_height_in = fontsize * 1.6 / 72.0
    title_in = 0.28

    height_budget_in = (_TABLE_Y_TOP - _TABLE_Y_BOTTOM) * fig_h
    rows_that_fit = max(6, int((height_budget_in - title_in) / line_height_in))

    n_col_groups = max(1, min(max_col_groups, math.ceil(n / rows_that_fit)))
    max_rows_per_col = math.ceil(n / n_col_groups)

    table_w = min(0.145 * n_col_groups + 0.045, 0.34)
    table_h = min((max_rows_per_col * line_height_in + title_in) / fig_h, _TABLE_Y_TOP - _TABLE_Y_BOTTOM)
    return table_w, table_h, n_col_groups, max_rows_per_col, fontsize


def _add_bp_location_partner_table(fig, chain_rows, bp_rearrangement_partner,
                                    table_w, table_h, n_col_groups, max_rows_per_col, fontsize):
    """
    Corner reference table: BP number, strand, genomic location, and its
    Rearrangement Partner BP — the one real edge every breakpoint always
    has in this branch (reference-partner status is a per-chain dangling/
    resolved question, reported instead in horizon_breakpoint_edges.csv).
    """
    rows = chain_rows[['Breakpoint number', 'chromosome', 'position', 'node_strand']].copy()
    rows['Breakpoint number'] = rows['Breakpoint number'].astype(int)
    rows = rows.sort_values('Breakpoint number')
    entries = [
        (int(r['Breakpoint number']), int(r['chromosome']), r['position'], r['node_strand'])
        for _, r in rows.iterrows()
    ]
    n = len(entries)

    ax_table = fig.add_axes(
        [1.0 - table_w + 0.01, max(_TABLE_Y_TOP - table_h, 0.02), table_w - 0.015, table_h]
    )
    ax_table.axis('off')
    ax_table.set_title(f'BP # → Location, Strand, Partner  (n={n})', fontsize=9.5,
                        fontweight='bold', loc='left', pad=3)

    col_width = 1.0 / n_col_groups
    for gi in range(n_col_groups):
        chunk = entries[gi * max_rows_per_col: (gi + 1) * max_rows_per_col]
        lines = [
            f'BP{num:<4} chr{chrom}:{pos:,.0f} ({strand})  partner=BP{bp_rearrangement_partner[num]}'
            for num, chrom, pos, strand in chunk
        ]
        ax_table.text(
            gi * col_width, 1.0, '\n'.join(lines),
            transform=ax_table.transAxes, ha='left', va='top',
            fontsize=fontsize, family='monospace', linespacing=1.6,
        )


def draw_horizon_chain(patient_id, chain_number, chain_rows, bp_to_node, node_to_bps,
                        ref_adj, bp_rearrangement_partner, save_dir):
    chr_to_pos = {
        int(chrom): sorted(grp['position'].unique())
        for chrom, grp in chain_rows.groupby('chromosome')
    }
    chr_to_rows = {int(chrom): grp for chrom, grp in chain_rows.groupby('chromosome')}
    sorted_chroms = sorted(chr_to_pos.keys())
    nmaps = {c: build_norm_map(chr_to_pos[c]) for c in sorted_chroms}
    chr_y = {c: (len(sorted_chroms) - 1 - i) * CHR_SEP for i, c in enumerate(sorted_chroms)}
    y_span = (len(sorted_chroms) - 1) * CHR_SEP

    # REARRANGEMENT edges — always present, one per Rearrangement number.
    rearr_arc_data = []
    for rearr_num, grp in chain_rows.groupby('Rearrangement number'):
        bps = grp['Breakpoint number'].astype(int).tolist()
        node_a, node_b = bp_to_node[bps[0]], bp_to_node[bps[1]]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya, yb = chr_y[node_a[0]], chr_y[node_b[0]]
        rearr_arc_data.append((xa, ya, xb, yb, node_a[0] != node_b[0]))

    # OUR OWN INFERRED REFERENCE edges — only drawn where the target node
    # also has a REAL recorded breakpoint in this chain (a dangling/open
    # end is deliberately left undrawn here, not connected to a shadow
    # position — its status is reported in horizon_breakpoint_edges.csv
    # instead). Deduped at the NODE-pair level (not bp-number) so a
    # colliding node with multiple breakpoints doesn't draw the same arc
    # more than once.
    seen_ref_node_pairs = set()
    ref_arc_data = []
    resolved_bp_numbers = set()
    for node, bps_here in node_to_bps.items():
        target = ref_adj.get(node)
        if target == 'TELOMERE' or target not in node_to_bps:
            continue
        resolved_bp_numbers.update(bps_here)
        pair_key = tuple(sorted([node, target]))
        if pair_key in seen_ref_node_pairs:
            continue
        seen_ref_node_pairs.add(pair_key)
        node_a, node_b = node, target
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya = chr_y[node_a[0]]
        ref_arc_data.append((xa, ya, xb, ya))

    rads = _assign_rads(rearr_arc_data, [a[4] for a in rearr_arc_data])

    b_i = {c: len(chr_to_pos[c]) for c in sorted_chroms}
    theta = f"Theta({len(sorted_chroms)},({','.join(str(b_i[c]) for c in sorted_chroms)}))"
    n_breakpoints = len(chain_rows)
    n_dangling = len(bp_to_node) - len(resolved_bp_numbers)
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

    table_w, table_h, n_col_groups, table_max_rows_per_col, table_fontsize = _bp_table_dims(n_breakpoints, fig_h)

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_facecolor('white')
    x_margin = (X_AXES_RANGE - 1.0) / 2
    ax.set_xlim(-x_margin, 1.0 + x_margin)
    ax.set_ylim(bottom_y - 2.0 - extra_data_bottom, top_y + 2.0 + extra_data_top)
    ax.axis('off')

    ax.set_title(
        f'{patient_id} — HORIZON Chain {chain_number}{erg_tag}\n'
        f'{theta}   n_breakpoints={n_breakpoints}   n_rearrangement_edges={len(rearr_arc_data)}'
        f'   n_reference_edges_resolved={len(ref_arc_data)}   n_dangling_ends={n_dangling}\n'
        f'Rearrangement + our own inferred reference edges ONLY — '
        f'adjacency and deletion-bridge edges NOT used (HORIZON branch)',
        fontsize=11.5, fontweight='bold', pad=14, loc='center',
    )

    for c in sorted_chroms:
        draw_chr_track(ax, c, chr_y[c], chr_to_pos[c], nmaps[c])
        draw_bp_number_labels(ax, c, chr_y[c], chr_to_rows[c], nmaps[c])

    for xa, ya, xb, yb in ref_arc_data:
        _arc(ax, xa, ya, xb, yb, -0.12, REF_EDGE_COLOR, 1.4, '-', 0.7, '-', 3)

    for (xa, ya, xb, yb, is_inter), rad in zip(rearr_arc_data, rads):
        col = COL['rearr_inter'] if is_inter else COL['rearr_intra']
        _arc(ax, xa, ya, xb, yb, rad, col, 2.0, '--', 0.9, '-', 5)

    fig.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Rearrangement edge — intrachromosomal'),
        mpatches.Patch(color=COL['rearr_inter'], label='Rearrangement edge — interchromosomal'),
        mpatches.Patch(color=REF_EDGE_COLOR,     label='Our inferred reference edge (resolved)'),
    ], loc='lower center', ncol=3, fontsize=10, framealpha=0.95, edgecolor='#CCCCCC',
       bbox_to_anchor=(0.5, 0.0))

    fig.subplots_adjust(left=0.04, right=max(0.96 - table_w, 0.6), top=0.88, bottom=0.08)
    _add_bp_location_partner_table(fig, chain_rows, bp_rearrangement_partner,
                                    table_w, table_h, n_col_groups, table_max_rows_per_col, table_fontsize)

    out_path = os.path.join(save_dir, f'{patient_id}_horizon_chain_{chain_number}.png')
    fig.savefig(out_path, dpi=150, facecolor='white')
    plt.close(fig)
    return out_path, n_dangling, len(ref_arc_data)


# ── Batch runner ─────────────────────────────────────────────────────────────

def run_all():
    os.makedirs(CSV_DIR, exist_ok=True)
    os.makedirs(CHAIN_IMG_DIR, exist_ok=True)

    df = load_mmc5_chains()
    clinical_df = pd.read_csv(BACA_CLINICAL_PHENO_FILE_PATH)
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])

    chain_summary_rows = []
    bp_edge_rows = []
    n_drawn, n_failed = 0, 0

    for (patient_id, chain_number), chain_rows in df.groupby(['Individual', 'Chain number']):
        try:
            bp_to_node, node_to_bps, rearr_adj, ref_adj, bp_rearr_partner, bp_ref_partner = compute_chain_edges(chain_rows)

            out_path, n_dangling, n_ref_resolved = draw_horizon_chain(
                patient_id, chain_number, chain_rows, bp_to_node, node_to_bps,
                ref_adj, bp_rearr_partner, CHAIN_IMG_DIR,
            )

            # OPEN ENDS — locked definition (confirmed with user 2026-08-24,
            # see CLAUDE.md HORIZON section): every breakpoint gap has 2
            # ends (a '-' end and a '+' end). The recorded end (the one in
            # Baca's data) always has a rearrangement edge, and every end
            # always has a reference edge (real neighbor, or a telomere
            # edge "with respect to this chain" if nothing further exists
            # in this chain's own data) — so the recorded end is NEVER
            # open. The breakpoint's OTHER end (same position, opposite
            # strand) was never independently recorded, so it has no
            # rearrangement edge — it is open, UNLESS that exact opposite
            # end happens to ALSO be independently recorded as its own
            # separate real breakpoint elsewhere in this chain (checked via
            # node_to_bps, not assumed).
            other_strand = {'+': '-', '-': '+'}
            open_bp_strand = []  # list of (bp_num, its_open_strand)
            for _, r in chain_rows.iterrows():
                bp_num = int(r['Breakpoint number'])
                chrom, pos, strand = int(r['chromosome']), r['position'], r['node_strand']
                opp_strand = other_strand[strand]
                if (chrom, pos, opp_strand) not in node_to_bps:
                    open_bp_strand.append((bp_num, opp_strand))
            open_bp_strand.sort(key=lambda t: t[0])
            open_ends_bp_strand_list = ", ".join(f'BP{bp}({s})' for bp, s in open_bp_strand)
            n_open_ends = len(open_bp_strand)

            b_i = chain_rows.groupby('chromosome').size().to_dict()
            chroms = sorted(b_i.keys())
            site_col = chain_rows['Site annotation'] if 'Site annotation' in chain_rows.columns else pd.Series([], dtype=str)
            contains_erg = bool(site_col.dropna().str.contains(ERG_PATTERN).any())
            genes = extract_genes(site_col.tolist())

            clin = clinical_df[clinical_df['Individual'] == patient_id]
            gleason = clin['Gleason Score'].values[0] if len(clin) > 0 else None
            stage = clin['Pathological stage'].values[0] if len(clin) > 0 else None

            chain_summary_rows.append({
                'patient_id': patient_id,
                'baca_chain_number': chain_number,
                'n_breakpoints': len(chain_rows),
                'n_rearrangements': chain_rows['Rearrangement number'].nunique(),
                'k_chromosomes': len(chroms),
                'chromosomes': ",".join(str(c) for c in chroms),
                'n_reference_edges_resolved': n_ref_resolved,
                'n_open_ends': n_open_ends,
                'open_ends_bp_strand_list': open_ends_bp_strand_list,
                'contains_ERG_or_TMPRSS2': contains_erg,
                'genes_in_chain': ",".join(genes),
                'ETS_status': ets_lookup.get(patient_id),
                'Gleason_Score': gleason,
                'Pathological_stage': stage,
            })

            open_lookup = dict(open_bp_strand)  # bp_num -> its open strand, only for open bp's
            for _, r in chain_rows.iterrows():
                bp_num = int(r['Breakpoint number'])
                strand = r['node_strand']
                is_open = bp_num in open_lookup
                bp_edge_rows.append({
                    'patient_id': patient_id,
                    'baca_chain_number': chain_number,
                    'breakpoint_number': bp_num,
                    'chromosome': int(r['chromosome']),
                    'position': r['position'],
                    'strand': strand,
                    'rearrangement_partner_bp': bp_rearr_partner[bp_num],
                    'reference_edge_target': bp_ref_partner[bp_num],
                    'other_end_strand': other_strand[strand],
                    'other_end_is_open': is_open,
                })

            n_drawn += 1
            print(f'{patient_id} chain {chain_number} -> {out_path}')
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    chain_summary_df = pd.DataFrame(chain_summary_rows).sort_values(['patient_id', 'baca_chain_number']).reset_index(drop=True)
    bp_edges_df = pd.DataFrame(bp_edge_rows).sort_values(['patient_id', 'baca_chain_number', 'breakpoint_number']).reset_index(drop=True)

    chain_summary_df.to_csv(CHAIN_SUMMARY_CSV, index=False)
    bp_edges_df.to_csv(BREAKPOINT_EDGES_CSV, index=False)

    print(f'\nChains drawn: {n_drawn}  Failed: {n_failed}  -> {CHAIN_IMG_DIR}')
    print(f'Wrote {len(chain_summary_df)} chains -> {CHAIN_SUMMARY_CSV}')
    print(f'Wrote {len(bp_edges_df)} breakpoints (their OTHER end tracked alongside) -> {BREAKPOINT_EDGES_CSV}')
    print(f'\nTotal open ends cohort-wide: {bp_edges_df["other_end_is_open"].sum()} / {len(bp_edges_df)} breakpoints have an open other end')


if __name__ == '__main__':
    run_all()
