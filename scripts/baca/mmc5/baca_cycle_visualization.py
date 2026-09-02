"""
baca_cycle_visualization.py — ONE IMAGE PER BACA CHAIN THAT HAS AT LEAST
ONE CLOSED CYCLE. Each closed cycle gets its own panel, but — per explicit
user request — the cycle is drawn ON TOP of the FULL chain (grayed out as
background context), not cropped into an isolated box showing only the
cycle's own chromosomes. A start/end marker is placed on the cycle since
it's a closed loop (returns to where it started).

Companion to baca_chain_visualization.py (draws the FULL chain, every real
edge, no cycle highlighting, all 366 chains, no fading). This script only
covers the 91/366 chains with >=1 closed cycle (per
results/mmc5_primary_chain_cycle_structures.csv /
baca_mmc5_chain_closure.py's has_any_closed_cycle flag).

Design, per user feedback on the first version of this script:
  - The full chain (ALL its breakpoints/chromosomes/edges) is always drawn
    as a faded gray backdrop, using the SAME coordinate layout in every
    panel for that chain — so you see exactly where in the whole chain
    each cycle sits, not an isolated, cropped fragment.
  - The cycle's own edges are drawn on top in full bright color (same
    color convention as baca_chain_visualization.py: rearrangement
    orange/blue, adjacency green, deletion bridge purple dotted).
  - A start/end marker (the cycle's lowest breakpoint number, arbitrary
    but deterministic choice) is placed and labelled 'START / END' —
    since an AMG cycle has no inherent direction, this just marks where
    the loop closes, not a true biological start.
  - No bounding box around panels — removed per explicit request.

A "closed cycle" here is a connected component of (real rearrangement
edge) UNION (real chosen adjacency-or-deletion-bridge reference edge)
where every node has degree exactly 2 — same definition, same algorithm,
as baca_mmc5_chain_closure.py / baca_mmc5_primary_chain_cycle_structures.py.
Nothing invented: every edge drawn here (faded background AND bright
highlight) is one of mmc5's 3 real edge types.

Output: results/baca_proposed_cycles/{patientID}_cycle_chain_{chainNumber}.png
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
from baca_mmc5_chain_closure import (
    _build_rearr_matching,
    _build_ref_graph,
    _cycle_structure_string,
    _parse_adjacent,
    _parse_deletion_bridge,
    MMC5_FILE_PATH,
)
from baca_chain_visualization import DELETION_BRIDGE_COLOR, STRAND_MAP
from baca_aberration_cycle_drawing import (
    build_norm_map,
    node_x,
    draw_chr_track,
    draw_bp_labels,
    _arc,
    _assign_rads,
    COL,
    CHR_SEP,
)

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PRIMARY_CSV = os.path.join(RESULTS_DIR, "mmc5_primary_chain_cycle_structures.csv")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "baca_proposed_cycles")

FADED_EDGE_COLOR = '#D5D5D5'


def load_mmc5_full():
    """Combines everything this script needs: chromosome (int), position
    (float), node_strand ('+'/'-'), adj_list, delb. Neither
    baca_mmc5_chain_closure.load_mmc5_table_s5a (no position/strand) nor
    baca_chain_visualization.load_mmc5 (no adj_list/delb) has both alone."""
    df = pd.ExcelFile(MMC5_FILE_PATH).parse('Table S5A')
    df.columns = [c.strip() for c in df.columns]
    df = df[df['Chromosome:position'].notna()].copy()
    df['chromosome'] = df['Chromosome:position'].str.split(':').str[0].str.replace('chr', '', regex=False).astype(int)
    df['position'] = df['Chromosome:position'].str.split(':').str[1].astype(float)
    df['node_strand'] = df['Strand'].map(STRAND_MAP)
    df['adj_list'] = df['Adjacent breakpoint(s)'].apply(_parse_adjacent)
    df['delb'] = df['Deletion bridge partner breakpoint'].apply(_parse_deletion_bridge)
    return df


# ── Identify closed-cycle components + their real edges (typed) ────────────

def _build_edge_type_maps(grp):
    """frozenset({a,b}) -> 'adjacency' or 'deletion_bridge' (a pair flagged
    as both keeps the deletion_bridge tag — the more specific real edge)."""
    adjacency_pairs = set()
    delb_pairs = set()
    for _, row in grp.iterrows():
        bp = row['Breakpoint number']
        for a in row['adj_list']:
            adjacency_pairs.add(frozenset((bp, a)))
        if row['delb'] is not None:
            delb_pairs.add(frozenset((bp, row['delb'])))
    return adjacency_pairs, delb_pairs


def find_closed_cycles(grp):
    """
    Returns a list of dicts, one per closed cycle found in this chain:
      {'breakpoints': [...], 'rearr_edges': [(a,b),...],
       'ref_edges': [(a,b,'adjacency'|'deletion_bridge'),...]}
    Same algorithm as baca_mmc5_chain_closure.compute_chain_closure, just
    also returning per-cycle membership instead of only aggregate counts.
    """
    nodes = set(grp['Breakpoint number'])
    rearr_adj = _build_rearr_matching(grp)
    ref_graph = _build_ref_graph(grp)
    matching = nx.algorithms.matching.max_weight_matching(ref_graph, maxcardinality=True)
    adjacency_pairs, delb_pairs = _build_edge_type_maps(grp)

    # MUST be a MultiGraph — see baca_mmc5_chain_closure.py docstring (a C2
    # cycle's rearrangement edge and reference edge can connect the exact
    # same two nodes; a plain Graph would collapse that to degree 1).
    combined = nx.MultiGraph()
    combined.add_nodes_from(nodes)
    seen_rearr_pairs = set()
    for a, b in rearr_adj.items():
        key = tuple(sorted((a, b)))
        if key not in seen_rearr_pairs:
            combined.add_edge(a, b)
            seen_rearr_pairs.add(key)
    for a, b in matching:
        combined.add_edge(a, b)

    cycles = []
    for component in nx.connected_components(combined):
        degrees = [combined.degree(n) for n in component]
        is_closed = len(component) >= 2 and all(d == 2 for d in degrees)
        if not is_closed:
            continue

        rearr_edges = [(a, b) for a, b in rearr_adj.items()
                       if a in component and b in component and a < b]
        ref_edges = []
        for a, b in matching:
            if a in component and b in component:
                pair = frozenset((a, b))
                edge_type = 'deletion_bridge' if pair in delb_pairs else 'adjacency'
                ref_edges.append((a, b, edge_type))

        cycles.append({
            'breakpoints': sorted(component),
            'rearr_edges': rearr_edges,
            'ref_edges': ref_edges,
        })

    # Largest cycle first — most interesting panel drawn first.
    cycles.sort(key=lambda c: -len(c['breakpoints']))
    return cycles


# ── Full-chain layout + faded background (shared across all panels of a chain) ──

def _build_chain_layout(grp):
    bp_to_node = {
        int(r['Breakpoint number']): (int(r['chromosome']), r['position'], r['node_strand'])
        for _, r in grp.iterrows()
    }
    chr_to_pos = {}
    for node in bp_to_node.values():
        chr_to_pos.setdefault(node[0], set()).add(node[1])
    chr_to_pos = {c: sorted(positions) for c, positions in chr_to_pos.items()}
    sorted_chroms = sorted(chr_to_pos.keys())
    nmaps = {c: build_norm_map(chr_to_pos[c]) for c in sorted_chroms}
    k = len(sorted_chroms)
    chr_y = {c: (k - 1 - i) * CHR_SEP for i, c in enumerate(sorted_chroms)}

    # Full chain's real edges (ALL of them, not just the ones the closure
    # algorithm happened to select) — same edge set baca_chain_visualization.py
    # draws for the full-chain image, used here as faded background context.
    rearr_arc_data = []
    n_unresolved_rearr = 0
    for rearr_num, sub in grp.groupby('Rearrangement number'):
        bps = sub['Breakpoint number'].astype(int).tolist()
        if len(bps) != 2:
            n_unresolved_rearr += len(bps)
            continue
        node_a, node_b = bp_to_node[bps[0]], bp_to_node[bps[1]]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya, yb = chr_y[node_a[0]], chr_y[node_b[0]]
        rearr_arc_data.append((xa, ya, xb, yb, node_a[0] != node_b[0]))
    rads = _assign_rads(rearr_arc_data, [a[4] for a in rearr_arc_data])

    seen_ref_pairs, seen_delb_pairs = set(), set()
    ref_arc_data, delb_arc_data = [], []
    for _, r in grp.iterrows():
        bp_num = int(r['Breakpoint number'])
        node_a = bp_to_node[bp_num]
        for adj_num in r['adj_list']:
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
        delb_num = r['delb']
        if delb_num is not None and delb_num in bp_to_node:
            pair_key = tuple(sorted([bp_num, delb_num]))
            if pair_key not in seen_delb_pairs:
                seen_delb_pairs.add(pair_key)
                node_b = bp_to_node[delb_num]
                xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
                xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
                ya = chr_y[node_a[0]]
                delb_arc_data.append((xa, ya, xb, ya))

    return {
        'bp_to_node': bp_to_node,
        'chr_to_pos': chr_to_pos,
        'sorted_chroms': sorted_chroms,
        'nmaps': nmaps,
        'chr_y': chr_y,
        'rearr_arc_data': rearr_arc_data,
        'rads': rads,
        'ref_arc_data': ref_arc_data,
        'delb_arc_data': delb_arc_data,
    }


def _draw_faded_full_chain(ax, layout):
    for c in layout['sorted_chroms']:
        draw_chr_track(ax, c, layout['chr_y'][c], layout['chr_to_pos'][c], layout['nmaps'][c], faded=True)
        draw_bp_labels(ax, c, layout['chr_y'][c], layout['chr_to_pos'][c], layout['nmaps'][c], {}, faded=True)
    for xa, ya, xb, yb in layout['ref_arc_data']:
        _arc(ax, xa, ya, xb, yb, -0.12, FADED_EDGE_COLOR, 1.2, '-', 0.55, '-', 2)
    for xa, ya, xb, yb in layout['delb_arc_data']:
        _arc(ax, xa, ya, xb, yb, 0.18, FADED_EDGE_COLOR, 1.2, ':', 0.55, '-', 2)
    for (xa, ya, xb, yb, is_inter), rad in zip(layout['rearr_arc_data'], layout['rads']):
        _arc(ax, xa, ya, xb, yb, rad, FADED_EDGE_COLOR, 1.4, '--', 0.55, '-', 2)


def _start_end_marker(ax, node, nmaps, chr_y, color):
    sx = node_x(node[1], node[2], nmaps[node[0]])
    sy = chr_y[node[0]]
    ax.plot(sx, sy, 'o', color=color, markersize=13,
            markeredgecolor='white', markeredgewidth=2.2, zorder=10)
    ax.text(sx, sy + 0.55, 'START / END', ha='center', va='bottom',
            fontsize=8, fontweight='bold', color=color, zorder=10)


# ── Draw one cycle, highlighted on top of the full faded chain ─────────────

def draw_cycle_panel(ax, cycle, layout, cycle_index):
    bp_to_node = layout['bp_to_node']
    nmaps, chr_y = layout['nmaps'], layout['chr_y']
    breakpoints = cycle['breakpoints']

    _draw_faded_full_chain(ax, layout)

    # Highlight: this cycle's own rearrangement edges, in full bright color,
    # drawn on top of (zorder > ) the faded background.
    cycle_rearr_arc_data = []
    for a, b in cycle['rearr_edges']:
        node_a, node_b = bp_to_node[a], bp_to_node[b]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya, yb = chr_y[node_a[0]], chr_y[node_b[0]]
        cycle_rearr_arc_data.append((xa, ya, xb, yb, node_a[0] != node_b[0]))
    cycle_rads = _assign_rads(cycle_rearr_arc_data, [a[4] for a in cycle_rearr_arc_data])

    for a, b, edge_type in cycle['ref_edges']:
        node_a, node_b = bp_to_node[a], bp_to_node[b]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya = chr_y[node_a[0]]
        if edge_type == 'deletion_bridge':
            _arc(ax, xa, ya, xb, ya, 0.18, DELETION_BRIDGE_COLOR, 2.2, ':', 0.95, '-', 7)
        else:
            _arc(ax, xa, ya, xb, ya, -0.12, COL['cycle_ref'], 2.0, '-', 0.95, '-', 6)

    for (xa, ya, xb, yb, is_inter), rad in zip(cycle_rearr_arc_data, cycle_rads):
        col = COL['rearr_inter'] if is_inter else COL['rearr_intra']
        _arc(ax, xa, ya, xb, yb, rad, col, 2.8, '--', 1.0, '-', 8)

    # Re-draw bp tick marks + (+/-) annotation for THIS cycle's own
    # breakpoints in full saturation, on top of the faded background ticks,
    # so the cycle's own nodes are visually distinguishable from the rest
    # of the chain even before considering edge color.
    cycle_active = set(bp_to_node[bp] for bp in breakpoints)
    for c in layout['sorted_chroms']:
        cycle_pos_here = sorted({bp_to_node[bp][1] for bp in breakpoints if bp_to_node[bp][0] == c})
        if cycle_pos_here:
            draw_bp_labels(ax, c, chr_y[c], cycle_pos_here, nmaps[c], cycle_active)

    # Start/end marker — the cycle's lowest breakpoint number (arbitrary
    # but deterministic; an AMG cycle has no inherent direction, this just
    # marks where the closed loop returns to itself).
    start_bp = min(breakpoints)
    _start_end_marker(ax, bp_to_node[start_bp], nmaps, chr_y, '#000000')

    n_breakpoints = len(breakpoints)
    cycle_chroms = sorted({bp_to_node[bp][0] for bp in breakpoints})
    ax.set_title(
        f'Cycle {cycle_index}: C{n_breakpoints}   (chromosomes in this cycle: {",".join(str(c) for c in cycle_chroms)})',
        fontsize=12, fontweight='bold', pad=10,
    )


# ── Draw all closed cycles for one chain into one figure ───────────────────

def draw_chain_cycles(patient_id, chain_number, grp, save_dir):
    cycles = find_closed_cycles(grp)
    if not cycles:
        return None

    layout = _build_chain_layout(grp)
    sorted_chroms = layout['sorted_chroms']
    k = len(sorted_chroms)
    n_cycles = len(cycles)

    fig_w_per_panel = 12.0
    fig_w = max(fig_w_per_panel * n_cycles, 14.0)
    fig_h0 = max(7.0, k * 3.0 + 3.0)

    # Row-aware arc-clearance margin (top-most / bottom-most row only —
    # see baca_chain_visualization.py for why), computed from the FULL
    # chain's edges since that's what's drawn (faded) in every panel.
    scale_x = fig_w_per_panel / 1.12
    y_span_layout = max((k - 1) * CHR_SEP, CHR_SEP)
    scale_y0 = fig_h0 / (y_span_layout * 1.05 + 4.0)
    top_y, bottom_y = max(layout['chr_y'].values()), min(layout['chr_y'].values())
    same_row_arcs = (
        [(ya, abs(rad) * abs(xb - xa) * scale_x) for (xa, ya, xb, yb, is_inter), rad in zip(layout['rearr_arc_data'], layout['rads']) if not is_inter]
        + [(ya, abs(-0.12) * abs(xb - xa) * scale_x) for xa, ya, xb, yb in layout['ref_arc_data']]
        + [(ya, abs(0.18) * abs(xb - xa) * scale_x) for xa, ya, xb, yb in layout['delb_arc_data']]
    )
    def _worst_at(row_y):
        return max((s for ya, s in same_row_arcs if ya == row_y), default=0.0)
    extra_top = _worst_at(top_y) * 1.15
    extra_bottom = _worst_at(bottom_y) * 1.15
    fig_h = fig_h0 + extra_top + extra_bottom

    fig, axes = plt.subplots(1, n_cycles, figsize=(fig_w, fig_h))
    if n_cycles == 1:
        axes = [axes]

    for ax in axes:
        ax.set_facecolor('white')
        ax.set_xlim(-0.06, 1.06)
        ax.set_ylim(bottom_y - 1.8 - extra_bottom / scale_y0, top_y + 1.8 + extra_top / scale_y0)
        ax.axis('off')

    for i, (ax, cycle) in enumerate(zip(axes, cycles), start=1):
        draw_cycle_panel(ax, cycle, layout, i)

    overall_structure = _cycle_structure_string(sorted((len(c['breakpoints']) for c in cycles), reverse=True))
    fig.suptitle(
        f'{patient_id} — Baca Chain {chain_number} — closed cycle(s) highlighted on full chain\n'
        f'cycle_structure = {overall_structure}',
        fontsize=14, fontweight='bold', y=0.99,
    )

    fig.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Rearrangement edge — intrachromosomal'),
        mpatches.Patch(color=COL['rearr_inter'], label='Rearrangement edge — interchromosomal'),
        mpatches.Patch(color=COL['cycle_ref'],   label='Adjacency reference edge'),
        mpatches.Patch(color=DELETION_BRIDGE_COLOR, label='Deletion bridge edge'),
        mpatches.Patch(color=FADED_EDGE_COLOR, label='Full chain (faded background)'),
    ], loc='lower center', ncol=5, fontsize=9.5, framealpha=0.95, edgecolor='#CCCCCC',
       bbox_to_anchor=(0.5, 0.0))

    fig.subplots_adjust(left=0.03, right=0.97, top=0.86, bottom=0.10, wspace=0.12)
    out_path = os.path.join(save_dir, f'{patient_id}_cycle_chain_{chain_number}.png')
    fig.savefig(out_path, dpi=150, facecolor='white')
    plt.close(fig)
    return out_path


def run_all(save_dir=OUTPUT_DIR):
    os.makedirs(save_dir, exist_ok=True)
    m5_df = load_mmc5_full()
    primary_df = pd.read_csv(PRIMARY_CSV)
    closed_chains = primary_df[primary_df['has_any_closed_cycle']]

    n_drawn, n_failed = 0, 0
    for _, row in closed_chains.iterrows():
        patient_id, chain_number = row['patient_id'], row['baca_chain_number']
        grp = m5_df[(m5_df['Individual'] == patient_id) & (m5_df['Chain number'] == chain_number)]
        try:
            out = draw_chain_cycles(patient_id, chain_number, grp, save_dir)
            n_drawn += 1
            print(f'{patient_id} chain {chain_number} -> {out}')
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    print(f'\nDrawn: {n_drawn}  Failed: {n_failed}  -> {save_dir}')


if __name__ == '__main__':
    run_all()
