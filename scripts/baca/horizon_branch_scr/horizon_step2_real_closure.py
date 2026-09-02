"""
horizon_step2_real_closure.py — HORIZON branch, Step 2.
Updated 2026-08-24 to read Step 1's corrected open-end/reference-edge
columns (horizon_breakpoint_edges.csv's 'reference_edge_target', not the
retired 'inferred_reference_partner'/'MISSING_REARR' — see CLAUDE.md's
HORIZON section for the locked definition this now reads).

Real cycle closure — "cycle already exists in chain" — using ONLY
rearrangement edges + reference edges whose target is ALSO an
independently-recorded breakpoint. No invented edges anywhere in this
script; that is Step 3 (obligate/probable), separate and later.

IMPORTANT TERMINOLOGY NOTE — this script's "open"/"closed" is a DIFFERENT
concept from Step 1's "open end":
  - Step 1 "open end" = does this specific end have a rearrangement edge?
    (a per-END fact, about one side of one breakpoint gap)
  - Step 2 "closed cycle" membership = does this breakpoint's own
    (always-rearrangement-bearing) node participate in a closed loop when
    you alternate rearrangement/reference hops? (a per-BREAKPOINT fact,
    about whether a whole loop closes)
  A breakpoint's own end is essentially NEVER "open" in the Step 1 sense
  (see CLAUDE.md), but it can absolutely fail to be part of any closed
  cycle in the Step 2 sense — these are not the same question, and this
  script never uses the word "open" to mean Step 1's definition.

Method mirrors the established mmc5/baca_mmc5_chain_closure.py pattern:
build a per-chain nx.MultiGraph from the rearrangement + reference edges,
then classify each connected component by degree — every node degree
exactly 2 = closed cycle; anything else (a degree-1 node) = not in a
closed cycle. A plain nx.Graph is NOT used (locked project convention) —
it would silently collapse a C2 cycle where the rearrangement edge and
the reference edge connect the exact same two breakpoints, wrongly
marking it unclosed.

This graph-degree method is mathematically equivalent to the alternating
rearrangement/reference HOP traversal already established elsewhere in
this codebase (_find_amg_cycles in core/baca_aberration_clinical_analysis.py):
a reference edge is only included here when its target is independently
real, which is exactly the condition under which a hop-based traversal
would be able to continue past it rather than dead-ending — confirmed
directly with the user 2026-08-24 before re-running this script, rather
than assumed.

No matching algorithm is needed here (unlike mmc5_chain_closure.py's
adjacency+deletion-bridge layer, which can offer a breakpoint more than
one real candidate partner): Step 1's reference edges are already an
unambiguous 1-to-1 pairing (verified: 0/4,440 breakpoints have an
ambiguous '/'-joined reference partner in this cohort).

Outputs:
  results/horizon_branch_res/horizon_csv_data/horizon_real_closure_summary.csv
    — 1 row per chain (366 rows)
  results/horizon_branch_res/horizon_chain_cycle_visualizations/cycles/
    {patient_id}_horizon_cycle_chain_{n}.png — 1 image per chain with >=1
    real closed cycle (full chain faded in the background, closed cycle(s)
    highlighted on top, START/END marker) — same visual style as the
    existing mmc5/baca_cycle_visualization.py, minus the deletion-bridge
    edge type (doesn't exist in this branch).
"""

import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
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

RESULTS_ROOT = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results/horizon_branch_res"
CSV_DIR = os.path.join(RESULTS_ROOT, "horizon_csv_data")
CYCLE_IMG_DIR = os.path.join(RESULTS_ROOT, "horizon_chain_cycle_visualizations", "cycles")

BREAKPOINT_EDGES_CSV = os.path.join(CSV_DIR, "horizon_breakpoint_edges.csv")
CHAIN_SUMMARY_CSV = os.path.join(CSV_DIR, "horizon_chain_structure_summary.csv")
REAL_CLOSURE_CSV = os.path.join(CSV_DIR, "horizon_real_closure_summary.csv")

FADED_EDGE_COLOR = '#D5D5D5'
REF_EDGE_COLOR = COL['cycle_ref']  # same green used throughout this codebase


# ── Closure computation (breakpoint-number graph, mirrors mmc5_chain_closure.py) ──

_PURE_INT_RE = re.compile(r'^\d+$')


def _real_ref_pairs(chain_bp):
    """
    Deduped (a,b) breakpoint-number pairs for every reference edge whose
    target IS a real, independently-recorded breakpoint in this chain
    (reads Step 1's 'reference_edge_target' column, post-2026-08-24 open-
    end fix: values are a plain bp-number string when real, or
    'TELOMERE (this chain)' / 'chr..(strand) [no separate BP# here]' when
    not — only plain bp-number values count here). This is exactly the
    traversal-continuability criterion: a reference hop only lets a cycle
    continue if it lands on a node that itself has a rearrangement edge —
    which is precisely what an independently-recorded target means.
    Reciprocal by construction (a's partner is b iff b's partner is a) —
    verified during Step 1.
    """
    seen = set()
    pairs = []
    for _, r in chain_bp.iterrows():
        target = str(r['reference_edge_target'])
        if not _PURE_INT_RE.match(target):
            continue
        a, b = int(r['breakpoint_number']), int(target)
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        pairs.append(key)
    return pairs


def _real_rearr_pairs(chain_bp):
    seen = set()
    pairs = []
    for _, r in chain_bp.iterrows():
        a, b = int(r['breakpoint_number']), int(r['rearrangement_partner_bp'])
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        pairs.append(key)
    return pairs


def _cycle_structure_string(cycle_lengths):
    if not cycle_lengths:
        return 'none'
    from collections import Counter
    counts = Counter(cycle_lengths)
    parts = [f"C{length}x{n}" for length, n in sorted(counts.items(), reverse=True)]
    return " + ".join(parts)


def find_closed_cycles(chain_bp):
    """
    Returns a list of dicts, one per closed cycle found in this chain:
      {'breakpoints': [...], 'rearr_edges': [(a,b),...], 'ref_edges': [(a,b),...]}
    """
    nodes = set(chain_bp['breakpoint_number'].astype(int))
    rearr_pairs = _real_rearr_pairs(chain_bp)
    ref_pairs = _real_ref_pairs(chain_bp)

    combined = nx.MultiGraph()
    combined.add_nodes_from(nodes)
    for a, b in rearr_pairs:
        combined.add_edge(a, b)
    for a, b in ref_pairs:
        combined.add_edge(a, b)

    cycles = []
    for component in nx.connected_components(combined):
        degrees = [combined.degree(n) for n in component]
        is_closed = len(component) >= 2 and all(d == 2 for d in degrees)
        if not is_closed:
            continue
        rearr_edges = [(a, b) for a, b in rearr_pairs if a in component and b in component]
        ref_edges = [(a, b) for a, b in ref_pairs if a in component and b in component]
        cycles.append({'breakpoints': sorted(component), 'rearr_edges': rearr_edges, 'ref_edges': ref_edges})

    cycles.sort(key=lambda c: -len(c['breakpoints']))
    return cycles


def compute_chain_closure(chain_bp):
    nodes = set(chain_bp['breakpoint_number'].astype(int))
    n_breakpoints = len(nodes)
    chrom_lookup = dict(zip(chain_bp['breakpoint_number'].astype(int), chain_bp['chromosome'].astype(int)))

    cycles = find_closed_cycles(chain_bp)
    cycle_lengths = [len(c['breakpoints']) for c in cycles]
    cycle_chrom_spans = [len({chrom_lookup[bp] for bp in c['breakpoints']}) for c in cycles]
    n_in_cycles = sum(cycle_lengths)
    n_open_breakpoints = n_breakpoints - n_in_cycles

    return {
        'cycles': cycles,
        'has_any_closed_cycle': len(cycle_lengths) > 0,
        'chain_fully_closed': len(cycle_lengths) > 0 and n_open_breakpoints == 0,
        'n_open_breakpoints': n_open_breakpoints,
        'pct_breakpoints_in_closed_cycles': round(100 * n_in_cycles / n_breakpoints, 1),
        'cycle_structure': _cycle_structure_string(cycle_lengths),
        'n_cycles': len(cycle_lengths),
        'max_cycle_chrom_span': max(cycle_chrom_spans) if cycle_chrom_spans else 0,
        'has_chromoplexy_strict': bool(cycle_chrom_spans) and max(cycle_chrom_spans) >= 3,
        'has_chromoplexy_loose': bool(cycle_chrom_spans) and max(cycle_chrom_spans) >= 2,
    }


# ── Drawing: full chain faded backdrop + closed cycle(s) highlighted ─────────

def _build_chain_layout(chain_bp):
    bp_to_node = {
        int(r['breakpoint_number']): (int(r['chromosome']), r['position'], r['strand'])
        for _, r in chain_bp.iterrows()
    }
    chr_to_pos = {}
    for node in bp_to_node.values():
        chr_to_pos.setdefault(node[0], set()).add(node[1])
    chr_to_pos = {c: sorted(positions) for c, positions in chr_to_pos.items()}
    sorted_chroms = sorted(chr_to_pos.keys())
    nmaps = {c: build_norm_map(chr_to_pos[c]) for c in sorted_chroms}
    k = len(sorted_chroms)
    chr_y = {c: (k - 1 - i) * CHR_SEP for i, c in enumerate(sorted_chroms)}

    rearr_pairs = _real_rearr_pairs(chain_bp)
    rearr_arc_data = []
    for a, b in rearr_pairs:
        node_a, node_b = bp_to_node[a], bp_to_node[b]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya, yb = chr_y[node_a[0]], chr_y[node_b[0]]
        rearr_arc_data.append((xa, ya, xb, yb, node_a[0] != node_b[0]))
    rads = _assign_rads(rearr_arc_data, [a[4] for a in rearr_arc_data])

    ref_pairs = _real_ref_pairs(chain_bp)
    ref_arc_data = []
    for a, b in ref_pairs:
        node_a, node_b = bp_to_node[a], bp_to_node[b]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya = chr_y[node_a[0]]
        ref_arc_data.append((xa, ya, xb, ya))

    return {
        'bp_to_node': bp_to_node, 'chr_to_pos': chr_to_pos, 'sorted_chroms': sorted_chroms,
        'nmaps': nmaps, 'chr_y': chr_y, 'rearr_arc_data': rearr_arc_data, 'rads': rads,
        'ref_arc_data': ref_arc_data,
    }


def _draw_faded_full_chain(ax, layout):
    for c in layout['sorted_chroms']:
        draw_chr_track(ax, c, layout['chr_y'][c], layout['chr_to_pos'][c], layout['nmaps'][c], faded=True)
        draw_bp_labels(ax, c, layout['chr_y'][c], layout['chr_to_pos'][c], layout['nmaps'][c], {}, faded=True)
    for xa, ya, xb, yb in layout['ref_arc_data']:
        _arc(ax, xa, ya, xb, yb, -0.12, FADED_EDGE_COLOR, 1.2, '-', 0.55, '-', 2)
    for (xa, ya, xb, yb, is_inter), rad in zip(layout['rearr_arc_data'], layout['rads']):
        _arc(ax, xa, ya, xb, yb, rad, FADED_EDGE_COLOR, 1.4, '--', 0.55, '-', 2)


def _start_end_marker(ax, node, nmaps, chr_y, color):
    sx = node_x(node[1], node[2], nmaps[node[0]])
    sy = chr_y[node[0]]
    ax.plot(sx, sy, 'o', color=color, markersize=13,
            markeredgecolor='white', markeredgewidth=2.2, zorder=10)
    ax.text(sx, sy + 0.55, 'START / END', ha='center', va='bottom',
            fontsize=8, fontweight='bold', color=color, zorder=10)


def draw_cycle_panel(ax, cycle, layout, cycle_index):
    bp_to_node = layout['bp_to_node']
    nmaps, chr_y = layout['nmaps'], layout['chr_y']
    breakpoints = cycle['breakpoints']

    _draw_faded_full_chain(ax, layout)

    cycle_rearr_arc_data = []
    for a, b in cycle['rearr_edges']:
        node_a, node_b = bp_to_node[a], bp_to_node[b]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya, yb = chr_y[node_a[0]], chr_y[node_b[0]]
        cycle_rearr_arc_data.append((xa, ya, xb, yb, node_a[0] != node_b[0]))
    cycle_rads = _assign_rads(cycle_rearr_arc_data, [a[4] for a in cycle_rearr_arc_data])

    for a, b in cycle['ref_edges']:
        node_a, node_b = bp_to_node[a], bp_to_node[b]
        xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
        xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
        ya = chr_y[node_a[0]]
        _arc(ax, xa, ya, xb, ya, -0.12, REF_EDGE_COLOR, 2.0, '-', 0.95, '-', 6)

    for (xa, ya, xb, yb, is_inter), rad in zip(cycle_rearr_arc_data, cycle_rads):
        col = COL['rearr_inter'] if is_inter else COL['rearr_intra']
        _arc(ax, xa, ya, xb, yb, rad, col, 2.8, '--', 1.0, '-', 8)

    cycle_active = set(bp_to_node[bp] for bp in breakpoints)
    for c in layout['sorted_chroms']:
        cycle_pos_here = sorted({bp_to_node[bp][1] for bp in breakpoints if bp_to_node[bp][0] == c})
        if cycle_pos_here:
            draw_bp_labels(ax, c, chr_y[c], cycle_pos_here, nmaps[c], cycle_active)

    start_bp = min(breakpoints)
    _start_end_marker(ax, bp_to_node[start_bp], nmaps, chr_y, '#000000')

    n_breakpoints = len(breakpoints)
    cycle_chroms = sorted({bp_to_node[bp][0] for bp in breakpoints})
    ax.set_title(
        f'Cycle {cycle_index}: C{n_breakpoints}   (chromosomes in this cycle: {",".join(str(c) for c in cycle_chroms)})',
        fontsize=12, fontweight='bold', pad=10,
    )


def draw_chain_cycles(patient_id, chain_number, chain_bp, cycles, save_dir):
    layout = _build_chain_layout(chain_bp)
    sorted_chroms = layout['sorted_chroms']
    k = len(sorted_chroms)
    n_cycles = len(cycles)

    fig_w_per_panel = 12.0
    fig_w = max(fig_w_per_panel * n_cycles, 14.0)
    fig_h0 = max(7.0, k * 3.0 + 3.0)

    scale_x = fig_w_per_panel / 1.12
    y_span_layout = max((k - 1) * CHR_SEP, CHR_SEP)
    scale_y0 = fig_h0 / (y_span_layout * 1.05 + 4.0)
    top_y, bottom_y = max(layout['chr_y'].values()), min(layout['chr_y'].values())
    same_row_arcs = (
        [(ya, abs(rad) * abs(xb - xa) * scale_x) for (xa, ya, xb, yb, is_inter), rad in zip(layout['rearr_arc_data'], layout['rads']) if not is_inter]
        + [(ya, abs(-0.12) * abs(xb - xa) * scale_x) for xa, ya, xb, yb in layout['ref_arc_data']]
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
        f'{patient_id} — HORIZON Chain {chain_number} — closed cycle(s) highlighted on full chain\n'
        f'cycle_structure = {overall_structure}   (rearrangement + our own inferred reference edges only)',
        fontsize=13.5, fontweight='bold', y=0.99,
    )

    fig.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Rearrangement edge — intrachromosomal'),
        mpatches.Patch(color=COL['rearr_inter'], label='Rearrangement edge — interchromosomal'),
        mpatches.Patch(color=REF_EDGE_COLOR,     label='Our inferred reference edge'),
        mpatches.Patch(color=FADED_EDGE_COLOR,   label='Full chain (faded background)'),
    ], loc='lower center', ncol=4, fontsize=9.5, framealpha=0.95, edgecolor='#CCCCCC',
       bbox_to_anchor=(0.5, 0.0))

    fig.subplots_adjust(left=0.03, right=0.97, top=0.86, bottom=0.10, wspace=0.12)
    out_path = os.path.join(save_dir, f'{patient_id}_horizon_cycle_chain_{chain_number}.png')
    fig.savefig(out_path, dpi=150, facecolor='white')
    plt.close(fig)
    return out_path


# ── Batch runner ─────────────────────────────────────────────────────────────

def run_all():
    os.makedirs(CYCLE_IMG_DIR, exist_ok=True)

    bp_df = pd.read_csv(BREAKPOINT_EDGES_CSV)
    chain_summary_df = pd.read_csv(CHAIN_SUMMARY_CSV)
    clinical_cols = chain_summary_df.set_index(['patient_id', 'baca_chain_number'])[
        ['n_breakpoints', 'n_rearrangements', 'k_chromosomes', 'chromosomes',
         'contains_ERG_or_TMPRSS2', 'genes_in_chain', 'ETS_status', 'Gleason_Score', 'Pathological_stage']
    ]

    rows = []
    n_drawn, n_failed, n_with_cycle = 0, 0, 0

    for (patient_id, chain_number), chain_bp in bp_df.groupby(['patient_id', 'baca_chain_number']):
        try:
            result = compute_chain_closure(chain_bp)
            clin = clinical_cols.loc[(patient_id, chain_number)]

            rows.append({
                'patient_id': patient_id,
                'baca_chain_number': chain_number,
                'n_breakpoints': clin['n_breakpoints'],
                'n_rearrangements': clin['n_rearrangements'],
                'k_chromosomes': clin['k_chromosomes'],
                'chromosomes': clin['chromosomes'],
                'has_any_closed_cycle': result['has_any_closed_cycle'],
                'chain_fully_closed': result['chain_fully_closed'],
                'n_open_breakpoints': result['n_open_breakpoints'],
                'pct_breakpoints_in_closed_cycles': result['pct_breakpoints_in_closed_cycles'],
                'cycle_structure': result['cycle_structure'],
                'n_cycles': result['n_cycles'],
                'max_cycle_chrom_span': result['max_cycle_chrom_span'],
                'has_chromoplexy_strict': result['has_chromoplexy_strict'],
                'has_chromoplexy_loose': result['has_chromoplexy_loose'],
                'contains_ERG_or_TMPRSS2': clin['contains_ERG_or_TMPRSS2'],
                'genes_in_chain': clin['genes_in_chain'],
                'ETS_status': clin['ETS_status'],
                'Gleason_Score': clin['Gleason_Score'],
                'Pathological_stage': clin['Pathological_stage'],
            })

            if result['has_any_closed_cycle']:
                n_with_cycle += 1
                out = draw_chain_cycles(patient_id, chain_number, chain_bp, result['cycles'], CYCLE_IMG_DIR)
                print(f'{patient_id} chain {chain_number}: {result["cycle_structure"]} -> {out}')
            n_drawn += 1
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    summary_df = pd.DataFrame(rows).sort_values(['patient_id', 'baca_chain_number']).reset_index(drop=True)
    summary_df.to_csv(REAL_CLOSURE_CSV, index=False)

    print(f'\nProcessed: {n_drawn}  Failed: {n_failed}')
    print(f'Chains with >=1 real closed cycle (images drawn): {n_with_cycle} / {len(summary_df)}')
    print(f'Chains FULLY closed: {summary_df["chain_fully_closed"].sum()} / {len(summary_df)}')
    print(f'Wrote {len(summary_df)} chains -> {REAL_CLOSURE_CSV}')
    print()
    print('Chromoplexy flags among chains with >=1 closed cycle:')
    any_closed = summary_df[summary_df['has_any_closed_cycle']]
    print(f'  strict (>=3 chrom in one cycle): {any_closed["has_chromoplexy_strict"].sum()} / {len(any_closed)}')
    print(f'  loose  (>=2 chrom in one cycle): {any_closed["has_chromoplexy_loose"].sum()} / {len(any_closed)}')


if __name__ == '__main__':
    run_all()
