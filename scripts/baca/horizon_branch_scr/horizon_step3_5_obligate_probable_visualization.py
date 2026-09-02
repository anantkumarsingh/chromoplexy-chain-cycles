"""
horizon_step3_5_obligate_probable_visualization.py — HORIZON branch,
Step 3.5.

Visualizes the obligate and most-probable completed structures reported
in horizon_cycle_completion.csv (Step 3) — one image per chain per lens,
drawn on the same full-chain-faded-backdrop style established in Step 2's
cycle images (horizon_step2_real_closure.py), with THREE edge types now
visually distinguished:
  - Real edges already part of a closed cycle (Step 2's result) — bright
    blue/orange (rearrangement) + green (reference), same colors as before.
  - Real edges that only become part of a cycle once the true loose ends
    are closed (the internal edges of an open path component) — same
    bright colors, since they ARE real data, just newly recognized as
    part of a larger loop.
  - The INVENTED connection(s) themselves — purple dashed, clearly
    distinct from both of the above, so a viewer never mistakes an
    invented edge for real data.

Since "obligate" / "probable" name a STRUCTURE (a cycle-length pattern),
not a single specific pairing, and more than one pairing of the true
loose ends can achieve the same structure, this script recomputes one
REPRESENTATIVE matching per chain per lens (the first one the same
enumeration+selection functions from Step 3 report as achieving that
structure) and draws exactly that — never an abstract structure with no
edges shown.

Scope: every enumerable chain (m<=12) — ALL 359, per explicit user
instruction, not just the 276 that need invention. For the 83 chains
already fully closed by real data (0 true loose ends), obligate/probable
trivially equal the real structure with 0 invented edges — still drawn,
so every enumerable chain has its own obligate and probable image.
Non-enumerable chains (7) have no obligate/probable to draw.

Outputs:
  results/horizon_branch_res/horizon_chain_cycle_visualizations/
    obligate_structures/{patient_id}_horizon_chain_{n}_obligate.png
    probable_structures/{patient_id}_horizon_chain_{n}_probable.png
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'phase2'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'mmc5'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))

from horizon_step2_real_closure import (
    _real_rearr_pairs, _real_ref_pairs, find_closed_cycles,
    _build_chain_layout, _draw_faded_full_chain,
)
from horizon_step3_obligate_probable import _build_combined_graph
from baca_phase2_obligate_probable_completion import (
    _closed_and_open_components, _enumerate_completions,
    _select_obligate, _select_probable,
)
from baca_mmc5_chain_closure import _cycle_structure_string
from baca_aberration_cycle_drawing import node_x, draw_bp_labels, _arc, _assign_rads, COL, CHR_SEP

RESULTS_ROOT = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results/horizon_branch_res"
CSV_DIR = os.path.join(RESULTS_ROOT, "horizon_csv_data")
VIZ_ROOT = os.path.join(RESULTS_ROOT, "horizon_chain_cycle_visualizations")
OBLIGATE_DIR = os.path.join(VIZ_ROOT, "obligate_structures")
PROBABLE_DIR = os.path.join(VIZ_ROOT, "probable_structures")

BREAKPOINT_EDGES_CSV = os.path.join(CSV_DIR, "horizon_breakpoint_edges.csv")
COMPLETION_CSV = os.path.join(CSV_DIR, "horizon_cycle_completion.csv")

INVENTED_COLOR = '#9C27B0'  # purple — distinct from real rearrangement (blue/orange) and real reference (green)
FADED_EDGE_COLOR = '#D5D5D5'


def _component_internal_pairs(component, pairs):
    return [(a, b) for a, b in pairs if a in component and b in component]


def draw_completion(patient_id, chain_number, chain_bp, real_cycles,
                     open_components, rearr_pairs, ref_pairs, matching,
                     structure_str, lens_name, save_dir):
    layout = _build_chain_layout(chain_bp)
    bp_to_node, nmaps, chr_y = layout['bp_to_node'], layout['nmaps'], layout['chr_y']
    sorted_chroms = layout['sorted_chroms']
    k = len(sorted_chroms)

    def _to_arc_data(pairs):
        out = []
        for a, b in pairs:
            node_a, node_b = bp_to_node[a], bp_to_node[b]
            xa = node_x(node_a[1], node_a[2], nmaps[node_a[0]])
            xb = node_x(node_b[1], node_b[2], nmaps[node_b[0]])
            ya, yb = chr_y[node_a[0]], chr_y[node_b[0]]
            out.append((xa, ya, xb, yb, node_a[0] != node_b[0]))
        return out

    # ── Real edges: already-closed cycles + open-path internal edges ──
    bright_rearr_pairs = []
    bright_ref_pairs = []
    for cyc in real_cycles:
        bright_rearr_pairs.extend(cyc['rearr_edges'])
        bright_ref_pairs.extend(cyc['ref_edges'])
    for comp in open_components:
        bright_rearr_pairs.extend(_component_internal_pairs(comp, rearr_pairs))
        bright_ref_pairs.extend(_component_internal_pairs(comp, ref_pairs))

    bright_rearr_arc_data = _to_arc_data(bright_rearr_pairs)
    bright_rads = _assign_rads(bright_rearr_arc_data, [a[4] for a in bright_rearr_arc_data])
    bright_ref_arc_data = _to_arc_data(bright_ref_pairs)

    # ── Invented edges — the specific pairing chosen for this lens ──
    invented_arc_data = _to_arc_data(matching)
    invented_rads = _assign_rads(invented_arc_data, [a[4] for a in invented_arc_data], rad_min=0.30, rad_max=0.85)

    # ── Sizing — compute AFTER every arc type is known, so invented edges
    # (which can be the widest/deepest arcs on the page) are never clipped.
    fig_w = 16.0
    fig_h0 = max(7.0, k * 3.0 + 3.0)
    y_span_layout = max((k - 1) * CHR_SEP, CHR_SEP)
    scale_x = fig_w / 1.12
    scale_y0 = fig_h0 / (y_span_layout * 1.05 + 4.0)
    top_y, bottom_y = max(chr_y.values()), min(chr_y.values())
    same_row_arcs = (
        [(ya, abs(rad) * abs(xb - xa) * scale_x) for (xa, ya, xb, yb, is_inter), rad in zip(layout['rearr_arc_data'], layout['rads']) if not is_inter]
        + [(ya, abs(-0.12) * abs(xb - xa) * scale_x) for xa, ya, xb, yb in layout['ref_arc_data']]
        + [(ya, abs(rad) * abs(xb - xa) * scale_x) for (xa, ya, xb, yb, is_inter), rad in zip(invented_arc_data, invented_rads) if not is_inter]
    )
    def _worst_at(row_y):
        return max((s for ya, s in same_row_arcs if ya == row_y), default=0.0)
    extra_top = _worst_at(top_y) * 1.15
    extra_bottom = _worst_at(bottom_y) * 1.15
    fig_h = fig_h0 + extra_top + extra_bottom

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_facecolor('white')
    ax.set_xlim(-0.06, 1.06)
    ax.set_ylim(bottom_y - 1.8 - extra_bottom / scale_y0, top_y + 1.8 + extra_top / scale_y0)
    ax.axis('off')

    _draw_faded_full_chain(ax, layout)

    for xa, ya, xb, yb, _ in bright_ref_arc_data:
        _arc(ax, xa, ya, xb, ya, -0.12, COL['cycle_ref'], 2.0, '-', 0.95, '-', 6)

    for (xa, ya, xb, yb, is_inter), rad in zip(bright_rearr_arc_data, bright_rads):
        col = COL['rearr_inter'] if is_inter else COL['rearr_intra']
        _arc(ax, xa, ya, xb, yb, rad, col, 2.6, '--', 1.0, '-', 8)

    for (xa, ya, xb, yb, is_inter), rad in zip(invented_arc_data, invented_rads):
        _arc(ax, xa, ya, xb, yb, rad, INVENTED_COLOR, 2.6, ':', 1.0, '-|>', 9)

    # Re-label the true-loose-end nodes (the invented edges' own endpoints)
    # in full saturation on top of the faded backdrop, so they're visually
    # easy to find.
    loose_bps = sorted({bp for pair in matching for bp in pair})
    for c in sorted_chroms:
        pos_here = sorted({bp_to_node[bp][1] for bp in loose_bps if bp_to_node[bp][0] == c})
        if pos_here:
            active = {bp_to_node[bp] for bp in loose_bps}
            draw_bp_labels(ax, c, chr_y[c], pos_here, nmaps[c], active)

    ax.set_title(
        f'{patient_id} — HORIZON Chain {chain_number} — {lens_name.upper()} completion\n'
        f'{lens_name}_cycle_structure = {structure_str}   '
        f'({len(matching)} invented edge{"s" if len(matching) != 1 else ""})',
        fontsize=13, fontweight='bold', pad=12,
    )

    fig.legend(handles=[
        mpatches.Patch(color=COL['rearr_intra'], label='Real rearrangement edge — intrachromosomal'),
        mpatches.Patch(color=COL['rearr_inter'], label='Real rearrangement edge — interchromosomal'),
        mpatches.Patch(color=COL['cycle_ref'],   label='Real reference edge'),
        mpatches.Patch(color=INVENTED_COLOR,     label=f'INVENTED edge ({lens_name})'),
        mpatches.Patch(color=FADED_EDGE_COLOR,   label='Full chain (faded background)'),
    ], loc='lower center', ncol=3, fontsize=9.5, framealpha=0.95, edgecolor='#CCCCCC',
       bbox_to_anchor=(0.5, 0.0))

    fig.subplots_adjust(left=0.04, right=0.96, top=0.86, bottom=0.14)
    out_path = os.path.join(save_dir, f'{patient_id}_horizon_chain_{chain_number}_{lens_name}.png')
    fig.savefig(out_path, dpi=150, facecolor='white')
    plt.close(fig)
    return out_path


def run_all():
    os.makedirs(OBLIGATE_DIR, exist_ok=True)
    os.makedirs(PROBABLE_DIR, exist_ok=True)

    bp_df = pd.read_csv(BREAKPOINT_EDGES_CSV)
    completion_df = pd.read_csv(COMPLETION_CSV)

    # ALL enumerable chains, including the 83 with 0 true loose ends (real
    # already fully closed) — per explicit user instruction: even though
    # obligate/probable trivially equal real for those chains with 0
    # invented edges, every enumerable chain gets its own obligate and
    # probable image, not just the ones requiring invention.
    target_chains = completion_df[completion_df['is_enumerable']]

    n_drawn, n_failed = 0, 0
    for _, row in target_chains.iterrows():
        patient_id, chain_number = row['patient_id'], row['baca_chain_number']
        chain_bp = bp_df[(bp_df.patient_id == patient_id) & (bp_df.baca_chain_number == chain_number)]
        try:
            chrom_lookup = dict(zip(chain_bp['breakpoint_number'].astype(int), chain_bp['chromosome'].astype(int)))
            rearr_pairs = _real_rearr_pairs(chain_bp)
            ref_pairs = _real_ref_pairs(chain_bp)

            combined = _build_combined_graph(chain_bp)
            closed_components, open_components, loose_ends = _closed_and_open_components(combined)
            real_cycles = find_closed_cycles(chain_bp)

            enum_results = _enumerate_completions(closed_components, open_components, loose_ends, chrom_lookup)

            obligate_structure, _, _, obligate_achieving = _select_obligate(enum_results)
            (probable_structure, _, _, _, _, probable_achieving) = _select_probable(enum_results)

            obligate_matching = obligate_achieving[0]['matching']
            probable_matching = probable_achieving[0]['matching']

            obligate_structure_str = _cycle_structure_string(list(obligate_structure))
            probable_structure_str = _cycle_structure_string(list(probable_structure))

            out_obl = draw_completion(
                patient_id, chain_number, chain_bp, real_cycles, open_components,
                rearr_pairs, ref_pairs, obligate_matching, obligate_structure_str,
                'obligate', OBLIGATE_DIR,
            )
            out_prob = draw_completion(
                patient_id, chain_number, chain_bp, real_cycles, open_components,
                rearr_pairs, ref_pairs, probable_matching, probable_structure_str,
                'probable', PROBABLE_DIR,
            )
            n_drawn += 1
            print(f'{patient_id} chain {chain_number}: obligate={obligate_structure_str} probable={probable_structure_str}')
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    print(f'\nChains processed: {n_drawn}  Failed: {n_failed}')
    print(f'Obligate images -> {OBLIGATE_DIR}')
    print(f'Probable images -> {PROBABLE_DIR}')


if __name__ == '__main__':
    run_all()
