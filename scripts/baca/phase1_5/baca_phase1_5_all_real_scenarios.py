"""
baca_phase1_5_all_real_scenarios.py — enumerate ALL real-data-only cycle
structure scenarios per patient under the NEW (mmc5-merged) reference-edge
model from baca_phase1_5_merged_closure.py.

Nothing invented anywhere in this script. Every candidate reference edge
considered — for both the forced ones and the ones inside a genuinely
ambiguous cluster — is either (a) our position-order fallback rule applied
to real recorded breakpoint positions, or (b) Baca's own real mmc5
adjacency/deletion-bridge edge. "Ambiguity" here means Baca's own data
gives MORE THAN ONE real, already-recorded candidate for the same
breakpoint (e.g. an mmc5 adjacency partner that disagrees with a
deletion-bridge partner, or with our fallback) — we don't know which one
the real biology used, so we enumerate every real option rather than
picking one arbitrarily.

Why this is feasible (verified empirically before writing this script):
52/57 patients have at least one node with multiple real reference-edge
candidates, which sounds like it could blow up combinatorially. It
doesn't. Repeatedly forcing every node that has only ONE remaining
candidate (which cascades — fixing one edge can force its neighbor's
neighbor, and so on) resolves almost the entire graph uniquely. What's
left, genome-wide across all 57 patients, is just 3 small residual
clusters where a genuine choice remains (two of 4 nodes, one of 6 nodes).
Enumerating every real combination across this handful of clusters is
cheap — at most a few dozen total combinations cohort-wide, not 2^87.

Method per patient:
  1. Build the same real rearrangement-edge matching + candidate
     reference-edge graph as baca_phase1_5_merged_closure.py (fallback
     always present, mmc5 adjacency/deletion-bridge added wherever that
     breakpoint is mmc5-covered).
  2. Repeatedly force any node with exactly one remaining candidate
     (fix that edge, remove both endpoints) until no more forced nodes
     remain.
  3. Whatever residual connected components remain (size >= 2 = genuine
     choice) get every maximum matching enumerated exhaustively — brute
     force, since components are tiny.
  4. For every combination across however many ambiguous components a
     patient has (Cartesian product), combine with the forced edges +
     real rearrangement matching, and classify closed/open exactly like
     baca_phase1_5_merged_closure.py (degree-based, nx.MultiGraph).
  5. One of these combinations, by construction, must reproduce EXACTLY
     the new_* result already in results/full_genome_cycle_phase1_5_comparison.csv
     (networkx's max_weight_matching always picks ONE valid maximum
     matching, which must appear somewhere in our enumeration of ALL
     valid maximum matchings) — flagged via is_original_new_result and
     verified against that file before trusting this script's output.

Output:
  results/full_genome_cycle_all_real_scenarios.csv — one row per
  (patient, scenario). 54/57 patients have exactly one real scenario (no
  ambiguity); the remaining ones get one row per distinct real
  combination found.
"""

import itertools
import os
import re
import sys
from collections import Counter

import networkx as nx
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'phase1'))
from baca_aberration_clinical_analysis import (
    load_data,
    _build_amg_nodes_and_edges,
    load_ets_status_from_mmc5,
    BACA_CHROM_ABER_FILE_PATH,
    BACA_CLINICAL_PHENO_FILE_PATH,
)
from baca_full_genome_cycles import _get_full_patient_rows
from baca_phase1_5_merged_closure import (
    load_mmc5_reference_edges,
    _cycle_structure_string,
)

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
ORIGINAL_CSV = os.path.join(RESULTS_DIR, "full_genome_cycle_phase1_5_comparison.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "full_genome_cycle_all_real_scenarios.csv")

ERG_PATTERN = re.compile(r'ERG|TMPRSS2')


# ── Force-propagation + small-component matching enumeration ───────────────

def _force_propagate(ref_graph):
    """
    Repeatedly fix any node with exactly one remaining real candidate edge,
    removing both endpoints, until no more such nodes exist. Returns
    (forced_edges, residual_graph). Residual graph nodes have degree 0
    (no real candidate left — a real open end) or degree >= 2 (genuine
    ambiguity — multiple real candidates that survive mutual constraints).
    """
    G = ref_graph.copy()
    forced_edges = []
    changed = True
    while changed:
        changed = False
        deg1_nodes = sorted([n for n in G.nodes() if G.degree(n) == 1])
        for n in deg1_nodes:
            if n not in G or G.degree(n) != 1:
                continue
            partner = sorted(G.neighbors(n))[0]
            forced_edges.append((n, partner))
            G.remove_node(n)
            if partner in G:
                G.remove_node(partner)
            changed = True
    return forced_edges, G


def _all_maximum_matchings(nodes, adj):
    """
    Brute-force enumeration of every matching of the GIVEN real edges
    (adj: node -> set of real neighbors within `nodes`) that achieves the
    maximum possible size. Only used on tiny residual components (<=6
    nodes cohort-wide, verified) — exhaustive search is cheap here.
    """
    all_matchings = []

    def backtrack(remaining, current):
        if not remaining:
            all_matchings.append(tuple(sorted(current)))
            return
        node = remaining[0]
        rest = remaining[1:]
        matched_here = False
        for nb in sorted(adj.get(node, [])):
            if nb in rest:
                new_rest = [n for n in rest if n != nb]
                backtrack(new_rest, current + [tuple(sorted((node, nb), key=str))])
                matched_here = True
        backtrack(rest, current)  # also try leaving `node` unmatched

    backtrack(sorted(nodes), [])
    max_size = max((len(m) for m in all_matchings), default=0)
    unique_max = sorted(set(m for m in all_matchings if len(m) == max_size))
    return unique_max


# ── Per-patient scenario enumeration ────────────────────────────────────────

def _chrom_set(chromosomes_str):
    return set(c for c in chromosomes_str.split(',') if c)


def _classify_combined(all_nodes, active_nodes, rearr_adj, full_ref_matching):
    """Same degree-based closure classification as baca_phase1_5_merged_closure.py."""
    combined = nx.MultiGraph()
    combined.add_nodes_from(sorted(all_nodes))
    seen_pairs = set()
    for u, v in rearr_adj.items():
        key = tuple(sorted((u, v), key=str))
        if key not in seen_pairs:
            combined.add_edge(u, v)
            seen_pairs.add(key)
    for u, v in full_ref_matching:
        combined.add_edge(u, v)

    cycle_lengths = []
    cycle_chrom_spans = []
    n_open_active_breakpoints = 0
    closed_cycle_dsb_positions = set()

    for component in nx.connected_components(combined):
        active_in_component = component & active_nodes
        if not active_in_component:
            continue
        degrees = [combined.degree(n) for n in component]
        is_closed = len(component) >= 2 and all(d == 2 for d in degrees)
        if is_closed:
            cycle_lengths.append(len(component))
            cycle_chrom_spans.append(len({n[0] for n in component}))
            closed_cycle_dsb_positions.update((n[0], n[1]) for n in component)
        else:
            n_open_active_breakpoints += len(active_in_component)

    b_i_cycles = Counter(chrom for chrom, _ in closed_cycle_dsb_positions)
    k_cycles = len(b_i_cycles)
    chromosomes_in_cycles = sorted(b_i_cycles.keys())
    theta_cycles = (
        f"Theta({k_cycles},({','.join(str(b_i_cycles[c]) for c in chromosomes_in_cycles)}))"
        if k_cycles else "Theta(0,())"
    )
    max_cycle_chrom_span = max(cycle_chrom_spans) if cycle_chrom_spans else 0

    return {
        'cycle_size_histogram': _cycle_structure_string(cycle_lengths),
        'n_closed_cycles': len(cycle_lengths),
        'rearr_edges_in_cycles': sum(cycle_lengths) // 2,
        'n_open_breakpoints': n_open_active_breakpoints,
        'theta_cycles': theta_cycles,
        'chromosomes_in_cycles': ",".join(str(int(c)) for c in chromosomes_in_cycles),
        'max_cycle_chrom_span': max_cycle_chrom_span,
        'has_chromoplexy_strict': max_cycle_chrom_span >= 3,
        'has_chromoplexy_loose': max_cycle_chrom_span >= 2,
    }


def enumerate_patient_scenarios(patient_id, chrom_df, mmc5_node_partners):
    patient_rows = _get_full_patient_rows(patient_id, chrom_df)
    if len(patient_rows) == 0:
        return []

    n_rejoins = len(patient_rows)
    n_dsbs = n_rejoins * 2
    all_nodes, rearr_adj, fallback_ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(patient_rows)
    k = len(chr_to_dsb_pos)
    max_b_i = max((len(p) for p in chr_to_dsb_pos.values()), default=0)
    active_nodes = set(rearr_adj.keys())

    ref_graph = nx.Graph()
    ref_graph.add_nodes_from(sorted(all_nodes))
    for u, v in fallback_ref_adj.items():
        if v != 'TELOMERE':
            ref_graph.add_edge(u, v)

    node_partners = mmc5_node_partners.get(patient_id, {})
    n_mmc5_covered = 0
    for node in sorted(active_nodes):
        if node in node_partners:
            n_mmc5_covered += 1
            for partner in sorted(node_partners[node]):
                if partner in all_nodes:
                    ref_graph.add_edge(node, partner)

    # Force-propagate on the FULL graph (active + passive nodes) — matching
    # baca_phase1_5_merged_closure.py exactly. Passive nodes always have
    # degree <= 1 here (mmc5 never covers them; they only ever get a single
    # fallback edge), so they always resolve in the first propagation pass
    # and a valid active-node fallback edge that happens to target a
    # passive node is preserved correctly. Restricting to active_nodes
    # BEFORE propagation (an earlier, buggy version of this function) wrongly
    # dropped exactly those edges, breaking the original-result cross-check
    # below — caught and fixed before trusting any output.
    forced_edges, residual = _force_propagate(ref_graph)

    ambiguous_components = [c for c in nx.connected_components(residual) if len(c) >= 2]
    component_matchings = []
    for comp in ambiguous_components:
        adj = {n: set(residual.neighbors(n)) & comp for n in comp}
        matchings = _all_maximum_matchings(sorted(comp), adj)
        component_matchings.append(matchings)

    n_scenarios = 1
    for m in component_matchings:
        n_scenarios *= max(len(m), 1)

    rows = []
    combos = itertools.product(*component_matchings) if component_matchings else [()]
    for scenario_idx, combo in enumerate(combos, start=1):
        full_ref_matching = list(forced_edges)
        for matching in combo:
            full_ref_matching.extend(matching)

        stats = _classify_combined(all_nodes, active_nodes, rearr_adj, full_ref_matching)
        frac_rejoins_in_cycles = stats['rearr_edges_in_cycles'] / n_rejoins if n_rejoins else 0.0

        rows.append({
            'patient_id': patient_id,
            'scenario_id': scenario_idx,
            'n_scenarios_total': n_scenarios,
            'n_ambiguous_components': len(ambiguous_components),
            'ambiguous_component_sizes': ",".join(str(len(c)) for c in ambiguous_components),
            'k': k,
            'n_rejoins': n_rejoins,
            'n_dsbs': n_dsbs,
            'max_b_i': max_b_i,
            'n_mmc5_covered_breakpoints': n_mmc5_covered,
            'n_closed_cycles': stats['n_closed_cycles'],
            'n_open_breakpoints': stats['n_open_breakpoints'],
            'frac_rejoins_in_cycles': round(frac_rejoins_in_cycles, 3),
            'cycle_size_histogram': stats['cycle_size_histogram'],
            'theta_cycles': stats['theta_cycles'],
            'chromosomes_in_cycles': stats['chromosomes_in_cycles'],
            'max_cycle_chrom_span': stats['max_cycle_chrom_span'],
            'has_chromoplexy_strict': stats['has_chromoplexy_strict'],
            'has_chromoplexy_loose': stats['has_chromoplexy_loose'],
        })

    return rows


# ── MAIN ─────────────────────────────────────────────────────────────────

def main():
    chrom_df, clinical_df = load_data(BACA_CHROM_ABER_FILE_PATH, BACA_CLINICAL_PHENO_FILE_PATH)
    mmc5_node_partners = load_mmc5_reference_edges()
    original_df = pd.read_csv(ORIGINAL_CSV)
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])

    all_rows = []
    for patient_id in sorted(chrom_df['Individual'].dropna().unique(), key=str):
        scenarios = enumerate_patient_scenarios(patient_id, chrom_df, mmc5_node_partners)
        if not scenarios:
            continue

        patient_rows = _get_full_patient_rows(patient_id, chrom_df)
        contains_erg = bool(
            patient_rows['Breakpoint 1 site'].dropna().str.contains(ERG_PATTERN).any() or
            patient_rows['Breakpoint 2 site'].dropna().str.contains(ERG_PATTERN).any() or
            patient_rows['Fusion'].dropna().str.contains(ERG_PATTERN).any()
        )
        ets_status = ets_lookup.get(patient_id)

        orig_row = original_df[original_df['patient_id'] == patient_id]
        orig_histogram = orig_row['new_cycle_size_histogram'].values[0] if len(orig_row) > 0 else None
        orig_theta = orig_row['new_theta_cycles'].values[0] if len(orig_row) > 0 else None
        orig_chroms = orig_row['new_chromosomes_in_cycles'].values[0] if len(orig_row) > 0 else None

        for row in scenarios:
            row['contains_ERG_or_TMPRSS2'] = contains_erg
            row['ETS_status'] = ets_status
            row['is_original_new_result'] = (
                row['cycle_size_histogram'] == orig_histogram and
                row['theta_cycles'] == orig_theta and
                row['chromosomes_in_cycles'] == (orig_chroms if pd.notna(orig_chroms) else '')
            )
            all_rows.append(row)

    out_df = pd.DataFrame(all_rows)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {len(out_df)} scenario-rows across {out_df['patient_id'].nunique()} patients -> {OUTPUT_CSV}")

    print()
    print(f"Patients with exactly 1 real scenario (no ambiguity): {(out_df.groupby('patient_id')['n_scenarios_total'].first() == 1).sum()} / {out_df['patient_id'].nunique()}")
    multi = out_df.groupby('patient_id')['n_scenarios_total'].first()
    multi = multi[multi > 1]
    print(f"Patients with >1 real scenario: {len(multi)}")
    for pid, n in multi.items():
        print(f"  {pid}: {n} scenarios")

    print()
    found_flags = out_df.groupby('patient_id')['is_original_new_result'].any()
    print(f"Original new_* result found among enumerated scenarios: {found_flags.sum()} / {len(found_flags)} patients")
    if found_flags.sum() < len(found_flags):
        print("  MISSING for:", list(found_flags[~found_flags].index))


if __name__ == '__main__':
    main()
