"""
horizon_step3_obligate_probable.py — HORIZON branch, Step 3.
COMPLETE REWRITE 2026-08-24 (superseding the first version, whose
"true loose end" definition was wrong — see CLAUDE.md's HORIZON section
for the full corrected write-up).

TRUE LOOSE END — locked definition, confirmed with the user:
  A bp-number-level open-path endpoint (e.g. BP6, found via the existing
  bp-number MultiGraph closure — same method as Step 2, unchanged) is
  NEVER itself the thing needing invented closure: it's a recorded end,
  always degree 2 (rearrangement + reference, per Step 1's locked rule).
  The true loose end is found by walking ONE HOP further, via that
  endpoint's OWN reference edge, to the actual open/shadow node it
  points at (e.g. BP6(+)'s reference edge points at BP18's un-recorded
  '-' side — "BP18(-)" is the true loose end, not BP6 itself) — UNLESS
  that reference edge points at TELOMERE (this chain's own boundary), in
  which case there is nothing to invent and this path tip contributes NO
  true loose end at all (explicitly excluded, not paired with anything).

  A true loose end is ALWAYS degree 1 by construction (asserted, not
  assumed): a bp-number-level path tip can only be degree-1 in the first
  place if its reference target is NOT an independently-recorded
  breakpoint (otherwise it would already be degree-2 in the bp-number
  graph and not a path tip at all) — so walking one hop always lands on a
  genuine shadow node, never a rare degree-2 collision case. Verified in
  code via an explicit assertion, per the user's instruction to make this
  a hard rule, not an assumption.

ODD COUNTS — per explicit user instruction: when a chain's number of true
loose ends is odd (a real, expected consequence of excluding telomere-
bound tips — not every open path necessarily contributes 2 valid
candidates), EVERY choice of which single end is left permanently open is
enumerated, combined with every perfect matching of the rest — the full
combined pool (not one arbitrary choice) is what obligate/probable
selection runs over.

CLOSURE PER CANDIDATE — computed by literally rebuilding the graph (real
rearrangement + real reference edges, restricted to the open portion of
the chain, PLUS this specific candidate's invented edges between true
loose ends) and decomposing it by degree, exactly as Step 2 already does
for the real-only case. This is NOT the same as the original Phase 2
script's fast union-find-over-components shortcut, which assumed a full
perfect matching always closes every open path completely — that
assumption breaks here, since a component can have only 1 (or 0) usable
true loose ends after telomere-exclusion, so "two components merge" no
longer implies "the merged thing closes." Rebuilding and checking degree
directly is slower but unambiguously correct.

Output: results/horizon_branch_res/horizon_csv_data/horizon_cycle_completion.csv
  — 1 row per chain (366 rows), same 30-column schema as before.
"""

import math
import os
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations

import networkx as nx
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'phase2'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'mmc5'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))

from horizon_step2_real_closure import _real_rearr_pairs, _real_ref_pairs
from baca_phase2_obligate_probable_completion import _select_obligate, _select_probable
from baca_mmc5_chain_closure import _cycle_structure_string
from baca_aberration_clinical_analysis import _all_perfect_matchings, _build_amg_nodes_and_edges

RESULTS_ROOT = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results/horizon_branch_res"
CSV_DIR = os.path.join(RESULTS_ROOT, "horizon_csv_data")

BREAKPOINT_EDGES_CSV = os.path.join(CSV_DIR, "horizon_breakpoint_edges.csv")
CHAIN_SUMMARY_CSV = os.path.join(CSV_DIR, "horizon_chain_structure_summary.csv")
REAL_CLOSURE_CSV = os.path.join(CSV_DIR, "horizon_real_closure_summary.csv")
OUTPUT_CSV = os.path.join(CSV_DIR, "horizon_cycle_completion.csv")

# Total candidate evaluations capped near the same order of magnitude as
# the original (m-1)!!<=10,395 precedent (m<=12). For odd m this is
# m * (m-2)!! (every leave-one-out choice x every matching of the rest);
# for even m it's (m-1)!! as before.
MAX_EVALUATIONS = 10_395
OTHER_STRAND = {'+': '-', '-': '+'}
_TELOMERE = 'TELOMERE'


def _closed_and_open_components(combined):
    """Same generic bp-number-level decomposition used everywhere else in
    this project (mirrors phase2's own version) — reproduced locally so
    this module has no dependency beyond what it actually reuses."""
    closed, open_, ends = [], [], []
    for comp in nx.connected_components(combined):
        degrees = {n: combined.degree(n) for n in comp}
        if all(d == 2 for d in degrees.values()):
            closed.append(frozenset(comp))
        else:
            endpoint_nodes = sorted(n for n, d in degrees.items() if d == 1)
            assert len(endpoint_nodes) == 2, f"Open component with {len(endpoint_nodes)} endpoints, expected 2: {comp}"
            assert all(d <= 2 for d in degrees.values()), f"Node with degree > 2: {comp}"
            open_.append(frozenset(comp))
            ends.extend(endpoint_nodes)
    return closed, open_, sorted(ends)


def _build_bp_number_graph(chain_bp):
    nodes = set(chain_bp['breakpoint_number'].astype(int))
    combined = nx.MultiGraph()
    combined.add_nodes_from(nodes)
    for a, b in _real_rearr_pairs(chain_bp):
        combined.add_edge(a, b)
    for a, b in _real_ref_pairs(chain_bp):
        combined.add_edge(a, b)
    return combined


def _build_rearr_rows_from_bp(chain_bp):
    """BP1/BP2-per-rearrangement rows, built directly from
    horizon_breakpoint_edges.csv (no mmc5 dependency) — same shape
    _build_amg_nodes_and_edges expects."""
    node_by_bp = {
        int(r['breakpoint_number']): (int(r['chromosome']), r['position'], r['strand'])
        for _, r in chain_bp.iterrows()
    }
    seen = set()
    out = []
    for _, r in chain_bp.iterrows():
        bp = int(r['breakpoint_number'])
        partner = int(r['rearrangement_partner_bp'])
        key = tuple(sorted((bp, partner)))
        if key in seen:
            continue
        seen.add(key)
        n1, n2 = node_by_bp[bp], node_by_bp[partner]
        out.append({
            'Breakpoint 1 chromosome': n1[0], 'Breakpoint 1 position': n1[1], 'Breakpoint 1 strand': n1[2],
            'Breakpoint 2 chromosome': n2[0], 'Breakpoint 2 position': n2[1], 'Breakpoint 2 strand': n2[2],
        })
    return pd.DataFrame(out)


def _resolve_true_loose_ends(chain_bp, path_tips, bp_to_node, ref_adj, node_to_bps):
    """
    For each bp-number path tip, walk its own reference edge one hop.
    Returns:
      tip_to_shadow: dict path_tip_bp -> (sibling_bp, shadow_strand) | None (telomere-excluded)
    """
    tip_to_shadow = {}
    for tip in path_tips:
        target = ref_adj[bp_to_node[tip]]
        if target == _TELOMERE:
            tip_to_shadow[tip] = None
            continue
        # target is a (chrom,pos,strand) node. Must be a genuine shadow —
        # never independently recorded — since a path tip can only be
        # degree-1 in the bp-number graph if its reference target ISN'T a
        # real recorded breakpoint (else it'd already be degree-2 there).
        assert target not in node_to_bps, (
            f"Path tip {tip}'s reference target {target} is independently "
            f"recorded ({node_to_bps[target]}) — should be impossible for a "
            f"bp-number-level degree-1 endpoint."
        )
        opp_strand = OTHER_STRAND[target[2]]
        recorded_side = (target[0], target[1], opp_strand)
        assert recorded_side in node_to_bps, (
            f"Path tip {tip}'s reference target {target} has no recorded "
            f"sibling at {recorded_side} — should be impossible."
        )
        sibling_bp = min(node_to_bps[recorded_side])
        tip_to_shadow[tip] = (sibling_bp, target[2])
    return tip_to_shadow


def _build_open_portion_base_graph(open_components, rearr_pairs, ref_pairs, tip_to_shadow):
    """
    Graph over bp-numbers (all members of open_components) PLUS shadow
    true-loose-end identities (one per enumerable tip) — edges are ALL
    real: rearrangement, real-real reference, and each enumerable tip's
    own real reference edge out to its shadow true loose end. No invented
    edges here — those are added per-candidate on a copy of this graph.
    """
    all_bp_nodes = set(n for comp in open_components for n in comp)
    G = nx.MultiGraph()
    G.add_nodes_from(all_bp_nodes)
    for a, b in rearr_pairs:
        if a in all_bp_nodes and b in all_bp_nodes:
            G.add_edge(a, b)
    for a, b in ref_pairs:
        if a in all_bp_nodes and b in all_bp_nodes:
            G.add_edge(a, b)
    for tip, shadow in tip_to_shadow.items():
        if shadow is not None:
            G.add_node(shadow)
            G.add_edge(tip, shadow)
    return G


def _bp_identity(node):
    """A bp-number node's identity is itself; a shadow true-loose-end
    node (the (bp, strand) tuple its own reference edge points at) is
    identified by that SAME bp number — its other side. Cn counts
    DISTINCT BREAKPOINTS in the closed structure, per explicit user
    instruction: a breakpoint whose two sides both end up swept into one
    cycle is counted ONCE, not twice, even though the graph itself (for
    degree-checking purposes) still tracks its two sides as separate
    nodes."""
    return node[0] if isinstance(node, tuple) else node


def _evaluate_candidate(base_graph, invented_pairs, chrom_lookup):
    G = base_graph.copy()
    for a, b in invented_pairs:
        G.add_edge(a, b)
    new_lengths, new_spans = [], []
    for comp in nx.connected_components(G):
        degrees = [G.degree(n) for n in comp]
        if len(comp) >= 2 and all(d == 2 for d in degrees):
            distinct_bps = {_bp_identity(n) for n in comp}
            new_lengths.append(len(distinct_bps))
            new_spans.append(len({chrom_lookup[n] for n in comp}))
    return new_lengths, new_spans


def _all_partial_matchings(elements, k):
    """
    Every way to choose exactly k disjoint pairs from `elements` — the
    remaining len(elements) - 2k elements are simply left untouched, not
    assigned to anything. This is the general form; a "full" matching is
    just the special case k = len(elements)//2.
    """
    if k == 0:
        yield []
        return
    n = len(elements)
    for chosen_idx in combinations(range(n), 2 * k):
        chosen = [elements[i] for i in chosen_idx]
        yield from _all_perfect_matchings(chosen)


def _enumerate_true_completions(base_graph, enumerable_shadows, chrom_lookup):
    """
    Returns list of dicts {structure, rule2_metric, max_span, matching}
    — same shape _select_obligate/_select_probable already expect
    (reused unchanged). `matching` here is a list of (shadow_a, shadow_b)
    true-loose-end pairs — NOT necessarily using every true loose end.

    LOCKED 2026-08-24 (explicit user instruction, corrected same session):
    obligate/probable must NOT be forced into a full perfect matching of
    every true loose end. Cornforth's own philosophy favors inventing as
    little as possible — if ONE invented edge already closes a cycle on
    its own, that is a valid, often more conservative, candidate than any
    solution requiring more invented edges elsewhere. So this enumerates
    EVERY invented-edge count from 1 up to floor(m/2): for each count k,
    every way to choose k disjoint pairs among the m true loose ends,
    leaving the rest of them simply untouched (not assigned to anything,
    not "excluded" as a special case — just not part of this particular
    candidate). This subsumes the earlier odd-count "leave-one-out"
    handling as the special case k = m//2 — no separate odd/even logic is
    needed here at all now.

    Real, obligate, and probable are three SEPARATE computations, never
    combined into one structure string. This function describes ONLY what
    a given invented completion itself produces — the chain's already-
    real closed cycles are NOT folded in here; they're reported
    separately as real_cycle_structure, computed once, unaffected by any
    of this.
    """
    results = []

    def _record(matching):
        new_lengths, new_spans = _evaluate_candidate(base_graph, matching, chrom_lookup)
        structure = tuple(sorted(new_lengths, reverse=True))
        rule2_metric = max(new_spans) if new_spans else 0
        max_span = max(new_spans) if new_spans else 0
        results.append({'structure': structure, 'rule2_metric': rule2_metric, 'max_span': max_span, 'matching': matching})

    m = len(enumerable_shadows)
    if m < 2:
        _record([])
        return results

    for k in range(1, m // 2 + 1):
        for matching in _all_partial_matchings(enumerable_shadows, k):
            _record(matching)

    return results


def compute_chain_completion(chain_bp):
    chrom_lookup = dict(zip(chain_bp['breakpoint_number'].astype(int), chain_bp['chromosome'].astype(int)))
    strand_lookup = dict(zip(chain_bp['breakpoint_number'].astype(int), chain_bp['strand']))

    bp_to_node = {
        int(r['breakpoint_number']): (int(r['chromosome']), r['position'], r['strand'])
        for _, r in chain_bp.iterrows()
    }
    node_to_bps = defaultdict(list)
    for bp, node in bp_to_node.items():
        node_to_bps[node].append(bp)

    rearr_rows = _build_rearr_rows_from_bp(chain_bp)
    all_nodes, rearr_adj_node, ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(rearr_rows)

    bp_graph = _build_bp_number_graph(chain_bp)
    closed_components, open_components, path_tips = _closed_and_open_components(bp_graph)

    real_structure_tuple = tuple(sorted((len(c) for c in closed_components), reverse=True))
    real_cycle_structure = _cycle_structure_string(list(real_structure_tuple))

    tip_to_shadow = _resolve_true_loose_ends(chain_bp, path_tips, bp_to_node, ref_adj, node_to_bps)
    enumerable_shadows = sorted(s for s in tip_to_shadow.values() if s is not None)
    m = len(enumerable_shadows)

    def _shadow_label(shadow):
        bp, strand = shadow
        return f'BP{bp}({strand})'
    true_loose_ends_bp_strand_list = ", ".join(_shadow_label(s) for s in enumerable_shadows)

    base = {
        'real_has_any_closed_cycle': len(closed_components) > 0,
        'real_cycle_structure': real_cycle_structure,
        'real_n_cycles': len(closed_components),
        'n_true_loose_ends': m,
        'true_loose_ends_bp_strand_list': true_loose_ends_bp_strand_list,
        'n_invented_edges': m // 2,  # MAXIMUM possible (a full pairing) — the
                                     # obligate/probable winners can each use
                                     # FEWER, see obligate_/probable_n_invented_edges
    }

    # Total candidate count across every invented-edge count k=1..floor(m/2):
    # sum of C(m,2k) * (2k-1)!! — generalizes the old (m-1)!!-only formula
    # now that partial (not just full) matchings are explored.
    def _dfact(n):
        r = 1
        while n > 1:
            r *= n
            n -= 2
        return r

    n_evaluations = sum(
        math.comb(m, 2 * k) * _dfact(2 * k - 1) for k in range(1, m // 2 + 1)
    ) if m >= 2 else 1

    if n_evaluations > MAX_EVALUATIONS:
        base.update({
            'is_enumerable': False,
            'n_matchings_enumerated': None,
            'obligate_cycle_structure': None, 'obligate_n_cycles': None,
            'obligate_max_cycle_chrom_span': None,
            'obligate_n_matchings_achieving_structure': None, 'obligate_n_invented_edges': None,
            'probable_cycle_structure': None, 'probable_n_cycles': None,
            'probable_max_cycle_chrom_span': None,
            'probable_n_matchings_achieving_structure': None,
            'probable_all_structures_enumerated': None,
            'probable_is_tied': None, 'probable_tied_alternatives': None, 'probable_n_invented_edges': None,
        })
        return base

    rearr_pairs = _real_rearr_pairs(chain_bp)
    ref_pairs = _real_ref_pairs(chain_bp)
    base_graph = _build_open_portion_base_graph(open_components, rearr_pairs, ref_pairs, tip_to_shadow)
    chrom_lookup_extended = dict(chrom_lookup)
    for shadow in enumerable_shadows:
        chrom_lookup_extended[shadow] = chrom_lookup[shadow[0]]

    enum_results = _enumerate_true_completions(base_graph, enumerable_shadows, chrom_lookup_extended)
    total_matchings = len(enum_results)

    obligate_structure, obligate_rule2, n_at_min, obligate_achieving = _select_obligate(enum_results)
    (probable_structure, probable_count, n_distinct,
     probable_is_tied, tied_alts, probable_achieving) = _select_probable(enum_results)

    obligate_max_span = min(r['max_span'] for r in obligate_achieving)
    probable_max_span = min(r['max_span'] for r in probable_achieving)
    # Among tied candidates achieving the winning structure, report the
    # FEWEST invented edges used to get there — consistent with obligate's
    # own least-coordination spirit, and simply the most informative single
    # number to show for probable too (several matchings can tie on
    # frequency using different edge counts).
    obligate_n_invented_edges = min(len(r['matching']) for r in obligate_achieving)
    probable_n_invented_edges = min(len(r['matching']) for r in probable_achieving)

    freq = Counter(r['structure'] for r in enum_results)
    breakdown_sorted = sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))
    probable_all_structures_enumerated = ", ".join(
        f'{_cycle_structure_string(list(s))}: {c}/{total_matchings}' for s, c in breakdown_sorted
    )

    base.update({
        'is_enumerable': True,
        'n_matchings_enumerated': total_matchings,

        'obligate_cycle_structure': _cycle_structure_string(list(obligate_structure)),
        'obligate_n_cycles': len(obligate_structure),
        'obligate_max_cycle_chrom_span': obligate_max_span,
        'obligate_n_matchings_achieving_structure': len(obligate_achieving),
        'obligate_n_invented_edges': obligate_n_invented_edges,

        'probable_cycle_structure': _cycle_structure_string(list(probable_structure)),
        'probable_n_cycles': len(probable_structure),
        'probable_max_cycle_chrom_span': probable_max_span,
        'probable_n_matchings_achieving_structure': len(probable_achieving),
        'probable_n_invented_edges': probable_n_invented_edges,
        'probable_all_structures_enumerated': probable_all_structures_enumerated,
        'probable_is_tied': probable_is_tied,
        'probable_tied_alternatives': "; ".join(
            _cycle_structure_string(list(s)) for s in tied_alts
        ) if tied_alts else "",
    })
    return base


COLUMN_ORDER = [
    'patient_id', 'baca_chain_number',
    'n_breakpoints', 'n_rearrangements', 'k_chromosomes', 'chromosomes',

    'real_has_any_closed_cycle', 'real_cycle_structure', 'real_n_cycles',

    'n_true_loose_ends', 'true_loose_ends_bp_strand_list', 'n_invented_edges',
    'is_enumerable', 'n_matchings_enumerated',

    'obligate_cycle_structure', 'obligate_n_cycles', 'obligate_max_cycle_chrom_span',
    'obligate_n_matchings_achieving_structure', 'obligate_n_invented_edges',

    'probable_cycle_structure', 'probable_n_cycles', 'probable_max_cycle_chrom_span',
    'probable_n_matchings_achieving_structure', 'probable_n_invented_edges',
    'probable_all_structures_enumerated', 'probable_is_tied', 'probable_tied_alternatives',

    'contains_ERG_or_TMPRSS2', 'genes_in_chain', 'ETS_status', 'Gleason_Score', 'Pathological_stage',
]


def run_all():
    bp_df = pd.read_csv(BREAKPOINT_EDGES_CSV)
    chain_summary_df = pd.read_csv(CHAIN_SUMMARY_CSV)
    real_closure_df = pd.read_csv(REAL_CLOSURE_CSV)

    clinical_cols = chain_summary_df.set_index(['patient_id', 'baca_chain_number'])[
        ['n_breakpoints', 'n_rearrangements', 'k_chromosomes', 'chromosomes',
         'contains_ERG_or_TMPRSS2', 'genes_in_chain', 'ETS_status', 'Gleason_Score', 'Pathological_stage']
    ]
    published_real = real_closure_df.set_index(['patient_id', 'baca_chain_number'])[
        ['cycle_structure', 'has_any_closed_cycle', 'n_cycles']
    ]

    rows = []
    n_failed = 0
    for (patient_id, chain_number), chain_bp in bp_df.groupby(['patient_id', 'baca_chain_number']):
        try:
            result = compute_chain_completion(chain_bp)

            pub = published_real.loc[(patient_id, chain_number)]
            assert result['real_cycle_structure'] == pub['cycle_structure'], (
                f"Real-structure mismatch for {patient_id} chain {chain_number}: "
                f"{result['real_cycle_structure']!r} != {pub['cycle_structure']!r}"
            )
            assert result['real_has_any_closed_cycle'] == bool(pub['has_any_closed_cycle'])
            assert result['real_n_cycles'] == pub['n_cycles']

            clin = clinical_cols.loc[(patient_id, chain_number)]
            row = {'patient_id': patient_id, 'baca_chain_number': chain_number}
            row.update(clin.to_dict())
            row.update(result)
            rows.append(row)
        except Exception as e:
            n_failed += 1
            print(f'FAILED {patient_id} chain {chain_number}: {e}')

    df = pd.DataFrame(rows).sort_values(['patient_id', 'baca_chain_number']).reset_index(drop=True)
    df = df[COLUMN_ORDER]
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {len(df)} chains ({df['patient_id'].nunique()} patients) -> {OUTPUT_CSV}  Failed: {n_failed}")

    print()
    print("=== Coverage ===")
    print(f"Chains with 0 true loose ends: {(df['n_true_loose_ends']==0).sum()} / {len(df)}")
    print(f"Chains enumerable: {df['is_enumerable'].sum()} / {len(df)}")
    not_enum = df[~df['is_enumerable']]
    print(f"Chains NOT enumerable: {len(not_enum)} -- "
          f"{list(zip(not_enum['patient_id'], not_enum['baca_chain_number'], not_enum['n_true_loose_ends']))}")

    print()
    print("=== Obligate vs probable, among enumerable chains that needed invention ===")
    enum_df = df[df['is_enumerable']]
    invented_df = enum_df[enum_df['n_invented_edges'] > 0]
    same = (invented_df['obligate_cycle_structure'] == invented_df['probable_cycle_structure']).sum()
    print(f"Of {len(invented_df)} chains requiring invention: obligate==probable for {same}, differ for {len(invented_df)-same}")
    print(f"Chains where probable_is_tied: {invented_df['probable_is_tied'].sum()} / {len(invented_df)}")


if __name__ == '__main__':
    run_all()
