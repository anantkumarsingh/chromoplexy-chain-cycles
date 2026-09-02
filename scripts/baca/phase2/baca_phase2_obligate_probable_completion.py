"""
baca_phase2_obligate_probable_completion.py — Phase 2: obligate and
most-probable cycle-structure completion for all 366 Baca mmc5 chains.

Full methodology, terminology, and the worked Cornforth 2001 examples this
implementation is built from are documented in
PHASE2_OBLIGATE_AND_PROBABLE_METHODOLOGY.md (project root). Summary:

  - Cornforth 2001 ("obligate cycle structure", Sheth's ref [16]) defines the
    obligate structure as whichever achievable cycle structure minimizes the
    largest cycle order (ties broken by the next-largest, etc.) -- i.e. the
    most conservative, least-coordination-assumed reading.
  - "Most probable" is simply whichever structure occurs most frequently
    among every combinatorially valid way of closing what's still open.
  - Cornforth's two-tiered pattern-closure rule (fewest new connections,
    then fewest additional chromosomes) is what decides WHICH open ends get
    paired before the obligate rule is even applied -- this directly
    answers the multi-open-path tie-break question this project had
    previously logged as unresolved.

WHAT THIS SCRIPT DOES, IN ONE SENTENCE: for every one of the 366 mmc5
chains, it starts from the exact same real (rearrangement + adjacency +
deletion-bridge) graph already validated in baca_mmc5_chain_closure.py,
finds whatever open ends remain (Cornforth's "one-way" ends), and computes
both the obligate and most-probable ways of closing them -- reporting both
as separate, clearly-labeled columns, never collapsed into one answer.

SCOPE DECISIONS, recorded explicitly to avoid any computational ambiguity:

  1. ALL 366 chains are processed, including the 91 that already have a
     real closed cycle and the 52 that are already fully closed (m=0) --
     per 2026-06-24 decision, to see what obligate/probable framing says
     even where real data already closes something. For m=0 chains,
     obligate == probable == the existing real structure, trivially.

  2. NO telomere/centromere distinction is applied to open ends in this
     version. Cornforth's original methodology distinguishes a true
     chromosome terminus (not missing data) from a genuinely unrecorded
     partner ("one-way" in the strict sense) -- but that requires
     centromere/telomere reference coordinates, which this project has
     explicitly deferred to a later, separate biological-constraints phase
     (see CENTROMERE_PQ_ARM_NOTES.md / CLAUDE.md). Every degree-1 node in
     the real combined graph is therefore treated as a closure candidate.
     This is a documented scope choice, not an oversight, and should be
     revisited if/when the centromere phase begins.

  3. Never invent a new breakpoint. Only ever connect two breakpoints that
     already exist in the real mmc5 data (Cornforth's explicit rejection of
     cryptic/invented breaks, see methodology doc Part 1.2).

  4. Cornforth's Rule 2 ("fewest additional chromosomes") is adapted, not
     applied literally: in his original radiation context, an unknown
     missing FRAGMENT could in principle be assumed to belong to any
     chromosome, so "fewest chromosomes" had a free choice to disambiguate.
     Here, every open end is an already-known, already-located real
     breakpoint -- pairing it with another known breakpoint never
     "introduces" a new chromosome to the dataset, so a literal reading of
     Rule 2 never discriminates between candidate pairings (verified: the
     UNION of chromosomes touched by all the open ends combined is fixed
     regardless of pairing choice). The metric actually used, OUR OWN
     ADAPTATION of the same conservative spirit, is: for each candidate
     pairing, take the MAXIMUM chromosome span among only the NEWLY-closed
     cycles it creates, and minimize that. This must be a MAX, not a SUM --
     an earlier version of this script summed chromosome spans across new
     cycles and was caught (by hand-tracing P01-28 chain 1 before trusting
     any aggregate result) to be systematically BIASED TOWARD merging:
     summing double-counts any chromosome shared by two open paths, so
     span(A)+span(B) > span(A union B) whenever A and B share a chromosome,
     making the merged (more coordinated) option look artificially cheaper
     -- exactly backwards from the conservative reading Cornforth's rule is
     supposed to produce. Using max() instead removes this bias entirely,
     since span(A union B) >= max(span(A), span(B)) always holds, so merging
     can never look cheaper than keeping paths separate. This is flagged
     here, and in the methodology doc, explicitly as our extension -- not a
     literal transcription of Cornforth's rule.

  5. Tractability cap: full enumeration of every way to pair up `m` open
     ends is (m-1)!! matchings. Capped at m <= 12 (<=10,395 matchings),
     the exact same numeric precedent already established for small-ERG-
     chain enumeration elsewhere in this codebase (MAX_REJOINS=6 ->
     n_DSBs<=12). Verified empirically before writing this script: 354/366
     chains have m <= 12; only 12 chains exceed it (m up to 24). Note this
     is a MUCH smaller fraction than earlier feared, because `m` here counts
     only the open ENDS needing invented closure, not the whole chain's
     breakpoint count -- most large chains are already mostly closed by
     real data, leaving comparatively few genuinely open ends.

CORRECTNESS, cross-checked before trusting any result:
  - Verified empirically (see scratch analysis run before writing this
    script) that across all 366 chains: m is always even, and no
    breakpoint ever has degree 0 in the real combined graph. Both are
    required for the algorithm below to be well-defined, and both were
    confirmed, not assumed.
  - Every node has degree <= 2 in the real combined graph (rearrangement
    contributes exactly 1 edge per node; the reference-layer matching
    contributes at most 1 more). A graph with max degree 2 decomposes only
    into disjoint simple paths and simple cycles -- so every non-closed
    component has EXACTLY 2 degree-1 endpoint nodes, never more. This is
    asserted in code, not just assumed.
  - Obligate selection uses Python's native tuple comparison over each
    candidate's sorted-descending cycle-length tuple. For partitions of a
    FIXED total (which is always true here -- every candidate closes the
    exact same set of open ends), tuple comparison is mathematically
    IDENTICAL to Cornforth's stated rule ("minimize the largest order, tie-
    break on the next-largest, etc.") -- this equivalence is explained in
    the inline comment on `_structure_tuple` below, not just asserted.
  - Reuses the already-validated `nx.MultiGraph()` construction and
    integer-node (not string-tuple) graph primitives from
    baca_mmc5_chain_closure.py verbatim -- both the MultiGraph
    requirement and the no-PYTHONHASHSEED-sensitivity property documented
    there carry over unchanged, since this script builds on the exact same
    primitives rather than reimplementing them.
  - Every per-chain result independently re-derives the REAL (pre-Phase-2)
    cycle structure via `compute_chain_cycle_structure`
    (baca_mmc5_primary_chain_cycle_structures.py, already published in
    results/mmc5_primary_chain_cycle_structures.csv) and asserts it
    matches this script's own from-scratch real-graph decomposition
    exactly, for every single chain, before adding any obligate/probable
    columns -- this catches any divergence between the two
    implementations immediately rather than silently trusting one.

Output: results/phase2_obligate_probable_completion.csv -- one row per
(patient, Baca chain number), all 366 chains, all 57 mmc5-covered patients.
"""

import os
import sys
from collections import Counter, defaultdict

import networkx as nx
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'mmc5'))
from baca_mmc5_chain_closure import (
    load_mmc5_table_s5a,
    extract_genes,
    _build_rearr_matching,
    _build_ref_graph,
    _cycle_structure_string,
    BACA_CLINICAL_PHENO_FILE_PATH,
    ERG_PATTERN,
)
from baca_mmc5_primary_chain_cycle_structures import compute_chain_cycle_structure
from baca_aberration_clinical_analysis import _all_perfect_matchings, load_ets_status_from_mmc5

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
OUTPUT_CSV = os.path.join(RESULTS_DIR, "phase2_obligate_probable_completion.csv")

# Same numeric precedent as MAX_REJOINS=6 elsewhere in this codebase
# (n_DSBs <= 12 -> (12-1)!! = 11!! = 10,395 matchings).
MAX_DANGLING_ENDS = 12


# ── Real-graph primitives (built on baca_mmc5_chain_closure's validated pieces) ──

def _build_real_combined_graph(grp):
    """The real (rearrangement + chosen reference-matching) MultiGraph --
    identical construction to compute_chain_closure/compute_chain_cycle_structure,
    duplicated here (not imported) only because we need the graph object
    itself, not just its summary statistics."""
    nodes = set(grp['Breakpoint number'])
    rearr_adj = _build_rearr_matching(grp)
    ref_graph = _build_ref_graph(grp)
    ref_matching = nx.algorithms.matching.max_weight_matching(ref_graph, maxcardinality=True)

    combined = nx.MultiGraph()
    combined.add_nodes_from(nodes)
    seen_rearr_pairs = set()
    for a, b in rearr_adj.items():
        key = tuple(sorted((a, b)))
        if key not in seen_rearr_pairs:
            combined.add_edge(a, b)
            seen_rearr_pairs.add(key)
    for a, b in ref_matching:
        combined.add_edge(a, b)
    return combined


def _closed_and_open_components(combined):
    """
    Split the real combined graph's connected components into:
      closed_components : list[frozenset] -- already a real closed cycle
                           (every node degree == 2)
      open_components    : list[frozenset] -- Cornforth's "one-way"/open
                           pattern (some node degree == 1)
      open_ends          : sorted list of every degree-1 node across all
                           open components (the Phase 2 closure candidates)

    Asserts, per component, that an "open" component has EXACTLY 2 degree-1
    nodes -- guaranteed by the max-degree-2 property of this graph (see
    module docstring), checked here rather than just assumed.
    """
    closed, open_, open_ends = [], [], []
    for comp in nx.connected_components(combined):
        degrees = {n: combined.degree(n) for n in comp}
        if all(d == 2 for d in degrees.values()):
            closed.append(frozenset(comp))
        else:
            endpoint_nodes = sorted(n for n, d in degrees.items() if d == 1)
            assert len(endpoint_nodes) == 2, (
                f"Open component with {len(endpoint_nodes)} endpoint(s), expected "
                f"exactly 2 (max-degree-2 graph property violated): {comp}"
            )
            assert all(d <= 2 for d in degrees.values()), f"Node with degree > 2 found: {comp}"
            open_.append(frozenset(comp))
            open_ends.extend(endpoint_nodes)
    return closed, open_, sorted(open_ends)


def _structure_tuple(cycle_lengths):
    """
    Sorted-descending tuple of cycle lengths, e.g. (6, 2, 2) means C6+C2+C2.

    WHY PYTHON TUPLE COMPARISON == CORNFORTH'S OBLIGATE RULE: Cornforth
    defines the obligate (smallest) structure as whichever minimizes the
    largest cycle order, tie-breaking on the next-largest, and so on. That
    is *exactly* what Python's built-in tuple comparison does when
    comparing two sorted-descending tuples element by element. The one
    subtlety -- what happens if the tuples have different lengths -- never
    actually arises here: every candidate structure compared against
    another closes the exact same fixed set of breakpoints (the chain's
    total node count never changes), so any two candidate partitions
    always sum to the same total and can never be a strict prefix of one
    another. This means `min()` over a list of these tuples is always a
    safe, exact implementation of Cornforth's rule -- not an approximation.
    """
    return tuple(sorted(cycle_lengths, reverse=True))


def _enumerate_completions(closed_components, open_components, open_ends, chrom_lookup):
    """
    Enumerate every perfect matching of `open_ends` (Cornforth: every way to
    bring the pattern to closure), and for each, compute, IN ONE PASS (no
    re-enumeration needed downstream):
      - structure  : the resulting FULL cycle-structure tuple (pre-existing
        real closed cycles, unchanged, plus whatever new cycle(s) form from
        merging open-path components via this matching's invented edges)
      - rule2_metric : the MAXIMUM chromosome span among only the newly-
        closed cycles this matching creates (see module docstring point 4
        for why this is our adapted reading of Cornforth's Rule 2, not a
        literal one -- and NOT a sum: summing is biased towards merging,
        see the in-line note below).
      - max_span   : the largest chromosome span among ALL cycles in the
        resulting structure (pre-existing closed cycles + new ones) --
        needed for the chromoplexy flags, computed here rather than by
        re-deriving a candidate's merge after the fact.

    Implemented via union-find over `open_components` rather than rebuilding
    a full graph per matching, for speed -- mathematically equivalent to
    adding each matching's invented edges to the real graph and re-running
    connected-components, because every open component has EXACTLY 2
    endpoint nodes (asserted in _closed_and_open_components), so "pairing
    two open ends" is exactly "merging the (at most two) components those
    ends belong to."

    Returns: list of dicts with keys structure, rule2_metric, max_span --
    one per matching, in the exact deterministic order
    _all_perfect_matchings yields them for a pre-sorted, plain-integer input
    list (not subject to PYTHONHASHSEED randomization -- see module
    docstring).
    """
    base_lengths = [len(c) for c in closed_components]
    base_spans = [len({chrom_lookup[n] for n in c}) for c in closed_components]

    end_to_component_idx = {}
    for idx, comp in enumerate(open_components):
        for n in comp:
            if n in open_ends:
                end_to_component_idx[n] = idx

    n_components = len(open_components)
    results = []
    for matching in _all_perfect_matchings(open_ends):
        parent = list(range(n_components))

        def find(i):
            while parent[i] != i:
                i = parent[i]
            return i

        for a, b in matching:
            ia, ib = find(end_to_component_idx[a]), find(end_to_component_idx[b])
            if ia != ib:
                parent[ia] = ib

        groups = defaultdict(list)
        for i in range(n_components):
            groups[find(i)].append(i)

        new_lengths, new_chrom_spans = [], []
        for idxs in groups.values():
            merged_nodes = set()
            for i in idxs:
                merged_nodes |= open_components[i]
            new_lengths.append(len(merged_nodes))
            new_chrom_spans.append(len({chrom_lookup[n] for n in merged_nodes}))

        structure = _structure_tuple(base_lengths + new_lengths)
        # MUST be max(), not sum(): summing chromosome spans across the new
        # cycles is biased TOWARD merging whenever two open paths share a
        # chromosome, because span(A union B) <= span(A) + span(B) -- i.e. a
        # "sum" metric makes merging look cheaper purely from removing
        # double-counted shared chromosomes, which is backwards (it would
        # reward the MORE coordinated answer, not the conservative one
        # Cornforth's rule is supposed to produce). max() does not have this
        # flaw: span(A union B) >= max(span(A), span(B)) always, so merging
        # never looks artificially cheaper -- caught by hand-tracing a real
        # chain (P01-28 chain 1) where the sum-based version wrongly favored
        # merging two open paths before this fix.
        rule2_metric = max(new_chrom_spans) if new_chrom_spans else 0
        all_spans = base_spans + new_chrom_spans
        max_span = max(all_spans) if all_spans else 0
        results.append({
            'structure': structure, 'rule2_metric': rule2_metric, 'max_span': max_span,
            'matching': matching,
        })

    return results


def _select_obligate(enum_results):
    """Cornforth's two-step selection: (1) keep only matchings tied at the
    minimum Rule-2 metric, (2) among those, pick the lexicographically
    smallest structure (== min() over the tuples, see _structure_tuple)."""
    min_rule2 = min(r['rule2_metric'] for r in enum_results)
    rule2_survivors = [r for r in enum_results if r['rule2_metric'] == min_rule2]
    obligate_structure = min(r['structure'] for r in rule2_survivors)
    n_matchings_at_minimum = len(rule2_survivors)
    achieving = [r for r in rule2_survivors if r['structure'] == obligate_structure]
    return obligate_structure, min_rule2, n_matchings_at_minimum, achieving


def _select_probable(enum_results):
    """Mode of the structure-frequency distribution. Ties are broken
    deterministically (smallest tuple, same total order as obligate uses)
    but flagged explicitly via probable_is_tied / probable_tied_alternatives
    -- never silently resolved, consistent with how every other ambiguity
    in this project has been handled."""
    freq = Counter(r['structure'] for r in enum_results)
    max_count = max(freq.values())
    tied = sorted(s for s, c in freq.items() if c == max_count)
    probable_structure = tied[0]
    is_tied = len(tied) > 1
    tied_alternatives = tied[1:] if is_tied else []
    achieving = [r for r in enum_results if r['structure'] == probable_structure]
    return probable_structure, max_count, len(freq), is_tied, tied_alternatives, achieving


def _resolve_span(achieving_results, label):
    """
    A cycle-LENGTH structure (e.g. "two C4s") does not uniquely determine
    chromosome span -- two different matchings can produce the same length
    partition while spanning different chromosomes per cycle. Among all
    matchings that achieve the selected structure (already filtered by
    whatever selection rule -- obligate's Rule-2-then-lexicographic, or
    probable's frequency), report the MINIMUM span found (the most
    conservative realization of that exact structure, consistent with the
    conservative spirit used throughout this script) and flag explicitly
    whether that span was actually unique across all such matchings or
    whether other realizations of the identical structure have a larger
    span -- never silently averaged or hidden.
    """
    spans = [r['max_span'] for r in achieving_results]
    min_span = min(spans)
    is_unique = len(set(spans)) == 1
    if not is_unique:
        print(
            f"  [note] {label}: structure achieved by {len(spans)} matching(s) with "
            f"varying chromosome spans {sorted(set(spans))} -- reporting minimum ({min_span})."
        )
    return min_span, is_unique


def _chromoplexy_flags(max_span):
    return max_span >= 3, max_span >= 2


# ── Per-chain Phase 2 computation ───────────────────────────────────────────

def compute_phase2_completion(grp):
    chrom_lookup = dict(zip(grp['Breakpoint number'], grp['chromosome']))

    combined = _build_real_combined_graph(grp)
    closed_components, open_components, open_ends = _closed_and_open_components(combined)

    m = len(open_ends)
    assert m % 2 == 0, f"Odd number of open ends ({m}) -- should be impossible (handshake lemma)."
    assert m == 2 * len(open_components), "Open-ends count must be exactly 2x open-components count."

    base = {
        'n_dangling_ends': m,
        'n_open_path_components': len(open_components),
        'uses_invented_edges': m > 0,
        'n_invented_edges': m // 2,
    }

    if m > MAX_DANGLING_ENDS:
        base.update({
            'phase2_enumerable': False,
            'n_matchings_enumerated': None,
            'obligate_cycle_structure': None, 'obligate_n_cycles': None,
            'obligate_max_cycle_chrom_span': None, 'obligate_max_span_is_unique': None,
            'obligate_has_chromoplexy_strict': None, 'obligate_has_chromoplexy_loose': None,
            'obligate_rule2_metric': None, 'obligate_n_matchings_at_rule2_minimum': None,
            'obligate_n_matchings_achieving_structure': None,
            'probable_cycle_structure': None, 'probable_n_cycles': None,
            'probable_max_cycle_chrom_span': None, 'probable_max_span_is_unique': None,
            'probable_has_chromoplexy_strict': None, 'probable_has_chromoplexy_loose': None,
            'probable_pct': None, 'probable_n_distinct_structures': None,
            'probable_is_tied': None, 'probable_tied_alternatives': None,
        })
        return base

    enum_results = _enumerate_completions(closed_components, open_components, open_ends, chrom_lookup)
    total_matchings = len(enum_results)

    obligate_structure, obligate_rule2, n_at_min, obligate_achieving = _select_obligate(enum_results)
    (probable_structure, probable_count, n_distinct,
     probable_is_tied, tied_alts, probable_achieving) = _select_probable(enum_results)

    obligate_max_span, obligate_span_unique = _resolve_span(obligate_achieving, 'obligate')
    probable_max_span, probable_span_unique = _resolve_span(probable_achieving, 'probable')

    obl_strict, obl_loose = _chromoplexy_flags(obligate_max_span)
    prob_strict, prob_loose = _chromoplexy_flags(probable_max_span)

    base.update({
        'phase2_enumerable': True,
        'n_matchings_enumerated': total_matchings,

        'obligate_cycle_structure': _cycle_structure_string(list(obligate_structure)),
        'obligate_n_cycles': len(obligate_structure),
        'obligate_max_cycle_chrom_span': obligate_max_span,
        'obligate_max_span_is_unique': obligate_span_unique,
        'obligate_has_chromoplexy_strict': obl_strict,
        'obligate_has_chromoplexy_loose': obl_loose,
        'obligate_rule2_metric': obligate_rule2,
        'obligate_n_matchings_at_rule2_minimum': n_at_min,
        'obligate_n_matchings_achieving_structure': len(obligate_achieving),

        'probable_cycle_structure': _cycle_structure_string(list(probable_structure)),
        'probable_n_cycles': len(probable_structure),
        'probable_max_cycle_chrom_span': probable_max_span,
        'probable_max_span_is_unique': probable_span_unique,
        'probable_has_chromoplexy_strict': prob_strict,
        'probable_has_chromoplexy_loose': prob_loose,
        'probable_pct': round(100 * probable_count / total_matchings, 2),
        'probable_n_distinct_structures': n_distinct,
        'probable_is_tied': probable_is_tied,
        'probable_tied_alternatives': "; ".join(
            _cycle_structure_string(list(s)) for s in tied_alts
        ) if tied_alts else "",
    })
    return base


# ── Driver ───────────────────────────────────────────────────────────────────

def run_phase2_analysis(m5_df, clinical_df):
    rows = []
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    for (patient_id, chain_number), grp in m5_df.groupby(['Individual', 'Chain number']):
        # Real (pre-Phase-2) baseline, reused verbatim from the already-
        # published primary-chain script -- never recomputed independently,
        # to guarantee this table's "real_*" columns match
        # results/mmc5_primary_chain_cycle_structures.csv exactly.
        real = compute_chain_cycle_structure(grp)

        # Cross-check: this script's own from-scratch real-graph decomposition
        # must agree with `real` exactly, for every chain, before anything
        # else is trusted.
        combined = _build_real_combined_graph(grp)
        closed_components, _, _ = _closed_and_open_components(combined)
        chrom_lookup = dict(zip(grp['Breakpoint number'], grp['chromosome']))
        independent_structure = _cycle_structure_string(
            sorted((len(c) for c in closed_components), reverse=True)
        )
        assert independent_structure == real['cycle_structure'], (
            f"Real-structure mismatch for {patient_id} chain {chain_number}: "
            f"{independent_structure!r} != {real['cycle_structure']!r}"
        )

        phase2 = compute_phase2_completion(grp)

        b_i_whole_chain = grp.groupby('chromosome').size().to_dict()
        chroms_whole_chain = sorted(b_i_whole_chain.keys(), key=lambda x: (len(x), x))
        site_col = grp['Site annotation'] if 'Site annotation' in grp.columns else pd.Series([], dtype=str)
        contains_erg = bool(site_col.dropna().str.contains(ERG_PATTERN).any())
        genes = extract_genes(site_col.tolist())

        clin = clinical_df[clinical_df['Individual'] == patient_id]
        ets_status = ets_lookup.get(patient_id)
        gleason = clin['Gleason Score'].values[0] if len(clin) > 0 else None
        stage = clin['Pathological stage'].values[0] if len(clin) > 0 else None

        row = {
            'patient_id': patient_id,
            'baca_chain_number': chain_number,
            'n_breakpoints': len(grp),
            'n_rearrangements': grp['Rearrangement number'].nunique(),
            'k_chromosomes_whole_chain': len(chroms_whole_chain),
            'chromosomes_whole_chain': ",".join(chroms_whole_chain),
            'contains_ERG_or_TMPRSS2': contains_erg,
            'genes_in_chain': ",".join(genes),
            'ETS_status': ets_status,
            'Gleason_Score': gleason,
            'Pathological_stage': stage,

            'real_has_any_closed_cycle': real['has_any_closed_cycle'],
            'real_chain_fully_closed': real['chain_fully_closed'],
            'real_n_nodes_in_open_paths': real['n_open_breakpoints'],
            'real_pct_breakpoints_in_closed_cycles': real['pct_breakpoints_in_closed_cycles'],
            'real_n_cycles': real['n_cycles'],
            'real_cycle_structure': real['cycle_structure'],
            'real_theta_cycles': real['theta_cycles'],
            'real_chromosomes_in_cycles': real['chromosomes_in_cycles'],
            'real_max_cycle_chrom_span': real['max_cycle_chrom_span'],
            'real_has_chromoplexy_strict': real['has_chromoplexy_strict'],
            'real_has_chromoplexy_loose': real['has_chromoplexy_loose'],
        }
        row.update(phase2)
        rows.append(row)

    return pd.DataFrame(rows).sort_values(['patient_id', 'baca_chain_number']).reset_index(drop=True)


COLUMN_ORDER = [
    'patient_id', 'baca_chain_number',
    'n_breakpoints', 'n_rearrangements', 'k_chromosomes_whole_chain', 'chromosomes_whole_chain',
    'contains_ERG_or_TMPRSS2', 'genes_in_chain', 'ETS_status', 'Gleason_Score', 'Pathological_stage',

    'real_has_any_closed_cycle', 'real_chain_fully_closed', 'real_n_nodes_in_open_paths',
    'real_pct_breakpoints_in_closed_cycles', 'real_n_cycles', 'real_cycle_structure',
    'real_theta_cycles', 'real_chromosomes_in_cycles', 'real_max_cycle_chrom_span',
    'real_has_chromoplexy_strict', 'real_has_chromoplexy_loose',

    'n_dangling_ends', 'n_open_path_components', 'uses_invented_edges', 'n_invented_edges',
    'phase2_enumerable', 'n_matchings_enumerated',

    'obligate_cycle_structure', 'obligate_n_cycles', 'obligate_max_cycle_chrom_span',
    'obligate_max_span_is_unique',
    'obligate_has_chromoplexy_strict', 'obligate_has_chromoplexy_loose',
    'obligate_rule2_metric', 'obligate_n_matchings_at_rule2_minimum',
    'obligate_n_matchings_achieving_structure',

    'probable_cycle_structure', 'probable_n_cycles', 'probable_max_cycle_chrom_span',
    'probable_max_span_is_unique',
    'probable_has_chromoplexy_strict', 'probable_has_chromoplexy_loose',
    'probable_pct', 'probable_n_distinct_structures', 'probable_is_tied',
    'probable_tied_alternatives',
]


def main():
    clinical_df = pd.read_csv(BACA_CLINICAL_PHENO_FILE_PATH)
    m5_df = load_mmc5_table_s5a()

    df = run_phase2_analysis(m5_df, clinical_df)
    df = df[COLUMN_ORDER]
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {len(df)} chains ({df['patient_id'].nunique()} patients) -> {OUTPUT_CSV}")

    print()
    print("=== Scope / enumeration coverage ===")
    print(f"Chains with m=0 (already fully closed, no invention needed): {(df['n_dangling_ends'] == 0).sum()} / {len(df)}")
    print(f"Chains requiring invented closure (m>0): {df['uses_invented_edges'].sum()} / {len(df)}")
    print(f"Chains enumerable (m<=12): {df['phase2_enumerable'].sum()} / {len(df)}")
    not_enum = df[~df['phase2_enumerable']]
    print(f"Chains NOT enumerable (m>12): {len(not_enum)} -- {list(zip(not_enum['patient_id'], not_enum['baca_chain_number'], not_enum['n_dangling_ends']))}")

    print()
    print("=== Obligate vs probable: how often do they differ? (enumerable chains only) ===")
    enum_df = df[df['phase2_enumerable']].copy()
    invented_df = enum_df[enum_df['uses_invented_edges']]
    same = (invented_df['obligate_cycle_structure'] == invented_df['probable_cycle_structure']).sum()
    print(f"Of {len(invented_df)} chains requiring invention: obligate==probable for {same}, differ for {len(invented_df) - same}")
    print(f"Chains where probable_is_tied: {invented_df['probable_is_tied'].sum()} / {len(invented_df)}")

    print()
    print("=== Chromoplexy flags, before vs after completion (strict, >=3 chrom) ===")
    print(f"real_has_chromoplexy_strict   : {df['real_has_chromoplexy_strict'].sum()} / {len(df)}")
    print(f"obligate_has_chromoplexy_strict: {enum_df['obligate_has_chromoplexy_strict'].sum()} / {len(enum_df)} (enumerable only)")
    print(f"probable_has_chromoplexy_strict: {enum_df['probable_has_chromoplexy_strict'].sum()} / {len(enum_df)} (enumerable only)")


if __name__ == '__main__':
    main()
