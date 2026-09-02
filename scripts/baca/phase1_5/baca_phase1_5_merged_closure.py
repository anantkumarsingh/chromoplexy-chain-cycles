"""
baca_phase1_5_merged_closure.py — Phase 1.5: merge Baca's own real reference
edges (mmc5.xlsx) into the full-genome AMG before any invented completion.

Why this phase exists (plain-language version, also in CLAUDE.md and memory):

  baca_full_genome_cycles.py builds the "normal connection" between two DNA
  breaks (the reference edge) with one simple rule: the next breakpoint
  after this one, in plain genomic position order, on the same chromosome.
  That rule is wrong whenever a deletion sits between two breakpoints — the
  truly adjacent DNA in the rearranged genome is no longer the positional
  neighbor. Baca's own ChainFinder already solved this and recorded the
  answer directly in mmc5.xlsx ('Adjacent breakpoint(s)' = statistically
  validated reference edge; 'Deletion bridge partner breakpoint' = a third
  edge type for breakpoints linked across two different fusions by a shared
  deleted region). We were not using either of mmc5's real edges in the
  full-genome closure — only our own position-order guess.

  This script does NOT invent anything. For every breakpoint that also
  appears in mmc5 (52/57 patients, and only the subset of each patient's
  rearrangements Baca's own algorithm clustered), it adds Baca's real
  adjacency/deletion-bridge connections as EXTRA candidate reference edges
  alongside our existing position-order guess — it never removes our
  fallback, it only gives the matching algorithm a better option to choose
  when one exists. Breakpoints absent from mmc5 keep using the old
  fallback rule unchanged. Any improvement in closure is therefore
  directly attributable to using more of Baca's own already-recorded data,
  not to a new modeling assumption.

  Obligate/minimal completion of whatever is STILL open after this (Phase
  2, inventing a closure) only makes sense once this real-data ceiling is
  known — otherwise we'd be inventing closures for gaps that were never
  actually real.

Method:
  For each patient, build the same node/edge primitives as
  baca_full_genome_cycles.py (_build_amg_nodes_and_edges), then build a
  reference-edge GRAPH (not a single-valued dict) per node:
    - always include the existing fallback partner (position order)
    - if this exact (chromosome, position, strand) breakpoint also exists
      in mmc5 for this patient, ALSO include mmc5's real adjacency and
      deletion-bridge partners as extra candidate edges
  Take the maximum matching of this graph, combine with the (always
  single-valued, real) rearrangement matching, and classify each connected
  component: every node at degree exactly 2 = closed cycle; anything else
  = open. Running this same algorithm with mmc5 edges OFF reproduces the
  existing baca_full_genome_cycles.py result exactly (verified against
  results/full_genome_cycle_compact_summary.csv) — so the only difference
  between the "old" and "new" columns in the output CSV is mmc5 edges
  being available or not.

Output:
  results/full_genome_cycle_phase1_5_comparison.csv — one row per patient,
  old (fallback-only) vs new (fallback+mmc5) cycle-closure stats side by
  side, plus delta/improvement columns.
"""

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

BACA_DATASET_FOLDER = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/baca_dataset"
MMC5_FILE_PATH = os.path.join(BACA_DATASET_FOLDER, "mmc5.xlsx")

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
COMPARISON_CSV = os.path.join(RESULTS_DIR, "full_genome_cycle_phase1_5_comparison.csv")

ERG_PATTERN = re.compile(r'ERG|TMPRSS2')
STRAND_MAP = {'Forward': '+', 'Reverse': '-'}


# ── LOAD mmc5 AND BUILD (chr,pos,strand)-KEYED REFERENCE EDGES ─────────────

def _parse_multi(value):
    """'Adjacent breakpoint(s)' can hold 0, 1, or 2 breakpoint numbers."""
    if pd.isna(value):
        return []
    return [int(x) for x in re.findall(r'\d+', str(value))]


def _parse_single(value):
    if pd.isna(value):
        return None
    s = str(value).strip()
    return int(float(s)) if s else None


def load_mmc5_reference_edges():
    """
    Returns:
        mmc5_node_partners : dict patient_id -> {node_tuple: set(partner_node_tuples)}
                              node_tuple = (chromosome_float, position_float, strand)
                              matching exactly the (chr,pos,strand) convention
                              used by _build_amg_nodes_and_edges on
                              chrom_aberrations_baca.csv.
    """
    df = pd.ExcelFile(MMC5_FILE_PATH).parse('Table S5A')
    df.columns = [c.strip() for c in df.columns]
    df = df[df['Chromosome:position'].notna()].copy()

    df['chromosome'] = (
        df['Chromosome:position'].str.split(':').str[0]
        .str.replace('chr', '', regex=False).astype(float)
    )
    df['position'] = df['Chromosome:position'].str.split(':').str[1].astype(float)
    df['strand_mapped'] = df['Strand'].map(STRAND_MAP)
    df['adj_list'] = df['Adjacent breakpoint(s)'].apply(_parse_multi)
    df['delb'] = df['Deletion bridge partner breakpoint'].apply(_parse_single)

    mmc5_node_partners = {}

    for patient_id, grp in df.groupby('Individual'):
        bp_to_node = {
            row['Breakpoint number']: (row['chromosome'], row['position'], row['strand_mapped'])
            for _, row in grp.iterrows()
        }

        node_partners = {}
        for _, row in grp.iterrows():
            node = bp_to_node[row['Breakpoint number']]
            partner_bps = set(row['adj_list'])
            if row['delb'] is not None:
                partner_bps.add(row['delb'])
            partner_nodes = {bp_to_node[p] for p in partner_bps if p in bp_to_node}
            node_partners[node] = partner_nodes

        mmc5_node_partners[patient_id] = node_partners

    return mmc5_node_partners


# ── PER-PATIENT CLOSURE (fallback-only OR fallback+mmc5) ───────────────────

def _chrom_set(chromosomes_str):
    return set(c for c in chromosomes_str.split(',') if c)


def _cycle_structure_string(cycle_lengths):
    if not cycle_lengths:
        return 'none'
    counts = Counter(cycle_lengths)
    return " + ".join(f"C{length}x{n}" for length, n in sorted(counts.items(), reverse=True))


def compute_patient_closure(patient_id, chrom_df, mmc5_node_partners, use_mmc5):
    """
    Computes cycle closure for one patient, genome-wide. If use_mmc5=False,
    this is methodologically identical to baca_full_genome_cycles.py's
    deterministic traversal (verified). If use_mmc5=True, mmc5's real
    adjacency/deletion-bridge edges are added as extra reference-edge
    candidates wherever available, on top of the same fallback.
    """
    patient_rows = _get_full_patient_rows(patient_id, chrom_df)
    if len(patient_rows) == 0:
        return None

    n_rejoins = len(patient_rows)
    n_dsbs = n_rejoins * 2

    all_nodes, rearr_adj, fallback_ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(patient_rows)
    b_i = {chrom: len(positions) for chrom, positions in chr_to_dsb_pos.items()}
    k = len(b_i)
    max_b_i = max(b_i.values()) if b_i else 0

    active_nodes = set(rearr_adj.keys())

    # Sort everything before building the graph. Python's hash randomization
    # (PYTHONHASHSEED, different every process launch) changes set/dict
    # iteration order for tuples containing strings — that changes edge
    # insertion order, which changes which maximum matching networkx picks
    # whenever multiple equally-maximal matchings exist (ties). Caught by
    # running this script 3x and getting 3 different "patients improved"
    # counts (5, 4, 5) before this fix — sorting makes every run identical.
    ref_graph = nx.Graph()
    ref_graph.add_nodes_from(sorted(all_nodes))
    for u, v in fallback_ref_adj.items():
        if v != 'TELOMERE':
            ref_graph.add_edge(u, v)

    n_mmc5_covered = 0
    if use_mmc5:
        node_partners = mmc5_node_partners.get(patient_id, {})
        for node in sorted(active_nodes):
            if node in node_partners:
                n_mmc5_covered += 1
                for partner in sorted(node_partners[node]):
                    if partner in all_nodes:
                        ref_graph.add_edge(node, partner)

    matching = nx.algorithms.matching.max_weight_matching(ref_graph, maxcardinality=True)

    # MUST be a MultiGraph: in a true C2 cycle the rearrangement edge and the
    # matched reference edge connect the EXACT same two nodes (e.g.
    # P05-1657's own ERG deletion). A plain nx.Graph collapses that into a
    # single edge, leaving both nodes at degree 1 instead of 2 and wrongly
    # marking the simplest, most common cycle as open — caught by comparing
    # against the existing full_genome_cycle_compact_summary.csv baseline
    # before trusting any result here.
    combined = nx.MultiGraph()
    combined.add_nodes_from(sorted(all_nodes))
    seen_pairs = set()
    for u, v in rearr_adj.items():
        key = tuple(sorted((u, v), key=str))
        if key not in seen_pairs:
            combined.add_edge(u, v)
            seen_pairs.add(key)
    for u, v in matching:
        combined.add_edge(u, v)

    cycle_lengths = []
    cycle_chrom_spans = []
    n_open_active_breakpoints = 0
    n_open_paths = 0
    closed_cycle_dsb_positions = set()  # (chrom, pos) pairs — DSBs actually inside a closed cycle

    for component in nx.connected_components(combined):
        active_in_component = component & active_nodes
        if not active_in_component:
            continue  # purely passive/free nodes never used by any active breakpoint
        degrees = [combined.degree(n) for n in component]
        is_closed = len(component) >= 2 and all(d == 2 for d in degrees)
        if is_closed:
            cycle_lengths.append(len(component))
            cycle_chrom_spans.append(len({n[0] for n in component}))
            closed_cycle_dsb_positions.update((n[0], n[1]) for n in component)
        else:
            n_open_active_breakpoints += len(active_in_component)
            n_open_paths += 1

    rearr_edges_in_cycles = sum(cycle_lengths) // 2
    frac_rejoins_in_cycles = rearr_edges_in_cycles / n_rejoins if n_rejoins else 0.0
    max_cycle_chrom_span = max(cycle_chrom_spans) if cycle_chrom_spans else 0

    # Theta notation and chromosome list SCOPED TO THE CLOSED-CYCLE SUBSET
    # ONLY — i.e. which chromosomes (and how many DSBs per chromosome) are
    # actually inside a real closed cycle, not the patient's whole genome
    # (which is identical between old/new since it's the same raw data —
    # only which DSBs CLOSE differs).
    b_i_cycles = Counter(chrom for chrom, _ in closed_cycle_dsb_positions)
    k_cycles = len(b_i_cycles)
    chromosomes_in_cycles = sorted(b_i_cycles.keys())
    theta_cycles = (
        f"Theta({k_cycles},({','.join(str(b_i_cycles[c]) for c in chromosomes_in_cycles)}))"
        if k_cycles else "Theta(0,())"
    )
    chromosomes_in_cycles_str = ",".join(str(int(c)) for c in chromosomes_in_cycles)

    return {
        'patient_id': patient_id,
        'k': k,
        'n_rejoins': n_rejoins,
        'n_dsbs': n_dsbs,
        'max_b_i': max_b_i,
        'n_mmc5_covered_breakpoints': n_mmc5_covered,
        'theta_cycles': theta_cycles,
        'chromosomes_in_cycles': chromosomes_in_cycles_str,
        'n_closed_cycles': len(cycle_lengths),
        'n_open_paths': n_open_paths,
        'n_open_breakpoints': n_open_active_breakpoints,
        'frac_rejoins_in_cycles': round(frac_rejoins_in_cycles, 3),
        'cycle_size_histogram': _cycle_structure_string(cycle_lengths),
        'max_cycle_chrom_span': max_cycle_chrom_span,
        'has_chromoplexy_strict': max_cycle_chrom_span >= 3,
        'has_chromoplexy_loose': max_cycle_chrom_span >= 2,
    }


# ── MAIN: run old vs new for all patients, write comparison CSV ────────────

def run_comparison(chrom_df, clinical_df, mmc5_node_partners):
    rows = []
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    for patient_id in sorted(chrom_df['Individual'].dropna().unique(), key=str):
        old = compute_patient_closure(patient_id, chrom_df, mmc5_node_partners, use_mmc5=False)
        new = compute_patient_closure(patient_id, chrom_df, mmc5_node_partners, use_mmc5=True)
        if old is None or new is None:
            continue

        patient_rows = _get_full_patient_rows(patient_id, chrom_df)
        contains_erg = bool(
            patient_rows['Breakpoint 1 site'].dropna().str.contains(ERG_PATTERN).any() or
            patient_rows['Breakpoint 2 site'].dropna().str.contains(ERG_PATTERN).any() or
            patient_rows['Fusion'].dropna().str.contains(ERG_PATTERN).any()
        )

        ets_status = ets_lookup.get(patient_id)

        rows.append({
            'patient_id': patient_id,
            'k': old['k'],
            'n_rejoins': old['n_rejoins'],
            'n_dsbs': old['n_dsbs'],
            'max_b_i': old['max_b_i'],
            'n_mmc5_covered_breakpoints': new['n_mmc5_covered_breakpoints'],
            'pct_dsbs_mmc5_covered': round(100 * new['n_mmc5_covered_breakpoints'] / old['n_dsbs'], 1) if old['n_dsbs'] else 0.0,

            'old_n_closed_cycles': old['n_closed_cycles'],
            'new_n_closed_cycles': new['n_closed_cycles'],
            'old_open_paths': old['n_open_paths'],
            'new_open_paths': new['n_open_paths'],
            'old_open_breakpoints': old['n_open_breakpoints'],
            'new_open_breakpoints': new['n_open_breakpoints'],
            'old_frac_rejoins_in_cycles': old['frac_rejoins_in_cycles'],
            'new_frac_rejoins_in_cycles': new['frac_rejoins_in_cycles'],
            'delta_frac_rejoins_in_cycles': round(new['frac_rejoins_in_cycles'] - old['frac_rejoins_in_cycles'], 3),
            'breakpoints_newly_closed': old['n_open_breakpoints'] - new['n_open_breakpoints'],

            'old_cycle_size_histogram': old['cycle_size_histogram'],
            'new_cycle_size_histogram': new['cycle_size_histogram'],

            'old_theta_cycles': old['theta_cycles'],
            'new_theta_cycles': new['theta_cycles'],
            'old_chromosomes_in_cycles': old['chromosomes_in_cycles'],
            'new_chromosomes_in_cycles': new['chromosomes_in_cycles'],
            'chromosomes_newly_in_cycles': ",".join(sorted(
                _chrom_set(new['chromosomes_in_cycles']) - _chrom_set(old['chromosomes_in_cycles']),
                key=int
            )),

            'old_max_cycle_chrom_span': old['max_cycle_chrom_span'],
            'new_max_cycle_chrom_span': new['max_cycle_chrom_span'],
            'old_has_chromoplexy_strict': old['has_chromoplexy_strict'],
            'new_has_chromoplexy_strict': new['has_chromoplexy_strict'],
            'old_has_chromoplexy_loose': old['has_chromoplexy_loose'],
            'new_has_chromoplexy_loose': new['has_chromoplexy_loose'],

            'improved_by_mmc5': new['frac_rejoins_in_cycles'] > old['frac_rejoins_in_cycles'],
            'newly_chromoplexy_strict': (not old['has_chromoplexy_strict']) and new['has_chromoplexy_strict'],
            'newly_chromoplexy_loose': (not old['has_chromoplexy_loose']) and new['has_chromoplexy_loose'],

            'contains_ERG_or_TMPRSS2': contains_erg,
            'ETS_status': ets_status,
        })

    return pd.DataFrame(rows)


def main():
    chrom_df, clinical_df = load_data(BACA_CHROM_ABER_FILE_PATH, BACA_CLINICAL_PHENO_FILE_PATH)
    mmc5_node_partners = load_mmc5_reference_edges()

    comparison_df = run_comparison(chrom_df, clinical_df, mmc5_node_partners)
    comparison_df.to_csv(COMPARISON_CSV, index=False)
    print(f"Wrote {len(comparison_df)} patients -> {COMPARISON_CSV}")

    print()
    print("=== Phase 1.5 summary: fallback-only vs fallback+mmc5 ===")
    print(f"Patients with ANY mmc5-covered breakpoints: {(comparison_df['n_mmc5_covered_breakpoints'] > 0).sum()} / {len(comparison_df)}")
    print(f"Patients improved by adding mmc5 edges: {comparison_df['improved_by_mmc5'].sum()} / {len(comparison_df)}")
    print(f"Total rearrangement-edges newly captured in closed cycles: {comparison_df['breakpoints_newly_closed'].sum() // 2}")
    print()
    print(f"Newly flagged chromoplexy (strict, >=3 chrom): {comparison_df['newly_chromoplexy_strict'].sum()}")
    print(f"Newly flagged chromoplexy (loose, >=2 chrom): {comparison_df['newly_chromoplexy_loose'].sum()}")
    print()
    print("Old vs new totals:")
    print(f"  Total closed cycles  -- old: {comparison_df['old_n_closed_cycles'].sum()}   new: {comparison_df['new_n_closed_cycles'].sum()}")
    print(f"  Total open paths     -- old: {comparison_df['old_open_paths'].sum()}   new: {comparison_df['new_open_paths'].sum()}")
    print(f"  Mean frac in cycles  -- old: {comparison_df['old_frac_rejoins_in_cycles'].mean():.3f}   new: {comparison_df['new_frac_rejoins_in_cycles'].mean():.3f}")
    print()
    print("By ETS status (improved_by_mmc5):")
    print(comparison_df.groupby('ETS_status')['improved_by_mmc5'].agg(['sum', 'count']))
    print()
    print("ERG/TMPRSS2-containing patients (improved_by_mmc5):")
    erg_only = comparison_df[comparison_df['contains_ERG_or_TMPRSS2']]
    print(f"  {erg_only['improved_by_mmc5'].sum()} / {len(erg_only)}")


if __name__ == '__main__':
    main()
