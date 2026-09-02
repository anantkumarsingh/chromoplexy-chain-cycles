"""
baca_mmc5_primary_chain_cycle_structures.py — the PRIMARY per-chain cycle
structure table for the publication.

Why this is now the primary substrate (full reasoning in CLAUDE.md
"FINAL DECISION: primary substrate is Baca's own mmc5 chains" and memory
cycle_completion_phase.md): the full-genome model (baca_full_genome_cycles.py,
baca_phase1_5_merged_closure.py) has two problems. (1) Its position-order
reference-edge rule is copy-number-blind — wherever a deletion bridge
exists between two breakpoints, that rule picks a biologically fictitious
"next breakpoint by position" connection, which can produce a FALSE
POSITIVE closed cycle or hide a real one (a FALSE NEGATIVE) — exactly the
professor's copy-number warning. (2) Merging in mmc5's real deletion-bridge
edges to fix that (baca_phase1_5_merged_closure.py) introduces a different
problem at full-genome scale: multiple structurally different maximum
matchings can exist for the same patient (proven on P04-1084 — a 138-pair
matching either way, but spanning 9 vs. 6 chromosomes depending on which
maximum matching is picked), and enumerating all of them is computationally
hard (#P-complete-adjacent) — not safely usable as-is.

Baca's own mmc5 chains (baca_mmc5_chain_closure.py) avoid both problems:
they already use all 3 real edge types (rearrangement, adjacency,
deletion bridge — no copy-number blindness), and they're scoped to much
smaller per-chain graphs (max 118 breakpoints vs. hundreds-to-thousands
genome-wide) where the same direct-matching-vs-force-propagation stability
check found only 1 disagreement across all 366 chains (P08-1042 chain 2,
size 18 vs 17 — a near-miss, not a structural problem). This is the
empirical basis for treating mmc5-chain closure as the primary,
trustworthy substrate, with the full-genome model kept only as a
secondary, broader-coverage (52/57 patients vs covers fewer rearrangements
per patient) cross-check.

This script extends baca_mmc5_chain_closure.py's per-chain closure result
with the same closed-cycle-SCOPED theta notation / chromosome list columns
already added to the full-genome Phase 1.5 comparison (theta_cycles,
chromosomes_in_cycles) — i.e. which chromosomes, and how many DSBs per
chromosome, are actually inside a real closed cycle for that chain, not
just the whole chain's chromosome membership.

Output:
  results/mmc5_primary_chain_cycle_structures.csv — one row per
  (patient, Baca chain number), all 366 chains, all 57 mmc5-covered
  patients.
"""

import os
import sys
from collections import Counter

import networkx as nx
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from baca_mmc5_chain_closure import (
    load_mmc5_table_s5a,
    extract_genes,
    _build_rearr_matching,
    _build_ref_graph,
    _cycle_structure_string,
    BACA_CLINICAL_PHENO_FILE_PATH,
    ERG_PATTERN,
)
from baca_aberration_clinical_analysis import load_ets_status_from_mmc5

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
OUTPUT_CSV = os.path.join(RESULTS_DIR, "mmc5_primary_chain_cycle_structures.csv")


def compute_chain_cycle_structure(grp):
    """
    Same closure computation as baca_mmc5_chain_closure.compute_chain_closure,
    extended to also return theta notation and chromosome list SCOPED TO
    THE CLOSED-CYCLE SUBSET ONLY (not the whole chain).
    """
    nodes = set(grp['Breakpoint number'])
    n_breakpoints = len(nodes)
    chrom_lookup = dict(zip(grp['Breakpoint number'], grp['chromosome']))

    rearr_adj = _build_rearr_matching(grp)
    ref_graph = _build_ref_graph(grp)
    matching = nx.algorithms.matching.max_weight_matching(ref_graph, maxcardinality=True)

    # MUST be a MultiGraph — see baca_mmc5_chain_closure.py docstring for why
    # (a C2 cycle's rearrangement edge and reference edge can connect the
    # exact same two nodes; a plain Graph would collapse that into degree 1).
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

    cycle_lengths = []
    cycle_chrom_spans = []
    n_open_breakpoints = 0
    closed_cycle_chrom_counts = Counter()  # chromosome (str) -> DSB count, closed-cycle nodes only

    for component in nx.connected_components(combined):
        degrees = [combined.degree(n) for n in component]
        is_closed_cycle = len(component) >= 2 and all(d == 2 for d in degrees)
        if is_closed_cycle:
            cycle_lengths.append(len(component))
            chroms_in_component = {chrom_lookup[n] for n in component}
            cycle_chrom_spans.append(len(chroms_in_component))
            for n in component:
                closed_cycle_chrom_counts[chrom_lookup[n]] += 1
        else:
            n_open_breakpoints += len(component)

    k_cycles = len(closed_cycle_chrom_counts)
    chromosomes_in_cycles = sorted(closed_cycle_chrom_counts.keys(), key=lambda x: (len(x), x))
    theta_cycles = (
        f"Theta({k_cycles},({','.join(str(closed_cycle_chrom_counts[c]) for c in chromosomes_in_cycles)}))"
        if k_cycles else "Theta(0,())"
    )

    return {
        'has_any_closed_cycle': len(cycle_lengths) > 0,
        'chain_fully_closed': len(cycle_lengths) > 0 and n_open_breakpoints == 0,
        'n_open_breakpoints': n_open_breakpoints,
        'pct_breakpoints_in_closed_cycles': round(100 * (n_breakpoints - n_open_breakpoints) / n_breakpoints, 1),
        'cycle_structure': _cycle_structure_string(cycle_lengths),
        'n_cycles': len(cycle_lengths),
        'theta_cycles': theta_cycles,
        'chromosomes_in_cycles': ",".join(chromosomes_in_cycles),
        'max_cycle_chrom_span': max(cycle_chrom_spans) if cycle_chrom_spans else 0,
        'has_chromoplexy_strict': bool(cycle_chrom_spans) and max(cycle_chrom_spans) >= 3,
        'has_chromoplexy_loose': bool(cycle_chrom_spans) and max(cycle_chrom_spans) >= 2,
    }


def run_primary_chain_analysis(m5_df, clinical_df):
    rows = []
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    for (patient_id, chain_number), grp in m5_df.groupby(['Individual', 'Chain number']):
        result = compute_chain_cycle_structure(grp)

        b_i_whole_chain = grp.groupby('chromosome').size().to_dict()
        chroms_whole_chain = sorted(b_i_whole_chain.keys(), key=lambda x: (len(x), x))
        site_col = grp['Site annotation'] if 'Site annotation' in grp.columns else pd.Series([], dtype=str)
        contains_erg = bool(site_col.dropna().str.contains(ERG_PATTERN).any())
        genes = extract_genes(site_col.tolist())

        clin = clinical_df[clinical_df['Individual'] == patient_id]
        ets_status = ets_lookup.get(patient_id)
        gleason = clin['Gleason Score'].values[0] if len(clin) > 0 else None
        stage = clin['Pathological stage'].values[0] if len(clin) > 0 else None

        rows.append({
            'patient_id': patient_id,
            'baca_chain_number': chain_number,
            'n_breakpoints': len(grp),
            'n_rearrangements': grp['Rearrangement number'].nunique(),
            'k_chromosomes_whole_chain': len(chroms_whole_chain),
            'chromosomes_whole_chain': ",".join(chroms_whole_chain),

            'has_any_closed_cycle': result['has_any_closed_cycle'],
            'chain_fully_closed': result['chain_fully_closed'],
            'n_open_breakpoints': result['n_open_breakpoints'],
            'pct_breakpoints_in_closed_cycles': result['pct_breakpoints_in_closed_cycles'],
            'n_cycles': result['n_cycles'],
            'cycle_structure': result['cycle_structure'],
            'theta_cycles': result['theta_cycles'],
            'chromosomes_in_cycles': result['chromosomes_in_cycles'],
            'max_cycle_chrom_span': result['max_cycle_chrom_span'],
            'has_chromoplexy_strict': result['has_chromoplexy_strict'],
            'has_chromoplexy_loose': result['has_chromoplexy_loose'],

            'contains_ERG_or_TMPRSS2': contains_erg,
            'genes_in_chain': ",".join(genes),
            'ETS_status': ets_status,
            'Gleason_Score': gleason,
            'Pathological_stage': stage,
        })

    return pd.DataFrame(rows).sort_values(['patient_id', 'baca_chain_number']).reset_index(drop=True)


def main():
    clinical_df = pd.read_csv(BACA_CLINICAL_PHENO_FILE_PATH)
    m5_df = load_mmc5_table_s5a()

    df = run_primary_chain_analysis(m5_df, clinical_df)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {len(df)} chains ({df['patient_id'].nunique()} patients) -> {OUTPUT_CSV}")

    print()
    print("=== Primary chain cycle-structure summary (Baca's own mmc5 chains) ===")
    print(f"Chains with >=1 closed cycle: {df['has_any_closed_cycle'].sum()} / {len(df)}")
    print(f"Chains FULLY closed: {df['chain_fully_closed'].sum()} / {len(df)}")
    print()
    print("By ETS status (has_any_closed_cycle):")
    print(df.groupby('ETS_status')['has_any_closed_cycle'].agg(['sum', 'count']))
    print()
    closed_only = df[df['has_any_closed_cycle']]
    print(f"Chromoplexy among closed-cycle chains -- strict (>=3 chrom): {closed_only['has_chromoplexy_strict'].sum()} / {len(closed_only)}")
    print(f"Chromoplexy among closed-cycle chains -- loose  (>=2 chrom): {closed_only['has_chromoplexy_loose'].sum()} / {len(closed_only)}")
    print()
    erg_df = df[df['contains_ERG_or_TMPRSS2']]
    print(f"ERG/TMPRSS2 chains: {len(erg_df)} ({erg_df['patient_id'].nunique()} patients)")
    print(f"  with >=1 closed cycle: {erg_df['has_any_closed_cycle'].sum()} / {len(erg_df)}")
    print(f"  fully closed: {erg_df['chain_fully_closed'].sum()} / {len(erg_df)}")


if __name__ == '__main__':
    main()
