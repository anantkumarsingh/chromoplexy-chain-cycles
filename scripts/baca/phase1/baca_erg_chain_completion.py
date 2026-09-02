"""
baca_erg_chain_completion.py  —  Phase 1b: complete the ERG chain open paths
using each patient's REAL full-genome rearrangement data (no hypothetical
matchings), and compare the resulting chromoplexy call against Baca's
15/26 heuristic.

Background: compute_erg_cycle_structures() in
baca_aberration_clinical_analysis.py builds its graph from
_get_erg_chain_rows — anchor rows (ERG/TMPRSS2) + a ONE-HOP inter_chr
expansion. Many of its "open paths" are an artifact of that restriction,
not a real biological gap: the dangling end's true rearrangement partner
is often recorded elsewhere in the same patient's full CSV, just outside
the one-hop chain.

Approach: don't restrict the EDGES, restrict the STARTING POINTS.
  1. Seed traversal from the original ERG chain's nodes (same anchor +
     one-hop set as before).
  2. But traverse using the FULL patient's real rearrangement/reference
     edges (built from every row for that patient, like in
     baca_full_genome_cycles.py) — not the chain-restricted edges.
  3. Because the traversal follows real edges wherever they actually lead,
     seeding only from ERG-chain nodes discovers exactly the real,
     data-grounded connected component reachable from the ERG/TMPRSS2
     locus — pulling in however many additional rows the real data
     requires, closing cycles only with edges that are actually in the
     data. This never wanders into totally unrelated clusters elsewhere
     in the genome, because nothing seeds a traversal there.
  4. Whatever is still open after exhausting all real data for that
     patient is tagged UNRESOLVED — its true partner genuinely doesn't
     appear anywhere in that patient's CSV. Sub-tagged by termination
     type (TELOMERE vs MISSING_REARR) since the later obligate /
     copy-number / proximity / centromere phase will likely treat those
     differently. UNRESOLVED is deferred to that phase, not guessed here.

This module reuses _get_erg_chain_rows / _build_amg_nodes_and_edges /
get_erg_chain_details_v2 from baca_aberration_clinical_analysis.py and
_get_full_patient_rows / classify_chromoplexy_from_cycles from
baca_full_genome_cycles.py — only the seeded-traversal logic and the
Baca-comparison reporting live here.

Output:
    results/erg_chain_completion_summary.csv — one row per ERG+ patient,
    chain-only ("before") vs real-data-completed ("after") stats, plus
    Baca's 15/26 heuristic call side by side with ours for direct comparison.
"""

import os
import sys
from collections import Counter

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
from baca_aberration_clinical_analysis import (
    load_data,
    _get_erg_chain_rows,
    _build_amg_nodes_and_edges,
    get_erg_chain_details_v2,
    BACA_CHROM_ABER_FILE_PATH,
    BACA_CLINICAL_PHENO_FILE_PATH,
)
from baca_full_genome_cycles import (
    _get_full_patient_rows,
    _find_amg_cycles_with_chroms,
    classify_chromoplexy_from_cycles,
)

RESULTS_DIR = (
    '/Users/anantkumarsingh/projects/prostate_cancer'
    '/nih-tcga-prad/results'
)
ERG_CHAIN_COMPLETION_CSV = os.path.join(RESULTS_DIR, 'erg_chain_completion_summary.csv')


# ── ERG+ patient list (same definition used in main(), reused directly) ────

def get_erg_positive_patients(clinical_df):
    """
    Same definition as main()'s erg_patients: clinically-confirmed ERG
    fusion by sequencing. 26 patients in the current dataset.
    """
    return clinical_df[
        clinical_df['ETS fusion detected by sequencing'].str.contains('ERG', na=False)
    ]['Individual'].tolist()


# ── Seeded traversal on the FULL real-data graph ────────────────────────────

def _find_amg_cycles_seeded(seed_nodes, full_all_nodes, rearr_adj_full, ref_adj_full):
    """
    Same alternating-edge traversal as _find_amg_cycles_with_chroms, but
    only STARTS new traversals from nodes in seed_nodes (the original ERG
    chain's nodes) — not every node in the patient's full genome graph.
    The edge dictionaries are still the FULL real-data graph, so a
    traversal can wander arbitrarily far from the seed along real edges;
    restricting the start set just keeps us from also discovering totally
    unrelated cycles elsewhere in the genome that have nothing to do with
    the ERG locus.

    Returns:
        cycles      : list of (cycle_length, set_of_chromosomes) — closed,
                      using only real data
        open_paths  : list of dicts {'path', 'termination', 'last_node'} —
                      still open even after using ALL real data for this
                      patient. termination is 'TELOMERE' or 'MISSING_REARR'.
        visited     : set of every node reached — the real-data-grounded
                      ERG-connected component (used to report how much
                      extra real data got pulled in beyond the one-hop chain)
    """
    visited = set()
    cycles = []
    open_paths = []

    candidates = sorted(n for n in seed_nodes if n in full_all_nodes)
    for start_node in candidates:
        if start_node in visited:
            continue
        if start_node not in rearr_adj_full:
            continue

        path = []
        current = start_node
        step = 'rearrangement'
        rearr_hops = 0
        is_open = False
        termination = None

        while current not in visited:
            visited.add(current)
            path.append(current)

            if step == 'rearrangement':
                nxt = rearr_adj_full.get(current)
                step = 'reference'
                if nxt is None:
                    is_open = True; termination = 'MISSING_REARR'; break
                rearr_hops += 1
            else:
                nxt = ref_adj_full.get(current)
                step = 'rearrangement'
                if nxt is None or nxt == 'TELOMERE':
                    is_open = True; termination = 'TELOMERE'; break

            current = nxt

        if (not is_open) and (current == start_node) and rearr_hops > 0:
            chroms = {n[0] for n in path}
            cycles.append((len(path), chroms))
        else:
            open_paths.append({'path': path, 'termination': termination, 'last_node': current})

    return cycles, open_paths, visited


# ── Per-patient before/after computation ────────────────────────────────────

def _cycle_histogram(cycle_lengths):
    counts = Counter(cycle_lengths)
    return (' + '.join(f'C{size}x{cnt}' for size, cnt in sorted(counts.items(), reverse=True))
            or 'none')


def complete_erg_chain_for_patient(patient_id, chrom_df, clinical_df, min_chromosomes=3):
    """
    Compare the original ERG chain (one-hop, possibly many open paths)
    against the real-data-completed version (seeded traversal on the full
    patient genome graph) for one patient. Also runs Baca's k_local /
    chain_size heuristic for a side-by-side comparison.
    """
    chain_rows = _get_erg_chain_rows(patient_id, chrom_df, clinical_df)
    if chain_rows.empty:
        return None

    # ---- BEFORE: original one-hop chain only ----
    all_nodes_chain, rearr_adj_chain, ref_adj_chain, chr_to_dsb_pos_chain = \
        _build_amg_nodes_and_edges(chain_rows)
    cycles_before, open_before = _find_amg_cycles_with_chroms(
        all_nodes_chain, rearr_adj_chain, ref_adj_chain)
    lengths_before = sorted((c[0] for c in cycles_before), reverse=True)

    n_rejoins_chain = len(chain_rows)
    n_dsbs_chain = n_rejoins_chain * 2
    k_local_chain = len(chr_to_dsb_pos_chain)

    # ---- AFTER: seeded traversal on the FULL real-data graph ----
    full_rows = _get_full_patient_rows(patient_id, chrom_df)
    all_nodes_full, rearr_adj_full, ref_adj_full, chr_to_dsb_pos_full = \
        _build_amg_nodes_and_edges(full_rows)

    cycles_after, open_after, visited = _find_amg_cycles_seeded(
        all_nodes_chain, all_nodes_full, rearr_adj_full, ref_adj_full)
    lengths_after = sorted((c[0] for c in cycles_after), reverse=True)

    chain_positions = {(n[0], n[1]) for n in all_nodes_chain}
    visited_positions = {(n[0], n[1]) for n in visited}
    n_extra_dsbs_pulled_in = len(visited_positions - chain_positions)

    n_unresolved_telomere = sum(1 for p in open_after if p['termination'] == 'TELOMERE')
    n_unresolved_missing  = sum(1 for p in open_after if p['termination'] == 'MISSING_REARR')

    chromoplexy_after = classify_chromoplexy_from_cycles(cycles_after, min_chromosomes)

    # ---- Baca's heuristic (k_local > 2 AND chain_size > 1) ----
    baca_details = get_erg_chain_details_v2(patient_id, chrom_df, clinical_df)
    baca_fusion_type = (
        'chromoplexy_embedded'
        if (baca_details['k_local'] or 0) > 2 and (baca_details['chain_size'] or 0) > 1
        else 'simple_fusion'
    )

    agrees_with_baca = (
        (baca_fusion_type == 'chromoplexy_embedded') == chromoplexy_after['has_chromoplexy']
    )

    return {
        'patient_id':                 patient_id,
        'k_local_chain':              k_local_chain,
        'n_rejoins_chain':            n_rejoins_chain,
        'n_dsbs_chain':               n_dsbs_chain,
        'open_paths_before':          len(open_before),
        'cycle_structure_before':     _cycle_histogram(lengths_before),
        'n_extra_dsbs_pulled_in':     n_extra_dsbs_pulled_in,
        'n_closed_cycles_after':      len(cycles_after),
        'cycle_structure_after':      _cycle_histogram(lengths_after),
        'max_cycle_chrom_span_after': chromoplexy_after['max_cycle_chrom_span'],
        'has_chromoplexy_after':      chromoplexy_after['has_chromoplexy'],
        'unresolved_total':           len(open_after),
        'unresolved_telomere':        n_unresolved_telomere,
        'unresolved_missing_rearr':   n_unresolved_missing,
        'baca_k_local':               baca_details['k_local'],
        'baca_chain_size':            baca_details['chain_size'],
        'baca_fusion_type':           baca_fusion_type,
        'agrees_with_baca':           agrees_with_baca,
    }


def run_erg_chain_completion(chrom_df, clinical_df, min_chromosomes=3, save_path=None):
    """
    Run complete_erg_chain_for_patient for every ERG+ patient. Saves and
    returns the summary DataFrame, and prints the Baca-vs-real-data
    comparison breakdown.
    """
    erg_patients = get_erg_positive_patients(clinical_df)
    results = []
    for pid in sorted(erg_patients, key=str):
        res = complete_erg_chain_for_patient(pid, chrom_df, clinical_df, min_chromosomes)
        if res is not None:
            results.append(res)

    summary_df = pd.DataFrame(results)
    if save_path:
        summary_df.to_csv(save_path, index=False)
        print(f"Saved ERG chain completion summary -> {save_path} ({len(summary_df)} patients)")

    print(f"\n{'='*70}")
    print(f"  ERG+ patients analyzed: {len(summary_df)} / {len(erg_patients)}")
    print(f"  Baca's heuristic        : "
          f"{(summary_df['baca_fusion_type'] == 'chromoplexy_embedded').sum()} chromoplexy_embedded, "
          f"{(summary_df['baca_fusion_type'] == 'simple_fusion').sum()} simple_fusion")
    print(f"  Real-data-completed cycles (min_chromosomes={min_chromosomes}): "
          f"{summary_df['has_chromoplexy_after'].sum()} chromoplexy, "
          f"{(~summary_df['has_chromoplexy_after']).sum()} no chromoplexy")
    print(f"  Agreement                : {summary_df['agrees_with_baca'].sum()} / {len(summary_df)}")
    print(f"{'='*70}")
    print(pd.crosstab(summary_df['baca_fusion_type'], summary_df['has_chromoplexy_after'],
                       rownames=['Baca heuristic'], colnames=['real-data cycle says chromoplexy']))
    print(f"{'='*70}\n")

    return summary_df


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    chrom_df, clinical_df = load_data(
        BACA_CHROM_ABER_FILE_PATH, BACA_CLINICAL_PHENO_FILE_PATH)
    run_erg_chain_completion(chrom_df, clinical_df, save_path=ERG_CHAIN_COMPLETION_CSV)
