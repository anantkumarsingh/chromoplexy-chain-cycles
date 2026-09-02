"""
baca_mmc5_chain_closure.py — reproduce Baca's own "closed chain" concept
using Baca's own chain boundaries and edges (mmc5.xlsx, Table S5A), and
compare the resulting closed-chain count against the paper's reported
121 closed chains (main text, Discussion section).

Why this is different from every other cycle-closure script in this repo:

  baca_full_genome_cycles.py / baca_erg_chain_completion.py both build the
  AMG from chrom_aberrations_baca.csv (the RAW rearrangement file) using
  only two edge types: rearrangement edges (real data) and reference edges
  WE infer from genomic position order. The chain boundary is also ours
  (full-genome connectivity, or one-hop ERG expansion).

  Here, every input is Baca's own: chain membership (Table S5A's
  'Chain number'), the rearrangement edge (two rows sharing the same
  'Rearrangement number'), AND a third edge type we were not previously
  using — the deletion bridge ('Deletion bridge partner breakpoint'),
  alongside the statistically-assigned adjacency edge
  ('Adjacent breakpoint(s)'). ChainFinder's methods section confirms these
  are exactly the three edge types it uses to build and statistically
  validate chains, and that it explicitly searches the resulting graph for
  closed cycles. We are not inferring a chain boundary or a reference edge
  here at all — we're using Baca's own.

Closure definition used here:
  IMPORTANT, verified empirically before settling on this design: not every
  breakpoint in a Baca chain has an adjacency or deletion-bridge partner —
  1,421 of 4,440 breakpoints cohort-wide (32%) have neither. So a chain
  cannot be required to close *entirely* (every single breakpoint forming
  one big loop) without throwing away real partial structure. Instead,
  for each chain we take the MAXIMUM (not necessarily perfect) matching of
  its ADJACENCY+DELETION-BRIDGE graph, combine it with the rearrangement
  matching (every breakpoint has exactly one rearrangement partner,
  verified globally), and inspect connected components of the combined
  graph:
    - a component where every node has degree exactly 2 is a genuine
      closed cycle (the two matchings' alternation returns to its start);
    - any other component (some node with degree 1) is an open path.
  Two columns capture both senses of "closed", since they answer different
  questions and Baca's paper doesn't specify which it means:
    - has_any_closed_cycle: chain contains >=1 closed cycle as a subset
      (matches the framing already used elsewhere in this repo, e.g.
      baca_erg_chain_completion.py's open/closed split within one chain)
    - chain_fully_closed: EVERY breakpoint in the chain is part of some
      closed cycle (the strict, whole-chain reading of "closed chain")

This is OUR reconstruction of Baca's unpublished closure criterion, built
only from columns Baca already computed — not the literal ChainFinder
source (not available to us). Flagged as such throughout.

Outputs (results/):
  mmc5_chain_closure_summary.csv      — all 366 chains, all 57 patients
  mmc5_chain_closure_erg_summary.csv  — subset: chains containing an
                                         ERG/TMPRSS2 site annotation
"""

import os
import re
import sys

import networkx as nx
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
from baca_aberration_clinical_analysis import BACA_CLINICAL_PHENO_FILE_PATH, load_ets_status_from_mmc5

BACA_DATASET_FOLDER = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/baca_dataset"
MMC5_FILE_PATH = os.path.join(BACA_DATASET_FOLDER, "mmc5.xlsx")

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
AGGREGATE_CSV = os.path.join(RESULTS_DIR, "mmc5_chain_closure_summary.csv")
ERG_CSV = os.path.join(RESULTS_DIR, "mmc5_chain_closure_erg_summary.csv")

GENE_PATTERN = re.compile(r'([A-Z][A-Z0-9\-]+)\([+-]\)')
ERG_PATTERN = re.compile(r'ERG|TMPRSS2')


# ── LOAD + CLEAN mmc5 Table S5A ─────────────────────────────────────────────

def _parse_adjacent(value):
    """'Adjacent breakpoint(s)' can hold 0, 1, or 2 breakpoint numbers,
    e.g. ' |', ' |37|', ' |40|174|'. Returns a list of ints."""
    if pd.isna(value):
        return []
    return [int(x) for x in re.findall(r'\d+', str(value))]


def _parse_deletion_bridge(value):
    if pd.isna(value):
        return None
    s = str(value).strip()
    if s == '':
        return None
    return int(float(s))


def load_mmc5_table_s5a():
    df = pd.ExcelFile(MMC5_FILE_PATH).parse('Table S5A')
    df.columns = [c.strip() for c in df.columns]
    df = df[df['Chromosome:position'].notna()].copy()
    df['chromosome'] = df['Chromosome:position'].str.split(':').str[0].str.replace('chr', '', regex=False)
    df['adj_list'] = df['Adjacent breakpoint(s)'].apply(_parse_adjacent)
    df['delb'] = df['Deletion bridge partner breakpoint'].apply(_parse_deletion_bridge)
    return df


def extract_genes(site_strings):
    genes = set()
    for s in site_strings:
        if pd.isna(s):
            continue
        genes.update(GENE_PATTERN.findall(str(s)))
    return sorted(genes)


# ── PER-CHAIN CLOSURE COMPUTATION ───────────────────────────────────────────

def _build_rearr_matching(grp):
    """Two rows sharing 'Rearrangement number' = one rearrangement edge.
    Verified globally: every group has exactly size 2, never crosses chains."""
    rearr_adj = {}
    for _, sub in grp.groupby('Rearrangement number'):
        bps = sub['Breakpoint number'].tolist()
        if len(bps) != 2:
            # Should never happen (verified globally) — skip defensively.
            continue
        a, b = bps
        rearr_adj[a] = b
        rearr_adj[b] = a
    return rearr_adj


def _build_ref_graph(grp):
    """Undirected graph of ADJACENCY + DELETION-BRIDGE edges, restricted to
    breakpoints within this chain (cross-chain references verified absent
    in mmc5, but defensively dropped here if encountered)."""
    chain_nodes = set(grp['Breakpoint number'])
    G = nx.Graph()
    G.add_nodes_from(chain_nodes)
    for _, row in grp.iterrows():
        bp = row['Breakpoint number']
        partners = set(row['adj_list'])
        if row['delb'] is not None:
            partners.add(row['delb'])
        for p in partners:
            if p in chain_nodes:
                G.add_edge(bp, p)
    return G


def _cycle_structure_string(cycle_lengths):
    if not cycle_lengths:
        return 'none'
    from collections import Counter
    counts = Counter(cycle_lengths)
    parts = [f"C{length}x{n}" for length, n in sorted(counts.items(), reverse=True)]
    return " + ".join(parts)


def compute_chain_closure(grp):
    """Returns a dict describing the closure result for one (Individual,
    Chain number) group. See module docstring for the has_any_closed_cycle
    vs chain_fully_closed distinction."""
    nodes = set(grp['Breakpoint number'])
    n_breakpoints = len(nodes)

    rearr_adj = _build_rearr_matching(grp)
    ref_graph = _build_ref_graph(grp)

    # Maximum (not necessarily perfect) matching of the reference layer.
    matching = nx.algorithms.matching.max_weight_matching(ref_graph, maxcardinality=True)

    # MUST be a MultiGraph: in a true C2 cycle the rearrangement edge and the
    # matched reference edge connect the EXACT same two nodes. A plain
    # nx.Graph collapses that into a single edge, leaving both nodes at
    # degree 1 instead of 2 and wrongly marking the simplest, most common
    # cycle as open — caught and fixed before trusting any result here.
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

    chrom_lookup = dict(zip(grp['Breakpoint number'], grp['chromosome']))
    cycle_lengths = []
    cycle_chrom_spans = []
    n_open_breakpoints = 0

    for component in nx.connected_components(combined):
        degrees = [combined.degree(n) for n in component]
        is_closed_cycle = len(component) >= 2 and all(d == 2 for d in degrees)
        if is_closed_cycle:
            cycle_lengths.append(len(component))
            cycle_chrom_spans.append(len({chrom_lookup[n] for n in component}))
        else:
            n_open_breakpoints += len(component)

    return {
        'has_any_closed_cycle': len(cycle_lengths) > 0,
        'chain_fully_closed': len(cycle_lengths) > 0 and n_open_breakpoints == 0,
        'n_open_breakpoints': n_open_breakpoints,
        'pct_breakpoints_in_closed_cycles': round(100 * (n_breakpoints - n_open_breakpoints) / n_breakpoints, 1),
        'cycle_structure': _cycle_structure_string(cycle_lengths),
        'n_cycles': len(cycle_lengths),
        'cycle_chrom_spans': cycle_chrom_spans,
    }


def run_closure_analysis(m5_df, clinical_df):
    rows = []
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    for (patient_id, chain_number), grp in m5_df.groupby(['Individual', 'Chain number']):
        result = compute_chain_closure(grp)

        b_i = grp.groupby('chromosome').size().to_dict()
        chroms = sorted(b_i.keys(), key=lambda x: (len(x), x))
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
            'k_chromosomes': len(chroms),
            'chromosomes': ",".join(chroms),
            'has_any_closed_cycle': result['has_any_closed_cycle'],
            'chain_fully_closed': result['chain_fully_closed'],
            'n_open_breakpoints': result['n_open_breakpoints'],
            'pct_breakpoints_in_closed_cycles': result['pct_breakpoints_in_closed_cycles'],
            'cycle_structure': result['cycle_structure'],
            'n_cycles': result['n_cycles'],
            'max_cycle_chrom_span': max(result['cycle_chrom_spans']) if result['cycle_chrom_spans'] else 0,
            'has_chromoplexy_strict': bool(result['cycle_chrom_spans']) and max(result['cycle_chrom_spans']) >= 3,
            'has_chromoplexy_loose': bool(result['cycle_chrom_spans']) and max(result['cycle_chrom_spans']) >= 2,
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

    summary_df = run_closure_analysis(m5_df, clinical_df)
    summary_df.to_csv(AGGREGATE_CSV, index=False)
    print(f"Wrote {len(summary_df)} chains ({summary_df['patient_id'].nunique()} patients) -> {AGGREGATE_CSV}")

    erg_df = summary_df[summary_df['contains_ERG_or_TMPRSS2']].copy()
    erg_df.to_csv(ERG_CSV, index=False)
    print(f"Wrote {len(erg_df)} ERG/TMPRSS2-containing chains ({erg_df['patient_id'].nunique()} patients) -> {ERG_CSV}")

    print()
    print("=== Aggregate closure results (all 366 chains, 57 patients) ===")
    print(f"Chains with >=1 closed cycle (subset OK): {summary_df['has_any_closed_cycle'].sum()} / {len(summary_df)}")
    print(f"Chains FULLY closed (every breakpoint in a cycle): {summary_df['chain_fully_closed'].sum()} / {len(summary_df)}")
    print("(Baca's paper reports 121 closed chains cohort-wide; exact definition not published)")
    print()
    print("By ETS status (has_any_closed_cycle):")
    print(summary_df.groupby('ETS_status')['has_any_closed_cycle'].agg(['sum', 'count']))
    print()
    print("Chromoplexy flags among chains with >=1 closed cycle:")
    any_closed = summary_df[summary_df['has_any_closed_cycle']]
    print(f"  strict (>=3 chrom in one cycle): {any_closed['has_chromoplexy_strict'].sum()} / {len(any_closed)}")
    print(f"  loose  (>=2 chrom in one cycle): {any_closed['has_chromoplexy_loose'].sum()} / {len(any_closed)}")

    print()
    print("=== ERG/TMPRSS2 chains only ===")
    print(f"Chains: {len(erg_df)}, patients: {erg_df['patient_id'].nunique()}")
    print(f"Chains with >=1 closed cycle: {erg_df['has_any_closed_cycle'].sum()} / {len(erg_df)}")
    print(f"Chains FULLY closed: {erg_df['chain_fully_closed'].sum()} / {len(erg_df)}")
    print("By ETS status (ERG chains, has_any_closed_cycle):")
    print(erg_df.groupby('ETS_status')['has_any_closed_cycle'].agg(['sum', 'count']))


if __name__ == '__main__':
    main()
