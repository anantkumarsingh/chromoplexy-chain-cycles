import pandas as pd
import warnings
from scipy import stats


warnings.filterwarnings('ignore')

## CONFIG
BACA_DATASET_FOLDER = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/baca_dataset"
BACA_CHROM_ABER_FILE_PATH = BACA_DATASET_FOLDER + "/chrom_aberrations_baca.csv"
BACA_CLINICAL_PHENO_FILE_PATH = BACA_DATASET_FOLDER + "/clinical_phenotypes.csv"
BACA_MMC5_XLSX_PATH = BACA_DATASET_FOLDER + "/mmc5.xlsx"

def load_data(CHROM_PATH, CLIN_PATH):
    chrom = pd.read_csv(CHROM_PATH)
    clinical = pd.read_csv(CLIN_PATH)
    return chrom, clinical

def load_ets_status_from_mmc5(mmc5_path=BACA_MMC5_XLSX_PATH, known_patients=None):
    """
    Baca's own dedicated ETS fusion call, from mmc5.xlsx Table S5B's
    'ETS fusion' column (POSITIVE/NEGATIVE) — this is the correct source.

    Previously every script in this project (including this one, via the
    now-removed classify_ets()) inferred ETS status from the free-text
    'ETS fusion detected by sequencing' column in clinical_phenotypes.csv
    with the rule "any non-'---' value -> ETS+". That rule is wrong: it
    flags ANY detected fusion, not specifically an ETS-family fusion.
    Confirmed on patient PR-4240, whose fusion text is "Protein fusion:
    mid-exon (NRF1-BRAF)" — a real fusion, but NRF1-BRAF is not an
    ETS-family gene, and S5B correctly calls this patient NEGATIVE while
    the old free-text rule wrongly called it ETS+. Cross-checked against
    all 57 patients: this was the only disagreement.

    Table S5B also contains other-cohort (BR-/LUAD-/HN_ prefixed) and
    simulated-control rows sharing the same 'Individual' column; pass
    known_patients (e.g. clinical_df['Individual']) to restrict the
    lookup to our own 57-patient cohort.
    """
    s5b = pd.read_excel(mmc5_path, sheet_name='Table S5B')
    if known_patients is not None:
        s5b = s5b[s5b['Individual'].isin(set(known_patients))]
    label_map = {'POSITIVE': 'ETS+', 'NEGATIVE': 'ETS-'}
    s5b = s5b.dropna(subset=['ETS fusion'])
    return dict(zip(s5b['Individual'], s5b['ETS fusion'].map(label_map)))


def get_erg_chain_k(patient_id, chrom_df):
    patient = chrom_df[chrom_df['Individual'] == patient_id]

    # Step 1 — find rows directly involving ERG or TMPRSS2
    erg_direct = patient[
        patient['Breakpoint 1 site'].str.contains('ERG|TMPRSS2', na=False) |
        patient['Breakpoint 2 site'].str.contains('ERG|TMPRSS2', na=False)
        ]

    if len(erg_direct) == 0:
        return None, 0

    # Step 2 — get chromosomes involved in those direct rows
    chr_involved = set(erg_direct['Breakpoint 1 chromosome'].dropna()) | \
                   set(erg_direct['Breakpoint 2 chromosome'].dropna())

    # Step 3 — find ALL inter_chr rearrangements that share
    # at least one of those chromosomes
    # These are likely part of the same chain
    inter_rows = patient[patient['Class'] == 'inter_chr']
    chain_rows = inter_rows[
        inter_rows['Breakpoint 1 chromosome'].isin(chr_involved) |
        inter_rows['Breakpoint 2 chromosome'].isin(chr_involved)
        ]

    # Step 4 — update chromosomes from chain rows
    chain_chr = set(chain_rows['Breakpoint 1 chromosome'].dropna()) | \
                set(chain_rows['Breakpoint 2 chromosome'].dropna())

    # k = number of chromosomes in the ERG fusion chain
    k_local = len(chain_chr)
    n_chain_rearrangements = len(chain_rows)

    return k_local, n_chain_rearrangements


def get_erg_chain_details(patient_id, chrom_df):
    patient = chrom_df[chrom_df['Individual'] == patient_id]

    # Step 1 — find rows directly involving ERG or TMPRSS2
    erg_direct = patient[
        patient['Breakpoint 1 site'].str.contains('ERG|TMPRSS2', na=False) |
        patient['Breakpoint 2 site'].str.contains('ERG|TMPRSS2', na=False)
        ]

    if len(erg_direct) == 0:
        return None

    # Step 2 — get chromosomes from direct rows
    chr_involved = set(erg_direct['Breakpoint 1 chromosome'].dropna()) | \
                   set(erg_direct['Breakpoint 2 chromosome'].dropna())

    # Step 3 — find all inter_chr rows sharing those chromosomes
    inter_rows = patient[patient['Class'] == 'inter_chr']
    chain_rows = inter_rows[
        inter_rows['Breakpoint 1 chromosome'].isin(chr_involved) |
        inter_rows['Breakpoint 2 chromosome'].isin(chr_involved)
        ]

    # Step 4 — collect all genes from site annotations in the chain
    all_sites = list(chain_rows['Breakpoint 1 site'].dropna()) + \
                list(chain_rows['Breakpoint 2 site'].dropna())

    # Step 5 — collect chromosomes in chain
    chain_chr = set(chain_rows['Breakpoint 1 chromosome'].dropna()) | \
                set(chain_rows['Breakpoint 2 chromosome'].dropna())

    return {
        'patient_id': patient_id,
        'k_local': len(chain_chr),
        'chain_size': len(chain_rows),
        'chromosomes': sorted(chain_chr),
        'genes_in_chain': all_sites,
        'chain_rows': chain_rows[['Breakpoint 1 chromosome',
                                  'Breakpoint 1 site',
                                  'Breakpoint 2 chromosome',
                                  'Breakpoint 2 site',
                                  'Class']]
    }


# def get_erg_chain_details_v2(patient_id, chrom_df, clinical_df):
#     patient = chrom_df[chrom_df['Individual'] == patient_id]
#
#     # First try direct site annotation search
#     erg_direct = patient[
#         patient['Breakpoint 1 site'].str.contains('ERG|TMPRSS2', na=False) |
#         patient['Breakpoint 2 site'].str.contains('ERG|TMPRSS2', na=False) |
#         patient['Fusion'].str.contains('ERG|TMPRSS2', na=False)
#         ]
#
#     # If nothing found — check Fusion column more broadly
#     if len(erg_direct) == 0:
#         erg_direct = patient[
#             patient['Fusion'].str.contains('ERG|TMPRSS2', na=False)
#         ]
#
#     # If still nothing — patient is clinically ERG+
#     # so use chromosome 21 as anchor
#     # (TMPRSS2 and ERG are both on chr21)
#     if len(erg_direct) == 0:
#         ets_status = clinical_df[
#             clinical_df['Individual'] == patient_id
#             ]['ETS fusion detected by sequencing'].values[0]
#
#         if pd.notna(ets_status) and 'ERG' in str(ets_status):
#             # Use chr21 inter_chr rows as the chain
#             erg_direct = patient[
#                 (patient['Class'] == 'inter_chr') &
#                 (
#                         (patient['Breakpoint 1 chromosome'] == 21) |
#                         (patient['Breakpoint 2 chromosome'] == 21)
#                 )
#                 ]
#
#     if len(erg_direct) == 0:
#         return {
#             'patient_id': patient_id,
#             'k_local': 0,
#             'chain_size': 0,
#             'chromosomes': [],
#             'genes_in_chain': [],
#             'source': 'not_found'
#         }
#
#     # Get chromosomes from direct rows
#     chr_involved = set(erg_direct['Breakpoint 1 chromosome'].dropna()) | \
#                    set(erg_direct['Breakpoint 2 chromosome'].dropna())
#
#     # Expand to all inter_chr rows sharing those chromosomes
#     inter_rows = patient[patient['Class'] == 'inter_chr']
#     chain_rows = inter_rows[
#         inter_rows['Breakpoint 1 chromosome'].isin(chr_involved) |
#         inter_rows['Breakpoint 2 chromosome'].isin(chr_involved)
#         ]
#
#     chain_chr = set(chain_rows['Breakpoint 1 chromosome'].dropna()) | \
#                 set(chain_rows['Breakpoint 2 chromosome'].dropna())
#
#     all_sites = list(chain_rows['Breakpoint 1 site'].dropna()) + \
#                 list(chain_rows['Breakpoint 2 site'].dropna())
#
#     return {
#         'patient_id': patient_id,
#         'k_local': len(chain_chr),
#         'chain_size': len(chain_rows),
#         'chromosomes': sorted(chain_chr),
#         'genes_in_chain': all_sites,
#         'source': 'site_annotation' if len(erg_direct) > 0 else 'chr21_fallback'
#     }

def get_erg_chain_details_v2(patient_id, chrom_df, clinical_df):
    patient = chrom_df[chrom_df['Individual'] == patient_id]

    # Step 1 — direct site annotation search for ERG or TMPRSS2
    erg_direct = patient[
        patient['Breakpoint 1 site'].str.contains('ERG|TMPRSS2', na=False) |
        patient['Breakpoint 2 site'].str.contains('ERG|TMPRSS2', na=False)
        ]

    # Step 2 — check Fusion column if still empty
    if len(erg_direct) == 0:
        erg_direct = patient[
            patient['Fusion'].str.contains('ERG|TMPRSS2', na=False)
        ]

    # Step 3 — fallback: patient is clinically ERG+ so use chr21 inter_chr rows
    source = 'site_annotation'
    if len(erg_direct) == 0:
        ets_val = clinical_df[
            clinical_df['Individual'] == patient_id
            ]['ETS fusion detected by sequencing'].values

        if len(ets_val) > 0 and pd.notna(ets_val[0]) and 'ERG' in str(ets_val[0]):
            erg_direct = patient[
                (patient['Class'] == 'inter_chr') &
                (
                        (patient['Breakpoint 1 chromosome'] == 21) |
                        (patient['Breakpoint 2 chromosome'] == 21)
                )
                ]
            source = 'chr21_fallback'

    # Nothing found at all
    if len(erg_direct) == 0:
        return {
            'patient_id': patient_id,
            'k_local': 0,
            'chain_size': 0,
            'chromosomes': [],
            'genes_in_chain': [],
            'source': 'not_found'
        }

    # Get chromosomes from anchor rows
    chr_involved = set(erg_direct['Breakpoint 1 chromosome'].dropna()) | \
                   set(erg_direct['Breakpoint 2 chromosome'].dropna())

    # Expand to all inter_chr rows sharing those chromosomes
    inter_rows = patient[patient['Class'] == 'inter_chr']
    chain_rows = inter_rows[
        inter_rows['Breakpoint 1 chromosome'].isin(chr_involved) |
        inter_rows['Breakpoint 2 chromosome'].isin(chr_involved)
        ]

    # If no inter_chr chain found — fusion is intrachromosomal (simple)
    # k_local = number of chromosomes in the anchor rows only
    if len(chain_rows) == 0:
        chain_chr = chr_involved
        all_sites = list(erg_direct['Breakpoint 1 site'].dropna()) + \
                    list(erg_direct['Breakpoint 2 site'].dropna())
        return {
            'patient_id': patient_id,
            'k_local': len(chain_chr),
            'chain_size': 0,
            'chromosomes': sorted(chain_chr),
            'genes_in_chain': all_sites,
            'source': source + '_intra_only'
        }

    chain_chr = set(chain_rows['Breakpoint 1 chromosome'].dropna()) | \
                set(chain_rows['Breakpoint 2 chromosome'].dropna())

    all_sites = list(chain_rows['Breakpoint 1 site'].dropna()) + \
                list(chain_rows['Breakpoint 2 site'].dropna())

    return {
        'patient_id': patient_id,
        'k_local': len(chain_chr),
        'chain_size': len(chain_rows),
        'chromosomes': sorted(chain_chr),
        'genes_in_chain': all_sites,
        'source': source
    }


# =============================================================================
# AMG CYCLE STRUCTURE COMPUTATION — Sheth 2026 Framework
# =============================================================================
#
# GLOSSARY (Sheth notation):
#   n          = total DSBs in chain = chain_size × 2  (each row = 1 rejoin = 2 DSBs)
#   b_i        = DSBs on chromosome i within the chain
#   k_local    = number of chromosomes in the ERG chain
#   Θ(k,(b_i)) = configuration notation
#   C_m        = cycle of length m (m rearrangement edges + m reference edges)
#   rank(C)    = n - |C|  where |C| = number of cycles  [Theorem 25]
#
# DATA CONVENTION:
#   Each CSV row = one rejoin = one rearrangement edge
#   Strand '+' at a breakpoint = RIGHT side of the DSB (higher genomic coordinate)
#   Strand '-' at a breakpoint = LEFT  side of the DSB (lower  genomic coordinate)
#
# REFERENCE EDGE RULE (fixed by genome structure):
#   For consecutive DSBs at positions p_i < p_{i+1} on the same chromosome:
#       reference edge: (chr, p_i, '+') <--> (chr, p_{i+1}, '-')
#   This connects the right end of the left break to the left end of the right break.
#   Open ends (leftmost '-' and rightmost '+' on each chr) connect to telomeres.
#
# CYCLE TRAVERSAL:
#   Every AMG node has degree 2: one rearrangement edge + one reference edge.
#   So the graph decomposes uniquely into disjoint cycles (and open paths if
#   data is incomplete — path hits a telomere or a missing rearrangement).
#   Traversal: start at node → rearrangement → node → reference → node → ...
#   A closed loop = a cycle. Cycle length m = number of rearrangement hops.
#
# ENUMERATION (Table 2 style):
#   The reference matching is FIXED. We enumerate all (2n-1)!! possible
#   rearrangement matchings on the 2n endpoints, compute the cycle structure
#   for each, and count how many matchings produce each cycle structure.
#   Tractable only for small n: n=2 → 1, n=4 → 3, n=6 → 15, n=8 → 105,
#   n=10 → 945 (Sheth's P05-1657 result), n=12 → 10395.
#   For n > 12, use Sheth's theoretical bounds only.
# =============================================================================

from collections import defaultdict
import math
import itertools


def _get_erg_chain_rows(patient_id, chrom_df, clinical_df):
    """
    Return all rows belonging to the ERG fusion chain for this patient.

    Includes:
      - Anchor rows: rows directly annotating ERG or TMPRSS2 in site columns
      - Inter-chr chain rows: inter_chr rows sharing chromosomes with the anchor
      - For intra-only patients: all same-chromosome rows on chr21

    Returns a DataFrame subset of chrom_df.
    """
    patient = chrom_df[chrom_df['Individual'] == patient_id].copy()

    # Find anchor rows (direct ERG/TMPRSS2 annotation)
    anchor = patient[
        patient['Breakpoint 1 site'].str.contains('ERG|TMPRSS2', na=False) |
        patient['Breakpoint 2 site'].str.contains('ERG|TMPRSS2', na=False) |
        patient['Fusion'].str.contains('ERG|TMPRSS2', na=False)
    ]

    if len(anchor) == 0:
        return pd.DataFrame()

    anchor_chromosomes = (
        set(anchor['Breakpoint 1 chromosome'].dropna()) |
        set(anchor['Breakpoint 2 chromosome'].dropna())
    )

    # Expand to all inter_chr rows sharing anchor chromosomes (one hop)
    inter = patient[patient['Class'] == 'inter_chr']
    chain_inter = inter[
        inter['Breakpoint 1 chromosome'].isin(anchor_chromosomes) |
        inter['Breakpoint 2 chromosome'].isin(anchor_chromosomes)
    ]

    # Combine anchor + chain inter-chr rows (union, no duplicates)
    chain_rows = pd.concat([anchor, chain_inter]).drop_duplicates(subset='Number')
    return chain_rows


def _build_amg_nodes_and_edges(chain_rows):
    """
    Build the AMG node set and edge dictionaries from observed rejoin rows.

    Each DSB position creates TWO nodes: (chr, pos, '+') and (chr, pos, '-').
    Both nodes exist even if only one appears as a rearrangement endpoint in the data.

    Returns:
        all_nodes      : set of (chr, pos, strand) tuples
        rearr_adj      : dict (chr,pos,strand) -> (chr,pos,strand)  [rearrangement edges]
        ref_adj        : dict (chr,pos,strand) -> (chr,pos,strand) or 'TELOMERE'
        chr_to_dsb_pos : dict chr -> sorted list of positions
    """
    # ---- STEP 1: Collect all DSB positions from chain rows ----
    # Each row contributes 2 DSB positions (BP1 and BP2)
    dsb_positions = defaultdict(list)   # chr -> [pos, ...]

    for _, row in chain_rows.iterrows():
        c1 = row['Breakpoint 1 chromosome']
        p1 = row['Breakpoint 1 position']
        c2 = row['Breakpoint 2 chromosome']
        p2 = row['Breakpoint 2 position']
        dsb_positions[c1].append(p1)
        dsb_positions[c2].append(p2)

    # Deduplicate and sort positions per chromosome
    chr_to_dsb_pos = {c: sorted(set(ps)) for c, ps in dsb_positions.items()}

    # ---- STEP 2: Build full node set (both strands for every DSB position) ----
    all_nodes = set()
    for chrom, positions in chr_to_dsb_pos.items():
        for pos in positions:
            all_nodes.add((chrom, pos, '+'))
            all_nodes.add((chrom, pos, '-'))

    # ---- STEP 3: Rearrangement edges from observed CSV rows ----
    # Each row: (chr1, pos1, strand1) <--> (chr2, pos2, strand2)
    rearr_adj = {}  # node -> node  (symmetric)

    for _, row in chain_rows.iterrows():
        u = (row['Breakpoint 1 chromosome'],
             row['Breakpoint 1 position'],
             row['Breakpoint 1 strand'])
        v = (row['Breakpoint 2 chromosome'],
             row['Breakpoint 2 position'],
             row['Breakpoint 2 strand'])
        # PSEUDOCODE NOTE: In a proper AMG each endpoint appears in exactly one
        # rearrangement edge. If the same endpoint appears in two rows (data error),
        # the last assignment wins — flag this in production code.
        rearr_adj[u] = v
        rearr_adj[v] = u

    # ---- STEP 4: Reference edges from genomic ordering ----
    # Rule: for consecutive DSBs p_i < p_{i+1} on same chromosome:
    #   reference edge: (chr, p_i, '+') <--> (chr, p_{i+1}, '-')
    # Open ends (leftmost '-' and rightmost '+') -> 'TELOMERE'
    ref_adj = {}

    for chrom, sorted_positions in chr_to_dsb_pos.items():
        for rank, pos in enumerate(sorted_positions):
            node_plus  = (chrom, pos, '+')
            node_minus = (chrom, pos, '-')

            # '+' node: connect RIGHT to the next DSB's '-' node
            if rank + 1 < len(sorted_positions):
                next_pos = sorted_positions[rank + 1]
                ref_adj[node_plus] = (chrom, next_pos, '-')
            else:
                ref_adj[node_plus] = 'TELOMERE'   # rightmost DSB on chr

            # '-' node: connect LEFT to the previous DSB's '+' node
            if rank - 1 >= 0:
                prev_pos = sorted_positions[rank - 1]
                ref_adj[node_minus] = (chrom, prev_pos, '+')
            else:
                ref_adj[node_minus] = 'TELOMERE'  # leftmost DSB on chr

    return all_nodes, rearr_adj, ref_adj, chr_to_dsb_pos


def _find_amg_cycles(all_nodes, rearr_adj, ref_adj):
    """
    Find all cycles in the AMG by alternating-edge traversal.

    Traversal: start at a node -> follow rearrangement edge -> follow reference edge -> ...
    Closed loop = cycle. Open path = hits TELOMERE or missing rearrangement (incomplete data).

    Cycle length m = number of rearrangement hops (= Sheth's C_m notation).

    Returns:
        cycles     : list of int (cycle lengths, e.g. [3, 2, 1] = C3+C2+C1)
        open_paths : list of lists (paths that didn't close — incomplete data)
    """
    visited = set()
    cycles = []
    open_paths = []

    for start_node in sorted(all_nodes):  # sorted for deterministic order
        if start_node in visited:
            continue
        if start_node not in rearr_adj:
            # This endpoint has no observed rearrangement partner.
            # It will appear as part of an open path when reached via reference.
            continue

        path = []
        current = start_node
        step = 'rearrangement'  # always start with rearrangement step
        rearr_hops = 0
        is_open = False

        while current not in visited:
            visited.add(current)
            path.append(current)

            if step == 'rearrangement':
                nxt = rearr_adj.get(current)
                step = 'reference'
                if nxt is None:
                    is_open = True; break  # missing rejoin = open path
                rearr_hops += 1
            else:
                nxt = ref_adj.get(current)
                step = 'rearrangement'
                if nxt is None or nxt == 'TELOMERE':
                    is_open = True; break  # telomere connection = open path

            current = nxt

        if (not is_open) and (current == start_node) and rearr_hops > 0:
            # Sheth convention: C_n where n = number of DSBs in the cycle.
            # Each cycle visit alternates rearrangement/reference hops, visiting
            # 2*rearr_hops nodes total — each node is the rearrangement endpoint
            # of a distinct DSB. So n_DSBs_in_cycle = len(path) = 2 * rearr_hops.
            # e.g. Θ(1,(2)) simple deletion → 2 DSBs → C2, not C1.
            cycles.append(len(path))
        else:
            open_paths.append(path)

    return cycles, open_paths


def _all_perfect_matchings(elements):
    """
    Generator: yield every perfect matching of `elements` as a list of (a, b) pairs.
    Total count = (n-1)!! for n elements (n must be even).

    n=2  →       1  matching
    n=4  →       3
    n=6  →      15
    n=8  →     105
    n=10 →     945   (Sheth's P05-1657 result)
    n=12 →  10 395
    """
    if not elements:
        yield []
        return
    first = elements[0]
    rest  = elements[1:]
    for i, partner in enumerate(rest):
        remaining = rest[:i] + rest[i + 1:]
        for sub in _all_perfect_matchings(remaining):
            yield [(first, partner)] + sub


def _classify_matching(matching, endpoint_to_chr):
    """
    Return 'ALL_INTER', 'ALL_INTRA', or 'MIXED' for one matching.

    ALL_INTER: every pair connects endpoints on different chromosomes.
    ALL_INTRA: every pair connects endpoints on the same chromosome.
    MIXED    : some pairs are inter, some are intra.

    ALL_INTRA is impossible whenever any chromosome has an odd number of
    rearrangement endpoints — odd-count chromosomes cannot self-pair completely,
    so at least one endpoint must cross to another chromosome.
    """
    n_inter = sum(
        1 for (u, v) in matching
        if endpoint_to_chr[u] != endpoint_to_chr[v]
    )
    n_total = len(matching)
    if n_inter == n_total:  return 'ALL_INTER'
    if n_inter == 0:        return 'ALL_INTRA'
    return 'MIXED'


def _cycle_struct_from_matching(matching, ref_adj, all_nodes):
    """
    Build a temporary rearr_adj from this matching, run cycle-finding,
    return the cycle structure as a sorted-descending tuple of ints.
    e.g. (3, 2, 1) means C3 + C2 + C1.
    """
    temp_rearr = {}
    for (u, v) in matching:
        temp_rearr[u] = v
        temp_rearr[v] = u
    cycles, _ = _find_amg_cycles(all_nodes, temp_rearr, ref_adj)
    cycles.sort(reverse=True)
    return tuple(cycles)


def _enumerate_cycle_structures(rearr_endpoints, endpoint_to_chr, ref_adj, all_nodes):
    """
    Enumerate all (2n-1)!! perfect matchings of rearr_endpoints.
    For each matching: compute cycle structure and classify as
    ALL_INTER, ALL_INTRA, or MIXED.

    Returns:
        counts : dict with keys 'ALL_INTER', 'ALL_INTRA', 'MIXED'
                 each value is a dict: cycle_struct_tuple -> int count
        total  : total number of matchings enumerated (= (2n-1)!!)
    """
    counts = {
        'ALL_INTER': defaultdict(int),
        'ALL_INTRA': defaultdict(int),
        'MIXED':     defaultdict(int),
    }

    for matching in _all_perfect_matchings(list(rearr_endpoints)):
        category = _classify_matching(matching, endpoint_to_chr)
        struct   = _cycle_struct_from_matching(matching, ref_adj, all_nodes)
        counts[category][struct] += 1

    total = sum(sum(c.values()) for c in counts.values())
    return counts, total


def compute_erg_cycle_structures(patient_id, chrom_df, clinical_df):
    """
    Compute cycle structures for the ERG fusion chain — Sheth Table 2 style.

    Output includes:
      - Observed cycle structure from the actual rejoin data
      - Sheth theoretical bounds (Proposition 4, Theorem 6, Theorem 25)
      - Full enumeration of all possible cycle structures (for small chains)
      - Table 2 formatted output: cycle structure | count | % of pathways
    """

    # ---- STEP 1: Extract ERG chain rows ----
    chain_rows = _get_erg_chain_rows(patient_id, chrom_df, clinical_df)

    if len(chain_rows) == 0:
        print(f"[{patient_id}] No ERG chain rows found.")
        return None

    n_rejoins = len(chain_rows)
    n_dsbs = n_rejoins * 2   # each row = 2 DSB positions

    # ---- STEP 2: Build AMG ----
    all_nodes, rearr_adj, ref_adj, chr_to_dsb_pos = _build_amg_nodes_and_edges(chain_rows)

    # b_i = DSBs per chromosome IN THE CHAIN (not whole genome)
    b_i = {chrom: len(positions) for chrom, positions in chr_to_dsb_pos.items()}
    k_local = len(b_i)
    max_b_i = max(b_i.values()) if b_i else 0

    # ---- STEP 3: Find actual cycles from observed rejoins ----
    observed_cycles, open_paths = _find_amg_cycles(all_nodes, rearr_adj, ref_adj)

    observed_cycles.sort(reverse=True)
    observed_structure = tuple(observed_cycles)
    observed_str = ' + '.join(f'C{m}' for m in observed_cycles) or 'none (all open paths)'
    # Theorem 25: rank = n - |C|. Only valid when |C| >= 1 (proper AMG).
    # If |C| = 0, the chain has no closed cycles — data is incomplete subset of full AMG.
    observed_rank = (n_dsbs - len(observed_cycles)) if observed_cycles else None

    # ---- STEP 4: Sheth theoretical bounds ----
    # Proposition 4: cycle count in [1, min(max_b_i, n//2)]
    min_cycles_bound = 1
    max_cycles_bound = min(max_b_i, n_dsbs // 2)

    # Theorem 6: most probable cycle count ≈ log(n) for large n
    probable_cycles = max(1, round(math.log(n_dsbs))) if n_dsbs > 0 else 1

    # Lemma 14: if every chromosome in the chain has exactly 1 DSB (max_b_i == 1)
    # then the only possible cycle structure is C_n (rank = n-1, maximum complexity)
    lemma14_applies = (max_b_i == 1)

    # ---- STEP 5: Enumerate all possible cycle structures ----
    # We have 2*n_rejoins rearrangement endpoints (2 per observed row).
    # Total matchings = (2*n_rejoins - 1)!!
    #   n_rejoins=1 →       1,  n_rejoins=3 →      15
    #   n_rejoins=5 →     945,  n_rejoins=6 →   10395
    # Tractable threshold: n_rejoins <= 6
    TRACTABLE_THRESHOLD = 6

    # Map each rearrangement endpoint to its chromosome (for inter/intra classification)
    rearr_endpoints  = list(rearr_adj.keys())          # 2 * n_rejoins endpoints
    endpoint_to_chr  = {ep: ep[0] for ep in rearr_endpoints}

    # Which chromosomes have odd endpoint counts → ALL_INTRA structurally impossible
    chr_endpoint_counts = defaultdict(int)
    for ep in rearr_endpoints:
        chr_endpoint_counts[ep[0]] += 1
    odd_chr = [c for c, n in chr_endpoint_counts.items() if n % 2 != 0]
    all_intra_possible = (len(odd_chr) == 0)

    structure_counts = None
    total_pathways   = None

    if n_rejoins <= TRACTABLE_THRESHOLD:
        structure_counts, total_pathways = _enumerate_cycle_structures(
            rearr_endpoints, endpoint_to_chr, ref_adj, all_nodes
        )

    # ---- STEP 6: Print output ----
    W = 67  # line width

    def _print_table(category_label, cat_counts, total_cat, observed_structure):
        """Print one Table-2-style block for a single matching category."""
        if total_cat == 0:
            print(f"  (no matchings in this category)")
            return
        print(f"  {'Cycle structure':<28} {'Count':>7}  {'%':>7}")
        print(f"  {'-'*46}")
        for struct, cnt in sorted(cat_counts.items(), key=lambda x: -x[1]):
            pct  = 100.0 * cnt / total_cat
            sstr = ' + '.join(f'C{m}' for m in struct) if struct else '(none — open paths only)'
            mark = '  ← actual data' if struct == observed_structure else ''
            print(f"  {sstr:<28} {cnt:>7}  {pct:>6.1f}%{mark}")
        print(f"  {'─'*46}")
        print(f"  {'Total':<28} {total_cat:>7}  {'100.0%':>7}")

    print(f"\n{'='*W}")
    print(f"  Patient : {patient_id}")
    print(f"  Θ({k_local}, ({', '.join(str(b_i[c]) for c in sorted(b_i))}))")
    print(f"  n_rejoins={n_rejoins}  n_DSBs={n_dsbs}  k_local={k_local}")
    print(f"  NOTE: n_rejoins/b_i count only rows in the ERG chain")
    print(f"        (anchor rows + inter_chr rows sharing anchor chromosomes).")
    print(f"        Each chromosome likely has more DSBs outside this chain.")
    print(f"  b_i (DSBs per chr, chain only) : {dict(sorted(b_i.items()))}")
    print(f"  Open paths (DSBs whose rearr. partner is outside chain) : {len(open_paths)}")
    rank_str = str(observed_rank) if observed_rank is not None else 'N/A (no closed cycles in chain subset)'
    print(f"\n  CYCLE STRUCTURE FROM ACTUAL DATA")
    print(f"  (computed by AMG traversal on the real sequencing-detected rejoins)")
    print(f"  Result : {observed_str}")
    print(f"  n_cycles={len(observed_cycles)}  rank={rank_str}")
    if observed_cycles:
        print(f"  (max rank in proper AMG = n-1 = {n_dsbs - 1})")
    print(f"\n  SHETH BOUNDS:")
    print(f"    Prop 4  — cycle count ∈ [{min_cycles_bound}, {max_cycles_bound}]")
    print(f"    Thm 6   — probable cycles ≈ log(n) = log({n_dsbs}) ≈ {math.log(n_dsbs):.2f} → {probable_cycles}")
    print(f"    Thm 25  — rank(actual data) = {rank_str}")
    print(f"    Lem 14  — Lemma 14 regime (all b_i=1) : {lemma14_applies}")
    print(f"    ALL_INTRA possible : {all_intra_possible}", end='')
    if not all_intra_possible:
        odd_strs = [f"chr{int(c)} (n={chr_endpoint_counts[c]}, odd)" for c in sorted(odd_chr)]
        print(f"  — impossible because {', '.join(odd_strs)}")
    else:
        print()

    if structure_counts is not None:
        total_inter = sum(structure_counts['ALL_INTER'].values())
        total_intra = sum(structure_counts['ALL_INTRA'].values())
        total_mixed = sum(structure_counts['MIXED'].values())

        print(f"\n  COUNTERFACTUAL ENUMERATION")
        print(f"  ─────────────────────────────────────────────────────────────")
        print(f"  The tables below ask: given the SAME {n_dsbs} DSB positions and")
        print(f"  strands, what cycle structure would result if the repair")
        print(f"  machinery had rejoined the endpoints differently?")
        print(f"  Each of the {total_pathways} rows is one possible rejoining scenario.")
        print(f"  Only ONE scenario matches the actual sequencing data — marked '← actual data'.")
        print(f"  All others are counterfactual (did not happen in this tumor).")
        print(f"  ─────────────────────────────────────────────────────────────")
        print(f"  TOTAL COUNTERFACTUAL PATHWAYS : {total_pathways}")
        print(f"  (ALL_INTER={total_inter}, ALL_INTRA={total_intra}, MIXED={total_mixed})")

        print(f"\n  ── TABLE A : ALL_INTER matchings ({total_inter} / {total_pathways})")
        print(f"     [every hypothetical rejoin crosses chromosome boundaries]")
        _print_table('ALL_INTER', structure_counts['ALL_INTER'], total_inter, observed_structure)

        print(f"\n  ── TABLE B : ALL_INTRA matchings ({total_intra} / {total_pathways})")
        print(f"     [every hypothetical rejoin stays within one chromosome]")
        if not all_intra_possible:
            print(f"  (structurally impossible — see odd-count chromosomes above)")
        else:
            _print_table('ALL_INTRA', structure_counts['ALL_INTRA'], total_intra, observed_structure)

        print(f"\n  ── TABLE C : MIXED matchings ({total_mixed} / {total_pathways})")
        print(f"     [hypothetical rejoins include both inter- and intra-chromosomal pairs]")
        _print_table('MIXED', structure_counts['MIXED'], total_mixed, observed_structure)

    else:
        print(f"\n  COUNTERFACTUAL ENUMERATION : not computed")
        print(f"  Chain has {n_rejoins} rejoins → (2×{n_rejoins}-1)!! matchings exceeds")
        print(f"  tractable threshold ({TRACTABLE_THRESHOLD} rejoins / 10,395 matchings).")
        print(f"  Use Sheth theoretical bounds above for analysis.")

    print(f"{'='*W}")

    return {
        'patient_id':        patient_id,
        'k_local':           k_local,
        'n_rejoins':         n_rejoins,
        'n_dsbs':            n_dsbs,
        'b_i':               b_i,
        'observed_structure':observed_str,
        'observed_n_cycles': len(observed_cycles),
        'observed_rank':     observed_rank,
        'min_cycles_bound':  min_cycles_bound,
        'max_cycles_bound':  max_cycles_bound,
        'probable_cycles':   probable_cycles,
        'lemma14_applies':   lemma14_applies,
        'all_intra_possible':all_intra_possible,
        'open_paths':        len(open_paths),
        'structure_counts':  structure_counts,
        'total_pathways':    total_pathways,
    }


def main():
    baca_chrom, baca_clinical = load_data(BACA_CHROM_ABER_FILE_PATH, BACA_CLINICAL_PHENO_FILE_PATH)


    # ---- CHROM ANALYSIS ----
    patients = baca_chrom.groupby('Individual')
    summary_results = []
    chromosome_results = []
    for patient_id, df in patients:
        df = df.sort_values(by='Breakpoint 1 chromosome') # Sorts by Ascending order of Breakpoint 1 Chromosome Number
        total_rearrangements = len(df) # Number of rows for each patients equals num of rearrangements
        inter_chrom = len(df[df['Class'] == 'inter_chr'])
        intra_classes = ['long_range', 'inversion',
                         'potential deletion', 'tandem_dup']
        intra_chrom = len(df[df['Class'].isin(intra_classes)])
        pct_inter = inter_chrom / total_rearrangements * 100
        bp_1_chromosomes = df['Breakpoint 1 chromosome']
        bp_2_chromosomes = df['Breakpoint 2 chromosome']
        total_chrom = pd.concat([bp_1_chromosomes, bp_2_chromosomes])
        unique_chrom = total_chrom.nunique() # this is our k from the mathematical notation
        tot_dsb = total_rearrangements * 2
        dbs_per_chromosome = total_chrom.value_counts().sort_index()

        summary_results.append({
            'Individual': patient_id,
            'total_rearrangements': total_rearrangements,
            'inter_chr_count': inter_chrom,
            'intra_chr_count': intra_chrom,
            'pct_inter': pct_inter,
            'k': unique_chrom,
            'total_DSBs':tot_dsb,
        })

        for chrom, count in dbs_per_chromosome.items():
            chromosome_results.append({
                'Individual': patient_id,
                'chromosome': chrom,
                'dsb_count': count
            })

    summary_df = pd.DataFrame(summary_results)
    chromosome_df = pd.DataFrame(chromosome_results)

    #print(summary_df.shape)
    #print(chromosome_df.shape)

    # ---- CLINICAL ANALYSIS ----
    clinical_df = pd.DataFrame(baca_clinical)
    clinical_df = clinical_df.drop(columns = ['Disease', 'ERG FISH result', 'Cohort',
                     'ABSOLUTE purity estimation', 'WGS-based purity estimation',
                      'Exome-sequenced in Barbieri et al. 2012'])

    #print(clinical_df.head())

    # Joining the 2 Dataframes together
    baca_df = pd.merge(summary_df, clinical_df, on='Individual')
    #print(baca_df.head())

    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    baca_df['ETS_status'] = baca_df['Individual'].map(ets_lookup)

    baca_df['has_tertiary_pattern'] = baca_df['Gleason Score'].str.contains(';').fillna(False)

    baca_df['gleason_primary'] = baca_df['Gleason Score'].str.split('+').str[0]
    baca_df['gleason_secondary'] = (baca_df['Gleason Score']
    .str.split('+').str[1]
    .str.split(';').str[0])

    baca_df['gleason_primary'] = pd.to_numeric(baca_df['gleason_primary'], errors='coerce')
    baca_df['gleason_secondary'] = pd.to_numeric(baca_df['gleason_secondary'], errors='coerce')
    baca_df['gleason_total'] = baca_df['gleason_primary'] + baca_df['gleason_secondary']

    stage_table = baca_df.groupby('Pathological stage').agg(
        n_patients=('Individual', 'count'),
        median_total_R=('total_rearrangements', 'median'),
        median_inter=('inter_chr_count', 'median'),
        median_intra=('intra_chr_count', 'median'),
        median_pct_inter=('pct_inter', 'median'),
        median_k=('k', 'median'),
        median_total_DSBs=('total_DSBs', 'median'),
        ets_positive=('ETS_status', lambda x: (x == 'ETS+').sum())
    ).reset_index()

    #print(stage_table.to_string())
    gleason_table = baca_df.groupby('Gleason Score').agg(
        n_patients=('Individual', 'count'),
        median_total_R=('total_rearrangements', 'median'),
        median_inter=('inter_chr_count', 'median'),
        median_pct_inter=('pct_inter', 'median'),
        median_k=('k', 'median'),
        median_total_DSBs=('total_DSBs', 'median'),
        ets_positive=('ETS_status', lambda x: (x == 'ETS+').sum())
    ).reset_index()

    #print(gleason_table.to_string())
    erg_patients = clinical_df[
        clinical_df['ETS fusion detected by sequencing']
        .str.contains('ERG', na=False)
    ]['Individual'].tolist()

    print(f"ERG fusion patients: {len(erg_patients)}")

    erg_rows = baca_chrom[
        baca_chrom['Breakpoint 1 site'].str.contains('ERG|TMPRSS2', na=False) |
        baca_chrom['Breakpoint 2 site'].str.contains('ERG|TMPRSS2', na=False)
        ]

    # Patients with ERG/TMPRSS2 in breakpoints
    erg_breakpoint_patients = set(erg_rows['Individual'].unique())

    # Get AMG metrics for all ERG fusion patients (from clinical ETS column)
    erg_clinical_patients = set(erg_patients)

    # Combine both — union of patients identified by either method
    # all_erg_patients = erg_clinical_patients.union(erg_breakpoint_patients)

    # Get their AMG summary
    erg_amg = baca_df[baca_df['Individual'].isin(erg_clinical_patients)][
        ['Individual', 'k', 'total_DSBs', 'pct_inter',
         'inter_chr_count', 'intra_chr_count',
         'total_rearrangements', 'ETS_status',
         'Gleason Score', 'Pathological stage']
    ].copy()
    #
    # # Classify as chromoplexy-embedded vs simple fusion
    # # Chromoplexy = multiple chromosomes involved
    # # A simple isolated TMPRSS2-ERG fusion involves:
    # #   chromosome 21 (ERG) + chromosome 21 (TMPRSS2)
    # #   both genes are on chromosome 21 — so k could be as low as 1
    # #   OR TMPRSS2 deletion bridge involves nearby chromosomes
    #
    # # Use k > 2 as the threshold
    # # Simple fusion: k <= 2 (just chr21 rearranging with itself or one partner)
    # # Chromoplexy embedded: k > 2 (multiple chromosomes involved)
    # erg_amg['fusion_type'] = erg_amg['k'].apply(
    #     lambda x: 'chromoplexy_embedded' if x > 2 else 'simple_fusion'
    # )
    #
    # print("=== ERG FUSION PATIENTS AMG COMPLEXITY ===")
    # print(erg_amg.sort_values('k', ascending=False).to_string())
    #
    # print(f"\n=== FUSION TYPE COUNTS ===")
    # print(erg_amg['fusion_type'].value_counts())
    #
    # print(f"\n=== CHROMOPLEXY EMBEDDED — AMG SUMMARY ===")
    # chrom_embedded = erg_amg[erg_amg['fusion_type'] == 'chromoplexy_embedded']
    # print(f"n = {len(chrom_embedded)}")
    # print(f"median k           : {chrom_embedded['k'].median()}")
    # print(f"median total_DSBs  : {chrom_embedded['total_DSBs'].median()}")
    # print(f"median pct_inter   : {chrom_embedded['pct_inter'].median():.1f}%")
    # print(f"k > 2 count        : {(chrom_embedded['k'] > 2).sum()}")
    # print(f"k > 3 count        : {(chrom_embedded['k'] > 3).sum()}")
    #
    # print(f"\n=== SIMPLE FUSION — AMG SUMMARY ===")
    # simple = erg_amg[erg_amg['fusion_type'] == 'simple_fusion']
    # print(f"n = {len(simple)}")
    # print(f"median k           : {simple['k'].median()}")
    # print(f"median total_DSBs  : {simple['total_DSBs'].median()}")
    # print(f"median pct_inter   : {simple['pct_inter'].median():.1f}%")
    #
    # # Statistical test
    # from scipy import stats
    # stat, p = stats.mannwhitneyu(
    #     chrom_embedded['k'],
    #     simple['k'],
    #     alternative='greater'
    # )
    # print(f"\n=== MANN-WHITNEY: k chromoplexy_embedded vs simple ===")
    # print(f"U statistic : {stat:.3f}")
    # print(f"P-value     : {p:.4f}")

    # Apply to all 26 ERG patients


    erg_amg['k_local'] = None
    erg_amg['chain_size'] = None
    erg_amg['source'] = ''

    for idx, row in erg_amg.iterrows():
        result = get_erg_chain_details_v2(row['Individual'], baca_chrom, clinical_df)
        erg_amg.at[idx, 'k_local'] = result['k_local']
        erg_amg.at[idx, 'chain_size'] = result['chain_size']
        erg_amg.at[idx, 'source'] = result['source']

    erg_amg['k_local'] = pd.to_numeric(erg_amg['k_local'], errors='coerce')
    erg_amg['chain_size'] = pd.to_numeric(erg_amg['chain_size'], errors='coerce')

    # Classify using k_local threshold
    erg_amg['fusion_type'] = erg_amg.apply(
        lambda row: 'chromoplexy_embedded'
        if row['k_local'] > 2 and row['chain_size'] > 1
        else 'simple_fusion',
        axis=1
    )

    # ---- PRINT RESULTS ----
    print("\n=== ERG FUSION PATIENTS — AMG COMPLEXITY ===")
    print(erg_amg[['Individual', 'k', 'k_local', 'chain_size',
                   'fusion_type', 'source', 'Gleason Score',
                   'Pathological stage']].sort_values('k_local', ascending=False).to_string())

    print(f"\n=== FUSION TYPE COUNTS ===")
    print(erg_amg['fusion_type'].value_counts())

    print(f"\n=== SOURCE BREAKDOWN ===")
    print(erg_amg['source'].value_counts())

    print(f"\n=== CHROMOPLEXY EMBEDDED — AMG SUMMARY ===")
    chrom_embedded = erg_amg[erg_amg['fusion_type'] == 'chromoplexy_embedded']
    print(f"n                  : {len(chrom_embedded)}")
    print(f"median k_local     : {chrom_embedded['k_local'].median()}")
    print(f"median chain_size  : {chrom_embedded['chain_size'].median()}")
    print(f"median total_DSBs  : {chrom_embedded['total_DSBs'].median()}")
    print(f"median pct_inter   : {chrom_embedded['pct_inter'].median():.1f}%")
    print(f"k_local > 2 count  : {(chrom_embedded['k_local'] > 2).sum()}")
    print(f"k_local > 3 count  : {(chrom_embedded['k_local'] > 3).sum()}")

    print(f"\n=== SIMPLE FUSION — AMG SUMMARY ===")
    simple = erg_amg[erg_amg['fusion_type'] == 'simple_fusion']
    print(f"n                  : {len(simple)}")
    print(f"median k_local     : {simple['k_local'].median()}")
    print(f"median chain_size  : {simple['chain_size'].median()}")
    print(f"median total_DSBs  : {simple['total_DSBs'].median()}")
    print(f"median pct_inter   : {simple['pct_inter'].median():.1f}%")

    # Statistical test — only if both groups have values
    # if len(chrom_embedded) > 0 and len(simple) > 0:
    #     stat, p = stats.mannwhitneyu(
    #         chrom_embedded['k_local'],
    #         simple['k_local'],
    #         alternative='greater'
    #     )
        # print(f"\n=== MANN-WHITNEY: k_local chromoplexy_embedded vs simple ===")
        # print(f"U statistic : {stat:.3f}")
        # print(f"P-value     : {p:.4f}")

    # ---- CHAIN DETAILS PER PATIENT ----
    print("\n=== ERG FUSION CHAIN DETAILS PER PATIENT ===\n")
    for patient_id in erg_amg['Individual'].tolist():
        result = get_erg_chain_details_v2(patient_id, baca_chrom, clinical_df)
        fusion_type = erg_amg[erg_amg['Individual'] == patient_id]['fusion_type'].values[0]
        print(f"Patient: {result['patient_id']}  [{fusion_type}]  source={result['source']}")
        print(f"  k_local    : {result['k_local']}")
        print(f"  chain_size : {result['chain_size']} rearrangements")
        print(f"  chromosomes: {result['chromosomes']}")
        if result['genes_in_chain']:
            print(f"  genes at breakpoints:")
            for site in result['genes_in_chain']:
                print(f"      {site}")
        print()

    # ---- AMG CYCLE STRUCTURES (Sheth Table 2 style) ----
    # Runs for all 26 ERG fusion patients.
    # Patients with n_rejoins <= 6 get full enumeration (Tables A/B/C).
    # Patients with n_rejoins > 6 get Sheth theoretical bounds only.
    print("\n" + "#" * 67)
    print("# AMG CYCLE STRUCTURES — Sheth 2026 Framework")
    print("# Enumeration threshold: n_rejoins <= 6  ((2n-1)!! <= 10395)")
    print("#" * 67)

    all_cycle_results = []
    erg_patient_list = erg_amg.sort_values('k_local', ascending=False)['Individual'].tolist()

    for patient_id in erg_patient_list:
        result = compute_erg_cycle_structures(patient_id, baca_chrom, clinical_df)
        if result is not None:
            fusion_type = erg_amg[erg_amg['Individual'] == patient_id]['fusion_type'].values[0]
            result['fusion_type'] = fusion_type
            all_cycle_results.append(result)

    # ---- SUMMARY TABLE across all patients ----
    print("\n\n" + "=" * 67)
    print("  CYCLE STRUCTURE SUMMARY — all 26 ERG patients")
    print("=" * 67)
    print(f"  {'Patient':<22} {'Type':<22} {'Θ':<18} {'Observed':<22} {'Paths':<8} {'Rank'}")
    print(f"  {'-'*100}")
    for r in all_cycle_results:
        theta = f"Θ({r['k_local']},({','.join(str(r['b_i'][c]) for c in sorted(r['b_i']))}))"
        rank  = str(r['observed_rank']) if r['observed_rank'] is not None else 'N/A'
        paths = str(r['total_pathways']) if r['total_pathways'] is not None else '>10395'
        print(f"  {r['patient_id']:<22} {r['fusion_type']:<22} {theta:<18} "
              f"{r['observed_structure']:<22} {paths:<8} {rank}")

if __name__ == "__main__":
    main()






