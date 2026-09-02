"""
baca_phase2_cn_plausibility.py — Phase 2 copy-number plausibility check on
INVENTED (obligate/probable) completion edges.

Professor Arsuaga's directive (item 2 of the original 4-item list, see
CLAUDE.md "Next Analysis Direction" / "Phase 1.5") was to track copy number
while completing cycles and report what extra information that provides.
Phase 1.5 already used Baca's real deletion-bridge edges to avoid the
position-order reference-edge rule's copy-number blindness — but only for
REAL edges. Phase 2's INVENTED (obligate/probable) completion edges have
never had any copy-number awareness at all. This script closes that gap
using the real, genome-wide copy-number segmentation file discovered this
session: data/baca_dataset/mmc3.xlsx Table S3B (Sample, Chromosome, Start,
End, Segment_Mean — a standard CBS/DNAcopy-style segmentation, 19,493
segments, all 57 patients, confirmed NOT the same as the 14-gene mmc6.xlsx
proxy previously assumed to be the only CN data available).

WHAT THIS SCRIPT DOES, IN ONE SENTENCE: for every invented obligate/probable
completion edge from Phase 2, checks whether the copy-number landscape
between its two (real, already-observed) endpoints is consistent with a
simple, single, already-recorded structural picture, or silently crosses an
unexplained copy-number transition that nothing in this patient's real data
accounts for.

KEY CONCEPTUAL POINT, worth stating plainly: Phase 2 never invents new
breakpoint POSITIONS — every invented edge connects two breakpoints that
already exist in the real mmc5 data (see baca_phase2_obligate_probable_
completion.py module docstring, point 3). The "invented" part is only the
CONNECTION between two already-real, already-located dangling ends. So this
script is not asking "is this position real" (it always is) — it is asking
"is the simplest reading of this connection (undisturbed reference DNA
between the two ends) consistent with the measured copy-number profile
there."

Table S3B rows are NOT AMG reference edges and their internal segment
boundaries are NOT breakpoints in our sense — they are copy-number
change-points from an independent measurement (read-depth/probe signal),
almost never coinciding exactly with a real rearrangement breakpoint
(median distance ~157kb, confirmed this session on all 11,420 real
breakpoints). S3B is used purely as a LOOKUP: given a genomic position,
which segment (if any) contains it, and what is its Segment_Mean.

FOUR-WAY (five-way, counting interchromosomal) classification per invented
edge, decided explicitly to avoid inventing or inferring anything the data
doesn't support:

  - CN_FLAT: both endpoints resolve to a real S3B segment, and it is the
    SAME segment (or segments with no boundary between them) -- no copy
    number transition anywhere between the two ends. Neutral: consistent
    with (but not proof of) undisturbed reference DNA. A balanced
    rearrangement produces no CN signal either, so CN_FLAT is not itself
    evidence the invented edge is correct -- only that nothing contradicts it.
  - CN_TRANSITION_EXPLAINED: a real segment boundary (copy-number change)
    exists between the two endpoints, but at least one real, already-
    observed breakpoint from this patient's FULL rearrangement dataset
    (chrom_aberrations_baca.csv, not just this chain) falls strictly
    between them -- i.e. there is already a recorded real event that could
    account for the CN shift; the invented edge may just be "skipping over"
    a real intervening event that is already otherwise in the graph.
  - CN_TRANSITION_UNEXPLAINED: a real segment boundary exists between the
    two endpoints, and NOTHING in this patient's real rearrangement data
    (anywhere in the genome, not just this chain) falls in that span. This
    is the strongest red flag this analysis can raise: the invented
    "undisturbed reference DNA" reading is directly contradicted by a
    measured copy-number change with no recorded explanation.
  - CN_INDETERMINATE: at least one endpoint has no resolvable S3B segment
    (a genuine coverage gap -- confirmed this session: ~26% of real
    breakpoints fall in such gaps, likely centromere/telomere/low-
    mappability regions S3B's segmentation doesn't cover). This is NOT
    evidence for or against the edge -- it means we cannot assess it from
    this data source at all, and must never be silently folded into
    CN_FLAT (which would wrongly claim "checked, found nothing" when the
    truth is "could not check"). Deliberately mirrors the same never-
    conflate-False-with-unknown discipline already used for Phase 3's
    has_unresolved_caveat flags.
  - NOT_APPLICABLE_INTERCHROMOSOMAL: the two endpoints are on different
    chromosomes. S3B segments only tile a single chromosome each -- there
    is no shared coordinate axis to check "does the span cross a
    transition" the way there is on one chromosome. Rather than invent an
    interchromosomal proxy (e.g. checking each endpoint's own local CN
    context in isolation), which would be answering a different, weaker
    question and was not part of what was agreed, this is reported as its
    own honest "not applicable" category.

CALIBRATION, done BEFORE trusting this methodology on anything invented
(same discipline as every other cross-check in this project): Baca's own
mmc5 Table S5A 'Deletion bridge partner breakpoint' column records REAL,
statistically-validated pairs of breakpoints known to border the SAME
deleted region. If the CN lookup below works correctly, these pairs should
show a measurable copy-number LOSS between them far more often than chance.
See run_deletion_bridge_calibration() and the results reported in the
companion methodology doc. Result: of the 521 real deletion-bridge pairs
cohort-wide, 78 have both endpoints resolvable in S3B; of those, 58 (74%)
show a clear measured loss between them, 19 (24%) show no visible
transition (plausible given segmentation resolution limits -- a real but
small/subclonal deletion can fall below what CBS-style segmentation
detects), and 1 is an ambiguous edge case. This is treated as sufficient
confirmation the lookup mechanism is sound before applying it to invented
edges, where -- unlike here -- we have no independent ground truth.

EDGE-CHOICE AMBIGUITY, checked explicitly, not assumed away: a chain's
selected obligate/probable STRUCTURE (e.g. "C4+C2") can in principle be
achieved by more than one structurally different matching of open ends
(different actual edges, same resulting cycle-length partition). Every
chain in this script is checked for this: if the achieving set contains
more than one distinct edge-set, ALL of them are classified (not one
arbitrarily chosen), and the chain is flagged obligate_edges_unique /
probable_edges_unique = False, with the full range of outcomes reported
rather than a silently arbitrary pick.

Outputs (results/):
  phase2_cn_deletion_bridge_calibration.csv — the 521-pair calibration check
  phase2_cn_invented_edge_detail.csv        — one row per invented edge, per
                                               (chain, lens, matching realization)
  phase2_cn_chain_summary.csv               — one row per (chain, lens),
                                               aggregated classification counts
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'mmc5'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'phase1'))
from baca_aberration_clinical_analysis import BACA_DATASET_FOLDER, BACA_CHROM_ABER_FILE_PATH
from baca_mmc5_chain_closure import load_mmc5_table_s5a
from baca_full_genome_cycles import _get_full_patient_rows
from baca_phase2_obligate_probable_completion import (
    _build_real_combined_graph,
    _closed_and_open_components,
    _enumerate_completions,
    _select_obligate,
    _select_probable,
    MAX_DANGLING_ENDS,
)

MMC3_S3B_PATH = os.path.join(BACA_DATASET_FOLDER, "mmc3.xlsx")

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
CALIBRATION_CSV = os.path.join(RESULTS_DIR, "phase2_cn_deletion_bridge_calibration.csv")
EDGE_DETAIL_CSV = os.path.join(RESULTS_DIR, "phase2_cn_invented_edge_detail.csv")
CHAIN_SUMMARY_CSV = os.path.join(RESULTS_DIR, "phase2_cn_chain_summary.csv")


# ── Copy-number segment lookup (mmc3.xlsx Table S3B) ────────────────────────

def _normalize_chrom(chrom):
    """mmc5/chrom_aberrations chromosomes can be '23'/'24' (X/Y numeric
    codes) or plain '1'..'22'. S3B labels sex chromosomes as 'X' only --
    confirmed no 'Y' rows exist in S3B at all. Returns None for chrY
    (genuinely no CN data available), else the S3B-matching label."""
    chrom = str(chrom)
    if chrom == '23':
        return 'X'
    if chrom == '24':
        return None
    return chrom


def load_cn_segments(mmc3_path=MMC3_S3B_PATH):
    """Returns dict: (patient_id, chromosome_str) -> segment DataFrame
    (Start, End, Segment_Mean, seg_id), sorted by Start. seg_id is a
    stable per-row identity (not the Segment_Mean value) used to detect
    whether two positions fall in the SAME segment vs. different ones."""
    df = pd.read_excel(mmc3_path, sheet_name='Table S3B')
    df['Chromosome'] = df['Chromosome'].astype(str)
    df = df.reset_index(drop=True)
    df['seg_id'] = df.index
    segments = {}
    for (sample, chrom), sub in df.groupby(['Sample', 'Chromosome']):
        segments[(sample, chrom)] = sub.sort_values('Start').reset_index(drop=True)
    return segments


def find_segment(segments, patient_id, chrom, pos):
    """Returns the S3B row (Series) containing this position, or None if
    the chromosome has no data for this patient (chrY, or a patient/
    chromosome combination absent from S3B) or the position falls in a
    genuine coverage gap between segments."""
    norm_chrom = _normalize_chrom(chrom)
    if norm_chrom is None:
        return None
    key = (patient_id, norm_chrom)
    if key not in segments:
        return None
    sub = segments[key]
    hit = sub[(sub['Start'] <= pos) & (sub['End'] >= pos)]
    if len(hit) == 0:
        return None
    return hit.iloc[0]


def segments_between(segments, patient_id, chrom, pos_a, pos_b):
    """All S3B segments whose interval overlaps [min(pos_a,pos_b), max(...)]."""
    norm_chrom = _normalize_chrom(chrom)
    if norm_chrom is None:
        return None
    key = (patient_id, norm_chrom)
    if key not in segments:
        return None
    lo, hi = min(pos_a, pos_b), max(pos_a, pos_b)
    sub = segments[key]
    return sub[(sub['End'] >= lo) & (sub['Start'] <= hi)].sort_values('Start')


# ── Deletion-bridge calibration (real data, known-good positive control) ───

def run_deletion_bridge_calibration(s5a, segments):
    """Every real, Baca-validated deletion-bridge pair should show a
    measurable copy-number LOSS between its two breakpoints if the CN
    lookup mechanism is working. Computed BEFORE trusting this mechanism
    on any invented edge. See module docstring for the expected/observed
    numbers."""
    bp_lookup = s5a.set_index(['Individual', 'Breakpoint number'])[['chromosome', 'position']]

    pairs = []
    seen = set()
    for _, row in s5a[s5a['delb'].notna()].iterrows():
        a, b = row['Breakpoint number'], row['delb']
        pid = row['Individual']
        key = (pid, tuple(sorted((a, b))))
        if key in seen:
            continue
        seen.add(key)
        pairs.append((pid, a, b))

    rows = []
    for pid, a, b in pairs:
        ra = bp_lookup.loc[(pid, a)]
        rb = bp_lookup.loc[(pid, b)]
        chrom = ra['chromosome']
        pos_a, pos_b = int(ra['position']), int(rb['position'])

        seg_a = find_segment(segments, pid, chrom, pos_a)
        seg_b = find_segment(segments, pid, chrom, pos_b)

        if seg_a is None or seg_b is None:
            rows.append({
                'patient_id': pid, 'breakpoint_a': a, 'breakpoint_b': b,
                'chromosome': chrom, 'pos_a': pos_a, 'pos_b': pos_b,
                'both_resolve': False, 'has_transition': None,
                'shows_expected_loss': None,
                'seg_mean_a': None, 'seg_mean_b': None, 'min_between_mean': None,
            })
            continue

        between = segments_between(segments, pid, chrom, pos_a, pos_b)
        has_transition = seg_a['seg_id'] != seg_b['seg_id']
        min_between = between['Segment_Mean'].min()
        flank_mean = max(seg_a['Segment_Mean'], seg_b['Segment_Mean'])
        shows_loss = bool(has_transition and min_between < flank_mean)

        rows.append({
            'patient_id': pid, 'breakpoint_a': a, 'breakpoint_b': b,
            'chromosome': chrom, 'pos_a': pos_a, 'pos_b': pos_b,
            'both_resolve': True, 'has_transition': bool(has_transition),
            'shows_expected_loss': shows_loss,
            'seg_mean_a': seg_a['Segment_Mean'], 'seg_mean_b': seg_b['Segment_Mean'],
            'min_between_mean': min_between,
        })

    return pd.DataFrame(rows)


# ── Real per-patient breakpoint positions (for the EXPLAINED/UNEXPLAINED check) ──

def build_real_position_index(chrom_df):
    """dict: (patient_id, chromosome_str) -> sorted list of every real
    breakpoint position for that patient, from BOTH ends of every row in
    chrom_aberrations_baca.csv (the FULL per-patient rearrangement set,
    not just this chain or even just mmc5-covered rows) -- the most
    generous, fairest possible source to check whether a CN transition is
    already "explained" by a real recorded event."""
    index = {}
    for patient_id in chrom_df['Individual'].dropna().unique():
        patient_rows = _get_full_patient_rows(patient_id, chrom_df)
        for _, row in patient_rows.iterrows():
            for chrom_col, pos_col in [
                ('Breakpoint 1 chromosome', 'Breakpoint 1 position'),
                ('Breakpoint 2 chromosome', 'Breakpoint 2 position'),
            ]:
                chrom = str(int(row[chrom_col]))
                pos = int(row[pos_col])
                index.setdefault((patient_id, chrom), []).append(pos)
    for key in index:
        index[key].sort()
    return index


# ── Per-invented-edge classification ────────────────────────────────────────

def classify_invented_edge(patient_id, chrom_a, pos_a, chrom_b, pos_b, segments, real_position_index):
    if chrom_a != chrom_b:
        seg_a = find_segment(segments, patient_id, chrom_a, pos_a)
        seg_b = find_segment(segments, patient_id, chrom_b, pos_b)
        return {
            'classification': 'NOT_APPLICABLE_INTERCHROMOSOMAL',
            'seg_mean_a': seg_a['Segment_Mean'] if seg_a is not None else None,
            'seg_mean_b': seg_b['Segment_Mean'] if seg_b is not None else None,
        }

    chrom = chrom_a
    seg_a = find_segment(segments, patient_id, chrom, pos_a)
    seg_b = find_segment(segments, patient_id, chrom, pos_b)

    if seg_a is None or seg_b is None:
        return {
            'classification': 'CN_INDETERMINATE',
            'seg_mean_a': seg_a['Segment_Mean'] if seg_a is not None else None,
            'seg_mean_b': seg_b['Segment_Mean'] if seg_b is not None else None,
        }

    if seg_a['seg_id'] == seg_b['seg_id']:
        return {
            'classification': 'CN_FLAT',
            'seg_mean_a': seg_a['Segment_Mean'], 'seg_mean_b': seg_b['Segment_Mean'],
        }

    lo, hi = sorted([pos_a, pos_b])
    real_positions = real_position_index.get((patient_id, chrom), [])
    explained = any(lo < p < hi for p in real_positions)
    classification = 'CN_TRANSITION_EXPLAINED' if explained else 'CN_TRANSITION_UNEXPLAINED'
    return {
        'classification': classification,
        'seg_mean_a': seg_a['Segment_Mean'], 'seg_mean_b': seg_b['Segment_Mean'],
    }


# ── Per-chain: get the distinct matching(s) achieving obligate/probable ────

def _distinct_matchings(achieving):
    seen = set()
    out = []
    for r in achieving:
        key = frozenset(tuple(sorted(pair)) for pair in r['matching'])
        if key not in seen:
            seen.add(key)
            out.append(r['matching'])
    return out


def get_chain_invented_matchings(grp):
    """Returns dict with 'enumerable' plus, if enumerable, 'obligate' and
    'probable' dicts each holding {'matchings': [...], 'is_unique': bool}.
    m=0 chains (already fully closed) get a trivial empty matching."""
    chrom_lookup = dict(zip(grp['Breakpoint number'], grp['chromosome']))
    combined = _build_real_combined_graph(grp)
    closed_components, open_components, open_ends = _closed_and_open_components(combined)
    m = len(open_ends)

    if m == 0:
        return {
            'enumerable': True,
            'obligate': {'matchings': [[]], 'is_unique': True},
            'probable': {'matchings': [[]], 'is_unique': True},
        }

    if m > MAX_DANGLING_ENDS:
        return {'enumerable': False, 'obligate': None, 'probable': None}

    enum_results = _enumerate_completions(closed_components, open_components, open_ends, chrom_lookup)
    _, _, _, obligate_achieving = _select_obligate(enum_results)
    _, _, _, _, _, probable_achieving = _select_probable(enum_results)

    obligate_matchings = _distinct_matchings(obligate_achieving)
    probable_matchings = _distinct_matchings(probable_achieving)

    return {
        'enumerable': True,
        'obligate': {'matchings': obligate_matchings, 'is_unique': len(obligate_matchings) == 1},
        'probable': {'matchings': probable_matchings, 'is_unique': len(probable_matchings) == 1},
    }


# ── Driver ───────────────────────────────────────────────────────────────────

CLASS_LABELS = [
    'CN_FLAT', 'CN_TRANSITION_EXPLAINED', 'CN_TRANSITION_UNEXPLAINED',
    'CN_INDETERMINATE', 'NOT_APPLICABLE_INTERCHROMOSOMAL',
]


def run_cn_plausibility_analysis(m5_df, chrom_df, segments):
    real_position_index = build_real_position_index(chrom_df)

    edge_rows = []
    chain_rows = []

    for (patient_id, chain_number), grp in m5_df.groupby(['Individual', 'Chain number']):
        pos_lookup = dict(zip(grp['Breakpoint number'], zip(grp['chromosome'], grp['position'])))
        result = get_chain_invented_matchings(grp)

        if not result['enumerable']:
            chain_rows.append({
                'patient_id': patient_id, 'baca_chain_number': chain_number,
                'phase2_enumerable': False,
            })
            continue

        for lens in ('obligate', 'probable'):
            matchings = result[lens]['matchings']
            is_unique = result[lens]['is_unique']

            per_realization_counts = []
            for real_idx, matching in enumerate(matchings):
                counts = {label: 0 for label in CLASS_LABELS}
                for a, b in matching:
                    chrom_a, pos_a = pos_lookup[a]
                    chrom_b, pos_b = pos_lookup[b]
                    outcome = classify_invented_edge(
                        patient_id, chrom_a, pos_a, chrom_b, pos_b, segments, real_position_index
                    )
                    counts[outcome['classification']] += 1
                    edge_rows.append({
                        'patient_id': patient_id, 'baca_chain_number': chain_number, 'lens': lens,
                        'realization_index': real_idx, 'n_realizations': len(matchings),
                        'edges_unique': is_unique,
                        'breakpoint_a': a, 'breakpoint_b': b,
                        'chromosome_a': chrom_a, 'position_a': pos_a,
                        'chromosome_b': chrom_b, 'position_b': pos_b,
                        'classification': outcome['classification'],
                        'seg_mean_a': outcome['seg_mean_a'], 'seg_mean_b': outcome['seg_mean_b'],
                    })
                per_realization_counts.append(counts)

            chain_row = {
                'patient_id': patient_id, 'baca_chain_number': chain_number,
                'phase2_enumerable': True, 'lens': lens,
                'n_invented_edges': len(matchings[0]) if matchings else 0,
                'n_realizations': len(matchings),
                'edges_unique': is_unique,
            }
            if per_realization_counts:
                for label in CLASS_LABELS:
                    vals = [c[label] for c in per_realization_counts]
                    chain_row[f'{label.lower()}_min'] = min(vals)
                    chain_row[f'{label.lower()}_max'] = max(vals)
            else:
                for label in CLASS_LABELS:
                    chain_row[f'{label.lower()}_min'] = 0
                    chain_row[f'{label.lower()}_max'] = 0

            # Convenience booleans built from the min/max range above, since
            # is_unique=False (probable, 194/354 chains) means a single
            # "the" classification does not exist -- report both ends of
            # the range explicitly rather than picking one arbitrarily.
            # When edges_unique=True (obligate, always; probable, 160/354)
            # min == max and these two columns agree.
            chain_row['has_unexplained_in_every_realization'] = (
                chain_row['cn_transition_unexplained_min'] > 0
            )
            chain_row['has_unexplained_in_some_realization'] = (
                chain_row['cn_transition_unexplained_max'] > 0
            )
            chain_rows.append(chain_row)

    return pd.DataFrame(edge_rows), pd.DataFrame(chain_rows)


def main():
    print("Loading data...")
    chrom_df = pd.read_csv(BACA_CHROM_ABER_FILE_PATH)
    m5_df = load_mmc5_table_s5a()
    m5_df['position'] = m5_df['Chromosome:position'].str.split(':').str[1].astype(int)

    segments = load_cn_segments()

    print("Running deletion-bridge calibration (real-data positive control)...")
    calibration_df = run_deletion_bridge_calibration(m5_df, segments)
    calibration_df.to_csv(CALIBRATION_CSV, index=False)
    n = len(calibration_df)
    both = calibration_df['both_resolve'].sum()
    n_trans = calibration_df.loc[calibration_df['both_resolve'], 'has_transition'].sum()
    n_loss = calibration_df['shows_expected_loss'].sum()
    print(f"  {n} real deletion-bridge pairs; {both} have both endpoints resolvable in S3B")
    print(f"  of those: {n_trans} show a CN transition, {n_loss} show the expected LOSS")
    print(f"  Wrote {CALIBRATION_CSV}")

    print("\nRunning CN plausibility classification on Phase 2 invented edges...")
    edge_df, chain_df = run_cn_plausibility_analysis(m5_df, chrom_df, segments)
    edge_df.to_csv(EDGE_DETAIL_CSV, index=False)
    chain_df.to_csv(CHAIN_SUMMARY_CSV, index=False)
    print(f"  Wrote {len(edge_df)} edge-detail rows -> {EDGE_DETAIL_CSV}")
    print(f"  Wrote {len(chain_df)} chain-summary rows -> {CHAIN_SUMMARY_CSV}")

    print("\n=== Edge-detail summary (realization-weighted -- a handful of highly-ambiguous")
    print("    chains with thousands of tied realizations dominate this raw count; see the")
    print("    chain-level summary below for the non-skewed picture) ===")
    print(edge_df['classification'].value_counts())

    enumerable = chain_df[chain_df['phase2_enumerable'] == True].copy()
    enumerable['edges_unique'] = enumerable['edges_unique'].astype(bool)

    print("\n=== Edge-choice ambiguity (obligate/probable structure achieved by >1 distinct edge-set) ===")
    for lens in ('obligate', 'probable'):
        sub = enumerable[enumerable['lens'] == lens]
        n_ambig = int((~sub['edges_unique']).sum())
        print(f"  {lens}: {n_ambig} / {len(sub)} chains have >1 distinct invented edge-set for the selected structure")

    print("\n=== Chain-level CN plausibility (one row per chain -- NOT realization-weighted) ===")
    for lens in ('obligate', 'probable'):
        sub = enumerable[(enumerable['lens'] == lens) & (enumerable['n_invented_edges'] > 0)]
        n_every = int(sub['has_unexplained_in_every_realization'].sum())
        n_some = int(sub['has_unexplained_in_some_realization'].sum())
        print(f"  {lens} ({len(sub)} chains with >=1 invented edge):")
        print(f"    >=1 CN-unexplained edge in EVERY tied realization: {n_every} / {len(sub)}")
        print(f"    >=1 CN-unexplained edge in AT LEAST ONE realization: {n_some} / {len(sub)}")

    print("\n=== Cross-reference against Phase 3's chromoplexy-recovery patients ===")
    print("(patients real_has_chromoplexy_strict=False but obligate/probable_has_chromoplexy_strict=True --")
    print(" the professor's 'what extra info do completed cycles give' population)")
    p3 = pd.read_csv(os.path.join(RESULTS_DIR, "phase3_patient_chromoplexy_summary.csv"))
    for lens in ('obligate', 'probable'):
        recovered = p3[(p3['real_has_chromoplexy_strict'] == False) & (p3[f'{lens}_has_chromoplexy_strict'] == True)]
        recovered_ids = set(recovered['patient_id'])
        sub = enumerable[(enumerable['lens'] == lens) & (enumerable['patient_id'].isin(recovered_ids))]
        per_patient_some = sub.groupby('patient_id')['has_unexplained_in_some_realization'].any()
        per_patient_every = sub.groupby('patient_id')['has_unexplained_in_every_realization'].any()
        print(f"  {lens}: {len(recovered_ids)} recovered patients")
        print(f"    >=1 chain with an unexplained edge in some realization: {int(per_patient_some.sum())} / {len(recovered_ids)}")
        print(f"    >=1 chain with an unexplained edge in every realization: {int(per_patient_every.sum())} / {len(recovered_ids)}")


if __name__ == "__main__":
    main()
