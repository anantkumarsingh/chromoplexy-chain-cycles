"""
baca_phase3_patient_aggregation.py — Phase 3 per-patient chromoplexy rollup.

phase2_obligate_probable_completion.csv is per-CHAIN (366 rows, 52 patients
covered by mmc5). Baca's own headline stat ("50 of 57 tumors show
chromoplexy") is per-TUMOR. This script rolls chain-level flags up to
patient level via an OR rule: a patient counts as positive on a given
metric if AT LEAST ONE of their chains is positive -- this mirrors Baca's
own framing exactly (a tumor needs only one qualifying chain).

Two independent chromoplexy criteria are aggregated, kept as separate
columns, never collapsed into one "final" call:

1. baca_chain_has_chromoplexy -- Baca's OWN published criterion, based
   purely on chain SIZE (>=5 rearrangements/chain), no cycle closure
   required at all. Reproduces his headline stat, quoted directly:
   "Baca et al. inferred chromoplexy in 50 of 57 prostate tumors based
   on the presence of at least one ChainFinder-defined, statistically
   interdependent chain containing five or more rearrangements."
   NOTE: only 52/57 Baca patients have any mmc5 chain at all, so this
   script's baca_chain_has_chromoplexy count will land somewhat below
   50/57 -- the 5 patients absent from mmc5 entirely cannot be scored
   here, not a disagreement with Baca, just a coverage gap.

   TODO (deferred, not yet computed): Baca's SECOND, distinct threshold
   -- chains of >=3 rearrangements where the chain also disrupts a
   KEGG-listed prostate TSG (26/57 tumors, 46%) -- is NOT the same as a
   plain >=3-rearrangement size cutoff (a plain size-only cutoff
   saturates to ~52/52 in this cohort and isn't a meaningful
   comparison). Revisit once we decide how to check TSG disruption
   against the `genes_in_chain` column.

2. real_/obligate_/probable_has_chromoplexy_strict/loose -- OUR
   closure-based criterion (a closed AMG cycle spanning >=3 / >=2 distinct
   chromosomes), computed at three lenses:
     - real       = only what closes from Baca's own recorded edges
                    (rearrangement + adjacency + deletion bridge). The
                    only one of the three that supports a "confirmed"
                    chromoplexy claim (Cornforth's irreducible/reducible
                    distinction).
     - obligate   = real cycle(s) preserved as-is, PLUS the smallest
                    possible invented closure of whatever's still open
                    (least coordination assumed).
     - probable   = real cycle(s) preserved as-is, PLUS the most
                    frequent invented closure among all combinatorially
                    valid completions.
   strict and loose are NESTED (strict implies loose), not a true/false
   partition -- both are always reported, never collapsed to one flag.

obligate/probable depend on enumeration (phase2_enumerable), which fails
for 12/366 chains (>12 dangling ends). A patient whose ONLY negative
evidence for a given obligate/probable flag includes one of these
non-enumerable chains gets a companion
*_has_unresolved_caveat_strict/loose = True, so a "False" aggregated flag
is never silently conflated with "actually unknown."

Outputs (all in RESULTS_DIR):
  phase3_patient_chromoplexy_summary.csv   -- all 52 mmc5-covered patients,
                                               every flag above + counts.
  phase3_baca_positive_vs_ours.csv         -- THE final table: subset of
                                               patients positive on
                                               baca_chain_has_chromoplexy
                                               (>=5 rearr, Baca's headline
                                               50/57 population), with all
                                               OUR cycle-chromosome-span
                                               criteria columns alongside
                                               for direct agree/disagree
                                               inspection.
"""

import os

import pandas as pd

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PHASE2_CSV = os.path.join(RESULTS_DIR, "phase2_obligate_probable_completion.csv")
ERG_HEURISTIC_CSV = os.path.join(RESULTS_DIR, "erg_chain_comparison_ours_vs_baca.csv")

SUMMARY_CSV = os.path.join(RESULTS_DIR, "phase3_patient_chromoplexy_summary.csv")
BACA_POSITIVE_TABLE_CSV = os.path.join(RESULTS_DIR, "phase3_baca_positive_vs_ours.csv")

BACA_MIN_REARR = 5


def _or_ignore_nan(series):
    """OR across a column that may contain True/False/NaN (NaN from
    non-enumerable chains) -- NaN never contributes a positive."""
    return bool(series.fillna(False).any())


def aggregate_patients(phase2_df):
    records = []
    for patient_id, g in phase2_df.groupby("patient_id", sort=True):
        clinical_cols = ["ETS_status", "Gleason_Score", "Pathological_stage"]
        for col in clinical_cols:
            # nunique(dropna=True): 0 means all-NaN for this patient (e.g.
            # metastatic samples with no primary Gleason score) -- allowed,
            # not an inconsistency. >1 means genuinely conflicting values.
            assert g[col].nunique(dropna=True) <= 1, (
                f"{patient_id}: inconsistent {col} across chains -- "
                f"expected identical clinical values per patient"
            )

        n_chains_total = len(g)
        n_chains_enumerable = int(g["phase2_enumerable"].sum())
        n_chains_unresolved = n_chains_total - n_chains_enumerable
        has_unresolved_chain = n_chains_unresolved > 0

        n_chains_baca_threshold = int((g["n_rearrangements"] >= BACA_MIN_REARR).sum())
        baca_chain_has_chromoplexy = n_chains_baca_threshold > 0

        real_strict = _or_ignore_nan(g["real_has_chromoplexy_strict"])
        real_loose = _or_ignore_nan(g["real_has_chromoplexy_loose"])

        obligate_strict = _or_ignore_nan(g["obligate_has_chromoplexy_strict"])
        obligate_loose = _or_ignore_nan(g["obligate_has_chromoplexy_loose"])
        obligate_caveat_strict = (not obligate_strict) and has_unresolved_chain
        obligate_caveat_loose = (not obligate_loose) and has_unresolved_chain

        probable_strict = _or_ignore_nan(g["probable_has_chromoplexy_strict"])
        probable_loose = _or_ignore_nan(g["probable_has_chromoplexy_loose"])
        probable_caveat_strict = (not probable_strict) and has_unresolved_chain
        probable_caveat_loose = (not probable_loose) and has_unresolved_chain

        records.append(dict(
            patient_id=patient_id,
            ETS_status=g["ETS_status"].iloc[0],
            Gleason_Score=g["Gleason_Score"].iloc[0],
            Pathological_stage=g["Pathological_stage"].iloc[0],
            n_chains_total=n_chains_total,
            n_chains_enumerable=n_chains_enumerable,
            n_chains_unresolved=n_chains_unresolved,
            n_chains_baca_threshold=n_chains_baca_threshold,
            baca_chain_has_chromoplexy=baca_chain_has_chromoplexy,
            real_has_chromoplexy_strict=real_strict,
            real_has_chromoplexy_loose=real_loose,
            obligate_has_chromoplexy_strict=obligate_strict,
            obligate_has_chromoplexy_loose=obligate_loose,
            obligate_has_unresolved_caveat_strict=obligate_caveat_strict,
            obligate_has_unresolved_caveat_loose=obligate_caveat_loose,
            probable_has_chromoplexy_strict=probable_strict,
            probable_has_chromoplexy_loose=probable_loose,
            probable_has_unresolved_caveat_strict=probable_caveat_strict,
            probable_has_unresolved_caveat_loose=probable_caveat_loose,
        ))

    out = pd.DataFrame.from_records(records)

    for method in ("real", "obligate", "probable"):
        out[f"{method}_agrees_with_baca_strict"] = (
            out[f"{method}_has_chromoplexy_strict"] == out["baca_chain_has_chromoplexy"]
        )
        out[f"{method}_agrees_with_baca_loose"] = (
            out[f"{method}_has_chromoplexy_loose"] == out["baca_chain_has_chromoplexy"]
        )

    if os.path.exists(ERG_HEURISTIC_CSV):
        erg = pd.read_csv(ERG_HEURISTIC_CSV)[["patient_id", "baca_has_ERG_chain", "our_fusion_type"]]
        erg = erg.rename(columns={
            "baca_has_ERG_chain": "has_ERG_chain_in_mmc5",
            "our_fusion_type": "erg_fusion_type_heuristic",
        })
        out = out.merge(erg, on="patient_id", how="left")

    return out


def print_summary(out):
    n = len(out)
    print(f"Patients covered (>=1 mmc5 chain): {n}")
    print()
    print("Positive counts:")
    for col in [
        "baca_chain_has_chromoplexy",
        "real_has_chromoplexy_strict", "real_has_chromoplexy_loose",
        "obligate_has_chromoplexy_strict", "obligate_has_chromoplexy_loose",
        "probable_has_chromoplexy_strict", "probable_has_chromoplexy_loose",
    ]:
        print(f"  {col}: {int(out[col].sum())}/{n}")
    print()
    print("Unresolved-caveat counts (negative flag, but >=1 non-enumerable chain contributed):")
    for col in [
        "obligate_has_unresolved_caveat_strict", "obligate_has_unresolved_caveat_loose",
        "probable_has_unresolved_caveat_strict", "probable_has_unresolved_caveat_loose",
    ]:
        print(f"  {col}: {int(out[col].sum())}/{n}")
    print()

    pos = out[out["baca_chain_has_chromoplexy"]]
    print(f"=== Baca-positive patients (>=5 rearr, his headline 50/57 stat): {len(pos)} ===")
    for method in ("real", "obligate", "probable"):
        for thresh in ("strict", "loose"):
            flag_col = f"{method}_has_chromoplexy_{thresh}"
            agree_col = f"{method}_agrees_with_baca_{thresh}"
            n_ours_positive = int(pos[flag_col].sum())
            n_agree = int(pos[agree_col].sum())
            print(
                f"  {method:9s} {thresh:6s}: {n_ours_positive}/{len(pos)} also positive "
                f"({n_agree}/{len(pos)} agree with Baca)"
            )
    print()


def main():
    phase2_df = pd.read_csv(PHASE2_CSV)
    out = aggregate_patients(phase2_df)

    out.to_csv(SUMMARY_CSV, index=False)
    print(f"Wrote {SUMMARY_CSV} ({len(out)} patients)")

    baca_positive_table = out[out["baca_chain_has_chromoplexy"]].copy()
    baca_positive_table.to_csv(BACA_POSITIVE_TABLE_CSV, index=False)
    print(f"Wrote {BACA_POSITIVE_TABLE_CSV} ({len(baca_positive_table)} patients)")

    print()
    print_summary(out)

    return out


if __name__ == "__main__":
    main()
