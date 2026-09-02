"""
Stat test 1 (chain level, N=354 enumerable chains out of 366): does chain
size predict whether a chain qualifies for chromoplexy under the
obligate/probable completion rules, among chains that do NOT already
qualify from real (observed) data alone?

Full plain-language write-up (what/why/how/results/limitations) lives in
the companion file test1_chain_size_confound.md in this same folder — read
that first if you're new to this test. This docstring only covers the
mechanics.

Population for each comparison: chains with real_has_chromoplexy_<thr> ==
False AND phase2_enumerable == True (the 12 non-enumerable chains have no
obligate/probable columns to compare, so they're excluded, not imputed).

Group "qualifies"        : <lens>_has_chromoplexy_<thr> == True
Group "does not qualify" : <lens>_has_chromoplexy_<thr> == False

For each of {obligate, probable} x {strict, loose} (4 comparisons), and each
size variable in {n_dangling_ends, n_rearrangements, k_chromosomes_whole_chain}:
Mann-Whitney U test (two-sided) comparing the "qualifies" group vs the "does
not qualify" group on that variable, plus rank-biserial correlation as the
effect size.

Companion check: Spearman rho between n_dangling_ends and
obligate_/probable_max_cycle_chrom_span across all 354 enumerable chains.

Output: results/stat_test1_chain_size_confound.csv (plain, one row per test).
"""
import os
import sys

import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PHASE2_CSV = os.path.join(RESULTS_DIR, "phase2_obligate_probable_completion.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "stat_test1_chain_size_confound.csv")

SIZE_VARS = ["n_dangling_ends", "n_rearrangements", "k_chromosomes_whole_chain"]
LENSES = ["obligate", "probable"]
THRESHOLDS = ["strict", "loose"]


def rank_biserial(u_stat, n1, n2):
    """Effect size for Mann-Whitney U: r in [-1, 1], 0 = no effect."""
    return 1 - (2 * u_stat) / (n1 * n2)


def load_enumerable_chains():
    df = pd.read_csv(PHASE2_CSV)
    n_total = len(df)
    enum_df = df[df["phase2_enumerable"] == True].copy()
    print(f"Loaded {n_total} chains total; {len(enum_df)} enumerable "
          f"(phase2_enumerable=True), {n_total - len(enum_df)} excluded "
          f"(too many dangling ends to enumerate).")
    return enum_df


def run_qualification_tests(enum_df):
    rows = []
    for thr in THRESHOLDS:
        real_col = f"real_has_chromoplexy_{thr}"
        population = enum_df[enum_df[real_col] == False]
        print(f"\n=== Threshold: {thr} | chains NOT already chromoplexy-positive "
              f"from real data: {len(population)} / {len(enum_df)} enumerable chains ===")

        for lens in LENSES:
            lens_col = f"{lens}_has_chromoplexy_{thr}"
            qualifies = population[population[lens_col] == True]
            not_qualifies = population[population[lens_col] == False]
            n1, n2 = len(qualifies), len(not_qualifies)
            print(f"\n  {lens} lens: now qualifies={n1}  still does not qualify={n2}")

            if n1 == 0 or n2 == 0:
                print("    skipped (one group empty)")
                continue

            for var in SIZE_VARS:
                a = qualifies[var].values
                b = not_qualifies[var].values
                u_stat, p_val = mannwhitneyu(a, b, alternative="two-sided")
                effect = rank_biserial(u_stat, n1, n2)
                med_a, med_b = pd.Series(a).median(), pd.Series(b).median()
                print(f"    {var:28s} median(qualifies)={med_a:5.1f}  "
                      f"median(does not qualify)={med_b:5.1f}  "
                      f"U={u_stat:8.1f}  p={p_val:.4f}  "
                      f"rank-biserial r={effect:+.3f}")
                rows.append({
                    "test": "qualification_confound",
                    "threshold": thr,
                    "lens": lens,
                    "variable": var,
                    "n_qualifying": n1,
                    "n_not_qualifying": n2,
                    "median_qualifying": med_a,
                    "median_not_qualifying": med_b,
                    "U_statistic": u_stat,
                    "p_value": p_val,
                    "rank_biserial_r": effect,
                })
    return rows


def run_span_correlation(enum_df):
    print("\n=== Spearman correlation: n_dangling_ends vs max achievable "
          f"cycle chromosome span (N={len(enum_df)} enumerable chains) ===")
    rows = []
    for lens in LENSES:
        span_col = f"{lens}_max_cycle_chrom_span"
        sub = enum_df[["n_dangling_ends", span_col]].dropna()
        rho, p_val = spearmanr(sub["n_dangling_ends"], sub[span_col])
        print(f"  {lens}: rho={rho:+.3f}  p={p_val:.4f}  (n={len(sub)})")
        rows.append({
            "test": "span_correlation",
            "threshold": None,
            "lens": lens,
            "variable": "n_dangling_ends_vs_" + span_col,
            "n_qualifying": None,
            "n_not_qualifying": None,
            "median_qualifying": None,
            "median_not_qualifying": None,
            "U_statistic": None,
            "p_value": p_val,
            "rank_biserial_r": None,
            "spearman_rho": rho,
        })
    return rows


def main():
    enum_df = load_enumerable_chains()
    rows = run_qualification_tests(enum_df)
    rows += run_span_correlation(enum_df)
    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {len(out_df)} test rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
