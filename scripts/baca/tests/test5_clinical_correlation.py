"""
Stat test 5 (patient level, N=52): does Gleason grade or pathological
T-stage correlate with any of our chromoplexy calls, or with Baca's
criterion?

Full plain-language write-up in the companion file
test5_clinical_correlation.md — read that first if you're new to this
test, especially for the encoding decisions (agreed with the user
2026-08-05, see the .md file section 1).

Encoding:
  - Gleason_Score -> base ISUP Grade Group from primary+secondary pattern
    ONLY (tertiary pattern, where recorded, does NOT change the base
    Grade Group number) + a separate has_tertiary_pattern boolean flag
    (and the raw tertiary grade, for reference). 2 patients with NaN
    Gleason_Score are excluded pairwise from any test using Grade Group.
  - Pathological_stage -> ordinal T-stage (pT2a < pT2b < pT2c < pT3a <
    pT3b < pT4) with any " N1" suffix stripped into a separate
    has_nodal_involvement flag. The 2 "Metastatic" patients have no
    T-stage at all and are EXCLUDED from the T-stage ordinal test (per
    explicit user instruction: "leave them be, we will analyze them
    later") — flagged as is_metastatic=True and reported separately, not
    silently dropped.

For each of real/obligate/probable x strict/loose (6) plus
baca_chain_has_chromoplexy (1) = 7 flags, and each of the 2 ordinal
clinical variables (Grade Group, T-stage): Mann-Whitney U comparing the
clinical variable between the flag's positive and negative groups, plus
rank-biserial effect size. 14 tests total. Spearman correlation was
considered and deliberately omitted as redundant with the rank-biserial
effect size already reported for a binary-vs-ordinal comparison (see the
.md Limitations section).

Outputs:
  results/stat_test5_clinical_correlation.csv
  results/stat_test5_gleason_grade_group.png
  results/stat_test5_pathological_t_stage.png
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PHASE3_CSV = os.path.join(RESULTS_DIR, "phase3_patient_chromoplexy_summary.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "stat_test5_clinical_correlation.csv")
GLEASON_PNG = os.path.join(RESULTS_DIR, "stat_test5_gleason_grade_group.png")
STAGE_PNG = os.path.join(RESULTS_DIR, "stat_test5_pathological_t_stage.png")

LENSES = ["real", "obligate", "probable"]
THRESHOLDS = ["strict", "loose"]

GRADE_GROUP = {
    (3, 3): 1,
    (3, 4): 2,
    (4, 3): 3,
    (4, 4): 4, (3, 5): 4, (5, 3): 4,
    (4, 5): 5, (5, 4): 5, (5, 5): 5,
}

T_STAGE_ORDER = {"pT2a": 1, "pT2b": 2, "pT2c": 3, "pT3a": 4, "pT3b": 5, "pT4": 6}


def encode_gleason(raw):
    if pd.isna(raw):
        return pd.Series({"grade_group": None, "has_tertiary_pattern": None, "tertiary_grade": None})
    base, has_tertiary, tertiary_grade = raw, False, None
    if ";" in raw:
        base, tertiary_str = raw.split(";")
        has_tertiary = True
        tertiary_grade = int(tertiary_str)
    primary, secondary = (int(x) for x in base.split("+"))
    grade_group = GRADE_GROUP[(primary, secondary)]
    return pd.Series({
        "grade_group": grade_group,
        "has_tertiary_pattern": has_tertiary,
        "tertiary_grade": tertiary_grade,
    })


def encode_stage(raw):
    if pd.isna(raw):
        return pd.Series({"t_stage_ordinal": None, "has_nodal_involvement": None, "is_metastatic": None})
    if raw == "Metastatic":
        return pd.Series({"t_stage_ordinal": None, "has_nodal_involvement": False, "is_metastatic": True})
    has_nodal = bool(re.search(r"N1", raw))
    t_code = re.sub(r"\s*N1\s*", "", raw).strip()
    return pd.Series({
        "t_stage_ordinal": T_STAGE_ORDER[t_code],
        "has_nodal_involvement": has_nodal,
        "is_metastatic": False,
    })


def rank_biserial(u_stat, n1, n2):
    return 1 - (2 * u_stat) / (n1 * n2)


def run_tests(p3, clinical_var, var_label):
    rows = []
    flags = [(lens, thr, f"{lens}_has_chromoplexy_{thr}") for lens in LENSES for thr in THRESHOLDS]
    flags.append((None, None, "baca_chain_has_chromoplexy"))

    for lens, thr, col in flags:
        sub = p3[["patient_id", col, clinical_var]].dropna()
        pos = sub[sub[col] == True][clinical_var].values
        neg = sub[sub[col] == False][clinical_var].values
        n1, n2 = len(pos), len(neg)
        label = col if lens is None else f"{lens}/{thr}"
        if n1 == 0 or n2 == 0:
            print(f"  [{var_label}] {label:28s} skipped (one group empty)")
            continue
        u_stat, p_val = mannwhitneyu(pos, neg, alternative="two-sided")
        effect = rank_biserial(u_stat, n1, n2)
        med_pos, med_neg = pd.Series(pos).median(), pd.Series(neg).median()
        print(f"  [{var_label}] {label:28s} n_pos={n1:2d} n_neg={n2:2d}  "
              f"median_pos={med_pos:.1f}  median_neg={med_neg:.1f}  "
              f"U={u_stat:6.1f}  p={p_val:.4f}  r={effect:+.3f}")
        rows.append({
            "clinical_variable": var_label,
            "flag": col,
            "lens": lens,
            "threshold": thr,
            "n_flag_pos": n1,
            "n_flag_neg": n2,
            "median_flag_pos": med_pos,
            "median_flag_neg": med_neg,
            "U_statistic": u_stat,
            "p_value": p_val,
            "rank_biserial_r": effect,
        })
    return rows


def plot_by_flag(p3, clinical_var, ylabel, out_path):
    flags = [f"{lens}_has_chromoplexy_{thr}" for lens in LENSES for thr in THRESHOLDS]
    flags.append("baca_chain_has_chromoplexy")
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    for ax, col in zip(axes.flat, flags):
        sub = p3[[col, clinical_var]].dropna()
        pos = sub[sub[col] == True][clinical_var].values
        neg = sub[sub[col] == False][clinical_var].values
        ax.boxplot([pos, neg], tick_labels=["+", "-"])
        ax.set_title(col, fontsize=9)
        ax.set_ylabel(ylabel)
    for ax in axes.flat[len(flags):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    p3 = pd.read_csv(PHASE3_CSV)
    print(f"Loaded {len(p3)} patients.\n")

    p3 = pd.concat([p3, p3["Gleason_Score"].apply(encode_gleason)], axis=1)
    p3 = pd.concat([p3, p3["Pathological_stage"].apply(encode_stage)], axis=1)

    n_tertiary = p3["has_tertiary_pattern"].sum()
    n_nodal = p3["has_nodal_involvement"].sum()
    n_meta = p3["is_metastatic"].sum()
    print(f"Encoding notes: {int(n_tertiary)} patients have a recorded tertiary "
          f"Gleason pattern (flagged, not folded into grade_group); "
          f"{int(n_nodal)} have nodal involvement (N1) noted; "
          f"{int(n_meta)} are Metastatic-stage (excluded from the T-stage "
          f"ordinal test, set aside per explicit instruction, not analyzed here).\n")

    print("=== Gleason Grade Group vs. chromoplexy flags ===")
    rows_gleason = run_tests(p3, "grade_group", "gleason_grade_group")
    print("\n=== Pathological T-stage vs. chromoplexy flags ===")
    rows_stage = run_tests(p3, "t_stage_ordinal", "pathological_t_stage")

    out_df = pd.DataFrame(rows_gleason + rows_stage)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {len(out_df)} test rows to {OUTPUT_CSV}")

    plot_by_flag(p3, "grade_group", "Gleason Grade Group", GLEASON_PNG)
    plot_by_flag(p3, "t_stage_ordinal", "Pathological T-stage (ordinal)", STAGE_PNG)


if __name__ == "__main__":
    main()
