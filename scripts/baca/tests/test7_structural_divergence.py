"""
Stat test 7 (chain level, N=354 enumerable chains out of 366): beyond the
strict/loose chromoplexy flags (already covered in test 1), how often do
obligate and probable completion pick a DIFFERENT cycle structure string
for the same chain, and does chain complexity predict that divergence?

Full plain-language write-up in the companion file
test7_structural_divergence.md — read that first if you're new to this
test.

Divergent chain: obligate_cycle_structure != probable_cycle_structure
(e.g. one says "C4+C2", the other says "C6" for the same chain).

Same method as test 1: Mann-Whitney U comparing n_dangling_ends /
n_rearrangements / k_chromosomes_whole_chain between divergent and
concordant chains, plus rank-biserial effect size.

Output:
  results/stat_test7_structural_divergence.csv
  results/stat_test7_divergence_by_size.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PHASE2_CSV = os.path.join(RESULTS_DIR, "phase2_obligate_probable_completion.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "stat_test7_structural_divergence.csv")
PLOT_PNG = os.path.join(RESULTS_DIR, "stat_test7_divergence_by_size.png")

SIZE_VARS = ["n_dangling_ends", "n_rearrangements", "k_chromosomes_whole_chain"]


def rank_biserial(u_stat, n1, n2):
    return 1 - (2 * u_stat) / (n1 * n2)


def main():
    p2 = pd.read_csv(PHASE2_CSV)
    enum_df = p2[p2["phase2_enumerable"] == True].copy()
    enum_df["divergent"] = enum_df["obligate_cycle_structure"] != enum_df["probable_cycle_structure"]

    n_divergent = int(enum_df["divergent"].sum())
    n_concordant = len(enum_df) - n_divergent
    print(f"Loaded {len(p2)} chains total; {len(enum_df)} enumerable.")
    print(f"Divergent (obligate != probable structure): {n_divergent} / {len(enum_df)}")
    print(f"Concordant (obligate == probable structure): {n_concordant} / {len(enum_df)}\n")

    rows = [{"metric": "n_divergent", "value": n_divergent},
            {"metric": "n_concordant", "value": n_concordant},
            {"metric": "n_enumerable_total", "value": len(enum_df)}]

    divergent = enum_df[enum_df["divergent"]]
    concordant = enum_df[~enum_df["divergent"]]
    n1, n2 = len(divergent), len(concordant)

    print("=== Mann-Whitney U: chain size, divergent vs. concordant ===")
    for var in SIZE_VARS:
        a = divergent[var].values
        b = concordant[var].values
        u_stat, p_val = mannwhitneyu(a, b, alternative="two-sided")
        effect = rank_biserial(u_stat, n1, n2)
        med_a, med_b = pd.Series(a).median(), pd.Series(b).median()
        print(f"  {var:28s} median(divergent)={med_a:5.1f}  median(concordant)={med_b:5.1f}  "
              f"U={u_stat:8.1f}  p={p_val:.4f}  r={effect:+.3f}")
        rows.append({"metric": f"{var}_median_divergent", "value": med_a})
        rows.append({"metric": f"{var}_median_concordant", "value": med_b})
        rows.append({"metric": f"{var}_mannwhitney_U", "value": u_stat})
        rows.append({"metric": f"{var}_p_value", "value": p_val})
        rows.append({"metric": f"{var}_rank_biserial_r", "value": effect})

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {len(out_df)} metric rows to {OUTPUT_CSV}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, var in zip(axes, SIZE_VARS):
        ax.boxplot([divergent[var].values, concordant[var].values],
                   tick_labels=["divergent", "concordant"])
        ax.set_title(var)
        ax.set_ylabel(var)
    fig.suptitle(f"Chain size vs. obligate/probable structural divergence "
                 f"({n_divergent} divergent / {n_concordant} concordant)")
    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved {PLOT_PNG}")


if __name__ == "__main__":
    main()
