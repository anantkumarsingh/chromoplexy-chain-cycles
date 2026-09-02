"""
Stat test 6 (chain level, N=354 enumerable chains out of 366): how reliable
/ how confident is the Phase 2 obligate/probable completion machinery
itself, independent of any chromoplexy question?

Full plain-language write-up in the companion file
test6_method_reliability.md — read that first if you're new to this test.

Three diagnostics, all descriptive (no group comparison, no hypothesis
about chromoplexy):

1. Distribution of probable_pct (what share of all enumerated matchings
   produced the "most probable" structure that got reported) — how
   dominant is the mode, on average?
2. How often probable_is_tied=True (no unique mode) and
   obligate_max_span_is_unique=False (the lexicographic obligate rule hit
   a tie) actually occur.
3. Spearman correlation between n_dangling_ends (chain complexity) and
   both probable_pct and probable_n_distinct_structures — does completion
   confidence degrade, and does the space of possible structures grow, as
   chains get more complex?

Output: results/stat_test6_method_reliability.csv (a flat metric/value
table, not repeated per lens/threshold like the other tests — this test
diagnoses the method, not a chromoplexy hypothesis).

Images:
  results/stat_test6_probable_pct_distribution.png
  results/stat_test6_probable_pct_vs_dangling_ends.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PHASE2_CSV = os.path.join(RESULTS_DIR, "phase2_obligate_probable_completion.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "stat_test6_method_reliability.csv")
DIST_PNG = os.path.join(RESULTS_DIR, "stat_test6_probable_pct_distribution.png")
SCATTER_PNG = os.path.join(RESULTS_DIR, "stat_test6_probable_pct_vs_dangling_ends.png")


def main():
    p2 = pd.read_csv(PHASE2_CSV)
    enum_df = p2[p2["phase2_enumerable"] == True].copy()
    print(f"Loaded {len(p2)} chains total; {len(enum_df)} enumerable.\n")

    rows = []

    print("=== probable_pct distribution ===")
    desc = enum_df["probable_pct"].describe()
    for stat in ["mean", "std", "min", "25%", "50%", "75%", "max"]:
        print(f"  {stat:6s} {desc[stat]:.2f}")
        rows.append({"metric": f"probable_pct_{stat.replace('%', 'pct')}", "value": desc[stat]})

    print("\n=== Tie/uniqueness frequency ===")
    n_tied = int((enum_df["probable_is_tied"] == True).sum())
    n_not_tied = int((enum_df["probable_is_tied"] == False).sum())
    n_unique = int((enum_df["obligate_max_span_is_unique"] == True).sum())
    n_not_unique = int((enum_df["obligate_max_span_is_unique"] == False).sum())
    print(f"  probable_is_tied=True            : {n_tied} / {len(enum_df)}")
    print(f"  probable_is_tied=False           : {n_not_tied} / {len(enum_df)}")
    print(f"  obligate_max_span_is_unique=True : {n_unique} / {len(enum_df)}")
    print(f"  obligate_max_span_is_unique=False: {n_not_unique} / {len(enum_df)}")
    rows += [
        {"metric": "n_probable_is_tied_true", "value": n_tied},
        {"metric": "n_probable_is_tied_false", "value": n_not_tied},
        {"metric": "n_obligate_max_span_unique_true", "value": n_unique},
        {"metric": "n_obligate_max_span_unique_false", "value": n_not_unique},
    ]

    print("\n=== Spearman correlations vs. n_dangling_ends ===")
    rho1, p1 = spearmanr(enum_df["n_dangling_ends"], enum_df["probable_pct"])
    print(f"  probable_pct               rho={rho1:+.3f}  p={p1:.3e}")
    rho2, p2v = spearmanr(enum_df["n_dangling_ends"], enum_df["probable_n_distinct_structures"])
    print(f"  probable_n_distinct_structures  rho={rho2:+.3f}  p={p2v:.3e}")
    rows += [
        {"metric": "spearman_rho_probable_pct_vs_dangling_ends", "value": rho1},
        {"metric": "spearman_p_probable_pct_vs_dangling_ends", "value": p1},
        {"metric": "spearman_rho_n_distinct_structures_vs_dangling_ends", "value": rho2},
        {"metric": "spearman_p_n_distinct_structures_vs_dangling_ends", "value": p2v},
    ]

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {len(out_df)} metric rows to {OUTPUT_CSV}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(enum_df["probable_pct"], bins=20, edgecolor="black")
    ax.set_xlabel("probable_pct (% of matchings producing the reported structure)")
    ax.set_ylabel("number of chains")
    ax.set_title(f"Distribution of probable_pct (N={len(enum_df)} enumerable chains)")
    fig.tight_layout()
    fig.savefig(DIST_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved {DIST_PNG}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(enum_df["n_dangling_ends"], enum_df["probable_pct"], alpha=0.5)
    ax.set_xlabel("n_dangling_ends")
    ax.set_ylabel("probable_pct")
    ax.set_title(f"probable_pct vs. n_dangling_ends (Spearman rho={rho1:+.3f})")
    fig.tight_layout()
    fig.savefig(SCATTER_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved {SCATTER_PNG}")


if __name__ == "__main__":
    main()
