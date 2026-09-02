"""
Stat test 3: does ETS fusion status (ETS+/ETS-) associate with any of our
chromoplexy calls, or with Baca's own criterion?

Full plain-language write-up in the companion file
test3_ets_status_comparison.md — read that first if you're new to this test.

Part A (PRIMARY, patient level, N=52): Fisher's exact test for ETS_status
vs. each of real/obligate/probable x strict/loose (6 tests) plus
baca_chain_has_chromoplexy (1 test) = 7 tests total, each with odds ratio +
95% CI (Haldane-Anscombe corrected when any cell is 0).

Part B (SECONDARY, chain level, N=366): reproduces the earlier
mmc5_chain_closure_summary.csv finding (ETS+ 58/141 vs ETS- 33/225 chains
with >=1 real closed cycle) directly from the Phase 2 CSV as a sanity
check, then extends the same chain-level split to
obligate_/probable_has_chromoplexy_strict/loose.

No multiple-testing correction is applied within this file — a combined
Benjamini-Hochberg pass across tests 2-4's Fisher tests is planned as a
separate follow-up step once test 4 is done (see the .md Limitations
section).

Outputs:
  results/stat_test3_ets_status_comparison.csv
  results/stat_test3_ets_status_patient_level.png
  results/stat_test3_ets_status_chain_level.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PHASE2_CSV = os.path.join(RESULTS_DIR, "phase2_obligate_probable_completion.csv")
PHASE3_CSV = os.path.join(RESULTS_DIR, "phase3_patient_chromoplexy_summary.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "stat_test3_ets_status_comparison.csv")
PATIENT_PNG = os.path.join(RESULTS_DIR, "stat_test3_ets_status_patient_level.png")
CHAIN_PNG = os.path.join(RESULTS_DIR, "stat_test3_ets_status_chain_level.png")

LENSES = ["real", "obligate", "probable"]
THRESHOLDS = ["strict", "loose"]


def odds_ratio_ci(a, b, c, d, alpha=0.05):
    """95% CI for the odds ratio via the log-odds normal approximation.
    Haldane-Anscombe correction (+0.5 to every cell) when any cell is 0."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    log_or = np.log((a * d) / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    z = 1.96
    return np.exp(log_or - z * se), np.exp(log_or + z * se)


def run_patient_level(p3):
    print("=== Part A (primary, patient level, N=52): ETS_status vs. each "
          "chromoplexy flag ===")
    rows = []
    flags = [(lens, thr, f"{lens}_has_chromoplexy_{thr}") for lens in LENSES for thr in THRESHOLDS]
    flags.append((None, None, "baca_chain_has_chromoplexy"))

    for lens, thr, col in flags:
        a = int(((p3["ETS_status"] == "ETS+") & (p3[col] == True)).sum())
        b = int(((p3["ETS_status"] == "ETS+") & (p3[col] == False)).sum())
        c = int(((p3["ETS_status"] == "ETS-") & (p3[col] == True)).sum())
        d = int(((p3["ETS_status"] == "ETS-") & (p3[col] == False)).sum())
        odds_ratio, p_val = fisher_exact([[a, b], [c, d]])
        ci_lo, ci_hi = odds_ratio_ci(a, b, c, d)
        label = col if lens is None else f"{lens}/{thr}"
        print(f"  {label:24s} ETS+ {a:2d}/{a+b:2d}  ETS- {c:2d}/{c+d:2d}  "
              f"OR={odds_ratio:8.3f}  95%CI=({ci_lo:.3f},{ci_hi:.3f})  p={p_val:.4f}")
        rows.append({
            "part": "patient_level",
            "level": "patient",
            "flag": col,
            "lens": lens,
            "threshold": thr,
            "n_ets_pos_flag_pos": a,
            "n_ets_pos_flag_neg": b,
            "n_ets_neg_flag_pos": c,
            "n_ets_neg_flag_neg": d,
            "odds_ratio": odds_ratio,
            "ci_95_low": ci_lo,
            "ci_95_high": ci_hi,
            "p_value": p_val,
        })
    return rows


def run_chain_level(p2):
    print("\n=== Part B (secondary, chain level, N=366): ETS_status vs. "
          "chain-level closure flags ===")
    rows = []
    chain_flags = [("real", None, "real_has_any_closed_cycle")]
    for lens in ["obligate", "probable"]:
        for thr in THRESHOLDS:
            chain_flags.append((lens, thr, f"{lens}_has_chromoplexy_{thr}"))

    for lens, thr, col in chain_flags:
        a = int(((p2["ETS_status"] == "ETS+") & (p2[col] == True)).sum())
        b = int(((p2["ETS_status"] == "ETS+") & (p2[col] == False)).sum())
        c = int(((p2["ETS_status"] == "ETS-") & (p2[col] == True)).sum())
        d = int(((p2["ETS_status"] == "ETS-") & (p2[col] == False)).sum())
        odds_ratio, p_val = fisher_exact([[a, b], [c, d]])
        ci_lo, ci_hi = odds_ratio_ci(a, b, c, d)
        label = col if thr is None else f"{lens}/{thr}"
        print(f"  {label:32s} ETS+ {a:3d}/{a+b:3d}  ETS- {c:3d}/{c+d:3d}  "
              f"OR={odds_ratio:8.3f}  95%CI=({ci_lo:.3f},{ci_hi:.3f})  p={p_val:.4f}")
        rows.append({
            "part": "chain_level",
            "level": "chain",
            "flag": col,
            "lens": lens,
            "threshold": thr,
            "n_ets_pos_flag_pos": a,
            "n_ets_pos_flag_neg": b,
            "n_ets_neg_flag_pos": c,
            "n_ets_neg_flag_neg": d,
            "odds_ratio": odds_ratio,
            "ci_95_low": ci_lo,
            "ci_95_high": ci_hi,
            "p_value": p_val,
        })
    return rows


def plot_patient_level(rows):
    fig, ax = plt.subplots(figsize=(9, 5))
    labels = [r["flag"] for r in rows]
    ets_pos_rate = [r["n_ets_pos_flag_pos"] / (r["n_ets_pos_flag_pos"] + r["n_ets_pos_flag_neg"]) for r in rows]
    ets_neg_rate = [r["n_ets_neg_flag_pos"] / (r["n_ets_neg_flag_pos"] + r["n_ets_neg_flag_neg"]) for r in rows]
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, ets_pos_rate, width, label="ETS+")
    ax.bar(x + width / 2, ets_neg_rate, width, label="ETS-")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=40, ha="right")
    ax.set_ylabel("fraction of patients positive")
    ax.set_title("Patient-level chromoplexy rate by ETS status (N=52)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(PATIENT_PNG, dpi=150)
    plt.close(fig)
    print(f"\nSaved {PATIENT_PNG}")


def plot_chain_level(rows):
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = [r["flag"] for r in rows]
    ets_pos_rate = [r["n_ets_pos_flag_pos"] / (r["n_ets_pos_flag_pos"] + r["n_ets_pos_flag_neg"]) for r in rows]
    ets_neg_rate = [r["n_ets_neg_flag_pos"] / (r["n_ets_neg_flag_pos"] + r["n_ets_neg_flag_neg"]) for r in rows]
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, ets_pos_rate, width, label="ETS+")
    ax.bar(x + width / 2, ets_neg_rate, width, label="ETS-")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("fraction of chains positive")
    ax.set_title("Chain-level closure rate by ETS status (N=366)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(CHAIN_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved {CHAIN_PNG}")


def main():
    p2 = pd.read_csv(PHASE2_CSV)
    p3 = pd.read_csv(PHASE3_CSV)
    print(f"Loaded {len(p3)} patients (phase3), {len(p2)} chains (phase2).\n")

    rows_a = run_patient_level(p3)
    rows_b = run_chain_level(p2)

    out_df = pd.DataFrame(rows_a + rows_b)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {len(out_df)} test rows to {OUTPUT_CSV}")

    plot_patient_level(rows_a)
    plot_chain_level(rows_b)


if __name__ == "__main__":
    main()
