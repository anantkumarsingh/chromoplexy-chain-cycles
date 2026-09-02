"""
Stat test 4 (ERG+ subset, N=23): does the older connectivity-based ERG
heuristic (chromoplexy_embedded / simple_fusion, from
results/erg_chain_comparison_ours_vs_baca.csv, populated for 23/52
mmc5-covered patients) agree with our real/obligate/probable closure-based
chromoplexy calls, or with Baca's own size-based criterion?

Full plain-language write-up in the companion file
test4_erg_heuristic_comparison.md — read that first if you're new to this
test.

For each of real/obligate/probable x strict/loose (6) plus
baca_chain_has_chromoplexy (1) = 7 tests: Fisher's exact test (small N=23
subset, several small cells) treating erg_fusion_type_heuristic ==
"chromoplexy_embedded" as the positive class, plus Cohen's kappa as a
chance-corrected agreement statistic against each flag.

No multiple-testing correction applied within this file — see the .md
Limitations section for the planned combined pass across tests 2-4.

Outputs:
  results/stat_test4_erg_heuristic_comparison.csv
  results/stat_test4_erg_heuristic_agreement.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PHASE3_CSV = os.path.join(RESULTS_DIR, "phase3_patient_chromoplexy_summary.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "stat_test4_erg_heuristic_comparison.csv")
PLOT_PNG = os.path.join(RESULTS_DIR, "stat_test4_erg_heuristic_agreement.png")

LENSES = ["real", "obligate", "probable"]
THRESHOLDS = ["strict", "loose"]


def cohens_kappa(a, b, c, d):
    n = a + b + c + d
    po = (a + d) / n
    p_row_pos = (a + b) / n
    p_col_pos = (a + c) / n
    pe = p_row_pos * p_col_pos + (1 - p_row_pos) * (1 - p_col_pos)
    if pe == 1:
        return float("nan")
    return (po - pe) / (1 - pe)


def run_tests(sub):
    rows = []
    flags = [(lens, thr, f"{lens}_has_chromoplexy_{thr}") for lens in LENSES for thr in THRESHOLDS]
    flags.append((None, None, "baca_chain_has_chromoplexy"))

    is_embedded = sub["erg_fusion_type_heuristic"] == "chromoplexy_embedded"
    for lens, thr, col in flags:
        a = int((is_embedded & (sub[col] == True)).sum())
        b = int((is_embedded & (sub[col] == False)).sum())
        c = int((~is_embedded & (sub[col] == True)).sum())
        d = int((~is_embedded & (sub[col] == False)).sum())
        odds_ratio, p_val = fisher_exact([[a, b], [c, d]])
        kappa = cohens_kappa(a, b, c, d)
        label = col if lens is None else f"{lens}/{thr}"
        print(f"  {label:28s} embedded {a:2d}/{a+b:2d}  simple {c:2d}/{c+d:2d}  "
              f"OR={odds_ratio:8.3f}  p={p_val:.4f}  kappa={kappa:+.3f}")
        rows.append({
            "flag": col,
            "lens": lens,
            "threshold": thr,
            "n_embedded_flag_pos": a,
            "n_embedded_flag_neg": b,
            "n_simple_flag_pos": c,
            "n_simple_flag_neg": d,
            "odds_ratio": odds_ratio,
            "p_value": p_val,
            "cohen_kappa": kappa,
        })
    return rows


def plot_agreement(rows):
    fig, ax = plt.subplots(figsize=(9, 5))
    labels = [r["flag"] for r in rows]
    kappas = [r["cohen_kappa"] for r in rows]
    colors = ["tab:green" if k >= 0 else "tab:red" for k in kappas]
    ax.bar(labels, kappas, color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Cohen's kappa (vs. ERG heuristic)")
    ax.set_title("Agreement between ERG connectivity heuristic and each chromoplexy flag (N=23)")
    plt.xticks(rotation=40, ha="right")
    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150)
    plt.close(fig)
    print(f"\nSaved {PLOT_PNG}")


def main():
    p3 = pd.read_csv(PHASE3_CSV)
    sub = p3[p3["erg_fusion_type_heuristic"].notna()].copy()
    print(f"Loaded {len(p3)} patients total; {len(sub)} have the ERG heuristic "
          f"populated ({(sub['erg_fusion_type_heuristic']=='chromoplexy_embedded').sum()} "
          f"chromoplexy_embedded, {(sub['erg_fusion_type_heuristic']=='simple_fusion').sum()} simple_fusion).\n")

    rows = run_tests(sub)
    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {len(out_df)} test rows to {OUTPUT_CSV}")

    plot_agreement(rows)


if __name__ == "__main__":
    main()
