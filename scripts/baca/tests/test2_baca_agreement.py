"""
Stat test 2 (patient level, N=52): formal agreement between Baca's own
chromoplexy criterion (baca_chain_has_chromoplexy, >=1 chain with >=5
rearrangements) and each of our real/obligate/probable x strict/loose
chromoplexy calls.

Full plain-language write-up in the companion file
test2_baca_agreement.md — read that first if you're new to this test.

Part A — six 2x2 contingency tables (baca_chain_has_chromoplexy vs.
<lens>_has_chromoplexy_<threshold>, for lens in {real, obligate, probable}
and threshold in {strict, loose}): Fisher's exact test (small/zero cell
counts rule out the chi-square approximation) plus Cohen's kappa as a
chance-corrected agreement statistic.

Part B — dose-response check: Spearman rho between n_chains_baca_threshold
(how many of a patient's chains meet Baca's >=5-rearrangement criterion —
an ordinal proxy for "how strongly Baca-positive") and the largest
chromosome span any of that patient's chains actually achieved under each
lens (real/obligate/probable), aggregated from the chain-level Phase 2 CSV.

Outputs:
  results/stat_test2_baca_agreement.csv
  results/stat_test2_baca_agreement_contingency.png
  results/stat_test2_baca_agreement_kappa.png
  results/stat_test2_baca_agreement_dose_response.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
PHASE2_CSV = os.path.join(RESULTS_DIR, "phase2_obligate_probable_completion.csv")
PHASE3_CSV = os.path.join(RESULTS_DIR, "phase3_patient_chromoplexy_summary.csv")
OUTPUT_CSV = os.path.join(RESULTS_DIR, "stat_test2_baca_agreement.csv")
CONTINGENCY_PNG = os.path.join(RESULTS_DIR, "stat_test2_baca_agreement_contingency.png")
KAPPA_PNG = os.path.join(RESULTS_DIR, "stat_test2_baca_agreement_kappa.png")
DOSE_RESPONSE_PNG = os.path.join(RESULTS_DIR, "stat_test2_baca_agreement_dose_response.png")

LENSES = ["real", "obligate", "probable"]
THRESHOLDS = ["strict", "loose"]


def cohens_kappa(a, b, c, d):
    """2x2 table: a=baca+/ours+, b=baca+/ours-, c=baca-/ours+, d=baca-/ours-."""
    n = a + b + c + d
    po = (a + d) / n
    p_baca_pos = (a + b) / n
    p_ours_pos = (a + c) / n
    pe = p_baca_pos * p_ours_pos + (1 - p_baca_pos) * (1 - p_ours_pos)
    if pe == 1:
        return float("nan")
    return (po - pe) / (1 - pe)


def run_contingency_tests(p3):
    rows = []
    tables = {}
    for lens in LENSES:
        for thr in THRESHOLDS:
            col = f"{lens}_has_chromoplexy_{thr}"
            a = int(((p3["baca_chain_has_chromoplexy"] == True) & (p3[col] == True)).sum())
            b = int(((p3["baca_chain_has_chromoplexy"] == True) & (p3[col] == False)).sum())
            c = int(((p3["baca_chain_has_chromoplexy"] == False) & (p3[col] == True)).sum())
            d = int(((p3["baca_chain_has_chromoplexy"] == False) & (p3[col] == False)).sum())
            odds_ratio, p_val = fisher_exact([[a, b], [c, d]])
            kappa = cohens_kappa(a, b, c, d)
            print(f"  {lens:9s} {thr:6s}  a={a:2d} b={b:2d} c={c:2d} d={d:2d}  "
                  f"OR={odds_ratio:.3f}  p={p_val:.4f}  kappa={kappa:+.3f}")
            rows.append({
                "test": "baca_agreement_2x2",
                "threshold": thr,
                "lens": lens,
                "n_baca_pos_ours_pos": a,
                "n_baca_pos_ours_neg": b,
                "n_baca_neg_ours_pos": c,
                "n_baca_neg_ours_neg": d,
                "odds_ratio": odds_ratio,
                "p_value": p_val,
                "cohen_kappa": kappa,
                "spearman_rho": None,
                "n_patients": None,
            })
            tables[(lens, thr)] = (a, b, c, d)
    return rows, tables


def run_dose_response(p2, p3):
    print("\n=== Dose-response: n_chains_baca_threshold vs. best achieved "
          "chromosome span per patient ===")
    rows = []
    span_by_patient = {}
    for lens in LENSES:
        span_col = f"{lens}_max_cycle_chrom_span"
        agg = p2.groupby("patient_id")[span_col].max()
        span_by_patient[lens] = agg
        merged = p3[["patient_id", "n_chains_baca_threshold"]].merge(
            agg.rename("max_span"), left_on="patient_id", right_index=True, how="inner"
        ).dropna()
        rho, p_val = spearmanr(merged["n_chains_baca_threshold"], merged["max_span"])
        print(f"  {lens:9s}  rho={rho:+.3f}  p={p_val:.4f}  (n={len(merged)})")
        rows.append({
            "test": "dose_response_correlation",
            "threshold": None,
            "lens": lens,
            "n_baca_pos_ours_pos": None,
            "n_baca_pos_ours_neg": None,
            "n_baca_neg_ours_pos": None,
            "n_baca_neg_ours_neg": None,
            "odds_ratio": None,
            "p_value": p_val,
            "cohen_kappa": None,
            "spearman_rho": rho,
            "n_patients": len(merged),
        })
    return rows, span_by_patient


def plot_contingency(tables):
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for row_i, thr in enumerate(THRESHOLDS):
        for col_i, lens in enumerate(LENSES):
            ax = axes[row_i, col_i]
            a, b, c, d = tables[(lens, thr)]
            grid = np.array([[a, b], [c, d]])
            ax.imshow(grid, cmap="Blues", vmin=0, vmax=max(grid.max(), 1))
            for (i, j), val in np.ndenumerate(grid):
                ax.text(j, i, str(val), ha="center", va="center", fontsize=14,
                        color="white" if val > grid.max() / 2 else "black")
            ax.set_xticks([0, 1])
            ax.set_yticks([0, 1])
            ax.set_xticklabels(["ours +", "ours -"])
            ax.set_yticklabels(["baca +", "baca -"])
            ax.set_title(f"{lens} / {thr}")
    fig.suptitle("Baca vs. ours: 2x2 chromoplexy agreement (N=52 patients)")
    fig.tight_layout()
    fig.savefig(CONTINGENCY_PNG, dpi=150)
    plt.close(fig)
    print(f"\nSaved {CONTINGENCY_PNG}")


def plot_kappa(rows):
    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.arange(len(LENSES))
    width = 0.35
    for i, thr in enumerate(THRESHOLDS):
        vals = [r["cohen_kappa"] for r in rows if r["threshold"] == thr]
        ax.bar(x + (i - 0.5) * width, vals, width, label=thr)
    ax.set_xticks(x)
    ax.set_xticklabels(LENSES)
    ax.set_ylabel("Cohen's kappa (vs. Baca)")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_title("Chance-corrected agreement with Baca's criterion")
    ax.legend(title="threshold")
    fig.tight_layout()
    fig.savefig(KAPPA_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved {KAPPA_PNG}")


def plot_dose_response(p3, span_by_patient):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True)
    for ax, lens in zip(axes, LENSES):
        merged = p3[["patient_id", "n_chains_baca_threshold"]].merge(
            span_by_patient[lens].rename("max_span"), left_on="patient_id",
            right_index=True, how="inner"
        ).dropna()
        jitter = np.random.default_rng(0).uniform(-0.1, 0.1, size=len(merged))
        ax.scatter(merged["n_chains_baca_threshold"] + jitter, merged["max_span"], alpha=0.6)
        ax.set_xlabel("n_chains_baca_threshold")
        ax.set_ylabel(f"{lens}_max_cycle_chrom_span (best per patient)")
        ax.set_title(lens)
    fig.suptitle("Dose-response: Baca chain-size strength vs. best achieved cycle span")
    fig.tight_layout()
    fig.savefig(DOSE_RESPONSE_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved {DOSE_RESPONSE_PNG}")


def main():
    p2 = pd.read_csv(PHASE2_CSV)
    p3 = pd.read_csv(PHASE3_CSV)
    print(f"Loaded {len(p3)} patients (phase3), {len(p2)} chains (phase2).")

    print("\n=== Part A: 2x2 contingency tests (Fisher's exact + Cohen's kappa) ===")
    rows_a, tables = run_contingency_tests(p3)
    rows_b, span_by_patient = run_dose_response(p2, p3)

    out_df = pd.DataFrame(rows_a + rows_b)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {len(out_df)} test rows to {OUTPUT_CSV}")

    plot_contingency(tables)
    plot_kappa(rows_a)
    plot_dose_response(p3, span_by_patient)


if __name__ == "__main__":
    main()
