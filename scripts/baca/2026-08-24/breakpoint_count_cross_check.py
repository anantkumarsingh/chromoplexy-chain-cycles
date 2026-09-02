"""
Cross-checks mmc5.xlsx Table S5B's 'Total breakpoints' column against a
directly-computed whole-genome breakpoint count, for each of the 57
patients.

  A. total_bp_whole_genome  — 2 x number of rearrangement rows for that
     patient in chrom_aberrations_baca.csv (== mmc3 Table S3C), i.e.
     every breakpoint anywhere in the patient's genome, chained or not.
  B. total_bp_mmc5_reported — Table S5B's own 'Total breakpoints' column,
     copied as-is.

An earlier pass also computed a third count -- distinct breakpoints
across all of a patient's mmc5 chains combined (Table S5A) -- to test
whether S5B's 'Total breakpoints' was actually a chain-scoped count.
Result: it never matched (0/57), while the whole-genome count matches
40/57 exactly. That ruled out the chain-scoped hypothesis, so this
script now only tracks the whole-genome vs. mmc5-reported comparison,
which is where the real (17/57 patient) disagreement lives.

Output: results/2026-08-24/breakpoint_count_cross_check.csv
"""

import pandas as pd
import os

BACA_DATASET_FOLDER = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/baca_dataset"
RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results/2026-08-24"

CLINICAL_PATH = f"{BACA_DATASET_FOLDER}/clinical_phenotypes.csv"
CHROM_PATH = f"{BACA_DATASET_FOLDER}/chrom_aberrations_baca.csv"
MMC5_PATH = f"{BACA_DATASET_FOLDER}/mmc5.xlsx"


def load_patients():
    clinical = pd.read_csv(CLINICAL_PATH)
    return sorted(clinical['Individual'].unique().tolist())


def compute_whole_genome_totals(patients):
    chrom = pd.read_csv(CHROM_PATH)
    chrom = chrom[chrom['Individual'].isin(set(patients))]
    counts = chrom.groupby('Individual').size() * 2
    return counts.to_dict()


def load_mmc5_reported_totals(patients):
    s5b = pd.read_excel(MMC5_PATH, sheet_name='Table S5B')
    s5b = s5b[s5b['Individual'].isin(set(patients))]
    return dict(zip(s5b['Individual'], s5b['Total breakpoints']))


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    patients = load_patients()

    whole_genome = compute_whole_genome_totals(patients)
    mmc5_reported = load_mmc5_reported_totals(patients)

    rows = []
    for p in patients:
        a = int(whole_genome.get(p, 0))
        c = mmc5_reported.get(p)
        c = int(c) if pd.notna(c) else None

        rows.append({
            'patientID': p,
            'total_bp_whole_genome': a,
            'total_bp_mmc5_reported': c,
            'agree': (a == c),
            'delta_whole_genome_minus_mmc5': a - c if c is not None else None,
        })

    df = pd.DataFrame(rows)
    out_path = f"{RESULTS_DIR}/breakpoint_count_cross_check.csv"
    df.to_csv(out_path, index=False)

    n = len(df)
    n_agree = df['agree'].sum()
    n_disagree = n - n_agree
    disagree_df = df[~df['agree']].copy()
    n_mmc5_higher = (disagree_df['delta_whole_genome_minus_mmc5'] < 0).sum()
    n_raw_higher = (disagree_df['delta_whole_genome_minus_mmc5'] > 0).sum()

    print(f"Wrote {out_path} — {n} patients\n")
    print(f"Agree (whole_genome == mmc5_reported)   : {n_agree}/{n}")
    print(f"Disagree                                : {n_disagree}/{n}")
    print(f"  of which mmc5 > whole_genome          : {n_mmc5_higher}")
    print(f"  of which whole_genome > mmc5          : {n_raw_higher}")

    print("\n--- Disagreeing patients ---")
    print(disagree_df.sort_values('delta_whole_genome_minus_mmc5')[
        ['patientID', 'total_bp_whole_genome', 'total_bp_mmc5_reported', 'delta_whole_genome_minus_mmc5']
    ].to_string(index=False))


if __name__ == '__main__':
    main()
