"""
baca_chain_comparison.py — compare OUR ERG chain construction (built from the
raw rearrangement file, chrom_aberrations_baca.csv) against BACA'S OWN chain
assignment (mmc5.xlsx, Table S5A, "Chain number" field).

Background:
  chrom_aberrations_baca.csv is the full raw rearrangement file — every
  detected rearrangement per patient, no grouping. Our existing code
  (_get_erg_chain_rows / get_erg_chain_details_v2 in
  baca_aberration_clinical_analysis.py) builds an ERG-anchored chain from
  this file ourselves: anchor on ERG/TMPRSS2 site annotations, then expand
  one hop via inter_chr rows sharing a chromosome with the anchor.

  mmc5.xlsx (Table S5A) is a DIFFERENT, smaller file: one row per
  breakpoint (not per rearrangement), with a "Chain number" column that is
  BACA'S OWN clustering of rearrangements into coordinated chromoplexy
  chains. It only covers 52 of the 57 patients, and only the subset of
  each patient's rearrangements that got assigned to a chain (verified:
  for P01-28, 138/140 mmc5 breakpoints match exact positions in
  chrom_aberrations_baca.csv, but mmc5 only includes 70 of the patient's
  96 rearrangements — the rest are unclustered singletons in Baca's data).

This script computes, for every chain in BOTH sources, the same set of
descriptive stats (chromosomes, k, b_i, Sheth Θ(k,(b_i)) notation, genes,
ERG/TMPRSS2 presence) so the two can be compared side by side. No
visualization here — that's a deferred next phase.

Outputs (results/):
  baca_proposed_chains_mmc5.csv         — every chain in mmc5, all patients
  our_erg_chain_chrom_aberrations.csv   — our ERG chain, ERG+ patients only
  erg_chain_comparison_ours_vs_baca.csv — per-patient side-by-side comparison
"""

import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'core'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'phase1'))
from baca_aberration_clinical_analysis import (
    load_data,
    get_erg_chain_details_v2,
    _get_erg_chain_rows,
    load_ets_status_from_mmc5,
    BACA_CHROM_ABER_FILE_PATH,
    BACA_CLINICAL_PHENO_FILE_PATH,
)
from baca_full_genome_cycles import _get_full_patient_rows

BACA_DATASET_FOLDER = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/baca_dataset"
MMC5_FILE_PATH = os.path.join(BACA_DATASET_FOLDER, "mmc5.xlsx")

RESULTS_DIR = "/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/results"
BACA_CHAINS_CSV = os.path.join(RESULTS_DIR, "baca_proposed_chains_mmc5.csv")
BACA_CHR21_CHAINS_CSV = os.path.join(RESULTS_DIR, "baca_erg_chr21_chains_mmc5.csv")
OUR_ERG_CHAIN_CSV = os.path.join(RESULTS_DIR, "our_erg_chain_chrom_aberrations.csv")
OUR_FULL_CHAIN_CSV = os.path.join(RESULTS_DIR, "our_full_chain_chrom_aberrations.csv")
COMPARISON_CSV = os.path.join(RESULTS_DIR, "erg_chain_comparison_ours_vs_baca.csv")
ERG_SET_COMPARISON_CSV = os.path.join(RESULTS_DIR, "erg_set_comparison_baca_chr21_vs_ours.csv")

GENE_PATTERN = re.compile(r'([A-Z][A-Z0-9\-]+)\([+-]\)')
ERG_PATTERN = re.compile(r'ERG|TMPRSS2')


def extract_genes(site_strings):
    genes = set()
    for s in site_strings:
        if pd.isna(s):
            continue
        genes.update(GENE_PATTERN.findall(str(s)))
    return sorted(genes)


def theta_notation(b_i):
    """Θ(k,(b_1,...,b_k)) sorted by chromosome for reproducibility."""
    chroms = sorted(b_i.keys())
    k = len(chroms)
    b_str = ",".join(str(b_i[c]) for c in chroms)
    return f"Theta({k},({b_str}))"


# ── BACA'S OWN CHAINS (mmc5.xlsx) ───────────────────────────────────────────

def load_mmc5_chains():
    df = pd.ExcelFile(MMC5_FILE_PATH).parse('Table S5A')
    df.columns = [c.strip() for c in df.columns]
    df = df[df['Chromosome:position'].notna()].copy()
    df['chromosome'] = df['Chromosome:position'].str.split(':').str[0].str.replace('chr', '', regex=False)
    df['position'] = df['Chromosome:position'].str.split(':').str[1].astype(float)
    return df


def compute_baca_chains(m5_df, clinical_df):
    rows = []
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    for (patient_id, chain_number), grp in m5_df.groupby(['Individual', 'Chain number']):
        chroms = grp['chromosome'].tolist()
        b_i = grp.groupby('chromosome').size().to_dict()
        site_col = grp['Site annotation'] if 'Site annotation' in grp.columns else pd.Series([], dtype=str)
        genes = extract_genes(site_col.tolist())
        contains_erg = bool(site_col.dropna().str.contains(ERG_PATTERN).any())

        n_breakpoints = len(grp)
        n_rearrangements = grp['Rearrangement number'].nunique() if 'Rearrangement number' in grp.columns else None

        end_loss_left = grp['Potential chromosomal end loss (left)'].sum() if 'Potential chromosomal end loss (left)' in grp.columns else None
        end_loss_right = grp['Potential chromosomal end loss (right)'].sum() if 'Potential chromosomal end loss (right)' in grp.columns else None
        n_deletion_bridges = grp['Deletion bridge partner breakpoint'].notna().sum() if 'Deletion bridge partner breakpoint' in grp.columns else None
        total_reads = grp['Supporting reads'].sum() if 'Supporting reads' in grp.columns else None
        mean_reads = grp['Supporting reads'].mean() if 'Supporting reads' in grp.columns else None

        clin = clinical_df[clinical_df['Individual'] == patient_id]
        ets_status = ets_lookup.get(patient_id)
        gleason = clin['Gleason Score'].values[0] if len(clin) > 0 else None
        stage = clin['Pathological stage'].values[0] if len(clin) > 0 else None

        rows.append({
            'patient_id': patient_id,
            'baca_chain_number': chain_number,
            'n_breakpoints': n_breakpoints,
            'n_rearrangements': n_rearrangements,
            'k': len(set(chroms)),
            'chromosomes': ",".join(sorted(set(chroms), key=lambda x: (len(x), x))),
            'b_i': str(b_i),
            'theta': theta_notation(b_i),
            'contains_ERG_or_TMPRSS2': contains_erg,
            'genes_in_chain': ",".join(genes),
            'total_supporting_reads': total_reads,
            'mean_supporting_reads': mean_reads,
            'n_potential_end_loss_left': end_loss_left,
            'n_potential_end_loss_right': end_loss_right,
            'n_deletion_bridge_partners': n_deletion_bridges,
            'ETS_status': ets_status,
            'Gleason_Score': gleason,
            'Pathological_stage': stage,
        })

    out = pd.DataFrame(rows).sort_values(['patient_id', 'baca_chain_number'])
    return out


# ── BACA'S ERG CHAINS BY CHR21 MEMBERSHIP (broader than gene-annotation match) ──
#
# contains_ERG_or_TMPRSS2 (used in build_comparison) only flags a chain if
# the Site annotation text literally mentions "ERG" or "TMPRSS2". Since both
# genes physically sit on chr21, that flag is always a SUBSET of "chr21 is
# one of the chromosomes in this chain" — a chain can include chr21 breaks
# that land elsewhere on the chromosome (no ERG/TMPRSS2 gene annotation at
# that specific breakpoint) and still legitimately be "the ERG-region chain"
# for that patient. This selects on chromosome membership directly.

def select_chr21_chains(baca_chains_df):
    has_chr21 = baca_chains_df['chromosomes'].apply(
        lambda s: '21' in str(s).split(',') if pd.notna(s) else False
    )
    return baca_chains_df[has_chr21].copy().sort_values(['patient_id', 'baca_chain_number'])


# ── OUR ERG CHAIN (chrom_aberrations_baca.csv) ──────────────────────────────

def get_erg_positive_patients(clinical_df):
    return clinical_df[
        clinical_df['ETS fusion detected by sequencing'].str.contains('ERG', na=False)
    ]['Individual'].tolist()


def compute_our_erg_chains(chrom_df, clinical_df):
    rows = []
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    for patient_id in get_erg_positive_patients(clinical_df):
        details = get_erg_chain_details_v2(patient_id, chrom_df, clinical_df)
        chain_rows = _get_erg_chain_rows(patient_id, chrom_df, clinical_df)

        if len(chain_rows) == 0:
            b_i = {}
        else:
            chr1 = chain_rows['Breakpoint 1 chromosome'].dropna().astype(int).astype(str)
            chr2 = chain_rows['Breakpoint 2 chromosome'].dropna().astype(int).astype(str)
            b_i = pd.concat([chr1, chr2]).value_counts().to_dict()

        genes = extract_genes(
            list(chain_rows.get('Breakpoint 1 site', pd.Series(dtype=str)).dropna()) +
            list(chain_rows.get('Breakpoint 2 site', pd.Series(dtype=str)).dropna())
        ) if len(chain_rows) > 0 else []

        fusion_type = (
            'chromoplexy_embedded'
            if details['k_local'] > 2 and details['chain_size'] > 1
            else 'simple_fusion'
        )

        clin = clinical_df[clinical_df['Individual'] == patient_id]
        ets_status = ets_lookup.get(patient_id)
        gleason = clin['Gleason Score'].values[0] if len(clin) > 0 else None
        stage = clin['Pathological stage'].values[0] if len(clin) > 0 else None

        rows.append({
            'patient_id': patient_id,
            'k_local': details['k_local'],
            'chain_size': details['chain_size'],
            'n_DSBs': details['chain_size'] * 2,
            'chromosomes': ",".join(sorted(set(b_i.keys()), key=lambda x: (len(x), x))),
            'b_i': str(b_i),
            'theta': theta_notation(b_i) if b_i else 'Theta(0,())',
            'genes_in_chain': ",".join(genes),
            'source': details['source'],
            'fusion_type': fusion_type,
            'ETS_status': ets_status,
            'Gleason_Score': gleason,
            'Pathological_stage': stage,
        })

    return pd.DataFrame(rows).sort_values('patient_id')


# ── OUR FULL GENOME CHAIN (chrom_aberrations_baca.csv, no ERG anchor) ──────
#
# Unlike compute_our_erg_chains, this does not restrict to ERG/TMPRSS2
# anchoring or a one-hop expansion — it takes EVERY rearrangement row that
# patient has anywhere in chrom_aberrations_baca.csv. This is the same row
# set used by baca_full_genome_cycles.py's cycle computation, just summarized
# here as a flat chromosome/Θ-notation table for all 57 patients, not just
# the 26 ERG+ ones.

def compute_our_full_chains(chrom_df, clinical_df):
    rows = []
    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    for patient_id in chrom_df['Individual'].dropna().unique():
        patient_rows = _get_full_patient_rows(patient_id, chrom_df)
        if len(patient_rows) == 0:
            continue

        chr1 = patient_rows['Breakpoint 1 chromosome'].dropna().astype(int).astype(str)
        chr2 = patient_rows['Breakpoint 2 chromosome'].dropna().astype(int).astype(str)
        b_i = pd.concat([chr1, chr2]).value_counts().to_dict()

        genes = extract_genes(
            list(patient_rows.get('Breakpoint 1 site', pd.Series(dtype=str)).dropna()) +
            list(patient_rows.get('Breakpoint 2 site', pd.Series(dtype=str)).dropna())
        )

        n_inter = (patient_rows['Class'] == 'inter_chr').sum()
        n_intra = len(patient_rows) - n_inter

        contains_erg = bool(
            patient_rows['Breakpoint 1 site'].dropna().str.contains(ERG_PATTERN).any() or
            patient_rows['Breakpoint 2 site'].dropna().str.contains(ERG_PATTERN).any() or
            patient_rows['Fusion'].dropna().str.contains(ERG_PATTERN).any()
        )

        clin = clinical_df[clinical_df['Individual'] == patient_id]
        ets_status = ets_lookup.get(patient_id)
        gleason = clin['Gleason Score'].values[0] if len(clin) > 0 else None
        stage = clin['Pathological stage'].values[0] if len(clin) > 0 else None

        rows.append({
            'patient_id': patient_id,
            'k': len(b_i),
            'chromosomes': ",".join(sorted(b_i.keys(), key=lambda x: (len(x), x))),
            'b_i': str(b_i),
            'theta': theta_notation(b_i),
            'n_rearrangements': len(patient_rows),
            'n_DSBs': len(patient_rows) * 2,
            'n_inter_chr': n_inter,
            'n_intra_chr': n_intra,
            'pct_inter': round(100 * n_inter / len(patient_rows), 1),
            'max_b_i': max(b_i.values()),
            'dominant_chromosome': max(b_i, key=b_i.get),
            'n_genes_in_chain': len(genes),
            'genes_in_chain': ",".join(genes),
            'contains_ERG_or_TMPRSS2': contains_erg,
            'ETS_status': ets_status,
            'Gleason_Score': gleason,
            'Pathological_stage': stage,
        })

    return pd.DataFrame(rows).sort_values('patient_id')


# ── ERG+ PATIENT SET COMPARISON: Baca's chr21-chain set vs. our clinical set ──
#
# "ERG+ per Baca" here = patient has at least one chain in mmc5 that includes
# chr21 (select_chr21_chains). "ERG+ per ours" = clinical_phenotypes.csv's
# 'ETS fusion detected by sequencing' column mentions ERG. These are
# different definitions and can disagree in both directions:
#   - a patient can have a chr21 chain in mmc5 for a reason unrelated to
#     ERG/TMPRSS2 (chr21 has other genes) while NOT being clinically ERG+
#   - a patient can be clinically ERG+ but have no chain in mmc5 touching
#     chr21 at all (either missing from mmc5, or chr21 breaks stayed
#     unclustered/singleton and never made it into any chain)

def build_erg_set_comparison(baca_chr21_chains_df, our_erg_chain_df, clinical_df):
    baca_patients = set(baca_chr21_chains_df['patient_id'].unique())
    our_patients = set(our_erg_chain_df['patient_id'].unique())
    all_patients = sorted(baca_patients | our_patients)

    ets_lookup = load_ets_status_from_mmc5(known_patients=clinical_df['Individual'])
    rows = []
    for patient_id in all_patients:
        in_baca = patient_id in baca_patients
        in_ours = patient_id in our_patients

        if in_baca and in_ours:
            match_status = 'both'
        elif in_baca:
            match_status = 'baca_only'
        else:
            match_status = 'our_only'

        ets_status = ets_lookup.get(patient_id)

        # Baca side: union of chromosomes + chain numbers across all of this
        # patient's chr21-containing chains (handles the 3 patients with 2)
        baca_rows = baca_chr21_chains_df[baca_chr21_chains_df['patient_id'] == patient_id]
        if len(baca_rows) > 0:
            baca_chain_numbers = ";".join(str(c) for c in baca_rows['baca_chain_number'])
            baca_chrom_set = set()
            for chrom_str in baca_rows['chromosomes']:
                baca_chrom_set.update(str(chrom_str).split(','))
            baca_chromosomes = ",".join(sorted(baca_chrom_set, key=lambda x: (len(x), x)))
            baca_k = len(baca_chrom_set)
        else:
            baca_chain_numbers = None
            baca_chromosomes = None
            baca_k = None

        # Raw side: our ERG chain built from chrom_aberrations_baca.csv
        our_row = our_erg_chain_df[our_erg_chain_df['patient_id'] == patient_id]
        if len(our_row) > 0:
            r = our_row.iloc[0]
            raw_chromosomes = r['chromosomes']
            raw_k_local = r['k_local']
            raw_theta = r['theta']
        else:
            raw_chromosomes = None
            raw_k_local = None
            raw_theta = None

        chromosome_sets_match = None
        if baca_chromosomes is not None and raw_chromosomes is not None:
            chromosome_sets_match = set(baca_chromosomes.split(',')) == set(raw_chromosomes.split(','))

        rows.append({
            'patient_id': patient_id,
            'ETS_status': ets_status,
            'in_baca_chr21_chain_set': in_baca,
            'in_our_clinical_ERG_set': in_ours,
            'match_status': match_status,
            'baca_chain_number(s)': baca_chain_numbers,
            'baca_k': baca_k,
            'baca_chromosomes': baca_chromosomes,
            'raw_k_local': raw_k_local,
            'raw_chromosomes': raw_chromosomes,
            'raw_theta': raw_theta,
            'chromosome_sets_match': chromosome_sets_match,
        })

    return pd.DataFrame(rows).sort_values(['match_status', 'patient_id'])


# ── COMPARISON ───────────────────────────────────────────────────────────

def build_comparison(baca_chains_df, our_chains_df):
    erg_chains = baca_chains_df[baca_chains_df['contains_ERG_or_TMPRSS2']].copy()
    # In case a patient somehow has >1 chain flagged with ERG/TMPRSS2, keep the largest (by k)
    erg_chains = erg_chains.sort_values('k', ascending=False).drop_duplicates('patient_id')

    rows = []
    for _, our in our_chains_df.iterrows():
        patient_id = our['patient_id']
        baca_match = erg_chains[erg_chains['patient_id'] == patient_id]
        in_mmc5 = patient_id in baca_chains_df['patient_id'].values

        if len(baca_match) == 0:
            rows.append({
                'patient_id': patient_id,
                'in_mmc5': in_mmc5,
                'baca_has_ERG_chain': False,
                'baca_chain_number': None,
                'baca_k': None,
                'baca_chromosomes': None,
                'baca_theta': None,
                'baca_n_rearrangements': None,
                'our_k_local': our['k_local'],
                'our_chromosomes': our['chromosomes'],
                'our_theta': our['theta'],
                'our_chain_size': our['chain_size'],
                'our_fusion_type': our['fusion_type'],
                'k_match': None,
                'chromosome_sets_match': None,
            })
        else:
            b = baca_match.iloc[0]
            our_chr_set = set(our['chromosomes'].split(',')) if our['chromosomes'] else set()
            baca_chr_set = set(b['chromosomes'].split(',')) if b['chromosomes'] else set()
            rows.append({
                'patient_id': patient_id,
                'in_mmc5': in_mmc5,
                'baca_has_ERG_chain': True,
                'baca_chain_number': b['baca_chain_number'],
                'baca_k': b['k'],
                'baca_chromosomes': b['chromosomes'],
                'baca_theta': b['theta'],
                'baca_n_rearrangements': b['n_rearrangements'],
                'our_k_local': our['k_local'],
                'our_chromosomes': our['chromosomes'],
                'our_theta': our['theta'],
                'our_chain_size': our['chain_size'],
                'our_fusion_type': our['fusion_type'],
                'k_match': b['k'] == our['k_local'],
                'chromosome_sets_match': our_chr_set == baca_chr_set,
            })

    return pd.DataFrame(rows).sort_values('patient_id')


def main():
    chrom_df, clinical_df = load_data(BACA_CHROM_ABER_FILE_PATH, BACA_CLINICAL_PHENO_FILE_PATH)
    m5_df = load_mmc5_chains()

    baca_chains_df = compute_baca_chains(m5_df, clinical_df)
    baca_chains_df.to_csv(BACA_CHAINS_CSV, index=False)
    print(f"Wrote {len(baca_chains_df)} Baca chains ({baca_chains_df['patient_id'].nunique()} patients) -> {BACA_CHAINS_CSV}")

    baca_chr21_chains_df = select_chr21_chains(baca_chains_df)
    baca_chr21_chains_df.to_csv(BACA_CHR21_CHAINS_CSV, index=False)
    print(f"Wrote {len(baca_chr21_chains_df)} Baca chains containing chr21 ({baca_chr21_chains_df['patient_id'].nunique()} patients) -> {BACA_CHR21_CHAINS_CSV}")

    our_chains_df = compute_our_erg_chains(chrom_df, clinical_df)
    our_chains_df.to_csv(OUR_ERG_CHAIN_CSV, index=False)
    print(f"Wrote {len(our_chains_df)} of our ERG chains -> {OUR_ERG_CHAIN_CSV}")

    our_full_chains_df = compute_our_full_chains(chrom_df, clinical_df)
    our_full_chains_df.to_csv(OUR_FULL_CHAIN_CSV, index=False)
    print(f"Wrote {len(our_full_chains_df)} full-genome chains (all patients) -> {OUR_FULL_CHAIN_CSV}")

    comparison_df = build_comparison(baca_chains_df, our_chains_df)
    comparison_df.to_csv(COMPARISON_CSV, index=False)
    print(f"Wrote {len(comparison_df)}-row comparison -> {COMPARISON_CSV}")

    erg_set_comparison_df = build_erg_set_comparison(baca_chr21_chains_df, our_chains_df, clinical_df)
    erg_set_comparison_df.to_csv(ERG_SET_COMPARISON_CSV, index=False)
    print(f"Wrote {len(erg_set_comparison_df)}-row ERG-set comparison -> {ERG_SET_COMPARISON_CSV}")
    print(erg_set_comparison_df['match_status'].value_counts().to_string())

    print()
    print("Comparison summary:")
    print(f"  ERG+ patients (our set)         : {len(our_chains_df)}")
    print(f"  ...present in mmc5              : {comparison_df['in_mmc5'].sum()}")
    print(f"  ...with a Baca chain flagged ERG/TMPRSS2 : {comparison_df['baca_has_ERG_chain'].sum()}")
    print(f"  ...k_local == baca_k            : {comparison_df['k_match'].sum()}")
    print(f"  ...chromosome sets identical    : {comparison_df['chromosome_sets_match'].sum()}")


if __name__ == '__main__':
    main()
