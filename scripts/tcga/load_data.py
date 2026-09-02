import os
from operator import index
from pickle import FALSE
from traceback import print_tb

import pandas as pd
import gzip
import shutil

import requests
import json
import re

from pandas import read_csv
from tabulate import tabulate

# ── CONFIG ──────────────────────────────────────────────────────────────────
CLINICAL_DIR = "/data/clinical"
MAF_DIR = "/data/maf"
CNV_DIR      = "/data/cnv"
SOMATIC_MUT_DIR = "/data/somatic_mutations"
SOMATIC_MUT_RAW_DIR = "/data/somatic_mutations/raw_maf"
SOMATIC_MUT_CSV_DIR= "/data/somatic_mutations/csv"
MANIFEST_FILE_PATH = "/Users/anantkumarsingh/research_data/tgca_prad_data/gdc_manifest.2026-04-12.155357.txt"
# ────────────────────────────────────────────────────────────────────────────

pd.set_option('display.max_columns', None)   # show all columns
pd.set_option('display.max_rows', None)      # show all rows
pd.set_option('display.width', None)         # no line wrapping
pd.set_option('display.max_colwidth', None)  # full content in each cell


# Create Clean CSV file for each clinical txt file
def load_core_patient_data_file():
    core_df = pd.read_csv(CLINICAL_DIR+"/nationwidechildrens.org_clinical_patient_prad.txt", delimiter = "\t", skiprows = [0,2])

    not_in_use_cols = ["bcr_patient_uuid","form_completion_date", "icd_o_3_site","icd_10",
                       "days_to_initial_pathologic_diagnosis", "tumor_tissue_site",
                       "clinical_N", "clinical_stage", "disease_code", "extranodal_involvement",
                       "pathologic_M", "pathologic_stage", "project_code", "stage_other", "system_version"
                       ]

    core_df.replace({'[Not Available]': None,
                     '[Not Applicable]': None,
                     '[Unknown]': None,
                     '[Discrepancy]': None}, inplace=True)

    core_df = core_df.drop(columns=not_in_use_cols)

    print(f"Dropped {len(not_in_use_cols)} columns: {not_in_use_cols}")
    print(f"Remaining shape: {core_df.shape}")

    core_df.to_csv(CLINICAL_DIR+'/clinical_patient_prad.csv', index=False)

    return core_df

def load_follow_up_data():
    df = pd.read_csv(CLINICAL_DIR + "/nationwidechildrens.org_clinical_follow_up_v1.0_prad.txt", delimiter="\t",
                     skiprows=[0, 2])

    not_in_use_cols = ["bcr_patient_uuid","form_completion_date"]


    df.replace({'[Not Available]': None,
                '[Not Applicable]': None,
                '[Unknown]': None,
                '[Discrepancy]': None}, inplace=True)
    df = df.drop(columns=not_in_use_cols)

    df = df.rename(columns={
        'vital_status': 'vital_status_updated'
    })

    print(f"Dropped {len(not_in_use_cols)} columns: {not_in_use_cols}")
    print(f"Remaining shape: {df.shape}")

    # This file contains duplicate patients
    # We keep the last record the patient iff vital status exists
    # If vital status does not exist, fall back to previous row
    print(f"Follow-up table  : {len(df)} rows (before dedup)")
    print(f"Unique patients  : {df['bcr_patient_barcode'].nunique()}")
    print(f"Patients with 2+ rows: "
          f"{(df.groupby('bcr_patient_barcode').size() > 1).sum()}")

    # Convert follow-up days to numeric so we can sort properly
    df['days_to_last_followup'] = pd.to_numeric(
        df['days_to_last_followup'], errors='coerce')

    sorted_df = df.sort_values(
        ['bcr_patient_barcode', 'days_to_last_followup'], # Order of Rows
        ascending=[True, True],
        na_position='first'     # Pushing Missing days so that entry is of no use

    )
    # Most recent entry is the row with the HIGHEST NUMBER OF DAYS
    latest_df = sorted_df.drop_duplicates(subset='bcr_patient_barcode', keep='last')

    print(f"Follow-up table  : {len(latest_df)} rows (after dedup)")

    latest_df.to_csv(CLINICAL_DIR + '/clinical_follow_up_prad.csv', index=False)

    return latest_df

def load_new_tumor_event_patients():
    nte_df = pd.read_csv(CLINICAL_DIR+"/nationwidechildrens.org_clinical_nte_prad.txt", delimiter = "\t", skiprows = [0,2])

    #threshold = 60
    not_in_use_cols = ["bcr_patient_uuid"]

    # cols_to_drop = list(set(
    #     [
    #     col for col in nte_df.columns
    #     if (nte_df[col] == '[Not Available]').sum() >= threshold or
    #         (nte_df[col] == '[Not Applicable]').sum() >= threshold
    # ] + not_in_use_cols
    # ))

    nte_df.replace({'[Not Available]': None,
                     '[Not Applicable]': None,
                     '[Unknown]': None,
                     '[Discrepancy]': None}, inplace=True)

    nte_df = nte_df.drop(columns=not_in_use_cols, errors='ignore')

    print(f"Dropped {len(not_in_use_cols)} columns: {not_in_use_cols}")
    print(f"Remaining shape: {nte_df.shape}")

    nte_df.to_csv(CLINICAL_DIR+'/clinical_new_tumor_event_prad.csv', index=False)

    return nte_df

# merge relevant columns from all the clinical files
def load_clinical_table(core_clinical_df, followup_df, nte_df):
    # Columns
    follow_up_cols = ['bcr_patient_barcode',
                      'vital_status_updated',
                      'person_neoplasm_cancer_status',
                      'days_to_last_followup',
                      'days_to_death',
                      'days_to_first_biochemical_recurrence',
                      'radiation_therapy',
                      'postoperative_rx_tx',
                      'new_tumor_event_after_initial_treatment']

    follow_up_cols = [c for c in follow_up_cols if c in followup_df.columns]

    # Merge core_patient + follow-up (left join — keep all patients)
    merged = core_clinical_df.merge(
        followup_df[follow_up_cols],
        on='bcr_patient_barcode',
        how='left',
        suffixes=('_patient', '_followup')
    )

    nte_cols = ['bcr_patient_barcode',
                'new_neoplasm_event_type',
                'days_to_new_tumor_event_after_initial_treatment',
                'tumor_progression_post_ht',
                'additional_radiation_therapy',
                'additional_pharmaceutical_therapy']

    nte_cols = [c for c in nte_cols if c in nte_df.columns]

    # Merge in NTE events (left join — patients without events get NaN)
    merged = merged.merge(
        nte_df[nte_cols],
        on='bcr_patient_barcode',
        how='left'
    )

    print(f"\nFinal merged table: {merged.shape[0]} patients × "
          f"{merged.shape[1]} columns")

    return merged


def cleanup_merged_data(df):
    print(f"Before cleanup: {df.shape}")

    # ── 1. DROP COMPLETELY EMPTY COLUMNS ────────────────────────────────────────
    empty_cols = [
        'clinical_N', 'clinical_stage', 'disease_code', 'extranodal_involvement',
        'pathologic_M', 'pathologic_stage', 'project_code', 'stage_other',
        'system_version'
    ]
    df = df.drop(columns=empty_cols, errors='ignore')

    # ── 2. DROP SINGLE-VALUE COLUMNS ────────────────────────────────────────────
    single_value_cols = ['gender', 'informed_consent_verified']
    df = df.drop(columns=single_value_cols, errors='ignore')

    # ── 3. DROP NEAR-EMPTY UNINFORMATIVE COLUMNS ────────────────────────────────
    df = df.drop(columns=['histological_type_other'], errors='ignore')

    # ── 4. RESOLVE DUPLICATE PAIRS — keep best version, drop the other ──────────
    # days_to_last_followup → keep _patient (492 vs 430)
    df = df.rename(columns={'days_to_last_followup_patient': 'days_to_last_followup'})
    df = df.drop(columns=['days_to_last_followup_followup'], errors='ignore')

    # days_to_death → keep _patient, drop _followup (identical)
    df = df.rename(columns={'days_to_death_patient': 'days_to_death'})
    df = df.drop(columns=['days_to_death_followup'], errors='ignore')

    # person_neoplasm_cancer_status → keep _patient (409 vs 344)
    df = df.rename(columns={'person_neoplasm_cancer_status_patient': 'cancer_status'})
    df = df.drop(columns=['person_neoplasm_cancer_status_followup'], errors='ignore')

    # radiation_therapy → keep _followup (405 vs 291)
    df = df.rename(columns={'radiation_therapy_followup': 'radiation_therapy'})
    df = df.drop(columns=['radiation_therapy_patient'], errors='ignore')

    # postoperative_rx_tx → keep _followup (401 vs 292)
    df = df.rename(columns={'postoperative_rx_tx_followup': 'postop_treatment'})
    df = df.drop(columns=['postoperative_rx_tx_patient'], errors='ignore')

    # new_tumor_event → keep _followup (400 vs 324)
    df = df.rename(columns={
        'new_tumor_event_after_initial_treatment_followup': 'had_new_tumor_event'
    })
    df = df.drop(columns=['new_tumor_event_after_initial_treatment_patient'],
                 errors='ignore')

    # days_to_first_biochemical_recurrence → keep _patient (70 vs 59)
    df = df.rename(columns={
        'days_to_first_biochemical_recurrence_patient': 'days_to_psa_recurrence'
    })
    df = df.drop(columns=['days_to_first_biochemical_recurrence_followup'],
                 errors='ignore')

    # ── 5. RENAME FOR CLARITY ────────────────────────────────────────────────────
    df = df.rename(columns={
        'primary_pattern':           'gleason_primary',
        'secondary_pattern':         'gleason_secondary',
        'tertiary_pattern':          'gleason_tertiary',
        'new_neoplasm_event_type':   'recurrence_type',
        'days_to_new_tumor_event_after_initial_treatment': 'days_to_recurrence',
        'vital_status_updated':      'vital_status_followup'
    })

    # ── 6. FINAL SHAPE AND SUMMARY ───────────────────────────────────────────────
    print(f"After cleanup : {df.shape}")
    print(f"\nFinal columns ({len(df.columns)}):")
    for col in df.columns:
        filled = df[col].notna().sum()
        print(f"  {col:<55} {filled}/500 filled")

    return df

def load_maf_files_to_csv():
    if len(os.listdir(SOMATIC_MUT_RAW_DIR)) != 501:
        unzip_file_count = 0
        csv_file_count = 0
        with os.scandir(MAF_DIR) as somatic_mut_folder:
            for mut_file in somatic_mut_folder:
                if mut_file.is_file() and mut_file.name.endswith(".gz"):
                    # print(mut_files.name)
                    input_path = mut_file.path
                    output_filename = mut_file.name[:-3]  # remove ".gz"
                    output_path = os.path.join(SOMATIC_MUT_RAW_DIR, output_filename)

                    # unzip
                    with gzip.open(input_path, 'rb') as f_in:
                        with open(output_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)

                    print(f"Unzipped: {mut_file.name} -> {output_filename}")
                    unzip_file_count += 1
            print("Total Files Unzipped:", unzip_file_count)
            with os.scandir(SOMATIC_MUT_RAW_DIR) as somatic_mut_unzipped_folder:
                for mut_raw_file in somatic_mut_unzipped_folder:
                    if mut_raw_file.is_file() and mut_raw_file.name.endswith(".maf"):
                        df = read_csv(mut_raw_file.path, delimiter = "\t", skiprows = [0,1,2,3,4,5,6])

                        # [CODE TO FILTER IMPORTANT COLUMNS]
                        MAF_COLS_TO_KEEP = [
                            # Patient identity
                            'Tumor_Sample_Barcode',  # join key to clinical data
                            'Tumor_Sample_UUID',  # tumor aliquot UUID
                            'case_id',  # case-level UUID

                            # Gene identity
                            'Hugo_Symbol',  # gene name - SPOP, PTEN, TP53 etc
                            'Chromosome',  # chr1 through chrX

                            # Genomic location
                            'Start_Position',
                            'End_Position',

                            # Mutation type
                            'Variant_Classification',  # Missense, Nonsense, Frame_Shift etc
                            'Variant_Type',  # SNP, INS, DEL
                            'IMPACT',  # HIGH, MODERATE, LOW

                            # The actual change
                            'Reference_Allele',
                            'Tumor_Seq_Allele2',
                            'HGVSp_Short',  # protein change e.g. p.V822M
                            'HGVSc',  # DNA change e.g. c.2464G>A
                            'Codons',

                            # Transcript information
                            'Exon_Number',  # which exon
                            'CANONICAL',  # is this the canonical transcript
                            'all_effects',  # all transcript consequences

                            # Sequencing quality - needed for VAF and clonality
                            't_depth',  # total reads
                            't_alt_count',  # mutant reads
                            't_ref_count',  # reference reads
                            'n_depth',  # normal tissue depth — quality check

                            # Functional impact
                            'SIFT',
                            'PolyPhen',
                            'hotspot',  # known cancer hotspot YES/NO
                            'callers',  # mutation calling algorithms

                            # Known cancer variants
                            'COSMIC',

                            # Population frequency - to confirm truly somatic
                            'gnomAD_AF',
                            'MAX_AF',

                            # Quality confirmation
                            'Mutation_Status',  # should always be Somatic
                        ]

                        filtered_df = df[MAF_COLS_TO_KEEP]
                        # compute Variant Allele Frequency for clonality analysis
                        filtered_df.insert(16, 'VAF', filtered_df['t_alt_count'] / filtered_df['t_depth'])

                        filtered_df.to_csv(SOMATIC_MUT_CSV_DIR + "/" + mut_raw_file.name[:-4] + '.csv', index= False)
                        csv_file_count += 1
                        print(mut_raw_file.name, " --> ", mut_raw_file.name[:-4], ".csv")
                print("Total Files converted to CSV:", csv_file_count)

    else: # files already unzipped and exist in SOM_MUT FOlder
        print("All Files exist in", SOMATIC_MUT_RAW_DIR)

    return MAF_COLS_TO_KEEP


def load_somatic_mutations_table(MAF_COLS_TO_KEEP):
    all_patients = []
    file_count = 0
    error_count = 0

    FINAL_COLS = MAF_COLS_TO_KEEP + ['Total mutations', 'VAF']

    with os.scandir(SOMATIC_MUT_CSV_DIR) as somatic_mut_csv_folder:
        for csv_mut_file in somatic_mut_csv_folder:
            if csv_mut_file.is_file() and csv_mut_file.name.endswith(".csv"):
                # combine entries of each file using partial 'Tumor_Sample_Barcode' as JOIN KEY
                try:
                    # Step 1 — read the file
                    df = pd.read_csv(csv_mut_file.path, low_memory=False)

                    # Step 2 — keep only columns that exist in this file
                    # (some files may be missing optional columns)
                    cols = [c for c in FINAL_COLS if c in df.columns]
                    df = df[cols]

                    # Step 3 — add patient_id from barcode first 12 characters
                    # TCGA-VN-A88N-01A-11D... -> TCGA-VN-A88N
                    df['patient_id'] = df['Tumor_Sample_Barcode'].str[:12]

                    # Step 4 — add to collection
                    all_patients.append(df)
                    file_count += 1

                    # Progress update every 100 files
                    if file_count % 100 == 0:
                        print(f"  Processed {file_count} files...")

                except Exception as e:
                    print(f"  Warning — could not read {csv_mut_file.name}: {e}")
                    error_count += 1

    # Step 5 — combine all patients into one table
    print(f"\n  Files read successfully : {file_count}")
    print(f" Files with errors : {error_count}")

    mut_table = pd.concat(all_patients, ignore_index=True)

    print(f"  Avg mutations/patient   : {round(len(mut_table) / mut_table['patient_id'].nunique(), 1)}")

    mut_table.to_csv(SOMATIC_MUT_DIR + '/somatic_mutation_table.csv', index = False)


def build_cnv_patient_mapping(manifest_path):
    """
    One-time GDC API query.
    Returns a dictionary: {uuid_in_filename -> patient_id}
    """
    cache_path = CNV_DIR + '/cnv_patient_mapping.csv'
    if os.path.exists(cache_path):
        print(f"Loading cached mapping from {cache_path}")
        mapping = pd.read_csv(cache_path)
        return dict(zip(mapping['filename_uuid'], mapping['patient_id']))

    manifest = pd.read_csv(manifest_path, sep='\t')
    # Bug fix 1: regex=False so the dots aren't treated as wildcards; na=False skips NaN filenames
    cnv_files = manifest[
        manifest['filename'].str.contains('ascat3.gene_level_copy_number', regex=False, na=False)
    ].copy()

    # Extract the UUID embedded in each filename
    # TCGA-PRAD.94b2fe48-....ascat3... → 94b2fe48-...
    cnv_files['filename_uuid'] = cnv_files['filename'].str.extract(
        r'TCGA-PRAD\.([a-f0-9-]+)\.ascat3'
    )

    # Use manifest file IDs to query GDC
    file_ids  = cnv_files['id'].tolist()
    batch     = 200
    api_rows  = []

    print(f"Querying GDC API for {len(file_ids)} files...")

    for start in range(0, len(file_ids), batch):
        chunk = file_ids[start : start + batch]

        payload = {
            "filters": {
                "op": "in",
                "content": {"field": "file_id", "value": chunk}
            },
            "fields": "file_id,file_name,cases.submitter_id",
            "size": batch,
            "format": "JSON"
        }

        # POST avoids 414 URI-too-long when sending hundreds of UUIDs
        response = requests.post(
            "https://api.gdc.cancer.gov/files",
            json=payload,
            timeout=30
        )
        response.raise_for_status()
        hits = response.json()['data']['hits']

        for hit in hits:
            # Bug fix 3: guard against files with no associated case
            cases = hit.get('cases', [])
            if not cases:
                print(f"  Warning: no case for file_id {hit['file_id']} — skipping")
                continue
            api_rows.append({
                'manifest_file_id': hit['file_id'],
                'filename':         hit['file_name'],
                'patient_id':       cases[0]['submitter_id']
            })

        print(f"  {min(start+batch, len(file_ids))}/{len(file_ids)} done")

    # Bug fix 4: handle the case where nothing came back from the API
    if not api_rows:
        raise RuntimeError("GDC API returned no results — check file IDs and network access")

    api_df = pd.DataFrame(api_rows)

    # Now merge to get filename_uuid alongside patient_id
    mapping = cnv_files[['id', 'filename', 'filename_uuid']].merge(
        api_df[['manifest_file_id', 'patient_id']],
        left_on='id',
        right_on='manifest_file_id',
        how='left'
    )

    # filename_uuid → patient_id
    uuid_to_patient = dict(
        zip(mapping['filename_uuid'], mapping['patient_id'])
    )

    print(f"\nDictionary built: {len(uuid_to_patient)} entries")
    print(f"Sample entries:")
    for k, v in list(uuid_to_patient.items())[:3]:
        print(f"  {k} → {v}")

    # Save so you never query again
    mapping.to_csv(CNV_DIR + '/cnv_patient_mapping.csv', index=False)

    return uuid_to_patient


def load_cnv_table(cnv_dir, uuid_to_patient):
    # Only keep what's needed for analysis — dropping start/end/gene_id cuts
    # 29M-row memory footprint significantly

    cnv_table_cache = CNV_DIR + '/copy_number_changes_table_prad.csv'
    if os.path.exists(cnv_table_cache):
        print("CNV Table already exists, returning CNV Dataframe...")
        df = pd.read_csv(cnv_table_cache)
        return df

    COLS_TO_KEEP = ['gene_name', 'chromosome', 'copy_number', 'min_copy_number', 'max_copy_number']

    all_dfs = []

    for fname in os.listdir(cnv_dir):
        if not fname.endswith('.tsv'):
            continue

        match = re.search(r'TCGA-PRAD\.([a-f0-9-]+)\.ascat3', fname)
        if not match:
            continue

        filename_uuid = match.group(1)
        patient_id = uuid_to_patient.get(filename_uuid)

        if patient_id is None:
            print(f"  No mapping found for {fname}")
            continue

        df = pd.read_csv(
            os.path.join(cnv_dir, fname),
            sep='\t',
            usecols=COLS_TO_KEEP  # read only what you need
        )
        df['patient_id'] = patient_id

        all_dfs.append(df)

    # Bug fix 1: guard against empty list
    if not all_dfs:
        raise RuntimeError("No CNV files matched — check cnv_dir and uuid_to_patient")

    cnv = pd.concat(all_dfs, ignore_index=True)

    # Bug fix 2: use the parameter, not the global
    out_path = os.path.join(cnv_dir, 'copy_number_changes_table_prad.csv')
    cnv.to_csv(out_path, index=False)
    print(f"Saved {len(cnv):,} rows → {out_path}")

    return cnv

# ---------- MAIN ----------

# core_clinical_df = load_core_patient_data_file()
# followup_df = load_follow_up_data()
# nte_df = load_new_tumor_event_patients()
#
# merged_clinical_table = load_clinical_table(core_clinical_df, followup_df, nte_df)
# final_clinical_features_table = cleanup_merged_data(merged_clinical_table)
#
# final_clinical_features_table.to_csv(CLINICAL_DIR + '/patients_clinical_features.csv', index=False)
#
# print("✓ Saved to data/clinical/clinical_merged.csv")

# MAF_COLS_TO_KEEP = load_maf_files_to_csv()
# load_somatic_mutations_table(MAF_COLS_TO_KEEP)

uuid_to_patient_dict = build_cnv_patient_mapping(MANIFEST_FILE_PATH)

cnv_df = load_cnv_table(CNV_DIR, uuid_to_patient_dict)

one_patient = cnv_df[cnv_df['patient_id'] == cnv_df['patient_id'].iloc[0]]

total       = len(one_patient)
filled      = one_patient['copy_number'].notna().sum()
empty       = one_patient['copy_number'].isna().sum()

print(f"One patient — total genes : {total:,}")
print(f"  Copy number filled      : {filled:,}  ({filled/total*100:.1f}%)")
print(f"  Copy number empty       : {empty:,}  ({empty/total*100:.1f}%)")



