import os
import shutil
import gzip

# ── CONFIG ──────────────────────────────────────────────────────────────────
# Change this to wherever your downloaded UUID folders are
RAW_DIR      = "/data/raw"
CLINICAL_DIR = "/data/clinical"
MAF_DIR      = "/data/maf"
CNV_DIR      = "/data/cnv"
# ────────────────────────────────────────────────────────────────────────────

def setup_folders():
    for folder in [CLINICAL_DIR, MAF_DIR, CNV_DIR]:
        os.makedirs(folder, exist_ok=True)
    print("✓ Output folders ready")

def organize_files():
    clinical_count = 0
    maf_count      = 0
    cnv_count      = 0
    skipped        = 0

    # Walk through every UUID subfolder
    for uuid_folder in os.listdir(RAW_DIR):
        uuid_path = os.path.join(RAW_DIR, uuid_folder)

        if not os.path.isdir(uuid_path):
            continue

        # Look at every file inside the UUID folder
        for filename in os.listdir(uuid_path):
            src = os.path.join(uuid_path, filename)

            # ── Clinical files ──────────────────────────────────────────────
            if filename.endswith('.txt') and 'nationwidechildrens' in filename:
                dst = os.path.join(CLINICAL_DIR, filename)
                shutil.copy2(src, dst)
                clinical_count += 1

            # ── MAF files ───────────────────────────────────────────────────
            elif filename.endswith('.maf.gz'):
                dst = os.path.join(MAF_DIR, filename)
                shutil.copy2(src, dst)
                maf_count += 1

            # ── CNV files (ASCAT3 only) ──────────────────────────────────────
            elif 'ascat3.gene_level_copy_number' in filename and filename.endswith('.tsv'):
                dst = os.path.join(CNV_DIR, filename)
                shutil.copy2(src, dst)
                cnv_count += 1

            else:
                skipped += 1

    print(f"✓ Clinical files copied : {clinical_count}")
    print(f"✓ MAF files copied      : {maf_count}")
    print(f"✓ CNV files copied      : {cnv_count}")
    print(f"  Skipped (logs etc)    : {skipped}")

def verify():
    print("\n── Verification ────────────────────────────────")
    for label, folder in [("Clinical", CLINICAL_DIR),
                           ("MAF",      MAF_DIR),
                           ("CNV",      CNV_DIR)]:
        files = os.listdir(folder)
        print(f"  {label:10s}: {len(files)} files")

    print("\n── Clinical files ──────────────────────────────")
    for f in sorted(os.listdir(CLINICAL_DIR)):
        print(f"  {f}")

    print("\n── First 3 MAF files ───────────────────────────")
    for f in sorted(os.listdir(MAF_DIR))[:3]:
        print(f"  {f}")

    print("\n── First 3 CNV files ───────────────────────────")
    for f in sorted(os.listdir(CNV_DIR))[:3]:
        print(f"  {f}")

if __name__ == "__main__":
    print("Starting data organization...\n")
    setup_folders()
    organize_files()
    verify()
    print("\n✓ Done — data is organized and ready for analysis")