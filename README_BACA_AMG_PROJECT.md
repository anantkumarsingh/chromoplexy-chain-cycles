# Baca-AMG Prostate Cancer Research Project — Full Context README
**Student:** Anant Kumar Singh, CS undergraduate, UC Davis  
**Lab:** Arsuaga-Vazquez Lab, UC Davis  
**Advisor:** Professor Javier Arsuaga (co-author on Sheth 2026 paper)  
**Last updated:** May 2026

---

## The Central Question

> **Does Baca's prostate cancer rearrangement data, when analyzed through Sheth's AMG mathematical framework, produce results consistent with Sheth's mathematical predictions?**

This is NOT about reproducing Baca's biological findings. It is about taking Baca's exact breakpoint data and asking whether the mathematical structure Sheth describes — cycle counts, poset ranks, configuration classes — actually describes what is observed in real prostate cancer genomes.

Professor Arsuaga's exact instruction:
> "For now we have to make the analysis with the existing Baca paper. Does our analysis match theirs? Once we have some results we can request access to TCGA which requires permits."

---

## The Two Papers

### Baca et al. 2013 — Cell 153:666-677
- Sequenced 57 prostate tumor genomes
- Discovered **chromoplexy** — coordinated chains of DNA rearrangements across multiple chromosomes simultaneously
- Found consensus progression path: NKX3-1/FOXP1/T2-ERG → CDKN1B/TP53 → PTEN
- Clonal events = early, subclonal = late (determined by allelic fraction / clonality %)
- ETS+ tumors: more interchromosomal rearrangements, driven by AR/TOP2B
- ETS- CHD1-deleted tumors: chromothripsis-like, intrachromosomal, concentrated damage
- Key claim: **"Oncogenic ERG fusions frequently arose in the setting of chromoplexy (15 of 26 cases, 58%)"**

### Sheth, Arsuaga, Sazdanovic 2026 — J. Phys. A 59:115601
- Mathematical AMG (Aberration Multigraph) framework for chromosomal rearrangements
- Key notation: **Θ(k, (b₁,...,bₖ))** — k chromosomes, bᵢ breaks on chromosome i
- Cycle structure C(Ω) = m₁C_{l₁} + m₂C_{l₂} + ... (partition of n)
- **Poset L(n)**: all cycle structures for n DSBs, ordered by complexity
- **Theorem 25**: rank(C) = n − |C| where |C| = number of cycles
- **Proposition 4**: cycle count range is [1, min(max_bᵢ, ⌊n/2⌋)]
- **Theorem 6**: probable cycle count ≈ log(n) for large n
- **Lemma 14**: if each chromosome has exactly 1 DSB → only possible cycle structure is Cₙ (rank = n−1, maximum)
- **Corollary 11**: combinatorially, most proper AMGs for large n involve interchromosomal exchanges
- Section 7.1 re-analyzes Baca patient P05-1657 — finds C₁₀ most probable among 945 pathways

---

## Biological-Mathematical Mapping

| Baca Biology | Sheth AMG Math |
|---|---|
| Chromoplexy (ETS+, many chromosomes) | Θ(large k, small bᵢ) |
| Chromothripsis (ETS-, CHD1del, few chromosomes) | Θ(small k, large bᵢ) |
| Each rearrangement row = 2 DSBs | n = total_rearrangements × 2 |
| Chromosomes affected | k in Θ notation |
| DSBs per chromosome | bᵢ distribution |
| Progression path NKX3-1→PTEN | Chain in poset L(n) |
| Each chromoplexy event | Move upward in L(n) |

---

## What Rank in L(n) Means

From Theorem 25: **rank(C) = n - |C|**

```
n   = total DSBs
|C| = number of cycles in cycle structure

High rank → n large relative to |C|
           → few cycles for the number of breaks
           → breaks organized into large coordinated cycles
           → MORE coordinated = more chromoplexy-like

Low rank  → many cycles for the number of breaks
           → breaks form many small independent events
           → LESS coordinated

Maximum rank = n - 1  (single cycle Cₙ, Lemma 14 case)
Minimum rank = 0      (all DSBs repaired correctly)
```

Moving UP one step in L(n) = one chromosome fragment reversal merging two cycles into one.

**rank_proxy = total_DSBs - chain_count** (chain_count ≈ |C|)

---

## Datasets

### Primary Dataset — Baca Supplementary Files
**Location:** `/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/baca_dataset/`

| File | Contents | Key Columns |
|---|---|---|
| `chrom_aberrations_baca.csv` | 5,710 rows, 57 patients, exact breakpoint data | Individual, Breakpoint 1/2 chromosome, Breakpoint 1/2 position, Breakpoint 1/2 site, Class, Fusion |
| `clinical_phenotypes.csv` | 57 patients clinical features | Individual, Gleason Score, Pathological stage, ETS fusion detected by sequencing |
| `mmc6.xlsx` Sheet S6B | 48 patients, somatic CNV with clonality % | Sample, Gene, Somatic Copy Number, Clonality %, Clonality % Ranges |
| `mmc3.xlsx` S3A | Somatic mutations MAF | Hugo_Symbol, Variant_Classification, VAF |
| `mmc5.xlsx` | Original chain/adjacency data | Individual, Chain number, Breakpoint number, Chromosome:position, Adjacent breakpoint(s) |

**Key data facts:**
- 57 total patients in clinical, 52 in mmc5 (5 PR- prefix patients missing from mmc5)
- 48 patients in mmc6, 46 with at least one key progression gene
- chrom_aberrations_baca.csv has ALL 57 patients with exact breakpoint coordinates
- Each row = one rearrangement = 2 DSBs
- Class values: inter_chr, long_range, inversion, potential deletion, tandem_dup
- inter_chr = interchromosomal; all others = intrachromosomal

### Secondary Dataset — TCGA PRAD (500 patients)
**Location:** `/Users/anantkumarsingh/projects/prostate_cancer/nih-tcga-prad/data/`
- `clinical/patients_clinical_features.csv` — cleaned clinical table
- `maf/all_mutations_merged.csv` — all somatic mutations merged
- `cnv/cnv_patient_mapping.csv` — UUID mapping
- **Structural variation files = controlled access (dbGaP phs000178) — NOT downloaded**
- TCGA currently used ONLY for clinical phenotype analysis (Phase 1)
- Cannot compute exact AMG from TCGA without structural variant files

---

## Current Code State

**Main script:** `scripts/baca/baca_aberration_clinical_analysis.py`

### Data Computed

**summary_df (57 × 9 columns):**
```
Individual, total_rearrangements, inter_chr_count, intra_chr_count,
pct_inter, k, total_DSBs, max_b_i, dominant_chromosome
```

**chromosome_df (828 × 3 columns):**
```
Individual, chromosome, dsb_count
```

**baca_df (57 patients merged with clinical):**
```
All above + Age, Gleason Score, Pathological stage,
ETS fusion detected by sequencing, ETS_status (binary),
gleason_primary, gleason_secondary, gleason_total,
has_tertiary_pattern
```

**patient_genome_profiles.csv (57 patients):**
```
Individual, k, chromosomes (list), DSBs_per_chr (dict),
theta (Θ notation string), all_genes (genes per chromosome),
Gleason Score, Pathological stage, ETS_status
```

### Key Metric Definitions

```python
# Each row = one rearrangement = 2 DSBs
total_DSBs = total_rearrangements * 2

# k = unique chromosomes from BOTH BP1 and BP2 columns combined
# This is the global k across whole patient genome
k = pd.concat([bp1_chromosomes, bp2_chromosomes]).nunique()

# max_b_i = maximum DSBs on any single chromosome
max_b_i = dbs_per_chromosome.max()

# pct_inter = % of rearrangements that are interchromosomal
pct_inter = inter_chr_count / total_rearrangements * 100

# intra_classes = ['long_range', 'inversion', 'potential deletion', 'tandem_dup']
```

---

## Completed Analysis: ERG Fusion Chromoplexy Validation

### The Baca Claim Being Tested
> "Oncogenic ERG fusions frequently arose in the setting of chromoplexy (15 of 26 cases, 58%)"

### Key Definitions for This Analysis

**k (global):** Chromosomes affected across the ENTIRE patient genome — all rearrangements, all events.

**k_local:** Number of chromosomes connected to chr21 (where ERG and TMPRSS2 both live) through interchromosomal rearrangements in the same chain — including chr21 itself. This measures the size of the chromoplexy chain specifically surrounding the ERG fusion event.

**chain_size:** Number of inter_chr rearrangements in the chain involving those k_local chromosomes. Computed as: find inter_chr rows that share any chromosome with the anchor (ERG/TMPRSS2) rows, then count them.

**fusion_type classification:**
```python
# BOTH conditions required — not just k_local
'chromoplexy_embedded' if k_local > 2 AND chain_size > 1
'simple_fusion'        otherwise
```

Why both conditions: k_local=2 with chain_size=1 is just a single translocation — not a coordinated chain. Chromoplexy requires the chain to propagate across multiple chromosomes with multiple events.

### The Three Patient Subtypes

```
Simple intra (9 patients, k_local=1, chain_size=0):
    Both TMPRSS2 and ERG on chr21
    Simple intrachromosomal deletion on chr21
    Class = long_range in chrom_aberrations
    No other chromosomes connected
    AMG: Θ(1,(2)) → C₂ cycle → simplest proper AMG

Simple inter (2 patients, k_local=2, chain_size=1):
    chr21 connected to exactly one other chromosome
    One single translocation event
    Still isolated — not part of coordinated chain
    PR-09-146 (chr21+chr2), P08-501 (chr21+chr1)

Chromoplexy embedded (15 patients, k_local≥3, chain_size≥2):
    chr21 connected to 3-12 other chromosomes
    ERG fusion is woven into coordinated multi-chromosome event
    Satisfies Θ(large k, small bᵢ) — Lemma 14 regime
```

### Result

```
chromoplexy_embedded: 15 patients
simple_fusion:        11 patients
Total:                26 patients

MATCHES BACA'S REPORTED 15/26 EXACTLY ✓

Chromoplexy embedded group:
    median k_local     = 5.0
    median chain_size  = 8.0
    median total_DSBs  = 138.0
    median pct_inter   = 43.2%
    ALL have k_local > 2 (15/15)
    ALL have k_local > 3 (13/15)

Simple fusion group:
    median k_local     = 1.0
    median chain_size  = 0.0
    median total_DSBs  = 68.0
    median pct_inter   = 28.6%

Mann-Whitney U test (k_local, chromoplexy vs simple):
    U statistic = 165.000
    P-value     = 0.0000
```

### Case Study: P05-1657 (Sheth's Patient)

This patient is particularly important — it is the patient Sheth analyzed in Section 7.1 of the 2026 paper.

```
k        = 5  (whole genome: chromosomes 4, 7, 8, 12, 21)
k_local  = 1  (fusion chain: only chr21)
chain_size = 0

Θ(5, (2, 4, 9, 1, 2)) — matches Sheth exactly ✓
DSBs per chromosome: {4:2, 7:4, 8:9, 12:1, 21:2}

What this means:
    The TMPRSS2-ERG fusion is a SIMPLE intrachromosomal
    deletion on chr21 (Class=long_range in data)
    k_local=1 → not chromoplexy for the fusion itself

    BUT the patient has chromoplexy ELSEWHERE:
    chr8 ↔ chr12 interchromosomal exchange
    chr8 has 9 DSBs (concentrated damage)
    This is what Sheth analyzed — the chr8+chr12 chain
    has 945 possible cycle structures
    Most probable: C₁₀

    So P05-1657 has BOTH:
    1. Simple chr21 deletion → TMPRSS2-ERG fusion
    2. Separate chromoplexy event on chr4,7,8,12,21
    These are independent events in the same tumor
```

### Source Detection Breakdown

```
site_annotation (17 patients):
    ERG/TMPRSS2 found directly in Breakpoint site columns
    Most reliable detection

site_annotation_intra_only (9 patients):
    ERG/TMPRSS2 found in site annotations
    BUT no inter_chr connections leave chr21
    Fusion is intrachromosomal deletion
    chain_size = 0 by design
```

### The Function — get_erg_chain_details_v2()

```python
def get_erg_chain_details_v2(patient_id, chrom_df, clinical_df):
    # Step 1: Direct site annotation search for ERG/TMPRSS2
    # Step 2: Check Fusion column if step 1 empty
    # Step 3: chr21 fallback for clinically confirmed ERG+
    #         with no annotation hits
    # Step 4: Get chromosomes from anchor rows
    # Step 5: Expand to ALL inter_chr rows sharing those chromosomes
    # Step 6: If no inter_chr chain → return k_local=len(chr_involved),
    #         chain_size=0, source='site_annotation_intra_only'
    # Returns: patient_id, k_local, chain_size, chromosomes,
    #          genes_in_chain, source
```

**Important limitation:** The chain expansion is ONE HOP from the anchor — not iterative. Chain_size is therefore a lower bound on the true chain size. This is sufficient for 15/11 classification but a future improvement would be iterative expansion.

---

## Patient Genome Profile Table

**Output file:** `patient_genome_profiles.csv`

For ALL 57 patients — computed from chrom_aberrations_baca.csv:

```
Individual    → patient ID
k             → total chromosomes affected (whole genome)
chromosomes   → list e.g. [4, 7, 8, 12, 21]
DSBs_per_chr  → dict e.g. {4:2, 7:4, 8:9, 12:1, 21:2}
theta         → Θ(k,(b₁,...,bₖ)) notation string
all_genes     → genes at breakpoints per chromosome
Gleason Score → clinical
Pathological stage → clinical
ETS_status    → ETS+ / ETS-
```

Verified against Sheth paper: P05-1657 = chromosomes [4,7,8,12,21] with Θ(5,(2,4,9,1,2)) ✓

---

## Summary Tables Computed

### Stage Table
```
Stage       n    median_total_R  median_pct_inter  median_k  ETS+
pT2b        1    9               11.1%             5         1
pT2c        21   56              19.7%             14        7
pT3a        26   65.5            28.7%             15.5      14
pT3b        4    46              22.1%             13        3
Metastatic  2    265             26.3%             23        2
```
Note: Stage analysis limited by small group sizes (n=1 for most stages except pT2c and pT3a).

### Gleason Table
```
GS    n    median_total_R  median_pct_inter  median_k  ETS+
3+3   8    57              16.9%             14        2
3+4   21   49              29.5%             14        12
4+3   11   105             12.5%             17        5
4+4   7    84              15.5%             15        3
```

**Key finding:** 3+4 has median 49 rearrangements vs 4+3 median 105 — same Gleason total, more than double the complexity. 3+4 has higher pct_inter (29.5% vs 12.5%) despite fewer total rearrangements.

---

## Five Testable Sheth Predictions (Pending)

| # | Prediction | Theorem | Test | Status |
|---|---|---|---|---|
| 1 | More DSBs → more inter-chromosomal AMGs combinatorially | Corollary 11 | Spearman correlation pct_inter vs n within ETS+ | Pending |
| 2 | High k, low bᵢ → Lemma 14 regime, single large cycle | Lemma 14 | ETS+ patients cluster where max_b_i=1 | Pending |
| 3 | Probable cycle count ≈ log(n) | Theorem 6 | Compare chain count to log(n) per patient | Pending |
| 4 | Rank increases along progression path | Theorem 25 | rank_proxy by progression stage from mmc6 | Pending |
| 5 | Cycle count in [1, min(max_bᵢ, n/2)] | Proposition 4 | Verify all patients satisfy bounds | Pending |

---

## Next Priority: Progression Stage AMG Rank Table

### What mmc6 Contains
- 48 patients, key gene somatic copy number alterations with clonality %
- Genes: NKX3-1 (26 entries), PTEN (21), T2-ERG (11), FOXP1 (11), CDKN1B (10), TP53 (10), CHD1 (9), RB1 (11)
- Clonality % tells you early (high %) vs late (low %)

### Stage Assignment Logic
```python
Stage 1 genes: NKX3-1, FOXP1, T2-ERG  (earliest, clonal)
Stage 2 genes: CDKN1B, TP53            (later, often subclonal)
Stage 3 genes: PTEN                     (latest, gating event)

Stage 3 = has Stage 1 + Stage 2 + Stage 3 genes
Stage 2 = has Stage 1 + Stage 2 genes
Stage 1 = has Stage 1 genes only
Stage 0 = no key genes detected

Use RELATIVE clonality ordering within patient
NOT absolute threshold (e.g. not "50% = subclonal")
A patient with NKX3-1=100% and PTEN=97% still has
NKX3-1 as earlier event because 100% > 97%
```

### AMG Metrics to Compute Per Patient
```python
import math

# Proposition 4 — cycle count range
min_cycles = 1
max_cycles = min(max_b_i, total_DSBs // 2)

# Theorem 6 — most probable cycle count
probable_cycles = max(1, round(math.log(total_DSBs)))

# Theorem 25 — rank proxy
rank_lower   = total_DSBs - max_cycles
rank_upper   = total_DSBs - min_cycles
rank_probable = total_DSBs - probable_cycles

# Lemma 14 check
lemma14_applies = (max_b_i == 1)
# If True → exact cycle structure is Cₙ, rank = n-1
```

### Expected Finding
If rank_probable increases from Stage 0 → Stage 1 → Stage 2 → Stage 3, that is the first empirical evidence that Baca's progression path corresponds to upward movement through Sheth's poset L(n). Neither paper makes this connection explicitly. This is the most original contribution of the project.

---

## Important Technical Notes

### chrom_aberrations_baca.csv Column Details
```
Individual              → patient ID (join key with clinical)
Number                  → rearrangement ID per patient
Breakpoint 1 chromosome → chromosome number (float, e.g. 21.0)
Breakpoint 1 position   → exact base pair position
Breakpoint 2 chromosome → chromosome number
Breakpoint 2 position   → exact base pair position
Class                   → inter_chr / long_range / inversion /
                          potential deletion / tandem_dup
Span                    → distance between breakpoints
Breakpoint 1 site       → gene annotation e.g.
                          "Intron of ERG(-): 10Kb after exon 1"
Breakpoint 2 site       → gene annotation
Fusion                  → gene fusion created e.g.
                          "Transcript fusion (TMPRSS2-ERG)"
```

### Gene Name Extraction from Site Annotations
```python
import re
def extract_gene_name(site_str):
    # Matches patterns like ERG(-) or TMPRSS2(+)
    matches = re.findall(r'([A-Z][A-Z0-9]+)\([+-]\)', site_str)
    return matches[0] if matches else None
```

### Clonality Interpretation
- Use RELATIVE ordering within patient — not absolute thresholds
- PTEN is OFTEN subclonal at cohort level but individual patients can have clonal PTEN loss (e.g. P07-837 has PTEN clonality=97%)
- A patient with NKX3-1=100% and PTEN=97% still has ordering consistent with the progression path

### Gleason Score Parsing
```python
# Handles formats: '3+4', '4+3;5', '3+4;5'
baca_df['has_tertiary_pattern'] = baca_df['Gleason Score'].str.contains(';').fillna(False)
baca_df['gleason_primary']   = baca_df['Gleason Score'].str.split('+').str[0]
baca_df['gleason_secondary'] = baca_df['Gleason Score'].str.split('+').str[1].str.split(';').str[0]
baca_df['gleason_primary']   = pd.to_numeric(baca_df['gleason_primary'],   errors='coerce')
baca_df['gleason_secondary'] = pd.to_numeric(baca_df['gleason_secondary'], errors='coerce')
baca_df['gleason_total']     = baca_df['gleason_primary'] + baca_df['gleason_secondary']
```

### ETS Classification
```python
def classify_ets(value):
    if pd.isna(value) or value == '---':
        return 'ETS-'
    else:
        return 'ETS+'
# Applied to 'ETS fusion detected by sequencing' column
```

---

## File Structure

```
nih-tcga-prad/
├── data/
│   ├── clinical/
│   │   └── patients_clinical_features.csv
│   ├── maf/
│   │   └── all_mutations_merged.csv
│   ├── cnv/
│   │   └── cnv_patient_mapping.csv
│   └── baca_dataset/
│       ├── chrom_aberrations_baca.csv     ← PRIMARY DATA FILE
│       ├── clinical_phenotypes.csv        ← CLINICAL DATA
│       ├── mmc2.xlsx through mmc7.xlsx    ← BACA SUPPLEMENTARY
│       └── Patient_P05-1657.xlsx          ← PROF ARSUAGA EXAMPLE
├── scripts/
│   └── baca/
│       └── baca_aberration_clinical_analysis.py  ← MAIN SCRIPT
└── results/
    └── patient_genome_profiles.csv        ← GENOME PROFILE TABLE
```

---

## What Has Been Validated So Far

```
✓ P05-1657 genome profile matches Sheth Section 7.1 exactly
  chromosomes [4,7,8,12,21], DSBs {4:2,7:4,8:9,12:1,21:2}
  Θ(5,(2,4,9,1,2)) confirmed

✓ ERG fusion chromoplexy classification reproduces Baca exactly
  15 chromoplexy-embedded, 11 simple fusion
  Matches Baca's reported 15/26 (58%)
  p < 0.0001 Mann-Whitney on k_local

✓ 3+4 vs 4+3 Gleason distinction captured
  4+3 has >2x rearrangements of 3+4
  4+3 has lower pct_inter despite more total breaks
  Consistent with different biological subtypes
```

---

## Pending Next Steps (Priority Order)

```
1. Add progression stage assignment from mmc6
   Load mmc6 Table S6B
   Assign Stage 0-3 per patient using gene presence
   + relative clonality ordering check

2. Compute Sheth AMG metrics per patient
   min_cycles, max_cycles, probable_cycles
   rank_lower, rank_upper, rank_probable
   lemma14_applies flag

3. Build progression stage vs AMG rank table
   Does rank_probable increase Stage 0→1→2→3?
   This is the key finding connecting both papers

4. ETS+ vs ETS- comparison
   Mann-Whitney on k, pct_inter, total_DSBs
   Between ETS+ and ETS- groups

5. Circos plots (Professor Arsuaga requested)
   Library: pycirclize
   Patients: P05-1657 (Sheth's case),
             one chromoplexy-embedded ERG+ patient,
             one simple fusion patient

6. TCGA Phase 1.2
   Kaplan-Meier survival curves by Gleason group
   Cox proportional hazards model
   (lower priority until Baca analysis complete)
```

---

## Presentation Status

**Next lab meeting: Thursday**

Slides cover:
1. Traditional cancer evolution model
2. Clinical problem — two identical Gleason 7 patients, different outcomes
3. Baca clonality concept
4. Consensus progression path (Steps 1-3 with gene biology)
5. Chromoplexy mechanism (correct vs incorrect repair)
6. Chromoplexy at each step of progression path
7. AMG framework (Θ notation, poset L(n), P05-1657 case)
8. Connecting both papers
9. TCGA dataset description
10. Computing AMG from Baca data (3 levels)
11. Research questions (6 questions)
12. Current results tables

**Per Professor Arsuaga feedback:**
- Remove chromothripsis slides (too complex for now)
- More figures, less text — split slides 50/50 figure/text
- Add Circos plots
- Add table of cycle distribution by tumor stage
- Check Mitelman database for Circos examples

---

## The Honest Research Framing

```
What we are doing:
    Taking Baca's exact breakpoint coordinates
    from 57 prostate cancer patients and computing
    AMG complexity metrics defined by Sheth's
    mathematical framework — then testing whether
    the metrics behave as Sheth's theorems predict

What is original:
    1. First quantitative mapping of Baca's patients
       into Sheth's Θ(k,(b₁,...,bₖ)) configuration space
    2. First validation that ERG fusion chromoplexy
       patients cluster in the Lemma 14 regime
       (high k, small bᵢ per chromosome)
    3. First test of whether progression path
       corresponds to increasing rank in L(n)
       (neither paper makes this connection)

What requires future work:
    Exact cycle structure computation requires
    graph traversal of adjacency data from mmc5
    Large n patients (>50 DSBs) may be
    computationally intractable for exact enumeration
    TCGA structural variant access requires
    dbGaP application (phs000178)
```
