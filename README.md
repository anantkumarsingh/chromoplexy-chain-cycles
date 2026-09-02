# Chromoplexy Chain Cycles

Genomic analysis of prostate adenocarcinoma (PRAD), studying **chromoplexy** —
coordinated, multi-chromosome chains of DNA rearrangement — through a formal
mathematical cycle-completion framework.

Conducted in the Arsuaga-Vazquez Lab at UC Davis, under Professor Javier
Arsuaga, as an extension of the AMG (Abstract Matching Graph) framework
described in Sheth, Arsuaga & Sazdanovic (2026), applied to the whole-genome
rearrangement data published in Baca et al. (2013, *Cell* 153:666-677).

## What This Project Is About

Baca et al. sequenced 57 prostate tumor genomes and identified chromoplexy —
clusters of rearrangements ("chains") that appear to have arisen from a
single coordinated mutational event rather than a series of independent
ones — and published chromoplexy calls for their cohort based on chain size
(number of rearrangements per chain).

This project re-analyzes the same underlying breakpoint data using the AMG
framework, which models a chromosome rearrangement chain as a graph of
alternating **rearrangement edges** (real, sequenced fusion junctions) and
**reference edges** (the genome's normal, un-rearranged connectivity). A
chain "closes" into a cycle when a path through this graph returns to its
starting point — a stronger, topological definition of coordination than a
simple rearrangement-count threshold. The central question: does closure
under this framework agree with Baca's own chain-size-based chromoplexy
calls, and what does copy-number and clinical data add to that picture?

A second, independent dataset — TCGA-PRAD (NIH GDC portal, ~500 patients) —
provides a complementary clinical/mutation/CNV pipeline for the same
disease, built and available for future analysis alongside the Baca-AMG
work.

## Tech Stack

- **Python 3.13**, run from a local virtual environment — a research
  pipeline, not a packaged library (no build system, package manager
  manifest, or test suite)
- **pandas** — tabular data processing across all pipelines
- **NetworkX** — graph construction, maximum matching, and cycle/component
  detection for AMG closure
- **Matplotlib** — all chain, cycle, and clinical-correlation diagrams
- **openpyxl** — reading Baca's supplementary Excel tables (`mmc1`–`mmc7`)
- **NIH GDC API** — programmatic mapping of TCGA copy-number files to
  patient barcodes

## Data Sources

- **Baca et al. 2013** raw rearrangement calls and all seven official
  supplementary tables (breakpoint-level chain assignments, chromoplexy
  summaries, genome-wide copy-number segmentation, clonality estimates)
- **TCGA-PRAD** clinical, somatic mutation (MAF), and copy-number (ASCAT3)
  data from the NIH GDC portal

## What Was Computed

- **AMG cycle-structure computation** on Baca's chromosomal rearrangement
  chains, including full enumeration of every combinatorially valid
  alternative pairing of breakpoints for chains small enough to enumerate
  exhaustively
- **Real-data cycle closure** using Baca's own three edge types (fusion,
  statistical adjacency, deletion bridge) as the primary substrate,
  cross-validated against an independent full-genome reconstruction
- **Obligate and most-probable cycle completion** for chains that don't
  close on real data alone — the conservative (fewest, smallest cycles) and
  most statistically frequent completions, respectively, following the
  Cornforth (2001) and Sheth (2026) definitions
- **Copy-number plausibility checks** on every invented completion edge
  against real genome-wide segmentation data, calibrated against Baca's own
  validated deletion-bridge calls
- **Deletion-bridge structural analysis**, tracing the exact final
  rearranged DNA structure implied by every deletion-bridge pair in the
  cohort to determine whether the bridged material is truly lost, moved,
  or unchanged
- **Per-patient chromoplexy reclassification**, rolling per-chain results
  up to the patient level and comparing directly against Baca's published
  headline statistic
- **Statistical correlation analysis** across ETS fusion status, chain
  size, Gleason score, and pathological stage against every closure-based
  chromoplexy definition
- Full visualization suite: chain diagrams, individual closed-cycle
  diagrams, and completion-structure diagrams for every patient and every
  chain in the cohort

## Key Results

- Real, physically-closed AMG cycles were found in only **4 of the 50
  patients (8%)** Baca classifies as chromoplexy-positive by chain size —
  obligate and most-probable completion recover 21/50 and 24/50
  respectively, still under half even under the most permissive
  reconstruction
- **ETS fusion-positive status is strongly associated with AMG cycle
  closure** (odds ratios of 4–15 across multiple tests), while Baca's own
  chain-size criterion shows no such association — a result not reported
  in the original paper
- **Deletion bridges represent relocated, not deleted, DNA**: across all
  521 deletion-bridge pairs in the cohort, 0 were found to be genuinely
  unaccounted for anywhere in the data; 94% showed evidence of relocation
  or reorientation
- Invented completion edges are overwhelmingly consistent with real
  copy-number data — only 1.3% of all invented edges cohort-wide
  contradict the measured copy-number profile at that location
- No significant association was found between AMG closure and Gleason
  score or pathological tumor stage
- Baca's own published figure of 121 closed chromoplexy chains could not
  be reconciled with any per-chain or per-patient identifying information
  in the paper's own supplementary tables — confirmed by direct inspection
  of the raw data and every relevant figure panel

## Presentations

Slide decks presented on this project are in [`presentations/`](presentations/):

- [Baca-AMG-TCGA-Connection.pdf](presentations/Baca-AMG-TCGA-Connection.pdf)
- [Cycle-Structure-Generation.pdf](presentations/Cycle-Structure-Generation.pdf)
- [Completing Cycle Structures.pdf](presentations/Completing%20Cycle%20Structures.pdf)
- [Deletion Bridges.pdf](presentations/Deletion%20Bridges.pdf)

## Repository Layout

- `scripts/baca/` — all analysis code for the Baca-AMG track, organized by
  phase and by named sub-track
- `scripts/tcga/` — the TCGA-PRAD clinical/mutation/CNV pipeline
- `data/` — raw and organized input data for both tracks
- `results/` — every generated CSV output and diagram, organized to mirror
  the scripts that produced them
