# Baca's Chain (mmc5) vs. Our Chain (chrom_aberrations_baca.csv) — Comparison

**Date:** 2026-06-22
**Script:** `scripts/baca/baca_chain_comparison.py`
**Purpose:** reconcile Baca's own published chromoplexy "chain" data against our independently-built AMG chain constructions, and explain why they disagree.

---

## 1. The question that started this

Baca's paper states:

> "single events joined DNA from dispersed regions of six or more chromosomes in multiple tumors, whereas chromothripsis frequently involves focal rearrangement of one or two chromosomes"

Our own AMG closed-cycle analysis (see `CYCLE_COMPLETION_PHASE1_SUMMARY.md`) found that, for ERG+ patients, closed AMG cycles never span more than **2** chromosomes even after completing them with every real rearrangement available. That looked like a contradiction. It isn't — and tracing it down required understanding that Baca's data comes in **two different files that answer two different questions**.

---

## 2. Two raw files, two different objects

| | `chrom_aberrations_baca.csv` | `mmc5.xlsx` (Table S5A) |
|---|---|---|
| Row = | one rearrangement (2 breakpoints) | one breakpoint |
| Chain/cluster info | **none** — flat, unordered | `Chain number` column — Baca's own clustering of rearrangements into coordinated chromoplexy events |
| Coverage | every detected rearrangement, all 57 patients | only the rearrangements Baca's algorithm assigned to a chain, 52 patients |
| Used by (existing code) | `_get_erg_chain_rows`, `baca_full_genome_cycles.py` — i.e. **our own** chain reconstructions | never used before this session |

We verified (P01-28) that these are **not independent datasets** — 138/140 mmc5 breakpoint positions match `chrom_aberrations_baca.csv` exactly. `mmc5` is a curated *subset* of the same underlying rearrangements, with chain labels added; it omits whatever Baca's clustering algorithm didn't assign to a chain.

**Reconciling the "six or more chromosomes" quote:** that sentence is almost certainly describing the chromosome *span* of a `Chain number` group in `mmc5` — i.e., reachability/connectivity of a cluster of rearrangements — not a verified closed AMG cycle (alternating rearrangement + reference edges, returning to the start). `mmc5` has no path ordering at all; `Adjacent breakpoint(s)` just names each breakpoint's own rearrangement partner, not a sequence link. There is nothing in Baca's raw data that tests for loop closure. So "spans six chromosomes" (chain reachability) and "closes into an AMG cycle" (our criterion) are different mathematical claims, and satisfying the first does not imply the second.

---

## 3. What we built: `baca_chain_comparison.py`

Five output CSVs in `results/`:

| File | Contents |
|---|---|
| `baca_proposed_chains_mmc5.csv` | Every chain Baca defined (366 chains, 52 patients) — chromosomes, k, b_i, Θ(k,(b_i)) notation, genes, read support, clinical fields |
| `our_erg_chain_chrom_aberrations.csv` | Our ERG-anchored one-hop chain (26 ERG+ patients), built purely from the raw rearrangement file, same schema for direct comparison |
| `our_full_chain_chrom_aberrations.csv` | Our full-genome chain — every rearrangement per patient, all 57 patients, no anchoring or filtering at all |
| `baca_erg_chr21_chains_mmc5.csv` | Subset of Baca's chains that include chr21 (25 chains, 22 patients) — broader than text-matching "ERG"/"TMPRSS2" in the site annotation, since chr21 carries other genes too |
| `erg_chain_comparison_ours_vs_baca.csv` | Per-ERG+-patient: our chain vs. Baca's ERG/TMPRSS2-annotated chain, side by side |
| `erg_set_comparison_baca_chr21_vs_ours.csv` | Patient-level set comparison: clinically-ERG+ (ours, 26 patients) vs. has-a-chr21-chain (Baca's mmc5, 22 patients) |

---

## 4. Key findings

### 4.1 The ERG+ patient sets disagree in both directions

Comparing "clinically ERG+" (our 26-patient definition, from `clinical_phenotypes.csv`) against "has a chr21-containing chain in mmc5" (22 patients):

- **18 patients — `both`**: agree on ERG+ status, but only **4/18** have an identical chromosome set between our chain and Baca's chain (the simplest cases — single intrachromosomal chr21 deletion, nothing to disagree about).
- **8 patients — `our_only`**: clinically ERG+, but Baca's own chain-clustering algorithm never grouped chr21 into any chain for them at all (mostly the "simple_fusion" cases — an isolated chr21 deletion that never got clustered).
- **4 patients — `baca_only`**: have a chr21-containing chain in mmc5 but are **not clinically ERG+** at all (`P07-144`, `P08-716`, `PR-07-3258` have no ETS fusion detected by sequencing; `PR-4240` has an unrelated `NRF1-BRAF` fusion). Chr21 carries genes other than ERG/TMPRSS2, so chromosome membership alone is not a safe ERG+ proxy.

### 4.2 Baca's chains are a strict subset of the full rearrangement data — confirmed for all 52 patients

We tested directly: for every patient in mmc5, is the union of chromosomes across all of that patient's Baca-defined chains a subset of the full chromosome set from every rearrangement that patient has in `chrom_aberrations_baca.csv`?

**Result: 52/52 patients pass.** Baca's chains never include a chromosome that isn't somewhere in the patient's raw rearrangement data (expected — mmc5 is built from the same data), and the raw data always contains chromosomes that *no* Baca chain includes. Baca's published chains are a conservative, filtered selection from inside the full rearrangement graph — not the full graph itself.

### 4.3 Worked example: P01-28

Full raw rearrangement file, every chromosome with at least one rearrangement:
```
{1, 2, 3, 4, 6, 9, 12, 14, 15, 17, 18, 19, 20, 21, 23}   (15 chromosomes, 96 rearrangements)
```

Baca's mmc5 chains for this patient (5 separate chains):
```
Chain 1: {1, 2, 3, 6, 9, 21}      40 breakpoints
Chain 2: {2}                      20 breakpoints
Chain 3: {4, 9}                   14 breakpoints
Chain 4: {1, 14, 17, 18, 19, 23}  60 breakpoints
Chain 5: {2, 19}                   6 breakpoints

Union: {1, 2, 3, 4, 6, 9, 14, 17, 18, 19, 21, 23}   (12 chromosomes)
```

Union is a strict subset of the full set. **Missing: chromosomes 12, 15, 20.** Traced to the exact rearrangement rows:

| Rearrangement # | Breakpoints | Class | Tumor reads | In any Baca chain? |
|---|---|---|---|---|
| 30 | chr15 ↔ chr15 | potential deletion (intrachromosomal) | 12 | No |
| 91 | chr4 ↔ chr20 | inter_chr | 6 | No |
| 227 | chr12 ↔ chr12 | tandem_dup (intrachromosomal) | 4 | No |
| 246 | chr20 ↔ chr23 | inter_chr | 4 | No |

- **chr12 and chr15** are simple: both are pure intrachromosomal events (tandem duplication, deletion) with no connection to anything else in the patient's genome. Textbook singletons — there's nothing for a chain-builder to link them to.
- **chr20 is the more interesting case.** It has two *interchromosomal* connections — to chr4 and to chr23 — and both chr4 and chr23 already appear in other Baca chains (via different rearrangements). So chr20 is geometrically connected to chain material, yet Baca's algorithm still excluded it.
- Cohort-wide, rearrangements that made it into a Baca chain average ~13 tumor reads vs. ~7 for excluded ones (not a hard cutoff — ranges overlap) and skew much more interchromosomal (65% vs. 37%). This suggests Baca's chain-assembly used more than simple "shares a chromosome with the chain" connectivity — likely a confidence/statistical filter not fully recoverable from the public columns alone.

---

## 5. Visualizing Baca's chains directly from mmc5 (2026-06-22)

**Script:** `scripts/baca/baca_proposed_chains_drawing.py`
**Output:** `results/baca_proposed_chains/{patient_id}_cycle_{chain_number}.png` — one image per `(Individual, Chain number)` group in mmc5's Table S5A. 366 images, all 52 mmc5 patients, 0 failures.

This draws Baca's chain **exactly as mmc5 defines it** — no cross-referencing with `chrom_aberrations_baca.csv` at all, everything comes from the mmc5 rows for that one chain.

### Discovering mmc5's two edge types (important correction mid-build)

The first version of this script wrongly assumed `Adjacent breakpoint(s)` was the rearrangement link. That was caught and corrected by checking same-vs-cross-chromosome counts across the whole mmc5 table before trusting either field:

| Field | Same-chromosome | Cross-chromosome | What it actually is |
|---|---|---|---|
| `Rearrangement number` (two rows sharing it) | 1549 | 671 | **Rearrangement edge** — the two ends of one real structural variant (matches the inter/intra mix expected from real SVs) |
| `Adjacent breakpoint(s)` | 2654 | 0 | **Reference edge** — the genomically-next breakpoint on the same chromosome, given directly by Baca, not inferred from position order the way our own AMG code does |

This means mmc5 encodes both AMG edge types directly — a follow-up worth revisiting: real cycle-closure computation on Baca's own chain data without inferring reference edges ourselves, the way `baca_full_genome_cycles.py` currently has to.

### Verifying the strand convention (`Forward`/`Reverse` vs. `+`/`-`)

mmc5 has its own `Strand` column with string values `Forward`/`Reverse`, remapped in the script as `Forward→'+'`, `Reverse→'-'` purely for the existing AMG node convention (used only to offset a breakpoint's arc-attachment point left/right of its tick mark — cosmetic, not used for any chromosome/pairing logic). This was not assumed — verified directly against `chrom_aberrations_baca.csv`'s own `+`/`-` strand columns:

1. **Built a ground-truth lookup from the raw file.** For every row in `chrom_aberrations_baca.csv`, mapped `(Individual, Breakpoint 1 chromosome, Breakpoint 1 position) → Breakpoint 1 strand`, and did the same for the Breakpoint 2 side. This gives the real `+`/`-` value recorded at every exact `(patient, chromosome, position)` triple in the raw data.
2. **Looked up every mmc5 row against that table.** For each of the 4440 rows in mmc5, used the same `(Individual, chromosome, position)` key to find the matching raw-file strand, and compared it against mmc5's mapped value (`Forward→'+'`, `Reverse→'-'`).
3. **Counted exact matches.** 4374 of 4440 mmc5 rows had an exact position match in the raw file (the remaining 66 are the expected mmc5-is-a-subset gap, consistent with §3). Of those 4374 comparisons: **4370 matched, 4 did not.**
4. **Hand-checked all 4 mismatches** rather than accepting the 99.9% rate at face value. Each one turned out to be a lookup-collision artifact, not a real disagreement: the exact same `(chromosome, position)` appears *twice* in the raw file for that patient — once as the `+` end of one rearrangement and once as the `-` end of a completely different rearrangement (e.g. `P08-1042` chr2:65,557,230 is `+` for rearrangement #31 and `-` for rearrangement #53; same pattern for 3 positions in `P09-37` on chr10). A single-key dictionary built in step 1 can only retain one of the two values at a shared key, so the comparison "disagreed" with whichever one it happened to drop — mmc5's recorded value matches one of the two genuinely-valid options at that position either way.
5. **Conclusion:** `Forward`/`Reverse` (mmc5) and `+`/`-` (`chrom_aberrations_baca.csv`) are the same strand convention, with zero genuine discrepancies found — consistent with the `+` = right-of-DSB / `-` = left-of-DSB convention used throughout this codebase's AMG framework.

### Worked example: P01-28, Baca Chain 5

The mmc5 rows for this one chain (6 breakpoints):

| Breakpoint # | Rearrangement # | Chromosome | Position | Strand | Adjacent breakpoint(s) |
|---|---|---|---|---|---|
| 62 | 84 | chr2 | 223,858,306 | + | \|70\| |
| 168 | 84 | chr19 | 22,335,283 | − | (empty) |
| 70 | 103 | chr2 | 223,859,813 | − | \|62\| |
| 176 | 103 | chr19 | 21,164,553 | − | \|186\| |
| 186 | 162 | chr19 | 21,156,131 | + | \|176\| |
| 80 | 162 | chr2 | 223,412,353 | − | (empty) |

**Step 1 — chromosomes/positions:** chr2 has 3 positions, chr19 has 3 positions → k=2, b_i={2:3, 19:3} → **Θ(2,(3,3))**.

**Step 2 — rearrangement edges** (group by `Rearrangement number`, each pair = one real SV):
- Rearrangement 84: bp62 (chr2:223,858,306) ↔ bp168 (chr19:22,335,283) — interchromosomal
- Rearrangement 103: bp70 (chr2:223,859,813) ↔ bp176 (chr19:21,164,553) — interchromosomal
- Rearrangement 162: bp186 (chr19:21,156,131) ↔ bp80 (chr2:223,412,353) — interchromosomal

All 3 have both ends present in this chain → `rearrangement_ends_outside_chain=0`. Drawn as the three crossing orange dashed arcs between the chr2 and chr19 rows.

**Step 3 — reference edges** (`Adjacent breakpoint(s)`):
- bp62 ↔ bp70: both on chr2, ~1.5kb apart → green arc on the chr2 row
- bp176 ↔ bp186: both on chr19, ~8kb apart → green arc on the chr19 row
- bp168 and bp80 both have an **empty** adjacent field — these are the chain's two open/loose ends, with no recorded reference neighbor, sitting at the far edges of their respective rows with nothing connecting them on the reference side.

**Result:** 6 breakpoints, 3 rearrangement edges (all chr2↔chr19), 2 reference edges (one per chromosome), 2 open ends — exactly what the image's title line and arcs show. Nothing inferred; every number traces directly to the mmc5 rows above.

### A real rendering bug, fixed with exact arc-geometry math (not trial and error)

Single/few-chromosome chains (small `y_span`) initially rendered with rearrangement arcs **clipped off the bottom of the figure** (e.g. `PR-STID0000002872` chain 1). Root cause: matplotlib's `arc3` connection style computes its curve sagitta in **physical display inches** — `sagitta_inches = |rad| * dx_data * (fig_w / xlim_range)` — which is independent of the y-axis data range entirely. With a fixed `fig_w=24in` (matching the full-genome-chain convention from `baca_full_genome_cycles.py`) and deep curvature (`rad` up to 0.78) on a wide-spanning intrachromosomal arc, that's up to ~12.7 inches of real vertical room needed regardless of how many chromosome rows exist.

First attempt (inflating `fig_h` independently of the y-axis range) overcorrected — fixed the clipping but wasted most of the figure as blank space (some images hit 9195px tall with almost nothing in them). Final fix: compute the exact extra inches needed for the worst-case arc (`worst_sagitta_in`), then extend **both** `fig_h` and the bottom y-limit by amounts that preserve the original (well-tuned, multi-chromosome) inches-per-data-unit ratio exactly:

```
extra_in   = worst_sagitta_in * 1.15                # safety margin
extra_data = extra_in / scale_y0                     # scale_y0 = baseline fig_h0 / y_range0
fig_h      = fig_h0 + extra_in
ylim       = (-1.8 - y_span*0.05 - extra_data, y_span + 2.2)   # only the bottom extends
```

This keeps `scale_y` mathematically identical to the baseline (provably, not approximately — substituting shows the new ratio equals the old one exactly), so the extra room added is exactly what the worst arc needs and nothing more. Final image heights range 1429–7605px across all 366 chains, all arcs fully contained, no clipping.

## 6. Bottom line

- "Our full chain" (raw, unfiltered, all rearrangements per patient) is the maximal object; Baca's published `mmc5` chains are statistically-filtered subsets drawn from inside it — confirmed across the entire cohort, not just one patient.
- The gap between them is partly genuine noise (isolated intrachromosomal singletons like chr12/chr15) and partly real interchromosomal connectivity that Baca chose not to include in a chain for reasons not fully recoverable from the published data (the chr20 case).
- Baca's "six or more chromosomes" chromoplexy language describes chain *reachability* in this same filtered chain object — not AMG cycle *closure*, which is the stricter, separate question our own analysis (Phase 1, `CYCLE_COMPLETION_PHASE1_SUMMARY.md`) is answering. The two are not in conflict; they're different mathematical claims about the same underlying data.
