# Phase 2 Copy-Number Plausibility Check

**Script:** `baca_phase2_cn_plausibility.py`
**Outputs:**
`results/phase2_cn_deletion_bridge_calibration.csv`,
`results/phase2_cn_invented_edge_detail.csv`,
`results/phase2_cn_chain_summary.csv`

## 1. What this is, and why we're doing it

Prof. Arsuaga's directive for this whole cycle-completion project has four
parts (see `CLAUDE.md`, "Next Analysis Direction" / Phase 1.5): (1) complete
the cycles, obligate/minimal first, (2) **track copy number and report what
extra information the completed cycles provide**, (3) check swap-stability
(Phase 4, not yet started), (4) defer biological constraints to a later
paper. Item (1) is done (Phase 2: `baca_phase2_obligate_probable_
completion.py`). Item (2) has a real-edge half that was already done in
Phase 1.5 (using Baca's deletion-bridge edges to avoid a copy-number-blind
reference rule) — but Phase 2's **invented** edges (the obligate/probable
completions that close open chain ends) have never had any copy-number
check at all. This script is that missing half.

This became possible this session because `data/baca_dataset/mmc3.xlsx`
Table S3B was discovered to be a real, genome-wide copy-number segmentation
file (19,493 segments, all 57 patients, standard CBS/DNAcopy-style output:
`Sample, Chromosome, Start, End, Segment_Mean`) — previously the project
believed no such file existed and only had `mmc6.xlsx`'s 14-gene proxy.

**The question this answers:** when Phase 2 invents a connection between
two already-real, already-observed dangling breakpoint ends to close a
cycle, is that connection's simplest reading (undisturbed reference DNA
between the two ends) consistent with the measured copy-number profile
there, or does it silently cross a copy-number transition that nothing in
this patient's real data explains? That's the concrete, citable "extra
information" the professor asked for.

**Important framing point:** Phase 2 never invents new breakpoint
*positions* — every invented edge connects two breakpoints that already
exist in the real mmc5 data (see `baca_phase2_obligate_probable_
completion.py`'s own docstring, point 3). So this script is never asking
"is this position real" — it always is. It is asking whether the specific
*connection* between two real positions is plausible given real copy-number
evidence.

## 2. What this does

1. **Calibrate the copy-number lookup on a known-good positive control**
   before trusting it on anything invented: Baca's own `mmc5.xlsx` Table
   S5A `Deletion bridge partner breakpoint` column records real,
   statistically-validated pairs of breakpoints known to border the *same*
   deleted region. If the lookup mechanism works, these pairs should show a
   measured copy-number **loss** between them.
2. **Build a position → copy-number-segment lookup** from Table S3B, by
   interval containment (not exact match — confirmed earlier this session
   that real breakpoints essentially never land exactly on a segment
   boundary; median distance ~157kb, because CN segmentation and
   rearrangement-breakpoint detection are different measurements at very
   different resolutions).
3. **Re-derive, for every one of the 354 enumerable Phase 2 chains, every
   distinct set of invented edges that achieves the selected obligate/
   probable structure** (not just one arbitrarily chosen realization — Phase
   2's own enumeration can tie multiple different edge-sets to the same
   resulting cycle-length structure; this script checks for that explicitly
   rather than assuming it away).
4. **Classify every invented edge** into one of five categories (below),
   using only real, already-recorded data — never inventing a new
   breakpoint or a magnitude threshold not justified by the data.
5. **Aggregate to chain level**, and cross-reference against Phase 3's
   patients who are chromoplexy-negative under `real` but positive under
   `obligate`/`probable` — the specific population the professor's question
   is actually about.

### The five-way classification

- **CN_FLAT** — both endpoints resolve to a real segment, and it's the
  *same* segment (no copy-number transition anywhere between them).
  Neutral: consistent with undisturbed reference DNA, but not proof (a
  balanced rearrangement produces no CN signal either).
- **CN_TRANSITION_EXPLAINED** — a real segment boundary exists between the
  two endpoints, but at least one real, already-observed breakpoint from
  this patient's **full** rearrangement dataset (`chrom_aberrations_baca.csv`,
  not just this chain) falls strictly inside that span — there's already a
  recorded real event that could account for the shift.
- **CN_TRANSITION_UNEXPLAINED** — a real segment boundary exists between the
  two endpoints, and nothing in this patient's real data (anywhere in the
  genome) falls in that span. The strongest red flag this analysis can
  raise.
- **CN_INDETERMINATE** — at least one endpoint has no resolvable segment (a
  genuine S3B coverage gap — ~26% of real breakpoints fall in such gaps,
  confirmed earlier this session). Not evidence for or against the edge —
  we simply cannot assess it. Never folded into CN_FLAT.
- **NOT_APPLICABLE_INTERCHROMOSOMAL** — the two endpoints are on different
  chromosomes. S3B segments only tile one chromosome each; there is no
  shared coordinate span to check. Reported as its own honest category
  rather than inventing a weaker interchromosomal proxy that wasn't agreed.

## 3. CSV column reference

### `phase2_cn_deletion_bridge_calibration.csv` (521 rows — every real deletion-bridge pair)

| Column | Meaning |
|---|---|
| `patient_id`, `breakpoint_a`, `breakpoint_b` | the real, Baca-validated deletion-bridge pair (mmc5 Table S5A) |
| `chromosome`, `pos_a`, `pos_b` | their shared chromosome and positions (deletion bridges are always intrachromosomal — confirmed, 0/1042 directed links cross chromosomes) |
| `both_resolve` | True if both positions found a containing S3B segment |
| `has_transition` | True if the two positions fall in different S3B segments |
| `shows_expected_loss` | True if a transition exists AND the lowest segment mean strictly between them is below the higher of the two flanking means — i.e. a real measured loss |
| `seg_mean_a`, `seg_mean_b`, `min_between_mean` | the raw log2-ratio values, for inspection |

### `phase2_cn_invented_edge_detail.csv` (195,032 rows — one per invented edge, per chain, per lens, per tied realization)

| Column | Meaning |
|---|---|
| `patient_id`, `baca_chain_number` | which mmc5 chain |
| `lens` | `obligate` or `probable` |
| `realization_index`, `n_realizations` | which tied edge-set this row belongs to, out of how many equally-valid ones achieve the selected structure |
| `edges_unique` | True if only one distinct edge-set achieves the selected structure for this chain+lens (obligate: always True, 0/354; probable: 194/354 False) |
| `breakpoint_a`, `breakpoint_b` | the two real breakpoints this invented edge connects |
| `chromosome_a/b`, `position_a/b` | their locations |
| `classification` | one of the five categories above |
| `seg_mean_a`, `seg_mean_b` | the resolved segment means, where available |

**Caveat, important for reading this file**: it is **realization-weighted** — a chain with 3,840 tied realizations contributes 3,840× as many rows as a chain with 1. Do not read raw counts/percentages from this file as a per-chain average; use `phase2_cn_chain_summary.csv` for that, or filter to `realization_index == 0` for one representative row per chain+lens.

### `phase2_cn_chain_summary.csv` (720 rows — one per chain, or per chain×lens where enumerable)

| Column | Meaning |
|---|---|
| `patient_id`, `baca_chain_number` | which chain |
| `phase2_enumerable` | False for the 12 chains with >12 dangling ends (no Phase 2 completion was computed for these at all; no CN analysis possible) |
| `lens` | `obligate` or `probable` (only present when enumerable) |
| `n_invented_edges` | how many invented edges this chain's selected structure requires (0 for the 52 chains already fully closed by real data) |
| `n_realizations`, `edges_unique` | same meaning as in the edge-detail file |
| `cn_flat_min/max`, `cn_transition_explained_min/max`, `cn_transition_unexplained_min/max`, `cn_indeterminate_min/max`, `not_applicable_interchromosomal_min/max` | for each classification, the minimum and maximum count found across all tied realizations (min==max when `edges_unique` is True) |
| `has_unexplained_in_every_realization` | True only if `cn_transition_unexplained_min > 0` — every equally-valid way of realizing this structure includes an unexplained edge |
| `has_unexplained_in_some_realization` | True if `cn_transition_unexplained_max > 0` — at least one equally-valid realization includes an unexplained edge (the more permissive, "worst case" reading) |

## 4. Step-by-step: what the script does

1. Load `mmc3.xlsx` Table S3B, group into per-`(patient, chromosome)` sorted segment tables.
2. Load `mmc5.xlsx` Table S5A (breakpoints) and `chrom_aberrations_baca.csv` (every real rearrangement, for the "explained" check).
3. **Calibration**: for every real deletion-bridge pair, look up both endpoints' segments, check for a transition and whether it's a loss. Write `phase2_cn_deletion_bridge_calibration.csv`.
4. For each of the 366 mmc5 chains: rebuild the real combined graph and its open ends exactly as Phase 2 does (reusing `_build_real_combined_graph`, `_closed_and_open_components` — not reimplemented), then re-run Phase 2's own enumeration (`_enumerate_completions`, `_select_obligate`, `_select_probable`) — this time additionally keeping each candidate's actual matching (an additive, no-op change to `baca_phase2_obligate_probable_completion.py`, verified byte-identical output on the existing pipeline before trusting anything downstream).
5. For each lens, collapse the achieving set down to its **distinct** edge-sets (some tied candidates use the literal same edges, just listed in a different order — deduped via `frozenset`).
6. Classify every edge of every distinct realization via the five-way scheme.
7. Aggregate to chain level (min/max across realizations) and write `phase2_cn_invented_edge_detail.csv` / `phase2_cn_chain_summary.csv`.
8. Cross-reference against `phase3_patient_chromoplexy_summary.csv`'s `real_has_chromoplexy_strict` vs. `obligate_/probable_has_chromoplexy_strict` to identify the specific patients whose chromoplexy call only exists *because* of invented completion, and report how many of those rest on a CN-unexplained edge.

## 5. Limitations and notes

- **Calibration coverage is partial by necessity**: only 79/521 (15%) real deletion-bridge pairs have both endpoints resolvable in S3B — the same ~26%-of-breakpoints coverage gap already documented (centromere/telomere/low-mappability regions S3B's segmentation doesn't reach). The calibration is still a clean, strong signal (59/59 = 100% of the resolvable-and-transitioning pairs show the expected loss), but it cannot speak to whether the lookup would behave the same way in the un-resolvable 85%.
- **`CN_FLAT` is not proof of correctness, and `CN_TRANSITION_UNEXPLAINED` is not proof of error.** A balanced rearrangement (simple translocation/inversion, no net gain or loss) produces no CN signal at all, so a flat profile is consistent with — but doesn't confirm — the invented connection. Conversely, an unexplained transition means the *simplest* reading (undisturbed reference DNA) is contradicted; it does not rule out a real, just-not-yet-recorded structural event there.
- **`NOT_APPLICABLE_INTERCHROMOSOMAL` (19.6% of invented edges) is a real, structural limitation of this method, not a gap to be closed later** — there is no meaningful same-coordinate CN check across two different chromosomes, and no interchromosomal proxy was built, deliberately, to avoid inventing a weaker heuristic never agreed to.
- **Edge-choice ambiguity is real for `probable` (194/354 chains) and essentially absent for `obligate` (0/354)** — this mirrors Test 6's earlier finding that Phase 2's tie-breaking concerns rarely bite in this dataset, now extended to a new metric (specific edge identity, not just structure or span) with the same clean result for obligate.
- **The "explained" check uses this patient's entire real rearrangement dataset**, not just this chain or even just the mmc5-covered subset — the most generous possible test, giving every invented edge its best chance of being explained before being flagged.
- **No telomere/centromere or driver-gene weighting is applied here** — consistent with the rest of Phase 2, these biological constraints remain explicitly deferred to the later, separate biological-constraints paper.

## 6. Analysis Corner

1. **The calibration step passed cleanly**: of 79 real, Baca-validated deletion-bridge pairs with resolvable copy-number data, 59 show a measurable transition, and *all 59* show it in the expected direction (a loss). This is strong confirmation the lookup mechanism itself is sound before trusting it on anything invented.
2. **The headline result: CN-unexplained invented edges are rare.** Using one representative realization per chain (the fair, non-realization-weighted view): of 1,368 total invented edges across both lenses, only **18 (1.3%)** are CN_TRANSITION_UNEXPLAINED. The rest split roughly evenly between CN_FLAT (31%), CN_INDETERMINATE (27%), CN_TRANSITION_EXPLAINED (21%), and NOT_APPLICABLE_INTERCHROMOSOMAL (20%).
3. **Obligate completions are essentially always CN-plausible**: of the 302 chains that actually require an invented edge, only **4 (1.3%)** have any CN-unexplained edge — and because obligate's edge choice is always unique (0/354 ambiguous), this is a clean, deterministic number, not a range.
4. **Probable completions carry more genuine risk, but mostly in a "possible, not certain" sense**: 3/302 chains have an unexplained edge in *every* tied realization (a firm flag), but up to 35/302 (11.6%) *could* have one depending on which equally-frequent realization is picked — reflecting probable's much higher edge-choice ambiguity (194/354 chains) relative to obligate's essentially none.
5. **The number that actually answers the professor's question**: among Phase 3's chromoplexy "recoveries" — patients who are chromoplexy-negative under real closure but positive under obligate/probable — **16 of 17 (94%) obligate recoveries have no CN-unexplained edge anywhere in their completion**, and even for the more liberal probable lens, only 2 of 20 (10%) recoveries are unexplained in *every* tied realization (12/20 could be, in a worst-case realization, but 18/20 have at least one equally-valid, CN-consistent way to achieve the same result). **Bottom line: the extra information copy number provides is that Phase 2's chromoplexy recoveries — especially the conservative obligate ones — are overwhelmingly not contradicted by real copy-number evidence.** This doesn't prove they're correct (CN_FLAT and CN_INDETERMINATE aren't proof either), but it substantially narrows where genuine doubt should be concentrated: a small, specific, now-identifiable set of chains/patients, rather than the whole recovered population.
6. **Next steps**: (a) the ~39 specific (chain, lens) pairs carrying an unexplained edge (`phase2_cn_invented_edge_detail.csv`, filtered to `classification == 'CN_TRANSITION_UNEXPLAINED'`) are a concrete, small, worth-inspecting-by-hand list before the paper write-up — each is a candidate for either "genuinely doubtful completion" or "a real event Baca's own algorithms happened to miss entirely." (b) This is a natural, real-data-only stopping point for the copy-number side of item 2 — Phase 4 (swap-stability) is the next item on the professor's original list. (c) The `mmc6.xlsx` 14-driver-gene proxy was not used here (S3B supersedes it for this purpose) but remains available as a targeted cross-check specifically on driver genes (ERG, PTEN, CDKN1B, TP53, NKX3-1) if the eventual paper wants a gene-specific angle, per the separately-deferred biological-constraints track.
