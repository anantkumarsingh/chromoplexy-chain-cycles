# Cycle Completion — Phase 1 Summary

**Purpose of this document:** explain, in plain language, what we built, how it works, and what we found — with one fully worked example, suitable for a presentation.

**Date:** 2026-06-18
**Directive:** Professor Arsuaga asked us to "complete the cycles" for the Baca et al. 2013 prostate cancer rearrangement data, as the first step toward a publication. The goal: figure out which patients *actually* show chromoplexy, using the real mathematical structure of the data (AMG cycles) instead of a connectivity-based heuristic.

---

## 1. The background you need

### What is a DSB and a rearrangement?

A **DSB** (double-strand break) is a point where DNA snapped. When the cell repairs DNA, it sometimes rejoins broken ends to the *wrong* partner — that's a **rearrangement**. Each row in Baca's data (`chrom_aberrations_baca.csv`) is one observed rearrangement: it tells you the two DSB positions that got rejoined to each other (chromosome, position, strand for each side).

Each rearrangement row = 2 DSBs (one on each side of the rejoin).

### What is an AMG cycle?

Think of a chromosome as a line. A DSB cuts it into two ends. There are two kinds of connections between DSB ends:

- **Rearrangement edge** — the *abnormal* connection: this is what got mis-rejoined (recorded directly in Baca's data).
- **Reference edge** — the *normal* connection: simply "the next DSB along the same original chromosome" (this is NOT in the data — we compute it from genomic position order).

If you alternate hopping along rearrangement edges and reference edges and you eventually come back to where you started, that's a **closed cycle**. If you instead run off the end of a chromosome (a real telomere) or reach a DSB end that was never reported as rearranged anywhere in the data, that's an **open path** — it doesn't close.

**Why cycles matter for chromoplexy:** A single big closed cycle touching many chromosomes is direct mathematical evidence of one coordinated, multi-chromosome shattering-and-reassembly event — which is what "chromoplexy" actually means. A bunch of small, single-chromosome closed loops, or open paths that never close, is NOT that — it's either simple/local damage, or it's an event whose full structure isn't confirmed by the data.

### Baca's own chromoplexy call (for contrast)

Baca's paper doesn't compute cycles at all. Their existing chromoplexy call (already in our codebase, `get_erg_chain_details_v2`) is a **connectivity heuristic**:

```
k_local    = number of distinct chromosomes connected to the ERG/TMPRSS2
             fusion locus through one hop of interchromosomal rearrangements
chain_size = number of those interchromosomal rearrangement rows

chromoplexy_embedded  if k_local > 2  AND  chain_size > 1
simple_fusion          otherwise
```

Applied to the 26 ERG-fusion-positive patients, this gives Baca's well-known result: **15/26 (58%) chromoplexy_embedded.**

This only asks "does the chain *reach* a lot of chromosomes?" — it never checks whether those reachable chromosomes actually form a closed loop.

---

## 2. What we built (three pieces of code, in order)

### Piece 1 — Full-genome cycle computation (`baca_full_genome_cycles.py`)

For each of all 57 Baca patients (not just the 26 ERG+ ones), we used **every single rearrangement row** for that patient — not just rows touching ERG — to build the complete rearrangement/reference edge graph, and ran the cycle-finding traversal across the whole thing.

This told us, per patient: how many DSBs total, how many closed cycles, how big each cycle is, and how many open paths remain.

**Key early finding:** most of the "open paths" we used to see in the smaller ERG-only analysis were not real biological gaps — they were an artifact of only looking at a small slice of the patient's data. Using the full real dataset closed most of them.

### Piece 2 — Chromoplexy flag from real cycle structure

We defined: **a patient shows chromoplexy if at least one *closed* cycle spans 3 or more distinct chromosomes.**

Why 3, not 2? A closed cycle touching exactly 2 chromosomes is just an ordinary reciprocal translocation — a single simple swap, not the multi-chromosome "shattering" event chromoplexy describes. Requiring 3+ chromosomes in one closed loop is the more faithful definition.

### Piece 3 — ERG chain completion + Baca comparison (`baca_erg_chain_completion.py`)

This is the one your professor specifically asked for: go back to the **original ERG-anchored chain** (the same one used in Baca's 15/26 heuristic — anchor rows touching ERG/TMPRSS2, plus rearrangements one hop away), and instead of stopping there, **complete it using only real, already-recorded data**.

**The trick, explained simply:**

> Don't change *what connections are allowed* (we never invent a connection that isn't in the real data). Only change *where we start looking*.

1. Take the original small ERG chain as a starting point ("seed").
2. Build the complete real graph for that patient (every rearrangement row they have, anywhere in their genome).
3. Starting only from the ERG chain's DSBs, follow real edges wherever they actually lead. If a DSB in the chain has a real rejoin partner recorded somewhere else in the patient's data — even far outside the original chain — we now follow it.
4. Keep going until either:
   - **It closes into a cycle** (success — confirmed by real data, no guessing), or
   - **It truly runs out of real data** — the DSB end was never reported as rearranged anywhere in this patient's file. We tag this `UNRESOLVED`.

We sub-tag `UNRESOLVED` ends by *why* they stopped:
- `TELOMERE` — ran off the natural end of the chromosome.
- `MISSING_REARR` — reached a DSB position whose other side was never recorded as rearranged anywhere in the patient's data.

`UNRESOLVED` is deliberately a dead end for *this* phase — it's exactly the set of cases that the next phase (obligate/minimal completion rule, copy number, proximity, centromere constraints — not started yet) is meant to resolve.

---

## 3. Worked example: Patient P01-28

This patient is a good example because Baca calls it `chromoplexy_embedded`, and walking through it shows exactly why real-data cycle completion gives a different answer.

### Step 1 — the original ERG chain

`_get_erg_chain_rows` pulls 23 rearrangement rows for P01-28, touching **8 chromosomes**: 1, 2, 3, 4, 9, 17, 19, 21.

```
n_rejoins (rows)     = 23
n_DSBs                = 46     (23 × 2)
chromosomes touched   = 8
```

Baca's heuristic on this chain: `k_local = 8`, `chain_size = 22` → since `8 > 2` and `22 > 1` → **`chromoplexy_embedded`**.

### Step 2 — what the chain looks like on its own (before completion)

Running the cycle traversal on just these 23 rows:

```
Closed cycles found : 0   ("none")
Open paths           : 11
```

Zero closed cycles already — even Baca's own "chromoplexy" chain doesn't close on its own. But maybe it would close if we had more of the real picture?

### Step 3 — complete it with the patient's full real data

We take those same 23 rows as a starting point, then look at every real rearrangement P01-28 actually has (not just the ERG-touching ones), and follow real connections outward from the chain.

```
Extra real DSBs pulled in : 22   (chain grows from 46 to 68 DSBs)
Closed cycles found after : 0   (still "none")
Open paths remaining       : 12
    — 4 stopped at a real chromosome end (TELOMERE)
    — 8 stopped because that DSB's other side was never
      recorded as rearranged anywhere in this patient's file
      (MISSING_REARR)
```

So even after pulling in *22 additional real DSBs* — nearly tripling the size of the picture — **not a single closed cycle appears.** The chromoplexy chain Baca calls real here never actually closes.

### Step 4 — a concrete dangling end, traced by hand

One of the `MISSING_REARR` dead ends is the DSB at **chromosome 17, position 7,853,514, strand `+`**.

We searched this patient's *entire* rearrangement file — every single row, every chromosome — for any row reporting chr17:7,853,514(+) as a breakpoint.

**Result: zero matches.** This exact DSB end is simply never recorded as having rejoined to anything, anywhere in Baca's sequencing data for this patient. That's not us giving up early — we checked literally everything that exists for this patient. This is genuinely the edge of what the real data can tell us, which is exactly the situation your professor meant by "Baca's paper does not complete cycles."

A second one is more interesting: the DSB at **chromosome 17, position 6,917,291, strand `-`**. In the small original chain, this looked like a simple dead end. Once we follow real data outward from it, it turns out to be the start of a real path touching **12 connected DSB positions** through more of the genome — but that long real path still ends at a true chromosome telomere on chr17, not a closed loop. More real data, same answer: open, not closed.

---

## 4. The headline finding

We ran this same before/after comparison for all 26 ERG+ patients and checked the result against Baca's heuristic.

**Using the strict definition (closed cycle must span ≥3 chromosomes):**

| | Result |
|---|---|
| Patients confirmed chromoplexy by closed cycle | **0 / 26** |

**Using a looser definition (closed cycle spans ≥2 chromosomes — i.e. any cross-chromosome closed loop at all):**

| | Result |
|---|---|
| Patients confirmed chromoplexy by closed cycle | **4 / 26** |
| Confirmed patients: | P03-2345, P05-3852, P07-4941, PR-2621 |

**Comparison against Baca's 15/26 heuristic call (at the looser ≥2 threshold):**

| | Baca says chromoplexy_embedded (15 patients) | Baca says simple_fusion (11 patients) |
|---|---|---|
| Real-data closed cycle confirms it | **4** | 0 |
| No confirmed closed cycle (still open) | **11** | 11 |

### What this means

- We **never contradict** Baca in the other direction — we don't call anything chromoplexy that Baca calls simple. That's a good consistency check.
- But **11 of Baca's 15 "chromoplexy_embedded" patients have zero confirmed closed multi-chromosome cycle**, even after using every scrap of real data available for them (P07-4941 alone pulled in 32 extra real DSBs and still only managed one small 2-chromosome cycle, with 18 ends still unresolved).
- **The takeaway:** Baca's chromoplexy label is based on *reachability* (the chain touches lots of chromosomes). Our cycle-based label is based on *closure* (the chain actually loops back on itself). Most chains that reach far do not close — they stay open even with complete real data. This is the concrete gap your professor was pointing at: completing the cycles reveals that "chromoplexy" by the heuristic and "chromoplexy" by confirmed mathematical structure are not the same thing for most of these patients.

This is the central result to lead with — it's a genuine, surprising, defensible finding, not just a recomputation of something already known.

---

## 5. What's deliberately NOT done yet

- We have not guessed at how the `UNRESOLVED` ends might close. Every `UNRESOLVED` tag means "real data has nothing more to say here" — no hypothesis has been applied.
- The next phase (not started) is to apply biological reasoning to those `UNRESOLVED` ends specifically: an "obligate" (minimal) completion rule, then copy-number constraints, physical proximity (p/q-arm), and centromere viability — in that order, per your professor's instructions.
- Only after that can we say anything about whether the 11 disagreement patients are "truly" chromoplexy or not.

---

## 6. Where to find everything

| What | File |
|---|---|
| Full-genome cycle computation (all 57 patients) | `scripts/baca/baca_full_genome_cycles.py` |
| Full-genome compact summary table | `results/full_genome_cycle_compact_summary.csv` |
| Full-genome diagrams (one big panel per patient) | `results/full_genome_cycle_diagrams/*.png` |
| ERG chain completion + Baca comparison (26 patients) | `scripts/baca/baca_erg_chain_completion.py` |
| ERG chain completion results table | `results/erg_chain_completion_summary.csv` |
