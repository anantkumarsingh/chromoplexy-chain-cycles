# Phase 3: From Chains to Patients — Comparing Our Chromoplexy Calls to Baca's

This document explains, in plain language, what Phase 3 of the Baca-AMG
project did, why it was needed, how it works step by step, and what it
found. For the full technical/algorithmic detail behind Phases 1-2 (which
Phase 3 builds on), see `CLAUDE.md` and `PHASE2_OBLIGATE_AND_PROBABLE_METHODOLOGY.md`.

Script: `scripts/baca/baca_phase3_patient_aggregation.py`
Inputs: `results/phase2_obligate_probable_completion.csv`,
`results/erg_chain_comparison_ours_vs_baca.csv`
Outputs: `results/phase3_patient_chromoplexy_summary.csv` (52 patients),
`results/phase3_baca_positive_vs_ours.csv` (the final comparison table)

---

## 1. The problem Phase 3 solves

Baca's own published headline finding is about **tumors**: "50 of 57
prostate tumors show evidence of chromoplexy." One number per patient.

Everything our project had computed up through Phase 2, though, was one
row per **chain** — and a single patient can have several chains (a
"chain" here means a cluster of rearrangements that Baca's own algorithm,
ChainFinder, already grouped together as statistically related). Our
Phase 2 table has 366 chain-rows spread across only 52 patients.

So before we can honestly say "we agree/disagree with Baca on X patients,"
we first need to collapse those 366 chain-level rows down into one row per
patient — using a rule that's fair and comparable to how Baca himself
counted. That collapsing step is Phase 3.

## 2. A 60-second recap of what each chain-level row already told us

Each of the 366 chains in the Phase 2 table has three different answers to
"does this chain show chromoplexy," because there are three different ways
of reading the data:

- **Real** — only using rearrangement connections Baca's own algorithm
  actually recorded. This is the conservative, "provable from the data"
  answer. Often, a chain only partly closes into a loop this way, leaving
  some ends dangling.
- **Obligate** — takes whatever *did* close for real, and closes the
  leftover dangling ends using the *smallest possible* invented
  connections (the "assume the least coordination" reading).
- **Probable** — same idea, but instead invents whichever connection is
  statistically the *most likely* one out of all the ways the leftover
  ends could combine.

Each of these three is then checked against two thresholds:

- **Strict**: does the resulting closed loop touch **3 or more**
  chromosomes?
- **Loose**: does it touch **2 or more**?

Strict automatically implies loose (a loop touching 3+ chromosomes touches
2+ by definition) — so these are never a "pick one" choice, both are
always reported side by side.

## 3. Step-by-step: what the Phase 3 script actually does

**Step 1 — Group chains by patient.**
Take all of a patient's rows out of the 366-row Phase 2 table (anywhere
from 1 chain up to 10 chains for one patient in this dataset).

**Step 2 — Roll each chromoplexy flag up with an "OR" rule.**
For each of the six flags above (real/obligate/probable × strict/loose),
the patient is marked positive if **at least one** of their chains is
positive. This isn't an arbitrary choice — it mirrors exactly how Baca
phrased his own result: a tumor only needs *one* qualifying chain to count
as showing chromoplexy, not all of them.

**Step 3 — Track when a "negative" isn't fully trustworthy.**
Some chains are too tangled to fully analyze (more than 12 loose ends —
past the point where every possible way of connecting them can be checked
by computer in reasonable time). If a patient's only negative evidence for
obligate or probable includes one of these un-checkable chains, we can't
honestly say "confirmed no chromoplexy" — the truth might be hiding in the
part we couldn't check. Those patients get a separate `unresolved_caveat`
flag instead of just being silently marked negative.

**Step 4 — Compute Baca's own criterion, the same way.**
Independently of anything above, check whether any of the patient's chains
has 5 or more rearrangements — Baca's own, exact, published rule for
calling a tumor chromoplexy-positive. This becomes the column
`baca_chain_has_chromoplexy`.

**Step 5 — Put them side by side.**
Build one table with the patient ID, Baca's flag, and all six of our
flags, so you can see directly, patient by patient, where we agree and
where we don't.

## 4. Worked example — why "real" and "obligate/probable" can disagree

Take patient **P07-4941**, one of their chains (chain #9, on chromosomes 7
and 17, 10 breakpoints total).

- Using only real, recorded connections, 4 of those 10 breakpoints close
  into a loop — but that loop only touches **chromosome 17**. Not
  chromoplexy by either threshold.
- The other 6 breakpoints form a loose end that never fully closes on its
  own. To force it closed, the obligate/probable methods add exactly
  **one** invented connection — and because that connection has to reach
  across to chromosome 7, the resulting (bigger) loop now touches **two**
  chromosomes, which does cross the loose threshold.

Nothing new was "discovered" here — it's the same chain, same real data.
The obligate/probable methods are simply forced to close *the entire
chain*, not just the part that already closes on its own, and closing the
leftover piece happens to cross a chromosome boundary. This is why
obligate/probable numbers are consistently higher than real numbers across
the whole dataset — it's a predictable side effect of forcing full
closure, not new evidence of coordinated breakage.

**A second example, for the "unresolved" caveat:** patient **P08-501** has
4 chains. One of them has 13 rearrangements and 16 loose ends — too many
to fully check — and it happens to be the exact chain that gives this
patient a positive Baca score (13 rearrangements clears his 5+ bar
easily). The other 3 chains are small and fully checkable; two of them
close into loops spanning 2 chromosomes, enough to confirm the loose
threshold. But nothing confirms the *strict* (3-chromosome) threshold,
and the one chain big enough to maybe do so was never checkable. So this
patient is marked `obligate_has_unresolved_caveat_strict = True` — the
"no strict chromoplexy" answer for this patient is a "we don't know," not
a "confirmed no."

## 5. Final findings

Out of the 57 Baca patients, only **52 have any usable chain data** (the
other 5 aren't covered by Baca's own chain-grouping file at all, so they
can't be scored by any of this). All numbers below are out of 52.

| Criterion | Positive patients |
|---|---|
| Baca's own rule (≥5 rearrangements in some chain) | **50 / 52** |
| Real closed cycle, strict (≥3 chromosomes) | **4 / 52** |
| Real closed cycle, loose (≥2 chromosomes) | 21 / 52 |
| Obligate closure, strict | 21 / 52 |
| Obligate closure, loose | 49 / 52 |
| Probable closure, strict | 24 / 52 |
| Probable closure, loose | 49 / 52 |

**The headline finding:** among the 50 patients Baca himself calls
chromoplexy-positive, only **4 (8%)** have a real, provable closed cycle
spanning 3+ chromosomes. Even the most generous reading (probable, loose
threshold) only recovers 48/50. In other words, Baca's own criterion
(clusters of statistically-related breakpoints) and ours (an actual closed
mathematical loop spanning multiple chromosomes) usually don't agree —
that gap is the project's central, publishable finding. It mirrors an
earlier, smaller-scale result from the ERG-specific analysis, where 11 of
15 of Baca's "chromoplexy_embedded" calls also failed to close into a real
cycle.

**Two patients — P09-396 and PR-2525 — disagree with Baca under every
lens we have** (real, obligate, and probable all say no, at both
thresholds). However, both are flagged with the `unresolved_caveat` —
each has at least one chain too large to fully check — so they should be
described as "no evidence found," not as clean, settled counter-examples
to Baca's call.

**A secondary check:** Baca also reports that 63% of his 57 tumors have 2
or more qualifying chains, not just one. In our 52-patient subset, 36
(69%) have 2 or more — a comparable, though not identical, figure (the
denominators differ: his 57 vs. our 52 mmc5-covered patients).

## 6. What's intentionally not done yet

- **The TSG-disruption version of Baca's second threshold** (his
  ≥3-rearrangement stat, 26/57 tumors, requires the chain to *also*
  disrupt a known prostate tumor-suppressor gene) was tried as a
  plain rearrangement-count cutoff, found to be uselessly non-selective
  (52/52 patients would qualify) without the gene condition, and was
  removed. Left as a follow-up task, not computed.
- **Formal statistics** (does chain size predict these disagreements,
  does ETS status or tumor grade correlate with any of the six flags,
  etc.) — a full plan has been agreed but nothing has been computed yet.
  See `CONTEXT_LOG_001.md` Session 1 for the plan.
