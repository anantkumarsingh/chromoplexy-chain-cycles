# Phase 2 Methodology: Obligate and Most-Probable Cycle Completion

**Purpose of this document:** explain, in a form suitable for presenting to others, (1) everything relevant from Cornforth 2001 — the paper that originated the term "obligate cycle structure" — and (2) the exact, separated, step-by-step algorithms Phase 2 will use to compute the obligate and most-probable completed cycle structures for each of Baca's mmc5 chains.

**Status: IMPLEMENTED (2026-06-24).** Parts 1-4 below are the original design (kept as written, with two corrections noted inline where the actual implementation deviated from the original plan). **Part 5 is new** and documents the actual script, its output file, every column's exact meaning, all caveats discovered while testing, and a real bug that was found and fixed before trusting any result.

**Script:** `scripts/baca/baca_phase2_obligate_probable_completion.py`
**Output:** `results/phase2_obligate_probable_completion.csv` (366 rows, one per (patient, Baca chain number), 47 columns)

---

## Part 1 — What We Found in Cornforth 2001

**Source:** Cornforth, M.N. "Analyzing Radiation-Induced Complex Chromosome Rearrangements by Combinatorial Painting." *Radiation Research* 155(5):643-659 (2001). This is the paper Sheth et al. 2026 cites (reference [16]) as the original source of "obligate cycle structure." It comes from radiation cytogenetics (whole-chromosome painting after X-ray/gamma-ray exposure), not cancer genomics — but the mathematical structure is the same one Sheth and our project apply to Baca's sequencing data.

### 1.1 The core problem Cornforth is solving

When radiation breaks several chromosomes near each other in time and space, the broken ends can rejoin to the *wrong* partners, producing a complex, multi-chromosome rearrangement. Painting (and, in our case, sequencing) lets you **see** the rejoining pattern, but it does not always tell you, with certainty, **which exact set of breaks acted together**. Several different real mechanisms can produce identical-looking output. Cornforth's paper is about how to analyze that ambiguity honestly rather than guessing.

### 1.2 Key terms, in the order they matter for us

| Term | Meaning | Why it matters here |
|---|---|---|
| **CAB(initial) / CAB(actual) / CAB(observed)** | Three different things: the breaks at the instant of injury, the true final rejoining, and what you actually observe. Only `observed` is ever directly seen. | Our `chrom_aberrations_baca.csv` rows are `observed`. Everything we compute about chain closure is an attempt to infer `actual`. |
| **Pattern closure** | The state where every broken end has an identified, real partner — no dangling pieces. | This is exactly our "open path" problem in the mmc5 chains, just under a different name. |
| **One-way exchange** | A broken end whose partner fragment is too small/unrecorded to detect, leaving a dangling end. | This is precisely what an open mmc5 chain end is: a real DSB with no recorded real partner. Cornforth calls this the single largest source of ambiguity in the whole field — not a rare edge case. |
| **Disintegrated rejoining outcome** | A single true coordinated multi-break event can, purely by chance, produce a pattern that looks exactly like two separate independent simple events. | Quoted directly: *"no amount of additional information, at any level of resolution, will allow a distinction... whether [they] are part of the same recombinational process."* This is why we cannot just "figure out the truth" — some ambiguity is permanent given the data type. |
| **Sequential Exchange Complex (SEC)** | A cycle structure that is *reducible* — it can be broken down into two or more smaller independent cycles (e.g. 2×C2 instead of one C4). | Distinguishes "looks simple but might be complex" (reducible) from "must be complex" (irreducible) — see the asymmetry below. |
| **Irreducible vs. reducible** | A structure is irreducible if it cannot be split into smaller independent cycles (e.g. a true C4); reducible if it can (e.g. 2×C2). | **Key asymmetry, quoted directly:** *"irreducibility... does ensure contemporaneity of breaks"* (an irreducible closed cycle in our data IS solid proof of one coordinated event) but *"reducibility is no guarantee that the breaks involved... were not proximate"* (a reducible/small structure does **not** prove independence — it's simply unprovable either way). This tells us how to phrase any chromoplexy claim: an irreducible multi-chromosome closed cycle = confirmed; an open or small/reducible result = "not provable from available data," never "disproven." |
| **Obligate cycle structure** | The smallest defensible cycle structure consistent with the data — defined precisely below. | This is the conservative, minimum-coordination-assumed reading. The direct justification for using it as our default reported number is the disintegrated-outcomes finding above: since you usually cannot prove coordination from the pattern alone, the responsible default is to assume the least. |
| **Most probable cycle structure** | Whichever cycle structure occurs most often among *all* the ways the data could combinatorially have come together. | The opposite end of the spectrum from obligate — assumes nothing about which structure is "right," just reports which one a random/uniform model of rejoining would produce most often. |
| **Homologue/homologue interaction** | A rearrangement between the *two copies* of the same chromosome (e.g. one copy of chr1 to the other copy of chr1) looks identical to a normal intrachromosomal event in painting — and in unphased short-read sequencing, like Baca's. | **Documented limitation, not fixable with current data:** some of our `Class != inter_chr` ("intrachromosomal") calls could secretly be inter-homolog events, which are mechanistically closer to interchromosomal exchanges. Should be stated as a caveat, not corrected for. |
| **Clustered/cryptic breaks** | The (rejected) idea of postulating extra, invisible breakpoints just to force a simpler-looking cycle structure. | Cornforth argues against this. **Direct rule for us: never invent a new breakpoint. Only ever connect breakpoints that already exist in the real data.** |
| **Restitution as a "cycle of order 1"** | A DSB that simply rejoins itself, producing no net change. | Minor technical note; not applicable to Baca's data since restitutions leave no trace in sequencing. |

### 1.3 The precise obligate-structure rule (this is more specific than "smallest" or "fewest cycles")

Quoted directly from the paper:

> "...the smallest (overall lowest) possible cycle structure for an actual exchange [is] the cycle structure that minimizes the largest order. Thus 2c3 is defined to be smaller than c4+c2; c4+c2 is defined to be smaller than c6. If two different cycle structures under consideration have the same largest order, one can look to minimize the next-largest order."

In plain terms: **write each candidate cycle structure as a list of cycle lengths sorted from biggest to smallest. The obligate structure is whichever list is "smallest" when you compare them entry by entry, biggest entry first.** This is *usually* the same as "the structure with the most cycles," but not strictly always — the rule is about minimizing the worst single cycle, not maximizing the count, and any implementation needs to use this exact rule, not a cycle-count shortcut.

Worked examples straight from the paper, confirming the rule:

- Among `{c4+c2, 2c3}` for a 6-break exchange: `2c3` has max part 3, `c4+c2` has max part 4 → **obligate = 2c3**.
- Among `{c5+c3, 2c4}` for an 8-break exchange: `2c4` has max part 4, `c5+c3` has max part 5 → **obligate = 2c4**.
- Among `{c6, c4+c2, 2c3}` for the paper's own worked Fig. 5 example: `2c3` (max 3) beats `c4+c2` (max 4) and `c6` (max 6) → **obligate = 2c3**.

### 1.4 The two-tiered pattern-closure rule (this directly answers a question we had previously left open)

Before you can even ask "what's the obligate cycle structure," you first have to decide **how to close** a pattern that has dangling/one-way ends — i.e., which open ends get connected to which. Cornforth gives an explicit, two-step rule for this (section 7.2 of the paper):

> Rule 1: close the pattern using the **fewest total additional/initial breaks** necessary.
> Rule 2 (used only if Rule 1 still leaves more than one option): among those, prefer the closure that involves the **fewest additional chromosomes**.

Only once closure is fixed by these two rules do you separately apply the obligate-structure rule (1.3) to whatever candidates remain tied.

**Why this matters for us specifically:** earlier in this project, we had explicitly logged as "genuinely unresolved" the question of how to decide which open ends to connect to which, when a chain has several separate open paths. Cornforth's rule above answers that question directly — it is not something we need to invent ourselves or ask the professor to define; it is established, citable methodology.

---

## Part 2 — How We Will Compute the Obligate Structure

This section describes the algorithm only for the **obligate** structure. The most-probable algorithm is described separately in Part 3 — the two are computed differently and should never be merged into a single procedure.

**Starting point for every chain:** the real rearrangement edges (fixed, from `chrom_aberrations_baca.csv` / mmc5), and the real reference edges (also fixed — adjacency, deletion-bridge, and our position-order fallback, exactly as already built in the existing closure scripts). Some breakpoints in the chain will already be fully connected (degree 2 in the combined real graph) and may already sit inside a real closed cycle. Others will be **open ends** — breakpoints with no real reference partner at all (Cornforth's "one-way" ends).

### Step-by-step

1. **Build the real graph for the chain.** Same primitives already used in `baca_mmc5_chain_closure.py`: rearrangement edges (always real, always fixed) + reference edges (adjacency / deletion bridge / fallback, also real). Identify all breakpoints with degree < 2 in this combined graph — these are the open ends needing closure.

2. **Separate true telomeric ends from genuinely ambiguous open ends.** A true chromosome-end telomere is not "missing data" — it is real biology and should never be force-closed. Only breakpoints whose other side was simply never recorded as rearranged anywhere in the data (Cornforth's actual "one-way" case) are candidates for closure. This step prevents inventing a connection where none should exist, consistent with Cornforth's rejection of cryptic/invented breakpoints (section 1.2 above).
   **DEVIATION, implemented deliberately, not an oversight (see Part 5.6):** this step was NOT applied in the actual implementation. Distinguishing a true telomere from a genuinely-missing partner requires centromere/telomere reference coordinates, which this project has explicitly deferred to a later, separate biological-constraints phase (`CENTROMERE_PQ_ARM_NOTES.md`). The shipped script treats every degree-1 breakpoint as a closure candidate. Revisit this step specifically when the centromere/telomere phase begins.

3. **Apply Rule 1 of the pattern-closure rule: minimize the number of new connections.** Since each invented edge closes exactly two open ends, the minimum number of new edges is fixed once we know how many genuinely-ambiguous open ends exist (from Step 2) — this step's real job is making sure we are not closing more ends than necessary (e.g. accidentally treating a real telomere as needing closure).
   **CLARIFICATION found during implementation:** because Step 2 (telomere exclusion) isn't applied, Rule 1 turned out to be entirely non-binding in practice — pairing up ALL open ends always requires exactly `n_dangling_ends / 2` invented edges, a number fixed before any pairing choice is made, so there is never more than one possible answer to "how many edges." Rule 1 therefore never actually discriminates between candidate closures in this implementation; only Rule 2 and the obligate rule itself do. This is recorded here so it isn't mistaken for an unused/broken step — it's mathematically always satisfied trivially.

4. **Apply Rule 2: among all ways of pairing up the open ends (perfect matchings), keep only those that minimize the largest chromosome span among the newly-formed cycles.** If a chain has more than 2 open ends, there can be multiple distinct ways to pair them up. **Implementation note, added after catching a real bug during testing:** the metric must be the *maximum* chromosome span among the new cycles, not the *sum*. An earlier version summed spans across the new cycles and was caught, by hand-tracing a real chain, to be systematically biased toward merging open paths together whenever they shared a chromosome (summing double-counts a shared chromosome, making the merged/more-coordinated option look artificially cheaper than keeping paths separate — backwards from the conservative reading this rule is supposed to produce). Using the maximum instead removes the bias: merging two paths can never produce a smaller chromosome span than the larger of the two paths alone, so merging never looks artificially cheaper.

5. **For each surviving candidate pairing from Step 4, compute the resulting full cycle structure.** Combine the chain's already-real edges with this candidate's invented closure edges, and decompose into cycles using the same connected-component/degree method already validated in the existing closure scripts (`nx.MultiGraph`, sorted node iteration — see existing bug-fix notes; both apply here unchanged).

6. **Select the obligate structure using the precise lexicographic rule (Part 1.3).** Among the candidate cycle structures produced in Step 5, write each as a cycle-length list sorted largest-to-smallest, and pick the list that is lexicographically smallest (smallest largest part; ties broken by the next-largest part, and so on). This is the chain's obligate cycle structure.

7. **Flag the result as INVENTED-closure where applicable.** If a chain required any closure at all (Steps 3-4 added edges), the obligate result must be visibly tagged as containing an invented completion — never merged into the same column/number as a chain that closed entirely from real data. Chains with zero open ends (`n_dangling_ends = 0`, 52/366 in the actual data — note this is a narrower, stricter count than the 91/366 "has at least one real closed cycle" figure used elsewhere in this project, since a chain can have one real closed cycle AND separate open ends at the same time) simply report their real structure unchanged, with `uses_invented_edges = False`, for direct comparison.

8. **Repeat for all 366 chains**, not only the 275 currently open ones — per the decision to also see what the obligate framing produces for the chains that already close on their own.

---

## Part 3 — How We Will Compute the Most-Probable Structure

This is a **separate procedure** from Part 2. Obligate asks "what's the most defensible minimum"; most-probable asks "if every valid way the data could have come together were equally likely, which outcome would show up most often." They are computed independently and reported as two separate columns — never collapsed into one "the" answer (this mirrors exactly how Sheth's own worked example in section 7.1 reports both side by side).

### Step-by-step

1. **Start from the same real graph and the same set of open ends as Part 2, Steps 1–2.** Real edges stay fixed; only the open ends are candidates for re-pairing.

2. **Enumerate every valid way to pair up the open ends — not just the minimal ones.** This is the same style of exhaustive enumeration already implemented for small ERG chains (`find_one_example_per_structure`, Tables A/B/C): every possible perfect matching of the open-end set is a separate "pathway." Unlike Part 2, we do not discard pairings that touch extra chromosomes or use more connections than the bare minimum — every topologically valid pairing counts.

3. **For each enumerated pairing, compute the resulting full cycle structure**, using the identical combine-and-decompose method as Part 2 Step 5 (same real edges + this pairing's invented edges → cycle decomposition).

4. **Tag each resulting cycle structure and count how many of the enumerated pairings produce it.** This produces a frequency table per chain: structure → count → percentage of all valid pairings — exactly the Table A/B/C format already used for the small ERG chains.

5. **The most-probable structure is whichever has the highest count/percentage** in that table — the mode of the distribution, with no judgment applied about plausibility, just raw frequency.

6. **If two or more structures tie for the highest frequency, flag this explicitly** rather than arbitrarily picking one — consistent with how every other ambiguity in this project has been handled so far (e.g. the Phase 1.5 ambiguous-node findings).

7. **Apply a tractability cap, consistent with the existing precedent.** The small-ERG-chain enumeration was already capped at `n_rejoins ≤ 6` (≤10,395 matchings) for computational reasons. The number of open ends in a chain plays the same role here.
   **RESOLVED during implementation, better than expected:** the cap is on `n_dangling_ends` (open ends needing closure), not on chain size — and most large chains are already mostly closed by real data, leaving comparatively few genuinely open ends. Verified empirically before finalizing the cap: 354/366 chains have `n_dangling_ends ≤ 12` (the exact same `(m-1)!! ≤ 10,395` threshold already used elsewhere), only 12 chains exceed it (up to 24 open ends). Those 12 are flagged `phase2_enumerable = False` with every obligate/probable column left null — not approximated, not sampled. Full list and exact counts in Part 5.4.

8. **Repeat for all 366 chains**, same as Part 2.

---

## Part 4 — Why Both, and How to Report Them

- **Obligate** is the conservative, defensible-minimum number: the answer to "what is the least we are forced to conclude, given the real data plus the smallest necessary closure." It is what we should lead with when making any chromoplexy claim, because it cannot be accused of overclaiming.
- **Most-probable** is the combinatorial-likelihood number: the answer to "if the rejoining process were a random matching among everything still ambiguous, what would typically come out." It is informative as a contrast (it usually points toward a much larger, more coordinated structure than obligate does) but is **not** itself proof of anything — it describes a probability distribution over possibilities, not a measurement.
- **Per the irreducible/reducible asymmetry (Part 1.2):** when either column produces an *irreducible* multi-chromosome closed cycle, that result can be stated as confirmed chromoplexy. When a chain stays open, or its obligate/most-probable structure reduces to small independent cycles, the correct statement is "not provable from available data" — never "proven independent" or "not chromoplexy."
- **Swap-stability (Phase 4)** is a separate, later analysis on top of whatever structures come out of Parts 2 and 3 — not part of this methodology, and not a way to choose between obligate and most-probable.

---

## Part 5 — Implementation: Output File, Column Reference, Caveats, and the Bug That Was Caught

### 5.1 A real bug, found and fixed before trusting any output

The first version of Part 2 Step 4's Rule-2 metric summed chromosome spans across the newly-formed cycles. This was **wrong**: summing is biased TOWARD merging open paths together whenever they share a chromosome, because `span(A) + span(B) > span(A∪B)` whenever A and B share a chromosome (summing double-counts the shared chromosome). This made the merged/more-coordinated closure look artificially cheaper than keeping paths separate — exactly backwards from the conservative reading Cornforth's rule is supposed to produce.

**How it was caught:** before trusting any aggregate result, obligate and probable were compared chain-by-chain. The buggy version produced `obligate == probable` for **all 302** chains that needed any invented closure at all — a statistically implausible result (hand-tracing a real chain with only 2 open paths showed the *un*-merged, conservative option should usually be at least a candidate, not eliminated every single time). Hand-tracing `P01-28` chain 1 confirmed the sum-based metric was systematically favoring the merged option.

**Fix:** use the **maximum** chromosome span among the new cycles, not the sum. This has no merge bias, because `span(A∪B) ≥ max(span(A), span(B))` always holds — merging two paths can never produce a smaller span than the larger of the two paths alone.

**After the fix:** obligate and probable differ for 194/302 chains requiring invention (a believable, non-degenerate result), and obligate is consistently the more conservative of the two cohort-wide (`obligate_has_chromoplexy_strict` = 35/354 enumerable chains vs `probable_has_chromoplexy_strict` = 39/354) — exactly the expected relationship between a conservative and a likelihood-based metric.

### 5.2 Cross-checks performed (not just claimed)

- Every chain's from-scratch real-graph decomposition (built independently inside the Phase 2 script) is asserted, in code, to exactly match the already-published `compute_chain_cycle_structure` result from `baca_mmc5_primary_chain_cycle_structures.py` — **0 mismatches across all 366 chains.**
- `n_matchings_enumerated` exactly follows the `(m-1)!!` formula at every observed `m`: 1, 1, 3, 15, 105, 945, 10395 for m = 0, 2, 4, 6, 8, 10, 12 — verified directly against the data, not just assumed from the formula.
- `m` (number of dangling open ends) is always even, and no breakpoint ever has degree 0 in the real combined graph — both verified empirically across all 366 chains before writing any enumeration code (a graph with max degree 2 can only decompose into simple paths and simple cycles, so every open component has exactly 2 endpoints — this was asserted in code, not just assumed).
- For all 52 chains with `n_dangling_ends = 0`, `obligate_cycle_structure` and `probable_cycle_structure` both exactly equal `real_cycle_structure` — **0 mismatches.**
- `n_invented_edges == n_dangling_ends / 2` and `uses_invented_edges == (n_dangling_ends > 0)` hold for all 366 rows.

### 5.3 Deliberate scope decisions (documented, not oversights)

- **No telomere/centromere distinction.** Every degree-1 breakpoint is treated as a closure candidate; a true chromosome terminus is not distinguished from a genuinely-missing real partner, because that requires centromere/telomere reference coordinates this project has explicitly deferred to a later phase (see Part 2 Step 2 above, and `CENTROMERE_PQ_ARM_NOTES.md`).
- **Never invent a new breakpoint.** Only ever connects two breakpoints that already exist in the real mmc5 data — consistent with Cornforth's rejection of cryptic/invented breaks (Part 1.2).
- **Rule 2 is OUR adaptation, not Cornforth's literal rule.** Cornforth's original "fewest additional chromosomes" assumed an unknown missing fragment could be assigned to any chromosome; here, every open end is an already-known, already-located real breakpoint, so the literal rule never discriminates (the union of chromosomes touched by all open ends is fixed regardless of pairing choice — verified). The adapted metric (max chromosome span among new cycles, see 5.1) is explicitly our own extension of the same conservative spirit, not a transcription of his rule.

### 5.4 The 12 chains that exceed the enumeration cap (`phase2_enumerable = False`)

| Patient | Chain | n_dangling_ends |
|---|---|---|
| P01-28 | 4 | 14 |
| P07-4941 | 4 | 14 |
| P07-4941 | 5 | 24 |
| P08-492 | 3 | 20 |
| P08-501 | 1 | 16 |
| P09-396 | 6 | 14 |
| PR-02-1431 | 5 | 16 |
| PR-02-1431 | 13 | 16 |
| PR-06-1749 | 3 | 22 |
| PR-07-360 | 5 | 22 |
| PR-07-4814 | 3 | 16 |
| PR-2525 | 5 | 14 |

For these 12 rows, every `obligate_*` and `probable_*` column is left blank (null) rather than approximated or sampled. `n_matchings_enumerated` is also null for these rows specifically (not zero — enumeration was never attempted).

### 5.5 Two things that look like missing data in the CSV but aren't

- **`probable_tied_alternatives` is blank for every single row** in the current data — not because anything failed, but because no chain ever produced a genuine frequency tie between two structures (`probable_is_tied = False` everywhere). The column exists so a future re-run that does find a tie won't have it silently hidden.
- **`real_chromosomes_in_cycles` is blank for 275 rows** — those are exactly the chains with `real_cycle_structure = "none"` (no real closed cycle at all), so there is no chromosome list to report. Both of these render as blank/NaN when the CSV is reloaded because they were written as empty strings, not because anything is missing computationally.

### 5.6 Full column reference

**Identifiers & clinical context**

| Column | Meaning |
|---|---|
| `patient_id` | Baca patient ID |
| `baca_chain_number` | Which of that patient's mmc5 chains this row describes (Baca's own `Chain number`) |
| `n_breakpoints` / `n_rearrangements` | Size of the chain — total breakpoints, and total rearrangement events (each rearrangement = 2 breakpoints) |
| `k_chromosomes_whole_chain` / `chromosomes_whole_chain` | How many, and which, chromosomes the chain touches at all (regardless of cycle closure) |
| `contains_ERG_or_TMPRSS2` | Whether this chain includes an ERG/TMPRSS2 site annotation |
| `genes_in_chain` | All genes annotated at any breakpoint in the chain |
| `ETS_status`, `Gleason_Score`, `Pathological_stage` | Patient-level clinical fields, joined in from `clinical_phenotypes.csv` (Gleason_Score null for 15 patients lacking a recorded score — real missing clinical data, unrelated to Phase 2) |

**`real_*` — the pre-Phase-2 baseline (Baca's real data only, no invented edges; identical to `results/mmc5_primary_chain_cycle_structures.csv`, cross-checked to match exactly)**

| Column | Meaning |
|---|---|
| `real_has_any_closed_cycle` | Does the chain already contain ≥1 real closed cycle, using only real rearrangement + adjacency + deletion-bridge edges? |
| `real_chain_fully_closed` | Are *all* breakpoints in the chain part of a real closed cycle (no open ends anywhere)? Equivalent to `n_dangling_ends == 0`. |
| `real_n_nodes_in_open_paths` | Total breakpoints sitting inside open (unclosed) paths (NOT the same as `n_dangling_ends` — this counts every node in an open path, `n_dangling_ends` counts only the 2 loose ends of each open path) |
| `real_pct_breakpoints_in_closed_cycles` | % of the chain's breakpoints already in a closed cycle |
| `real_n_cycles` | Number of real closed cycles found |
| `real_cycle_structure` | The real cycle structure, e.g. `"C6x1 + C2x2"` (one 6-cycle, two 2-cycles); `"none"` if nothing closes |
| `real_theta_cycles` | Sheth's Θ(k,(b₁,...)) notation, scoped only to breakpoints inside real closed cycles |
| `real_chromosomes_in_cycles` | Which chromosomes are inside a real closed cycle (blank if none — see 5.5) |
| `real_max_cycle_chrom_span` | The largest number of distinct chromosomes spanned by any single real closed cycle |
| `real_has_chromoplexy_strict` / `real_has_chromoplexy_loose` | Chromoplexy flag from real data alone: strict = some closed cycle spans ≥3 chromosomes; loose = ≥2 |

**Scope / enumeration metadata — describes the Phase 2 problem itself, before any selection**

| Column | Meaning |
|---|---|
| `n_dangling_ends` | How many breakpoints are "open" (missing a real partner) and need an invented connection — this is `m` |
| `n_open_path_components` | How many separate open paths exist (`= n_dangling_ends / 2`, since every open path has exactly 2 loose ends) |
| `uses_invented_edges` | Whether this chain needed *any* invented connection at all (`n_dangling_ends > 0`) |
| `n_invented_edges` | How many invented connections were used to fully close the chain (`= n_dangling_ends / 2`) |
| `phase2_enumerable` | Whether the chain was small enough to enumerate every possible closure (`n_dangling_ends ≤ 12`) |
| `n_matchings_enumerated` | How many distinct ways of closing the chain were actually tried — follows `(m-1)!!` exactly; null if not enumerable |

**`obligate_*` — the conservative, "least coordination assumed" completion (Cornforth's rule, Part 1.3/2)**

| Column | Meaning |
|---|---|
| `obligate_cycle_structure` | The chosen conservative completed structure — among all closures tied at the minimum chromosome-span penalty (Rule 2), whichever has the smallest largest-cycle (Cornforth's exact rule) |
| `obligate_n_cycles` | Number of cycles in that structure |
| `obligate_max_cycle_chrom_span` | Largest chromosome span among the obligate structure's cycles |
| `obligate_max_span_is_unique` | Whether every matching that achieved this exact structure had the same chromosome span. False means two different specific closures realize the same cycle-LENGTH structure but with different chromosome spans — the minimum span is reported in that case, never silently averaged. A console `[note]` is printed during the run whenever this happens (none occurred in the actual run). |
| `obligate_has_chromoplexy_strict` / `_loose` | Same chromoplexy thresholds (≥3 / ≥2 chromosomes), applied to the obligate structure |
| `obligate_rule2_metric` | The minimized value itself: the smallest achievable "largest chromosome span among new cycles" |
| `obligate_n_matchings_at_rule2_minimum` | How many of the enumerated matchings tied at that minimum (the candidate pool obligate was chosen from) |
| `obligate_n_matchings_achieving_structure` | Of those, how many actually produced the exact reported structure (the structure itself is always unique by construction — see Part 1.3's note on why tuple comparison is a strict total order here — this column just shows how many different specific closures realize it) |

**`probable_*` — whichever completion is combinatorially most common, no conservatism assumed (Part 1.3/3)**

| Column | Meaning |
|---|---|
| `probable_cycle_structure` | The structure that occurred most often across every enumerated way of closing the chain |
| `probable_n_cycles` | Number of cycles in that structure |
| `probable_max_cycle_chrom_span` | Largest chromosome span among its cycles (same minimum-among-ties handling as obligate) |
| `probable_max_span_is_unique` | Same uniqueness flag as above, for the probable structure's achieving matchings |
| `probable_has_chromoplexy_strict` / `_loose` | Chromoplexy thresholds applied to the probable structure |
| `probable_pct` | What % of all enumerated matchings produced this exact structure |
| `probable_n_distinct_structures` | How many genuinely different structures appeared at all (a measure of how combinatorially ambiguous the chain is) |
| `probable_is_tied` | Whether two or more structures tied for most-frequent (none did, in the actual run — see 5.5) |
| `probable_tied_alternatives` | The other structures tied at that same frequency, if any (always blank in the actual run — see 5.5) |

### 5.7 Aggregate results from the actual run (cohort-wide, 366 chains / 52 patients)

```
Chains with n_dangling_ends = 0 (already fully closed, no invention needed): 52 / 366
Chains requiring invented closure (n_dangling_ends > 0):                   314 / 366
Chains enumerable (n_dangling_ends <= 12):                                 354 / 366
Chains NOT enumerable (> 12, see 5.4 for the list):                         12 / 366

Of 302 enumerable chains requiring invention: obligate == probable for 108, differ for 194
Chains where probable_is_tied: 0 / 302

real_has_chromoplexy_strict      : 4 / 366
obligate_has_chromoplexy_strict  : 35 / 354 (enumerable only)
probable_has_chromoplexy_strict  : 39 / 354 (enumerable only)
```

This is the first concrete, quantitative answer to "what extra information do completed cycles provide" (the professor's Phase 3 question) — even before the formal Phase 3 reclassification write-up, the raw flags already show completion roughly 9x's the strict-chromoplexy count under the conservative (obligate) reading, and the probable reading goes slightly further still.

---

## References

- Cornforth, M.N. "Analyzing Radiation-Induced Complex Chromosome Rearrangements by Combinatorial Painting." *Radiation Research* 155(5):643-659 (2001).
- Sheth, S., Arsuaga, J., Sazdanovic, R. (2026). *J. Phys. A: Math. Theor.* 59:115601.
- Baca, S.C. et al. (2013). "Punctuated Evolution of Prostate Cancer Genomes." *Cell* 153:666-677.
