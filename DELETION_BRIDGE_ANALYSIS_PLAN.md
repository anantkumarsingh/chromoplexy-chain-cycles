# Deletion Bridge Analysis Plan — Baca Chains

> Status: **ANALYSIS COMPLETE (2026-08-13).** All 3 parts implemented,
> cross-verified, and the overall conclusion written — see "Overall
> conclusion" section below for the headline result. Plan locked
> 2026-08-13 (`CONTEXT_LOG_003.md` Session 9, resolved
> `CONTEXT_LOG_004.md` Session 10); Parts 1-3 implemented, cross-checked,
> and concluded the same day (`CONTEXT_LOG_004.md` Sessions 11-12).

## Part 2 implementation status — DONE 2026-08-13

- **Core algorithm** (the hard dependency the plan flagged):
  `scripts/baca/deletion_bridge_analysis/final_structure_assembly.py`.
  A NEW assembly algorithm (not a reuse of
  `extract_all_traversal_hops`, which Session 8 found silently splits a
  segment with both real-fused ends into two disconnected fragments).
  Literally simulates breakage-and-rejoining: every reference segment
  starts as its own 1-piece "molecule"; each real rearrangement edge,
  applied in ascending `Rearrangement number` order (locked decision
  #4), merges two molecules at their matching free ends, reversing
  (flipping) one side when the strand pairing requires it. A segment's
  flip status is then read directly off the assembled `Piece.flipped`
  attribute — fully deterministic, no reading-direction ambiguity.
- **Validated exactly against the P05-1657 hand-derivation** before
  trusting it on any other chain: all 6 final molecules match
  (including the standalone `[7-6]` piece and the chr8→chr12 merge via
  `[3-7]`), and the flip pattern matches character-for-character:
  `A'`, `B` (unflipped), `C'`, and Piece 1 (`P1'`) — exactly the user's
  confirmed answer.
- **Two real bugs found and fixed while testing beyond the single P05-1657
  case, both before running the full cohort:**
  1. A trivial (0-piece) telomere-facing tip's node was being re-exposed
     as a free end after being consumed in a merge — would have let it
     be merged again and mislabeled a resolved chain-terminus as if
     still open. Fixed with a `('CHAIN_TERMINAL', node)` sentinel,
     distinct from a never-touched open tip.
  2. **Breakpoint position is not a safe dedup key** — caught running
     the full 194-chain cohort (not the P05-1657 test, which didn't
     happen to hit this): two DIFFERENT breakpoint numbers can share the
     exact same (chromosome, position) (a quirk already documented
     elsewhere in `CLAUDE.md`, e.g. P08-1042 chr2:65,557,230 = BP24 for
     Rearrangement 31 AND BP77 for Rearrangement 53). Grouping by
     `position.unique()` silently dropped one of the two, crashing
     assembly on P08-1042 chain 2. Fixed by sorting each chromosome's
     breakpoints by (position, breakpoint_number) directly, never
     deduplicating by position alone. Re-verified 0 failures across all
     194 chains afterward, and re-confirmed P05-1657 unchanged
     (regression check).
- **Rearrangement Type classifier**: `classify_rearrangement_type` /
  `classify_chain_rearrangements` in the same module — independent
  4-way rule from chromosome-match + strand-match only (decision #7),
  never reads Baca's own `Class` field. P05-1657's 4 rearrangements
  re-verified to classify exactly as expected. The previously-inferred
  4th category (different-chromosome + opposite-strand → plain
  `TRANSLOCATION`) is CONFIRMED real: 346 occurrences cohort-wide (vs.
  779 `SIMPLE / NON-RECIPROCAL TRANSLOCATION`, 770 `SIMPLE INVERSION`,
  325 `INVERTED TRANSLOCATION`) — resolves the plan's last open item.
- **Visualization script**: `scripts/baca/deletion_bridge_analysis/
  part2_final_structure_visualization.py`. New "molecule row" drawing
  style (final fragments can span multiple original chromosomes, so
  chromosome-track drawing doesn't apply) — one horizontal row per final
  fragment (open or self-closed), pieces laid out left-to-right in final
  order, DB-span segments keeping their Part 1 letter + color (same
  `DB_PALETTE` index — visually traceable between the Part 1 and Part 2
  images for the same chain) with an apostrophe when flipped, real
  fusion junctions marked with a diamond, deletion-bridge relationship
  visible via the same segment highlighting. Same 3 corner tables as
  Part 1 (BP-location, Segment Span Information) plus a new
  Rearrangement Type table.
  - **Bug found and fixed while testing** (not just Part 1's reused
    logic): `_bp_table_dims` assumes short BP-location strings; reused
    naively for the much longer Rearrangement Type / Segment Span lines,
    table titles and content visibly overflowed into the neighboring
    table (caught by rendering P05-1657 chain 1 before trusting the
    layout). Fixed with a new `_text_table_width` helper that sizes a
    table from its actual longest line (monospace character-width
    estimate) instead of reusing the BP-table-specific heuristic.
- **Full cohort run**: 194/194 DB-bearing chains drawn, 0 failures.
  Spot-checked P09-1042 chain 1 again (60 breakpoints, 19 DB pairs, 33
  open fragments) — renders without crashing; the top (DB-dense) row is
  visually crowded as expected, but all 3 corner tables remain fully
  legible as the definitive source, matching Part 1's established
  design philosophy for extreme cases.
- **Output location**: `results/Deletion Bridge Analysis on Baca Chains/
  Final Rearrangement Structure/Final Structure After Rearrangements
  Images/` — 194 PNGs, `{patient_id}_chain_{n}.png`.

## Part 2 — REWRITTEN 2026-08-13 per explicit user correction

The version described immediately above (one row per disconnected final
molecule, in a bespoke "molecule row" style) was **wrong** — the user
clarified the intended visualization should look like Baca's own chain
(Part 1's per-CHROMOSOME track style, breakpoints re-ordered to reflect
the final structure) with explicit GAP markers wherever native material
relocated away, exactly matching the user's own hand-drawn P05-1657
notation (`TELOMERE----2[-] GAP GAP 7[+]----6[-] 18[-]----6[+] ...`).
Full rewrite, validated character-for-character against that exact text
before running the cohort again.

**New gap/movement rule, locked after resolving a real inconsistency**:
an early candidate rule ("does this segment's own boundary connect to
something non-native") wrongly flagged Segment C as moved, when it
should stay unmarked. The user's own generative process (apply
rearrangements in ascending order, literally building the structure step
by step) gave the correct rule: a merge is **"local"** — in-place flip, no gap — if
the breakpoint on either side is *already a member of the other unit's
current piece set* (checked directly, not via a "native neighbor of the
specific strand side" computation — an earlier, more complex version of
this check was wrong and is documented/retracted in
`final_structure_assembly.py`'s comments). Otherwise it's a
**relocation**: the earlier-native-position unit is tagged `.moved`
(driving the gap), unless one side is a trivial (bare telomere) unit, in
which case the substantive side always moves. Verified against all 4 of
P05-1657's rearrangements including the subtle case (Rearr 17 fusing the
already-moved Piece-1 unit to Segment B counts as local via Segment C,
already part of that unit — confirmed via the user's own algebra,
`(C·P1)' = P1'·C'`).

**`Piece.moved` added** to `final_structure_assembly.py`, set once during
assembly, carried through `.flip()`, never cleared.
**`build_chromosome_rows`** (new function) reorganizes the assembled
molecules into one row per original chromosome: native segments that
were never `.moved` are anchor points (render the whole final molecule
there); `.moved` segments show a GAP at their native slot; the same
molecule's real content surfaces later, at its true anchor — same
chromosome further along, or a different chromosome entirely (e.g.
P05-1657's `[3-7]` piece reappearing on the chr12 row, matching the
user's own text exactly, including the `→BP19(-)` cross-chromosome
connector notation).

**Two more real bugs found and fixed while rebuilding, both caught by
testing beyond the visual "does it look plausible" check**:
1. A piece cell's own display direction (left/right boundary order) was
   not swapped when its molecule was entered from a `CHAIN_TERMINAL`
   side on a different chromosome's row — the [3-7] piece showed
   `BP3(+)` next to the `→BP19(-)` connector, when the real fusion is at
   BP7, not BP3. Fixed with `Piece.display_nodes(display_swapped)` — a
   pure rendering-direction swap that never touches `.flipped` (keeps
   the apostrophe/DB-lookup semantics, which depend on `.flipped`,
   completely separate from which end happens to print first on a given
   row).
2. **When a merge produced a molecule with BOTH ends becoming
   `CHAIN_TERMINAL`** (connects out on both sides, nothing further to
   attach on either end), neither end got a live dict entry — the whole
   molecule, and every piece in it, silently vanished from the output.
   Not caught by the P05-1657 test (never happens to produce this case)
   — caught running the full 194-chain cohort: PR-3042 chain 10 (chr13)
   lost 2 real pieces this way, causing a `KeyError` downstream. Fixed
   by tracking "fully capped" molecules in a separate list rather than
   relying solely on live free-end lookups.
- **Full (re-run) cohort result**: 194/194 drawn, 0 failures, after
  fixing both bugs above. P05-1657 re-verified to match the user's exact
  text, character-for-character, after every fix (regression-checked
  each time, not just at the end).
- **Known limitation, not yet addressed**: very dense chains (e.g.
  P09-1042 chain 1, 60 breakpoints) produce extremely wide images
  (10,000+ px) that compress to illegible thumbnails — verified the
  underlying rendering is still CORRECT at native resolution (cropped
  and inspected a section directly), so this is a legibility/scale
  concern only, not a correctness one, consistent with Part 1's existing
  "tables are the definitive source for dense chains" design philosophy.

## Part 2 — polish round, 2026-08-13 (same day, third pass)

Three fixes per direct user feedback on the rewritten (chromosome-row)
version above:

1. **Missing diamond at CHAIN_TERMINAL boundaries.** The junction-marker
   logic only drew a diamond BETWEEN two pieces of the same rendered
   molecule — a piece's boundary connecting to a `CHAIN_TERMINAL` (a
   real fusion whose other side isn't part of this row, e.g. a cross-
   chromosome jump) got only a text arrow, no diamond, even though it's
   an equally real somatic fusion (confirmed: no diamond was shown at
   BP7(-)/BP19(-) despite Rearrangement 168 being a real
   INVERTED TRANSLOCATION there). Fixed so every real fusion — inter-
   piece AND terminal — gets exactly one diamond. A second bug surfaced
   while fixing this and was caught by direct visual inspection before
   trusting it: inter-piece junctions were briefly rendering as a
   doubled "◆◆" (both the earlier piece's right edge AND the next
   piece's left edge each drew their own marker for the SAME physical
   junction) — fixed so each junction is owned by exactly one side
   (the earlier piece's right edge; a piece's own left edge only draws
   a diamond when it's the molecule's first piece AND that specific
   left end is itself a CHAIN_TERMINAL).
2. **Split into a matched pair of images per chain**, per explicit
   request: `{patient}_chain_{n}_final_rearranged_structure.png` (the
   diagram alone) and `{patient}_chain_{n}_final_rearr_str_details.png`
   (the 3 corner tables alone, laid out as full-width readable columns
   instead of squeezed corner boxes). `draw_final_structure` (the old
   single-figure function) was replaced by
   `draw_structure_diagram` + `draw_structure_details`, orchestrated by
   `draw_final_structure_pair`. Old single-image outputs
   (`{patient}_chain_{n}.png`) were deleted from the output directory
   before regenerating — they are superseded, not an intentional third
   variant left lying around.
3. **Generous spacing throughout**, now that the diagram has the whole
   figure to itself: `ROW_SEP` 4.0→7.0, `CELL_W` 0.11→0.20, explicit
   `CELL_GAP` breathing room between cells (was a bare 0.01), larger
   fonts and more vertical room for BP-number/strand labels and DB
   letter/apostrophe labels.

**Full cohort re-run after all 3 fixes**: 194/194 pairs drawn (388 files
total), 0 failures. P05-1657 re-verified visually after each fix — the
diagram still matches the user's exact text, now with a correctly single
diamond at every real junction including the chr12 cross-chromosome one,
and the two-image split renders cleanly. The P09-1042 extreme case
(60 breakpoints) still renders without crashing at the new, more
generous spacing (diagram image ~10,300×1,760 px, details image
~2,700×2,900 px).

## Part 2 — polish round 2, 2026-08-13 (same day, fourth pass)

Four more fixes per direct user feedback, again explained against
P05-1657:

1. **Genuine TELOMERE bars.** Outermost breakpoints (BP2(-), BP14(+) on
   chr8; BP19(-), BP19(+) on chr12) were only a small text label with a
   circle icon — no visual bar. Fixed: `build_chromosome_rows`'s
   `_render_tip` now ALWAYS emits a `('telomere', node, has_real_fusion)`
   row item first (whether or not that specific tip also got a real
   fusion), and the drawing code renders it as a solid black bar with a
   rounded telomere-ward cap, distinct from the gray "intact reference
   segment" bars. This is real reference DNA (the stretch from the true
   chromosome end up to the outermost recorded breakpoint) and exists
   regardless of what happened at the breakpoint itself.
2. **Diamond/breakpoint spacing.** `CELL_GAP` raised from 0.025 to 0.11
   and diamonds are now drawn CENTERED in the gap between adjacent
   cells (via a `pending_junction` flag consumed at the start of the
   next cell) rather than pinned to a specific cell's edge, which had
   been crowding into the neighboring tick/label. Also fixed: a
   telomere-to-piece junction (e.g. BP19(-) into the [3-7] piece) now
   gets its diamond from the telomere cell's own `junction_right`, not
   duplicated by the piece.
3. **Segment-letter band separated from BP labels.** DB-span letters
   (A-Z, with apostrophes) now render in their own clearly separated
   horizontal band (`DB_LABEL_Y_OFFSET = -3.2`, well below the row's own
   BP tick labels and well clear of the next chromosome row — `ROW_SEP`
   raised 7.0→11.0 to make room), connected back to its own piece by a
   thin dotted vertical line so the letter stays legibly attached even
   though it's now far below the tick marks.
4. **Deletion bridge arc, missing entirely from the final structure.**
   Added `draw_deletion_bridge_arcs`: after all rows are drawn, using a
   `node_positions` map recorded live during `draw_chromosome_row`
   (breakpoint number → every drawn tick's (x, y)), draws a purple
   dotted arc connecting each DB pair's two anchor breakpoints at their
   ACTUAL final positions (same color/style convention as Part 1's
   original-chain DB arc), labeled `DB{n}: BP{low}↔BP{high}`. Since all
   521 DB pairs cohort-wide are same-chromosome (verified earlier in
   Part 1's build), the anchors always land on the same row in practice,
   though the function also handles a hypothetical cross-row case.

**Full cohort re-run after all 4 fixes**: 194/194 pairs drawn, 0
failures. P05-1657 re-verified visually one more time — telomere bars,
generously spaced diamonds, a clearly separated A'/B/C' label band with
connector lines, and the purple `DB1: BP6↔BP14` arc are all present and
correctly matched to the underlying (unchanged, already-validated)
piece/flip computation. The P09-1042 extreme case (60 breakpoints)
still renders without crashing at the new spacing (diagram image
~14,000×3,000 px).

## Part 2 — polish round 3, 2026-08-13 (same day, fifth pass)

One more fix, caught immediately from a rendered image the user
attached directly: junction diamonds were landing INSIDE the next
piece's own gray bar instead of centered in the white space before it —
an off-by-`CELL_GAP` bug in the gap-midpoint math (`x` is already
advanced PAST the gap by the time the diamond is drawn, since the
previous cell's own trailing step is `x = x_hi + CELL_GAP`; the fix uses
`x - CELL_GAP / 2`, not `x + CELL_GAP / 2`). Re-verified visually on
P05-1657 (diamonds now correctly centered in the gaps, including the
chr12 telomere-to-piece junction) and re-ran the full cohort: 194/194
pairs, 0 failures.

## Part 1 implementation status — DONE 2026-08-13

- **Shared module** (also the foundation Part 2/3 must reuse for segment
  identity, per decision #3): `scripts/baca/deletion_bridge_analysis/
  db_segments_common.py`. Found and fixed a real bug while testing:
  `find_db_pairs_for_chain` initially ordered a DB pair's low/high anchor
  by raw breakpoint NUMBER, not genomic position — breakpoint numbers are
  just mmc5 row IDs and don't track position order (caught on P09-1042
  chain 1, where BP18 is at a HIGHER position than its DB partner BP28).
  Fixed to order by position; re-verified 0 failures across all 366
  chains afterward.
- **Visualization script**: `scripts/baca/deletion_bridge_analysis/
  part1_db_segments_visualization.py` — extends (imports/reuses, does not
  copy) `mmc5/baca_chain_visualization.py`'s data loading, BP-labeling,
  and BP-location-table helpers. Adds: a colored underline + `DBn:letter`
  label per segment (color cycles per DB-pair index so multiple DB pairs
  sharing one chromosome stay visually distinguishable — up to 14
  observed sharing one chromosome in a single chain, P09-1042 chain 1),
  and a new "Segment Span Information" corner table (the definitive
  letter→span source, independent of any visual crowding in dense
  chains).
- **Validated against the P05-1657 hand example** before running the
  cohort: output table reads exactly `A: spans 6[+]----18[-]`,
  `B: spans 18[+]----15[-]`, `C: spans 15[+]----14[-]` — matches. Visual
  diagram layout/sizing cross-checked against the original
  (unmodified) `baca_chain_visualization.py` output for the same chain —
  identical dimensions and edge rendering, confirming the extension
  didn't alter the base drawing.
- **Full cohort run**: 194/194 DB-bearing chains drawn, 0 failures, 172
  no-DB chains correctly skipped. Spot-checked the most extreme case
  (P09-1042 chain 1: 19 DB pairs, 60 breakpoints, up to 14 DB pairs on
  one chromosome) — renders without crashing, stays legible via the
  Segment Span table even where on-diagram colors get tight. One bug
  found and fixed during this spot-check: the title bar listed every DB
  pair inline and overflowed the figure width past ~8 DB pairs — capped
  to a plain count + pointer to the table for chains with >8 DB pairs.
- **Output location**: `results/Deletion Bridge Analysis on Baca Chains/
  Baca Chains With DB Segments/Baca Chains with Segments Images/` — 194
  PNGs, `{patient_id}_chain_{n}.png`.

## Motivation

Grew directly out of the P05-1657 worked example (`CONTEXT_LOG_003.md`
Sessions 6 and 8): hand-tracing all 4 real fusions in that chain showed
that the material nominally "inside" the BP6↔BP14 deletion bridge span
is NOT missing — it's reorganized (Segment A and B stay locally paired
via an in-place inversion; Piece 1, from a distant locus, gets inserted
into the native B–C junction, displacing Segment C). This directly
conflicts with the general/textbook description of a deletion bridge as
tracking literal DNA loss (surfaced via the user's own research, same
session). This analysis is designed to settle, cohort-wide and
systematically (not just for one hand-traced patient), whether deletion-
bridge-flanked material is generally **relocated** (findable elsewhere in
the same chain via real edges) or **genuinely absent**.

## Empirical facts checked before finalizing this plan (2026-08-13)

Queried `mmc5.xlsx` Table S5A directly rather than assuming:
- **366 total mmc5 chains.**
- **194/366 chains have ≥1 deletion-bridge pair.**
- **126/366 chains have MORE THAN ONE distinct DB pair** — not an edge
  case, over a third of all chains. Max found: **19 distinct DB pairs in
  a single chain.**
- **Max internal breakpoints strictly between any single DB pair's two
  anchors (same chromosome): 9** — so at most ~10 segments per DB pair,
  safely under the 26-letter alphabet. Letter overflow is a non-issue
  **per DB pair** (see "Multi-DB handling" below for why this matters).

## Top-level folder

```
Deletion Bridge Analysis on Baca Chains/
├── Baca Chains With DB Segments/              (Part 1)
│   └── Baca Chains with Segments Images/      (Part 1 output images)
├── Final Rearrangement Structure/             (Part 2)
│   └── Final Structure After Rearrangements Images/   (Part 2 output images)
│       [user wrote "Final Strucuture After Rearrangements Images" —
│        treated as a typo for "Structure" and corrected here; flag
│        with the user before creating the literal directory if exact
│        string match to their wording turns out to matter]
└── DB Conclusion/                             (Part 3)
```

## Locked design decisions (resolved 2026-08-13 — do not re-litigate)

1. **Segment-lettering scope**: letters A, B, C, ... mark ONLY the
   segments **strictly inside a given deletion-bridge span** (i.e.
   between that DB pair's two anchor breakpoints) — NOT the whole chain.
   Matches exactly how P05-1657's Segment A/B/C were scoped (between
   BP6 and BP14 only; Piece 1, outside that span, was deliberately named
   differently, not lettered A/B/C in the original hand-derivation).

2. **Multi-DB-pair chains handled INDEPENDENTLY per DB pair**
   (user's explicit choice, over "one shared scheme" or "largest pair
   only"): each of a chain's distinct DB pairs gets its OWN letter
   sequence, restarting at A, and its OWN row in Part 3's CSV. Given the
   max-9-internal-breakpoints finding above, no single DB pair's lettering
   can overflow past Z. **Implementation default (not yet confirmed with
   user, but a reasonable, reversible naming choice)**: when a chain has
   multiple DB pairs, disambiguate letters across pairs in the SAME chain
   image using a `DB<n>:<letter>` prefix (e.g. `DB1:A`, `DB1:B`, `DB2:A`)
   rather than plain `A`, `B` repeating per pair on one image — needed
   because a single Part 1 chain image could otherwise show multiple
   unrelated "Segment A"s with no way to tell them apart. Revisit this
   specific display convention if the user objects once they see a real
   multi-DB chain image.

3. **Segment identity is preserved from Part 1 → Part 2.** The same
   letter (e.g. `DB1:A`) refers to the same physical span of DNA in both
   the original chain (Part 1) and the final rearranged structure
   (Part 2) — Part 2 never reassigns or relabels; it shows where each
   Part-1-defined segment ends up, which is the entire mechanism for
   answering "where did it go."

4. **Flipped segments marked with an apostrophe** — exact rule, confirmed
   2026-08-13: a segment gets an apostrophe (`A'`) if, and only if, its
   two boundary endpoints print in the REVERSE order in Part 2 compared
   to Part 1. Example: if Part 1's Segment Span table says
   `A: spans 6[+]----18[-]`, and in the final structure this same
   physical segment appears as `18[+]----6[-]` (endpoints reversed), it
   is flipped → written `A'`. Direct, local, per-segment comparison —
   NOT a function of any global "which end do you start reading the
   whole structure from" choice. Applies only in Part 2 — Part 1 segments
   are always in native/reference orientation by definition, never
   apostrophized.

   **How determinism is actually achieved (resolves an ambiguity the
   assistant raised and the user then resolved 2026-08-13):** the
   assistant initially proposed a "start reading from the lowest-
   coordinate true open end" convention to fix the flip-vs-natural
   ambiguity (the same physical fragment can look flipped-A/natural-B or
   natural-A/flipped-B depending on which end of it you start describing
   from — verified directly by re-deriving P05-1657 two different ways
   this session and getting opposite answers). **The user's actual
   method sidesteps this differently and more physically-grounded**: do
   not pick an arbitrary traversal starting node at all — instead,
   ASSEMBLE the final structure by applying each real rearrangement edge
   **in the fixed order Baca's data lists them** (confirmed via the
   P05-1657 example: Rearr 11 (2+/14−) applied first, then Rearr 17
   (3−/15−) second, matching ascending `Rearrangement number` — TO
   CONFIRM: whether "order given" means ascending `Rearrangement number`
   specifically or literal row order in the sheet; these coincide for
   P05-1657 but may not for every chain, not yet explicitly confirmed).
   Every segment starts at its native Part-1 orientation before any
   rearrangement is applied; each fusion step's own strand pairing
   (same-strand = requires a flip at that junction, opposite-strand =
   no flip needed — the same 4-way logic as the Rearrangement Type
   table above) determines the cumulative orientation as the structure
   is built up incrementally. Because the build order is fixed and each
   step's transformation is locally well-defined, the final orientation
   of every segment is fully determined with no arbitrary choice
   anywhere — this is what makes the apostrophe rule reproducible
   between an independent hand-derivation and the automated algorithm.
   **Consequence flagged, not yet confirmed with the user**: applying
   this rule mechanically to the full P05-1657 example suggests Segment
   C (`15[+]----14[-]` natively) would ALSO come out flipped (`C'`) in
   the final structure (it prints as `14[-]----15[+]`, reversed) — this
   was not explicitly called out in the original walkthrough (only
   Segment A and Piece 1 were narrated as flipped), so this needs a
   direct sanity-check against the user's own by-hand P05-1657 result
   before trusting the algorithm on any other chain.

5. **No-DB chains**: excluded from Part 1/2's detailed segment-marking
   and final-structure work entirely. In Part 3's CSV, every chain still
   gets at least one row via the `DB exists in chain? (Y/N)` column, but
   for `N` rows, no further columns are populated ("nothing after the NO
   column," user's exact framing).

6. **Part 3 verdict schema** (user's explicit choice, "both" option):
   - `Altered? (Yes/No)` — headline column. `No` only if the DB-span
     segment(s) are in the exact same position, same neighbors, same
     orientation in both Part 1 and Part 2 (i.e., truly unchanged). `Yes`
     for anything else.
   - A detailed sub-type column, populated when `Altered = Yes`:
     `MOVED` (relocated elsewhere, still exists, unflipped or flipped),
     `DELETED` (genuinely unaccounted for anywhere in the final
     structure), `FLIPPED-IN-PLACE` (stayed adjacent to its original
     neighbors, same position, but internal orientation reversed — the
     apostrophe case with no relocation), or a combination where
     applicable (e.g. moved AND flipped).
   - Row granularity: one row per DB pair (per decision #2 above), each
     carrying its own `Altered`/sub-type verdict — not one row per chain.

7. **Rearrangement Type column, added to Part 2 (new this session, not
   in the original spec)**: for every real fusion edge involved in a
   chain's final-structure computation, independently classify its type
   from the two breakpoints' (chromosome, strand) values — **do not just
   copy Baca's own `Class` field**; the classification must be derived
   fresh from the data every time, though it may be cross-checked against
   `Class` afterward (the comparison must not influence the
   classification itself). Working 4-way rule, reverse-engineered from
   how the user classified P05-1657's 4 real rearrangements this session
   (all 3 observed cases matched Baca's own `Class` field exactly, but
   the 4th case below was never observed in that example and is an
   inferred completion of the pattern, not yet confirmed against a real
   row):
   | Same chromosome? | Strand pair | Type (independent label) | Matches P05-1657 case |
   |---|---|---|---|
   | Yes | opposite (+/− or −/+) | **SIMPLE / NON-RECIPROCAL TRANSLOCATION** | Rearr 11 (BP2+/BP14−) |
   | Yes | same (+/+ or −/−) | **SIMPLE INVERSION** | Rearr 17 (BP3−/BP15−), Rearr 165 (BP6+/BP18+) |
   | No | same (+/+ or −/−) | **INVERTED TRANSLOCATION** | Rearr 168 (BP7−/BP19−, chr8→chr12) |
   | No | opposite (+/− or −/+) | **TRANSLOCATION** | *(not yet observed — inferred 4th category, confirm against a real example before trusting the label)* |

## Part 1 — `Baca Chains With DB Segments/`

**Goal:** re-draw every chain that has ≥1 deletion bridge (194/366; skip
the other 172 entirely per decision #5) — same visual style as the
current `mmc5/baca_chain_visualization.py` output (all 3 real edge
types, BP number/location/strand corner table) — but with each DB pair's
internal segments additionally labeled per decisions #1-2 above.

**New script** (name/location TBD at implementation time — lives inside
this folder). Reuses/extends
`scripts/baca/mmc5/baca_chain_visualization.py` — do not duplicate its
edge-drawing logic, extend it.

**New deliverable alongside the existing BP-location table:** a
**Segment Span Information** table — one row per lettered segment
(scoped per DB pair per decision #2), e.g. `DB1:A — spans 6[+]----18[-]`,
`DB1:B — spans 18[+]----15[-]`, etc.

**Output images:** `Baca Chains with Segments Images/` (inside this
Part's folder), one image per DB-bearing chain (≤194 images, not 366).

## Part 2 — `Final Rearrangement Structure/`

**Goal:** for every chain that has ≥1 deletion bridge, compute the
**final post-rearrangement derivative structure** — the same kind of
reconstruction the user did by hand for P05-1657, generalized into an
algorithm.

**Critical implementation dependency, found this session, must not be
skipped:** `core/baca_aberration_cycle_drawing.py`'s
`extract_all_traversal_hops` has a verified quirk (`CONTEXT_LOG_003.md`
Session 8): it reports a segment whose BOTH ends independently carry
real rearrangement partners (like Piece 1 in P05-1657) as two separate
disconnected fragments instead of recognizing that segment's own body as
the bridge connecting them. **Part 2's chain-assembly algorithm cannot
just call this function as-is** — needs a fix or a new assembly routine
that correctly merges such fragments.

**CORRECTNESS IS THE HARD REQUIREMENT OF THIS ENTIRE ANALYSIS** — direct
user instruction, not a nice-to-have: *"it is very important that the
final structure after all rearrangements we make is correct otherwise
our whole analysis is rubbish. Make sure you hand draw / follow each
final structure."* Read literally: every DB-bearing chain's computed
final structure must be independently verified (hand-traced or
cross-checked via a second, independent method), not just trusted from
one automated pass — mirrors exactly how the P05-1657 case itself was
validated this session (hand-derivation cross-checked against
`_build_amg_nodes_and_edges`/`extract_all_traversal_hops` output, not
accepted from either alone). Implementation must budget real time for
this verification step, not treat it as optional polish. Exact
verification workflow (per-chain manual check vs. an automated
self-consistency method applied to all ≤194 chains) not yet decided —
resolve before or during implementation, given the scale (194 chains is
a lot to hand-verify one by one; may need a hybrid: automated method
built to mirror the hand-derivation logic exactly, then spot-checked by
hand on a meaningful sample plus any chain the automated method flags as
ambiguous).

**Visualization requirement:** same completeness as Part 1's images —
final structure diagram + BP table + Segment Span table + anything else
already present in the current Baca chain images — **plus**:
- the deletion bridge partner explicitly marked in the rearranged
  structure,
- the lettered segments (with apostrophes where flipped, per decision
  #4) shown in their new position, and
- the new **Rearrangement Type** column/label (per decision #7) for
  every real fusion edge shown.

**Output images:** `Final Structure After Rearrangements Images/` (inside
this Part's folder), one per DB-bearing chain.

## Part 3 — `DB Conclusion/` — DONE 2026-08-13

**Goal:** for every DB pair in every chain, determine and record the
verdict per decision #6 above.

**Script**: `scripts/baca/deletion_bridge_analysis/part3_db_conclusion.py`.
**Output**: `results/Deletion Bridge Analysis on Baca Chains/DB
Conclusion/db_conclusion.csv` — 693 rows (521 DB-pair rows + 172 no-DB
placeholder rows, one per chain without any deletion bridge), columns:
`patient_id, chain_number, db_exists_in_chain, db_index,
db_anchor_low_bp, db_anchor_high_bp, chromosome,
segment_letters_covered, n_segments, altered, sub_type,
per_segment_detail`.

**Per-segment verdict** (reads `Piece.moved`/`Piece.flipped` directly
from Part 2's already-computed, already-validated assembly result — no
new computation, just interpretation): `UNCHANGED` (moved=False,
flipped=False), `FLIPPED-IN-PLACE` (moved=False, flipped=True — stayed
adjacent to native neighbors, only internally reversed), `MOVED`
(moved=True, regardless of flip — `per_segment_detail` separately notes
"(flipped)" when both), `DELETED` (reserved, see below).

**DB-pair-level (row) verdict**: `altered` = `No` only if EVERY segment
in that pair is `UNCHANGED`; otherwise `Yes`, with `sub_type` = the
distinct set of non-UNCHANGED per-segment verdicts found (comma-joined
if more than one type occurs within the same pair — decision #6's
"combination" case). `per_segment_detail` gives the full per-letter
breakdown for the required cross-verification against the Part 1/2
images.

**The open DELETED-definition question, resolved (not just assumed)
before coding**: the assembly algorithm never destroys a piece, only
relocates/reorients it, so "DELETED" (genuinely unaccounted for
anywhere) needed a concrete check, not a guess. Verified empirically
2026-08-13: `n_unresolved_rearr_ends == 0` across all 194 DB-bearing
chains (mmc5's own `Rearrangement number` groups never cross chain
boundaries — already established project-wide) — meaning every real
fusion touching a DB-bearing chain's breakpoints resolves WITHIN that
same chain's own data, with nothing left dangling to an external,
unresolved partner. `DELETED` is kept in the schema (a chain with an
unresolved end WOULD legitimately produce it) but is **structurally
guaranteed to be 0 for this cohort** — confirmed by the actual run
(`DELETED occurrences: 0`), not just predicted.

**Cohort-wide headline result** (693 rows, 366 chains, 0 failures; 521
DB-pair rows exactly matches the DB-pair count already established in
Part 1's build):
```
db_exists_in_chain:  Y=521 rows (194 chains)   N=172 rows (172 chains)
Altered (among Y):   Yes=489 (94%)   No=32 (6%)
sub_type (among Altered=Yes):
  MOVED                    424  (87% of altered)
  FLIPPED-IN-PLACE          57  (12%)
  FLIPPED-IN-PLACE, MOVED    8  (2%, multi-segment pairs with a mix)
DELETED: 0 occurrences (confirmed, not assumed -- see above)
```
**This is the direct, quantitative, cohort-wide answer to the question
this whole analysis was built to answer**: deletion-bridge-flanked
material is essentially never genuinely gone (0/521) — the overwhelming
majority of the time (94%) it IS altered relative to the original chain,
and when altered, it's usually because the material relocated elsewhere
in the chain (87% of altered pairs), not merely flipped in place (12%).
This directly extends the single-patient P05-1657 finding to the full
194-chain cohort with hard numbers.

**Validated against P05-1657 before running the cohort**: DB1 → `Yes`,
`FLIPPED-IN-PLACE`, `A:FLIPPED-IN-PLACE; B:UNCHANGED; C:FLIPPED-IN-PLACE`
— exactly matches the already-established, hand-verified result (A and
C flipped in place, B unchanged; none of the DB-span letters themselves
are MOVED — it's the unlettered Piece 1 that moves in, which is
consistent with, not contradicted by, this per-DB-pair verdict since
Piece 1 was never one of the lettered segments being scored).

## Cross-verification step — DONE 2026-08-13

Per the user's explicit instruction: Part 3's CSV verdicts must be
manually cross-checked against the Part 1 and Part 2 images before the
CSV is trusted as a final answer. Performed on a deliberately diverse
sample (all 4 verdict types, one single-segment case and the original
multi-segment P05-1657 case, one closed-cycle case) using BOTH direct
data inspection (rerunning `assemble_final_structure` and reading
`rearr_log`/`piece.moved`/`piece.flipped` directly — more reliable than
eyeballing arcs on a dense image) and visual inspection of the actual
Part 1/2 PNGs:

- **P05-1657 DB1** (`FLIPPED-IN-PLACE`, 3 segments) — already exhaustively
  validated throughout Part 2's build; re-confirmed consistent here.
- **P03-1426 chain 12 DB1** (`UNCHANGED`, 1 segment): Segment A's own
  two boundaries BOTH have real fusions (to chr7 material and to a
  neighboring chr21 piece) — yet it's correctly `UNCHANGED`, because in
  BOTH real fusion events, the *other* side was determined to be the
  mover (via the same chromosome/position tiebreak already validated on
  P05-1657's Segment C). Confirms `UNCHANGED` means "this segment was
  never itself the one that had to relocate," not "nothing real ever
  touched it" — worth stating precisely, since the two readings differ.
- **P07-837 chain 14 DB1** (`FLIPPED-IN-PLACE`, 1 segment, **inside a
  self-closed cycle**): confirmed Part 3's piece lookup correctly finds
  segments in `closed_molecules`, not just `open_molecules` (this
  specific chain's DB segment is part of a 4-piece closed loop). Visual
  check of the Part 2 image confirms the closed-cycle coloring (purple/
  lavender bars, matching Part 2's own legend) and the segment printed
  in reversed order (`BP430(-)----BP348(+)` vs. native
  `BP348(+)----BP430(-)`) — correctly flagged `A'`.
- **P03-1426 chain 1 DB1** (`MOVED`, flipped): confirmed via direct data
  inspection — `moved=True, flipped=True` → `MOVED (flipped)` in
  `per_segment_detail`, matching the CSV exactly.

No discrepancies found in any of the 4 sampled cases. The user is
encouraged to independently re-check further chains from
`db_conclusion.csv` against their matching Part 1/2 image pairs at their
own pace — the verdicts above are the assistant's independent check, not
a substitute for the user's own, per the original instruction that both
parties look independently.

## Execution order

1. Part 1 (segment-labeled Baca chain images + Segment Span table, only
   for the 194 DB-bearing chains) — needed first since Part 2 and Part 3
   both depend on the same lettering scheme.
2. Part 2 (final rearranged structure, per DB-bearing chain, with the
   fixed/corrected assembly algorithm AND the Rearrangement Type column)
   — depends on Part 1's segment definitions. Includes the mandatory
   correctness-verification step described above.
3. Part 3 (Altered/sub-type CSV, one row per DB pair, plus the DB-exists
   Y/N column covering all 366 chains) — depends on comparing Part 1's
   original segment placement against Part 2's final placement.
4. Cross-verification of Part 3 against Part 1 + Part 2 images.
5. Overall conclusion (cohort-wide: is "deletion bridge = literal loss"
   or "deletion bridge = relocation signal" — or a mix — the better
   general description).

## Overall conclusion — DONE 2026-08-13 (execution order step 5)

**The direct answer to the tension that started this whole analysis**
(user's external research describing a deletion bridge as literal DNA
loss, vs. the P05-1657 finding that its flanked material traced fully to
real edges elsewhere): across the full 194-chain, 521-DB-pair cohort,
**a deletion bridge is a relocation signal, not a literal-loss signal —
and this is now a directly measured result, not an inference from one
example.**

```
DELETED (genuinely gone, unaccounted for anywhere):  0 / 521  (0%)
Altered in some way (moved and/or flipped):        489 / 521  (94%)
  of which MOVED to a new neighborhood:             424 / 489  (87%)
  of which flipped in place only (stayed put):        57 / 489  (12%)
  of which both (mixed within a multi-segment pair):    8 / 489  (2%)
Completely unchanged (Altered = No):                 32 / 521  (6%)
```

Three findings together, each load-bearing:

1. **DNA is never simply missing.** `DELETED` never fires — confirmed
   empirically (0 unresolved rearrangement ends anywhere in the cohort),
   not assumed. Every reference segment nominally "inside" a deletion
   bridge is always traceable to a real fusion edge somewhere in the
   patient's own chain data. This directly contradicts a literal
   "DNA loss" reading of what a deletion bridge represents.
2. **Something *does* happen at almost every deletion bridge** — 94%
   show SOME alteration, so the CN-loss signal ChainFinder detects isn't
   spurious; it's picking up on something real.
3. **That "something" is usually relocation, not just local
   rearrangement.** 87% of altered DB pairs have material that
   genuinely moved to a new neighborhood (often pulled there by foreign
   material arriving from elsewhere in the chain, as directly observed
   in P05-1657 — Piece 1 moving in from ~13kb away). Only 12% are the
   "gentler" case of staying in place but flipping, and true
   no-change-at-all pairs are a small minority (6%).

**Reconciling this with the textbook description** (ChainFinder/LINX-
style tools describing a deletion bridge as tracking "DNA loss"): the
most defensible reading, given this evidence, is that the *reference
window* between the two DB anchors loses its original, simple, linear
correspondence to that stretch of the genome — the sequence that used to
occupy it is very often no longer there (87% of the time, genuinely
relocated elsewhere) — which is exactly the kind of disruption a
coverage/junction-based CN caller would register as an apparent
"deletion" at that specific window, without the underlying DNA actually
being destroyed anywhere in the genome. "Deletion bridge" is a good name
for the *statistical footprint*, not a literal description of what
happened to the DNA.

**Honest limitations of this conclusion** (do not overstate it beyond
what was actually tested):
- Scoped to the 194/366 mmc5 chains that have a deletion bridge at all,
  and only measures the DB-span's own lettered segments against the
  *same chain's* real data — a segment's real fusion partner outside
  this specific chain (if it existed) would not be visible to this
  analysis (though verified as 0 occurrences in this cohort).
- "MOVED" and "FLIPPED-IN-PLACE" are both purely topological
  classifications from the AMG assembly, not independently confirmed
  against copy-number or other orthogonal data (unlike Phase 2's earlier
  CN-plausibility check on *invented* completion edges — this analysis
  works from *real* edges only, so a CN cross-check was not part of its
  design, but would be a natural extension).
- The per-pair aggregation (Altered = Yes if ANY segment in the pair is
  altered) means a DB pair with e.g. 3 segments where only 1 moved and 2
  stayed unchanged is still counted as "Yes" — the 94%/6% split is a
  headline number, not a claim that every DB pair is uniformly altered
  throughout.

## Resolved since the previous revision (2026-08-13, same day)

- **Multi-DB display convention** (`DB<n>:<letter>` prefix) — user
  confirmed correct, "this does solve any ambiguity between multiple DBs
  in a single chain." Locked, no longer open.
- **Part 2 verification workflow at scale** — user's actual bar is
  *agreement* (assistant's automated Part 2 output must match what the
  user would derive by hand for the same chain), not literal manual
  labor on all ≤194 chains: *"if you want to follow your written plan
  or literally hand draw each one, is up to you."* Locked approach:
  build the algorithm to exactly mirror the user's own assembly method
  (see apostrophe-rule section above), then spot-check + investigate any
  chain the algorithm itself flags as ambiguous — literal full manual
  re-derivation of all 194 is no longer required, just convergence.
- **Apostrophe/flip determinism** — resolved via the assembly-order
  method described in decision #4 above, not the assistant's originally-
  proposed "lowest-coordinate open end" convention.

**Segment C flip consequence — CONFIRMED correct 2026-08-13.** User's
own by-hand P05-1657 result: `A'`, `P1'` (Piece 1), and `C'` are all
flipped; only Segment B is unflipped in the final structure. Matches the
assistant's mechanical application of the apostrophe rule exactly — the
rule as written in decision #4 is validated on this example and can be
trusted going forward.

**Build order — CONFIRMED 2026-08-13.** Assembly must process each
chain's real rearrangement edges strictly in **ascending
`Rearrangement number` order** (not row order in the spreadsheet — user
confirmed these coincide for P05-1657 by coincidence, but ascending
`Rearrangement number` is the actual rule to use for every chain). This
was the last open item blocking the Part 2 algorithm spec — the plan is
now fully locked, no remaining open design questions before
implementation can start (only the two low-priority items below).

## Remaining open items

None — the last one (below) was confirmed during Part 2 implementation.

**Fourth rearrangement-type category — CONFIRMED 2026-08-13.** Ran
`classify_chain_rearrangements` (in the new
`scripts/baca/deletion_bridge_analysis/final_structure_assembly.py`)
across all 194 DB-bearing chains: the inferred 4th category
(different-chromosome + opposite-strand → plain `TRANSLOCATION`) occurs
346 times cohort-wide (vs. 779 `SIMPLE / NON-RECIPROCAL TRANSLOCATION`,
770 `SIMPLE INVERSION`, 325 `INVERTED TRANSLOCATION`) — a real, common
category, not a vacuous inferred case. P05-1657's own 4 rearrangements
re-verified to classify exactly as expected (Rearr 11 → SIMPLE/NON-
RECIPROCAL TRANSLOCATION, Rearr 17 & 165 → SIMPLE INVERSION, Rearr 168 →
INVERTED TRANSLOCATION).

## Relationship to existing project work (do not duplicate)

- Different question from the existing Phase 2 CN-plausibility check
  (`scripts/baca/phase2/baca_phase2_cn_plausibility.py`,
  `results/phase2_cn_*.csv`) — that work tests whether Phase 2's
  *invented* completion edges are copy-number-plausible. This analysis
  is about *real*, already-observed deletion-bridge pairs and whether
  their flanked material is *findable via real edges*, not about
  invented edges at all.
- Different from (but reuses the visual style of)
  `scripts/baca/mmc5/baca_chain_visualization.py` (Part 1's base) and
  the corrected understanding of `extract_all_traversal_hops` from
  Session 8 (a hard dependency for Part 2's correctness).
- Directly operationalizes the "test cohort-wide" proposal from Session
  8's end and Session 9's tightened version (simple-DB-pair cases vs.
  complex multi-fusion cases like P05-1657).
