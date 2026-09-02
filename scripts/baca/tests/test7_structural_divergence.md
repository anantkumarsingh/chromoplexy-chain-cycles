# Test 7: Obligate-vs-Probable Structural Divergence

**Script:** `test7_structural_divergence.py`
**Output:** `results/stat_test7_structural_divergence.csv`
**Image:** `results/stat_test7_divergence_by_size.png`

## 1. What the test is

Test 1 checked whether obligate/probable *agree with each other on the
strict/loose chromoplexy flag* (they mostly do — test 1 found obligate's
qualifying chains are always a subset of probable's). This test asks a
finer-grained question: even when they agree on the flag, do they agree on
the actual reported cycle structure? Two chains could both be flagged
chromoplexy-positive while one lens reports "C4+C2" and the other reports
"C6" for the exact same chain — that's a real, substantive disagreement
about the coordination story, even if the binary flag happens to match.

We used the same **Mann-Whitney U + rank-biserial** approach as test 1, for
the same reason (small non-negative integer chain-size counts, unequal
group sizes) — this test is essentially test 1's method applied to a
different outcome variable (structural divergence instead of chromoplexy
qualification).

## 2. What the test does

For every enumerable chain (354/366), compares
`obligate_cycle_structure` against `probable_cycle_structure` as plain
strings — if they differ at all, the chain is "divergent." Then compares
`n_dangling_ends`, `n_rearrangements`, and `k_chromosomes_whole_chain`
between divergent and concordant chains with Mann-Whitney U.

## 3. CSV column reference

This test's CSV is a flat metric/value table (same style as test 6).

| Row (`metric`) | Meaning | Source |
|---|---|---|
| `n_divergent` | chains where `obligate_cycle_structure != probable_cycle_structure` | `phase2_obligate_probable_completion.csv` |
| `n_concordant` | chains where the two structure strings are identical | same |
| `n_enumerable_total` | 354, the denominator for both counts above | same |
| `<var>_median_divergent` / `_median_concordant` | median of that size variable in each group, for `var` in `n_dangling_ends`, `n_rearrangements`, `k_chromosomes_whole_chain` | derived |
| `<var>_mannwhitney_U` | Mann-Whitney U statistic for that variable | `scipy.stats.mannwhitneyu` |
| `<var>_p_value` | two-sided p-value | `scipy.stats.mannwhitneyu` |
| `<var>_rank_biserial_r` | effect size, -1 to +1 | computed as `1 - 2U/(n1*n2)` |

## 4. Step-by-step: what the test does

1. Load `results/phase2_obligate_probable_completion.csv`, keep the 354 `phase2_enumerable` chains.
2. Flag each chain `divergent` if `obligate_cycle_structure != probable_cycle_structure`.
3. Run Mann-Whitney U for each of the 3 size variables, divergent group vs. concordant group.
4. Write 18 metric rows to the output CSV.
5. Render 3 side-by-side box plots (one per size variable) comparing the divergent and concordant groups.

## 5. Limitations and Notes

- This test only compares the *reported* structure strings — it doesn't weigh how different two structures are (e.g. "C4+C2" vs. "C6" is counted the same as "C4+C2" vs. "C4+C2+C2+C2..." even though the latter pair is a much bigger structural gap). A finer divergence-magnitude metric was not built here.
- No multiple-testing correction is applied within this file (3 Mann-Whitney tests). Small relative to the running total from tests 1-6, but still uncorrected.
- Divergence can only be evaluated for the 354 enumerable chains — the 12 non-enumerable chains have no `obligate_cycle_structure`/`probable_cycle_structure` to compare.
- Same caution as test 1's `k_chromosomes_whole_chain` result: `n_dangling_ends`'s relationship to divergence here is (mostly) mechanical, not a free-standing empirical discovery — see Analysis Corner point 1.

## 6. Analysis Corner

Step-by-step, what the results actually mean:

1. **`n_dangling_ends` shows perfect separation (rank-biserial r=-1.000, p<0.0001) — and this part is mechanical, not a real finding.** A direct cross-tab (`n_dangling_ends` × `divergent`) shows: every chain with `n_dangling_ends` ∈ {0, 2} is concordant (160/160), and every chain with `n_dangling_ends` ≥ 4 is divergent (194/194) — a perfect, deterministic split. This is expected combinatorially: with 0 dangling ends there's nothing to complete; with exactly 2 dangling ends there's exactly `(2-1)!!=1` possible matching, so obligate and probable are FORCED to pick the same (only) option; divergence only becomes *possible* once `n_dangling_ends ≥ 4` (3+ matchings to choose between).
2. **But the fact that divergence is realized 100% of the time whenever it's possible (all 194 chains with `n_dangling_ends ≥ 4` diverge, not just some of them) IS a genuine, non-tautological empirical finding worth reporting** — it wasn't guaranteed that obligate's conservative rule and probable's frequency-mode rule would land on different answers every single time they had the chance to; in this dataset, they always do.
3. **`n_rearrangements` shows a real, moderate, non-mechanical effect** (r=-0.257, p<0.0001) — divergent chains skew toward more rearrangements (median 4, but with a long tail up to 36) than concordant chains (median also 4, but with a much tighter, smaller-tailed distribution, max 19). This mirrors test 1's finding that chain size matters, but weaker here — divergence is dominated by the `n_dangling_ends` structural effect described in point 1.
4. **`k_chromosomes_whole_chain` shows a small effect in the OPPOSITE direction** (r=+0.169, p=0.002) — concordant chains have a slightly *higher* median chromosome span (2) than divergent chains (1). This is the one genuinely counter-intuitive result in this test: chains that touch more chromosomes are, if anything, slightly *less* likely to have obligate and probable disagree on structure, not more.
5. **Bottom line:** whether obligate and probable report the same structure for a chain is almost entirely a function of how many dangling ends that chain has (mechanical, and once ≥4, divergence is a near-certainty in this data) rather than the chain's chromosome span or even its full rearrangement count. This is a useful caveat for anyone reading a Phase 2 obligate/probable structure column directly: any chain with 4+ dangling ends should be assumed to have obligate and probable disagreeing on structure specifics (194/194 do), even where they still happen to agree on the coarser strict/loose chromoplexy flag.
