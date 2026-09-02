# Test 6: Method-Reliability Diagnostics

**Script:** `test6_method_reliability.py`
**Output:** `results/stat_test6_method_reliability.csv`
**Images:** `results/stat_test6_probable_pct_distribution.png`, `results/stat_test6_probable_pct_vs_dangling_ends.png`

## 1. What the test is

This test is different in kind from tests 1-5 — it doesn't test a
chromoplexy hypothesis at all. It diagnoses the Phase 2 obligate/probable
completion *method itself*: when we report "the probable structure for
this chain is X," how confident should anyone reading that be? `X` is just
whichever structure happened to be the most frequent among however many
valid matchings existed for that chain — if the runner-up structure is
almost as frequent, "probable" is a much weaker claim than if the reported
structure dominates. Likewise, the Phase 2 implementation notes flagged
two possible failure modes at the time (no unique most-frequent structure;
no unique lexicographic-minimum obligate structure) but never checked how
often those failure modes actually occur in the real 366-chain dataset.

No inferential test is chosen here on purpose — this is a descriptive
audit (distributions and frequency counts), plus one Spearman correlation
to quantify a relationship (does confidence degrade with chain complexity)
rather than to test a yes/no hypothesis.

## 2. What the test does

Restricted to the 354/366 `phase2_enumerable` chains:

1. Summarizes the distribution of `probable_pct` (mean, std, min, quartiles, max).
2. Counts how often `probable_is_tied=True` (no unique most-frequent structure) and `obligate_max_span_is_unique=False` (the lexicographic obligate rule hit a tie) actually occur.
3. Correlates `n_dangling_ends` against `probable_pct` and `probable_n_distinct_structures` with Spearman's rho.

## 3. CSV column reference

This test's CSV is a flat metric/value table (not repeated per lens/threshold like the other tests, since it diagnoses the method rather than testing a chromoplexy hypothesis).

| Row (`metric`) | Meaning | Source |
|---|---|---|
| `probable_pct_mean` / `_std` / `_min` / `_25pct` / `_50pct` / `_75pct` / `_max` | summary statistics of `probable_pct` across the 354 enumerable chains | `probable_pct` in `phase2_obligate_probable_completion.csv` |
| `n_probable_is_tied_true` / `_false` | count of chains where the most-frequent structure was/wasn't tied with another structure at the same frequency | `probable_is_tied` in same file |
| `n_obligate_max_span_unique_true` / `_false` | count of chains where the obligate lexicographic rule did/didn't land on a unique answer | `obligate_max_span_is_unique` in same file |
| `spearman_rho_probable_pct_vs_dangling_ends` / `spearman_p_...` | correlation between chain openness and how dominant the reported "probable" structure is | computed, `n_dangling_ends` + `probable_pct` |
| `spearman_rho_n_distinct_structures_vs_dangling_ends` / `spearman_p_...` | correlation between chain openness and how many distinct cycle structures are even possible | computed, `n_dangling_ends` + `probable_n_distinct_structures` |

## 4. Step-by-step: what the test does

1. Load `results/phase2_obligate_probable_completion.csv` (366 chains), keep only the 354 with `phase2_enumerable == True`.
2. Compute `.describe()` on `probable_pct`.
3. Count `probable_is_tied` and `obligate_max_span_is_unique` value frequencies.
4. Run `scipy.stats.spearmanr` for `n_dangling_ends` vs. `probable_pct`, and separately vs. `probable_n_distinct_structures`.
5. Write all 15 metric rows to the output CSV.
6. Render a histogram of `probable_pct` and a scatter plot of `probable_pct` vs. `n_dangling_ends`.

## 5. Limitations and Notes

- This is a descriptive/diagnostic test, not a hypothesis test — there's no p-value to interpret as "significant" except for the two Spearman correlations, and even those are reported for their effect size (rho), not primarily for significance (with N=354 and this strong a relationship, significance was never really in question).
- Restricted to the 354/366 enumerable chains, same as test 1 — the 12 non-enumerable chains have no `probable_pct`/`obligate_max_span_is_unique` values to diagnose.
- The zero-ties results (point 1 in Analysis Corner below) describe this specific 366-chain dataset only — they are not a mathematical proof that ties can never occur in this framework, just an empirical fact about the data actually observed here.

## 6. Analysis Corner

Step-by-step, what the results actually mean:

1. **Neither flagged failure mode ever actually occurs in this dataset.** `probable_is_tied` is `False` for all 354/354 enumerable chains, and `obligate_max_span_is_unique` is `True` for all 354/354 — the tie-breaking concerns raised during Phase 2 development (why the precise lexicographic and frequency-mode rules mattered) turned out not to bite in practice for any real chain in this cohort. Every "probable" and "obligate" structure reported anywhere in this project is a clean, unambiguous answer, not an arbitrary tie-break.
2. **`probable_pct` is bimodal, not smoothly distributed** — see the histogram: a large cluster at exactly 100% (160/354 chains, the ones with only 1-3 total matchings, where the reported structure is trivially "the" structure) and a second broad cluster spanning roughly 37-67% for chains with more matchings. Median is 66.7%, but the IQR (53.3%-100%) shows this splits into two very different regimes rather than one typical value.
3. **`probable_pct` is almost entirely determined by `n_dangling_ends` alone** (Spearman rho=-0.967, p≈3×10⁻²¹²) — the scatter plot shows an almost step-function relationship: chains with the same number of dangling ends cluster tightly at nearly the same `probable_pct`, regardless of which specific chromosomes or topology are involved. This makes sense combinatorially (the total number of matchings, `(m-1)!!`, depends only on `m = n_dangling_ends`), but it's a useful, now-quantified fact: confidence in the "probable" completion is predictable almost purely from how open a chain is, not from anything else about it.
4. **`probable_n_distinct_structures` grows just as predictably** (rho=+0.965) — more open ends means more possible outcomes, exactly mirroring point 3 from the opposite direction (more possibilities dilutes any single structure's share, which is exactly what point 3 shows).
5. **Bottom line, methods-section framing:** the obligate/probable machinery is reliable in the narrow sense that it never hits an ambiguous tie in this dataset, but the *confidence* behind any individual "probable" call should be read alongside `n_dangling_ends` — a chain with many dangling ends can have its reported "probable" structure represent as little as ~37% of the combinatorial space, which is a real caveat for how strongly any single probable-lens chromoplexy call (used throughout tests 1-5) should be trusted for the more complex chains specifically.
