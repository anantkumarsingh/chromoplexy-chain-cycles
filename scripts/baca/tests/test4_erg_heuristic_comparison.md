# Test 4: ERG-Heuristic 3-Way Comparison

**Script:** `test4_erg_heuristic_comparison.py`
**Output:** `results/stat_test4_erg_heuristic_comparison.csv`
**Image:** `results/stat_test4_erg_heuristic_agreement.png`

## 1. What the test is

Before any of the mmc5/Phase 1-3 closure work existed, this project already
had a simpler connectivity-based heuristic for ERG+ patients:
`chromoplexy_embedded` (k_local>2 AND chain_size>1) vs. `simple_fusion` —
built specifically to reproduce Baca's own published "15 of 26 ERG fusions
arose in the setting of chromoplexy" statistic, and it does reproduce that
number exactly. This test asks whether that older, simpler heuristic
actually agrees with the newer, more rigorous closure-based criteria
(real/obligate/probable) or with Baca's own size criterion, now that we
can check formally.

We used the same **Fisher's exact test + Cohen's kappa** pairing as tests 2
and 3, for the same reason: this is an even smaller population (N=23, the
subset of patients where the ERG heuristic happens to be populated), so
several cells are single digits or zero, ruling out chi-square, and raw
percent agreement would again be misleading since `chromoplexy_embedded`
already dominates the group (15 of 23).

## 2. What the test does

Restricts to the 23 patients where `erg_fusion_type_heuristic` is not null
(from `results/erg_chain_comparison_ours_vs_baca.csv`, joined into
`results/phase3_patient_chromoplexy_summary.csv`). Treats
`chromoplexy_embedded` as the "positive" class. For each of `real`/
`obligate`/`probable` × `strict`/`loose` (6 flags) plus
`baca_chain_has_chromoplexy` (1 flag) = 7 tests: builds a 2x2 table against
the heuristic, runs Fisher's exact test, computes Cohen's kappa.

## 3. CSV column reference

| Column | Meaning | Source |
|---|---|---|
| `flag` | the chromoplexy flag column being compared against the ERG heuristic | `phase3_patient_chromoplexy_summary.csv` |
| `lens` | `real`, `obligate`, `probable`, or blank (for `baca_chain_has_chromoplexy`) | same |
| `threshold` | `strict`/`loose`, or blank where not applicable | same |
| `n_embedded_flag_pos` | count where heuristic=`chromoplexy_embedded` AND flag is True | derived |
| `n_embedded_flag_neg` | heuristic=`chromoplexy_embedded`, flag False | derived |
| `n_simple_flag_pos` | heuristic=`simple_fusion`, flag True | derived |
| `n_simple_flag_neg` | heuristic=`simple_fusion`, flag False | derived |
| `odds_ratio` | Fisher's exact odds ratio | `scipy.stats.fisher_exact` |
| `p_value` | two-sided Fisher's exact p-value | `scipy.stats.fisher_exact` |
| `cohen_kappa` | chance-corrected agreement, -1 to +1 | computed as `(po - pe) / (1 - pe)` |

## 4. Step-by-step: what the test does

1. Load `results/phase3_patient_chromoplexy_summary.csv` (52 patients).
2. Restrict to the 23 rows where `erg_fusion_type_heuristic` is populated (15 `chromoplexy_embedded`, 8 `simple_fusion`).
3. For each of the 7 flags, build the heuristic × flag 2x2 table, run Fisher's exact test, compute Cohen's kappa.
4. Write all 7 rows to the output CSV.
5. Render a bar chart of the 7 kappa values.

## 5. Limitations and Notes

- N=23 is small, and several cells are 0 (e.g. `real/strict`: 0 `simple_fusion` patients positive) — several odds ratios are technically infinite and should be read alongside the raw counts in the CSV, not in isolation.
- No multiple-testing correction is applied within this file (7 Fisher's exact tests here, on top of 6 from test 2 and 12 from test 3 = 25 tests across tests 2-4 combined). The combined Benjamini-Hochberg pass across tests 2-4 is planned as the next step now that all three are done — see `test2`/`test3`'s Limitations sections for the same note.
- This heuristic (`chromoplexy_embedded`/`simple_fusion`) was built specifically to reproduce Baca's own "15 of 26" ERG statistic — so a high agreement with `baca_chain_has_chromoplexy` here is partly circular (it was tuned to match a Baca number, not derived independently), not a fully independent validation.
- Only 23/52 patients (and only 23/26 clinically-ERG+ patients) have this heuristic populated — see `results/erg_chain_comparison_ours_vs_baca.csv` for why the other 3 ERG+ patients are missing from that file.

## 6. Analysis Corner

Step-by-step, what the results actually mean:

1. **Agreement is modest across the board — never above 0.51 kappa** ("moderate" at best on the conventional Landis & Koch scale, most values "fair" or worse). The best agreement is `obligate/strict` and `probable/strict` (kappa=0.506, p=0.026, the only two statistically significant comparisons in this test) — the heuristic aligns best with the strict-threshold invented-completion lenses, not with real data alone.
2. **`real/strict` has the lowest agreement of the strict-threshold comparisons** (kappa=0.202, not significant, p=0.257) — the connectivity heuristic (which only asks whether a chain *reaches* multiple chromosomes) and real closed-cycle evidence (which requires the cycle to actually *close*) diverge the most here, consistent with the project's earlier, separate finding that 11/15 `chromoplexy_embedded` calls never close into a real cycle.
3. **The loose-threshold comparisons (`obligate/loose`, `probable/loose`) have the weakest kappa of all** (0.157, not significant) — again a ceiling effect: nearly every patient in this ERG+ subset qualifies at the loose threshold regardless of heuristic label (15/15 embedded, 7/8 simple both positive), leaving little room for the two labels to disagree usefully.
4. **`baca_chain_has_chromoplexy` agreement is moderate but not significant** (kappa=0.303, p=0.111) — some real disagreement exists here (6/8 `simple_fusion` patients still meet Baca's size criterion), interesting because it shows the heuristic and Baca's own size criterion aren't simply restating each other either, despite the heuristic being tuned to reproduce a different Baca statistic.
5. **Bottom line:** the older connectivity-based ERG heuristic is at best moderately related to any of the newer closure-based criteria, and the relationship is strongest with the strict-threshold invented-completion lenses (obligate/probable), not with real data. This supports treating `chromoplexy_embedded`/`simple_fusion` as a distinct, connectivity-only signal rather than a proxy for cycle-closure-based chromoplexy — the two are correlated but clearly not interchangeable, echoing what test 2 already found for Baca's criterion more broadly.
