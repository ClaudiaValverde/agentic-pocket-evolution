# Part A v3: orthology O2 (real tier-2 placement, O2_r1 SOGs, D20 test half). Snapshot 20261030T085136Z

**Status: development subset, preliminary; not yet verified.** Written 2026-10-30 17:00 CET (`date`). The v2.1 results
(O1) are not modified. Code commits: `7e550db`, `d106854`, `95b3490`. All bulk data is in `/tmp/pe_big/partA/20261030T085136Z_O2/`;
this directory is 1.3 MB. prereg sha256 `52f94618…0d83`. Numbers come from `results/numbers.json` → `partA_O2` and from
the CSVs here. The comparison table is `o2_vs_o1.csv`, the class moves are `class_moves.csv` / `class_moves_summary.json`.

## What changed vs v2.1
- **Orthology.** `data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz` (157 species, 57,332 root OGs).
- **Tier-1 ID gap.** Resolved with `idgap_mapping.parquet`: 27,087 proteins now map to their gene's O2 member; 1,121
  were dropped because the member was a duplicate or clashed.
- **Tier-2 placement.** Uses `placement/tier2`, the DIAMOND placement into O2 HOGs. This **replaces the PANTHER proxy
  for 365,882 tier-2 proteins.** A further 7,889 Swiss-Prot proteins are O2 members themselves.
- **Unplaced tier-2 proteins** (DIAMOND no-hit) keep the flagged fallback: PANTHER proxy for 16,728, singleton for
  18,217. This is conservative, because grouping more proteins means fewer N pairs.
- **Same-site split.** Uses O2_r1 SOGs at τ = 0.2 (frozen, D8). Tier-2 placed proteins are not SOG members, so their H12
  pairs still use the mapped-Jaccard placeholder (8.2 M placeholder calls vs 1.3 M SOG calls among kept pairs).
  Same-site fractions: H1 0.124, H2 0.130, H12 0.134.
- **Test half (D20).** Taken from the SOG agent's split table `results/sogs/…/O2_r1_split/og_split_o2.parquet`, which
  was present, so it was **used, not recomputed**. Crossover families and non-HOG groups are excluded. Sites: test
  629,206, tune 577,790, excluded 2,014,706.
- **Null and significance.** Recomputed on the O2 context; the null masks homologs by O2 groups. Result: 9,778,355
  significant pairs, against 9,776,291 under O1. The extended calibration here uses 2 × 10⁸ pairs.

## O2 v3 vs O1 v2.1 (pooled tail, E < 0.01 in both directions; 200 bootstrap reps, halves 100)

| metric | O1 v2.1, all | **O2 v3, all** | O1 v2.1, test half (hash split) | O2 v3, test half (D20) |
|---|---|---|---|---|
| F_anc, per family pair | 0.176 [0.169, 0.185] | **0.060 [0.057, 0.063]** | 0.175 | 0.083 [0.077, 0.091] |
| F_analog, per family pair | 0.819 [0.810, 0.826] | **0.935 [0.932, 0.937]** | 0.819 | 0.910 [0.900, 0.916] |
| F_analog, per pair | 0.367 [0.323, 0.416] | 0.367 [0.325, 0.416] | 0.285 | 0.136 [0.104, 0.181] |
| r_adj H1s | 0.0100 [0.0086, 0.0119] | 0.0132 [0.0116, 0.0149] | 0.0107 | 0.0173 [0.013, 0.025] |
| r_adj H2s | 0.0057 [0.0021, 0.0091] | 0.0065 [0.0038, 0.0082] | 0.0019 | 0.0077 |
| r_adj H12s | 0.022 [0.017, 0.029] | 0.023 [0.019, 0.028] | 0.019 | 0.024 |
| r_adj H3j | 2.0e-5 | 1.9e-5 | 1.7e-5 | 7.2e-6 |
| r_adj H4 | 1.0e-5 | 1.1e-5 | 1.5e-5 | 8.4e-6 |
| r_adj N | 7.0e-7 [6.8e-7, 7.2e-7] | 7.0e-7 [6.7e-7, 7.3e-7] | 6.8e-7 | 6.8e-7 |
| F_conv,func, per family pair, all labels (headline) | 0.081 [0.007, 0.178] | 0.095 [0.021, 0.204] (n = 191) | 0.127 (n = 41) | 0.39 (**n = 2**, not interpretable) |
| F_conv,func, specific ligands only (tokens, no metal) | 0.021 [0, 0.062] | 0.025 [0, 0.089] | 0.021 | 0 (n = 2) |
| Species-Henikoff reweighting (sensitivity): F_analog per family pair / per pair | 0.700 / 0.258 | 0.912 / 0.242 | not computed | not computed |
| descent M1, ortholog / paralog (qm90, capped) | 0.774 / 0.780 | 0.932 / 0.936 (**not comparable**, see below) | 0.775 / 0.784 | 0.928 / 0.931 |
| descent: ortholog − paralog at matched identity | −0.76 pp [−1.34, −0.12] | −0.38 pp [−0.59, −0.15] | −0.94 pp | −0.23 pp [−0.62, +0.20] |

How to read it:
- **The rates are stable.** r_N, r_H1s, r_H12s and per-pair F_analog are unchanged within their CIs.
  - Positional orthologs are about 19,000 times more often significantly similar than unrelated pairs (0.013 vs
    7.0e-7).
  - The functional excess (F_conv,func) and its specific-ligand sensitivity are unchanged (0.095 and 0.025).
- **The per-family-pair split moves: F_analog 0.82 → 0.93.** This comes from the weighting, not from the pairs.
  - O2 places tier-2 proteins into fewer, larger root groups: 89,051 groups with sites, against 148,511 under O1.
  - Significant H3j pairs now fall into 137,714 distinct group cells instead of 494,695. Their pair count barely
    changes (1.39 M → 1.30 M).
  - Each cell weighs 1, so ancestry cells lose weight. Within-group cells are about equal: 12,641 under O1, 12,945
    under O2.
  - The headline per-family-pair weighting is therefore **sensitive to how finely the orthology partitions
    families**. The verifier flagged this; it is now quantified. The per-pair rates r_c are the robust Part A
    statement.
- **Descent M1 on O2 is not comparable to O1.** The O2_r1 states apply D11: an A call needs a footprint spanning more
  than 10 residues, while P has no such requirement. This is the asymmetry the verifier predicted, and it removes
  informative non-matches, which inflates M1 (0.93).
  - O2 M1 still shows no ortholog advantage (−0.4 pp), so "ortholog conjecture not supported for presence" holds.
  - The O2 descent tables come from the descent agent's builder, unmodified, run through
    `scripts/partA/17_descent_o2.py`. The wrapper restricts M2 site matching to tier-1/2 universe sites, because
    O2 SOGs contain tier-3/L sites without fingerprints. M1 is unaffected.
  - The fix (give P the same > 10-residue rule) belongs to the SOG and descent agents.
- **F_conv,func on the D20 test half:** only 2 annotated test pairs survive the crossover exclusion, so it is not
  interpretable.

## Pairs that move class because of real placement (`class_moves_summary.json`)
- **Overall.** 9,757,299 pairs are significant under both orthologies. **1,529,928 (15.7 %) change class**, and
  932,692 of those involve a tier-2 protein.
- **Of the 3,970,709 tier-2 pairs that were H12 under the O1 PANTHER proxy:**

  | O2 class | pairs | comment |
  |---|---|---|
  | H12s | 3,509,735 | the proxy and the placement agree |
  | **H3j** | **386,685** | placed into different O2 HOGs; still homologous by joins |
  | H12x | 66,261 | |
  | H1s, H2s, H1x, H2x | 3,548 + 586 + 140 + 37 | now O2 members with an ortholog/paralog call |
  | H3s / H4 | 2,727 / 927 | |
  | **N** | **63** | |

- **Into O2 H12:** 527,054 pairs that were H3j under O1, now joined into one O2 HOG.
- **Net:** real placement mostly reshuffles mass between H12 and H3j. Both count as ancestry, so per-pair F_anc does
  not change. Almost nothing moves into or out of N (63 pairs out and 69 in).

## Inverse-clade-size reweighting (sensitivity; point estimates only)
- **Weighting.** Each site gets w = 1 / (number of universe proteins of the same taxid in its root group); the pair
  weight is w_i · w_j. This is applied jointly to K_c, to the exact totals (weighted `exact_class_totals` and
  `LabelSubsets`, checked against brute force in `test_weighted_exact_totals_match_bruteforce`) and to the
  denominator sample.
- **Limitation.** This is a species-level Henikoff weighting that down-weights in-paralog and strain redundancy. It is
  not a tree-based GSC weighting; there was no time for that.
- **Results:**

  | | O1 | O2 |
  |---|---|---|
  | F_analog, per family pair | 0.819 → 0.700 | 0.935 → 0.912 |
  | F_analog, per pair | 0.367 → 0.258 | 0.367 → 0.242 |
  | F_anc, per pair | 0.615 → 0.727 | → 0.742 |

  Down-weighting redundancy shifts mass toward ancestry at the pair level. No conclusion changes.
- Files: `estimands_cos_henikoff.csv` (O2) and `O1_ref/estimands_cos_henikoff.csv` (O1).

## Deviations (PA ids; for the ledger)
- **PA23:** unplaced tier-2 proteins (DIAMOND no-hit) keep the PANTHER-proxy or singleton fallback (flagged). Placements
  are used whether or not they are `compatible`.
- **PA24:** test half = D20 split table, used as provided. Non-HOG groups are excluded along with crossover groups.
- **PA25:** the O2 null extended calibration uses 2 × 10⁸ pairs (O1 used 10⁹).
- **PA26:** descent on O2 restricts M2 to universe sites (wrapper). O2 M1 is flagged as not comparable to O1 (D11
  asymmetry).
- **PA27:** the species-Henikoff weighting stands in for the inverse-clade-size reweighting (point estimates only).
- **Changes a headline estimate:** the per-family-pair F on O2 (0.935), which comes from the group-granularity effect
  above. **Sent to the PI.**

## Runtime (64-core node at load 230; GPU 3)
| stage | wall |
|---|---|
| context | 200 s |
| null (2 × 10⁸ extended) | 4.5 min |
| significance | 60 s |
| classify | 374 s |
| estimands (200 + 8 × 100 replicates) | about 50 min (CPU contention) |
| functional, descent, figures, numbers | about 5 min |
| Henikoff, 2 runs | 3 min each |
| comparison | 30 s |
| **total elapsed** | 15:43–17:00 |
