# Part A v2: similarity decomposition, functional and directional convergence. Snapshot 20261030T085136Z

**Status: development subset, preliminary.** Not yet verified (docs/VERIFY_CHECKLIST.md). It covers tiers 1–2 and
orthology O1. **v2 changes from v1** (`../PARTA.md`, `../v1/`):
- the real O1 SOGs give the same-site/offsite split (τ = 0.2 from the SOG run's tune half);
- new: functional convergence F_conv,func, the cross-fitted directional test F_conv,dir, the undetected-homology
  check, CIs for strata;
- the acceptance-stat bug is fixed.

Stages 1–5 (universe, kNN, null, significance) are unchanged from v1. Their summaries are copied here.

Code commits: `3f744ad` (stages 6, 7, 8, 14), `3f11d4c` (stages 12, 13, csls 6/7, 10), plus later provenance-only commits for 04, 04b (sidecars regenerated from the cached v1 artifacts; calibration tables are byte-identical) and 11 (`f913f54`). prereg sha256
`52f94618…0d83`. Every number below comes from `numbers_partA.json` (= `results/numbers.json` → `partA`) or from
a CSV in this directory. Tests: `tests/test_partA.py`, 16 pass (`pytest_partA.log`), including the new directional
toy worlds.

## 1. Decomposition (pooled tail, E < 0.01 both directions; FP nominal × λ = 1.72; 200 two-way bootstrap reps)
9,776,291 significant pairs. Expected false positives: 27,276 in N (q0 = 5.34e-9 per pair).

| weighting | F_anc | F_offsite | F_H4 | F_analog (N) |
|---|---|---|---|---|
| **per family pair (headline)** | **0.176 [0.169, 0.185]** | 0.0001 | 0.005 [0.004, 0.006] | **0.819 [0.810, 0.826]** |
| per site | 0.288 [0.278, 0.298] | 0.008 [0.007, 0.009] | 0.004 | 0.700 [0.689, 0.710] |
| per pair (descriptive) | 0.615 [0.566, 0.661] | 0.011 [0.008, 0.015] | 0.007 | 0.367 [0.323, 0.416] |

- **Hidden homology barely changes F_analog.** F_analog corrected for undetected homology (§4, ε = 0.005) is
  **0.815 [0.806, 0.822]** per family pair and 0.365 per pair.
- **Same-site split (SOG at τ = 0.2).**
  - 1.46 M significant H1/H2 pairs were split by SOG identity. 8.06 M pairs, mostly H12 (tier-2 proxy-placed
    proteins, which have no SOG), fell back to the mapped-Jaccard placeholder at the same τ.
  - The denominator split was estimated from 200,000 within-group pairs: same-site fractions H1 0.120, H2 0.128,
    H12 0.120. In that sample, 99.7 % of H1 and 99.5 % of H2 pairs were called by SOG.
- **Rates.**

  | class | r_adj [95 % CI] |
  |---|---|
  | H1s (same-site orthologs) | **0.010 [0.0086, 0.012]** |
  | H2s | 0.0057 [0.0021, 0.0091] |
  | H12s | 0.022 [0.017, 0.029] |
  | offsite H1x / H2x / H12x | 5.5e-5 / 4.4e-5 / 6.1e-5 |
  | N | **7.0e-7 [6.8e-7, 7.2e-7]** |

  Positional orthologs are therefore about 14,000 times more often significantly similar than unrelated pairs.
  Offsite pairs are about 80 times above N but about 180 times below same-site orthologs: being in the same family,
  without being in the same position, adds little.
- **Sensitivity** (`sensitivity_grid.csv`; F_analog per family pair / per pair):
  - E < 1e-4: 0.800 / 0.330;
  - E < 1: 0.796 / 0.381;
  - per-query tail: 0.921 / 0.638;
  - CSLS score: 0.905 / 0.581;
  - review FP definition (bound): 0 / 0.
- **Strata with CIs** (`strata_cos.csv`; 100 bootstrap reps per stratum, both sites in the stratum):

| stratum | F_analog per family pair | F_analog per pair | F_anc per family pair |
|---|---|---|---|
| pLDDT ≥ 90 | **0.588 [0.571, 0.604]** | 0.158 [0.134, 0.181] | 0.401 [0.386, 0.418] |
| both tier 1 | 0.886 [0.879, 0.892] | 0.514 [0.396, 0.651] | 0.106 [0.100, 0.112] |
| both tier 2 (Swiss-Prot) | 0.390 [0.362, 0.412] | 0.091 [0.074, 0.111] | 0.608 [0.586, 0.636] |
| volume small / mid / large | 0.653 / 0.671 / 0.753 (CIs ±0.02) | 0.184 / 0.208 / 0.241 | 0.344 / 0.324 / 0.233 |

## 2. Calibration (unchanged from v1; `fig_null_calibration_cos.png`)
- **Bulk:** calibrated. On held-out 10^6 pairs, observed/expected is 0.98–1.04 for p from 0.1 to 1e-3.
- **Deep tail:** the excess (1.7× at 1e-4, 38× at 1e-6) is entirely N by search.
- **At the threshold:** on 10^9 random pairs, the rate of significance in both directions at E < 0.01 is 9.9e-7,
  which equals K_N/n_N = 7.0e-7 within the k = 1000 limit.
- **Agreed interpretation (PI):** geometric N similarity at the threshold matches the random-pair rate. That is
  consistent with pocket-space degeneracy, so convergence has to be shown by functional excess (§3) and by the
  directional test (§5).

## 3. Functional convergence F_conv,func (`functional_convergence.csv`)
- **Labels** (experimental only; a label attaches to a site when ≥ 50 % of the label's residues lie in it):
  - BioLiP2 CCD ligands, regular and metal;
  - UniProt binding sites with ECO:0000269/ECO:0007744, as ChEBI ids, with common ligands mapped to CCD
    (`config/partA_ligand_map.tsv`);
  - M-CSA reference catalytic residues, as 4-digit EC.

  This gives 23,540 annotated sites (v1: 19,081; BioLiP2 is most of it).
- **Shared class:** a common token, or ECFP4 Tanimoto ≥ 0.5 between non-metal CCDs.
- **Definitions:**
  - s_sig = shared share among annotated significant N pairs. pocketcompare train/val proteins are excluded.
  - s_ctl = shared share among 53,417 random annotated N pairs (seeded; MMseqs2/Foldseek-negative), reweighted to
    the deciles of covariate-predicted similarity μ_ij of the test pairs.
  - δ = s_sig − s_ctl.
  - **F_conv,func = F_analog × max(δ, 0)**: the share of all significant similarity mass carried by N pairs whose
    functional sharing is above the covariate-matched background. It assumes the annotated subset represents all
    significant N pairs; this is stated, not verified.
  - AF = δ/s_sig, the attributable fraction.
  - Category rows are additionally scaled by the share of annotated test pairs in which both sites carry that
    category.
- **CIs:** two-way bootstrap with the same Poisson group multiplicities as the estimand bootstrap, so F_analog and
  δ are resampled jointly.

| subset | category | n sig N pairs | s_sig | s_ctl (matched) | δ [95 % CI] | enrichment | F_conv,func per family pair [CI] | per pair |
|---|---|---|---|---|---|---|---|---|
| all | any | 192 | 0.156 | 0.061 | 0.095 [0.004, 0.207] | 2.55 | **0.078 [0.004, 0.169]** | 0.035 [0.002, 0.078] |
| all | metal | 19 | 0.53 | 0.28 | 0.24 [−0.12, 0.65] | 1.87 | 0.020 [0, 0.053] | 0.009 |
| all | cofactor | 4 | 0.75 | 0.07 | 0.68 (n too small for CI) | 10.8 | 0.012 | 0.005 |
| all | small molecule | 121 | 0.124 | 0.050 | 0.074 [−0.02, 0.19] | 2.50 | 0.038 [0, 0.098] | 0.017 |
| all | catalytic EC | 0 | none | 0.000 | none | none | none | none |
| pLDDT ≥ 90 | any | 100 | 0.120 | 0.049 | 0.071 [−0.02, 0.19] | 2.44 | 0.042 [0, 0.095] | 0.011 [0, 0.027] |

- **Reference shared-class rates** (significant pairs, annotated subset): H1s 0.79 (n = 222), H2s 0.74, H12s 0.80,
  H3j 0.76, H3s 0.80 (n = 5), H4 0.75 (n = 4), H1x 0.50 (n = 6). Homologous significant pairs share function about
  5 times as often as significant N pairs.
- **Without the training exclusion:** n = 210, δ = 0.101, F_conv,func = 0.083.
- **Reading.** Significant analog pairs are enriched for shared experimental ligands at about 2.5×. The excess
  corresponds to at most a few percent of the similarity mass. The CI barely excludes 0 overall and includes 0 at
  pLDDT ≥ 90. **n is the binding constraint:** UniProt experimental sites mostly overlap BioLiP2, and M-CSA
  reference residues fall in only 839 site attachments, giving 0 doubly-annotated N pairs.

## 4. Undetected-homology check (`hidden_homology.csv`)
- **Method:** 200 seeded significant N pairs, re-aligned with Foldseek TM-align mode (`--alignment-type 1
  --exact-tmscore 1`, full chains). HHsearch is not installed; prereg allows "or Foldseek TM-align mode".
- **Result:** 1 of 200 has max(qTM, tTM) ≥ 0.5 (lenient), so **ε = 0.005 (Wilson 0.0009–0.028)**; 0 of 200 meet
  min ≥ 0.5 (strict). The median max TM is 0.10.
- **Reference:** 200 random N pairs give ε = 0.010.
- **Sanity check:** 50 significant H1s pairs have median TM 0.89, and H4 pairs 0.38.
- **Correction:** F_analog × (1 − ε), with the difference moved to H4 (above).

## 5. Directional convergence test (`directional_by_class.csv`, `directional_summary.json`)
- **Method** (`src/pocketevo/directional.py`): for each significant N/H3j/H3s/H4 pair, each orientation (a, b) in
  which b lies in an O1 SOG with ≥ 3 species is tested.
  - The derived clade is the smallest species-tree clade around b's species that holds ≥ 4 other SOG members from
    ≥ 2 species. Members outside it are the sister/outgroup proxy, and at least 2 are required.
  - Derived members are split (seeded) into a selection half that contains b and a test half.
  - A pair is selected if Δ_sel > 0. The test statistic is **Δ = mean z(a, test half) − mean z(a, outside)**, so
    b's own value never enters.
- **Nulls:**

  | null | n | mean | sd | q95 |
  |---|---|---|---|---|
  | (i) significant same-SOG ortholog pairs (a's protein excluded) | 247,901 | 0.24 | 0.63 | 1.34 |
  | (ii) random non-homologous partners of the tested b sites | 296,338 | 0.005 | 0.45 | 0.74 |

- **Decision:** Δ above the larger q95, and BH (q = 0.05) on the conservative p (the max over both nulls; normal-fit
  tail beyond the null's range, PA15).
- **Toy calibration** (tests): the no-convergence BM world reads F_conv,dir = 0 (0 false passes out of 6,400). A
  planted derived-clade convergence is detected with power 1.
- **Real data:**
  - 5,107,143 significant N/H3/H4 pairs, of which **2,627,453 are testable (51 %)**.
  - Tested Δ has mean 0.20 and q95 1.21. That is the same as the ortholog null (0.24 / 1.34): the shift reflects
    lineage similarity carried by selection, not convergence. 3.7 % of selected orientations exceed q95 (below the
    5 % expected under the null).
  - **3 pairs pass**, all N. **F_conv,dir = 1.1e-6 per pair [0, 4e-6]**, 2e-6 per family pair [0, 7e-6]. H3j, H3s
    and H4 have 0 passes.
  - Sensitivity, random null only (anti-conservative: ignores the selection-induced shift): 11,350 orientations pass.
- **Candidates:** `convergence_candidates.parquet` holds the 3 primary passes, all random-null-only passes and the
  top 500 selected orientations by p, with flags `passed_primary` / `passed_random_null_only`. Columns: proteins,
  sites, SOGs, s, z, Δ, p and ligand tokens.
  - Primary passes: Q24298:s0 ↔ E1JGJ8:s0 (SOG HOG0035362:g0, Δ 5.2); P46723:s3 ↔ C0NED7:s2 (HOG0036631:g2, Δ 4.8);
    Q22937:s1 ↔ Q6R5P0:s1 (HOG0029116:g8, Δ 4.7).
  - None of the three has experimental ligand labels. Of the listed candidates with labels on both sides, 5 share a
    token and 2 do not.
- **Conclusion for Part A:** pocket similarity across unrelated families at the significance threshold equals the
  random-pair rate (degeneracy). It carries a modest functional excess (F_conv,func ≈ 0.08 [0.004, 0.17] per family
  pair). It shows **no detectable lineage-directional convergence** (F_conv,dir ≈ 1e-6). Convergence is therefore
  rare at the resolution of this fingerprint and of O1 sampling.

## Deviations (PA ids; mapping to the project ledger DESIGN §12: PA1 (pooled tail) = ledger D12; PA3 (expected FP) = ledger D13)
- **PA1–PA13:** as in `../PARTA.md` (v1), unchanged, with these changes:
  - PA4 is replaced by PA14.
  - PA5 is still active (no O2 placement table).
  - PA13 items are now done, except inverse-clade-size reweighting.
- **PA14 Same-site split.**
  - SOG identity at τ = 0.2 (`results/sogs/…/O1/tau_choice.json`, tuned on the tune half; prereg grid). It replaces
    the prereg default 0.3 used by the v1 placeholder.
  - When either site has no SOG (OG with < 3 species, tier-2 proxy placement), the mapped-Jaccard placeholder at
    τ = 0.2 is used. The source of each call is recorded in `split_source`.
  - Effect: r_H1s 0.031 → 0.010. The F estimands are unchanged to 3 decimals.
- **PA15 Directional p-values.**
  - Normal-fit upper tail of each null beyond its range, because empirical p cannot fall below 1/(n_null + 1) and
    would make BH powerless (shown in the toy). p is the max over the two nulls.
  - Selection on Δ_sel > 0. A pair passes if any orientation passes.
  - Member sets are capped at 64 (seeded).
  - **This determines F_conv,dir. Sent to the PI for the ledger.**
- **PA16 F_conv,func definition.** As in §3 (PI request). Category rows are scaled by category share. Labels add
  UniProt ECO-experimental binding sites and M-CSA reference EC; there is no cross-source ChEBI↔CCD matching
  beyond the map file. The control is a fresh seeded draw of 60,000 (53,417 N after search). **Determines
  F_conv,func; sent to the PI.**
- **PA17 Hidden homology:** Foldseek TM-align mode instead of HHsearch (prereg alternative). Lenient rule. ε is
  applied as a point correction; its own CI is listed separately.
- **PA18 Strata CIs:** 100 replicates per stratum (headline uses 200).
- **Fix:** `rest_sample_accept` now equals accepted/drawn (0.986).

## Open issues
1. O2 placement for tier 2 (H12 is still a proxy, and 85 % of significant same-group pairs are H12).
2. The F_conv,func n is small (192 test pairs). The full atlas plus better cross-source ligand matching would raise
   it.
3. Directional testability is 51 % (SOGs exist only for tier-1 OGs with ≥ 3 species). Tier L will raise the
   species density.
4. Inverse-clade-size reweighting was not done.

## Runtime (v2 stages; 64-core node shared with other agents, GPU 3)
| stage | wall |
|---|---|
| 06 classify | 769 s (SOG split + 200,000-pair within-group sample) |
| 07 estimands | 2,170 s (200 reps primary + 6 × 100 strata) |
| 12 functional | 120 s |
| 13 directional | 217 s (4.17 M orientations, 3.13 M testable) |
| 14 TM-align | 10 s |
| 06/07 csls | about 6 min |
