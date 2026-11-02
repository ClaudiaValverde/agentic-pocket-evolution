# Part A v2.1: test halves, τ sensitivity, functional and directional convergence corrections. Snapshot 20261030T085136Z

**Status: development subset, preliminary.** Written 2026-10-30 15:41 CET (`date`). v2.1 answers the verifier report
`reviews/verify/verify_partA_v2.md` (commit 7327c7a) and PI request v2.1 items 1–6. **v2 (`../v2/`) is not modified.**
- Code commits: `b2bb8f7` (all v2.1 stages) and `0c966a5` (strata-bootstrap filename and functional cell-count fixes,
  before the reruns of stages 12, 15 and 11).
- Unchanged and copied from v2 (clean commits): universe, kNN, null, significance, classification of the primary
  τ = 0.2 run, CSLS.
- prereg sha256 `52f94618…0d83`. Every number comes from `numbers_partA.json` (= `results/numbers.json` → `partA`) or
  from a CSV here.
- Tests: `tests/test_partA.py` and `tests/test_descent.py`. In the faiss env, 26 pass. Under
  `uv run --project . pytest …`, 23 pass and 3 are skipped with a message because they need torch (`pytest_partA.log`;
  docs/ENV.md).

## 1. Headline estimates on the test half and with τ = 0.3 (`headline_halves_tau.csv`)
- **Split.** prereg tune/test: the key is the smallest accession among root-HOG members (all forest genes of the
  root; for proxy and singleton groups, the group's proteins), assigned by sha1 % 2. Control-list families are forced
  into the test half.
- **Half subsets.** A pair belongs to a half when both of its groups are in that half, so cross-half pairs appear in
  neither.
- **τ = 0.3.** SOGs come from `results/sogs/…/O1_r2/sensitivity_tau0.3`. The same-site split, r_c, F_conv,func and
  F_conv,dir are recomputed. F_anc and F_analog do not depend on τ except through the s/x split inside F_anc.
- **Uncertainty.** All CIs are 95 % two-way bootstrap over root groups: 200 replicates for the primary run on all
  pairs, 100 replicates for the halves, the strata and τ = 0.3.

| metric | τ 0.2, all | τ 0.2, **test half** | τ 0.2, tune half | τ 0.3, all | τ 0.3, test half |
|---|---|---|---|---|---|
| F_analog, per family pair | **0.819 [0.810, 0.826]** | 0.819 [0.805, 0.830] | 0.811 [0.798, 0.823] | 0.819 [0.810, 0.826] | 0.819 [0.805, 0.830] |
| F_anc, per family pair | 0.176 [0.169, 0.185] | 0.175 [0.163, 0.189] | 0.184 [0.172, 0.198] | 0.176 [0.169, 0.184] | 0.175 [0.163, 0.188] |
| F_analog, per pair | 0.367 [0.323, 0.416] | 0.285 [0.245, 0.340] | 0.219 [0.173, 0.281] | 0.367 [0.321, 0.413] | 0.285 [0.245, 0.340] |
| F_anc, per pair | 0.615 [0.566, 0.661] | 0.695 [0.636, 0.737] | 0.762 [0.705, 0.805] | 0.601 [0.555, 0.642] | 0.679 [0.621, 0.725] |
| r_adj H1s | 0.0100 [0.0086, 0.0119] | 0.0107 [0.0088, 0.0137] | 0.0094 [0.0080, 0.0112] | 0.0139 [0.0121, 0.0164] | 0.0149 [0.0121, 0.0191] |
| r_adj H2s | 0.0057 [0.0021, 0.0091] | 0.0019 [0.0013, 0.0036] | 0.0087 [0.0027, 0.0117] | 0.0087 [0.0032, 0.0139] | 0.0029 [0.0019, 0.0054] |
| r_adj H12s | 0.022 [0.017, 0.029] | 0.019 [0.014, 0.028] | 0.025 [0.019, 0.034] | 0.030 [0.023, 0.038] | 0.026 [0.019, 0.038] |
| r_adj N | 7.0e-7 [6.8e-7, 7.2e-7] | 6.8e-7 [6.6e-7, 7.2e-7] | 7.2e-7 [6.8e-7, 7.5e-7] | 7.0e-7 | 6.8e-7 |
| F_conv,func, per family pair (§3) | 0.081 [0.007, 0.178] | 0.127 [0, 0.317] (n = 41) | not computed | 0.081 [0.004, 0.155] | 0.127 (n = 41) |
| F_conv,dir, per pair (§4) | 1.1e-6 [0, 3.9e-6]; 3 passes | **0 (0 of 632,569 testable)** | not computed | 8.5e-7 [0, 2.9e-6]; 2 passes | 0 (0 of 559,090) |
| descent M1, ortholog / paralog (§5) | 0.774 / 0.780 | 0.775 / 0.784 | 0.773 / 0.777 | 0.739 / 0.739 | 0.741 / 0.743 |

How to read the table:
- **No conclusion depends on the half or on τ.** The family-pair F values are identical within their CIs.
- The per-pair F values differ between halves because cross-half pairs are excluded. This is a definitional
  difference, not a disagreement.
- **τ = 0.3 changes rates, not conclusions.** It lowers the same-site fraction of the denominators (H1 0.120 → 0.086),
  so r_H1s rises to 0.014. The ancestry-to-N rate ratio stays about 14,000–20,000×.
- F_conv,dir is 0 on both test halves.

## 2. numbers.json keys (verifier item 3)
- `11_numbers.py` now uses `to_dict("records")`. F_conv,func appears under named keys:
  - `functional[*]["F_conv_func.per_family_pair"]`, `….ci_lo`, `….ci_hi`;
  - `delta.family_pair`, `delta.pair`, and so on.
- CIs are separate numeric columns. They are no longer strings.

## 3. Functional convergence: family-pair δ, label sensitivity (`functional_convergence.csv`)
- **Headline weighting (v2.1).** F_conv,func = per-family-pair F_analog × per-family-pair δ. Each root-group cell of
  the test set or of the matched control set gets weight 1. The v2 form (per-pair δ) is kept as
  `F_conv_func.per_family_pair_x_pair_delta_v2`, which gives 0.078.
- **Mechanism of the sharing calls.** The 192 test pairs fall into 180 family cells. Of the 30 pairs that share a class, the mechanism breaks down as:

  | mechanism | share |
  |---|---|
  | ECFP4-Tanimoto only (mostly nucleotide ↔ nucleotide) | 47 % |
  | common metal token only | 33 % |
  | common non-metal token | 20 % |

- **Label sensitivity** (all subsets; pocketcompare training proteins excluded; δ is per family pair):

| labels used | n test pairs (cells) | δ, per family pair [CI] | F_conv,func, per family pair [CI] | F_conv,func, per pair |
|---|---|---|---|---|
| full (primary) | 192 (180) | 0.100 [0.008, 0.217] | **0.081 [0.007, 0.178]** | 0.035 [0.002, 0.078] |
| no Tanimoto (common tokens only) | 192 (180) | 0.054 [−0.006, 0.124] | 0.044 [0, 0.100] | 0.020 |
| no metal tokens | 149 (138) | 0.081 [−0.017, 0.184] | 0.066 [0, 0.150] | 0.028 |
| tokens only, no metal | 149 (138) | 0.026 [−0.012, 0.076] | **0.021 [0, 0.062]** | 0.011 |
| pLDDT ≥ 90, full | 100 (96) | 0.076 [−0.017, 0.193] | 0.044 [0, 0.097] | 0.011 |
| pLDDT ≥ 90, tokens only, no metal | 79 (75) | −0.009 [−0.012, −0.007] | 0 | 0 |
| test half, full | 41 (37) | 0.156 [−0.053, 0.414] | 0.127 [0, 0.317] | 0.038 |

- **Reading.** The functional excess is carried by generic nucleotide-Tanimoto and metal matches.
  - With specific non-metal ligand identity only, it falls to 0.021 [0, 0.062].
  - With specific non-metal ligand identity in high-confidence models, it is 0.
  - So F_conv,func ≈ 0.08 is an upper bound on "functional convergence", and it consists mostly of
    nucleotide/metal-binding site sharing.
  - For reference, homologous significant pairs share a class at 0.79 (H1s), 0.74 (H2s) and 0.76 (H3j) with full
    labels, and 0.66, 0.61 and 0.60 with tokens only and no metal.
- **Per-class rows (metal, cofactor, small molecule, EC) are EXPLORATORY.** The category scaling `frac_c` was added
  after the first functional output had been seen. This is the ledger D18 correction, as the verifier established
  from the commit times.
- **Residual covariate imbalance after μ matching** (verifier): aromatic fraction 0.175 vs 0.152, polar 0.504 vs
  0.476. This is not corrected; it is listed as a limitation.
- Selection-shift numbers in §4 (11–18 % of selected tests above q95; mean shifts) are from the verifier's toy
  (`verify_partA_v2.md` §7) and were not recomputed here.

## 4. Directional test: F_conv,dir ≈ 0; 3 passes, all flagged artifact-like (`directional_*`, `convergence_candidates.parquet`)
**The test, the nulls and the cutoff are unchanged from v2** (PA15 / ledger). Two things were added:
- **Toy no-convergence worlds over a fixed seed list** (`directional_toy_seeds.csv`; seeds 0–49, the prereg seeds
  [0, 1, 2] first; production toy: 20 families × 32 species):

  | world | result |
  |---|---|
  | no convergence | false-pass rate **4.4e-5** (Wilson 2.6e-5–7.3e-5) over 320,000 tested pairs; 5 of 50 seeds have ≥ 1 false pass; prereg seeds 0–2: 0 |
  | planted, frac 0.5 | power 0.96 |
  | planted, frac 0.2 | power 0.50 |
  | planted, frac 0.1 | 0: the planted pairs never reach significance |

  The real-data rate (1.1e-6) is **below** the toy false-pass rate, so the passes are consistent with chance.
- **Artifact checks** on every listed candidate. The flag rule was declared in code (`b2bb8f7`) before these runs. A
  candidate is flagged if any of the following holds:
  - minimum site volume < 100 Å³;
  - volume ratio > 5;
  - pLDDT < 80, or a fraction below 70 that exceeds 0.2, on either side;
  - SOG footprint ≤ 6 columns;
  - fewer than 8 residues on either side;
  - burrow score < 0.1 on either side.

| pair (τ 0.2) | volumes (Å³) | pLDDT, a / b (frac < 70 on b) | footprint | residues | burrow | flags |
|---|---|---|---|---|---|---|
| Q24298:s0 ↔ E1JGJ8:s0 | 233 / 276 | 84 / 87 | 8 | 11 / **6** | 0.72 / **0.09** | few residues, low burrow |
| P46723:s3 ↔ C0NED7:s2 | **46** / 835 | 98 / **76** (0.38) | 12 | 18 / 16 | **0.04** / 0.10 | small volume, 18× ratio, low pLDDT, low burrow |
| Q22937:s1 ↔ Q6R5P0:s1 | **53** / 890 | 96 / **74** (0.42) | **5** | 12 / 50 | **0.04** / 0.21 | small volume, 17× ratio, low pLDDT, small footprint, low burrow |

- **Result: F_conv,dir ≈ 0.** All 3 passes are flagged artifact-like. 85 % of all listed candidates (the top 500 by
  p plus the random-null-only passes) are flagged.
- At τ = 0.3 there are 2 passes, and both are flagged. On the test halves there are 0 passes.
- **Limitation: the nulls are not conditioned on selection (Δ_sel > 0).** In the toy, selection shifts mean Δ from
  0.022 to 0.037 and the random null to 0.032. Under the null, about 11–18 % of selected tests exceed q95, against a
  nominal 5 %, so the q95 criterion is anti-conservative.
  - This works against a false negative. Despite the anti-conservative criterion, F_conv,dir is still about 0.
  - Power was shown only for strong planted effects (frac ≥ 0.2). Weak lineage-level convergence would not be
    detected. Not changed here, per PI instruction.

## 5. Descent M1: ortholog vs paralog at matched identity AND divergence time (`descent_2d_summary.csv`, `descent_2d_strata.csv`)
**Method.**
- Inputs: the descent pair tables from `scripts/partA_descent/build_descent_pairs.py`, unmodified. The O1 τ 0.2 run is
  the verified DESCENT.md run. O1_r2 was rebuilt here at τ 0.2 and at τ 0.3 (34 s each).
- Primary variant qm90, capped sampling.
- Divergence time from TimeTree (`analyze_descent.divergence_table`); same-species pairs get 0.
- Standardised ortholog − paralog difference over cells where both classes have ≥ 20 informative observations; 200 OG
  bootstrap replicates.

| control | O1 τ 0.2, all | test half | tune half | O1_r2 τ 0.2 | O1_r2 τ 0.3 |
|---|---|---|---|---|---|
| identity only (existing) | −0.76 pp [−1.34, −0.12] | −0.94 [−1.87, −0.05] | −0.58 [−1.45, +0.14] | −0.72 [−1.21, −0.07] | −0.15 [−0.78, +0.56] |
| time only | −0.02 pp [−0.83, +0.75] | −0.30 [−1.41, +0.84] | +0.28 [−0.80, +1.44] | +0.05 [−0.72, +0.79] | +1.42 [+0.66, +2.28] |
| **identity × time (2-D)** | **−3.55 pp [−4.23, −2.79]** | −3.77 [−4.62, −2.62] | −3.28 [−4.30, −2.37] | −3.45 [−4.07, −2.67] | −3.02 [−3.66, −2.18] |
| identity × time, cross-species | identical to the 2-D row: the same-species bin holds paralogs only, so it drops out of the matched cells | | | | |

Reading:
- **The ortholog conjecture is not supported for pocket presence.**
- At matched identity and matched species-pair divergence time, ortholog pairs keep the pocket position no more often
  than paralog pairs. The point estimate is lower (−3.5 pp).
- This should **not** be read as "orthologs are less conserved". For paralogs, "time" is the divergence of the two
  species that carry them, not the age of the duplication, which is unknown. 30 % of paralog pairs are within one
  species, so time-matched contrasts are ill-posed, as DESCENT.md and the verifier note.
- The time-only contrast is about 0 at τ 0.2 and +1.4 pp at τ 0.3.
- Cell-level M1 is in `descent_2d_strata.csv` (identity decile × time bin × relation, with n). Example: identity
  0.3–0.4 at 100–300 Mya gives ortholog 0.698 vs paralog 0.786 (n 32k / 38k).

## 6. Tests and environment (verifier item 12)
- In the torch-dependent tests, `pytest.importorskip("torch", reason="torch not installed in this env (GPU path
  optional; docs/ENV.md)")` is used consistently.
- `test_covmodel_recovers_size_effect` is split so that its numpy part always runs. The block-form check moved to
  `test_covmodel_block_form_equals_pair_form`, which is skipped without torch.
- The toy world moved to `pocketevo.directional.toy_world`, shared by the tests and stage 13b.
- `docs/ENV.md` documents both environments.

## Deviations (PA ids; ledger mapping PA1 = D12, PA3 = D13; PA15–PA16 are ledger D17–D18)
- **PA19 Family-pair δ.** F_conv,func now uses the family-pair δ (headline 0.0815, was 0.0778). This changes a
  headline estimate and is **sent to the PI for the ledger**. Label-variant rows are sensitivity analyses; category
  rows are exploratory (D18 correction).
- **PA20 Artifact flag rule.** Declared in code before the v2.1 runs (§4), descriptive only. The test is unchanged.
- **PA21 Half membership of pairs.** A pair belongs to a half only if both of its root groups do. For non-HOG groups,
  the key is the smallest accession of the group's proteins.
- **PA22 Bulk-file rename.** The strata bootstrap file of the τ 0.3 run was written as
  `bootstrap_strata_cos_t0.3_t0.3.parquet` (filename bug, fixed in `0c966a5`) and renamed in bulk. Its contents are
  unchanged.

## Runtime (2026-10-30, GPU 3, 64-core node shared)
| stage | wall |
|---|---|
| 07 estimands, τ 0.2 (200 primary + 8 strata × 100 replicates) | 2,107 s |
| 06 classify, τ 0.3 | 350 s |
| 07 estimands, τ 0.3 (100 + 8 × 100 replicates) | 1,811 s |
| 12 functional, ×2 | 2 × about 130 s |
| 13 directional, ×2 | 240 s and 250 s |
| 13b toy, 80 worlds on 16 processes | 16 s |
| 16 descent 2-D (3 datasets × 3 halves × 200 replicates) | 64 s |
| descent pair rebuilds | 2 × 34 s |
| 14, 08, 10, 15, 11 | < 1 min in total |
| **total elapsed** | 14:35–15:41 (66 min) |
