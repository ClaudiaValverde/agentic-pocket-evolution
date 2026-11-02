# Part A: similarity decomposition (ancestry, analogy, convergence). Snapshot 20261030T085136Z

**Status: development subset, preliminary.** Tiers 1–2 only. Orthology is O1 (46 proteomes). The same-site split
uses a placeholder because SOGs are not ready. Producer: Part A lead analyst agent. These results have not been
verified (docs/VERIFY_CHECKLIST.md is still to run).
prereg sha256 `52f94618…0d83` (checked by every stage). Every number below comes from `numbers_partA.json`
(also `results/numbers.json` → `partA`) or from a CSV in this directory.

## What was computed
Universe: **3,221,702 sites** (1,779,029 tier 1 and 1,442,673 tier 2) on 861,985 proteins. Each site is its rep
fingerprint. The pair unit is an unordered site pair on two different proteins; there are 5.19 × 10^12 pairs.

| stage | script | what it does | wall |
|---|---|---|---|
| 1 | `01_universe.py` | concatenates the site metadata and fp16 fingerprints | 46 s |
| 2 | `02_knn.py` | exact inner-product kNN, FAISS GPU flat with fp32 storage, k = 1000, one RTX 4090 | **422 s** (371 s search) |
| 3 | `03_context.py` | root groups, region labels (Pfam clan and TED CATH-H at ≥ 50 % of site residues), PANTHER, covariates | 181 s |
| 4 | `04_null.py` | covariate regression (10^6 pairs), per-query GPD, held-out calibration (10^6), extended calibration (10^9) | 877 s (calibration on 10^9 pairs: 845 s) |
| 4b | `04b_tail_audit.py` | searches the held-out pairs in the deep tail | 10 s |
| 5 | `05_sig.py` | E-values for every kNN entry, both directions, then symmetrised | 74 s |
| 6 | `06_classify.py` | joins, MMseqs2 and Foldseek on the residual pairs and the 10^5 sample, same-site split | 479 s |
| 7 | `07_estimands.py` | r_c, the F estimands, 200-replicate two-way bootstrap, sensitivity grid, strata | 1,069 s |
| 8 | `08_topk_baselines.py` | top-k composition, size-only and covariate-only baselines | 30 s |
| 9 | `09_ligands.py` | BioLiP2 shared-ligand enrichment | 60 s |
| 10/11 | `10_figures.py`, `11_numbers.py` | calibration figure, numbers.json | <1 min |

Toy worlds ran before the real data: `tests/test_partA.py`, 12 tests, all pass (`pytest_partA.log`). They cover:
- LCA ortholog/paralog calls and exact forest counts against brute force;
- exact H3j/rest totals against enumeration of all pairs, with and without bootstrap weights and site masks;
- GPD tail calibration in exponential and bounded worlds;
- recovery of the CovModel coefficients;
- recovery of r_c and F in a synthetic pair world with planted analog pairs;
- the no-convergence world: F_analog is above 0 raw and about 0 after FP subtraction;
- the weighting schemes.

## Calibration (required figure: `fig_null_calibration_cos.png`)
- **Covariate model.** Mean and log-scale are symmetric quadratic in both pockets' covariates: log site volume,
  burrow score, lining polar fraction, lining aromatic fraction, pLDDT and model source. The fit used 10^6 random
  non-homologous pairs (seed 0). Covariates explain **R² = 0.19**; size alone explains **R² = 0.027**. The
  residual z has sd 0.97.
- **Tail.** The pooled shape is **ξ = −0.221** (per-query ξ: 5–95 % range −0.35 to −0.11). That is a bounded tail,
  as expected for a cosine.
- **Held-out 10^6 pairs (seed 1), one direction.** Observed/expected is 0.99 at p = 0.1, 0.98 at 1e-2 and 1.01 at
  1e-3, so **the null is calibrated in the bulk**. The ratio rises to 1.72 at 1e-4, 6.8 at 1e-5 and 38 at 1e-6.
- **Extended 10^9 pairs (seed 2).** At the significance threshold (E = 0.01, p = 3.1e-9), random non-homologous
  pairs pass in both directions at a rate of **9.9 × 10⁻⁷**. That is about 190 times the nominal 5.3 × 10⁻⁹.
- **The deep-tail excess is not hidden homology.** The 4b audit sent all ≈307 held-out pairs with p < 1e-4 in either direction, plus
  3,000 random pairs, to MMseqs2 and Foldseek, and **100 % of the tail pairs resolved to N** (`null_cos_tail_audit.csv`).
  The excess is therefore real pocket similarity between unrelated proteins: a heavier upper tail in the N population
  than the GPD extrapolated from its own bulk. It cannot be told apart from GPD misfit. This is the Gao–Skolnick
  degeneracy the design predicted, and it is the same thing F_analog measures.
- **Variants.** The per-query shape (the literal prereg reading) is worse in the tail: 13,600× at E = 0.01. With the
  shape clipped at ≥ 0, the test is over-conservative: **0 pairs reach E < 0.01**, which shows that significance at
  this threshold depends entirely on the bounded-tail extrapolation.
- **Cross-check.** Random-pair rate × n_N ≈ 5.1 × 10^6 expected significant N pairs. The kNN lists give K_N = 3.6 × 10^6.
  The two are consistent given k = 1000 (0.42 % of queries saturated, which is below the 1 % flag).

## Headline (pooled tail, E < 0.01 in both directions, expected FP nominal × λ = 1.72; 200 two-way bootstrap reps)
9,776,291 significant site pairs. Expected false positives: q0 = 5.34 × 10⁻⁹ per pair, **27,276 in N** and about
28,000 in total, against K_N = 3,605,586.

| weighting | F_anc | F_offsite (H1x/H2x/H12x) | F_H4 | F_analog (N) |
|---|---|---|---|---|
| **per family pair (headline)** | **0.176 [0.169, 0.185]** | 0.0003 [0.0003, 0.0003] | 0.005 [0.004, 0.006] | **0.819 [0.810, 0.826]** |
| per site | 0.280 [0.271, 0.291] | 0.016 [0.015, 0.018] | 0.004 [0.003, 0.006] | 0.700 [0.689, 0.710] |
| per pair (descriptive) | 0.601 [0.552, 0.646] | 0.025 [0.020, 0.033] | 0.007 [0.004, 0.012] | 0.367 [0.323, 0.416] |

Per-class rates r_c = (K_c − E0_c)/n_c. Brackets are 95 % CIs; E0 is the expected FP count.

| class | K | n_c | E0 | r_adj |
|---|---|---|---|---|
| H1s ortholog, same site | 141,883 | 4.5e6 | 0.02 | **0.031** [0.027, 0.037] |
| H2s paralog, same site | 456,864 | 1.9e7 | 0.1 | 0.024 [0.008, 0.040] |
| H12s same group, proxy-placed, same site | 3,823,450 | 1.2e8 | 0.7 | 0.031 [0.024, 0.041] |
| H1x / H2x / H12x offsite | 9,855 / 26,984 / 210,112 | 1.2e8 / 6.2e8 / 1.4e9 | 0.6 / 3.3 / 7.3 | 8.5e-5 / 4.4e-5 / 1.5e-4 |
| H3j joins (Pfam clan / TED SF / PANTHER fam) | 1,393,868 | 7.0e10 | 375 | 2.0e-5 [1.7e-5, 2.4e-5] |
| H3s MMseqs2 E < 1e-3 | 44,138 | 2.0e9 | 11 | 2.2e-5 |
| H4 Foldseek only | 63,551 | 6.2e9 | 33 | 1.0e-5 [5e-6, 1.9e-5] |
| **N** | 3,605,586 | 5.1e12 | 27,276 | **7.0e-7** [6.8e-7, 7.2e-7] |
| unknown | 0 | 0 | 0 | none (all proteins had a model and a sequence) |

How to read it:
- **Pockets descend by ancestry when they are positional homologs.** A same-site ortholog pair is significantly
  similar at a rate about 45,000 times higher than an unrelated pair (0.031 vs 7 × 10⁻⁷).
- **Within-family offsite pairs** (homologous proteins, different sites) are no more ancestry-like than H3/H4 pairs.
- **Most significant pairs involve no detectable homology.** This holds per site and per family pair, because the
  N class has about 70 times more pairs than all homologous classes together (5.1e12 vs 7.2e10).
- **F_analog = "similarity without detectable homology".** It is an upper bound on convergence, not an estimate of
  it: the directional test was not run, and the N tail looks like pocket-space degeneracy (see Calibration).
- **The F_anc vs F_analog split depends on the weighting.** This is the redundancy effect the review anticipated
  (Swiss-Prot and model-organism paralogs).

**Sensitivity** (`sensitivity_grid.csv`; per family pair / per pair F_analog):
- E < 1e-4: 0.800 / 0.330.
- E < 1: 0.796 / 0.381.
- Per-query-shape tail: 0.921 / 0.638 (25.4 M significant pairs).
- CSLS score (complete pipeline rerun with s′ = 2s − r_i − r_j): 0.905 / 0.581.
- The review's literal FP definition, Ê0 = observed rate among random non-homologous pairs × n_c, gives
  **F_analog = 0** in every cell. That is true by construction, because the random pairs are the N population; see
  PA3.

**Strata** (`strata_cos.csv`, both sites in the stratum; F_analog per family pair / per pair):
- **pLDDT ≥ 90: 0.588 / 0.158.** Moving to high-confidence models roughly halves the analog share at the family
  level and cuts it by more than half at the pair level, so a large part of the N signal comes from lower-confidence
  models.
- Both tier 1 (v6 proteomes, lower mean pLDDT): 0.886 / 0.514.
- Both tier 2 (Swiss-Prot): 0.390 / 0.091.
- Site volume tertiles, small / mid / large: 0.653 / 0.671 / 0.753 per family pair. Larger pockets are more
  analog-heavy.

**Raw top-k composition** (`topk_composition.csv`, 100k seeded queries, join classes only):

| k | H1 | H2 | H12 | H3j | rest (H3s/H4/N, not searched) |
|---|---|---|---|---|---|
| 10 | 2.4 % | 2.1 % | 22.4 % | 9.0 % | 64.2 % |
| 50 | 0.9 % | 1.1 % | 15.5 % | 7.9 % | 74.5 % |
| 200 | 0.4 % | 0.7 % | 9.7 % | 6.7 % | 82.6 % |

22.5 % of top-10 entries pass E < 0.01 in one direction. **Baselines.** Ranking neighbours by similarity predicted
from covariates alone gives top-10 entries that are 97.9 % "rest" with 0.04 % H1/H2/H12. Size alone gives 98.4 %
"rest", and random pairs 98.6 %. Covariates and size therefore carry essentially none of the homology signal.

## Functional validation (`ligand_enrichment.csv`)
- **Test set.** BioLiP2 regular ligands were mapped to sites (≥ 50 % of the BioLiP residues inside the site), giving
  19,081 annotated sites. Proteins in pocketcompare's training and validation data were excluded (1,943
  accessions). That leaves 139 significant N pairs with both sites annotated.
- **Result.** **13.7 % share a ligand class** (same CCD, or ECFP4 Tanimoto ≥ 0.5; Wilson CI 8.9–20.4 %). The control
  is 26,399 random annotated N pairs (MMseqs2/Foldseek-negative) matched on deciles of covariate-predicted
  similarity: 5.7 % share a class. **Enrichment is 2.4× (bootstrap CI 1.4–3.5).**
- **Without the training exclusion:** 155 pairs, 2.5× [1.6, 3.6].
- **PPV-type references in the annotated subset:** significant H1s 80 %, H2s 73 %, H3j 77 %.
- **Interpretation.** Significant analog pairs are modestly enriched for shared chemistry, far less than homologs.
  The enrichment is reported uncorrected.

## Deviations and choices not fixed by prereg (conservative option, logged 2026-10-30)
Mapping to the project ledger (DESIGN §12): PA1 (pooled tail) = ledger D12; PA3 (expected FP) = ledger D13.
**Superseded by v2/PARTA.md.**
- **PA1 Tail shape.**
  - Primary = pooled GPD shape (median of per-query PWM ξ) with per-query threshold and scale. Declared after the
    toy test (`test_gpd_tail_calibration`) showed that per-query shape noise makes extrapolated tails too light, and
    before any real-data result.
  - The prereg-literal per-query shape and the ξ ≥ 0 variant are reported throughout.
  - The declared fallback rule ("use the conservative variant if the primary fails calibration") would select ξ ≥ 0,
    which gives zero significant pairs at E < 0.01. I report that as the finding rather than as an estimand, because
    it measures the threshold's dependence on extrapolation and not the decomposition.
- **PA2 Background.** One seeded set of 12,000 random sites, with homologs masked per query (same protein, root
  group, or any shared Pfam clan/family, TED-H or PANTHER family label), replaces per-query draws. 981 queries
  (0.03 %) have fewer than 10,000 usable background sites (minimum 9,345).
- **PA3 Expected false positives.**
  - Headline: q0 = λ·E_thr/N_db, where λ = max(1, held-out observed/expected at the deepest threshold with ≥ 100
    expected) = 1.72, measured at p = 1e-4.
  - The review's definition (observed held-out null rate × n_c) is circular here. The held-out pairs are 99.9 % N
    by search, so that definition subtracts the whole N signal (F_analog ≡ 0). It is reported as the "bound" rows.
  - Under the weighted schemes, each expected FP pair is subtracted with weight 1 (its maximum).
- **PA4 Same-site split: PLACEHOLDER.**
  - Method: MMseqs2 local alignment of the two proteins, site residues mapped through it, mapped-residue Jaccard ≥ τ
    = 0.3 (the untuned prereg default). Applied to significant H1/H2/H12 pairs.
  - The denominator split comes from a 20,000-pair seeded within-group sample: same-site fractions H1 3.7 %, H2 3.0 %,
    H12 8.3 %.
  - Recompute with `06_classify.py --sogs <site_id→sog_id parquet>` when `results/sogs/…` is final; the hook is
    implemented. Until then, the H1s/H1x split and r_H1s are provisional.
- **PA5 Tier 2 / ID-gap proteins.** No O2 placement table exists yet. Each such protein takes the root HOG holding
  most O1 proteins of its PANTHER subfamily (363,756 proteins), else a PANTHER-subfamily pseudo-group (21,884), else
  a singleton (107,568). Within-group pairs involving a proxy-placed protein form the class H12 (flagged): no
  ortholog/paralog call. `03_context.py --placement` will use O2 when it arrives.
- **PA6 H4 TM rule.** The rule uses max(qTM, tTM) of the full-chain Foldseek alignment, not the region domain. It is
  more lenient, so more pairs become H4, which is conservative for N. The E-values of the MMseqs2/Foldseek pair
  alignments are relative to databases of all 861,985 universe sequences / 1.10 M AFDB v6 F1 models. The sequences
  are the model sequences.
- **PA7 Tiered searches.**
  - The primary set (pooled E < 0.01) is searched exhaustively: 3.71 M residual site pairs.
  - The other kept pairs (E < 1 under any variant; 24.9 M) used a seeded 1.2 % sample (299,533) with
    Horvitz–Thompson weights. This affects only the E = 1 and per-query-shape sensitivity cells.
- **PA8 Covariates.**
  - 44,740 sites (1.4 %) have imputed polarity/aromatic values because the UniProt sequence was missing (7,028
    proteins) or its length differed from the model (3,658 proteins). These sites are excluded from the
    fit/background/calibration pools.
  - "Model source" = tier (proteome v6 tar vs Swiss-Prot v6 tar); both are v6.
  - Polar = DEHKNQRSTY, aromatic = FWY.
- **PA9 Bootstrap.**
  - Poisson(1) group multiplicities, with seeds seed·1000 + r.
  - Weights: cross-group pairs w_g·w_h, within-group pairs w_g. Both are unbiased; w_g² would double the
    within-group mass.
  - Exact totals and the denominator sample are reweighted jointly in each replicate.
  - Per-site and per-family base weights are not recomputed inside replicates.
  - 200 replicates (amendment 3), not the 1,000 in the prereg `uncertainty` block.
- **PA10 Added, not in prereg.**
  - Extended 10^9-pair calibration (seed 2) and the 4b tail audit.
  - CSLS r_i = mean cosine of the 10 nearest neighbours, homologs included. The CSLS extended calibration uses
    2 × 10^8 pairs.
- **PA11 Functional control.** Drawn separately (30,000 random pairs of annotated sites, seed 20, resolved by the same
  searches), because the 10^5 denominator sample holds too few doubly annotated pairs. Matching is on deciles of
  μ_ij.
- **PA12 Labels.**
  - Pfam region = envelope (seq_start–seq_end); label = clan, or the family if unclustered.
  - TED labels are used only at CATH H level (superfamily); 1.37 M sites have no TED-H region label.
  - PANTHER family is protein-level (prereg).
- **PA13 Not done** (open):
  - HHsearch/TM-align undetected-homology correction on 200 N pairs.
  - The directional convergence test F_conv,dir, which needs SOGs.
  - Inverse-clade-size reweighting.
  - Stratum CIs (point estimates only).

## Open issues
1. Replace the PA4 placeholder with SOGs and the PA5 proxy with O2 placement, then rerun stages 6–7 and 11. Stage
   1–5 outputs stay valid.
2. Undetected-homology correction for N (HHsearch on 200 pairs). F_conv,dir.
3. The deep-tail excess of the N null is real non-homologous similarity. E < 0.01 relies on bounded-tail
   extrapolation, so treat F_analog as "significantly similar beyond the covariate bulk", not as convergence.
4. The pLDDT ≥ 90 stratum shows that F_analog depends on model quality. The full atlas (tier 3, v4 Foldcomp) needs
   the model-source covariate re-fit.
5. Minor: the `rest_sample_accept` field in cls_*_summary.json is computed wrongly (it reports n/drawn for one
   oversized draw). Random pairs are 98.4 % "rest" (`topk_composition.csv`, random row).

## Files
Bulk (node-local, regenerable): `/tmp/pe_big/partA/20261030T085136Z/` holds:
- `universe.parquet`, `fp.f16.npy`;
- `knn/{I.i32,D.f16}.npy` (19 GB);
- `null_*/`, `sig_*/pairs.parquet`, `cls_*/pairs_cls.parquet`, `sample_rest.parquet`, `sample_within.parquet`;
- `bootstrap_cos.parquet`.

Foldseek DB of all AFDB v6 F1 models: `/tmp/pe_big/partA/foldseek/afdb_v6`. Faiss/torch env:
`/tmp/pe_big/envs/faiss`.

Rerun: `cd scripts/partA && for s in 01_universe 02_knn 03_context 04_null 04b_tail_audit 05_sig 06_classify
07_estimands 08_topk_baselines 09_ligands 10_figures 11_numbers; do CUDA_VISIBLE_DEVICES=3
/tmp/pe_big/envs/faiss/bin/python $s.py; done`. CSLS is run with `--score csls` for 04/05/06/07. Stages are
resumable or cached and parameterized by `--snapshot/--sites/--orthology/--annotations`.
