# provenance checks 1243 (2026-10-30T11:43:04Z)

Read-only run of the three check scripts on the current repo. No outputs modified, no commit. REPORT.md does not exist yet; check_report_numbers was also run on results/partA/20261030T085136Z/PARTA.md as a demonstration.

## check_prereg.py
```
| check | result | detail |
|---|---|---|
| prereg.yaml hash | PASS | 52f946189098... |
| json outputs scanned | 997 | ok=387 differ=0 missing=610 |
| ledger ids | 8 | C1, D10, D11, D12, D13, D7, D8, D9 |
| referenced ids | 13 | unknown=['D1', 'D2', 'D3', 'D4', 'D5', 'D6'] |

MISSING prereg hash (610), first 15:
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/ca2_zn_azm_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-1_sample-0/ca2_zn_azm_seed-1_sample-0_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-1_sample-1/ca2_zn_azm_seed-1_sample-1_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-1_sample-2/ca2_zn_azm_seed-1_sample-2_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-1_sample-3/ca2_zn_azm_seed-1_sample-3_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-1_sample-4/ca2_zn_azm_seed-1_sample-4_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-2_sample-0/ca2_zn_azm_seed-2_sample-0_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-2_sample-1/ca2_zn_azm_seed-2_sample-1_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-2_sample-2/ca2_zn_azm_seed-2_sample-2_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-2_sample-3/ca2_zn_azm_seed-2_sample-3_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-2_sample-4/ca2_zn_azm_seed-2_sample-4_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-3_sample-0/ca2_zn_azm_seed-3_sample-0_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-3_sample-1/ca2_zn_azm_seed-3_sample-1_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-3_sample-2/ca2_zn_azm_seed-3_sample-2_summary_confidences.json
   results/cofold/bench/ca2_zn_azm/ca2_zn_azm/seed-3_sample-3/ca2_zn_azm_seed-3_sample-3_summary_confidences.json

unknown ledger id D3: results/partA/20261030T085136Z/PARTA.md:100, results/partA/20261030T085136Z/PARTA.md:147, results/partA/20261030T085136Z/v1/PARTA.md:100, results/partA/20261030T085136Z/v1/PARTA.

unknown ledger id D1: results/partA/20261030T085136Z/PARTA.md:136, results/partA/20261030T085136Z/v1/PARTA.md:136

unknown ledger id D2: results/partA/20261030T085136Z/PARTA.md:144, results/partA/20261030T085136Z/v1/PARTA.md:144

unknown ledger id D4: results/partA/20261030T085136Z/PARTA.md:153, results/partA/20261030T085136Z/PARTA.md:203, results/partA/20261030T085136Z/v1/PARTA.md:153, results/partA/20261030T085136Z/v1/PARTA.

unknown ledger id D5: results/partA/20261030T085136Z/PARTA.md:160, results/partA/20261030T085136Z/PARTA.md:203, results/partA/20261030T085136Z/v1/PARTA.md:160, results/partA/20261030T085136Z/v1/PARTA.

unknown ledger id D6: results/partA/20261030T085136Z/PARTA.md:164, results/partA/20261030T085136Z/v1/PARTA.md:164

RESULT: FAIL
 - 0 json outputs with differing and 610 with missing prereg hash
 - deviation ids referenced but not in ledger: ['D1', 'D2', 'D3', 'D4', 'D5', 'D6']
exit code: 1
```

## check_provenance.py
```
| check | count | detail |
|---|---|---|
| distinct sidecar commits | 4 | missing from git: 0 |
| sidecars w/o commit field | 8 | |
| sidecars citing unknown commit | 0 | |
| dirty-tree sidecars | 374 | commits: ['0e112016ea', 'bab66c35d5', 'c325bfaa7a', 'cbe924346b'] |
| tracked results modified vs HEAD | 0 | |
| untracked results paths | 12 | (info) |
| manifest results/sites/20261030T085136Z/checksums.sha256 | ok | verified=0 bad=0; bulk-missing (node-local): 646 |
| manifest results/sogs/20261030T085136Z/O1/BULK_MANIFEST.sha256 | ok | verified=0 bad=0; bulk-missing (node-local): 10563 |
| manifest results/sogs/20261030T085136Z/O1/checksums.sha256 | ok | verified=545 bad=0 |

DIRTY (374), first 12:
   results/partA/20261030T085136Z/estimands_cos.csv.json bab66c35d5
   results/partA/20261030T085136Z/estimands_csls.csv.json bab66c35d5
   results/partA/20261030T085136Z/fig_null_calibration_cos.png.json 0e112016ea
   results/partA/20261030T085136Z/ligand_enrichment.csv.json bab66c35d5
   results/partA/20261030T085136Z/null_cos_tail_audit.csv.json 0e112016ea
   results/partA/20261030T085136Z/null_csls_summary.json 0e112016ea
   results/partA/20261030T085136Z/numbers_partA.prov.json bab66c35d5
   results/partA/20261030T085136Z/topk_composition.csv.json bab66c35d5
   results/partB/events/20261030T085136Z/O1/tau0.2/primary/candidates.parquet.json bab66c35d5
   results/partB/events/20261030T085136Z/O1/tau0.2/primary/events.parquet.json bab66c35d5
   results/partB/events/20261030T085136Z/O1/tau0.2/primary/funnel.csv.json bab66c35d5
   results/partB/events/20261030T085136Z/O1/tau0.2/primary/manifest.json bab66c35d5

NO COMMIT FIELD (8), first 12:
   results/numbers.json
   results/partA/20261030T085136Z/numbers_partA.json
   results/partA/20261030T085136Z/v1/numbers_partA.json
   results/partB/sim/measured/primary_known_gt10/thresholds.json
   results/partB/sim/measured/primary_known_le10/thresholds.json
   results/partB/sim/measured/primary_redetect_gt10/thresholds.json
   results/partB/sim/measured/primary_redetect_le10/thresholds.json
   results/partB/sim/thresholds.json

UNTRACKED (info) (12), first 12:
   cases/
   results/cofold/controls/
   results/cofold/controls_ca2/
   results/cofold/controls_cdk2/
   results/cofold/controls_pygl/
   results/cofold/controls_yadh1/
   results/cofold/msatest/
   results/cofold/msatest_cdk2_cacheunpaired/
   results/cofold/msatest_cdk2_cf/
   results/partA/20261030T085136Z/v1/
   results/partB/events/
   results/partB/sim/measured/

RESULT: FAIL
 - 374 sidecars were produced from a dirty tree
exit code: 1
```

## check_report_numbers.py
```
FAIL: report <repo>/REPORT.md not found
exit code: 2
```

## check_report_numbers.py results/partA/20261030T085136Z/PARTA.md
```
results/partA/20261030T085136Z/PARTA.md: 210 notable numbers, 743 values in numbers.json, 56 unmatched
  line 10: 1,779,029    | Universe: **3,221,702 sites** (1,779,029 tier 1 and 1,442,673 tier 2) on 861,985 proteins. Each site
  line 10: 1,442,673    | Universe: **3,221,702 sites** (1,779,029 tier 1 and 1,442,673 tier 2) on 861,985 proteins. Each site
  line 10: 861,985    | Universe: **3,221,702 sites** (1,779,029 tier 1 and 1,442,673 tier 2) on 861,985 proteins. Each site
  line 11: 5.19    | fingerprint. The pair unit is an unordered site pair on two different proteins; there are 5.19 × 10^
  line 11: 12    | fingerprint. The pair unit is an unordered site pair on two different proteins; there are 5.19 × 10^
  line 16: 4090    | | 2 | `02_knn.py` | exact inner-product kNN, FAISS GPU flat with fp32 storage, k = 1000, one RTX 409
  line 16: 422    | | 2 | `02_knn.py` | exact inner-product kNN, FAISS GPU flat with fp32 storage, k = 1000, one RTX 409
  line 16: 371    | | 2 | `02_knn.py` | exact inner-product kNN, FAISS GPU flat with fp32 storage, k = 1000, one RTX 409
  line 17: 181    | | 3 | `03_context.py` | root groups, region labels (Pfam clan and TED CATH-H at ≥ 50 % of site resid
  line 18: 877    | | 4 | `04_null.py` | covariate regression (10^6 pairs), per-query GPD, held-out calibration (10^6), 
  line 18: 845    | | 4 | `04_null.py` | covariate regression (10^6 pairs), per-query GPD, held-out calibration (10^6), 
  line 20: 74    | | 5 | `05_sig.py` | E-values for every kNN entry, both directions, then symmetrised | 74 s |
  line 21: 479    | | 6 | `06_classify.py` | joins, MMseqs2 and Foldseek on the residual pairs and the 10^5 sample, same
  line 22: 1,069    | | 7 | `07_estimands.py` | r_c, the F estimands, 200-replicate two-way bootstrap, sensitivity grid, s
  line 27: 12    | Toy worlds ran before the real data: `tests/test_partA.py`, 12 tests, all pass (`pytest_partA.log`).
  line 40: 0.97    | residual z has sd 0.97.
  line 41: −0.221    | - **Tail.** The pooled shape is **ξ = −0.221** (per-query ξ: 5–95 % range −0.35 to −0.11). That is a
  line 41: 95 %    | - **Tail.** The pooled shape is **ξ = −0.221** (per-query ξ: 5–95 % range −0.35 to −0.11). That is a
  line 46: 9.9    | pairs pass in both directions at a rate of **9.9 × 10⁻⁷**. That is about 190 times the nominal 5.3 ×
  line 46: 190    | pairs pass in both directions at a rate of **9.9 × 10⁻⁷**. That is about 190 times the nominal 5.3 ×
  line 46: 5.3    | pairs pass in both directions at a rate of **9.9 × 10⁻⁷**. That is about 190 times the nominal 5.3 ×
  line 47: 307    | - **The deep-tail excess is not hidden homology.** The 4b audit sent all ≈307 held-out pairs with p 
  line 52: 13,600    | - **Variants.** The per-query shape (the literal prereg reading) is worse in the tail: 13,600× at E 
  line 55: 3.6    | - **Cross-check.** Random-pair rate × n_N ≈ 5.1 × 10^6 expected significant N pairs. The kNN lists g
  line 59: 5.34    | 9,776,291 significant site pairs. Expected false positives: q0 = 5.34 × 10⁻⁹ per pair, **27,276 in N
  line 60: 28,000    | 28,000 in total, against K_N = 3,605,586.
  line 68: 95 %    | Per-class rates r_c = (K_c − E0_c)/n_c. Brackets are 95 % CIs; E0 is the expected FP count.
  line 84: 45,000    | similar at a rate about 45,000 times higher than an unrelated pair (0.031 vs 7 × 10⁻⁷).
  line 94: 0.330    | - E < 1e-4: 0.800 / 0.330.
  line 95: 0.381    | - E < 1: 0.796 / 0.381.
  line 96: 0.921    | - Per-query-shape tail: 0.921 / 0.638 (25.4 M significant pairs).
  line 96: 0.638    | - Per-query-shape tail: 0.921 / 0.638 (25.4 M significant pairs).
  line 96: 25.4    | - Per-query-shape tail: 0.921 / 0.638 (25.4 M significant pairs).
  line 125: 19,081    | 19,081 annotated sites. Proteins in pocketcompare's training and validation data were excluded (1,94
  line 125: 1,943    | 19,081 annotated sites. Proteins in pocketcompare's training and validation data were excluded (1,94
  line 127: 8.9    | - **Result.** **13.7 % share a ligand class** (same CCD, or ECFP4 Tanimoto ≥ 0.5; Wilson CI 8.9–20.4
  line 130: 155    | - **Without the training exclusion:** 155 pairs, 2.5× [1.6, 3.6].
  line 130: 2.5    | - **Without the training exclusion:** 155 pairs, 2.5× [1.6, 3.6].
  line 130: 1.6    | - **Without the training exclusion:** 155 pairs, 2.5× [1.6, 3.6].
  line 130: 3.6    | - **Without the training exclusion:** 155 pairs, 2.5× [1.6, 3.6].
  line 144: 12,000    | - **D2 Background.** One seeded set of 12,000 random sites, with homologs masked per query (same pro
  line 145: 981    | group, or any shared Pfam clan/family, TED-H or PANTHER family label), replaces per-query draws. 981
  line 146: 9,345    | (0.03 %) have fewer than 10,000 usable background sites (minimum 9,345).
  line 150: 99.9 %    | - The review's definition (observed held-out null rate × n_c) is circular here. The held-out pairs a
  line 161: 363,756    | most O1 proteins of its PANTHER subfamily (363,756 proteins), else a PANTHER-subfamily pseudo-group 
  line 161: 21,884    | most O1 proteins of its PANTHER subfamily (363,756 proteins), else a PANTHER-subfamily pseudo-group 
  line 162: 107,568    | a singleton (107,568). Within-group pairs involving a proxy-placed protein form the class H12 (flagg
  line 166: 861,985    | alignments are relative to databases of all 861,985 universe sequences / 1.10 M AFDB v6 F1 models. T
  line 166: 1.10    | alignments are relative to databases of all 861,985 universe sequences / 1.10 M AFDB v6 F1 models. T
  line 169: 3.71    | - The primary set (pooled E < 0.01) is searched exhaustively: 3.71 M residual site pairs.
  line 170: 24.9    | - The other kept pairs (E < 1 under any variant; 24.9 M) used a seeded 1.2 % sample (299,533) with
  line 170: 299,533    | - The other kept pairs (E < 1 under any variant; 24.9 M) used a seeded 1.2 % sample (299,533) with
  line 173: 44,740    | - 44,740 sites (1.4 %) have imputed polarity/aromatic values because the UniProt sequence was missin
  line 173: 7,028    | - 44,740 sites (1.4 %) have imputed polarity/aromatic values because the UniProt sequence was missin
  line 174: 3,658    | proteins) or its length differed from the model (3,658 proteins). These sites are excluded from the
  line 194: 1.37    | - TED labels are used only at CATH H level (superfamily); 1.37 M sites have no TED-H region label.
RESULT: FAIL
exit code: 1
```
