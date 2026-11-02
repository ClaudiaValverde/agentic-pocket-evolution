# REPORT v0: pocket evolution across the AlphaFold DB pocket atlas

**Status: final, written 2026-10-31 21:03 CET (`date`), after the 09:45 deadline (see §4 item 15). Development subset, preliminary.** Snapshot
`20261030T085136Z`, prereg sha256 `52f94618…0d83` (DESIGN.md v2, frozen at G1). Branch `pocket-evolution`.
The G5 second-snapshot comparison (§5) was completed last. Every number below comes from
[results/numbers.json](../results/numbers.json) or from the committed JSON/CSV linked next to it
(`scripts/check_report_numbers.py` checks this). Numbers we need but do not yet have in a machine-readable output are
listed under [TODO numbers](#todo-numbers), not typed in.

**What to push:** branch `pocket-evolution` (commits only; no push credential in the container).

## 0. Summary for reviewers

1. **Ancestry vs analogy (Part A).** Pocket similarity, counted as significant site pairs, is mostly not explained by
   detectable homology when each pair of families gets weight 1. The **preregistered headline** (per family pair) is
   F_analog = 0.819 [0.810, 0.826] on orthology O1. That number depends on how finely the orthology cuts families: it
   is 0.935 on O2 (D29). The per-pair estimate, which agrees between O1 and O2 on all pairs (post hoc framing, D29), is F_analog = 0.367 [0.323, 0.416] on O1 and
   0.367 [0.325, 0.416] on O2. F_analog is an upper bound on convergence, not convergence.
2. **Convergence is weak.** Functional excess among non-homologous similar pairs is F_conv,func = 0.081 [0.007, 0.178]
   (per family pair), and 0.021 [0, 0.062] when only specific non-metal ligands count. **Directional convergence ≈ 0**
   (3 passes among 2,627,453 testable pairs, all artifact-like; real rate below the toy false-pass rate). Non-homologous
   similarity is consistent with pocket-space degeneracy, with at most a small functional excess.
3. **Descent.** Pocket positions are retained in most informative ortholog pairs and decay with divergence, but the
   pocket "ortholog conjecture" is **not supported** for presence: orthologs keep the site no more often than paralogs
   at matched identity, on both orthologies.
4. **Part B events: 0 convincing cases in every category.** Gains pass a simulated FDR of 0.05, but the **empirical
   null FDR is higher** (curveball null 0.199 on O1 on 175 pre-rescue calls, 0.267 on O2 on rescued states, not directly comparable;
   [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json)). **Losses are not
   FDR-calibrated** (D21): null inputs give as many loss calls as real data. **Specificity changes: 0** pass BH
   (preregistered), but BH cannot fire at the null's resolution, and at the funnel stage real sister pairs exceed the
   pseudo-clade null (513 vs 71–88 on O1, 536 vs 66–84 on O2; about 1.4–2× attributable to the fingerprint shift once the residue gate is passed;
   empirical FDR ≈ 0.15; D34, post hoc, 3 null seeds): a population-level signal without nameable cases. Most
   gain candidates fail the artifact cascade, mainly on orthology (on O1, outgroups that O2 places in a different root HOG, 79.5 %; on O2,
   comparisons against another paralog subfamily of the same HOG, D33) and present-side fragments.
5. **Known events: recall 0/7.** The annotation-first validation recovers none of the 7 literature events. The causes
   are mixed (C-D31): some events lie between orthogroups, which a within-orthogroup pipeline cannot see (MR/GR on O1,
   IDH, GID1/HSL); the others are within one orthogroup but lack confident absence calls, show the pocket in both
   states, or have no significant specificity test. A cross-orthogroup extension also recovered 0/7 (D32). The
   validation negative controls (uricase, GULO, M. leprae truncations) produce no spurious losses; the separate
   false-loss co-fold controls are in item 6.
6. **Co-folding** is calibrated (AF3; ligand pLDDT AUC 0.79 [0.68, 0.89],
   [threshold.json](../results/cofold/calibration/threshold.json)); the D23 composite is specific (0.96) but insensitive
   (0.42), so a "binds" call is informative and a "does not bind" call is weak. **On cases:** no gain or loss is
   co-fold-consistent with a cognate ligand, because in the 3 cases that pass the positive-control gate the cognate
   ligand is placed 12–34 Å (per-member medians) from the gained footprint, in a different pocket (BURROW rank
   1–3), and no cognate job passes D23. Ligand-free refolds reproduce 3
   of 4 testable gains (the footprint pocket appears only on the present side), but confirm only 1 of 4 testable
   false-loss controls. Known-answer deep dive: AF3 resolves the MR/GR steroid specificity difference in 3/3 seeds
   (small effect); LDH/MDH is uninterpretable (gate failed). Scoring was **unblinded** (§3.8).

## 1. Approach

The pipeline runs one command per stage on a frozen atlas snapshot manifest, with every threshold in
`config/prereg.yaml` and every output carrying `{snapshot_id, prereg_sha256, git_commit, seed}`. An agent that did not
produce a result verified each headline before use (reviews/verify/).

**Data and the development subset** (results/coverage/20261030T085136Z/COVERAGE.md). Dataset sources, versions,
checksums and licences: [DATA.md](../DATA.md).
- Tier 1: 46 AFDB v6 reference proteomes, 555,846 structures, 1,788,750 kept pockets. Tier 2: Swiss-Prot, 436,886
  structures. Tier 3 (AFDB v4 via Foldcomp): 7,841 work units done and 72,463 pending, 22,962,643 structures
  ([tier_summary.csv](../results/coverage/20261030T085136Z/tier_summary.csv),
  [tier3_done_vs_pending_units.csv](../results/coverage/20261030T085136Z/tier3_done_vs_pending_units.csv)).
- Tier L: pockets computed here for 39 eukaryotic proteomes (637,503 structures, 1,837,831 pockets) so that
  persistence across species is testable near human, fly, nematodes and plants.

**Orthology.** O1 = Ancestra HOGs on the 46 tier-1 proteomes. O2 = O1 plus tier L and context prokaryotes/protists:
157 species, 1,549,554 genes, 57,332 root groups
([O2 manifest](../data/orthology/ancestra_O2/accurate/manifest.json)). Tier-2 proteins are placed into O2 HOGs by DIAMOND
(placed fraction 0.897, [placement.json](../data/orthology/ancestra_O2/placement/tier2/placement.json)). Check against
PANTHER: 0.158 of same-subfamily human–mouse pairs fall in different O2 root groups
([validation.json](../data/orthology/ancestra_O2/validation/validation.json)).

**Sites and site orthogroups (SOGs).** A site is a fragment-deduplicated BURROW pocket (3,221,702 sites in tiers 1–2).
Per orthogroup with ≥ 3 species, sites are mapped onto a FAMSA MSA; two sites are positional homologs if column Jaccard
≥ τ; SOGs are Leiden communities. τ = 0.2 was tuned on the tune half of a benchmark of experimentally evidenced
ortholog site pairs (D10: grid edge, τ = 0.3 reported as sensitivity). O1: 260,991 SOGs; O2: 463,014 SOGs
([O1_r3 funnel](../results/sogs/20261030T085136Z/O1_r3/funnel.json), [O2_r1 funnel](../results/sogs/20261030T085136Z/O2_r1/funnel.json)).

**Observation states.** Each SOG × member is P (present), W (weak present), A (absent observed: footprint pLDDT ≥ 90,
aligned, length-normal, no signal peptide/TM, footprint > 10 residues by D11), U (uninformative) or G (no gene). A calls
are rescued (BURROW with relaxed settings and other seeds; k = 1 of 6 settings turns A into W, D9). Error rates were
measured on known answers before any event was called (§3.1).

**Part A.** All tier-1/2 site pairs above a significance threshold (GPU kNN k = 1000; covariate-regression null plus a
pooled generalised-Pareto tail, E < 0.01 in both directions, D12), split into homology classes (ortholog/paralog
same-site, sequence homolog, remote structural, none = N), with expected false positives subtracted. Estimands:
F_anc, F_H4, F_analog under three weightings (per family pair = preregistered headline; per pair; per site), the
functional convergence excess F_conv,func, and a cross-fitted directional convergence test F_conv,dir. Descent: M1
(site retained in the other member, given informative), M2 (fingerprint similarity), M3 (tree consistency vs curveball
null).

**Part B.** Sankoff parsimony per SOG on the HOG gene tree (gain cost 2, loss 1; 1:1 and 3:1 as sensitivity).
Thresholds were set to 10 % FDR on a simulated no-event null that uses the measured error rates. Every candidate goes
through null-input pipelines (curveball-permuted states, random pseudo-clades; 20 seeds each), then an artifact
cascade (footprint pLDDT, disorder, TED/PAE proxy, fragment on either side, orthology via O2 and PANTHER, a DIAMOND
reciprocal-best-hit check, species-tree lineages), then three evidence filters (persistence, experimental,
computational). Specificity changes use a residue-first test (≥ 2 clade-consistent BLOSUM62 < 0 substitutions plus a
fingerprint shift above the noise floor, empirical matched-depth null, BH). An annotation-first validation track traces
known events and negative controls through the same outputs.

**Co-folding.** AlphaFold3 inference on 4× RTX-4090, cached AF3 MSAs where the sequence matches, ColabFold MSAs
otherwise, `--no-templates` primary. Readouts were frozen on a BioLiP2 calibration panel before any case was folded
(D23). Details in docs/COFOLD.md.

**Web app.** FastAPI plus vendored 3Dmol.js. It reads index files and streams structures from the AFDB tars and the
Foldcomp DB in place (no AFDB copy). Views: protein → pockets in 3D; orthogroup → SOG × species state matrix (O1_r2,
O2_r1); cases with evidence and filter status (webapp/README.md). The state matrix shows O1_r2; Part B results use O1_r3 (same SOGs). The
tier-3 protein index lives in /tmp and must be rebuilt (`python -m webapp build-index`) after a node restart.

### Runtime of the main stages

| stage | wall time | source |
|---|---|---|
| Orthology O1: families / scores / merge stages | 266 s / 407 s / 520 s | [manifest](../data/orthology/ancestra_tier1_v0/accurate/manifest.json) |
| Orthology O2: families / scores / merge stages | 1,587 s / 2,290 s / 1,965 s | [manifest](../data/orthology/ancestra_O2/accurate/manifest.json) |
| Tier-2 DIAMOND placement into O2 | 188.5 s | [placement.json](../data/orthology/ancestra_O2/placement/tier2/placement.json) |
| Tier L pockets (39 proteomes) | 5,246.6 s | [numbers.json](../results/numbers.json) `observation.tierL.run_wall_s` |
| SOG states + export, O1_r3 / O2_r1 | 71.9 s / 144.4 s | [O1_r3](../results/sogs/20261030T085136Z/O1_r3/funnel.json), [O2_r1](../results/sogs/20261030T085136Z/O2_r1/funnel.json) |
| Error rates: re-detection / false-A on known sites | 28.2 s / 1,032.1 s | [redetection.json](../results/observation/20261030T085136Z/redetection.json), [false_a.json](../results/observation/20261030T085136Z/false_a.json) |
| Part A directional test / descent analysis | 225.8 s / 354.8 s | [numbers.json](../results/numbers.json) |
| Validation track (annotation-first) | 26.9 s | [validation_table.parquet.json](../results/cases/20261030T085136Z/validation_table.parquet.json) |
| AF3 job, 3 seeds: CDK2 + ATP + Mg / PYGL + PLP | 79.5 s / 312.2 s | [runtime.tsv](../results/cofold/bench/runtime.tsv) |
| Co-fold calibration panel: median job / total GPU | 236 s / 6,051.7 GPU-s (24 jobs) | [threshold.json](../results/cofold/calibration/threshold.json) |

Part B wall times: O1 events 51.2 s, rescue 61.7 s, nulls 244.5 s, artifact cascade 17.2 s, specificity 1,876.3 s; O2
events 101.6 s, nulls 379.4 s, cascade 20.1 s; DIAMOND RBH 293.6 s ([numbers.json](../results/numbers.json)
`partB.wall_times`). Part A end-to-end wall time was not exported.

## 2. Best methods (what we chose and why)

- **Absence calls only where the model is good and the footprint is large.** A requires footprint pLDDT ≥ 90 and
  > 10 footprint residues (D11), then rescue k = 1 (D9). This is the most conservative option on the grid; no k met the
  preregistered 5 % false-A target (§3.1).
- **Pooled GPD tail for Part A significance (D12).** Per-query tail shapes were too light-tailed on a toy with known
  answer; the pooled tail was chosen before real results. Per-query results are reported alongside.
- **Per-pair weighting as the robust Part A statement (D29).** The preregistered per-family-pair headline is kept, but
  it moves with orthology granularity; per pair agrees on all pairs (0.367 vs 0.367) but not on the
  test halves (0.285 vs 0.136, different splits).
- **O2 for case finding (D28).** O1 families are over-merged: most O1 gains rest on outgroup absences that O2 puts in a
  different family (§3.6).
- **Symmetric quality gates in descent (D15, D29b).** The same pLDDT and footprint-size rules apply to the P and A side;
  the asymmetric version inflated M1 on O2.
- **Co-folding concordance rule (D23).** A state difference counts only if ligand pLDDT and a geometry readout
  (centroid distance or site recall) agree in ≥ 2 of 3 seeds. ipTM is not used.

## 3. Test-set results

### 3.1 Observation error rates (known answers)

| quantity | value | source |
|---|---|---|
| Same-protein site re-detection, v6 atlas vs Foldcomp v4 re-description | 0.653 (0.716 at pLDDT ≥ 90) | [numbers.json](../results/numbers.json) `observation.redetection` |
| Re-detection across BURROW seeds | 0.807 | [numbers.json](../results/numbers.json) |
| False-A on 73,107 known sites (BioLiP2, UniProt ECO, LIGYSIS), primary A rule | 0.399 | [numbers.json](../results/numbers.json) `observation.false_A` |
| False-A by ligand class: metal / small molecule / cofactor | 0.609 / 0.356 / 0.225 | [numbers.json](../results/numbers.json) |
| False-A after rescue k = 1 (tune / test half) | 0.129 / 0.142 | [numbers.json](../results/numbers.json) `observation.rescue` |
| False-A after rescue for sites > 10 residues (tune half) | 0.054 | [rescue_breakdown.csv](../results/observation/20261030T085136Z/rescue_breakdown.csv) |
| True absences wrongly rescued to W (k = 1, tune half) | 0.348 | [rescue_breakdown.csv](../results/observation/20261030T085136Z/rescue_breakdown.csv) |

Reading: even a high-quality AFDB model misses a known site often, so a single absence is weak evidence. This is why
losses cannot reach the FDR target (D21).

### 3.2 Site matching across orthologs (test half)

| orthology | τ | F1 [95 % CI] | recall (both sites detected) / precision | source |
|---|---|---|---|---|
| O1 | 0.2 | 0.681 [0.604, 0.714] | 0.591 / 0.803 | [tau_benchmark.csv](../results/sogs/20261030T085136Z/O1/tau_benchmark.csv) |
| O2, leakage-clean split (D20) | 0.2 | 0.695 [0.618, 0.775] | 0.633 / 0.771 | [tau_benchmark_o2_split.csv](../results/sogs/20261030T085136Z/O2_r1_split/tau_benchmark_o2_split.csv) |
| O2, leakage-clean split (D20) | 0.3 | 0.699 [0.638, 0.774] | 0.603 / 0.831 | [tau_benchmark_o2_split.csv](../results/sogs/20261030T085136Z/O2_r1_split/tau_benchmark_o2_split.csv) |

### 3.3 Part A: what fraction of pocket similarity is ancestry? (D29 framing)

The first two value columns are the **preregistered headline** (per family pair); the last two are the
**granularity-robust per-pair** estimate. 95 % CIs from a two-way bootstrap over root groups. N significant pairs:
9,776,291 on O1 ([numbers.json](../results/numbers.json)).

| run | subset | F_anc, per family pair | F_analog, per family pair (prereg) | F_anc, per pair | F_analog, per pair (robust) | source |
|---|---|---|---|---|---|---|
| O1 | all | 0.176 [0.169, 0.185] | **0.819 [0.810, 0.826]** | 0.615 [0.566, 0.661] | **0.367 [0.323, 0.416]** | [numbers.json](../results/numbers.json) `partA` |
| O2 | all | 0.060 [0.057, 0.063] | **0.935 [0.932, 0.937]** | 0.616 [0.567, 0.658] | **0.367 [0.325, 0.416]** | [numbers.json](../results/numbers.json) `partA_O2` |
| O1 | test half | 0.175 [0.163, 0.189] | 0.819 [0.805, 0.830] | 0.695 [0.636, 0.737] | 0.285 [0.245, 0.340] | [numbers.json](../results/numbers.json) `partA.strata.test_half` |
| O2 | test half (D20) | 0.083 [0.077, 0.091] | 0.909 [0.900, 0.916] | 0.849 [0.802, 0.881] | 0.136 [0.104, 0.181] | [numbers.json](../results/numbers.json) `partA_O2.strata.test_half` |
| O1 | both pLDDT ≥ 90 | 0.401 [0.386, 0.418] | 0.588 [0.571, 0.604] | 0.820 [0.794, 0.840] | 0.158 [0.134, 0.181] | [numbers.json](../results/numbers.json) `partA.strata` |
| O2 | both pLDDT ≥ 90 | 0.116 [0.109, 0.123] | 0.870 [0.862, 0.878] | 0.821 [0.796, 0.844] | 0.158 [0.135, 0.178] | [numbers.json](../results/numbers.json) `partA_O2.strata` |

- Per-class significance rates are weighting-free: same-site orthologs r_H1s = 0.010 vs non-homologous r_N = 7.0e-7 on
  O1 ([numbers.json](../results/numbers.json) `partA.per_class`).
- Hidden-homology correction (Foldseek check on an N sample): F_analog per family pair 0.815 on O1.
- Why O1 and O2 differ per family pair: O2 groups are coarser for all proteins, so homolog pairs fall into fewer
  family-pair cells (the verifier also sees the shift within tier-1-only pairs). Only the cell partition changes:
  1,529,928 of 9,757,299 pairs significant in both runs change class, mostly between homologous classes
  ([class_moves_summary.json](../results/partA/20261030T085136Z/v3/class_moves_summary.json)).
- On the per-pair halves, cross-half pairs drop out, so half rows are not directly comparable with the "all" rows.
  τ = 0.3 changes no conclusion (results/partA/20261030T085136Z/v2_1/PARTA.md §1).

### 3.4 Convergence

| estimand | value | source |
|---|---|---|
| F_conv,func, per family pair, full labels (n = 192 significant N pairs) | 0.081 [0.007, 0.178]; shared ligand class 2.62× matched random | [numbers.json](../results/numbers.json) `partA.functional` |
| F_conv,func, specific non-metal ligand tokens only (n = 149) | 0.021 [0, 0.062] | [numbers.json](../results/numbers.json) `partA.functional` |
| F_conv,dir, per pair (O1, all) | 1.1e-6 [0, 3.9e-6]; 3 passes of 2,627,453 testable | [numbers.json](../results/numbers.json) `partA.directional` |
| F_conv,dir, test half | 0 of 632,569 testable | [numbers.json](../results/numbers.json) |
| Toy no-convergence world: false-pass rate | 4.4e-5 (power 0.96 for strong planted effects) | [numbers.json](../results/numbers.json) `partA.directional_toy_seeds` |

The real directional pass rate is below the toy false-pass rate, and all 3 passes carry artifact flags (small volume,
low pLDDT, low BURROW score). The functional excess is carried by nucleotide (Tanimoto) and metal matches; with specific
ligands in pLDDT ≥ 90 models it is 0. Convergence positives are detection-limited: catalytic residues of the
chymotrypsin and carbonic-anhydrase controls are mostly outside BURROW pockets (validation table).

### 3.5 Descent

| quantity | ortholog | paralog | ortholog − paralog at matched identity [CI] | source |
|---|---|---|---|---|
| M1, O1, prereg-literal qm90 (all) | 0.774 | 0.780 | −0.0076 [−0.0134, −0.0012] | [numbers.json](../results/numbers.json) `partA.descent_2d_summary` |
| M1, O1, symmetric D11 (test half) | 0.889 | 0.900 | −0.0110 [−0.0176, −0.0048] | [descent_O1_qm90d11_summary.csv](../results/partA/20261030T085136Z/v3_1/descent_O1_qm90d11_summary.csv) |
| M1, O2, symmetric D11 (test half) | 0.885 | 0.888 | −0.0022 [−0.0083, +0.0051] | [descent_O2_qm90d11_summary.csv](../results/partA/20261030T085136Z/v3_1/descent_O2_qm90d11_summary.csv) |
| M1 at matched identity × time, O1 | | | −0.0355 [−0.0423, −0.0279] | [numbers.json](../results/numbers.json) |

- Retention is relative to a re-detection ceiling below 1 (0.869 for > 10-residue sites at pLDDT ≥ 90,
  [numbers.json](../results/numbers.json) `partA_descent.ceiling`).
- Ortholog conjecture for presence: **not supported** on either orthology. This is not a claim that orthologs are less
  conserved: paralog "time" is the species split, not the (unknown) duplication age.
- M3 tree consistency: 0.933 vs curveball null 0.919 [0.918, 0.919], p = 0.001 over 1,000 replicates, and 0.942 under
  simulated noisy-vertical evolution ([numbers.json](../results/numbers.json) `partA_descent.M3`). Presence tracks the tree
  beyond row/column margins, modestly.

### 3.6 Part B: event candidates, nulls and the artifact cascade

| orthology | type | called (O1 pre-rescue; O2 on rescued states) | empirical null FDR: curveball / pseudo-clade | source |
|---|---|---|---|---|
| O1 | gain | 175 | 0.199 / 0.080 (simulated FDR 0.051) | [real_vs_null_same_cut.csv](../results/partB/20261030T085136Z/O1_r3/real_vs_null_same_cut.csv) |
| O1 | loss | 61 | 1.000 / 0.407 (not FDR-calibrated) | [real_vs_null_same_cut.csv](../results/partB/20261030T085136Z/O1_r3/real_vs_null_same_cut.csv) |
| O2 | gain | 104 | 0.267 / 0.093 | [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) |
| O2 | loss | 66 | 1.000 / 0.674 (not FDR-calibrated) | [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) |

| orthology | type | primary calls | pass artifact cascade (strict) | paralog-branch flag | source |
|---|---|---|---|---|---|
| O1 (O1_r3, post-rescue) | gain | 50 | 7 (4) | 32 | [artifact_summary.json](../results/partB/20261030T085136Z/O1_r3/artifact_summary.json) |
| O1 | loss | 17 | 1 (0) | 11 | [artifact_summary.json](../results/partB/20261030T085136Z/O1_r3/artifact_summary.json) |
| O2 (O2_r1 + rescue) | gain | 104 | 17 (13) | 57 | [artifact_summary.json](../results/partB/20261030T085136Z/O2_r1/artifact_summary.json) |
| O2 | loss | 66 | 6 (1) | 23 | [artifact_summary.json](../results/partB/20261030T085136Z/O2_r1/artifact_summary.json) |

- **Gains are statistically real in aggregate but carry a higher error than simulated.** The empirical curveball FDR
  exceeds the 10 % target (D24); no threshold was changed after seeing this.
- **Losses carry no statistical signal.** The curveball null calls as many losses as the real data. They go forward
  only at audit level and must stand on case evidence (D21).
- **O1 gains are mostly orthology artifacts.** Of 50 O1 gains, 2 are re-called on O2 (0 of the 32 paralog-branch
  gains); most outgroup A members behind them are not in the derived clade's O2 SOG
  ([o1_o2_gain_agreement.json](../results/partB/20261030T085136Z/O1_r3/o1_o2_gain_agreement.json); a root-HOG recount by
  the verifier gives 79.5 %, reviews/verify/verify_G4b.md). Where the RBH check fails, the outgroup's best hit is
  usually in the same root HOG but outside the derived clade, i.e. the call compares against another paralog
  subfamily (D33).
- **Specificity: 0 preregistered BH calls, uninformative by construction** (D34): the smallest attainable empirical p
  (3.5e-5 O1, 1.6e-5 O2) is above the single-call BH threshold (4.7e-6, 2.6e-6). O1: 222,489 branch pairs, 513 pass the
  residue gate + shift vs 71–88 under the pseudo-clade null (3 seeds); O2 (fingerprints read from atlas rows): 743,349
  pairs, 536 vs 66–84 ([funnel_real_vs_null.csv](../results/partB/20261030T085136Z/O1_r3/specificity/tau0.2/funnel_real_vs_null.csv)).
  Empirical FDR ≈ 0.15 at the funnel and 0.17–0.19 after the cascade (200 vs 30–44 on O1, 230 vs 32–45 on O2; 4 cascade
  checks are not evaluable on null pairs). Given the residue gate, the fingerprint shift adds only 1.4–2.0×. A
  population-level signal, post hoc; no individual pair is supported.
- **Duplication-associated gains** (exploratory, not preregistered): derived paralog copy P vs sister copy A in ≥ 3
  species. Some pass on O1 and O2, none clears the artifact cascade on O1, and the selection-conditioned null is
  near-empty by construction. Hypothesis-generating only (8 of 32 on O1, 18 of 55 on O2; unselected background with ≥ 3 species:
  0 of 2,597 and 4 of 5,276; [numbers.json](../results/numbers.json) `partB.duplication_gains`).
- **False-loss negative controls** (required by the prompt): losses called on pre-quality-filter states whose
  derived-clade A members are all poor models (low footprint pLDDT or ≤ 10-residue footprint) disappear once the
  primary A rule and D11 apply; 5 are carried as case controls (pools on O1_r3: 319 bad-model, 48 rescue-instability).
- **Robustness of O1 gains:** 98 % are also called at cost 1:1 and 84 % at τ = 0.3, but only **34 % at 3:1**
  ([numbers.json](../results/numbers.json) `partB.O1_primary_gains`).
- **RBH orthology:** reciprocal best hits of the derived clade among outgroup supporters: 0.095 (O1), 0.197 (O2);
  Pfam-architecture match 0.731 / 0.835 (`partB.rbh_orthology`).

### 3.7 Cases and the evidence filters

Headline rule: **convincing** = persistence + experimental + computational filters + artifact cascade. Counts from
[case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) `counts` (frozen at G4, D30):

| category | candidates | f1 persistence | f2 experimental | f3 computational | artifacts pass | supported | **convincing** | source |
|---|---|---|---|---|---|---|---|---|
| O1 gain | 50 | 42 | 1 | 47 | 7 | 39 | **0** | [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) |
| O1 loss | 17 | 8 | 2 | 16 | 1 | 8 | **0** | [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) |
| O2 gain | 104 | 76 | 0 | 61 | 17 | 47 | **0** | [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) |
| O2 loss | 66 | 19 | 3 | 28 | 6 | 7 | **0** | [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) |
| O1 specificity (exploratory top 5) | 5 | 2 | 0 | 0 | 0 | 0 | **0** | [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) |
| O1 convergence (directional passes) | 3 | 3 | 0 | 0 | 0 | 0 | **0** | [case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json) |

- **0 convincing cases**, also with the post hoc literature line for f2 (D25b). f2 fails because no experimental record
  lies at the footprint in the state that has the pocket; a verifier's own UniProt/PDBe/RCSB search found nothing
  missed (reviews/verify/verify_G4.md B5).
- **Co-fold pool** (rule-based, "supported" + artifacts pass, top-5 by S + 1 seeded audit per category): 6 gains, 1
  loss, 5 bad-model false-loss controls, 2 exploratory duplication gains. Verifier flags on pool members: VAC14
  present-side pLDDT below 90, citrate synthase thin conservation margin and possible hexamer interface, secernin
  indel/gene-model risk, connexin footprint in the TM bundle in contact with lipid.
- **Annotation-first validation** ([validation_table.parquet.json](../results/cases/20261030T085136Z/validation_table.parquet.json)):
  known events recovered **0 of 7** (steroid-receptor LBP, NR4A/Nurr1 loss, PfLDH vs MDH, IDH NADP→NAD, GID1 vs HSL,
  D14 vs KAI2, lysozyme vs α-lactalbumin) on both O1 and O2. Main reason: derived and ancestral proteins fall into
  different orthogroups for MR/GR (O1), IDH and GID1/HSL; the remaining items are within one orthogroup and are missed
  for lack of confident absences, a pocket detected in both states, or a non-significant specificity test (C-D31).
  Negative controls: uricase is G (no gene) in human and GULO is G on O1 (on O2 the human rows are P/U paralogs);
  M. leprae truncations produce no M. leprae-only loss calls.

### 3.8 Co-folding

Calibration panel: 24 BioLiP2 holo positives (6 each cofactor, drug-like, metal, nucleotide), 72 negatives (48
property-matched in-class decoys, 24 site ablations); job statistic = median over 3 seeds × 5 samples
([threshold.json](../results/cofold/calibration/threshold.json), frozen before any case was folded).

| readout | AUC [95 % CI] | source |
|---|---|---|
| ligand pLDDT | 0.79 [0.68, 0.89] | [threshold.json](../results/cofold/calibration/threshold.json) |
| min protein–ligand PAE | 0.76 [0.62, 0.88] | [threshold.json](../results/cofold/calibration/threshold.json) |
| site recall / centroid distance | 0.77 [0.69, 0.85] / 0.73 [0.62, 0.83] | [threshold.json](../results/cofold/calibration/threshold.json) |
| ligand pLDDT, metal ligands only | 0.58 [0.31, 0.84] | [threshold.json](../results/cofold/calibration/threshold.json) |
| D23 composite (pLDDT ≥ 92 AND (centroid ≤ 2 Å OR recall ≥ 0.9), ≥ 2/3 seeds) | sensitivity 0.42, specificity 0.96 | [threshold.json](../results/cofold/calibration/threshold.json) |

Metal sites cannot be called by co-folding (pLDDT at chance). Panel decoys are hard in-class near-misses, so these AUCs
are a lower bound relative to the easier decoys in docs/COFOLD.md §4.3. ipTM is recorded but not used (D23), because
decoys and site-destroyed variants keep high ipTM in the positive-control runs.

**Case results** ([numbers.json](../results/numbers.json) `cofold`; per-case readouts in cases/<id>/case.yaml;
total GPU time for cases and deep dives: TODO number, not yet exported). Scoring is **unblinded**: a post-processing bug re-minted the blinding salt, and the
state mapping was recovered by sequence matching before the final scoring pass; the folds themselves ran blind.

| case | test | result | verdict |
|---|---|---|---|
| citrate synthase gain (O1) | cognate OAA, gate 1.34 Å | ligand at a rank-1 pocket on both sides, not at the footprint | uninformative (ligand not at footprint) |
| MetE gain (O1) | cognate MET, gate 1.14 Å | ligand off-footprint on both sides; footprint pocket re-detected 3/3 present vs 0/3 absent | uninformative for ligand; gain reproduced in fresh models |
| choline kinase (exploratory duplication gain, O2) | cognate ADP, gate 0.37 Å | ligand off-footprint; footprint 2/3 present vs 0/3 absent | uninformative for ligand; gain reproduced |
| RNA pol II loss (O1) | cognate nucleotide | positive-control gate failed (7.52 Å) | uninterpretable |
| gain HOG0034532 (O2) | ligand-free refold, site recall present vs absent | 1.00 vs 0.00 | consistent with gain |
| gain HOG0056135 (O2, audit case) | ligand-free refold | 0.923 vs 0.077 | consistent with gain |
| gain HOG0026737 (O1) | ligand-free refold | 0.75 vs 0.00 | consistent with gain |
| gain HOG0025964 (O1) | ligand-free refold | present side 0.222 | uninformative |
| duplication gain HOG0046679 (O2) | ligand-free refold | 0.706 vs 0.353 | inconclusive |
| 5 false-loss controls | ligand-free refold of the bad-model absent side (prediction: pocket reappears) | reappears in 1 (0.60), not in 3 (0.25, 0.00, 0.00), 1 untestable | **prediction mostly refuted** |

Reading: the gained pockets that survive are *secondary sites*, reproducible in independent AF3 models but not where
the family's known ligand binds, so their function is unknown and co-folding cannot supply a ligand-level
confirmation. The false-loss controls were chosen as losses whose absent-side AFDB models are poor; for 3 of them a
fresh ligand-free AF3 model still lacks the pocket, so a fresh model does not restore the pocket. Whether these are real losses discarded by
the quality filter or errors shared across AF2/AF3 is unresolved; the control label is unsupported for them.

**Deep dives with known answers** (success criterion pre-committed in 21b5c9e before results):
- **MR vs GR** (aldosterone vs cortisol): gate passed (0.31 Å). All four cells bind in the site (centroid 0.5–1.5 Å,
  recall 1.00). The state × ligand interaction in ligand pLDDT is +2.96 / +2.68 / +2.80 over 3 seeds, the expected
  sign every time. This meets the pre-committed "resolved" criterion, but the effect (~3 pLDDT) is small; all four
  cells, including GR + aldosterone, pass D23, and receptor–steroid complexes are likely in AF3's training data, so
  this is weak evidence of ligand-level discrimination.
- **LDH vs MDH** (pyruvate vs oxaloacetate): interaction +18.9 / +9.9 / +17.5, the expected sign, but the
  positive-control gate failed (2.66 Å on symmetric malonate without symmetry-corrected RMSD), so it is reported
  **uninterpretable** and not overridden. An independent rerun (results/cofold/deep_dives/, from the parallel session's agent)
  reproduces the steroid numbers exactly (same inputs: reproducibility, not replication) and uses an LDHA holo gate that
  passes, but in that rerun all four PfLDH/MDH cells place the ligand about 10–12 Å from the site ([deep_dives.tsv](../results/cofold/deep_dives/deep_dives.tsv); none passes D23), and in
  the primary run LDH + pyruvate also fails D23, so neither run supports an LDH/MDH resolution; switching gates after seeing results
  would be post hoc, so the primary verdict stays "uninterpretable", and that rerun is noted only.

## 4. Limitations (stated plainly)

1. **0 convincing cases.** No gain, loss, specificity or convergence case passes all three filters. The report has no
   headline case.
2. **Recall 0/7 on known events**, with mixed causes (C-D31). MR/GR (O1), IDH and GID1/HSL lie between orthogroups,
   which the within-orthogroup pipeline cannot see; D14/KAI2, lysozyme/α-lactalbumin, Nurr1 and ER vs AR/GR (O1) and
   PfLDH/MDH are within one orthogroup but are missed for lack of confident absences, a pocket detected in both
   states, or a missing or non-significant specificity test. The cross-OG extension also recovers 0/7 (D32).
3. **Losses are not FDR-calibrated** (D21). Under the measured noise no loss threshold reaches 10 % simulated FDR with
   non-zero recall, and the curveball null produces as many losses as the real data.
4. **Gains: the empirical null FDR (0.199 O1, 0.267 O2, curveball) is well above the simulated 0.05**
   ([case_list_v2.parquet.json](../results/cases/20261030T085136Z/case_list_v2.parquet.json)). Gains are a population
   signal with a large false fraction, not individually trustworthy calls.
5. **Specificity: 0** BH-significant branch pairs on O1 and O2; BH cannot fire at the null's resolution (3 null
   seeds); the population-level excess (D34, post hoc) has no individually supported pair.
6. **Directional convergence ≈ 0.** The directional null is not conditioned on selection, which makes the q95 cut
   anti-conservative, so the null result is not caused by a strict cut; weak lineage-level convergence would still be
   missed. Convergence controls (catalytic triads, CA classes) are detection-limited.
7. **Part A headline depends on orthology granularity** (per family pair 0.819 on O1 vs 0.935 on O2); the per-pair
   estimate is robust ([numbers.json](../results/numbers.json)). F_analog is an upper bound on convergence. Residual
   covariate imbalance remains after matching in the functional test.
8. **Development-subset biases** (results/coverage/20261030T085136Z/COVERAGE.md): tier 3 has 7,841 work units done
   and 72,463 pending ([tier3_done_vs_pending_units.csv](../results/coverage/20261030T085136Z/tier3_done_vs_pending_units.csv));
   Archaea and unclassified/metagenomic taxa are over-represented (enrichment 9.3 vs AFDB), Bacteria (0.72) and
   Eukaryota (0.53) under-represented ([tier3_superkingdom.csv](../results/coverage/20261030T085136Z/tier3_superkingdom.csv));
   within a taxon the done slice is contiguous, not random; tier-1 eukaryote models have mean pLDDT 73.7
   ([quality_by_superkingdom.csv](../results/coverage/20261030T085136Z/quality_by_superkingdom.csv)), so many absences are
   U. Model source (v6 file vs Foldcomp v4) is a covariate. All findings are preliminary until rerun on the full atlas.
9. **Absence is weak evidence** even in good models (false-A 0.399 before rescue, 0.129 after; metal sites 0.609;
   [numbers.json](../results/numbers.json)).
10. **Orthology:** deep (Eukaryota-level) splits are uncertain; duplication levels are unreliable for near-identical
    families; 0.158 of same-PANTHER-subfamily human–mouse pairs are split across O2 root groups
    ([validation.json](../data/orthology/ancestra_O2/validation/validation.json)).
11. **Descent** ortholog−paralog contrasts are age-confounded (duplication age unknown).
12. **Co-folding** is relative evidence only: specificity 0.96 but sensitivity 0.42
    ([threshold.json](../results/cofold/calibration/threshold.json)); metal ligands uninterpretable; positive controls
    are in the AF3 training set.
13. **The 3D viewer was not browser-tested.** Headless Chromium lacks system libraries in the container and the
    alternative install route was denied; API endpoints and server-rendered pages are covered by tests. A human needs to
    open the viewer once.
14. **Process:** many post-freeze deviations (Appendix A), each with a stated bias direction; several outputs carry
    `+dirty` code hashes in their sidecars (reviews/verify/provenance_1243.md).
15. **Schedule and oversight.** No commits 18:07–22:56 on 30 Oct (co-folding stalled ~18:40–22:58), and the PI was idle
    from 23:58 (30 Oct) to 20:39 (31 Oct). The 09:45 delivery deadline was missed; this report was finalised afterwards.
16. **Co-fold case scoring is unblinded.** A post-processing bug re-minted the blinding salt, and the state mapping was
    recovered by sequence before the final scoring pass. The folds ran blind; the deep dives stayed blinded.
17. **Two PI sessions ran in parallel** after a client restart (~23:12 on 30 Oct) until the user stopped one. Its
    "PI review" commits 48f713e/4d54c6b were reverted; O2 specificity was run twice (the D34 rerun is bit-identical:
    reproducibility, not replication); the LDHA-gated deep-dive rerun also came from that session.
18. **Re-runnability on a new snapshot failed at first.** On snapshot 20261031T020009Z, partB_O1, partB_O2 and the final
    Part A stage failed on snapshot-specific inputs that were not carried over (rescue_k.json, the O1 events path, the
    TimeTree pairwise table), and the G5 Part A stage overwrote the snapshot-1 `partA` numbers (restored in 531f95f).
    See §5 for the fixed rerun.
19. **Post hoc choices.** The O2 RBH orthology rule was chosen after PANTHER/O2-root came back "unknown" for all gains
    (D33: 17 pass instead of 34); also post hoc: the D34 specificity population signal and tail fits, the D29 per-pair
    framing, the D25b literature f2 line, duplication gains, and the cross-OG extension (D31/D32). None enters the
    convincing headline.
20. **Unrescued absences.** Roughly a hundred thousand O2 A calls with no filtered pockets were never rescued; they make up 23 % of the O2
    gain outgroup-A supporters (verify_G4b A5).
21. **False-loss controls:** for 3 of 5 bad-model controls a fresh ligand-free AF3 model also lacks the pocket, so the
    A-quality filter may be discarding some real losses as well as artifacts.

## 5. Exploratory extension and second-snapshot comparison

- **Cross-orthogroup extension** (exploratory, D31/D32; method committed before any output in 77db16c; CROSSOG.md):
  structure-aligned site matching between OGs of the same TED/Pfam family recovers **0/7** validation items, so its
  ≥ 3/7 gate failed and no genome-wide run was made. Site linking itself is specific (610 linked SOG pairs vs 1 under a
  shuffled correspondence in a 2,000-pair sample;
  [O1_funnel.json](../results/partB_crossog/20261030T085136Z/characterization_sample2000/O1_funnel.json)), but the specificity test is too conservative under a global
  matched-identity null, and the presence test fires more under the shuffle (10 vs 33).
- **G5 second-snapshot comparison** (D27; snapshot 20261031T020009Z, taken when tier 3 was complete:
  80,304 work units, 689,760,778 pockets; [snapshot_manifest.json](../results/pipeline/G5/20261031T020009Z/snapshot_manifest.json),
  [G5_STATUS.md](../results/pipeline/G5_STATUS.md)). The first unattended run (02:00–05:20) failed in three stages on
  snapshot-specific inputs (§4 item 18). After fixes (b83541b, 97e031f: frozen inputs looked up under the prereg
  snapshot, the O1→O2 agreement step ordered after both lanes, plus three further defects found on the way) every stage
  finishes with rc 0 (partB_O1 65 min with BURROW recomputed from scratch, partB_O2 16 min, Part A 50 min).
  - **O1 (identical inputs): reproduces exactly.** SOGs are byte-identical; F_analog per pair / per family pair,
    F_conv,func, M1, and the 50 gains / 17 losses (7 / 1 after the cascade) are identical
    ([g5__numbers_partA.json](../results/pipeline/G5/20261031T020009Z/g5__numbers_partA.json), [O1 artifact summary](../results/pipeline/G5/20261031T020009Z/O1_g5__artifact_summary.json)).
    Residual nondeterminism: M1 bootstrap CIs in the 4th decimal (row order) and 2 RBH fractions (DIAMOND threading),
    with no change in any pass/fail ([reproducibility_notes.md](../results/pipeline/G5/20261031T020009Z/reproducibility_notes.md)).
  - **O2 (new tier-3 data): stable.** 463,318 SOGs vs 463,014 ([snapshot 2](../results/pipeline/G5/20261031T020009Z/O2_g5__funnel.json),
    [snapshot 1](../results/sogs/20261030T085136Z/O2_r1/funnel.json)); clean-test F1 at τ = 0.2 essentially unchanged
    ([O2 benchmark](../results/pipeline/G5/20261031T020009Z/O2_g5__tau_benchmark.csv)); O2 gains 110 after rescue, 17 after the artifact cascade (vs 104 → 17)
    ([O2 artifact summary](../results/pipeline/G5/20261031T020009Z/O2_g5__artifact_summary.json)). Most candidates recur, but not all, so individual O2
    candidates are snapshot-sensitive while the counts and the conclusions are not.
  - Reading: the pipeline now reruns unchanged on a new snapshot; Part A and the O1 results are reproducible, and the
    O2 results are stable at the level of counts and headline numbers. Gains of tier-3 coverage (10× more pockets)
    did not change any conclusion, though O1/O2 orthology (and therefore SOG membership) were not recomputed for the
    new tier-3 species (D27).

## 6. Follow-ups and remaining ideas

1. **Sister-paralog comparison instead of outgroup absence** for paralog-branch gains (possible neofunctionalisation),
   with an unselected background null; the exploratory duplication-gain test is the starting point.
2. **Cross-orthogroup event calling** (D31) so that MR/GR, LDH/MDH, IDH and GID1/HSL-type events become visible;
   validate on the 7 known events before any genome-wide run.
3. **≥ 20 pseudo-clade null seeds and a preregistered tail model for specificity**, so that individual specificity
   calls become attainable.
4. **Selection-conditioned null for the directional test**, plus power for weak lineage-level convergence.
5. **Measured false-P rate** (currently a placeholder in the Part B noise model) and a likelihood (Mk) observation
   model in place of simulation-calibrated parsimony.
6. **PAE from the AFDB API** instead of the TED-domain proxy; oligomeric assemblies for interface footprints.
7. **Ancestral sequence reconstruction** for the deep dives (IQ-TREE `-asr` plus posterior samples), cut for time.
8. **Full-atlas rerun** when tier 3 completes: same commands, new snapshot id (§7); report Part A/B stability.
9. Inverse-clade-size reweighting for Part A; HHsearch-level hidden-homology correction on a larger N sample.

## 7. Reproduce

Everything reads a snapshot manifest under `data/snapshot/<id>/` and `config/prereg.yaml`. Bulk intermediates go to
`/tmp/pe_big` (node-local, regenerable).

```bash
# whole pipeline on a new atlas snapshot (print commands first, then run)
.venv/bin/python scripts/run_pipeline.py --snapshot new --run-id g5 --dry-run
.venv/bin/python scripts/run_pipeline.py --snapshot new --run-id g5

# Part B alone, reproducing the committed results (O1 -> O1_r3, O2 -> O2_r1 by default)
bash scripts/partB/run_partB.sh --orthology O1
bash scripts/partB/run_partB.sh --orthology O2

# research browser (single command, from the repo root)
uv run --project . python -m webapp --port 8765      # open http://127.0.0.1:8765/

# tests and report checks
.venv/bin/python -m pytest -q tests/
.venv/bin/python scripts/check_report_numbers.py REPORT.md
```

## Numbers not exported

- Part A end-to-end wall time (kNN, null, significance, classification, bootstrap).
- Total GPU time of case co-folding and deep dives.
- Web app: largest-OG matrix latency and test count (reported in agent notes only).

## Appendix A. Deviation ledger summary (DESIGN.md §12; full text there)

Bias direction in the last column. D0–D6 and D14 are decisions made before or at the freeze (RESEARCH_TREE.md decision
log), not deviations.

| id | deviation (short) | bias |
|---|---|---|
| D7 | Loss cut raised and counted after persistence; superseded by C-D7 | conservative |
| D8 | SOGs with footprint ≤ 3 columns excluded from events | fewer events |
| D9 | Rescue k = 1 although no k meets the 5 % false-A target; residual propagated | fewer absences |
| D10 | τ = 0.2 at the grid edge; τ = 0.3 sensitivity for every headline | none; reported |
| D11 | A requires a footprint > 10 residues | fewer absences |
| D12 | Pooled GPD tail (per-query shape reported) | per-query gives higher F_analog |
| D13 | Expected-FP subtraction reported as upper and lower bound | bounds bracket truth |
| C1 | Timestamps corrected to commit times | none |
| D15 | M1 primary "qm90": P informative only if it passes the A quality rules | lower M1 |
| D16 | M1 also normalised with a rule-matched ceiling | none |
| D17 | Directional p-values: normal tail beyond the empirical null, max over two nulls | conservative |
| D18 | F_conv,func scaling and label union defined; per-class rows post hoc (exploratory) | sets F_conv,func |
| D19 | Same-site split by SOG identity at τ = 0.2 | changes r_H1s, not F |
| D20 | O2 tune/test split excludes crossover families (leakage fix) | smaller test set |
| D21 | Losses not FDR-supported: audit-level only, labelled in every report | weaker loss support, stated |
| D22 | Curveball swap count based on min(#1, #0) instead of #1; mixing checked at half and double | none expected |
| D23 | Co-fold concordance rule; ipTM dropped; ablation negatives in the panel | fewer confirmations |
| D24 | Gains reported with simulated and empirical FDR; paralog-branch and fragment checks added; false-loss controls must be bad models | gains less certain, stated |
| C-D7 | Loss cut re-derived under measured noise; about zero simulated recall | none |
| D25 | 0 convincing in v1: disjointness over-exclusion checked (bug fixed with a test); post hoc literature f2 reported separately; co-fold the supported pool | can only add cases |
| D26 | F_conv,func uses family-pair δ; specific-ligand sensitivity reported | headline nearly unchanged |
| D27 | G5 scope: rerun only stages that see new data; storage guard | O1 results get a reproducibility check only |
| D28 | Symmetric misassigned-orthology check (O2 for O1, PANTHER for O2); pool shifts to O2 | removes most O1 gains |
| D29 | Part A headline stays per family pair, always shown with per pair; symmetric D11 in descent | framing; lowers inflated O2 M1 |
| D30 | G4 freeze; PANTHER side-family check report-only; validation track built | none on discovery pool |
| D31 | Validation recall 0/7 from within-OG scope; exploratory cross-OG extension with its own null | can only add candidates |
| D32 | Cross-OG extension 0/7, gate failed, no genome-wide run | none; negative result |
| D33 | O2 artifact orthology check = DIAMOND RBH (post hoc, after O2/PANTHER were 'unknown'); 17 vs 34 O2 gains under D28 as written; RBH failures = other paralog subfamily of the same HOG | fewer cases |
| C-D31 | 0/7 causes are mixed, not only between-OG | wording; 0/7 stands |
| D34 | Specificity: BH unattainable; population-level real-vs-null excess and tail-fit sensitivity reported post hoc | can only add; outside headline |
