# Round 1 — Statistics review (phylogenetic comparative methods, large-scale similarity)

Scope: DESIGN.md v0 §3.3, §4, §5.1–5.2, §8.3. Verdict: the Part A pipeline is buildable, but F_anc/F_conv as written is
not a well-defined estimand (depends on k, sampling density and null construction), two tests are biased by
construction (directional Δ, specificity-change permutation), and the per-query significance rule cannot be met.
Changes below are ordered by importance; cost = extra build time.

## A. Similarity decomposition (§4.2)
1. **The significance rule is unattainable and the wrong scale.** With 1,000 null draws the smallest empirical p is
   1/1001 > 1e-3, so no pair passes. Worse, a per-draw p ignores database size: with N_db ≈ 3.2 M targets, p = 1e-3
   is expected ~3,200 times per query by chance. Change: fit a parametric tail to the per-query null (GPD over the
   top 5 % of ≥ 10,000 size-matched draws, or a Gumbel fit as in BLAST) and use **E = N_db · p_tail**; call pairs
   with E < 0.01 (BH q = 0.05 over the E-derived p is an equivalent alternative; state q explicitly). Cost ~1 h.
2. **Do not estimate F from kNN lists.** Top-50 lists are a sample biased toward homologs (75 % kinase→kinase),
   so H-class mass depends on k and on how many orthologs exist (a family with n members gives O(n²) H1 pairs;
   Swiss-Prot gives many human/model paralogs). Define the estimand on the **full pair population**, estimated by
   stratified sampling: per homology class c, draw pairs uniformly (≥ 10⁶ for N, all/most for H1–H3 which are
   enumerable from HOGs/Pfam), estimate the excess rate r_c = P(E < 0.01 | c) − r_null, and population counts n_c.
   Then **F_anc = Σ_{c∈H1–H4} n_c r_c / Σ_c n_c r_c**, F_analog = n_N r_N / Σ. kNN is kept for discovery only.
   Report the estimand at three preregistered weightings: per pair, per query site (each site weight 1), and per
   family pair (each root-HOG × root-HOG cell weight 1, removes O(n²) sampling inflation). Cost ~2 h.
3. **The null removes the thing it should detect.** Random size-matched non-homologous pockets *are* the N
   population, so N excess over that null is ~α by construction and F_conv is a tautology. Use a **covariate
   ("physics") null**: quantile regression of s on the covariates of *both* pockets (volume, burial, pLDDT, residue
   class composition, polarity, AF-mode flag, model source), fitted on random N pairs; excess = s above the
   conditional 95th/99.9th quantile. Convergence-like excess is then similarity *not explained by generic pocket
   physicochemistry* (the Gao–Skolnick point). Rename F_conv → **F_analog** (static similarity cannot establish
   convergence); reserve "convergence" for pairs passing the directional test (#7). Cost ~2 h.
4. **Homology misclassification biases F_analog upward.** N contains undetected homologs. Estimate contamination
   on a random 2,000-pair sample of significant N pairs with HHblits/HHsearch (or Foldseek 3Di+AA at E < 1) → rate ε,
   report F_analog × (1 − ε) with CI. Report H4 separately and as a sensitivity (assigned to anc vs analog), since
   TM ≥ 0.5 can itself be convergent at fold level. Exclude H0 (same-protein) pairs from both numerators; count
   "unknown" pairs and give bounds with them assigned to either side. Cost ~2 h CPU.
5. **Uncertainty is dyadic.** Pairs share two clusters; a one-way bootstrap over root orthogroups understates
   variance. Use a **pigeonhole (two-way) bootstrap** over root HOGs (resample HOGs with multiplicity, keep a pair
   with weight m_i·m_j); 200 replicates suffice for SE. "Inverse clade size" needs a definition: use per-HOG
   Henikoff/GSC sequence weights on the species tree, pair weight w_i·w_j. Cost negligible.
6. **Hubness and the PPV correction.** Hubs generate spurious N pairs: symmetrise scores (mean of query-side and
   target-side z), and report the result with targets of k-occurrence > p99 removed. The F × PPV correction needs
   assumptions stated (true analog pairs share ligand class with prob. ≈ 1; random pairs at the background rate;
   annotated subset representative). The annotated subset is human/enzyme-biased, so report PPV per stratum and as
   a bound, not a point correction. Cofactor-site sharing (ATP/NAD/FAD) should be its own stratum.

## B. Directional convergence (§4.2)
7. **Δ > 0 is guaranteed by selection (winner's curse).** Pairs are selected because s(a, b_derived) is high, so
   any sister comparison regresses to the mean; the role-swap null doesn't remove this. Change: **cross-fit**:
   choose the pair using one derived representative, compute Δ with *held-out* derived species vs sister species
   (pLDDT- and source-matched), and compare the Δ distribution with Δ for (i) homologous H1 pairs selected the
   same way and (ii) random below-threshold N pairs. Require ≥ 2 held-out derived species. Cost ~1 h.

## C. Phylogenetic signal and tree consistency (§4.1)
8. **Fritz–Purvis D under these conditions.** D needs branch lengths (TimeTree; fall back to Grafen, flagged) and
   a resolved tree (resolve polytomies randomly, 10 resolutions, report spread). It is unstable at extreme
   prevalence; restrict to 0.1 ≤ prevalence ≤ 0.9 and **report the excluded fraction** (it excludes the most
   conserved pockets, i.e. the strongest evidence for descent; report the prevalence distribution too). Observation
   error (recall ~0.55–0.7: only 55 % of AF-mode pockets have a non-AF counterpart, 60 % for Foldcomp v4) pushes D
   toward 1. Calibrate by simulation: threshold-BM trait on our tree + fitted emission errors + the real U pattern
   → expected D under "perfect descent observed noisily"; report D relative to that. U is not missing at random
   (clustered in low-pLDDT eukaryotes); include the sensitivity with U imputed both ways. Cost ~1 h.
9. **M4 null.** A permutation preserving only per-species rates does not preserve SOG prevalence, and Dollo counts
   depend strongly on prevalence. Use a **fixed-margins (curveball) null** preserving both per-species and per-SOG
   totals, and the same U mask. Also report Dollo counts with false negatives allowed (Dollo + error, below).

## D. Observation model and reconstruction (§3.3, §5.1)
10. **Emission probabilities must be site-specific and estimated, not fixed at 0.7.** The 0.7 is ligand-site
    recall in benchmarks, not P(observe P | true P) for an A-eligible footprint (pLDDT ≥ 70). Estimate
    P(obs A | true P) as a function of footprint pLDDT, burrow score and size from replicate observations of the
    same protein: tier-1 v6 vs tier-2/tier-3 copies, overlapping AFDB fragments, and BURROW seed/preset reruns on a
    random 1,000-site sample. False-positive rate from the same replicates. This is also required because recall
    and loss rate are **not jointly identifiable** from tips alone (scattered A = losses or misses).
11. **Mk fitting.** Per-SOG rates from ~20 tips are poorly determined; pool rates hierarchically (empirical Bayes
    over SOGs in strata of size/score) and use marginal posteriors under those rates. Treat "posterior ≥ 0.9" as a
    local fdr and control the **Bayesian FDR** of the called set (mean(1 − posterior) ≤ 0.1). Run Sankoff with
    gain:loss 1:1, 2:1, 3:1 and the ML-estimated ratio; keep events stable across all (report, don't tune).

## E. Event calls and multiple testing (§5.2, §5.3)
12. **Replace the label-permutation null for event counts with a parametric bootstrap**: simulate under the fitted
    no-event model (ancestral presence everywhere + fitted emission errors + real U/G masks), run the full calling
    pipeline, and count calls → empirical FDR per event type and per clade size. Use the same simulator to set the
    ≥ 3 species / 80 % / ≥ 2 sister thresholds to a target FDR of 10 % before seeing real calls (§8.3 already has
    the toy machinery). Permuted labels test "no tree structure", which is not the relevant null.
13. **Specificity change test is confounded with drift.** Under any Brownian descent, between-clade cosine <
    within-clade, so a tip permutation (non-independent tips) flags phylogenetic signal, not a shift. Change: test
    a **mean shift on a branch** against a BM (or OU) null fitted on the tree (e.g. phylogenetic ANOVA with
    simulated null, or l1ou-style shift on the first 2–3 PCs of the fingerprint), and require the shift to exceed
    the replicate-noise floor (same-site AF vs non-AF cosine median 0.886, p10 0.69). BH (q = 0.05) across all tested
    SOG × branch hypotheses; report the number tested.
14. **Report the discovery funnel with denominators.** For each event type: SOGs tested → candidates (FDR-
    controlled) → artifact-filtered → filter 1/2/3. Filters 2–3 are evidence, not tests: don't feed case counts back
    into prevalence estimates (annotated proteins are a biased sample). Fragment-pair negative control (§5.3 iii) is
    a direct FDR check: report its call rate with exact binomial CI.

## F. Minor
- M1/M2: model decay with a logistic GLMM (random effects for root HOG and species pair) rather than binned curves;
  the ortholog-vs-paralog test at matched identity needs identity as a smooth covariate, not strata.
- Size control by decile is coarse; use both pockets' volumes (2-D bins or the #3 regression).
- State the seeds for bootstrap/simulation replicates beyond [0,1,2] (e.g. seed s → replicate stream s).

None of these changes alters project scope; all reuse planned components (§8.3 simulator, HOGs, FAISS).
