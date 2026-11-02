# Round 1 — External PI review of DESIGN.md v0
Reviewer: external PI (structural/evolutionary biology). Read: PROMPT.md, DESIGN.md v0, LITERATURE.md, docs/tools_pocket.md,
docs/tools_ancestra.md, COVERAGE.md (snapshot 20261030T085136Z), docs/COFOLD.md (skim).

## Verdict
The framework is sound: site → positional homology via the MSA → P/A/U/G states → reconstruction, with absence treated
as an observation that has a known error rate. That is the right core. Two problems remain. (i) As written, the Part A
estimand (F_anc/F_conv) mostly reflects how the samples were drawn and the k=50 truncation, not biology. (ii) The plan
covers about 2× what 22 h allows, and Part B candidates are generated "data-first", which makes filter 2 (experimental
support) the last and least likely gate. Changes are listed by priority, each with a rationale and a cost.

## Changes (must-do = 1–7; should-do = 8–12)
1. **Redefine the Part A estimand so it does not depend on sampling.** "Share of excess mass in H1–H4" grows with the
   number of homologs each query happens to have. Swiss-Prot (tier 2) holds many orthologs of well-studied enzymes, so the
   share is inflated by redundancy. Report two quantities. (a) A per-pair rate: P(significant | class) for each class, with
   the denominator being all candidate pairs of that class. (b) The current mass share, after redundancy weighting
   (cluster sequences at 50 % identity in tiers 1–2, or weight by 1/|HOG members per species|). Report F_conv only as
   "similarity without detectable homology (upper bound on convergence)". Cost: ~1 h of analysis code, negligible CPU.
2. **Replace the kNN k=50 universe with a threshold (range) search, or report the result as a function of k.** Homologs fill
   the top-50 of most neighbour lists (75 % of kinase top-10 are kinases), so with k=50 the N pairs are truncated and F_conv
   is biased low. The bias size depends on family size. Use a FAISS range search at a size-conditioned threshold (s_null99.9
   per volume decile), capped per query, and report F at k ∈ {10, 50, 200} as a sensitivity check. Cost: 1–2 h CPU.
3. **Fix the empirical-p resolution.** 1,000 null pockets per query cannot give p < 1e-3. Fit a parametric tail (e.g. GPD or
   skew-normal on the volume-decile null, pooled per decile from about 10^5 samples) or use 10^5 shared nulls per decile.
   Without this, the significance and BH steps are degenerate. Cost: minutes.
4. **Measure a replicate ceiling for pocket detection and normalise M1 by it.** Ortholog conservation of 0.6 is meaningless
   without a same-protein reproducibility figure. Use accessions present both as v6 (tiers 1–2) and as Foldcomp v4 (tier 3),
   plus the BURROW seed runs s1/s2, to estimate P(site re-detected | same sequence). Report M1 as observed/ceiling. The same
   numbers parameterise the error layer in §5.1 (r, false-positive rate) empirically instead of using r = 0.7 from a
   different benchmark. Cost: ~1 h, cheap compute (overlap ids + 2k re-runs).
5. **Use the atlas's unfiltered pockets (`pockets_all`) in the state calls.** A footprint with a pocket that was dropped by
   the score < 0.03 or lining pLDDT < 70 filter must be P(weak) or U, never A. This is cheap, applies globally (not only to
   candidates), and removes the biggest structural source of false "losses" before the expensive rescue step. Cost: about
   30 min, if `pockets_all` is readable per unit (check POCKETOME.md).
6. **Add an annotation-first track for Part B, run in parallel with the data-first scan.** Start from footprints with
   experimental evidence: UniProt ECO:0000269 binding/active sites, BioLiP2 or PDB holo sites, M-CSA. Project them through
   SOGs and look for P→A or chemistry shifts. Include known-answer positives with expected outcomes: pseudoenzymes, i.e.
   pseudokinases or catalytically dead paralogs with a lost or altered ATP/active site (expected: loss/specificity change at
   a conserved position); steroid-receptor ligand-binding-domain paralogs in human/zebrafish (specificity change after
   duplication); and the M-CSA convergent pairs (Riziotis 2025). If the pipeline cannot recover these, its de novo cases are
   not credible. Cost: 2–3 h of agent time; it replaces some of the open-ended scanning.
7. **Cut scope to fit 22 h.** Cut or defer: (a) O2 (~150 context proteomes). Tier-3 near-complete taxa are mostly
   archaea/bacteria (archaea 9.3× enriched; large eukaryotes only as 4,000-entry slices), so O2 adds little to the eukaryotic
   clades where most tier-1 events will be. Do O2 only if O1 persistence filters fail, and then only for targeted clades.
   (b) M3 Fritz–Purvis D (redundant with M4, and the polytomy-heavy guide tree is weak). (c) The ML Mk2 + error layer:
   Sankoff with the empirical error rates (change 4) plus the permutation null is enough for v1; Count/GLOOME is a follow-up.
   (d) The QfO re-benchmark (Ancestra already has one; cite it). (e) Foldseek structure-aware MSAs unless benchmark 8.1 shows
   that sub-30 % identity orthogroups fail. Saves an estimated 5–7 h.
8. **Make the permutation null the explicit FDR for each candidate class.** Report the expected number of false
   gains/losses at the candidate thresholds (the label-permutation null preserves per-species and per-SOG rates) next to
   the observed counts. Thousands of SOGs are scanned, so tens of "gains" will arise from noise. Cost: already designed;
   it only needs to be reported.
9. **Fix the uncertainty for cross-family pairs.** An N pair spans two orthogroups, so a one-way cluster bootstrap over
   root orthogroups understates the variance. Use a two-way (pairs-cluster) bootstrap or resample families and keep the
   pairs whose two families are both drawn. Report the main stratified F also by superkingdom and by tier. Cost: small.
10. **Choose co-folding cases where both states are cheap to fold.** Non-human, non-yeast proteins need rate-limited MSAs
    (1–5 min each, server). Prefer human paralog pairs (neo/subfunctionalisation) or human/yeast orthologs for the deep
    dives, and give each case a fixed cap of ≥ 2 proteins × (ligand, decoy, apo/ablation) × 3 seeds. Once the AF3 per-job
    runtime from the benchmark is known, fix the case list by ~T+12 h. Cost: none; this is planning.
11. **Keep threshold selection separate from case discovery.** τ, the Jaccard merge cut-off and the significance thresholds
    are chosen on benchmark 8.1. Benchmark families must be flagged, and any case drawn from them must be marked
    "benchmark-family" in the case table. Leiden on the column-Jaccard graph has a resolution parameter: preregister it and
    report SOG size and split statistics (chaining from union-find at Jaccard 0.3 can create mega-sites; cap or report).
12. **Add hard time gates to §10.** T+6 h: sites/SOGs and benchmark 8.1 numbers done, τ frozen. T+10 h: Part A v1 numbers with
    CIs on tier 1 only (tier 2 added if time permits). T+12 h: candidate list frozen for co-folding. T+18 h: web app feature
    freeze. Missing a gate triggers the cut list, not overtime.

## Things that would make the results misleading if left unaddressed
- Calling the N class "convergence". Without the directional Δ test it is only "no detectable homology" (remote homology
  leaks in, and small greasy pockets match by physics, as Gao & Skolnick showed). Keep the language in REPORT.md strict.
- Mixing v4 and v6 models inside a SOG. Model version is correlated with tier and therefore with taxon, which can produce
  clade-consistent artificial "changes". Re-describing case members from one source is good; also require the M1–M4
  results to be stratified by model version.
- Size leakage in the fingerprint (TOUGH-M1 L-mode at 0.95). Every similarity claim needs the size-only baseline beside it,
  and the specificity-change test must condition on volume (already designed; keep it mandatory).
- Eukaryotic tier-1 pLDDT (38 % of lining residues < 70). Absences in plants, kinetoplastids and Plasmodium will be mostly U.
  Report the per-clade informative fraction so that "few losses in X" is not read as biology.
- The fingerprint was trained on BioLiP ligand-sharing pairs, which may include homologs. The "shared ligand" PPV
  validation of N pairs is therefore not fully independent of the training data. State this, and if possible exclude
  training-set proteins from the PPV estimate.

No change is needed to the goals, the P/A/U/G state logic, the artifact negative controls in §5.3, or the co-folding
control design (ablation plus decoy plus holo positive control). These are the design's strongest parts.
