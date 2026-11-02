# DESIGN — How do binding pockets evolve across the protein universe?

Version: **v2 — FROZEN at G1** (2026-10-30 10:47 CET; config/prereg.yaml sha256 52f94618…0d83). v1 (10:39) was revised after review round 1 (six reviewers, reviews/round1/). The change log
mapping each review point to a decision is §12. v0 is in git history (commit "DESIGN v0 draft").
Inputs: PROMPT.md, docs/tools_pocket.md, docs/tools_ancestra.md, docs/COFOLD.md, LITERATURE.md,
results/coverage/20261030T085136Z/COVERAGE.md, POCKETOME.md.

## 0. Principles and safeguards
- **One command per stage**, parameterized by a snapshot manifest (`data/snapshot/<id>/manifest.parquet`) and
  `config/prereg.yaml`. Full-atlas re-run = new snapshot + same commands. Streaming per orthogroup / per file / FAISS;
  nothing loads all fingerprints.
- **Preregistration enforced**: every threshold, metric, case-ranking score and co-folding readout lives in
  `config/prereg.yaml`; its sha256 is recorded in RESEARCH_TREE.md when frozen (gate G1). Every output carries
  `{snapshot_id, prereg_sha256, git_commit, seed}` in a sidecar `.json`. `scripts/check_prereg.py` fails if an output's
  prereg hash differs from the frozen one without a dated deviation entry in §12.
- **Tune/test split**: thresholds that are tuned (τ, alignment method, absence pLDDT cut) are tuned only on the
  *tune* half of benchmark orthogroups (hash split by root HOG id: `sha1(hog_id) % 2 == 0`); results are reported on
  the *test* half and with a sensitivity grid.
- **Seeds**: `seeds: [0, 1, 2]` (stochastic steps); bootstrap seeds derived as `seed*1000 + replicate`.
- **No number without a producer**: every number in REPORT.md comes from `results/numbers.json` (written by
  scripts) or a sidecar; `scripts/check_report_numbers.py` checks this. Each stage logs counts dropped per reason
  (funnel tables). No bare `except:`; failures are counted, not swallowed.
- **Independent verification**: an agent that did not produce a result reviews it against a checklist
  (`docs/VERIFY_CHECKLIST.md`) from committed artifacts only, re-runs a slice and compares hashes before the PI relies
  on it. Toy data with known answers (including "no event" and "no convergence" worlds) run before real data.
- **Development subset**: results on snapshot `20261030T085136Z` are preliminary hypotheses.

## 1. Data
| set | content | structures | pockets | role |
|---|---|---|---|---|
| tier 1 | 46 AFDB v6 reference proteomes | 555,846 | 1.79 M | orthology backbone, presence/absence |
| tier 2 | Swiss-Prot v6 (≈11k taxa, median 2/taxon) | 436,886 | 1.44 M | presence-only evidence, Part A pair universe, annotations |
| tier 3 (atlas) | AFDB v4 via Foldcomp, 9.8 % of units | 22.96 M | 68.1 M | context taxa ≥ 90 % present (mostly prokaryotes/protists) |
| **tier L (local)** | ~30–40 eukaryotic reference proteomes near tier-1 species, pockets computed **here** from the Foldcomp v4 DB with the atlas worker (same method record as tier 3) | ~0.4–0.6 M | ~1.5 M | densify vertebrates/insects/plants/fungi so persistence across species is testable |

Biases (COVERAGE.md): tier-3 done set over-represents Archaea and small phyla; no near-complete vertebrates, insects
or land plants in tier 3 (hence tier L); tier-1 eukaryote pLDDT low (mean 73.7). Model source (v6 file vs Foldcomp
v4) is a covariate everywhere; case studies re-describe all members from one source (Foldcomp v4, local, which has
every AFDB v4 model including human; same method record as tier 3).

## 2. Orthology
- **O1** = Ancestra HOGs on the 46 tier-1 proteomes (done, 36,736 root OGs). **O2** = O1 + tier-L proteomes + ≤100
  prokaryote/protist taxa ≥ 90 % present in the snapshot (deterministic selection script). ~1.5 M proteins, ≈1.3 h.
- ID gap (7 % of tier-1 atlas accessions not in gene-centric proteomes) mapped to canonical gene members; unmapped
  counted and excluded.
- Tier 2 and other tier-3 proteins: DIAMOND placement into O2 HOGs → presence-only evidence.
- Species tree: NCBI guide tree (polytomies). Divergence times: TimeTree node ages for tier-1/L species (table
  fetched once; DATA.md), rank proxy otherwise (flagged).
- Independent check: PANTHER family agreement for every case and for a 2,000-pair sample (reported).
- Deep (Eukaryota-level) splits are flagged "orthology-uncertain" (Ancestra recall is low there).

## 3. Matching pockets across the tree
### 3.1 Sites
Per protein, BURROW merged pockets (fragments converted to UniProt numbering) with residue Jaccard ≥ 0.3 → sites by
union-find; we record the largest-site size distribution to detect percolation (sites > 60 residues flagged) and run
the analysis also on raw pockets as a sensitivity check.

### 3.2 Site orthogroups (SOGs)
Per root OG with ≥ 3 species: FAMSA MSA of model sequences (identity-checked against the atlas model); columns;
site → column set; inverted index column → sites; two sites in different members are positionally homologous if
column Jaccard ≥ τ (τ tuned on the tune half of the benchmark in §8.1 from {0.2, 0.3, 0.4, 0.5}); SOGs = Leiden
communities on that graph. Footprint = columns in ≥ 50 % of member sites; insertion columns are allowed in the
footprint (lineage-specific loops, e.g. GID1 lid). Low-identity OGs (mean pairwise identity < 25 %): FoldMason
structural MSA as alternative, compared on the benchmark; used only if it wins on the tune half.

### 3.3 Observation states per SOG × member
- **P** present: a site overlapping the footprint (Jaccard ≥ τ).
- **W** weak-present: no kept site, but an atlas-filtered pocket (`n_pockets_all` > kept; re-described with
  `--min-score none`) or a rescue-setting site overlaps the footprint. W is never called absence.
- **A** absent-observed (primary): footprint ≥ 80 % aligned non-gap; footprint mean pLDDT ≥ 90 and frac(< 70) ≤ 0.2
  (Akdel 2022); protein length 0.7–1.3 × OG median, footprint not in a terminal truncation; no signal-peptide/TM
  overlap of the footprint; status ok. **Sensitivity stratum**: pLDDT ≥ 70.
- **U** uninformative: everything else, including footprints flagged as likely oligomeric interface (UniProt
  subunit/PDB assembly evidence for case members) and metal-only sites.
- **G** gene absent in the species (no HOG member in a complete proteome).
- **Rescue** (symmetric, applied to every A call used by a gain *or* loss candidate): BURROW with `--min-score none`,
  no AF mode, seeds s1/s2, fast/accurate presets; a footprint site in ≥ 2/6 settings turns A → W.
- **Measured error rates** (replace the borrowed 70 % recall): (i) re-detection rate of a site in the same protein
  across model sources (v6 atlas vs Foldcomp v4 re-description of a 2,000-protein tier-1 sample) and BURROW seeds;
  (ii) false-A rate on experimentally evidenced sites (BioLiP2 / LIGYSIS sites in tier-1/2 proteins) by pLDDT band and
  organism. These rates feed the simulations of §5.

## 4. Part A — global patterns (preregistered)
### 4.1 Do pockets descend by orthology?
Framed as turnover: conservation is expected; we measure its rate and exceptions.
- **M1** P(site has a positional match in the other protein | footprint informative there) for ortholog vs paralog
  pairs, vs divergence time and sequence identity, normalised by the same-protein re-detection ceiling (§3.3).
  Ortholog-vs-paralog difference at matched identity = "pocket ortholog conjecture" test.
- **M2** fingerprint similarity of positionally matched sites vs divergence time, relative to the technical
  noise floor (same site, different source/seed) and the covariate-matched background (§4.2).
- **M3** tree consistency: fraction of SOGs explained by ≤ 1 gain/loss (Dollo parsimony), against a null that keeps
  per-species and per-SOG totals (curveball swap) and against simulated noisy-but-vertical evolution.
- Deferred (only if time): Fritz–Purvis D with simulation calibration.

### 4.2 Similarity decomposition: ancestry vs analogy vs convergence
- **Pair universe**: tier-1/2 sites (rep fingerprints). **Discovery**: GPU FAISS exact kNN k = 1000; for the estimand
  we use *all pairs above a significance threshold*, not a fixed top-k (reported also at k = 10/50/200).
- **Significance**: per-query, covariate-matched null: similarity to ≥ 10,000 pockets from non-homologous proteins,
  stratified jointly on volume, burrow score, lining polarity/aromatic fraction and mean pLDDT (both pockets'
  covariates), with a generalised-Pareto tail fit; E = N_db × p; significant if E < 0.01. A size-only and a covariate-
  only baseline (similarity predicted from covariates alone) are reported.
- **Homology classes** (first match wins): H0 same protein; H1 ortholog; H2 paralog (same root HOG); H3 sequence-
  homologous other (shared Pfam clan at the pocket region / PANTHER family / TED CATH superfamily / MMseqs2 E < 1e-3);
  H4 remote-structural only (Foldseek E < 1e-3 or TM ≥ 0.5 on the region domain; reported separately as
  "ancestry-uncertain"); N no detectable homology (Foldseek and MMseqs2 searched and negative). Pfam/PANTHER absence is
  not required for N (coverage 58–62 % in non-model proteomes); undetected-homology rate estimated by HHsearch-
  level check on a 200-pair N sample (or Foldseek TM-align mode) and reported as a correction.
- **Estimands** (all pairs, sampled by class):
  - per-class rate r_c = P(pair significant | class c), with n_c the number of pairs in the class;
  - **F_anc** = Σ_{c∈H1–H3} n_c r_c / Σ_c n_c r_c; **F_H4** likewise; **F_analog** = N share. Three preset weightings:
    per pair, per site (each query weight 1), per family pair (each pair of root OGs weight 1 — removes redundancy).
    F_analog is "similarity without detectable homology" — an upper bound on convergence, not convergence.
  - **F_conv,dir**: share of significant pairs (any non-orthologous class N, H3, H4) with a significant *directional*
    convergence signal (below).
  - **Functional validation**: enrichment of shared experimental ligand class (BioLiP2 chemotype, ECFP4 Tanimoto ≥
    0.5, or shared M-CSA/EC) among significant N pairs vs covariate-matched random N pairs; reported uncorrected and as
    PPV separately (PPV comes from the annotated subset; stated as such).
- **Directional convergence test (cross-fitted)**: for a significant pair (site a in family F; site b in family G with
  a SOG reconstructed in G), split G's derived-lineage members into selection and test halves; select pairs on the
  selection half; on the test half compute Δ = s(b_test, a) − s(b_sister, a) where b_sister are orthologous sites in
  G's sister lineages / reconstructed ancestral-proxy (outgroup) members. Calibrated against Δ computed identically for
  ortholog pairs and on a toy no-convergence world. Convergence ⇔ Δ > 0 beyond the calibration null (BH across pairs).
  Prediction stated in advance: F_analog > 0 even with no convergence (Gao & Skolnick pocket-space degeneracy).
- **Uncertainty**: two-way (crossed) bootstrap over root OGs (1,000 reps); stratification by model source (v4/v6),
  pocket size, pLDDT ≥ 90; inverse-clade-size reweighting.

## 5. Part B — events
### 5.1 Reconstruction
Per SOG on the HOG tree (nested speciation/duplication nodes): Sankoff parsimony, gain cost 2, loss cost 1, U/W
handled as missing/present-ambiguous. Event-calling thresholds are set to a 10 % FDR using a **simulated no-event
null** (same tree, observed per-species error rates from §3.3, no true events) and checked on simulated planted
events (recall reported). Mk/ML layer deferred.

### 5.2 Candidate definitions
- **Gain**: SOG P in ≥ 3 species of clade C (≥ 80 % of informative members), A (rescued) in ≥ 2 informative members
  from ≥ 2 sister/outgroup lineages; parsimony places one gain on C's stem.
- **Loss**: SOG P in ≥ 2 outgroup lineages; A (rescued) in ≥ 2 species of derived clade D (≥ 80 % of informative
  members); one loss on D's stem.
- **Specificity change** (residue-first): SOG P throughout; ≥ 2 clade-consistent substitutions at footprint
  residues (BLOSUM62 < 0) between clade D and its sister; supported by a fingerprint mean-shift on the D branch that
  exceeds the technical noise floor (same-site replicate cosines across seeds/sources) after regressing out volume,
  pLDDT and model source; clade-level permutation; BH across all tested SOG × branch pairs. Paralog branches after
  duplication included. Sensitivity floor measured with a gatekeeper-type single-substitution control.
- **Convergence**: pairs passing the directional test, with shared experimental ligand class.

### 5.3 Two tracks
- **Discovery track**: candidates from the SOG reconstruction, ranked by the preregistered score
  `S = persistence (n species) × min(state confidence) × evidence count`.
- **Annotation-first track** (positives with known answers): experimentally annotated sites (UniProt ECO:0000269/
  ECO:0007744 binding/active sites, BioLiP2, M-CSA) traced through SOGs; plus the round-1 control list
  (reviews/round1/prior_work.md): steroid-receptor LBPs (specificity), NR4A/Nurr1 pocket loss, PfLDH vs MDH, E. coli
  IDH NADP→NAD, GID1 vs HSL, D14/KAI2, lysozyme vs α-lactalbumin; detection-limited convergence (serine-protease triads,
  CA classes, SODs) reported as such; negative controls uricase/GULO (should be G), M. leprae vs M. tuberculosis
  (truncations must not become losses). Pipeline recovery of these is reported as a validation table.

### 5.4 Artifact cascade (required before any gain/loss call)
footprint pLDDT (primary ≥ 90) and disorder; PAE (AFDB API if reachable, else TED domain boundaries as proxy, stated);
fragment/truncation and gene-model outliers (length, start Met, terminal truncation); signal peptide/TM; oligomeric
interface; misassigned orthology (PANTHER disagreement, Ancestra deep-level flag); pocket-finder instability (rescue);
model source (re-describe from one source); apo/holo state (AHoJ-DB lookup where the family has PDB entries).
**Negative controls**: (i) the whole Part B pipeline run on null inputs (curveball-permuted P/A, random pseudo-clades,
same-protein fragments) — funnel table for real vs null; (ii) simulated losses (masked pLDDT, truncation) must give U;
(iii) uricase/GULO/M. leprae controls.

### 5.5 Evidence filters and the bar
1. Persistence: change in ≥ 3 species of the derived lineage, single tree-consistent event.
2. Experimental: PDB holo ligand at the site (BioLiP2, residue overlap ≥ 50 %), UniProt experimental ECO binding/
   active site, M-CSA, ChEMBL pChEMBL ≥ 6 on the protein, cited mutagenesis/biochemistry. Evidence must be
   *about the change* (e.g. ligand bound in the state that has the pocket) or is recorded as "protein-level only".
3. Computational: footprint conservation/co-variation, rescue-setting consistency, single-source re-description,
   co-folding (§6), ancestral reconstruction fold (§6).
Labels: **convincing** = all three filters pass; **supported** = filter 1 + one of 2/3. Both reported; only
"convincing" cases are headline cases. Every identifier (PDB, DOI, ChEMBL) is checked against its API and a short
supporting quote/record is cached in `cases/<id>/evidence/`.
Case selection for deep dives: top-k by S per category + 1 randomly sampled passing case (seeded) as audit;
literature-famous examples count as "rediscovered positives", not discoveries.

## 6. Co-folding (AF3; docs/COFOLD.md)
- Comparable inputs: one MSA pipeline for all members compared within a case (ColabFold MMseqs2 server MSAs for all,
  including human/yeast, unless all members have cached AF3 MSAs); `--no-templates` primary (templates as secondary);
  Neff recorded; no MSA reuse across derived/ancestral boundary.
- Per case: ≥ 2 derived-state and ≥ 2 ancestral/sister-state proteins × relevant ligand × property-matched decoy ×
  3 seeds; positive control = known holo PDB complex of the family (ligand RMSD after pocket superposition);
  **site-ablation control** (footprint → Gly) benchmarked on the four positive controls first; **ligand-free refold**
  + BURROW re-description (separates "AFDB model missed it" from "sequence doesn't encode it").
- Preregistered readouts: ligand pLDDT, min protein–ligand PAE, ligand centroid distance to footprint, PoseBusters
  validity; ipTM recorded but not used alone (decoys reach 0.72–0.86). Decision rule: derived vs ancestral difference
  in ligand pLDDT and min-PAE must exceed the decoy-vs-ligand difference within the same protein.
- **Ancestral sequence reconstruction** for 2–4 deep dives: IQ-TREE `-asr` on the OG MSA/tree; fold the marginal ML
  ancestor + ≥ 2 posterior samples (AltAll-style); require a consistent pocket call; consensus-bias caveat.
- Budget: ≈ 80 s per 3-seed job at ~300 aa; a 300-job programme ≈ 7 GPU-h ≈ 2 h on 4 GPUs.

## 7. Web app
`uv run python -m webapp` (README): FastAPI + static HTML, 3Dmol.js vendored. Reads index files (parquet) from
results/; structures streamed from AFDB tars via the atlas `tar_index/` and from Foldcomp (no AFDB copy). Views:
protein → pockets in 3D; orthogroup → SOG × species matrix on the tree (P/W/A/U/G); case → evidence, filter status,
co-folding models, figures, literature links.

## 8. Validation & benchmarks
1. **Site matching across orthologs** (tunes τ / alignment): pairs of experimentally evidenced sites (UniProt ECO
   binding features, BioLiP2) for the same ligand in orthologs of different species → SOG recall/precision; tune/test
   split by HOG hash. Atlas site recall on LIGYSIS (shipped with burrow).
2. **Convergence detection**: M-CSA convergent pairs *restricted to sites BURROW detects* (detection rate reported);
   cofactor-site convergence positives (adenine/NAD/FAD; Barelier pairs from ProSPECCTs if obtainable); TOUGH-M1
   pairs on AF models as fingerprint sanity check.
3. **Toy worlds**: simulated trees with planted gains/losses + measured observation noise → event precision/recall;
   synthetic pair sets with known class composition → estimand recovery; a no-convergence world → Δ null.
4. **Orthology**: PANTHER agreement sample; QfO overlap reported only.
Every stage reports runtime.

## 9. Compute and storage
CPU: O2 ≈ 1.3 h; tier-L pockets ≈ 70 min (64 cores); FAMSA MSAs ≈ 1 h; homology labels by joins (HOG, TED, Pfam,
PANTHER) + MMseqs2/Foldseek on residual significant pairs (< 1 h). GPU: faiss-gpu kNN (k = 1000, 3.2 M pockets,
≈ 20 min on one 4090); co-folding ≈ 2–7 GPU-h. Storage: repo holds code + compact results (home quota ~50 GB, 39 GB
used); bulk in /tmp/pe_big (node-local, regenerable; recipes in docs/).

## 10. Timeline and gates (T0 = 09:45)
- **G1 (T+4, 13:45)** DESIGN v2 frozen; prereg.yaml hashed.
- **G2 (T+7)** O2 + tier L + sites + SOGs on a 5-family end-to-end slice; τ frozen on tune half.
- **G3 (T+10)** Part A v1 numbers on tiers 1–2; artifact cascade + null pipeline running.
- **G4 (T+12)** co-folding case list frozen (favour human/yeast members with cached MSAs where comparable).
- **G5 (T+18)** web app frozen; **G6 (T+22)** REPORT.md + verifier sign-off.

## 11. Risks
Deep-tree orthology recall; detection-limited catalytic sites (convergence positives); fingerprint early precision
(8 % top-100) → convergence via enrichment + directional test; tier-3 incompleteness → absences only in complete
proteomes (tier 1, tier L, ≥ 90 % taxa); interface and metal sites invisible in monomers.

## 12. Review log
### Round 1 (reviews/round1/*.md) — decisions
Accepted:
- External PI: estimand redefined (rates per class, redundancy weighting, F_analog naming); threshold-based pair set
  not top-k; measured re-detection ceiling; W state; annotation-first track with known positives; time gates.
- Statistics: tail-fitted null with E-values; all-pairs estimand with three weightings; covariate null; F_analog vs
  F_conv,dir; cross-fitted Δ; residue-first + noise-floor specificity test with BH; simulated no-event null sets
  FDR thresholds; curveball null for M3; two-way bootstrap.
- Data/algorithms/compute: GPU kNN k = 1000; homology by joins + targeted MMseqs2/Foldseek; Foldcomp local single
  source for cases; tier L local pockets; FoldMason fallback; inverted-index SOG building; 5-family slice gate.
- Prior work: field definitions (analogy/convergence/parallelism); H4 separate; joint covariate null; detection check
  for control sites; interface → U; residue-first specificity; insertion columns in footprints; control list.
- Structural/co-folding: primary A at pLDDT ≥ 90; measured false-A rate on known sites; technical noise floor;
  comparable MSAs, no-templates primary, preregistered co-folding readouts; ligand-free refold; ablation benchmark;
  ASR for deep dives; AHoJ-DB; signal peptide/TM filter; symmetric rescue.
- AI failure modes: prereg.yaml with hash check; null-input pipeline + funnel; rule-based case selection with audit
  sample; numbers.json + report check; independent verifier; tune/test split; "convincing" = all three filters.
Deferred (time): ML/Mk observation-error model (replaced by simulation-calibrated parsimony); Fritz–Purvis D;
QfO re-benchmark. Rejected: dropping context species entirely (PI review) — instead re-targeted to eukaryotic tier L,
because without it the persistence filter cannot be met near human/mouse/fly/plants.

### Round 2 (reviews/round2/pi_statistics.md, ai_structural.md) — decisions (v2)
All thresholds now live in `config/prereg.yaml` (42 reviewer-proposed values accepted, 2 deferred); where this
document and prereg.yaml disagree, prereg.yaml wins. Changes relative to v1:
1. **Pair unit and ancestry**: unit = site pair on different proteins. Ortholog/paralog pairs are split into
   *same-SOG* (positional homologs → ancestry) and *different-SOG* (homologous proteins, non-homologous sites → counted
   with analogy, reported separately as "H1/H2-offsite"). Significance requires E < 0.01 in both query directions.
2. **Null**: covariate regression of similarity on both pockets' covariates (volume, burrow score, polarity,
   aromatic fraction, pLDDT, model source) fitted on 10^6 random non-homologous pairs (seeded), then a per-query GPD
   tail on residuals; validated on a held-out 10^6-pair set (calibration plot is a required output).
3. **Expected false positives subtracted**: r_c = (significant_c − expected FP_c) / n_c; expected FP from the
   calibrated null; FP count reported next to every rate. Class totals: exact for H1/H2/H3-by-join; H4/N totals from a
   seeded 10^5 protein-pair sample. CI: two-way bootstrap over root OGs, 200 replicates (1,000 if time).
4. **Headline weighting = per family pair** (Swiss-Prot redundancy); per-pair, per-site also reported. Hubness:
   pairs counted once (symmetrised); CSLS variant reported as sensitivity.
5. **Specificity null**: empirical null from other sister-clade pairs of matched depth (permutation secondary).
6. **Tune/test split key** = smallest member accession of the root HOG (stable O1→O2); control-list families forced
   into the test half; annotations used for tuning τ are excluded from evidence filter 2 for the same SOG.
7. **Rescue k** tuned to a false-A target of 5 % on known-present sites (tune half), grid {1..4}; rescue also applied
   to P calls of candidates (symmetric). Measured ceilings stratified by pLDDT band and model source.
8. **Interface rule** applied to all A calls entering event candidates and to case members (not to M1 bulk rates —
   no proteome-wide assembly data; stated caveat). TED-domain rule: footprint spanning ≥ 2 TED domains needs PAE or → U.
9. **Co-folding**: primary readout ligand pLDDT; secondary min P–L PAE, centroid distance, PoseBusters; decision
   rule must hold in all 3 seeds; blinded scoring (hashed labels); positive-control gate (RMSD ≤ 2 Å) or case
   co-fold = "uninterpretable"; cofactor present in both states; homodimer fold for oligomeric case SOGs.
10. **ASR**: footprint mean posterior ≥ 0.8 required; shuffled-ancestor and extant-under-same-MSA controls.
11. **Yield preflight** (pre-declared): if < 5 gain+loss candidates pass at pLDDT ≥ 90 on tier 1+L, report the
    count and run the pLDDT ≥ 70 stratum as the declared secondary analysis.
12. **Cut order if behind** (from round 2): ASR → FoldMason → AHoJ/PAE lookups → extra prokaryote taxa in O2 →
    M2 → extra weightings → HHsearch correction → ligand-free refold. Never cut: rescue, null-input pipeline,
    simulated FDR, decoy controls, null calibration, verifier sign-off.
13. Producers never mark their own RESEARCH_TREE node done; the verifier (docs/VERIFY_CHECKLIST.md) does.
No third review round: remaining open items are implementation details tracked in RESEARCH_TREE.md.

### Deviation ledger (post-freeze; each entry dated, with direction of bias)
| id | date | deviation from prereg | reason | expected bias direction |
|---|---|---|---|---|
| D7 | 2026-10-30 11:47 | loss cut raised to ≥3 species & ≥3 outgroup lineages; loss candidates counted only after persistence filter | simulation calibration at 10 % FDR (prereg procedure) — not a deviation in rule, recorded for transparency | fewer losses (conservative) |
| D8 | 2026-10-30 12:13 | SOGs with footprint ≤ 3 columns flagged `small_footprint`, excluded from event candidates | Leiden chaining of small heterogeneous sites (4,507 SOGs) | fewer events (conservative) |
| D9 | 2026-10-30 12:14 | rescue k = 1 although no k meets the 5 % false-A target (prereg had no fallback) | k=1 minimises false-A (0.13 residual); residual propagated into simulations | fewer absences → fewer gains/losses (conservative) |
| D10 | 2026-10-30 12:15 | τ = 0.2 is at the grid edge (tune F1 still rising); τ = 0.3 sensitivity reported for every headline metric | literal prereg rule; grid not extended post hoc | none on choice; uncertainty reported |
| D11 | 2026-10-30 12:36 | A (absent-observed) additionally requires the member's footprint to span > 10 residues; smaller footprints → U | measured on tune half: false-A after rescue k=1 is 0.054 for > 10-residue sites vs 0.148 for ≤ 10 (results/observation/.../rescue_breakdown.csv); test half reported | fewer absences → fewer gains/losses (conservative) |
| D12 | 2026-10-30 12:37 | Part A significance uses a pooled GPD tail shape (per-query shape reported alongside) | toy test showed per-query shape estimates too light-tailed; decided before real results | per-query: higher F_analog (0.92 vs 0.82 per family pair) |
| D13 | 2026-10-30 12:37 | expected-FP subtraction reported two ways: nominal × λ (headline, upper bound on F_analog) and observed random-pair rate (lower bound; ≈0 because the random-pair population *is* the N population) | deep-tail null is 186× miscalibrated vs nominal at E<0.01; the circularity was anticipated in round-1 statistics review | bounds bracket the truth |
| C1 | 2026-10-30 12:40 | timestamp correction: decision times D2–D13 and the G1 time were initially logged from an estimated clock (13:30–17:00); corrected to git commit times (G1 = 10:47, G2 = 12:13). prereg.yaml `pi_freeze.date` still reads 13:30 and is left unchanged to keep its hash; the true freeze time is the commit time of d586f1a (10:47). | process error, no effect on results | none |
| D15 | 2026-10-30 12:55 | M1 primary uses "qm90": a P member counts as informative only if its footprint also passes the A-primary quality criteria (symmetric with A); prereg-literal (any P informative) reported alongside (+8 pp) | compare like with like: literal rule lets low-quality P but not low-quality A into the denominator | lower M1 (conservative for "descent") |
| D16 | 2026-10-30 12:55 | M1 ceiling normalisation additionally reported with a rule-matched ceiling (partial overlaps excluded, as in SOG states); prereg ceiling (error_rates.parquet) reported too | prereg ceiling not on the M1 scale (ratio > 1 for close orthologs) | none (both shown) |
| D17 | 2026-10-30 14:16 | (Part A PA15) directional-convergence p-values use a normal-fit tail beyond the empirical null range; p = max over the two nulls (ortholog-Δ and random) | empirical p cannot reach BH levels with the null sample size; max-p is the conservative combination | sets F_conv,dir (1.1e-6 per pair, 3 passes) |
| D18 | 2026-10-30 14:16 | (Part A PA16) F_conv,func = per-family-pair excess of shared experimental ligand class over covariate-matched random N pairs, scaled by the significant-N fraction; label sources = BioLiP2 + UniProt experimental sites + M-CSA | prereg named the estimand but not its scaling or label union; definition and label union fixed before the result; CORRECTION (verifier 7327c7a): the per-class scaling was committed 4 min AFTER the first functional output, so per-class rows are post hoc and reported as exploratory | sets F_conv,func (0.078 [0.004,0.169]) |
| D19 | 2026-10-30 14:16 | (Part A PA14) same-site split uses SOG identity at τ = 0.2; sites without a SOG fall back to the footprint-overlap placeholder | real SOGs now available | changes r_H1s (0.010), not F |
| D20 | 2026-10-30 14:20 | O2 tune/test split: an O2 family is 'test' only if none of its members were in an O1 tune family; families whose members span halves are 'crossover' and excluded from tune and test reporting | recomputing the split from O2 membership moved O1-tune families (e.g. CDK2) into O2 test → leakage of tuned choices into test | smaller O2 test set; no effect on point estimates on all data |
| D21 | 2026-10-30 14:23 | Loss calls: no threshold reaches 10% simulated FDR with non-trivial recall under the measured noise, so losses are not FDR-supported. They go forward only as audit-level candidates that must pass the full artifact cascade and all three evidence filters at case level, and they are labelled 'not FDR-calibrated' in every report. Gains keep the prereg minima (sim FDR 0.05, recall ~0.1%). G3 slips to ~15:45 for the null funnel | the prompt requires lost-pocket cases and a false-loss negative control; dropping losses would remove that, and relaxing the FDR target would be tuning to results | loss cases carry weaker statistical support; stated, not hidden |
| D22 | 2026-10-30 14:58 | Curveball null uses masked checkerboard swaps at 100 × min(#1, #0) instead of the prereg 100 × #1 | the prereg count stalls on the sparse, masked presence matrix (the swap chain cannot finish in budget); the swap count is a mixing parameter, not a statistic | none expected on the null distribution once mixed; mixing should be checked by comparing 50× and 200× |
| D23 | 2026-10-30 15:06 | Co-folding case call needs concordance: a state difference counts only if ligand pLDDT AND (centroid distance to the predicted site OR site-residue recall) move in the same direction in ≥2 of 3 seeds, with thresholds frozen on the calibration panel. ipTM is not used as a readout; PoseBusters is reported only. The calibration panel includes site-ablation negatives as well as decoys, because cases compare protein variants, not ligands | controls a6c3030: ablation ≈ decoy cost in ligand pLDDT; PYGL ablation was caught only by centroid/recall; ipTM moves 0.08–0.22 under total site destruction | fewer co-folding confirmations (conservative) |
| D24 | 2026-10-30 15:14 | Part B gains are reported with BOTH error estimates: simulated FDR 0.05 (prereg) and the empirical curveball-null FDR (~16–20% before rescue, ≥20 seeds). No post hoc threshold change. Paralog-branch gains (gain just below a duplication; 33/51) and fragment-length flags enter the artifact cascade as 'orthology/duplication' and 'fragment' checks. The false-loss negative controls must be losses whose A side is a demonstrably bad model (low footprint pLDDT, fragment, high PAE) in the pre-D11/pre-rescue funnel, not rescue-unstable high-quality models | verifier e70b420: empirical null FDR exceeds the 10% target; rescue-overturned controls show instability, not bad predictions; the prompt requires a false loss that is really a bad prediction | gains labelled less certain than the simulated FDR suggests; stated, not hidden |
| D25 | 2026-10-30 15:18 | Case list v1 has 0 convincing cases (f2 = 0 for gains). (a) Correctness check, not a relaxation: annotation_disjointness excludes only sites actually in the τ-tuning set (tune-half families, §8.1); confirm the builder doesn't also exclude test-half or non-tuning annotations. (b) Published biochemistry or mutagenesis found by a literature search (PubMed IDs, about the change, reviewed by a second agent) is added as an f2 line: the prompt lists it but the prereg f2 list omits it. Convincing counts are reported with and without (b), and (b) is labelled post hoc. (c) With no convincing cases, co-folding (the structural-confirmation step) runs on the top-5 'supported' cases per category by S plus 1 seeded random audit case (prereg case_selection applied to the 'supported' pool); headline counts stay convincing_only | prereg rule correct in intent (no circularity) but possibly over-applied; prompt explicitly names published biochemistry/mutagenesis; co-folding needs cases to confirm | (b) can only add cases → reported separately; (a) may add cases if over-exclusion is found |
| C-D7 | 2026-10-30 15:24 | Correction to D7: the loss cut in D7 (>=3 species, >=3 outgroup lineages) came from the placeholder-noise simulation and is superseded. Under measured noise (D9 rescue k=1, D11, false-W 0.35) the 2:1 cuts are gain = prereg minima (3 species, 2 lineages, simulated FDR 0.05) and loss = 2 species, margin 2, 3 lineages, which has ~0 simulated recall (D21: losses not FDR-supported, audit only). Source: results/partB/sim/measured/d11_primary_known_gt10/thresholds.json | the earlier D7 cut is out of date after the measured-noise re-derivation | none (record correction); losses already audit-level |
| D26 | 2026-10-30 15:43 | (Part A PA19) F_conv,func headline uses family-pair δ, matching the headline weighting: 0.081 [0.007, 0.178]. Sensitivity with specific non-metal, non-Tanimoto ligand matches only: 0.021 [0, 0.062]; 0 at pLDDT ≥ 90. 0.08 is reported as an upper bound dominated by nucleotide/metal sharing. PA20–22 are descriptive (artifact-flag rule written in code before the runs; a pair belongs to a half only if both groups do; filename fix) | verifier 7327c7a: δ was per pair, not per family pair | headline nearly unchanged; the sensitivity is much lower and is stated |
| D27 | 2026-10-30 16:03 | G5 scope (03:00): rerun only the stages that see new data (snapshot, sites, web index, O2 SOGs reusing the MSA cache, O2 benchmark, D20 split, O2 rescue pool, Part B on O2), plus a reproducibility check of O1/Part A/Part B O1 on identical inputs if time allows. Tier 3 should be complete around 03:00 (about 226M structures). Storage guard: generated data is now 210 GB (/tmp/pe_big 202 GB + repo 7.6 GB); stop and prune regenerable work/ dirs if the total passes 450 GB | a full rerun takes about 5.8 h serial / 3.5 h parallel and would end after 06:00; only tier 3 changes | the second snapshot tests O2-based results only; O1-based results get a reproducibility check, not new data |
| D28 | 2026-10-30 16:13 | Misassigned-orthology check made explicit and symmetric. An O1 candidate fails it if most (>50%) of the outgroup A members it rests on are in a different O2 root HOG from the derived clade. An O2 candidate fails if its outgroups are split from the derived clade by PANTHER family (independent of Ancestra). The co-fold pool (D25) draws from O1 and O2 candidates that pass the full artifact cascade, deduplicated by protein set. G4 freeze moves to 18:30 so that O2 candidates get the cascade | targeted rescue showed that 79% of outgroup absences behind O1 gains are not in the derived clade's O2 family, so O1 gains largely reflect over-merged O1 families; recovery on O2 is 2/50 | removes most O1 gains; the case pool shifts to O2 (more proteomes, finer families) |
| D29 | 2026-10-30 17:02 | (a) Part A headline is still per family pair, as preregistered, but every report gives it next to the per-pair estimate and states that it depends on orthology granularity: F_analog per family pair is 0.819 on O1 and 0.935 on O2, while per pair it is 0.367 on both. Per pair is called the granularity-robust statement; the headline is not switched post hoc. (b) Descent M1 applies D11 symmetrically: a P member counts as informative only if its footprint is > 10 residues, matching A, on both O1 and O2. Both the asymmetric and the symmetric versions are reported | O2's larger groups pack homolog pairs into 3.6× fewer family-pair cells, which moves the weighted fraction; the verifier flagged the M1 asymmetry under D11 | (a) none on the estimates, it is a framing choice; (b) lowers O2 M1 from its inflated 0.93 |
| D30 | 2026-10-30 17:11 | G4 freeze at 17:12 (ahead of 18:30, because the inputs are final): the co-fold pool is case_list_v2 (5e5a3ac), 6 gains (5 top-k + 1 audit) + 1 loss + 5 bad-model false-loss controls; 2 exploratory duplication gains, outside the headline. Headline: 0 convincing. The PANTHER side-family check added by the case builder after seeing results is reported only and does not enter selection. The annotation-first validation track (DESIGN §5.3 control list: steroid receptors, NR4A, LDH/MDH, IDH NADP→NAD, GID1/HSL, D14/KAI2, lysozyme/α-lactalbumin, uricase/GULO, M. leprae truncations) was missed: it is now built as a recall/validation table and supplies the co-folding deep dives with known answers | inputs final; DESIGN §5.3 requires the validation table; known-answer cases give co-folding deep dives with an experimental ground truth | none on the discovery pool; validation shows recall, not new discoveries |
| D31 | 2026-10-30 17:24 | The annotation-first validation recovers 0 of 7 known events, and all 5 negative controls behave. The dominant cause is design scope: SOGs and events are within one orthogroup, but most documented specificity and loss events (MR/GR, PfLDH/MDH, IDH, GID1/HSL, Nurr1 on O2) are between paralogous orthogroups. This is reported as the headline limitation of Part B. An exploratory cross-orthogroup extension (structure-aligned site matching between OGs in the same TED/Pfam family) is built and tested on the validation items first, then genome-wide if it works; it is labelled not preregistered | validation_table 75f2cf5 | the extension can only add candidates; it is reported separately with its own null |
| D32 | 2026-10-30 18:07 | The cross-OG extension, with its method committed before any output (77db16c), recovers 0/7 known events on O1 and O2, so the ≥3/7 gate fails and there is no genome-wide run. Site linking across OGs is specific (610 real vs 1 shuffled links in 2,000 same-family pairs), but the specificity test (a) is too conservative under the global matched-identity null, and the presence test (b) fires more under shuffle (10 vs 33). No v2 is tuned on the 7 items, since that would overfit the validation set. The site-specific statistic (D divided by the larger within-OG spread; ER vs ketosteroid receptors in the top 1.5%) is recorded as a follow-up to preregister for the full-atlas run | a gate failed as preregistered | none; negative result reported |
| D33 | 2026-10-30 23:22 | (Missing entry, found by verifier 7ead1de.) On O2, orthology_split inside pass_artifacts is the DIAMOND RBH rule (0cf2502), not D28's PANTHER rule, which can never fire on O2 because all gains are 'unknown'; PANTHER is reported only. As written, D28 would pass 34 O2 gains; RBH passes 17. Interpretation corrected: in 88% of RBH failures, the outgroup's best hit is in the SAME root HOG but outside the derived clade. So the check detects comparisons against a different paralog subfamily (gains on paralog branches), not misassigned HOGs. O1's 'over-merged families' reading stands only for the root-HOG recount (734/923 = 79.5% of outgroup supporters are in a different O2 root HOG) | RBH rule chosen after O2 orthology came back 'unknown' (post hoc, with bias toward fewer cases); wording corrected per verifier | O2 pool 17 vs 34 gains: conservative |
| C-D31 | 2026-10-30 23:22 | Correction to D31: the 0/7 validation misses do NOT all come from events lying between paralogous OGs. Per item (verify_G4b B): OG split explains MR/GR on O1 (same O2 root HOG HOG0035184), IDH and GID1/HSL. D14/KAI2, lysozyme/α-lactalbumin, Nurr1 (O1), ER vs AR/GR (O1) and PfLDH vs PfMDH are within one OG, and are missed for other reasons: no A calls on the absent side (quality and length rules), the pocket detected in both states, or no specificity test on a derived-only branch or one that is not significant. On O2, the human GULO rows are P/U paralogs, not G | verifier 7ead1de | REPORT wording changes; the conclusion that recall is 0/7 stands |
| D34 | 2026-10-30 23:24 | Specificity: the preregistered BH call gives 0 on O1 and O2, but at the residue-gate + fingerprint-shift stage real data exceed the pseudo-clade null 6.5× on O1 (513 vs 71–88) and 7× on O2 (536 vs 66–84), an empirical null FDR of about 0.15. That is the same estimator used for gains (D24). BH is unattainable by construction: the smallest empirical p (1.8e-5) is above the BH threshold for a single call (about 0.05/10,584 tested ≈ 4.7e-6 on O1). Reported: (a) prereg 0 calls (headline); (b) a population-level specificity signal with empirical FDR; (c) a sensitivity using a GPD/normal tail fit beyond the empirical range, as in D17, labelled post hoc. Funnel-stage pairs go through the artifact cascade as exploratory specificity candidates | the specificity null had too few draws to reach BH; the population-level enrichment is the same real-vs-null estimator as for gains | (b)/(c) can only add; they are labelled post hoc and kept out of the convincing headline |
| D35 | 2026-10-31 22:04 | G5 re-runnability fixes (b83541b, 97e031f): frozen and species-level inputs (rescue_k for D9, TimeTree, observation calibrations) are looked up under the prereg snapshot, with provenance, instead of the run snapshot; the O1→O2 agreement step runs after both Part B lanes; Part A numbers are snapshot-keyed and 11_numbers.py refuses to overwrite partA for another snapshot. Also fixed: the stage-04b sequence DB dropped 7,028 proteins on new snapshots; rbh was missing from the new-snapshot O1 stage list; run_artifacts.py failed on rerun | the first unattended G5 run failed on snapshot-specific inputs and overwrote snapshot-1 numbers | none on snapshot-1 results; enables an unchanged rerun |
