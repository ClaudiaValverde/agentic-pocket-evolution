# Round-1 review — structural biology / AF artifacts / co-folding

Read: PROMPT.md, DESIGN.md v0, docs/tools_pocket.md, docs/COFOLD.md + results/cofold/bench/{runtime,placement_summary}.tsv,
LITERATURE.md §1/§5/§7. Twelve changes, all in scope, with GPU cost. None needs new data or new tools.

**Verdict.** The skeleton is sound and unusually well defended for a draft: site-level (not pocket-id) matching, MSA-column
footprints, an explicit U state, a rescue step, and a §6 that already knows co-folding is relative evidence only. The AF3 setup
is real and its positive controls reproduce crystal poses (CA2/AZM 0.28 Å, CDK2/ATP 0.74, PYGL/PLP 1.22, yADH1/NAD 0.41). Not
yet defensible: the absence call rests on one confidence threshold with no measured false-loss rate; the specificity test has no
technical noise floor; and derived-vs-ancestral co-folding is confounded by MSA provenance and templates.

## A. Site definition and matching
**1. Guard the Jaccard-0.3 union-find against percolation.** pocketcompare already merges sub-pockets at residue Jaccard ≥ 0.3
(METHOD §1); §3.1 adds a *second* transitive 0.3 merge. Single-linkage percolates in large proteins: a 900-residue cleft becomes
one "site", the sph8 fingerprint saturates (tools_pocket §9), and footprint and specificity readout are both destroyed. Require
overlap coefficient ≥ 0.5 *and* Jaccard ≥ 0.3 (sub-pocket inclusion is what you mean to merge); cap component size (> ~3× median
pocket residue count → split by Leiden as in §3.2); report site-size and max-site/length as QC; tune on §8.1 as for τ. CPU-minutes.

**2. Quantify MSA-method sensitivity instead of sampling it.** §3.2 treats FAMSA vs MAFFT as a "cross-check on a sample", but
below 30 % identity the MSA *is* the positional-homology call. Compute the P/A/U flip rate between FAMSA, MAFFT and the
Foldseek/US-align mapping for *all* case candidates; demote any flipped state to U. CPU ~1 h.

## B. Observation states and the absence call
**3. pLDDT ≥ 70 for an A call is too permissive, and asymmetric the wrong way.** Akdel 2022 — your own cited precedent — shows AF
pockets match experimental performance only above mean pocket pLDDT **90**; below 90 they are larger and generally miss known
sites. Tier-1 eukaryotes sit at mean 73.7 with 38 % of lining residues < 70 (COVERAGE.md): exactly the unreliable regime. And a
*present* call had to clear the atlas filter (lining pLDDT ≥ 70, score ≥ 0.03) while an *absent* call clears nothing. Make primary
A require footprint mean pLDDT ≥ 90 **and** frac(< 70) ≤ 0.2, with ≥ 70 as a declared sensitivity stratum. Expect candidate losses
to drop sharply — report by how much; that attrition is itself the §7 "artifact cascade" result.

**4. Measure the false-loss rate; the simulated-loss control does not.** §5.3(i) checks the pipeline reports U when you destroy a
model — a unit test, not an error rate. Build a held-out known-present set (BioLiP2/LIGYSIS experimental sites in tier-1/2
proteins), run the state caller blind, report false-A rate overall and by pLDDT band and organism. That number turns "we ruled out
artifacts" into "our loss calls have FDR ≤ x". CPU ~2 h, reuses §8.1 infrastructure.

**5. PAE is not available in bulk here — use TED domains as the scale proxy.** /cluster/work/beltrao/shared/alphafold/v6 tars hold
only cif/pdb (verified), so "PAE from the API for case members only" is right but leaves inter-domain footprints unscreened at
scale. Flag footprints spanning ≥ 2 TED (Lau 2024) domains as PAE-required: no such footprint may be called A without a fetched
PAE showing low inter-domain error. Closes the "your pocket is an inter-domain placement artifact" objection.

**6. Add signal peptides, propeptides and TM segments as a hard filter — absent from §5.3.** AFDB models include the signal
peptide, which packs as a spurious hydrophobic cavity, and TM proteins are modelled without a bilayer so lipid-facing grooves
score as pockets. Lineage-specific SP presence is also a classic gene-model artifact that will masquerade as a clean,
tree-consistent gain. Exclude/flag SOGs whose footprint overlaps UniProt SIGNAL/PROPEP/TRANSMEM (DeepTMHMM where unannotated) by
> 30 %. CPU-minutes.

**7. Monomer-only atlas: make the interface and cofactor caveats filters, not prose.** The atlas has no interface pockets
(tools_pocket §9) and heme/NAD/metal sites are described as empty cavities with small metal sites dropped. A "gain" that is really
an interface, and a "loss" that is really a cofactor site, are both live. Annotate each case SOG with UniProt SUBUNIT /
ComplexPortal state plus an AF3 `--copies 2` homodimer fold, and with UniProt COFACTOR / metal BINDING; co-fold cofactor sites
*with* the cofactor present in both states. Budget from the updated runtime.tsv: a 2×847-residue dimer is ~352 s at 23.3 GB,
i.e. at the 24 GB ceiling — anything larger needs `-m unified` and runs ~3× slower (3×847: 1141 s), so cap dimer controls at
~1,700 tokens or accept the unified-memory cost.

**8. Calibrate the rescue rule and apply it symmetrically.** "≥ 2 of 6 settings" is arbitrary and the 6 settings are far from
independent (HMGCR top-5 Jaccard 0.30–1.0 across seeds/presets; only 55 % of AF-mode pockets have a non-AF-mode counterpart at
Jaccard ≥ 0.5). Pick k as the smallest value meeting a preregistered false-A target on the change-4 set, and rescue P calls too (a
P found in 1 of 6 settings is P(weak)). Rescuing only A biases the pipeline toward discovering losses.

## C. Specificity-change claims
**9. Give the specificity test a technical noise floor — the likeliest single rejection.** §5.2 asks whether between-clade cosine
is below within-clade. But the *same site* re-described under a different setting already gives cosine median 0.886 (p10 0.69),
and v4-Foldcomp vs v6 of the same site gives median 0.95 (p5 0.76) (tools_pocket §7, §9). Clade correlates with model version and
with pLDDT here, so a 0.1 between-clade drop is fully explainable technically. (a) Estimate a technical Δcosine null by
re-describing a few hundred proteins across model source, AF-mode and preset, and require the between-clade effect to exceed its
95th percentile; (b) regress out site_volume_A3, footprint pLDDT and model version before the permutation test; (c) make the §5.3
single-source re-description a *precondition* of the test, not a follow-up; (d) strengthen the lining-residue requirement to a
clade-consistent chemistry change at positionally conserved sites. CPU ~1 h.

## D. Co-folding protocol
**10. Preregister the readouts that discriminate — your benchmark already shows which do not.** In placement_summary.tsv the CDK2
decoys (AZM, caffeine) reach ipTM 0.85–0.86 with `frac_in_site` 1.0: ipTM and "is it in the site" do **not** separate binders from
non-binders, exactly as Kim/Fraser 2026 report. What separates cleanly is ligand pLDDT (true 92.1–97.4 vs decoy 58.2–74.6),
minimum interface PAE (0.78–0.89 vs 1.54–4.02) and centroid distance (0.10–0.49 Å vs 4.1–7.7 Å). Fix the case readout to
(lig_pLDDT, min interchain PAE, centroid distance, PoseBusters) with thresholds set from these distributions *before* cases run;
keep ipTM descriptive only. Widen the decoy panel from 4 systems/2 decoys to ~10 decoys matched on MW, logP, HBD/HBA and
rotatable bonds, so the decoy tail is usable. ~10 jobs × 80 s ≈ 15 min on one GPU.

**11. Match MSA provenance and drop templates across every derived/ancestral comparison — the top confound.** AF3 confidence
tracks MSA depth, and as written you would compare a human protein with a cached jackhmmer MSA *and up to 4 PDB templates*
against a plant or fungal ortholog with a ColabFold MSA and none. Any Δ in ligand pLDDT or PAE is then uninterpretable, and if the
derived protein carries a holo template the claim is circular. (a) Build MSAs for **both** members with the same pipeline
(ColabFold for all, even where a cache exists; 1–5 min/sequence), record and report Neff per member; (b) run case comparisons
`--no-templates` as primary — cdk2_atp_mg_notmpl is indistinguishable from the templated run (ipTM 0.96 vs 0.95) — templated as
secondary; (c) never use `--reuse-msa-for-variant` across the derived/ancestral boundary: it imports the extant family's
coevolution signal into the ancestor.

**12. Two missing controls: a ligand-free re-fold, and the ablation you specified but have not run.** (a) *Apo re-fold* — fold
every case member protein-only (3 seeds × 5 samples) and re-run BURROW on the models. This separates "the AFDB v4/v6 model missed
the pocket" from "the sequence does not encode the pocket", the strongest single defence of a loss claim, and the multi-sample
ensemble partly answers the cryptic-pocket objection (§7.7). ~80 s per 300-residue member. (b) *Site ablation* — §6 specifies the
Masters 2025 poly-Gly control but the benchmark never ran it. Run it on the four positive controls **now**, before any case: if
AZM still lands in ablated CA2 at ligand pLDDT 95, that is the exchange rate for every case-level Δ you later report. 4 × 80 s ≈
6 min. Make both standing columns of the case table.

## E. Ancestral sequence reconstruction + folding the ancestor?
**Yes — adopt it for the 2–4 deep dives only, with three guards.** It answers the strongest objection to Part B ("your ancestral
state is a parsimony label on a pocket predictor's output, never a structure"), and it is the Thornton/Ortlund design LITERATURE
§1 already names as the template. Cost is small: Ancestra gives HOGs with no branch lengths and discards its local trees
(tools_ancestra), so build a per-family FAMSA + IQ-TREE tree and run IQ-TREE `--ancestral` (CPU minutes to ~1 h per family);
folding an ancestor is an ordinary ~80–300 s AF3 job. Total ≤ 1–2 GPU-h for 4 families × (ML ancestor + AltAll + 2 posterior
samples) × (cognate + decoy) × 3 seeds. Guards: (i) *ambiguity* — fold AltAll and ≥ 2 posterior samples, accept only if the pocket
call is consistent, refuse ancestors with mean posterior < ~0.8 over footprint positions (Eick 2017); (ii) *consensus bias* —
reconstructed ancestors are consensus-like and hyper-stable, so AF3 folds them confidently regardless; control with a
shuffled-ancestor and with the extant descendant folded under the identical MSA, and never compare an ancestor's absolute
confidence to an extant protein's; (iii) *MSA matching* — use the family MSA with the ancestor as query row for ancestor and
comparators alike, so change 11 holds. Do not extend ASR beyond deep dives: a gene tree per SOG at scale is out of budget and
would not survive the ambiguity guard.

## F. Rejection risks vs coverage
| Objection a reviewer will raise | Covered by |
|---|---|
| "Absence is just your predictor's 30 % miss rate" | rescue §3.3 + change 4 (measured false-A) + change 8 |
| "Low-confidence region" / "inter-domain artifact" | changes 3, 5 |
| "Fragment or broken gene model" / "misassigned ortholog" | §5.3 length + OMArk-like test; PANTHER §2 (eggNOG would be stronger) |
| "One ortholog modelled apo, the other holo" | **not covered** — add AHoJ-DB apo/holo class per case + change 12a |
| "Cryptic pocket, not a lost one" | change 12a ensemble; name PocketMiner as the explicit screen |
| "Needs the dimer / needs the cofactor" | change 7 |
| "That is a signal peptide or a lipid-facing groove" | change 6 |
| "Your cosine shift is model version, not evolution" | change 9 |
| "Co-folding confidence ≠ binding" / "the two folds were not comparable" | §6 framing + changes 10, 12b, 11 |

The one real gap in that table: **AHoJ-DB apo/holo labelling per case pocket** is cited in LITERATURE §1 and §7.6 as the control
for state-driven fake pocket turnover but never enters the §5.3 cascade. Add it — a table join, not compute.
