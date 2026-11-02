# Co-folding negative controls: the exchange rate

What a case-level Δ in a co-folding readout is *worth*, measured on the four
positive-control systems before any real case is folded (DESIGN.md §6; round-1 review
changes 10 and 12; `config/prereg.yaml` → `cofold`).

Scripts: `scripts/cofold/site_ablation.py`, `run_controls.py`, `cofold_readouts.py`,
`apo_describe.py`. Tables: `results/cofold/controls/{summary,exchange_rate,apo_burrow}.tsv`.
Full AF3 outputs: `/tmp/pe_big/cofold/controls/` (node-local, regenerable).

## 1. Design

| | |
|---|---|
| Systems | CDK2 + ATP/Mg (P24941, 1HCK), CA2 + Zn/AZM (P00918, 3HS4), PYGL + PLP (P06737, 1FA9), yeast ADH1 + NAD/2 Zn (P00330, 4W6Z) |
| Site | every holo polymer residue with a heavy atom within **4.5 Å** of the crystal ligand, mapped onto the UniProt sequence by pairwise alignment (100 % identity over the observed chain in all four; 13/22/14/36 residues, all mapped) |
| Ablations | all site residues → **Gly** (Masters 2025 style, strips the lining) and → **Phe** (fills the cavity with bulk) |
| Decoy | caffeine (CCD `CFF`) replacing the cognate organic ligand, structural metals kept |
| Apo | wild type, protein only, then `pocketcompare describe --af-mode` on the refolded model |
| Settings | `--no-templates`, 3 seeds × 5 samples, 1 RTX 4090 per system |

**MSA.** WT, Gly and Phe share one MSA *by construction*: the cached AF3 jackhmmer MSA
with the query row replaced (`make_af3_input.py --reuse-msa-for-variant`). The ablation
Δ is therefore purely sequence-driven and cannot be an MSA-depth artifact. This is a
deliberate, recorded departure from `cofold.msa`, which bars query-row reuse *across the
derived/ancestral boundary*: an ablation is not an ancestral state, and holding the MSA
fixed is the **stricter** test, since the MSA still carries the intact family's
coevolution signal. The toy case in §5 repeats the CDK2 Gly ablation with an independent
ColabFold MSA and reaches the same verdict.

**Readouts** (preregistered; primary = ligand pLDDT). `site_centroid_dist` here is
model-internal — the distance from the predicted ligand centroid to the centroid of the
site residues' Cα atoms *in the same model* — so it is defined for variants and needs no
holo superposition. It is **not** on the same scale as the holo-superposed centroid
distance in `results/cofold/bench/` (a perfectly placed ligand sits ~1–4 Å from the Cα
centroid, not ~0 Å). Within a system it is still directly comparable across rows.

## 2. Results

Top-ranked sample of 15; `sd` over all 15. `recall` = fraction of site residues within
5 Å of the predicted ligand. PB = PoseBusters `dock` validity.

| System | Variant | Ligands | lig pLDDT | P–L PAE min | centroid (Å) | recall | ipTM | PB |
|---|---|---|---|---|---|---|---|---|
| **CDK2** | WT | ATP+MG | **92.1** ±0.4 | **0.84** | 1.8 | 1.00 | 0.96 | ✓ |
| | site→Gly | ATP+MG | 49.4 ±2.2 | 2.68 | 1.44 | 0.68 | 0.76 | ✗ clash |
| | site→Phe | ATP+MG | 63.3 ±4.6 | 2.62 | 14.21 | 0.14 | 0.80 | ✗ clash |
| | decoy | MG+CFF | 72.0 ±2.0 | 1.55 | 3.56 | 0.55 | 0.86 | ✓ |
| **CA2** | WT | ZN+AZM | **97.3** ±0.2 | **0.78** | 3.69 | 1.00 | 0.98 | ✓ |
| | site→Gly | ZN+AZM | 66.7 ±6.1 | 3.81 | 21.20 | 0.00 | 0.79 | ✓ |
| | site→Phe | ZN+AZM | 55.3 ±2.7 | 2.96 | 7.29 | 0.38 | 0.76 | ✗ clash |
| | decoy | ZN+CFF | 79.4 ±4.1 | 2.07 | 7.17 | 0.46 | 0.83 | ✓ |
| **PYGL** | WT | PLP | **97.5** ±0.1 | **0.90** | 1.82 | 1.00 | 0.99 | ✗ (covalent) |
| | site→Gly | PLP | 77.9 ±4.7 | 3.24 | 31.19 | 0.00 | 0.88 | ✗ clash |
| | site→Phe | PLP | 77.7 ±3.3 | 2.47 | 30.76 | 0.00 | 0.91 | ✗ clash |
| | decoy | CFF | 55.7 ±8.3 | 4.71 | 20.23 | 0.00 | 0.80 | ✓ |
| **yADH1** | WT | NAD+2ZN | **96.8** ±0.1 | **0.85** | 1.06 | 1.00 | 0.97 | ✓ |
| | site→Gly | NAD+2ZN | 57.7 ±3.6 | 3.14 | 6.53 | 0.44 | 0.77 | ✓ |
| | site→Phe | NAD+2ZN | 49.6 ±4.7 | 4.16 | 8.57 | 0.39 | 0.75 | ✗ clash |
| | decoy | 2ZN+CFF | 71.7 ±2.4 | 2.74 | 8.76 | 0.31 | 0.84 | ✓ |

### The exchange rate (Δ vs that system's own WT + cognate run)

| System | Perturbation | Δ lig pLDDT | Δ PAE min | Δ centroid (Å) | Δ ipTM | Δ recall |
|---|---|---|---|---|---|---|
| CDK2 | →Gly | **−42.7** | +1.84 | −0.4 | −0.20 | −0.32 |
| CDK2 | →Phe | −28.8 | +1.78 | +12.4 | −0.16 | −0.86 |
| CDK2 | decoy | −20.1 | +0.71 | +1.8 | −0.10 | −0.45 |
| CA2 | →Gly | −30.6 | +3.03 | +17.5 | −0.19 | −1.00 |
| CA2 | →Phe | **−42.0** | +2.18 | +3.6 | −0.22 | −0.62 |
| CA2 | decoy | −17.9 | +1.29 | +3.5 | −0.15 | −0.54 |
| PYGL | →Gly | −19.6 | +2.34 | **+29.4** | −0.11 | −1.00 |
| PYGL | →Phe | −19.8 | +1.57 | **+28.9** | −0.08 | −1.00 |
| PYGL | decoy | **−41.8** | +3.81 | +18.4 | −0.19 | −1.00 |
| yADH1 | →Gly | −39.1 | +2.29 | +5.5 | −0.20 | −0.56 |
| yADH1 | →Phe | **−47.2** | +3.31 | +7.5 | −0.22 | −0.61 |
| yADH1 | decoy | −25.1 | +1.89 | +7.7 | −0.13 | −0.69 |

## 3. What this means for case-level claims

1. **Destroying the site is worth 20–47 ligand-pLDDT points and 1.6–4.2 PAE units.**
   That is the scale a genuine presence/absence difference should reach. A derived-vs-
   ancestral Δ of a few pLDDT points is noise: seed-to-seed SD on a WT control is
   0.1–0.4, but on a perturbed system it rises to 2–8, so **report the SD over all 15
   samples, not just the top-ranked value**.
2. **The decoy baseline is 18–42 pLDDT and 0.7–3.8 PAE — the same order as the
   ablation.** This is exactly why the prereg decision rule is relative: a state Δ only
   counts if it *exceeds* the within-protein cognate-vs-decoy gap. On three of four
   systems ablation > decoy, so the rule has headroom; on PYGL it does not (below).
3. **ipTM is useless here, again.** Every perturbation, including total site
   destruction, keeps ipTM ≥ 0.75 — Δ only 0.08–0.22. Descriptive only, as frozen.
4. **Placement does not discriminate, and can invert.** CDK2+ATP with all 22 site
   residues mutated to Gly puts ATP *closer* to the site centroid than the WT run
   (1.44 vs 1.80 Å) at ligand pLDDT 49. AF3 will dock a ligand into an open cleft
   whatever lines it. Never read "the ligand is in the pocket" as evidence.
5. **PoseBusters is a genuinely orthogonal signal.** All four WT runs and all four
   decoys are valid (PYGL excepted, see below); 5 of 8 ablations fail, every one on
   `minimum_distance_to_protein` — the model pushes the ligand into the backbone the
   side chains no longer shield. So PB failure flags a *structurally impossible* pose,
   not a weak binder, and it separates ablations from decoys in the opposite direction
   to pLDDT. Keep both.
6. **PYGL is the cautionary case, and it argues for the full readout set.** The Gly and
   Phe ablations cost only ~20 pLDDT — less than half the caffeine decoy's 42 — so on
   the primary readout alone the ablation would look *milder* than a non-binder. The
   geometry readouts catch it unambiguously: the ablations throw PLP 29–31 Å off the
   site with site recall 0.00, while the WT sits at 1.8 Å with recall 1.00. PLP is a
   small ligand in a very large protein (847 aa), and the 14-residue footprint leaves
   plenty of surface elsewhere, so AF3 relocates it confidently rather than binding it
   badly. **Consequence: a case may not be called on ligand pLDDT alone.** The prereg
   rule already requires pLDDT *and* min-PAE jointly; PYGL shows centroid distance and
   site recall must be reported beside them, and a large centroid shift at unchanged
   pLDDT is itself evidence of pocket loss.
   (PYGL WT fails PoseBusters because PLP is covalent — a Schiff base to Lys680 — but is
   folded non-covalently; `bondedAtomPairs` is not yet exposed. Known caveat, §5 of
   docs/COFOLD.md. PB is therefore uninformative for PYGL in either direction.)

## 4. Ligand-free refold + BURROW re-description

Each WT protein refolded with no ligand (3 seeds × 5 samples), top-ranked model passed
to `pocketcompare describe --af-mode`, every pocket scored against the known holo site.

| System | n pockets | site res | best pocket | BURROW rank | recall | precision | Jaccard | pocket pLDDT |
|---|---|---|---|---|---|---|---|---|
| CA2 | 4 | 13 | `m0` | 1 | **1.00** | 0.52 | 0.52 | 98.6 |
| CDK2 | 4 | 22 | `m0` | 1 | 0.59 | 0.72 | 0.48 | 93.0 |
| PYGL | 6 | 14 | `m0` | 1 | 0.86 | 0.29 | 0.27 | 98.3 |
| yADH1 | 5 | 36 | `m1` | 2 | 0.44 | 0.62 | 0.35 | 95.5 |

Three of four recover the known site as BURROW's **rank-1** pocket from a ligand-free
AF3 model. This is the control that separates "the AFDB model missed the pocket" from
"the sequence does not encode the pocket", and it works.

Two calibration notes for the loss calls in Part B:

* **Recall is not 1.0 even on true positives.** CDK2 0.59, yADH1 0.44. A 4.5 Å
  crystal-ligand footprint and a BURROW cavity are different objects: BURROW finds the
  concave core, so it misses rim residues that touch the ligand but line no cavity
  (precision 0.62–0.72 is the better-behaved number). **Do not set a "site detected"
  threshold at recall ≥ 0.8**; on this panel ≥ 0.4 recall at rank ≤ 2 is what a true
  positive looks like.
* **yADH1 is at rank 2, not 1.** Its 36-residue footprint spans the NAD cleft *and* both
  Zn sites, and BURROW splits them. Any site-detection rule must consider the top few
  pockets, not only rank 1 — consistent with the atlas rescue rule (DESIGN §3.3).
* PYGL's precision of 0.29 is the mirror image: a 42-residue BURROW pocket against a
  14-residue PLP footprint. Report recall and precision, never Jaccard alone.

## 5. Independent-MSA check

`cases/cofold/cdk2_toy.yaml` re-runs human CDK2 (derived) against the CDK2 Gly ablation
(as a stand-in "ancestral" state) with ATP/Mg and a caffeine decoy, using **ColabFold
MSAs fetched separately for both sequences** instead of the shared cached MSA — i.e. the
route a real case will take. MMseqs2 still returns 16,847 sequences for the ablated
query (vs 20,013 for WT), so the ablation does not starve the MSA.

Verdict **`gain_in_derived`**, positive-control gate passed (1HCK ATP RMSD 1.14 Å ≤ 2.0),
and the rule held in all three seeds with a wide margin:

| seed | derived pLDDT | ancestral pLDDT | D pLDDT | decoy gap X | D PAE | decoy gap X | rule met |
|---|---|---|---|---|---|---|---|
| 1 | 90.8 | 46.8 | 44.1 | 11.5 | 2.69 | 1.79 | ✓ |
| 2 | 90.6 | 47.6 | 43.0 | 11.7 | 2.14 | 1.75 | ✓ |
| 3 | 90.8 | 47.4 | 43.4 | 10.4 | 2.35 | 1.74 | ✓ |

So the ablation Δ survives an independent MSA (43 pLDDT here vs 42.7 with the shared
cached MSA), confirming §1: holding the MSA fixed neither inflates nor deflates the
exchange rate for CDK2. Tables:
`results/cofold/cases/cdk2_toy/{case_table,decision_per_seed,msa_source}.tsv`,
`decision.json`.

## 6. Cost

| Job set | Jobs | Bucket | GPU wall |
|---|---|---|---|
| CDK2 (wt, gly, phe, decoy, apo) | 5 | 384 | 294 s |
| CA2 | 5 | 384 | 275 s |
| yADH1 | 5 | 512 | 382 s |
| PYGL | 5 | 1024 | 1397 s |
| **Total** | **20** | | **2348 s ≈ 0.65 GPU-h** (≈ 23 min wall on 3 GPUs) |

The whole control programme is cheap enough to rerun whenever the readout definitions
change, and it should be.

## 7. Post-hoc check on the calibration panel's decoys (recorded, not acted on)

After the threshold was frozen, the case-level decoy selection revealed a flaw in pure
property matching: it returns decoys that are the *same molecule* as the cognate. For the
case ligands it proposed MET -> MED (D- vs L-methionine, ECFP4 Tanimoto 1.00) and
APC -> ACP (both AMP-PCP, 0.79). Those are not non-binders.

Re-checking the frozen panel, **21 of its 48 decoys exceed Tanimoto 0.4** to their
cognate, including VCA/4I1 and GAL/MAN at 1.00, ATP/AGS at 0.88, NDP/NAI at 0.88 and
SAM/SAH at 0.59 — SAH being SAM's own demethylated product, which occupies the same site.

**Direction of the bias.** Some panel "negatives" are really binders. A composite that
calls them binding is scored as a false positive, so the measured specificity of 0.96 is
a *lower* bound on the truth, and the measured sensitivity of 0.42 and the AUCs are
biased *down*. The frozen threshold is therefore conservative in the direction that
matters: a "binds" call is at least as trustworthy as the frozen number claims.

**Not acted on.** Per the standing instruction not to tune further, the panel was not re-run and
`threshold.json` is unchanged. `scripts/cofold/pick_decoys.py` applies the structural
floor (ECFP4 Tanimoto <= 0.4 after property matching, near-cognates logged as rejected)
and is used for all case decoys; `calibration_panel.py` carries the same helper for any
future panel. Any re-run of the panel with the floor applied should be expected to raise
the AUCs, and must be frozen as a new, separately versioned threshold file.
