# Cross-orthogroup site matching — method (frozen before any validation output)

**Status: exploratory, not preregistered (DESIGN §12 D31).** Written 2026-10-30 17:24–17:31 CET (`date`; first commit 77db16c), before any cross-OG
output was computed; committed before the validation run so the commit order shows it. Any later change to this
file is a deviation and is listed in CROSSOG.md with its time and reason. Code: `scripts/partB/crossog.py`,
toy tests `tests/test_crossog.py`. Results: CROSSOG.md, `results/partB_crossog/20261030T085136Z/`.

## 0. Inputs (read unchanged)
- SOGs and states at the frozen τ = 0.2: O1 = `O1_r3` (D8, D11), O2 = `O2_r1` (D8, D11). Bulk:
  `/tmp/pe_big/sogs/20261030T085136Z/{O1_r3,O2_r1}/{msa,sogs,states_t0.2,inputs}`; flags from `results/sogs/.../web/sogs.parquet`.
- Sites: `inputs/shard_sites` of each run (residues in UniProt numbering = MSA row numbering).
- Fingerprints: pocketcompare 384-d site fingerprints, `/tmp/pe_big/sites/20261030T085136Z/site_fp/tier*/` (rep pocket of the site).
- Structures: the Foldseek DB of AFDB v6 models `/tmp/pe_big/partA/foldseek/afdb_v6` (tier 1/2, 1.10 M entries),
  Foldseek 10.941cd33 (`envs/ancestra/bin/foldseek`), same binary as Part A.
- Families: TED domains (`data/annotations/ted_domains.parquet`, CATH H-level labels only) and Pfam 38.2
  (`data/annotations/pfam_regions.parquet`, family + clan).
- Prereg state parameters from `config/prereg.yaml` via `pocketevo.sog.prereg_params`, plus D11 (A needs > 10 footprint residues).

## 1. Eligible SOGs and representatives
- Eligible SOG: not `small_footprint`, not `empty_footprint`, and P in ≥ 3 species (post-D11 states).
- Representative of a SOG = its P member with the highest footprint pLDDT whose model is in the Foldseek DB
  (ties: lexicographically smallest accession). No such member → SOG not linkable (counted in the funnel).
- Representative of an OG (target of projections, §4) = member with `seq_ok` whose model is in the Foldseek DB and
  has the highest mean pLDDT over the OG MSA columns it covers × (fraction of OG columns covered) — i.e. highest
  sum of pLDDT over covered columns / n_cols; ties by accession.

## 2. Same-family OG pairs
- OG family labels: a CATH H-level TED label is an OG label if ≥ 50 % of the OG's members carry a TED domain with
  that label. If an OG has **no** member with any H-level TED label, its labels are the Pfam families and Pfam
  clans carried by ≥ 50 % of its members (`PF…`, `CL…` labels).
- A pair of OGs (A, B), A ≠ B, is a same-family pair if they share ≥ 1 label (TED with TED, or Pfam/clan with
  Pfam/clan; an OG with TED labels is compared to a Pfam-only OG through Pfam/clan labels computed for both).
- Both OGs need ≥ 1 eligible SOG for links / test (a); for test (b) only the present side needs one.

## 3. Structural alignment and footprint mapping
- Every needed (rep_x, rep_y) pair is aligned with Foldseek `structurealign` on exactly that pair (pairs written
  as a prefilter DB with `tsv2db`, so no prefilter can drop a pair — equivalent to `--exhaustive-search` for
  these pairs), default 3Di+AA local alignment, `-e 10`, output `fident, alntmscore, qtmscore, ttmscore, cigar,
  qstart, tstart`. Pairs without an alignment are unlinkable.
- Structural support: a rep pair is usable only if `alntmscore ≥ 0.5`.
- Mapping SOG footprint F_x (OG_x MSA columns) into OG_y: columns → rep_x residues (columns where rep_x has a
  gap are *unmapped*) → aligned rep_y residues via the CIGAR (unaligned residues are *unmapped*) → OG_y MSA
  columns via rep_y's row. Result: mapped column set M_xy plus n_unmapped.
- Directional column Jaccard J_xy = |M_xy ∩ F_y| / (|M_xy ∪ F_y| + n_unmapped) (same convention as
  `homology.mapped_jaccard`). Link score J = mean(J_xy, J_yx), each direction using that direction's
  representatives (rep(s_x) → rep(s_y) for J_xy, rep(s_y) → rep(s_x) for J_yx = the same alignment read the other way).
- SOGs s_x ∈ OG_x and s_y ∈ OG_y are **linked** if J ≥ τ = 0.2 (the frozen within-OG τ; not re-tuned).
  Tested links: for each eligible SOG, its best partner (max J, ties by sog_id) in the other OG; the union of
  best-partner links from both sides.

## 4. Test (a) — specificity divergence at a linked site
- Members per SOG: up to 10 P members, one per species (highest footprint pLDDT per species), with a fingerprint.
  Needs ≥ 2 members on each side.
- d(i, j) = 1 − cosine(fp_i, fp_j). Cross statistic D = median of d over all cross pairs (i ∈ s_x, j ∈ s_y).
  Within statistics W_x, W_y = median d over within-SOG cross-species pairs of each side.
- Identity of the link = Foldseek `fident` of the two SOG representatives.
- Null: within-OG ortholog pairs = cross-species pairs of P members inside the same eligible SOG (same member
  selection), over all eligible SOGs of the same orthology version; up to 200,000 pairs, seed 0. Pair identity =
  identical / both-non-gap columns of the OG MSA. Identity bins: [0, .2), [.2, .3), …, [.9, 1.0]. A bin with
  < 500 null pairs is merged into the next lower-identity bin (the lowest bin merges upward).
- p_a = (1 + #{null pairs in the link's identity bin with d ≥ D}) / (1 + n_bin). One-sided (divergence).
- Call (a): BH q ≤ 0.05 over all links tested in the run (validation run: its own BH family; genome-wide run:
  its own) **and** D > max(W_x, W_y) (the shift must exceed both OGs' own variability at that site).
- Reported, not gating: residue-level chemistry — number of mapped footprint positions where each OG's members
  are ≥ 80 % identical (column consensus) and the two consensus residues differ.

## 5. Test (b) — presence/absence across paralog OGs
- For each eligible SOG s_x in OG_x and each same-family OG_y, project F_x into OG_y with rep(s_x) → rep(OG_y)
  (§3, alntmscore ≥ 0.5 required). Testable if ≥ 80 % of F_x columns are mapped.
- States of every OG_y member at the projected footprint F' are computed with `pocketevo.states.compute_states`
  (prereg A_primary rules, τ = 0.2, D11 > 10 residues; same code as the SOG pipeline; G not needed).
- Rescue (prereg, D9 k = 1, hit-Jaccard 0.3, the 6 BURROW settings, `pocketevo.rescue`) is run on the A rows of
  every pair that passes the pre-rescue count below; A → W if any setting hits. Rows that cannot be rescued
  (model not extractable / multi-fragment) → U.
- Call (b): OG_x has P in ≥ 3 species at s_x, and OG_y at F' has **no** P or W member and A (after rescue) in
  ≥ 3 species. Direction reported as "present in OG_x, absent in OG_y" (no gain/loss polarity without outgroups).

## 6. Validation (run first; DESIGN §5.3 items in `scripts/cases/validation_controls.yaml`, role `known`, 7 items)
- For each item and orthology version (O1, O2): derived OGs = OGs of the `derived` proteins, ancestral OGs =
  OGs of the `ancestral` proteins (from validation_table.parquet `ogs`). Tested OG pairs = (derived OG, ancestral OG)
  with derived OG ≠ ancestral OG. **No such pair → item out of scope for the cross-OG method in that version**
  (expected for D14/KAI2 and lysozyme/α-lactalbumin, both within one OG); it counts as not recovered.
  Pairs are tested whether or not they pass the §2 family rule (the family rule is reported).
- Anchor = the validation table's site SOG (`sog_id`, chosen in 75f2cf5 by overlap with experimental residues of
  the `reference` protein, without any cross-OG information).
- Recovered (per item, per version):
  - expected specificity: test (a) called on a tested link that contains the anchor SOG (primary). Lenient
    (reported only): (a) called on any link between a derived and an ancestral OG.
  - expected loss: test (b) called with absent side = a derived OG and present side = an ancestral OG (any SOG).
  - expected gain: test (b) called with present side = a derived OG and absent side = an ancestral OG.
  - items with several expected types are recovered if any one is.
- Item recovered = recovered in O1 or O2 (both reported). **Recall k/7.** Threshold to proceed genome-wide: k ≥ 3.
- No parameter in this file is changed after the validation result. If k < 3 the genome-wide run is not made
  (stated in CROSSOG.md), except the null/funnel on a small pair sample for method characterisation.

## 7. Genome-wide run (O1, only if k ≥ 3) and null
- All O1 same-family OG pairs (§2). Both tests as above; (a) BH over all tested links.
- Rescue in the genome-wide run is applied to the A rows of (b) pre-rescue passes (as in §5).
- Null: a seeded (seed 0) random sample of 2,000 same-family OG pairs (all if fewer), run through the identical
  pipeline with the site correspondence shuffled: every rep→rep aligned residue pair (q, t) is replaced by
  (q, ((t − 1 + o) mod L_t) + 1), with one offset o per alignment drawn uniformly from [L_t/4, 3 L_t/4]
  (destroys the structural correspondence, keeps footprint sizes, contiguity and site densities). The (a) cut is
  the real run's BH p threshold (same cut); (b) uses the same counts (rescue applied identically on null passes,
  or, if time does not allow, both real and null reported at the pre-rescue cut).
- Reported: funnel (pairs → aligned → structurally supported → links → tested → called) for real and null; per-pair
  call rates real vs null on the sampled pairs (same pairs, real vs shuffled); empirical FDR ≈ null rate / real rate.
- Caveat stated up front: the shuffled null estimates how often the tests fire without positional homology; for
  (b) it is not a detection-error null (absence at a random position is common), so (b) calls are audit-level.

## 8. Known limitations (stated before results)
- Interface rule and TED two-domain rule (prereg states) are not applied to the new A calls; P-side rescue (P→W)
  not run. Tier-3/L models are not in the Foldseek DB: O2 reps are restricted to tier 1/2 models.
- One structural alignment per representative pair; no MSA-of-MSAs. Paralogous OGs from over-merged or split
  orthology (D28) can make the same biology appear as within- or cross-OG.
