# τ benchmark (DESIGN §8.1), snapshot 20261030T085136Z, orthology O1 (Ancestra tier-1 v0, accurate)

Development subset, preliminary. prereg sha256 `52f94618…9d83` (config/prereg.yaml, frozen). Produced by
`scripts/sog_tau_benchmark.py` on 2026-10-30 (553 s wall, 1 core). Numbers come from `tau_benchmark.csv`,
`tau_benchmark_funnel.json` and `tau_choice.json` in this directory.

## Choice
**τ = 0.2**, the maximum F1 on the TUNE half (prereg `sogs.tau_column_jaccard.objective`). The TEST half was not
used for the choice. MSA method: FAMSA only (see "Not done"). Leiden resolution 1.0 (frozen).

| τ | tune recall | tune precision | **tune F1** [95 % CI] | test recall | test precision | **test F1** [95 % CI] |
|---|---|---|---|---|---|---|
| **0.2** | 0.693 | 0.745 | **0.718** [0.661, 0.774] | 0.591 | 0.803 | **0.681** [0.604, 0.714] |
| 0.3 | 0.621 | 0.769 | 0.687 [0.613, 0.756] | 0.574 | 0.833 | 0.680 [0.589, 0.716] |
| 0.4 | 0.552 | 0.810 | 0.656 [0.588, 0.735] | 0.537 | 0.868 | 0.663 [0.556, 0.709] |
| 0.5 | 0.457 | 0.825 | 0.588 [0.536, 0.642] | 0.442 | 0.878 | 0.588 [0.513, 0.629] |

Primary subset: metal-only pairs excluded, both annotation sources pooled, recall conditional on detection.
The CIs come from a bootstrap over root OGs (200 replicates, seed rule `seed*1000+replicate`, seed 0).
Denominators, primary subset: tune 2,170 detected positive pairs (5,242 including undetected) and 3,677 negative
pairs in 529 OGs; test 5,991 detected positives (9,852) and 5,897 negatives in 580 OGs.

**Caveats for the PI**
1. τ = 0.2 is the **lowest value in the grid**, and tune F1 falls monotonically with τ. The optimum may lie below
   the grid. Extending the grid would be a prereg deviation, so I did not do it.
2. On the test half, F1 is flat between τ = 0.2 and τ = 0.3 (0.681 vs 0.680, overlapping CIs). Going from 0.3 to 0.2
   buys recall at the cost of precision: the test false-merge rate rises from 0.117 to 0.148.
3. The CIs are wide, and the tune-half differences between neighbouring τ lie within them.

## Full grid by subset (full table: `tau_benchmark.csv`)
| subset | half | F1 τ=0.2 | 0.3 | 0.4 | 0.5 | n pos detected / neg |
|---|---|---|---|---|---|---|
| primary (no metal-only) | tune | 0.718 | 0.687 | 0.656 | 0.588 | 2,170 / 3,677 |
| primary | test | 0.681 | 0.680 | 0.663 | 0.588 | 5,991 / 5,897 |
| with metal-only pairs | tune | 0.666 | 0.642 | 0.617 | 0.552 | 2,918 / 5,298 |
| with metal-only pairs | test | 0.633 | 0.634 | 0.638 | 0.562 | 9,021 / 8,518 |
| UniProt ECO only | tune | 0.667 | 0.665 | 0.687 | 0.609 | 538 / 1,786 |
| UniProt ECO only | test | 0.593 | 0.620 | 0.702 | 0.609 | 1,511 / 2,331 |
| BioLiP2 only | tune | 0.706 | 0.677 | 0.648 | 0.580 | 2,001 / 3,671 |
| BioLiP2 only | test | 0.672 | 0.672 | 0.654 | 0.575 | 5,487 / 5,883 |
| identity [0, 0.3) | tune | 0.686 | 0.654 | 0.625 | 0.552 | 1,026 / 2,258 |
| identity [0, 0.3) | test | 0.687 | 0.692 | 0.681 | 0.606 | 4,224 / 3,911 |
| identity [0.3, 0.5) | tune | 0.738 | 0.713 | 0.683 | 0.613 | 930 / 1,297 |
| identity [0.3, 0.5) | test | 0.620 | 0.600 | 0.580 | 0.519 | 1,274 / 1,739 |
| identity [0.5, 1] | tune | 0.801 | 0.744 | 0.691 | 0.653 | 214 / 122 |
| identity [0.5, 1] | test | 0.791 | 0.777 | 0.719 | 0.602 | 493 / 247 |

The UniProt-only subset points to τ = 0.4 on both halves. BioLiP2 contributes about 90 % of the pairs and
therefore drives the pooled optimum.

## Definitions (implemented in `scripts/sog_tau_benchmark.py`)
- **Annotated sites (experimental only).**
  - UniProt `Binding site` features with ECO:0000269 or ECO:0007744 and a ChEBI id, grouped per (accession, ChEBI,
    ligand label). UniProt release 2026_03.
  - BioLiP2 (2025-11-12) sites grouped per (accession, CCD id). PDB entries are clustered by residue overlap
    coefficient ≥ 0.5, and a cluster's site is made of the residues seen in at least half of its entries.
  - Peptide, DNA and RNA ligands are excluded.
  - Annotations on proteins whose model sequence ≠ the current UniProt sequence are dropped (591).
  - In total: 68,660 annotated sites (BioLiP2 62,498, UniProt 6,162) in OGs that went through the pipeline.
- **Detection** (prereg detectability preflight, overlap ≥ 0.5): 42,602 / 68,660 (62.0 %) of annotated sites are
  covered by an atlas site, 61.8 % for BioLiP2 and 64.6 % for UniProt. Recall is reported conditional on detection.
  Unconditional recall at τ = 0.2 is 0.287 on tune and 0.359 on test.
- **Pairs.** The unit is a pair of distinct atlas sites in members of different species within one root OG.
  - Positive: the two sites carry annotations of the same ligand (same ChEBI, or the same CCD id) from the same
    source.
  - Negative: the two sites carry annotations of different ligands whose annotated residues map to disjoint MSA
    column sets.
  - TP: a positive pair whose two sites fall in the same SOG (Leiden membership). FP: a negative pair in the same
    SOG.
  - precision = TP / (TP + FP), recall = TP / detected positives.
- **Split.** prereg `tune_test`: the key is the smallest member accession, and sha1 % 2 decides the half.
  Control-list families are forced into test (13 OGs).

## Annotation disjointness (prereg tune_test.annotation_disjointness)
`tau_tuning_annotations.parquet` lists the 33,528 annotated sites in TUNE-half OGs that were used for tuning
(BioLiP2 30,517, UniProt 3,011; 20,176 detected). Evidence filter 2 must exclude these for the same SOG. The
overlap count belongs in the Part B funnel.

## Not done / open
- **FoldMason** (prereg `msa_choice_rule`) is not installed and was not compared. That decision is still open for
  the 1,633 low-identity OGs (identity < 0.25). FoldMason is the second item in the DESIGN §12 cut order.
- **Leiden resolution** is frozen at 1.0, so the grid [0.5, 1.0, 1.5] was not run.
- **Benchmark scope.** Only the tier-1 (46 proteomes) annotations count; Swiss-Prot tier 2 is not in O1. Re-run on
  O2 by passing `--sogs /tmp/pe_big/sogs/<snapshot>/O2`.
