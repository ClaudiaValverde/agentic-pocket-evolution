# Verification G4: Part A v3 / v3.1 (O2 orthology) and case list v2

Verifier: independent agent (did not produce these results). Written 2026-10-30 17:27 CET (`date`).
Branch `pocket-evolution`, HEAD `98fdb04`. Snapshot 20261030T085136Z. I read the inputs only. My scripts are in
`/tmp/verify_g4/`, and the cached network responses (UniProt REST, PDBe graph API, RCSB) are in `/tmp/verify_g4/net/`.
**Runtime:** 17:11:58 to 17:27, about 15 min wall time. The heaviest step was loading `pairs_cls` (~6 s, 7 GB RSS).

## Verdicts

| # | Claim | Verdict |
|---|---|---|
| A1 | F_anc / F_analog per pair and per family pair, recomputed from the stored tables | **PASS** |
| A2 | The family-pair shift 0.819 → 0.935 is a weighting effect; per-pair values are unchanged | **PASS-with-notes** |
| A3 | 15.7 % of pairs change class | **PASS-with-notes** |
| A4 | qm90d11: P is gated by footprint > 10 exactly like A | **PASS** |
| A5 | The O2 test half uses O2_r1_split with crossover families excluded | **PASS-with-notes** |
| B1 | Selection = top-5 by S among pass_artifacts per category + 1 seeded audit, deduplicated | **PASS** |
| B2 | Nothing was backfilled | **PASS** |
| B3 | f1/f2/f3 labels re-derived independently on 3 cases | **PASS** |
| B4 | The annotation_disjointness fix excludes only tune-half benchmark annotations | **PASS-with-notes** |
| B5 | f2 fails for all 7 pool cases | **PASS** (no missed evidence found) |
| B6 | Hand inspection of the 7 pool cases for artifacts | **PASS-with-notes** (4 cases carry flags) |
| T | Tests | **PASS**: 49 passed, 3 skipped (torch not installed) |

## A. Part A v3 / v3.1

**A1. Recomputation.** `/tmp/verify_g4/recompute_F.py` rebuilds F from K_c and E0_c as Σ(K−E0)+ over class groups. It
does not read the stored F columns. I ran it on `v3/estimands_cos.csv`, `v2_1/estimands_cos.csv`, both `strata_cos.csv`
and both Henikoff CSVs. Every F_anc and F_analog matches the stored value to 4 decimals:
- per pair: O1 0.3671 and O2 0.3671;
- per family pair: O1 0.8186 and O2 0.9347;
- test half: O1 0.2847 / 0.8190 and O2 0.1357 / 0.9095;
- pLDDT ≥ 90 and Henikoff rows: every value in the PARTA tables.

I also recomputed K per family pair directly from the significant-pair tables, with my own cell weights
(`/tmp/verify_g4/cells.py`). The result matches the CSV: K_H3j is 487,849 on O1 and 130,264 on O2, and K_N is
2,378,954 on O1 and 2,307,844 on O2.

**A2. Weighting effect.** These are independent counts from `pairs_cls.parquet` and `context.npz`, and they match every
number in the PARTA text:
- H3j cells: 494,695 on O1 → 137,714 on O2 (3.59×);
- within-group cells: 12,641 → 12,945;
- groups with sites: 148,511 → 89,051;
- significant H3j pairs: 1.394 M → 1.305 M.

Per-pair F (0.3671) and per-site F (0.6997) are identical on O1 and O2. Only 51,162 of the 9,757,299 pairs significant
in both runs (0.52 %) switch between ancestry and non-ancestry classes. The shift is therefore the cell partition, not
the pairs.

*Note:* the per-family-pair estimate also moves within tier-1-only pairs (both_tier1: 0.886 → 0.915). So the cause is
O2's coarser groups in general, not only tier-2 placement.

**A3. 15.7 % class changes.** From `class_moves.csv`, excluding `not_sig`: 1,529,928 of 9,757,299 pairs change class
(15.68 %). This matches the summary JSON. The merge on (i, j) is valid because `universe.parquet` is the same file for
both runs (symlink). Two notes:
1. **Misattribution.** The section is titled "move class because of real placement", but 597,236 of the moves (39 %)
   are tier-1-only pairs, which involve no tier-2 protein. Among tier-1-only pairs, 22 % change class, against 13.6 %
   of pairs with a tier-2 protein. Those moves come from the O1 → O2 HOG change, not from placement.
2. **Mislabelled N counts.** "63 pairs out and 69 in" are H12 ↔ N moves only. All N moves among pairs significant in
   both runs are 165 out and 190 in. That is still negligible, so the conclusion stands.

**A4. D11 symmetric (qm90d11).**
- **Rule in code.** `descent.informative_masks` applies `fp_res > 10` to every informative Y row, P or A.
  `build_descent_pairs.py` computes `fp_res = round(footprint_size × aligned_frac)`. This is the same quantity as
  `n_al` in the state rule (`states.py:153`, `n_al > A_min_footprint_residues`).
- **Check on O2_r1 states** (64/1024 shards, then every 16th shard):
  - 0 A rows have fp_res ≤ 10 (minimum 11), so the reconstruction reproduces D11 exactly.
  - 53 % of P rows (87,414 of 165,175) have fp_res ≤ 10. These are the rows qm90d11 now removes.
  - No P or A row has a missing fp_res, so no row falls back to `inf` (ungated).
- **Numbers.** My own M1 from `pairs.parquet` gives O1 0.8885 / 0.8972 and O2 0.8901 / 0.8924 under qm90d11, and O2
  qm90 0.9323 / 0.9363. All match the v3.1 tables.
- **Toy test.** It passes.

**A5. O2 test half.** I rebuilt each site's half independently, from O2 HOG membership, the placement tables and
idgap, against `og_split_o2.parquet`.
- Agreement with `site_half.npy` is 98.99 %. Test sites from crossover OGs: **0**. The D20 table was used and
  crossover is excluded.
- Descent test half: `og_split_o2` reproduces 259,712 / 86,619 pairs and M1 0.8847 / 0.8881 exactly. The raw hash
  without the crossover exclusion would give 495k pairs.
- *Note:* all 32,654 disagreeing sites are `proxy_hog` proteins (PANTHER-proxy fallback) assigned into an O2 HOG root
  group. They inherit that root's half. This is defensible because the group is the bootstrap unit, but the PARTA
  wording "non-HOG groups excluded" does not mention it. It covers 2.7 % of half sites.

## B. Case list v2

**B1 and B2. Selection and backfill.** `/tmp/verify_g4/sel.py` reimplements the selection.
- **Eligible set.** I defined it as: label supported or convincing, `artifacts == pass`, not exploratory. That gives
  14 rows, identical to the stage `eligible` set: O1 6 (5 gains + 1 loss) and O2 8 gains.
- **Deduplication.** No pair of eligible rows shares a present-side protein (Jaccard 0), so dedup removes nothing.
- **Gains.** Top 5 by S, ties broken by case_id: 34532 (S 7), 35156 (6), 25964 (5), 26737 (4), 35403 (4). The seeded
  audit (`default_rng(0)` over the sorted remainder) picks 56135.
- **Losses.** The only eligible loss is 27153, so there is no audit.
- **Result.** This matches the parquet exactly.
- **Backfill.** None. Every selected row has `artifacts == pass` and `supported`. The higher-S rows that are excluded
  (28478 S 40 convincing, 34468 S 66) all fail artifacts. All primary candidates enter the list: O1 50 gains / 17
  losses, O2 104 / 66.
- *Hygiene note:* `cases/` still holds three v1 directories with only a `manual.yaml` (25451, 35363_g0, 35363_g6).
  They are not in v2.

**B3. Independent re-derivation** (`/tmp/verify_g4/rederive.py`). For each case I rebuilt the clade from the HOG tree,
read post-rescue states, wrote my own conservation code over the MSA, and ran an online f2 search. All labels match.

| Case | My re-derivation | Stored labels |
|---|---|---|
| gain_35403 | f1 pass (4 species, purity 1.0, 9 outgroup); conservation 6/11 columns specific vs background 0.514 → pass | f1 pass, f2 fail, f3 pass |
| gain_34532 | f1 pass (7 species, 0.90, 3 outgroup); conservation 2/11 vs 0.136 → pass | f1 pass, f2 fail, f3 pass |
| loss_27153 | f1 pass (3 species, 1.0, 6 outgroup); conservation fail (8/14 vs 0.62); f3 passes via rescue consistency | f1 pass, f2 fail, f3 pass |

**B4. annotation_disjointness.**
- `used_in_tuning` keeps 12,929 of 33,528 tuning-list rows, and all of them are tune-half OGs.
- All 12,929 keys appear in tune-half benchmark pairs, and 0 appear in test-half pairs.
- In the whole list, only gain_35363_g6 flips f2 because of the exclusion, and it fails artifacts anyway. No pool case
  has any excluded hit (`f2_diag_no_exclusion` is fail for all 7).
- *Note:* 2,448 used annotations with no `site_id` (undetected sites) are dropped from the per-SOG keys. They are
  therefore never excluded in primary mode, only in the OG-wide sensitivity. This is a possible under-exclusion
  (leak), with no effect here.

**B5. f2: my own search for experimental evidence** (`f2_local.py`, `f2_online.py`, `rcsb_seq.py`, `lig_dist.py`)
- **Sources.** For every member (P, A and U) I queried current UniProt features with ECO:0000269/0007744 and the PDBe
  `ligand_sites` endpoint (all PDB entries, including lipids and detergents). For one derived member per case I ran an
  RCSB sequence search at ≥ 60 % identity. No case yields site-level evidence at the footprint. The closest
  near-misses:
  - **34532 (GJB2/Cx26).** Phosphatidylethanolamine in cryo-EM entries 7QEQ, 7QEU, 7QEO and 8QA2 contacts 2/11
    footprint residues (76, 150). That is 0.05 of the lipid site, far below the 0.5 overlap threshold. No UniProt
    experimental site lies in the footprint.
  - **35403 (E. coli CS, P0ABH7).** NADH, CMC, OAA and SO4 in 9 E. coli entries all lie outside residues 320–368. The
    outgroup P. aeruginosa structure 6ZU0 has no ligand.
  - **27153 (POLRMT, loss).** IMT1B (7A8P, R4Q) is the source of the protein-level ChEMBL line. Its nearest footprint
    residue is 10.5 Å away (A706), and the rest are 11–22 Å. That is a different site.
  - **35156, 25964, 26737, 56135.** No PDB entry at ≥ 60 % identity for any derived member. Outgroup structures
    (O50008, P82610, the IQGAPs, Q08AM6) have no ligand at the footprint.
- **Verdict.** "f2 fails everywhere" is correct, and the builder missed nothing I could find.

**B6. Hand inspection** (`case_members.py`, `trunc.py`). For all 7 cases: no UniProt fragment flags, all members
resolve, and all orthology checks (O2/RBH/PANTHER/side-family) pass. The flags:
- **gain_26737 (VAC14, Spirurina).** Present-side footprint pLDDT is 81–84, while the absent outgroups are 92–95. The
  pocket exists only in lower-confidence models, which is a quality-asymmetry risk. Conservation fails; f3 rests on
  rescue consistency only.
- **gain_35403 (citrate synthase, Enterobacteriaceae).** Present-side pLDDT is 87 (the flag `present_plddt_lt90` is
  set), while the outgroups are 90–96. The enzyme is a type II hexamer, so an interface pocket is possible.
  Conservation margin is thin (0.55 vs 0.51).
- **gain_56135 (secernin, Protostomia; the audit case).** Three of the four derived P members are 0.83–0.91 of the
  median length. They lack about 40–70 residues upstream of the footprint (the C-termini align). An indel or
  gene-model effect could create the cavity. V4CJ56 (Lottia, 0.97 of median length) is also P, so this is not
  conclusive. The Trichuris member is A.
- **gain_34532 (connexin, Amniota).**
  - Location: the footprint lies in the TM2/TM3/TM4 bundle, in contact with lipid in the cryo-EM entries.
  - Paralog inconsistency: the GJB6 paralogs are mostly U (partial overlap: mouse P70689 0.125), and human GJB6
    O95452 is P only at Jaccard 0.20, exactly at τ. The gain is mainly a GJB2 signal.
  - Within-species split: zebra finch H0ZM33 is A and H0ZM34 is P, with identical footprint pLDDT (93.2955).
- **gain_25964 (IQGAP, leotiomyceta).** Only 3 distant outgroup species (rat, S. pombe, Candida). Derived conservation
  is 0.67 and conservation fails.
- **gain_35156 (MetE)** and **loss_27153 (POLRMT)** look clean. For the loss, the derived members are 3 filarial
  species with pLDDT 93–95.

## Commands
```
python3 /tmp/verify_g4/recompute_F.py O2 v3/estimands_cos.csv O1 v2_1/estimands_cos.csv ...   # A1
uv run --project . python /tmp/verify_g4/cells.py          # A2 (pairs_cls, context.npz)
uv run --project . python /tmp/verify_g4/d11b.py <states> <sogs.parquet>; .../m1all.py; .../desc_half.py   # A4, A5
uv run --project . python /tmp/verify_g4/half.py           # A5 site halves
uv run --project . python /tmp/verify_g4/sel.py            # B1, B2
uv run --project . python /tmp/verify_g4/rederive.py       # B3
.venv/bin/python -I /tmp/verify_g4/f2_online.py; .../rcsb_seq.py; .../lig_dist.py net/pdb/7a8p.cif.gz R4Q ...   # B5
uv run --project . pytest -q tests/test_partA.py tests/test_descent.py tests/test_evidence.py   # 49 passed, 3 skipped, 6.4 s
```
