# verify_G2 — independent verification of the τ benchmark, observation error rates and sites

Verifier: independent agent, fresh context. Date 2026-10-30 (wall clock ~11:15–11:45 UTC). Snapshot 20261030T085136Z.
I edited no producer code and made no commit. Scratch scripts and re-run outputs are in
`<agent-scratchpad>/v/`
(split.py, ktab.py, ctrl.py, numbering.py, posneg.py, negj.py, states_spot.py, frag.py, inv.py, taurerun/).

## Verdicts

| item | verdict |
|---|---|
| (a) SOG τ benchmark and states | **PASS-with-notes** |
| (b) Observation error rates | **PASS-with-notes** for the numbers. **FAIL** for the report: OBSERVATION.md is stale, and the results are not committed. |
| (c) Sites | **PASS** |
| Provenance (checklist §1) | **FAIL**, details under "Provenance" below |

---

## Provenance (checklist §1–2)

- **Prereg hash.** `sha256sum config/prereg.yaml` gives `52f94618…0d83`, which equals RESEARCH_TREE.md and every
  sidecar I checked. PASS.
- **Checks missing.** `scripts/check_prereg.py`, `scripts/check_report_numbers.py` and `results/numbers.json` do not
  exist. OBSERVATION.md cites `results/numbers.json["observation"]` and `OBSERVATION_numbers.md`, and neither file
  exists. FAIL (checklist items 3 and 18).
- **Observation results are not committed.** `git status` shows `?? results/observation/`, and
  `scripts/observation/`, `observation.py` and `rescue.py` are tracked. FAIL ("committed artifacts").
- **Code version of the SOG outputs.** The SOG/τ sidecars carry `git_commit 0e11201…-dirty` and `code_sha256 7416036f…`.
  - That hash matches neither the committing commit ed79fdc/HEAD (`23d99c69…`) nor the current working tree
    (`2a86a66a…`; `scripts/build_sogs.py` has +170/−51 uncommitted lines).
  - So the code that produced the committed SOG/state outputs is not recoverable from git.
- **Bulk overwritten during verification.** While I was verifying, the producer re-ran
  `build_sogs.py --stage inputs,sogs,states,export` in place on `/tmp/pe_big/sogs/20261030T085136Z/O1`. It used the
  uncommitted code: PID 3210287, started 12:18:32 CET, new `sogs_run.json` at 12:26 CET.
  - `sha256sum -c BULK_MANIFEST.sha256` now fails 4,103 of 10,563 files (inputs/, sogs/, states_t0.2/).
  - The committed G2 bulk therefore no longer exists on disk.
  - **Fix:** write new builds to a new run directory (e.g. `O1_v2/`) and never overwrite a manifested bulk. Commit
    the code that produced every committed output (code_sha256 must match a commit).
- **Timestamps.** The decision log (D7–D10 at 15:00–16:15) and the `pi_freeze.date` (13:30 CET) are later than the
  wall clock while I was verifying (~12:40 CET). Every artifact, prereg.yaml included (mtime 10:47 CET), predates
  the stated freeze time. The content hash is unchanged, so no tuning leak follows from this. The ledger times are
  still not real times and should be corrected.

---

## (a) τ benchmark and states — PASS-with-notes

### Tune/test split
I re-derived the split independently (`split.py`, my own orthoXML parse):

- key = min accession over all root-group geneRefs; rule = `sha1 % 2`.
- Over 36,736 OGs there are **0 mismatches** with `inputs/ogs.parquet`.
- All 22 forced accessions were found; they lie in 13 OGs. All 13 are flagged `forced_test` and are in the test
  half.
- Eligible OGs split 12,683 tune / 12,664 test, matching `inputs_run.json`.

The code (`sog.split_key/split_half`, build_sogs.py:281-287) implements the prereg key and rule exactly. PASS.

### τ chosen on the tune half only
- `main()` picks τ from `R[(half=="tune") & (subset=="primary")]`, sorted by f1, then precision, then smaller τ.
- `tau_tuning_annotations.parquet` contains 33,528 rows: 100 % tune-half OGs and 0 forced-test OGs.
- No site appears in pairs from both halves: 0 shared sites, and the half is per-OG by construction.

PASS.

### Numbers reproduce
- I re-ran the full benchmark into `scratchpad/v/taurerun` (667 s).
- `tau_benchmark.csv` matches with **max abs diff 0.0**, bootstrap CIs included.
- `tau_benchmark_pairs.parquet` is `DataFrame.equals` → True, and `tau_tuning_annotations.parquet` has the same md5.
- τ choice = 0.2.
- Recomputing TP/FP/F1 by hand from the committed pairs gives tune F1 0.7179/0.6869/0.6564/0.5880 and test
  0.6808/0.6796/0.6634/0.5878, exactly as in TAU_BENCHMARK.md.

PASS.

### Pair construction
There is no leakage: pairs are within one OG, across different species, and with distinct accessions; there are no
duplicate pairs, and no pair is both pos and neg. There are, however, construct-validity notes:

1. **Positives are not position-checked.**
   - A positive = same ligand key + same source; MSA overlap of the annotated residues is not required.
   - In a seeded slice (60 shards, 47 OGs, 309 detected positive site pairs; `posneg.py`), **10 %** of positives
     have disjoint annotated MSA columns. These are the same ligand at non-homologous positions (SAH, SO4, NDP,
     TTP…), and their recall is 0.03.
   - Recall on positionally overlapping positives is 0.69 (all halves) and 0.80 at τ=0.2 on tune.
   - Effect: recall is understated by roughly a constant amount; the τ ranking is unchanged in the slice.
   - Crystallisation additives (SO4 and others) count as positives.
   - **Fix:** require annotated-column Jaccard > 0 for positives, or report that subset.
2. **Negatives are defined on annotation residues, not on the atlas sites being scored.**
   - In a seeded slice of 150 OGs (`negj.py`), 11.7 % of negatives join two atlas sites whose own column Jaccard
     is ≥ 0.2, i.e. the same pocket with two ligands. These account for **58 % of FPs** at τ=0.2.
   - Precision is therefore understated, and most so at low τ. This would not move τ off 0.2, which is already the
     grid edge (D10).
   - **Fix:** exclude negatives whose atlas sites have column J ≥ min(τ grid), or report a sensitivity row.
3. UniProt (ChEBI) and BioLiP (CCD) keys never match each other. So cross-source pairs can only be negatives, and
   only with disjoint columns. This is harmless.

### State rules vs prereg (states.py)
- A requires: aligned non-gap ≥ 0.8, mean pLDDT ≥ 90, frac < 70 ≤ 0.2, length ratio in [0.7, 1.3] × OG median
  (median over members with a model), no terminal truncation, SP/TM overlap ≤ 0, model status ok, seq_ok, no
  partial site overlap.
- P requires footprint Jaccard ≥ τ and seq_ok.
- This matches the prereg. Notes:
  - **pLDDT quantisation.** plddt_q = floor(2·pLDDT) and the mean is rebuilt with +0.25, so the reported footprint
    mean has error in (−0.25, +0.25].
    - 12,687 of 748,649 A calls (1.7 %) have reported mean in [90, 90.25). Their true mean can be as low as 89.75,
      i.e. below the frozen cut.
    - In my spot checks the deviation was up to 0.12 (e.g. raw 90.36 vs reported 90.45).
    - **Fix:** apply the cut on the lower bound, `sum(q)/n/2 >= 90`, or store pLDDT at 0.01 resolution.
  - SP/TM is only checked where UniProt annotates it. The prereg says "DeepTMHMM where unannotated", and that is
    not implemented (`sptm_checked` records the gap).
  - 1,189 multi-member SOGs at τ=0.2 have an **empty footprint**; all their members are U(aligned<0.8). This is
    conservative, but should be counted in the funnel.

### 20 state calls spot-checked against the raw AFDB v6 mmCIF and the MSA (`states_spot.py`)
- Sample: seed 20261030, 8 shards; 8 A, 6 U, 6 P.
- Each model was extracted from the AFDB tar and read with gemmi (CA B-factor). I recomputed the footprint →
  UniProt positions from the MSA, aligned fraction, mean pLDDT, frac < 70, length ratio, truncation and best site
  Jaccard.
- **20/20 agree** with the stored state and reason: 8 A, 6 P, U for length_ratio, partial_site_overlap, aligned<0.8
  ×3 (one with an empty footprint), and plddt<90.
- Raw and stored sequences were identical in all 20.

---

## (b) Observation error rates — numbers PASS-with-notes; report FAIL

### Re-detection sampling
- I re-ran `sample_redetection.py` into scratch. The accession list and seed subset are **identical** (2,000 / 300).
- Allocation is equal per kingdom × pLDDT band: 245 each, Archaea <70 takes all 39, and the remainder goes to
  Archaea 70–90 (246).
- The seeds subset has 33–34 proteins per stratum.

PASS.

### Method-record hash
The hash check is real:

- `parquet_method_hash` reads the `pocketcompare` schema metadata from the files themselves: atlas pockets file,
  `v4_describe.parquet` and `seed_s*_part*.parquet`. All hash to b110dd…19b1, which equals the hash of the
  hard-coded atlas record.
- The seed files carry `pocketevo_burrow_model` sha256 values that match the actual s1/s2 model files (`8fbafb…`,
  `abff3d…`).
- Caveat: the 11 fields hashed cannot show a seed or model swap (this is stated in OBSERVATION.md).

PASS.

### False-A definition
In `false_a_sites.parquet` I recomputed `false_A == A_primary & best_atlas_cov ≥ 0.5` and
`A_primary == (mean ≥ 90) & (frac70 ≤ 0.2) & no skip & no SP/TM`. Both hold for every row (21,127 false-A of 52,949
A-eligible → 0.399). PASS. Note:

- The state caller does not use coverage. It calls U for any partial overlap and A only when no site overlaps at all.
- The false-A analogous to the state caller's rule is therefore the "any overlap" row, **0.245**, not 0.399.
- The simulations should use the rate that matches the caller's rule.

### Numbering spot checks (`numbering.py`, seed 20261030)
- **BioLiP2, 10/10 correct.** For every site, all residues_pdb one-letter codes equal the v6 model residue at the
  mapped UniProt position. This includes non-trivial offsets: P55776 PDB 36 → UniProt 56; P06005 and P02185 +1.
- **LIGYSIS af_set, 10/10 consistent.** Site residues lie within 5 Å of the transferred ligand atoms in the v6 model
  frame (median 2.2–4.5 Å; ≥ 0.86 of residues within 5 Å in 9/10, 0.50 in 1/10).
- Region pLDDT recomputed from the raw CA B-factors matched the stored values to 2 dp in all 20.

### Rescue control regions (`ctrl.py`)
I regenerated the controls with the producer's seeded procedure:

- The case sample is identical: 2,500.
- The control set is identical: 2,454 ids.
- Recomputed n_hits from the cached BURROW outputs: **0 mismatches**.

On non-site status:

- Controls overlap 0 atlas sites, 0 known sites used and 0 LIGYSIS sites.
- **5.4 % (133) overlap a BioLiP regular or metal site** (1–3 residues) that was excluded upstream (<3 residues,
  n_mapped < 0.8, or deduplicated).
- 11 % overlap any BioLiP entry, including peptide, DNA and RNA interfaces; 0.8 % overlap UniProt binding or active
  sites.
- Effect is small: tune k=1 control-rescued is 0.348 overall and 0.344 on clean controls only.
- **Fix:** forbid all raw BioLiP and UniProt binding-site residues when picking patches.

### k table
I recomputed `rescue_k_table.csv` from `false_a_sites.parquet` and `rescue_hits.parquet` for both rules, all halves
and k = 1..4. **Max abs diff 0.** I also recomputed case n_hits from the cached CSVs on a seeded slice: 0/26
mismatches. Tune k=1: false-A after rescue 0.1285, controls rescued 0.348. No k meets 0.05, so `k = null` and D9
(k=1) is correctly reported as a PI decision.

Notes:

- The case sample is proportional, not equal, allocation. It also stratifies on source, which OBSERVATION.md does
  not mention.
- The rescue hit rule uses Jaccard 0.3; the chosen τ is 0.2.

### Report — FAIL
OBSERVATION.md (11:45) predates the current outputs (12:13) and describes the superseded
`ligysis_only_v1/` run:

- It says false-A is 0.656 with "BioLiP2/UniProt not yet available". The current `false_a.json` gives 0.399 on
  73,107 sites from BioLiP2 + UniProt_ECO + LIGYSIS.
- Its k table (0.175/0.328, …) is not the committed `rescue_k.json` (0.129/0.348, …).
- It reports 1,500 cases / 1,392 controls; the current run has 2,500 / 2,454.
- RESEARCH_TREE mixes the two runs: D9 uses 0.129/35 %, while D7 and the "[!] Interim" line use 0.65/66 % (the old
  LIGYSIS-only number).

Fix: regenerate OBSERVATION.md from the current JSON, label `ligysis_only_v1` as superseded, and correct D7 and the
tree line.

---

## (c) Sites — PASS

- **Fragment rule, independent check.**
  - I extracted 6 seeded accessions (25 fragments) from the AFDB tar and read `_struct_ref_seq` with gemmi. Every
    db_align_beg − auth_beg = 200(k−1).
  - Residue names in UniProt numbering agree 100 % across overlapping fragments (F1 vs F2: 1,200 positions; F1 vs
    F4: 800).
  - All site residues are ≤ the protein length.
  - The committed `fragment_offsets.tsv` has 3,251 rows with 0 offset and 0 window failures.
- **Checksums.** `sha256sum -c checksums.sha256` in `/tmp/pe_big/sites/20261030T085136Z`: all 646 files OK.
- **Invariants** (`inv.py`, 9 units: 8 seeded tier-1/2 plus unit 0 with fragments):
  - the sites partition the pockets;
  - site residues = union of member pockets' residues after the +200(k−1) offset;
  - rep = the max-score member;
  - 0 cross-site pocket pairs with J ≥ 0.3;
  - site_index is ordered by score.

  All hold.
- pytest.log: 17 passed.
