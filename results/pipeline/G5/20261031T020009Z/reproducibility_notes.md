### Reproducibility notes (O1 side: identical inputs)

Rerun of the 3 failed stages and the stages downstream of them, 2026-10-31 19:51–21:01 UTC, code b83541b plus the
follow-up fix commit (logs `rerun_*.jsonl`). The stage table above still lists the first-run failures (lanes laneA/laneB); the reruns are
the `rerun_*` lanes.

- **O1 SOGs:** every file in `results/sogs/<snap>/O1_*` (web states, SOGs, site→SOG table) is byte-identical to
  O1_r3 (`checksums.sha256`). Only the timing columns of `og_stats.parquet` (`msa_s`, `graph_s`, `og_s`) differ.
- **Part B O1:** the pre-rescue events (called 898, candidates 7,173) and the rescue are identical: k = 1 (frozen,
  D9), 548 SOGs, A→W 3,338, P→W 131. BURROW was recomputed from scratch here (1,601 s) and reproduces the reused
  run exactly. The candidates are identical (932 rows / 560 events, Jaccard 1.0), and so is the artifact cascade
  (7 gains + 1 loss pass). The first attempt gave 9 + 2. That came from a driver defect, not from
  nondeterminism: the new-snapshot O1 stage list left out `rbh`, so `orthology_split` (D28) was `unknown` for
  all 67 primary candidates. This is fixed: `rbh` now runs before `artifacts`.
  Residual: the DIAMOND RBH fractions differ in 2 of 154 rows (`frac_outgroup_rbh`) and 7 of 149 rows
  (`frac_arch_match`), with no change in pass/fail. Snapshot 1 ran `rbh_check.py` with 16 threads over O1_r3 and
  O2_r1 queries together. G5 ran it with 32 threads over O1 only. DIAMOND `--max-target-seqs 1` hits depend on
  the query batch and the thread count (tie-breaking).
- **Part A:** F_analog, F_conv,func and every other `numbers_partA.json` value equal snapshot 1. The first G5
  attempt did not: F_analog per pair was 0.3613 against 0.3671, and per family pair 0.8139 against 0.8186. The
  cause was a stage-order defect. `04b_tail_audit.py` ran before `06_classify.py` and built the shared
  `<work>/mmseqs/seqdb` from UniProt sequences only, which dropped 7,028 proteins without a UniProt sequence.
  06 then reused that DB: 87,796 significant pairs became `unknown`. Fixed with one shared builder
  (`_common.build_seqdb`). The rebuilt `all.fasta` is byte-identical to the one from snapshot 1.
- **M1:** the reported snapshot-1 M1 (0.7736 / 0.7803) uses the pre-D11 `O1` states. G5 builds O1 with the O1_r3
  rules, so its M1 (0.9278 / 0.9344) is compared with snapshot 1 recomputed on O1_r3 states. That gives the same
  point estimates, and the descent pair table has identical content. The bootstrap CIs differ in the 4th decimal
  because `build_descent_pairs.py` writes rows (and `pair_id`) in worker completion order, and
  `16_descent_2d.py` draws OG bootstrap weights over OG codes in table order. This is row-order dependence, not a
  data difference. The same script on the same table reproduces the reported snapshot-1 M1 and CI exactly.
- **O2:** the differences are expected because tier 3 has new data (80,304 vs 7,841 units). The O1→O2 agreement
  now uses this snapshot's O1 run: 2 of the 50 O1 gains are re-called on O2, the same as snapshot 1.
- **cases:** rc 0, 0 cases selected (not compared with snapshot 1; that run used different false-loss and convergence inputs). Case dirs are written under the bulk root.
  False-loss controls are off on a new snapshot (they need a pre-D11 SOG release).
