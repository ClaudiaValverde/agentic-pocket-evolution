# G5 status: snapshot 20261031T020009Z vs 20261030T085136Z

Updated Sat Oct 31 22:03:29 CET 2026 by `scripts/g5_compare.py` (launcher `scripts/g5_launch.sh`, logs `/tmp/pe_big/logs/g5/20261031T020009Z/`). Development subset, preliminary.

Rerun of the failed stages finished Sat Oct 31 22:03:29 CET 2026 (see Reproducibility notes).

Tier 3: 7,841 -> 80,304 units; 68,104,408 -> 689,760,778 pockets.

| lane | stage | status | wall s | note |
|---|---|---|---|---|
| final | og_index | ok | 9 |  |
| laneA | sogs_O2 | ok | 0 |  |
| laneA | sogs_O2 | ok | 2980 |  |
| laneA | sogs_O2 | ok | 253 |  |
| laneA | sogs_O2 | ok | 121 |  |
| laneA | sogs_O2 | ok | 142 |  |
| laneA | bench_O2 | ok | 1263 |  |
| laneA | split_O2 | ok | 62 |  |
| laneA | rescue_O2 | ok | 828 |  |
| laneA | partB_O2 | FAILED rc=1 | 599 | Error(path) FileNotFoundError: <repo>/results/partB/events/20261031T020009Z/O1_r3/tau0.2/primary/rescue/events/events.parquet  |
| laneB | sogs_O1 | ok | 0 |  |
| laneB | sogs_O1 | ok | 473 |  |
| laneB | sogs_O1 | ok | 113 |  |
| laneB | sogs_O1 | ok | 46 |  |
| laneB | sogs_O1 | ok | 73 |  |
| laneB | sogs_O1 | ok | 19 |  |
| laneB | bench_O1 | ok | 825 |  |
| laneB | partB_O1 | FAILED rc=1 | 108 | ^^^^^^^^^ FileNotFoundError: [Errno 2] No such file or directory: '<repo>/results/observation/20261031T020009Z/rescue_k.json'  |
| laneB | partA | ok | 90 |  |
| laneB | partA | ok | 519 |  |
| laneB | partA | ok | 329 |  |
| laneB | partA | ok | 927 |  |
| laneB | partA | ok | 19 |  |
| laneB | partA | ok | 56 |  |
| laneB | partA | ok | 375 |  |
| laneB | partA | ok | 2155 |  |
| laneB | partA | ok | 29 |  |
| laneB | partA | ok | 36 |  |
| laneB | partA | ok | 2 |  |
| laneB | partA | ok | 0 |  |
| laneB | partA | ok | 123 |  |
| laneB | partA | ok | 204 |  |
| laneB | partA | ok | 10 |  |
| laneB | partA | ok | 0 |  |
| laneB | partA | FAILED rc=1 | 1 | ^^^^^^^^^^^^^ FileNotFoundError: [Errno 2] No such file or directory: '<repo>/data/timetree/20261031T020009Z_O1_pairwise.csv'  |
| rerun_cases | cases | ok | 9 |  |
| rerun_cases | cases | ok | 9 |  |
| rerun_partA | partA | ok | 16 |  |
| rerun_partA | partA | ok | 376 |  |
| rerun_partA | partA | ok | 2179 |  |
| rerun_partA | partA | ok | 29 |  |
| rerun_partA | partA | ok | 38 |  |
| rerun_partA | partA | ok | 2 |  |
| rerun_partA | partA | ok | 0 |  |
| rerun_partA | partA | ok | 124 |  |
| rerun_partA | partA | ok | 219 |  |
| rerun_partA | partA | ok | 10 |  |
| rerun_partA | partA | ok | 0 |  |
| rerun_partA | partA | ok | 28 |  |
| rerun_partA | partA | ok | 25 |  |
| rerun_partA | partA | ok | 0 |  |
| rerun_partB_O1 | partB_O1 | ok | 3916 |  |
| rerun_partB_O1O2 | partB_O1O2 | ok | 40 |  |
| rerun_partB_O1_rbh | partB_O1 rbh,artifacts (run_partB.sh direct) | FAILED rc=1 | 126 | run_artifacts.py: ValueError Duplicate column names (cascade columns appended twice; fixed: idempotent) |
| rerun_partB_O1_rbh | partB_O1 artifacts (run_partB.sh direct, x2 idempotency check) | ok | 37 |  |
| rerun_partB_O2 | partB_O2 | ok | 931 |  |
| shared | snapshot | ok | 148 |  |
| shared | sites | ok | 7 |  |
| shared | sites | ok | 754 |  |
| shared | sites | ok | 4 |  |
| shared | webindex | ok | 946 |  |

## Snapshot 2 vs snapshot 1

| quantity | snapshot 1 | snapshot 2 |
|---|---|---|
| O1 sogs | 260991 | 260991 |
| O1 sogs_multispecies | 116583 | 116583 |
| O1 members_with_model | 382896 | 382896 |
| O1 sites_in_run_ogs | 1431296 | 1431296 |
| O1 state P | 1078376 | 1078376 |
| O1 state A | 188247 | 188247 |
| O1 state U | 13391766 | 13391766 |
| O1 state G | 2185119 | 2185119 |
| O1 §8.1 test F1 (τ 0.2) [CI] | (0.681, 0.604, 0.714) | (0.681, 0.604, 0.714) |
| O2 sogs | 463014 | 463318 |
| O2 sogs_multispecies | 192317 | 192473 |
| O2 members_with_model | 1121082 | 1123060 |
| O2 sites_in_run_ogs | 3908023 | 3914880 |
| O2 state P | 2780379 | 2787322 |
| O2 state A | 548061 | 546547 |
| O2 state U | 40745746 | 40845276 |
| O2 state G | 22894474 | 22905561 |
| O2 §8.1 test F1 (τ 0.2) [CI] | (0.703, 0.668, 0.738) | (0.705, 0.67, 0.74) |
| O2 D20 clean-test F1 (τ 0.2) | (0.695, 0.618, 0.775) | (0.698, 0.632, 0.77) |
| O2 rescue pool rows -> W / A | 100967 -> {'A': 50659, 'W': 50308} | 100883 -> {'W': 50703, 'A': 50180} |
| Part A F_analog per pair [95% CI] | 0.3671 [0.3229, 0.4164] | 0.3671 [0.3229, 0.4164] |
| Part A F_analog per family pair [95% CI] | 0.8186 [0.8102, 0.8257] | 0.8186 [0.8102, 0.8257] |
| Part A F_conv,func per family pair [95% CI] | 0.0815 [0.0066, 0.1783] | 0.0815 [0.0066, 0.1783] |
| Part A F_conv,func per pair [95% CI] | 0.0349 [0.0017, 0.0782] | 0.0349 [0.0017, 0.0782] |
| Part A M1 ortholog [95% CI] | 0.7736 [0.7696, 0.7780] | 0.9278 [0.9258, 0.9300] |
| Part A M1 paralog [95% CI] | 0.7803 [0.7749, 0.7866] | 0.9344 [0.9314, 0.9372] |
| Part A M1 ortholog - paralog (identity x time) [95% CI] | -0.0355 [-0.0423, -0.0279] | -0.0202 [-0.0238, -0.0164] |
| Part A M1 ortholog, like-for-like (snapshot 1 on O1_r3 states) | 0.9278 [0.9252, 0.9301] | 0.9278 [0.9258, 0.9300] |
| Part A M1 paralog, like-for-like (snapshot 1 on O1_r3 states) | 0.9344 [0.9312, 0.9371] | 0.9344 [0.9314, 0.9372] |
| Part A M1 ortholog - paralog (identity x time), like-for-like (snapshot 1 on O1_r3 states) | -0.0202 [-0.0237, -0.0178] | -0.0202 [-0.0238, -0.0164] |
| O1 gains after rescue (2:1, τ 0.2) | 50 | 50 (in both: 50, Jaccard 1.0) |
| O1 gains after artifact cascade (2:1, τ 0.2) | 7 | 7 (in both: 7, Jaccard 1.0) |
| O1 losses after rescue (2:1, τ 0.2) | 17 | 17 (in both: 17, Jaccard 1.0) |
| O1 losses after artifact cascade (2:1, τ 0.2) | 1 | 1 (in both: 1, Jaccard 1.0) |
| O1 all FDR/audit candidate rows, all costs and τ (unique events) | 932 (560) | 932 (560; in both 560, Jaccard 1.0) |
| O2 gains after rescue (2:1, τ 0.2) | 104 | 110 (in both: 90, Jaccard 0.726) |
| O2 gains after artifact cascade (2:1, τ 0.2) | 17 | 17 (in both: 16, Jaccard 0.889) |
| O2 losses after rescue (2:1, τ 0.2) | 66 | 67 (in both: 55, Jaccard 0.705) |
| O2 losses after artifact cascade (2:1, τ 0.2) | 6 | 7 (in both: 4, Jaccard 0.444) |
| O2 all FDR/audit candidate rows, all costs and τ (unique events) | 666 (456) | 685 (465; in both 398, Jaccard 0.761) |

Overlap key: sog_id+type+node (same event = same SOG, same type, same tree node). Snapshot-1 comparators: O1_r3, O2_r1 (same rules: D8, D11, verifier); Part A results/numbers.json[partA] vs results/partA/20261031T020009Z/g5.

## Reproducibility of the O1 side

Tiers 1-2 manifest rows (units, file sizes, mtimes, row counts, method config) between the snapshots: tier 1: 162 units, identical, tier 2: 110 units, identical. O1 SOGs, Part A (tiers 1-2 universe) and Part B O1 therefore have identical inputs: every difference in their rows above measures reproducibility (code version, seeds, nondeterminism), not new data. Only tier 3 (O2 members, web index) sees new data.

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


n/a = the stage did not run or was blocked (see the table above).
