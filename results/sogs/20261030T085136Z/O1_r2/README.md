# results/sogs/20261030T085136Z/O1_r2: rebuild of O1 (supersedes nothing; O1 = committed G2 run, unchanged)

Same inputs, prereg (sha256 52f94618…9d83) and frozen τ = 0.2 (D8; `../O1/tau_choice.json`) as `../O1`. Code: commit
recorded in every sidecar (`git_commit`, `code_version`). Development subset, preliminary.

Changes vs O1 (PI decision D8 + verifier notes):
1. `sogs.parquet` / `sogs_all_tau.parquet`: `small_footprint` (footprint <= 3 columns, D8: excluded from event candidates)
   and `empty_footprint` (>= 2 sites, footprint empty: excluded from candidates). τ=0.2: 4,507 multi-species SOGs
   small_footprint (6,595 any), 1,189 multi-site SOGs empty_footprint.
2. A thresholds use the conservative pLDDT lower bound (stored floor(2·pLDDT)/2; true means in [90.0, 90.25) can no
   longer pass): A 748,649 -> 735,962 (12,687 -> U `plddt<90`); `footprint_plddt` still reports the midpoint.
   SOGs, P and G are byte-for-byte the same as O1 (P 1,078,376; G 2,185,119; 260,991 SOGs).
3. State rows carry `model_accession` / `model_tier` (O2 ID-gap support; = accession / 1 for O1).
4. `sensitivity_tau0.3/`: τ = 0.3 sensitivity (D8): per-SOG counts (`web/sogs.parquet`) and funnel; full state table in
   bulk `/tmp/pe_big/sogs/20261030T085136Z/O1_r2/states_t0.3/`. τ=0.3: 372,699 SOGs (131,010 multi-species),
   P 1,163,181, A 1,575,138 (s70 2,377,964), U 34,111,765, G 3,382,696.
G2 slice at τ = 0.2: `slice_G2_*.csv` (same SOGs as O1; ESR1/AR SOG A 8 -> 7).
Bulk: `/tmp/pe_big/sogs/20261030T085136Z/O1_r2/` (`BULK_MANIFEST.sha256`). Method/caveats: `../O1/SOGS.md`.
Web app: `SOG_DIRS=results/sogs/20261030T085136Z/O1_r2/web`.
