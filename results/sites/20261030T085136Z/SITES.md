# Sites — snapshot 20261030T085136Z (DESIGN.md §3.1)

Built 2026-10-30 by `scripts/build_sites.py` (code: `src/pocketevo/sites.py`, schema `sites-v1`,
Jaccard threshold 0.3, code sha256 `74f576eb…`). Tests: `tests/test_sites.py` (17 pass).
Stats: `scripts/site_stats.py` (CSVs in this directory). All numbers below come from those CSVs.

## Where the data is (not in git: 2.9 GB)
`/tmp/pe_big/sites/20261030T085136Z/` (node-local /tmp; it disappears with the job, so regenerate it with the
commands below, ~10 s for tiers 1–2). `checksums.sha256` here lists sha256 for all 644 parquet files.
- `sites/tier<T>/<unit>.parquet`: one row per site, metadata only (≈ 86 B/site; tiers 1–2: 274 MB).
- `site_fp/tier<T>/<unit>.parquet`: `site_id, rep_pocket_id, fingerprint` (float16[384], copied unchanged
  from the atlas row of the rep pocket; ≈ 730 B/site; tiers 1–2: 2.3 GB).
- `state/<unit>.json`: input key, counts, timings for each unit; `run_*.json`: run summaries.
- `split_accessions.tiers*.parquet`: accessions whose fragments are spread over more than one unit file.
  In this snapshot there is one: P98161, units 1–2. It is processed by unit 1.

```bash
export UV_CACHE_DIR=/tmp/pe_big/uvcache
.venv/bin/python scripts/build_sites.py --manifest data/snapshot/20261030T085136Z/manifest.parquet \
    --tiers 1,2 --out /tmp/pe_big/sites/20261030T085136Z --workers 32
.venv/bin/python scripts/build_sites.py ... --tiers 3 --sample 50 --seed 0      # timing sample
.venv/bin/python scripts/site_stats.py --sites /tmp/pe_big/sites/20261030T085136Z \
    --manifest data/snapshot/20261030T085136Z/manifest.parquet --out results/sites/20261030T085136Z
```

## Site record
`site_id` (`<accession>:s<k>`, where k counts from 0 in order of decreasing `max_burrow_score`, with ties broken
by `rep_pocket_id`), `accession, site_index, taxid, organism, tier, unit_id` (the unit file holding the site),
`model_version, source` (filled for tiers 1–2), `n_pockets, member_pocket_ids` (sorted by decreasing score),
`fragments, n_fragments, residues` (sorted UniProt positions, int32), `n_res, res_min, res_max,
max_burrow_score, mean_burrow_score, plddt_mean_rep/min/max` (over the members' `plddt_mean`), `plddt_min_min,
plddt_frac_lt70_max, volume_sum, volume_max, n_points_sum` (sums double-count duplicates that come from
overlapping fragments), `rep_pocket_id` (highest score), `rep_fragment, rep_site_volume_A3`.

## Fragment numbering rule (verified)
For fragment k, UniProt position = author resnum + 200·(k−1). Fragments are 1400-residue windows with a
200-residue step, and the last window ends at the protein length. `scripts/verify_fragment_offsets.py` read
`_struct_ref_seq` from every fragmented model in the AFDB v6 human tar (seeking through the atlas tar_index):
**3,251 fragments of 215 accessions, 0 offset mismatches, 0 window mismatches** (`fragment_offsets.tsv`). The test
suite also checks P04275 (VWF, 9 fragments, F9 = 1601–2813) directly against its mmCIF. Fragments occur only in
tier 1 (human).

## Results
| tier | units | proteins ok | with ≥1 site | pockets | sites | sites/protein (mean / median / p90 / p99 / max) | multi-pocket sites | pockets in multi-pocket sites |
|---|---|---|---|---|---|---|---|---|
| 1 | 162 (all) | 552,810 | 453,269 (82.0 %) | 1,788,750 | 1,779,029 | 3.92 / 3 / 7 / 14 / 140 | 0.200 % (3,559) | 0.74 % |
| 2 | 110 (all) | 436,886 | 408,716 (93.6 %) | 1,442,678 | 1,442,673 | 3.53 / 3 / 6 / 11 / 31 | 0.0003 % (5) | 0.0007 % |
| 3 (sample: 50 of 7,841 units, seed 0) | 50 | 151,905 | 129,480 (85.2 %) | 454,608 | 454,606 | 3.51 / 3 / 6 / 11 / 38 | 0.0004 % (2) | 0.0009 % |

The sites/protein column counts only proteins that have at least one site. Including proteins with no site, the
means are 3.22 (tier 1), 3.30 (tier 2) and 2.99 (tier 3). Pockets per site, tier 1: 1 → 1,775,470; 2 → 1,334;
3 → 760; 4 → 502; 5 → 348; ≥6 → 615 (max 21). Tiers 2–3 have at most 2. Histograms are in
`pockets_per_site_hist.csv` and `sites_per_protein_hist.csv`.

**Key finding: within one model, merging at Jaccard 0.3 does almost nothing, by construction.**
- pocketcompare already merged sub-pockets at residue Jaccard ≥ 0.3 (`PC/docs/METHOD.md`), so atlas pockets of
  the same model essentially never overlap that much.
- `threshold_sensitivity.csv` takes 40 units per tier, about 1.4 M pockets in total. Its within-model pairs:
  0.25 ≤ J < 0.3: 3,907 pairs; J ≥ 0.3: 8 pairs (maximum 0.333); J ≥ 0.5: none.
- Nearly all of the merging in tier 1 is cross-fragment deduplication. Of the 3,559 multi-pocket sites, 3,549
  span ≥ 2 fragments. The 192 fragmented proteins that have sites go from 17,970 pockets to 8,259 sites.
- So at J = 0.3, a site is a fragment-deduplicated pocket. It does **not** regroup a binding site that BURROW
  split into several pockets (POCKETOME §6 pitfall 6). That would need a lower threshold or a spatial criterion.
- The threshold sensitivity on the same sample gives:

| tier | pockets/site at J=0.1 | J=0.2 | J=0.3 | J=0.5 |
|---|---|---|---|---|
| 1 | 1.182 | 1.031 | 1.015 | 1.011 |
| 2 | 1.178 | 1.014 | 1.000 | 1.000 |
| 3 | 1.184 | 1.015 | 1.000 | 1.000 |

  Going down to J = 0.1 merges about 15 % of pockets, and does so by chaining. Nothing here validates that
  choice: picking the threshold belongs to the benchmark in DESIGN §8, not to these counts.

## Runtime (64-CPU node, 32 workers, atlas on /cluster/scratch)
| run | units | pockets | wall | CPU (sum over units) |
|---|---|---|---|---|
| tiers 1+2, full | 272 | 3.23 M | 5.3 s total (4.5 s processing; pre-pass < 1.1 s; peak RSS 0.43 GB) | 137 s (42.5 s per 1 M pockets) |
| tier 3, 50-unit sample | 50 | 0.45 M | 1.6 s total (0.9 s processing) | 19.9 s (0.40 s/unit; 43.7 s per 1 M pockets) |

Tier-3 estimates from the sample, at 43.7 CPU-s per 1 M pockets (0.40 CPU-s per unit):
- **Current snapshot** (7,841 units, 68.1 M pockets): ≈ 3,100 CPU-s, about 2 min wall.
- **Complete atlas** (80,304 tier-3 units, ≈ 0.8 B pockets): ≈ 35,000 CPU-s ≈ 10 CPU-h, about 20 min wall on
  32 workers. Reading ~7 MB of pocket file per unit (≈ 0.56 TB in total from scratch) probably makes it
  I/O-bound, so 20–40 min is realistic.
- **Storage**: site metadata ≈ 86 B/site, about 70 GB for 0.8 B sites. Rep fingerprints would add about
  560 GB, more than the whole 500 GB budget. For tier 3, run with `--no-fingerprints` and read the rep
  fingerprint from the atlas row `rep_pocket_id`; it is identical by construction.

## Scaling and correctness notes
- Each worker reads one unit file (≈ 1e4 pockets) at a time, fingerprints included. Edges come from a sparse
  pocket × (accession, position) incidence product, so pockets of different accessions can never be joined.
- Accessions split across unit files: a pre-pass reads only parquet footer statistics (fragment max, accession
  min/max) and caches the result. The lowest unit_id that contains the accession is its owner. The owner also
  reads that accession's rows from the other units, and the other units drop them.
- Runs are resumable. The input key covers the manifest's size and mtime for each input file (rechecked on
  disk, so a mismatch fails), the split assignment, the threshold, the schema version and the sha256 of
  sites.py. A rerun with unchanged inputs skips all 272 units in 0.16 s.
- Invariants are checked at build time and in the tests: every pocket belongs to exactly one site;
  Σ n_pockets = number of input pockets; the site's residues are the union of its members' residues; the rep
  is the highest-scoring member; no pair of pockets in different sites of an accession reaches the threshold.
  Two runs give identical tables, including after shuffling the input rows.
- Empty unit files (6 in tier 3) have a truncated schema and are handled: they produce an empty site table.
