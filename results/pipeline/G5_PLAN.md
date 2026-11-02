# G5 rerun plan (written 2026-10-30 16:02 CET, `date`)

Driver: `scripts/run_pipeline.py --snapshot new --run-id g5 [--msa-seed /tmp/pe_big/sogs]` (`--dry-run` prints every
command). Test on the G2 families (current snapshot, /tmp/pe_big/pipeline_test) ran sites, SOGs O1 + O2, the
benchmarks, the D20 split, the O2 rescue pool and the Part B O1 events with rc = 0. A committed run id is refused.

## Atlas growth (atlas_snapshot.py, 2026-10-30 14:47 UTC)
| | 20261030T085136Z | now | 03:00 CET (projected) |
|---|---|---|---|
| tier-3 units | 7,841 | 32,811 (of 80,304) | about 80,300 (complete about 02:05 UTC at about 4,200 units/h) |
| tier-3 structures / pockets | 23.0 M / 68.1 M | 92.5 M / 292 M | about 226 M / 0.71 B (about 580 GB of pocket files) |
Tiers 1–2 are unchanged.

## Runtime estimate (measured timings, scaled)
| stage | estimate |
|---|---|
| snapshot | 3 min |
| sites tier 3 (0.71 B pockets) | 15–35 min (I/O bound) |
| web index tier 3 | about 16 min |
| SOGs O1 (MSA cache seeded) + τ 0.2/0.3 + benchmark | about 25 min |
| SOGs O2 (MSA cache seeded) + benchmark | about 45 min |
| O2 rescue pool | 74 min (minutes if `--work` reuses the earlier BURROW outputs) |
| Part A, GPU 3 | about 2 h |
| Part B O1 | about 40 min (BLOCKED for a new snapshot) |
Run serially this is about 5.8 h. With an O1/Part A lane and an O2 lane in parallel it is about 3.5 h, which ends
around 06:30, later than 06:00.

## Proposal
Tiers 1–2, and so O1, the Part A universe and Part B O1, have identical inputs on the new snapshot. Only tier-3 sites,
the web index and the tier-3 members of O2 see new data.
- Run at 03:00: snapshot, sites (tier 3 + L), web index, SOGs O2 (seeded), benchmark O2, D20 split, and O2 rescue with
  BURROW reuse. That is about 2–2.5 h and should be done by about 05:30.
- O1, Part A and Part B O1 are reproducibility checks. Run them in parallel only if CPU allows: Part A on GPU 3,
  O1 with the MSA cache seeded.

## Blockers (not patched)
- Part B: `run_o1_pipeline.sh` hard-codes S, and 8 Part B scripts set `SNAP=20261030T085136Z` with no `--snapshot`
  option. There is no Part B driver for O2.
- `build_cases.py` uses the SNAP constant in its input defaults.
- `webapp_build_og_index.py` has no input root.
