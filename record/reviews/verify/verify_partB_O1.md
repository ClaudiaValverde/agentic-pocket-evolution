# Independent verification: Part B on O1 (gain/loss calls, null funnel, rescue, controls)

Verifier: independent agent (I did not produce these results). Written 2026-10-30 15:13 CET (`date`).
Branch `pocket-evolution`, HEAD `f97f705` (target commit `029e0d7` plus the ledger-only commit `f97f705`). I edited no
code, result, threshold or test file. `scripts/partB` and `results/partB` were read only.
Check scripts are in `/tmp/verify_partB/` (scratch, not committed). The parsimony, the rule evaluation, the OrthoXML
parser, the species-tree lineage check and the nulls are my own code (`recompute.py`, `null_check.py`). They do not
import `pocketevo`. Inputs are the stored state tables:
- pre-rescue D11 states: `/tmp/pe_big/partB/20261030T085136Z/O1/states_d11_t0.2/` (+ raw `/tmp/pe_big/sogs/.../O1/states_t0.2`, `msa/`, `sogs/`);
- post-rescue states: `results/partB/events/20261030T085136Z/O1/d11/tau0.2/primary/rescue/sog_states/bucket=all.parquet`;
- gene trees: `data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz`; species tree `tier1_guide_tree.nwk`.

**Wall time of my checks: about 8 min of compute (recompute pre 51 s, post 4 s, D11/r2 5 s, null replication 370 s, inspection and controls ~30 s); about 15 min elapsed (14:58 to 15:13).**

## Verdicts

| # | claim / check | verdict |
|---|---|---|
| 1a | Calls reproduce with independent Sankoff 2:1 + D11 (30 gains, 18 losses, 30 non-called; in fact all SOGs) | **PASS** (exact: 179/63 pre-rescue, 51/18 post-rescue, 0 discrepancies) |
| 1b | Each gain has ≥ 3 derived species and ≥ 2 outgroup lineages | **PASS** in the gene tree. **PASS-with-notes** against the species tree: 5/51 gains have < 2 species-tree outgroup lineages, and 33/51 sit just below a duplication node |
| 1c | Input states are the stated SOG release (O1_r2) | **FAIL (non-blocking)**: the run used O1, not O1_r2. Under O1_r2 states the result is 50 gains + 17 losses |
| 2 | Null funnel: same code, same cuts, like for like; rescue symmetry | **PASS-with-notes** |
| 3 | Rescue applied to every member of the candidate SOGs; how W is used | **PASS-with-notes** |
| 4 | Gain cut = prereg minima; losses `fdr_supported = false` | **PASS** |
| 5 | Hand inspection of 10 random gains | **PASS-with-notes**: the absent side is clean. The present side is often low-pLDDT or length-divergent. 34/51 gains carry at least one artifact flag |
| 6 | 5 false-loss controls are artifacts | **PASS-with-notes**: they are rescue-unstable, not low quality |
| 7 | `tests/test_recon.py` | **PASS** (18/18, uv `.venv`, both at `029e0d7` and in the working tree). The faiss env cannot run it (no numba) |

## 1. Independent recomputation

Command: `.venv/bin/python -I /tmp/verify_partB/recompute.py {pre,post}`. Then `sample_check.py`.

Method:
- Own OrthoXML parser, giving one tree per root HOG.
- Leaf costs: P and W forbid absence, A forbids presence, U is free. Cost matrix gain 2, loss 1.
- Down pass plus outside (up) pass. A node is fixed if only one state reaches the optimum.
- An event is an edge whose two endpoints are fixed in opposite states, so it does not depend on tie-breaking.
- Margin = |total(0) − total(1)| at both endpoints (min).
- Outgroup lineages = sibling subtrees of every ancestor of the clade node.
- Rules: prereg minima plus the chosen 2:1 cuts. Gain: ≥ 3 P species, (P+W)/inf ≥ 0.8, ≥ 2 outgroup A members in ≥ 2 lineages, margin ≥ 1. Loss: ≥ 2 A species, A/inf ≥ 0.8, ≥ 3 outgroup P lineages, margin ≥ 2.

| input | SOGs | theirs gain / loss | mine gain / loss | in both | only theirs | only mine |
|---|---|---|---|---|---|---|
| pre-rescue, D11, D8, ≥2P & ≥2A prefilter | 9,323 (= their count) | 179 / 63 | 179 / 63 | 242 | 0 | 0 |
| post-rescue (560 rescued SOGs) | 560 | 51 / 18 | 51 / 18 | 69 | 0 | 0 |

Seeded samples (seed 20261030):
- 30 of the 51 primary gains and all 18 losses agree on every field: n_species_P/A, n_P/W/A, out_A members and lineages, out_P lineages, margin, clade, out-A and out-P species. There are 0 mismatches in 12 fields.
- None of 30 random non-called SOGs is called by my code.

**D11 is correctly applied** (`d11_and_r2.py`). For all 7,294 A rows of the 560 rescued SOGs, I recounted the member's non-gap residues at the SOG footprint columns from the FAMSA MSA. My count equals `n_fp_res` in 7,294/7,294 rows, and state = A iff the count is > 10 in 7,294/7,294. No non-A row was changed.

**Species-tree check of "≥ 3 derived species, ≥ 2 outgroup lineages"** (guide tree = NCBI taxonomy used by the pipeline; `data/timetree` has only pairwise ages):
- All 51 gains have ≥ 3 P species, and in the gene tree they all have ≥ 3 outgroup A lineages (the minimum is 3).
- Mapped onto the species tree, the A members outside the gene clade fall in ≥ 2 distinct lineages outside the clade's species LCA for **46/51**. The other 5 are listed below.

| SOG | species-tree LCA of clade | species-tree outgroup lineages | outgroup A species inside the clade LCA |
|---|---|---|---|
| HOG0027482:g15 | Bilateria | 0 | all 6 (paralogs in fly and nematodes) |
| HOG0033466:g4 | Enterobacteriaceae | 1 | 2 |
| HOG0034186:g10 | Mycobacteriales | 1 | 2 |
| HOG0034899:g8 | Pseudomonadati | 1 | 7 |
| HOG0035533:g7 | Enterobacteriaceae | 1 | 4 |

- 37/51 gains have at least some outgroup A members in species inside the clade's species LCA. These are paralogs.
- 33/51 gain edges have a **duplication parent**, so the "gain" is a pocket present in one paralog copy and absent in its sister copy.
- Pre-rescue, 28/179 gains have < 2 species-tree lineages.

This is a legitimate gene-tree reading of the prereg ("A in ≥ 2 members from ≥ 2 sister/outgroup lineages"), and DESIGN §5.2 includes paralog branches. However, the candidate table labels each gain with `clade_lca_species_tree` (for example "Enterobacteriaceae"), which a reader will take as a gain on that species-tree stem. **Recommendation:** add `parent_kind` and a species-tree outgroup-lineage count to CANDIDATES.md, and report paralog-branch gains separately from species-stem gains.

**Input release (1c).** The task named `results/sogs/.../O1_r2` as the SOG input. The sidecars show the run used O1, not O1_r2:
- `run_events.py` ran with `--sogs results/sogs/20261030T085136Z/O1/web --sog-meta results/sogs/.../O1/sogs_all_tau.parquet`.
- The states came from `/tmp/pe_big/sogs/.../O1/states_t0.2`, through `apply_d11.py`, where BULK is hard-coded to O1.

O1_r2 differs only by its conservative pLDDT lower bound. D8 is unaffected, because every `empty_footprint` SOG also has footprint < 4 columns, so the 6,595 dropped SOGs are the same.

In the 560 rescued SOGs:
- 86 A rows that pass in O1 fail the O1_r2 rule (A→U).
- 83 of them survive D11, and 22 of those are in primary SOGs.
- I set those rows to U in the post-rescue states and recomputed (`post_r2.py`). Primary calls become **50 gains + 17 losses**. HOG0026417:g18 (gain) and HOG0028295:g7 (loss) drop out, and no new call appears.

The effect is small, but the provenance should say "O1 (pre-r2) states". RESEARCH_TREE already lists a rerun on O1_r3 as pending.

Also note: `run_events.py` and `run_rescue.py` have uncommitted edits in the working tree (the D8 flag read and `--bulk`), and the sidecars record git `4fe605f+dirty` / `3f11d4c+dirty`. The exact code that produced the events is therefore not pinned by a clean commit. My reimplementation reproduces the outputs, so this is a provenance note, not a correctness problem.

## 2. Null funnel

Code path. Verified by reading `scripts/partB/run_nulls.py` at `029e0d7`, plus the sidecar argv:
- Real and null inputs go through the same `run()`, which calls `events.run_assembled` with the same `settings` (primary 2:1) and the same `thresholds.json`.
- Both use the same `--prefilter`, D8 meta, D11 state shards, τ 0.2 and the same 9,323 SOGs.
- The real row in `funnel_real_vs_null.csv` (179 / 63) equals `run_events.py` and my recomputation.
- Curveball: masked checkerboard swaps on the SOG×species matrix (strongest member state), with U and no-member cells fixed, remapped to members.
- Pseudo-clade: per SOG, member states are permuted over all leaves of the OG tree.

**Rescue is not applied to either side** of the 179-vs-null comparison. Both are pre-rescue D11 states, so the comparison is like for like.

Caveats:
- (i) The null ratio describes the **pre-rescue** call set. RESEARCH_TREE line 6 attaches "~6× over nulls" to the 51 post-rescue gains. That ratio was never computed post-rescue, and rescue removes 73 % of real pre-rescue gains (130/179). Its effect on null calls is unknown, because pseudo-clade nulls move states onto proteins that have no rescue result.
- (ii) Read as an empirical FDR, the curveball null gives 29.3/179 = **16 %** (max 32/179 = 18 %). That is above the 10 % target and above the simulated FDR of 5 %. The pseudo-clade null gives 8.8 %. `thresholds.json` names a `real_data_rule` (empirical FDR re-check on real SOGs), but I found no output that applies it.
- (iii) The replicate null was not run (no `replicate_states/`), as stated.
- (iv) D22 (swap count 100 × min(#1, #0)) was logged at 14:58, two minutes after the null funnel was written (13:56Z = 14:56 CET). The 50×/200× mixing check it asks for was not in the outputs. My own check is below.
- (v) Losses: curveball 65.7 vs real 63. The real data show no loss signal over the null, so D21 is right.

Independent null replication (`null_check.py`, own code, own seeds):

| null (my code) | swaps (multiplier × min(#1,#0) = 40,046) | swaps done | cells changed | gains called | losses called |
|---|---|---|---|---|---|
| pseudo-clade s11 / s12 / s13 | – | – | – | 13 / 13 / 16 | 20 / 22 / 27 |
| checkerboard 50× s11 | 2.0 M target | 128,246 (attempt cap) | 27,428 | 36 | 73 |
| checkerboard 100× s11 / s12 / s13 | 4.0 M target | ~258,000 (attempt cap) | 27,862 / 27,812 / 27,724 | 36 / 25 / **57** | 69 / 69 / 77 |
| checkerboard 200× s11 | 8.0 M target | 518,899 | 27,742 | 26 | 59 |
| *theirs* curveball s0–s2 / pseudo-clade s0–s2 | 4.0 M | 4,004,600 | ~27,850 | 32 / 27 / 29 ; 14 / 21 / 12 | 65 / 70 / 62 ; 23 / 25 / 23 |
| real (mine = theirs) | | | | 179 | 63 |

- The pseudo-clade null reproduces: 13–16 here vs 12–21 in theirs.
- The checkerboard null agrees in level, but between-seed spread is large: 25–57 gains in my 5 runs vs 27–32 in their 3.
- On mixing: the number of changed cells plateaus at ~27.4–27.9 k from 50× to 200× in my chain and at 100× in theirs. So mixing is adequate, even though my sampler accepts fewer swaps.
- With 3 seeds, the null mean is poorly estimated. Pooling my 5 checkerboard runs gives a mean of 36 gains, an empirical FDR of about 20 % (max 57/179 = 32 %).
- The real gain excess over the null is robust: 179 vs ≤ 57.
- Neither null supports a loss signal. Real losses (63) are at null level (59–77).
- **Recommendation:** run ≥ 20 null seeds and report the null as an empirical FDR next to the simulated 5 %.

## 3. Rescue

`rescue_summary.json`, `rescue_hits.parquet`, `rescue/sog_states`, and `run_rescue.py` / `states.apply_rescue` at `029e0d7`:
- The rescued SOGs are exactly the SOGs with a pre-rescue call at any cost setting: 560 = 560.
- **Every P and A member** of those SOGs (15,196/15,196 rows) has a rescue row, and every row has all 6 settings run.
- There are no multi-fragment or no-model members and 0 failures.
- The rule is symmetric and complete: A→W iff ≥ 1/6 settings hit (Jaccard ≥ 0.3); P→W iff 0/6 hit. Crosstab: A 7,035 → A 3,660 + W 3,375; P 8,161 → P 8,030 + W 131.
- There is no selective rescue inside a SOG.

How W is used. I verified this by reimplementing it, and the reimplementation reproduces the calls exactly:
- W is "present" in the Sankoff leaf costs (absence is forbidden).
- W counts toward (P+W)/informative for gains, and against the loss fraction.
- W never counts toward n_species_P, and has state confidence 0.5.
- An outgroup A rescued to W therefore removes gain support.

Notes:
- (a) Rescue is applied only to SOGs that were called before rescue. It is not purely subtractive: **7 post-rescue calls (2 gains, 5 of the 18 losses) were not called pre-rescue at 2:1.** They come from SOGs that were rescued because of another call, where A→W reshaped the MPR. Equivalent SOGs that had no pre-rescue call were never rescued, so such rescue-created calls are missed there. This is a small selection asymmetry. It is worth stating, because 5/18 losses exist only because of rescue.
- (b) k = 1 rescues 35 % of true-negative control regions (D9), so W is a weak signal.

## 4. Thresholds and labels

- `config/prereg.yaml` sha256 = `52f94618…9d83`, unchanged since `d586f1a`, and equal to the hash in `thresholds.json` and every sidecar.
- The chosen 2:1 gain cut is min_species 3 and min_out_lineages 2, the same as the prereg minima. The frozen frac 0.8 and the ≥ 2 outgroup A members are also applied.
- The extra `min_margin = 1` does not bind: with integer costs, an unambiguous event has margin ≥ 1. Pre-rescue, candidate = called = 179.
- `Rules.check_not_relaxed` guards against relaxing the prereg minima.
- Labels: all 108 loss rows have `fdr_supported = False` / "not FDR-calibrated (audit)", and all 845 gain rows have True.
- The applied loss cut is min_species 2, out-P lineages ≥ 3, margin ≥ 2. My recomputation needs exactly these values to reproduce the 63 / 18.
- Note: ledger D7 still states "loss ≥ 3 species & ≥ 3 lineages". That was superseded by the measured-noise re-derivation (species 2), but D7 was not updated.
- The gain "FDR-calibrated" label rests on recall ≈ 0.1 % in simulation. This is stated in CANDIDATES.md.

## 5. Hand inspection of 10 random gains

Seed 20261030; `inspect_gains.py`, tables `/tmp/verify_partB/inspect_gains{,_summary}.csv`. For all 51 gains, see `gain_quality.py`. The fragment check is length vs the OG median, recomputed from the MSA sequences. It agrees with `len_ratio`.

**Absent side (outgroup A): clean in all 51 gains.**
- Footprint pLDDT is ≥ 90.0 (min 90.02).
- 0 members have a length outside 0.7–1.3 × the OG median.
- The aligned footprint fraction is ≥ 0.82, and there are ≥ 11 footprint residues.

These are guaranteed by the A rule and confirmed independently. One caveat: many A members sit right at the D11 floor of 11 residues, and footprints are only 11–22 columns.

**Present side (derived P): no quality rule applies, and it is where the artifacts are.**

| SOG | fp size | P pLDDT (min / med) | P length / OG median | flag |
|---|---|---|---|---|
| HOG0026417:g18 | 20 | 76.8 / 77.9, frac<70 up to 0.25 | 0.93–0.95 | pLDDT asymmetry (P ~78 vs A ≥ 90); also drops under O1_r2 |
| HOG0029217:g0 | 22 | 82.8 / 84.8 | **1.27–1.56** | nematode P copies longer (insertion) + low pLDDT |
| HOG0029246:g3 | 13 | 93.5 / 96.3 | 0.95–0.99 | clean |
| HOG0030712:g13 | 11 | 95.4 / 97.5 | 1.10–1.30 | clean-ish; outgroup A mostly in-species paralogs; fp 11 columns |
| HOG0031031:g3 | 11 | 91.9 / 94.5 | **0.58–0.67** | P proteins ~180 aa vs 309 median (short paralog in 4 vertebrates): possible truncation-exposed cavity |
| HOG0034782:g3 | 11 | 88.2 / 92.8 | 0.87–1.00 | plausible; deep (parent Eukaryota); two gains in one SOG (f1 false) |
| HOG0035256:g19 | 12 | 87.0 / 90.2 | **0.16–0.30** | **likely artifact**: P members are 70–130 aa fragments of a 428-aa OG |
| HOG0035363:g0 | 14 | 69.2 / 80.0, frac<70 up to 0.5 | 0.93–1.18 | SLC1/EAAT transporters (P43003, P43005, Q15758…): P pLDDT 70–83 vs A ≥ 96; Jaccard down to 0.20 |
| HOG0036328:g8 | 19 | 90.0 / 90.8 | **1.41–1.59** | P copies 40–60 % longer (gene model / extra domain) |
| HOG0036512:g4 | 12 | 96.8 / 97.1 | 1.06–1.35 | clean |

Across all 51 gains:
- 17 have derived-P median footprint pLDDT < 90 (4 below 80).
- 18 have most P members outside 0.7–1.3 × the OG median length.
- 5 fail the species-tree lineage check.
- **34/51 carry at least one flag, and 17 are clean.**

This is a systematic asymmetry: the A rule screens length and pLDDT on the absent side only. Gains therefore enrich for derived clades whose models are lower confidence or length-divergent. CANDIDATES.md says the artifact cascade (fragment/gene-model, PAE/TED, …) is not yet applied. These flags show it must be applied before any gain is used.

## 6. False-loss negative controls

`controls.py`; 5 designated controls in `false_loss_controls.parquet`:
- They are not low-quality artifacts. All clade A members have footprint pLDDT 93.8–98.3, length 0.85–1.12 × median, aligned fraction 1.0, and 11–12 footprint residues.
- They are rescue artifacts: every clade A member was re-detected (A→W) by ≥ 1 of the 6 settings. Most of the hits come from `min_score_none`, a pocket that exists but falls below the score floor.

| control | clade A | members hit by ≥ 2 settings | ≥ 3 settings |
|---|---|---|---|
| HOG0035844:g7 (Eukaryota) | 16 | 88 % | 38 % |
| HOG0035864:g3 | 11 | 82 % | 45 % |
| HOG0035379:g7 | 9 | 100 % | 78 % |
| HOG0036716:g3 | 6 | 67 % | 33 % |
| HOG0028153:g6 | 5 | 40 % | 20 % |

They are fair examples of "pocket-finder instability turns a loss into W". However, "truly artifacts" is not shown:
- k = 1 also rescues 35 % of true-absence control regions.
- At k = 2, HOG0028153:g6 would keep 3/5 A members.

Describe them as "rescue-unstable losses", not as "proven false". The five were chosen by the largest rescued fraction, which is fine for illustration.

## 7. Tests

- `cd /tmp/verify_partB/at029 && PYTHONPATH=src <repo>/.venv/bin/python -m pytest tests/test_recon.py -q` gives **18 passed** (6.4 s). This is the `git archive 029e0d7` tree.
- The same file in the working tree also gives 18 passed. (`src/pocketevo/events.py` was modified in the working tree during my run. The change is a specificity p-value speed-up and does not touch gain/loss.)
- `/tmp/pe_big/envs/faiss/bin/python` fails at collection: `ModuleNotFoundError: numba`.

## Summary for the PI

- The gain/loss calls are computed as specified. An independent parsimony reproduces 179/63 and 51/18 exactly.
- Thresholds and labels match the prereg and D21.
- Rescue is complete and symmetric.
- The null comparison is like for like, but only before rescue.

Before relying on the 51 gains:
1. Rerun on O1_r2/r3 states (50 gains + 17 losses on r2).
2. Separate the 33 paralog-branch gains from the species-stem gains, and add a species-tree lineage count (5 gains fail it).
3. Apply the fragment/gene-model and pLDDT-asymmetry parts of the cascade (34/51 flagged; HOG0035256:g19 is a likely fragment artifact).
4. Report the curveball null as an empirical FDR (16 % from their seeds, ~20 % from mine, pre-rescue, using ≥ 20 seeds).
5. Do not quote "6× over null" for the post-rescue 51.
