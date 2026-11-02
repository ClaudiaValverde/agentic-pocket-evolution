# Independent verification: D34 specificity population signal and co-folding follow-up

- Verifier: independent verifier agent (Claude Opus 5.5). I did not produce any of the results checked here. Not a PI review.
- Date: 2026-10-30T22:48:15Z (`date -u`); repo HEAD at review time 4ff737e, branch pocket-evolution.
- Commits checked: 75d1df0, 8eba86a, 8e8d456 (specificity); ab4de9d, 464a0b9, 001be32 (co-folding).
- Scripts (all mine, CPU only, read-only on project files): `/tmp/verify_d34/verify_spec.py`, `/tmp/verify_d34/spot.py`,
  `/tmp/verify_d34/na.py`, `/tmp/verify_d34/verify_cofold.py`. Outputs: `/tmp/verify_d34/spec_out.json`,
  `/tmp/verify_d34/cofold_allsamples.tsv` (270 sample models), `/tmp/verify_d34/cofold_pockets.tsv`.
- Python: `<repo>/.venv/bin/python`. Co-fold geometry was recomputed with my own
  gemmi/numpy code, not by importing `scripts/cofold`. BURROW was re-run with `pocketcompare describe --af-mode --device cpu`
  on apo structures I wrote myself.

## Verdict summary

| # | Claim | Verdict |
|---|---|---|
| A1 | Funnel and cascade counts, real vs null, and empirical FDR (O1 513 vs 88/77/71, FDR 0.153; O2 536 vs 84/66/80, FDR 0.143; cascade 200 vs 37/44/30, FDR 0.185; 230 vs 45/32/43, FDR 0.174; strict 55 / 50; tail counts) | PASS |
| A2 | The prereg BH cannot fire (minimum attainable p is above alpha/n_tested) | PASS |
| A3 | The null goes through the identical code path and cascade | PASS-with-notes |
| A4 | The two O2 runs (O2_r1_spec2 and O2_r1_specD34) agree | PASS (bit-identical; this is reproducibility, not independent replication) |
| A5 | "Real exceeds null about 6.5-7x", and "the excess survives the full artifact cascade at a similar ratio" | PASS-with-notes (6.5x/7.0x is the funnel point estimate; after the cascade it is 5.4x/5.75x; the 3-seed CIs are wide) |
| A6 | Pseudoreplication (whether passing pairs concentrate in a few families) | PASS (no concentration) |
| A7 | Spot-check of 5 real pairs | PASS-with-notes (no outright artifact, but 4 cascade checks are never evaluated) |
| B1 | The 3 gate-passing cases: site recall, centroid and ligand pocket recomputed from the models | PASS (matches the stored tables exactly) |
| B2 | "Ligand 11-34 A off the footprint" | PASS-with-notes (all 270 sample models: 12.0-49.0 A; per-member medians 12.2-34.0 A; recall at most 0.08) |
| B3 | "Footprint pocket detected 3/3 present vs 0/3 absent (MetE), 2/3 vs 0/3 (choline kinase)" | PASS-with-notes (reproduced exactly; n is tiny and only the top model was used) |
| B4 | The decision code now requires the D23 composite | PASS-with-notes (some fields are stale) |
| B5 | The blinding status is recorded as UNBLINDED | PASS-with-notes (recorded in all 4 case.yaml files; not in decision.json or docs/COFOLD.md) |

---

## A. Specificity, D34

Inputs: the bulk rows `/tmp/pe_big/partB/specD34/{O1_r3,O2_r1}/rows/rows_{real,null0,null1,null2}.parquet` and the cascades
`/tmp/pe_big/partB/specD34/*/cascade_*.parquet`. I checked these against the committed `funnel_real_vs_null.csv`,
`specificity_exploratory.parquet`, `SPECIFICITY_*.md` and `results/numbers.json` (`partB.specificity_D34_O1/O2`, `partB.specificity_O2`).

Command: `.venv/bin/python /tmp/verify_d34/verify_spec.py`

### A1. Counts and FDR, recomputed

| release | stage | real | null per seed | null mean | FDR (mean/real) | real/null |
|---|---|---|---|---|---|---|
| O1_r3 | residue gate | 10,584 | 2230/2260/2243 | 2244 | 0.212 | 4.72 |
| O1_r3 | gate + shift (funnel) | 513 | 88/77/71 | 78.7 | 0.153 | 6.52 |
| O1_r3 | + cascade | 200 | 37/44/30 | 37.0 | 0.185 | 5.41 |
| O1_r3 | + cascade, strict | 55 | 15/17/6 | 12.7 | 0.230 | 4.34 |
| O1_r3 | tail-normal BH | 3 | 3/1/1 | 1.67 | 0.556 | 1.8 |
| O1_r3 | tail-GPD BH | 1 | 1/0/0 | 0.33 | 0.333 | 3 |
| O2_r1 | residue gate | 19,462 | 5596/5295/5493 | 5461 | 0.281 | 3.56 |
| O2_r1 | funnel | 536 | 84/66/80 | 76.7 | 0.143 | 6.99 |
| O2_r1 | + cascade | 230 | 45/32/43 | 40.0 | 0.174 | 5.75 |
| O2_r1 | + cascade, strict | 50 | 9/9/7 | 8.33 | 0.167 | 6.0 |
| O2_r1 | tail-normal BH | 6 | 1/0/1 | 0.67 | 0.111 | 9 |
| O2_r1 | tail-GPD BH | 1 | 0/0/0 | 0 | 0 | n/a |

Every number matches the reports, the CSVs, `cascade_summary.json` and numbers.json. The committed exploratory parquets have
513 and 536 rows, with 200/55 and 230/50 passing (cascade and strict); 93 and 129 of the passing pairs have at least 3 species.
Prereg BH calls are 0 in all 8 inputs.

There is one small numbers.json issue. `specificity_D34_O1.funnel[*].min_p_tail_normal` is stored as `0.0`, and
`tail_called_real.p_tail_normal` as `0.0` / `q 1e-10`. The CSV holds 7.2e-15 and 7.6e-11. These values were rounded to 10
decimals, so the stored zeros are an artifact of the rounding. Cosmetic, but a p of 0.0 should not appear in numbers.json.

### A2. Prereg BH is unattainable

Recomputed min over tested pairs of 1/(1+n_pool): 3.46e-5 for O1 and 1.60e-5 for O2. alpha/n_tested is 4.72e-6 and 2.62e-6.
So no single pair can reach q <= 0.05. **PASS.** This is a property of the null resolution, as the report says. The headline
"0 calls" is therefore uninformative by construction; it is not a negative result.

### A3. The null uses the identical code path

- In `scripts/partB/run_specificity.py`, real and null inputs go through the same loop
  (`for seed in [None] + list(range(n_null_seeds))`). That loop calls the same `_og` → `events.specificity_tests`, then
  `events.call_specificity`, then `d34_columns`. The only difference is in `_og`: when `seed` is set, it permutes the
  residue, fingerprint and covariate dicts among the P/W leaves of each SOG tree.
- The (sog_id, node, d_node) pair universe is identical between real and each null seed (`pair_universe_identical: true`),
  with 222,489 pairs for O1 and 743,349 for O2.
- `scripts/partB/spec_cascade_d34.py` handles every `rows_*.parquet` with the same `run_artifacts.py` subprocess. For every
  input, the set of cascade rows equals the set of funnel rows (`cascade_rows_eq_funnel: true`).
- **Note 1.** The cascade evaluates the proteins at the tree positions. Under the null, the fingerprint/residue data at each
  position comes from a different, permuted member. So for null pairs the fragment and orthology checks are not looking at
  the proteins that generated the signal. The code path is identical, but the cascade is not semantically equivalent for
  the null. The report says this in parentheses. It is one reason the null pass rate is higher (see A5).
- **Note 2.** "Full artifact cascade" overstates what ran. In every real and null pair, `a_footprint_plddt`, `a_disorder`
  and `a_pae_ted` are `na`, and `a_instability` is `na` in 510/513 (O1) and 536/536 (O2). `a_orthology_panther` is `na`
  in 18 (O1) and 162 (O2). In practice the effective cascade is fragment, fragment_present, orthology_o2 (O1 only) and
  panther. The report only mentions instability, species-tree and PAE as `na`.
- **Note 3.** The pseudo-clade null destroys clade-consistent substitutions by construction. The residue gate alone gives
  real/null of 4.7x (O1) and 3.6x (O2). Conditional on passing the residue gate, the shift pass rate is:
  - O1: real 4.85% vs null 3.51%, ratio 1.38 (Fisher OR 1.40, p = 2e-5);
  - O2: real 2.75% vs null 1.40%, ratio 1.96 (OR 1.99, p = 4e-19).

  So the part of the 6.5-7x that can be attributed to the fingerprint shift is about 1.4-2x. The report flags the funnel
  ratio as an "upper bound". That is correct, and this is the number that bound refers to.

### A4. O2_r1_spec2 vs O2_r1_specD34

I compared `O2_r1_spec2/.../specificity_tested.parquet` with the D34 real rows. All 19,102 tested pairs match on
(sog_id, node, d_node), none are left over on either side, max |Δshift| = 0, max |Δp| = 0, n_subst is equal everywhere, and
both runs have 536 funnel pairs. The funnel CSVs are identical for all 4 inputs. The same holds for O1: the original
`O1_r3/specificity` against the O1_r3_specD34 rerun, 10,584/10,584, Δ = 0.

**PASS.** These are the same deterministic code on the same inputs, so this confirms reproducibility, not independent
replication. The spec2 "top-5 exploratory cascade 2/5 pass, 0 strict" is a subset view and does not conflict with D34.

### A5. How robust is "6.5-7x" to 3 seeds

Real/null ratio with 95% intervals. The t interval uses df = 2 on the mean of the 3 null seeds. The Poisson interval is exact,
on the pooled null count divided by 3. Null variance/mean is 0.9-1.3 at the funnel and cascade stages, so Poisson is a
reasonable model.

| release | stage | point | seed-wise range | t95 (3 seeds) | Poisson 95% |
|---|---|---|---|---|---|
| O1 | funnel | 6.52 | 5.8-7.2 | 5.1-9.0 | 5.7-7.4 |
| O1 | cascade | 5.41 | 4.5-6.7 | 3.7-10.2 | 4.5-6.6 |
| O1 | strict | 4.34 | 3.2-9.2 | 2.0-∞ | 3.2-6.1 |
| O2 | funnel | 6.99 | 6.4-8.1 | 5.4-10.1 | 6.1-8.0 |
| O2 | cascade | 5.75 | 5.1-7.2 | 4.0-10.2 | 4.8-6.9 |
| O2 | strict | 6.0 | 5.6-7.1 | 4.5-9.1 | 4.1-9.3 |

Funnel FDR t95 is 0.11-0.20 (O1) and 0.10-0.19 (O2). Cascade FDR t95 is 0.10-0.27 (O1) and 0.10-0.25 (O2).

- That real exceeds null is robust: the lower bound is at least 3.7x at the cascade stage under either interval.
- The specific "6.5-7x" is a point estimate at the funnel. After the cascade the ratio is 5.4x/5.75x, not "a similar ratio".
- The null passes the cascade more often than real pairs do:
  - O1: 47% null vs 39% real, Fisher p = 0.046;
  - O2: 52% null vs 43% real, p = 0.022.

  "Pass rates similar" therefore slightly understates a modest depletion of real pairs by the cascade.
- D24 asks for at least 20 seeds and only 3 were run. Tail-stage FDRs rest on 0-3 null calls and carry no information
  (O1 tail-normal FDR t95 runs from 0 to 1.5).

### A6. Pseudoreplication

Cascade-passing real pairs collapsed to root HOG (`sog_id` prefix):

- O1: 200 pairs cover 189 SOGs and **185 distinct root HOGs**. The largest family has 2 pairs, 170 roots have a single pair,
  the top 10 roots hold 10% of pairs, and the effective number of roots is 174.
- O2: 230 pairs cover 213 SOGs and **210 distinct root HOGs**. The largest family has 3 pairs (HOG0002719), the top 10 hold
  9%, and the effective number of roots is 194.
- At root level, real vs null is still about 5-6x: O1 funnel 451 vs 85/72/70, cascade 185 vs 35/41/30; O2 funnel 471 vs
  82/64/80, cascade 210 vs 45/32/43.

**PASS:** the excess does not come from a few families.

### A7. Spot-check of 5 real pairs

Five O1 cascade-passing pairs were drawn at random (`random_state=7`, `/tmp/verify_d34/spot.py`):

| sog_id | level, kind | n_d / n_s / species_d | substitutions | shift | p_emp | issues |
|---|---|---|---|---|---|---|
| HOG0028540:g1 | Eukaryota, speciation | 16 / 1 / 4 | 0:N>V;8:C>M;10:P>A;13:P>S;14:G>P | 0.71 | 0.0034 | sister is a single sequence; present-side min footprint pLDDT 84 |
| HOG0029830:g9 | Trypanosomatidae, duplication | 4 / 1 / 3 | 4:F>H;5:M>D | 0.42 | 0.022 | paralog branch (flag_paralog_branch); single-sequence sister; pLDDT 83 |
| HOG0019324:g0 | Rhabditida, speciation | 2 / 4 / 1 | 5:R>S;6:T>P | 0.48 | 0.011 | derived clade is one species; panther `na`; pLDDT 80 |
| HOG0016220:g0 | leotiomyceta, speciation | 2 / 2 / 2 | 19:E>G;21:I>R;22:G>L | 0.60 | 0.0046 | panther `na`; pLDDT 79 |
| HOG0035703:g6 | Eurotiomycetes, speciation | 2 / 2 / 2 | 0:K>Y;1:K>H;4:D>H;6:P>A | 0.41 | 0.011 | none (strict pass); pLDDT 94 |

None of the five has an outright artifact: no X or gap substitutions anywhere among the passing pairs, and O2 consistency
is 0.94-1.0. But they are weak, as a group:

- individual p values are 0.003-0.02, which is null-like;
- 4/5 have present-side footprint pLDDT < 90;
- 2/5 are paralog or single-species comparisons;
- 4 cascade checks were never evaluated.

Across all passing pairs:

- present-side footprint pLDDT < 90: 114/200 (O1) and 146/230 (O2);
- paralog-branch: 72/200 and 93/230;
- single-sequence sister clade (n_s = 1): 68/200 and 101/230;
- single-species derived clade: 25/200 and 29/230.

**PASS-with-notes.** This is a population-level excess. No individual pair is supported.

---

## B. Co-folding Tier 1 follow-up

Inputs: AF3 outputs `/tmp/pe_big/cofold/cases/<case>/af3/*/<job>/` (top model plus 3 seeds × 5 samples), `blind_key.json` and
`jobs.json` (footprint residues per job). For each of the 18 cognate jobs I computed:

- **ligand centroid distance**: |ligand heavy-atom centroid − footprint CA centroid|;
- **footprint site recall**: share of footprint residues with any heavy atom within 5 Å of the ligand;
- **ligand site**: protein residues within 4.5 Å of the ligand.

I then stripped the ligand, ran BURROW (`pocketcompare describe --af-mode`, CPU) and found the pocket that best covers the
ligand site and the pocket that best covers the footprint. "Detected" means recall >= 0.3.

Command: `CUDA_VISIBLE_DEVICES="" .venv/bin/python /tmp/verify_d34/verify_cofold.py`

### B1. Recomputed values match the stored tables

Top-model centroids and recalls match `case_table.tsv` to 0.01 Å. Ligand pocket ids and ranks, footprint pocket ids, ranks and
recalls, and the ligand-site/footprint overlap (0 residues in all 18) match `ligand_vs_burrow_pockets.tsv` exactly.

| case | ligand pocket rank (my run) | footprint detected, present | footprint detected, absent |
|---|---|---|---|
| MetE (gain_HOG0035156) | MET: rank 2 ×5, rank 3 ×1 (6/6) | 3/3 (ranks 3, 5, 6; recall 1.0/0.92/0.85) | 0/3 (recall 0, 0.077, 0) |
| choline kinase (dupgain_HOG0054071) | ADP: rank 1 ×4, 2 ×1, 5 ×1 | 2/3 (ranks 3, 8; recall 0.75, 1.0; third 0.083) | 0/3 (0.083, 0.083, 0) |
| citrate synthase (gain_HOG0035403) | OAA: rank 1 ×6 | 0/3 | 1/3 (K7MR36, rank 34, recall 1.0) |

In no member is the ligand pocket the same as the footprint pocket.

### B2. "Ligand 11-34 Å off the footprint"

Over all 270 cognate sample models, the centroid ranges and maximum recalls are:

| case | centroid range | max recall |
|---|---|---|
| choline kinase | 12.0-20.0 Å | 0 |
| MetE | 12.0-38.7 Å | 0.08 (1 of 13 residues) |
| citrate synthase | 19.0-49.0 Å | 0 |

The 49.0 Å value is a single-sample outlier in Q0P7U8, and O68883 has one at 39.3 Å.

- Per-member medians are 12.2-34.0 Å and top models 12.2-33.9 Å, so "11-34 Å" fits the medians and top models, not every sample.
- No sample comes within 2 Å or reaches recall 0.9, so the qualitative claim holds.
- The 464a0b9 commit message gives MetE as "centroid 21-34 A". MetE cognate per-member medians are actually 32.9-34.1 Å, and
  some samples are as close as 12 Å. The 21-25 Å values are decoy (9X6) jobs. The sub-range is mislabelled.
- "site recall 0.000 in all 12 member jobs" is approximately true: the MetE P9WK07 median is 0.016, and some samples reach 0.08.

**PASS-with-notes.**

### B3. "Detected 3/3 vs 0/3 (MetE), 2/3 vs 0/3 (choline kinase)"

Reproduced exactly. Caveats:

1. **Small n.** Fisher p is 0.10 for MetE alone and 0.015 for the pooled 5/6 vs 0/6.
2. **Top model only.** Only 1 of 15 samples per job was described.
3. **Two absent-side choline-kinase members have a single BURROW pocket.** P46560 and P35790 have n_pockets = 1, so a
   footprint pocket could hardly be found there.
4. **Not independent of the atlas structures.** These are AF3 models of the same sequences, so "BURROW independently
   reproduces the gain" holds for a fresh structure predictor, not for fresh data.
5. **The footprint residues differ across the event.** For MetE, present-side `..TPYGAV` vs absent-side `..ALP.A..`. This is
   consistent with a real local sequence change, but it also means the "footprint" on the absent side is a mapped alignment
   position.

### B4. The decision code requires the D23 composite

In `scripts/cofold/case_cofold.py` (HEAD), `d23_binds()` implements pLDDT >= 92 AND (centroid <= 2.0 OR recall >= 0.9) in >= 2
of 3 seeds, with thresholds read from `threshold.json` (t_plddt 92, t_centroid 2.0, t_recall 0.9).

`decide()` emits:

- `gain_in_derived` only if all(passes_gain) AND derived binds AND NOT ancestral binds;
- `loss_in_derived` symmetrically;
- `contradicts_call` only when ancestral binds and derived does not;
- `uninformative_ligand_not_at_footprint` when no cognate job binds, min median centroid > 2 Å and max recall < 0.1.

All three `decision.json` files say `uninformative_ligand_not_at_footprint`. My per-sample geometry (no sample within 12 Å, recall
at most 0.08) means no aggregation could make any job pass. Notes:

- **(a)** The `decide()` change is in `case_cofold.py`, whose last commit is bf803ee, not in 464a0b9; 464a0b9 changed
  `run_pool.py`. The code is correct at HEAD.
- **(b)** Stale text remains in every Tier 1 `case.yaml` cofold block: `decision_rule: "all seeds must meet D > X for BOTH
  ligand pLDDT and min P-L PAE"`, together with `seeds_loss: [true, true, true]`. The `decision.json` `"rule"` string has the
  same text. Both still describe the old rule and could be misread as a loss signal.
- **(c)** `d23_call.json` writes `recall: null` for every seed. `d23_call.py` reads `site_recall`/`recall` from a different key
  layout. The verdict is unaffected because the centroid term fails, but the per-seed recall record is missing.

**PASS-with-notes.**

### B5. Blinding recorded as UNBLINDED

`blinding.status: UNBLINDED` is present in the cofold block of all 4 Tier 1 `case.yaml` files, with an explanation of the salt
bug, and it is reflected in RESEARCH_TREE.md (line 85, read-only). It is not in `results/cofold/cases/*/decision.json`,
`cofold_spec.yaml` or `docs/COFOLD.md`. **PASS-with-notes:** the label should also go in the decision.json files and in COFOLD.md.

## Bottom line

- **Specificity.**
  - The D34 numbers are reproduced exactly, and the null uses the same code.
  - Real exceeds pseudo-clade null robustly (≥ 3.7x at the cascade stage at 95%), and the excess is not driven by a few families.
  - But "6.5-7x" is a funnel point estimate, and after the cascade it is 5.4-5.75x.
  - The part attributable to the fingerprint shift (given the residue gate) is only 1.4-2.0x.
  - The cascade runs with 4 checks `na`.
  - It is a population-level signal with no individually supported pair.
- **Co-folding.** The geometry and pocket tables are reproduced exactly from the models. Fix the "21-34 Å" MetE sub-range, the
  stale `decision_rule`/`seeds_loss` fields, the null recall in `d23_call.json`, and add the UNBLINDED label to decision.json.
