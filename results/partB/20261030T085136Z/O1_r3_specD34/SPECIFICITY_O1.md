# Specificity change on O1 (O1_r3, run O1_r3_specD34): D34 analysis

Generated 2026-10-30 23:40 CET by `scripts/partB/write_spec_d34.py` (git 4cb08f1). Producer: Part B specificity agent; **not yet verified**. Runs: `run_specificity.py --orthology O1 --release O1_r3 --run-id O1_r3_specD34 --tau 0.2 --workers 16 --n-null-seeds 3 --rows-dir /tmp/pe_big/partB/specD34/O1_r3/rows` (runtime 180 s, CPU) and `spec_cascade_d34.py` (35 s). Prereg sha256 `52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83`; frozen thresholds unchanged (noise floor 0.3942, >= 2 BLOSUM62<0 clade-consistent substitutions, BH alpha 0.05). Everything below except the prereg BH count is **post hoc (D34)** and is kept out of the convincing headline.

## 1. Prereg BH (headline): 0 calls, and why BH cannot fire

| input | seed | residue gate | gate + shift | BH tested | BH called | min empirical p | min attainable p | BH threshold for one call | tail-normal called | tail-GPD called |
|---|---|---|---|---|---|---|---|---|---|---|
| real |  | 10,584 | 513 | 10,584 | 0 | 0.00013 | 3.5e-05 | 4.7e-06 | 3 | 1 |
| null_pseudo_clade | 0 | 2,230 | 88 | 2,230 | 0 | 0.00013 | 3.3e-05 | 2.2e-05 | 3 | 1 |
| null_pseudo_clade | 1 | 2,260 | 77 | 2,260 | 0 | 0.00011 | 3.3e-05 | 2.2e-05 | 1 | 0 |
| null_pseudo_clade | 2 | 2,243 | 71 | 2,243 | 0 | 0.00013 | 3.3e-05 | 2.2e-05 | 1 | 0 |

The smallest empirical p any real pair can reach is 1/(1 + n_pool) of its matched comparator pool; the best case over tested real pairs is 3.5e-05, above the BH threshold for a single call, alpha/n_tested = 0.05/10,584 = 4.7e-06. So the prereg BH call is unattainable by construction at this null resolution, whatever the data (observed minimum 0.00013).

## 2. Population-level signal (empirical null FDR, the D24 estimator)

| stage | real | null per seed | null mean | empirical FDR |
|---|---|---|---|---|
| residue gate + shift (funnel) | 513 | [88, 77, 71] | 78.67 | 0.153 |
| funnel + artifact cascade | 200 | [37, 44, 30] | 37.00 | 0.185 |
| funnel + cascade, strict (no quality flags) | 55 | [15, 17, 6] | 12.67 | 0.230 |
| tail-normal BH (D17 rule) | 3 | [3, 1, 1] | 1.67 | 0.556 |
| tail-GPD BH | 1 | [1, 0, 0] | 0.33 | 0.333 |

Real sister pairs exceed the pseudo-clade null about 6.5-7x at the funnel stage, and the excess survives the full artifact cascade at a similar ratio. The cascade is applied identically to null pairs (it evaluates the real proteins of the derived/sister clades; the pseudo-clade null permutes member data, not the tree). The pass rates are similar for real and null pairs, so the cascade does not enrich the signal; it removes mostly present-side fragment / length-ratio cases. Caveat: most of the excess is in the residue gate (real clades share substitutions, pseudo-clades do not, by construction), so the funnel ratio is an upper bound on a fingerprint-shift signal.

## 3. Tail-fit sensitivity (post hoc, D34(c))

Same matched comparator pool as the prereg empirical p (depth x log2 clade-size bins, widened to >= 20). Normal (D17 rule, `directional.decide`): normal-fit p beyond the pool maximum, else max(empirical, normal). GPD: PWM fit to exceedances over the pool's q95; beyond the pool maximum the GPD p, else max(empirical, GPD) above q95. BH over the same tested set; same call rule (gate & shift above floor & q <= 0.05). Applied to the null seeds identically.

- Normal tail: real 3 vs null [3, 1, 1]. The normal fit is a poor model for the bounded, right-skewed shift (1 - cos) distribution: it gives p down to ~1e-39 and fires in null seeds too, so it is anti-conservative and not usable as a call rule.
- GPD tail: real 1 vs null [1, 0, 0].
- No null-seed tail call passes the artifact cascade.

Tail-called real pairs:

| sog_id | level | n_d | n_s | species_d | n_subst | substitutions | shift | pool max | n_pool | p_emp | p_normal | q_normal | p_gpd | q_gpd | pass_artifacts | strict | fails |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HOG0034501:g3 | saccharomyceta | 5 | 1 | 4 | 2 | 3:R>T;5:W>L | 0.841 | 0.804 | 5.76e+03 | 0.000174 | 7.18e-15 | 7.6e-11 | 8.42e-05 | 0.272 | True | False |  |
| HOG0034036:g2 | Gammaproteobacteria | 4 | 1 | 4 | 4 | 3:R>A;5:V>W;9:R>T;11:L>W | 0.826 | 0.804 | 5.76e+03 | 0.000174 | 2.32e-14 | 1.23e-10 | 0.000105 | 0.272 | False | False | a_fragment_present |
| HOG0029488:g0 | Spirurina | 4 | 1 | 4 | 4 | 4:G>C;6:A>I;9:S>L;20:G>E | 0.562 | 0.399 | 66 | 0.0149 | 9.08e-08 | 0.00032 | 9.08e-08 | 0.000961 | False | False | a_fragment_present |

## 4. Exploratory table

`specificity/tau0.2/specificity_exploratory.parquet`: all 513 funnel-stage real pairs with the cascade columns and flags `flag_bh_prereg_called`, `flag_tail_normal_called`, `flag_tail_gpd_called`, `flag_pass_artifacts`, `flag_pass_artifacts_strict`, `flag_ge3_species_d`, `flag_tail_and_cascade`. 200 pass the cascade (55 strict); 93 of those have >= 3 derived species. Cascade fail counts: {"a_fragment_present": 209, "a_orthology_o2": 130, "a_orthology_panther": 86, "a_fragment": 52}. All rows are labelled exploratory, not significant. Null rows and per-input cascades are bulk: /tmp/pe_big/partB/specD34/O1_r3/.

## Caveats

- 3 null seeds (as the earlier specificity runs); empirical FDRs on null means of ~8-80 pairs are rough.
- Candidates files (`candidates.parquet`) were not read for writing; the cascade ran in scratch results roots with an emptied candidate table.
- O2 instability / species-tree lineage / PAE checks are `na` for specificity; O2 W members without an atlas site have no fingerprint.
- Design scope (D31): within-orthogroup sister pairs only.
