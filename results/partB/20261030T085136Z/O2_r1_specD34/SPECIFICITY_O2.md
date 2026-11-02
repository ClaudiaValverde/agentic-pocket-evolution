# Specificity change on O2 (O2_r1, run O2_r1_specD34): D34 analysis

Generated 2026-10-30 23:40 CET by `scripts/partB/write_spec_d34.py` (git 4cb08f1). Producer: Part B specificity agent; **not yet verified**. Runs: `run_specificity.py --orthology O2 --release O2_r1 --run-id O2_r1_specD34 --tau 0.2 --workers 16 --n-null-seeds 3 --rows-dir /tmp/pe_big/partB/specD34/O2_r1/rows` (runtime 318 s, CPU) and `spec_cascade_d34.py` (74 s). Prereg sha256 `52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83`; frozen thresholds unchanged (noise floor 0.3942, >= 2 BLOSUM62<0 clade-consistent substitutions, BH alpha 0.05). Everything below except the prereg BH count is **post hoc (D34)** and is kept out of the convincing headline.

## 1. Prereg BH (headline): 0 calls, and why BH cannot fire

| input | seed | residue gate | gate + shift | BH tested | BH called | min empirical p | min attainable p | BH threshold for one call | tail-normal called | tail-GPD called |
|---|---|---|---|---|---|---|---|---|---|---|
| real |  | 19,462 | 536 | 19,102 | 0 | 1.8e-05 | 1.6e-05 | 2.6e-06 | 6 | 1 |
| null_pseudo_clade | 0 | 5,596 | 84 | 5,461 | 0 | 4.9e-05 | 1.6e-05 | 9.2e-06 | 1 | 0 |
| null_pseudo_clade | 1 | 5,295 | 66 | 5,199 | 0 | 5.2e-05 | 1.6e-05 | 9.6e-06 | 0 | 0 |
| null_pseudo_clade | 2 | 5,493 | 80 | 5,394 | 0 | 3.3e-05 | 1.6e-05 | 9.3e-06 | 1 | 0 |

The smallest empirical p any real pair can reach is 1/(1 + n_pool) of its matched comparator pool; the best case over tested real pairs is 1.6e-05, above the BH threshold for a single call, alpha/n_tested = 0.05/19,102 = 2.6e-06. So the prereg BH call is unattainable by construction at this null resolution, whatever the data (observed minimum 1.8e-05).

## 2. Population-level signal (empirical null FDR, the D24 estimator)

| stage | real | null per seed | null mean | empirical FDR |
|---|---|---|---|---|
| residue gate + shift (funnel) | 536 | [84, 66, 80] | 76.67 | 0.143 |
| funnel + artifact cascade | 230 | [45, 32, 43] | 40.00 | 0.174 |
| funnel + cascade, strict (no quality flags) | 50 | [9, 9, 7] | 8.33 | 0.167 |
| tail-normal BH (D17 rule) | 6 | [1, 0, 1] | 0.67 | 0.111 |
| tail-GPD BH | 1 | [0, 0, 0] | 0.00 | 0.000 |

Real sister pairs exceed the pseudo-clade null about 6.5-7x at the funnel stage, and the excess survives the full artifact cascade at a similar ratio. The cascade is applied identically to null pairs (it evaluates the real proteins of the derived/sister clades; the pseudo-clade null permutes member data, not the tree). The pass rates are similar for real and null pairs, so the cascade does not enrich the signal; it removes mostly present-side fragment / length-ratio cases. Caveat: most of the excess is in the residue gate (real clades share substitutions, pseudo-clades do not, by construction), so the funnel ratio is an upper bound on a fingerprint-shift signal.

## 3. Tail-fit sensitivity (post hoc, D34(c))

Same matched comparator pool as the prereg empirical p (depth x log2 clade-size bins, widened to >= 20). Normal (D17 rule, `directional.decide`): normal-fit p beyond the pool maximum, else max(empirical, normal). GPD: PWM fit to exceedances over the pool's q95; beyond the pool maximum the GPD p, else max(empirical, GPD) above q95. BH over the same tested set; same call rule (gate & shift above floor & q <= 0.05). Applied to the null seeds identically.

- Normal tail: real 6 vs null [1, 0, 1]. The normal fit is a poor model for the bounded, right-skewed shift (1 - cos) distribution: it gives p down to ~1e-39 and fires in null seeds too, so it is anti-conservative and not usable as a call rule.
- GPD tail: real 1 vs null [0, 0, 0].
- No null-seed tail call passes the artifact cascade.

Tail-called real pairs:

| sog_id | level | n_d | n_s | species_d | n_subst | substitutions | shift | pool max | n_pool | p_emp | p_normal | q_normal | p_gpd | q_gpd | pass_artifacts | strict | fails |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HOG0043018:g0 | Ecdysozoa | 3 | 3 | 3 | 3 | 13:Q>A;17:G>V;18:Q>Y | 1.17 | 1.15 | 5.45e+04 | 1.83e-05 | 1.39e-39 | 2.65e-35 | 5.54e-06 | 0.053 | False | False | a_fragment_present;a_orthology_panther |
| HOG0044697:g11 | Rhabditida | 3 | 1 | 3 | 5 | 0:F>D;1:P>R;4:C>S;6:F>V;8:P>H | 1.07 | 0.956 | 3.25e+04 | 3.08e-05 | 2.92e-34 | 1.86e-30 | 2.47e-05 | 0.157 | False | False | a_fragment_present;a_orthology_panther |
| HOG0049595:g1 | Sordariomycetes | 2 | 8 | 2 | 2 | 6:G>T;7:C>G | 0.92 | 0.863 | 2.27e+04 | 4.4e-05 | 6.7e-37 | 6.4e-33 | 7.75e-05 | 0.172 | True | False |  |
| HOG0046460:g2 | Gnathostomata | 12 | 1 | 10 | 3 | 2:C>D;3:C>L;4:A>E | 0.879 | 0.851 | 1.06e+04 | 9.42e-05 | 1.26e-25 | 4.81e-22 | 8.06e-05 | 0.172 | True | False |  |
| HOG0043054:g2 | Bilateria | 21 | 1 | 13 | 2 | 2:H>C;3:A>L | 0.863 | 0.73 | 3e+03 | 0.000333 | 1.18e-28 | 5.66e-25 | 4.57e-05 | 0.172 | False | False | a_fragment_present |
| HOG0056747:g4 | Aparaglossata | 10 | 3 | 3 | 2 | 3:R>I;9:Y>A | 0.47 | 0.419 | 214 | 0.00465 | 1.1e-07 | 0.00035 | 1.1e-07 | 0.0021 | True | False |  |

## 4. Exploratory table

`specificity/tau0.2/specificity_exploratory.parquet`: all 536 funnel-stage real pairs with the cascade columns and flags `flag_bh_prereg_called`, `flag_tail_normal_called`, `flag_tail_gpd_called`, `flag_pass_artifacts`, `flag_pass_artifacts_strict`, `flag_ge3_species_d`, `flag_tail_and_cascade`. 230 pass the cascade (50 strict); 129 of those have >= 3 derived species. Cascade fail counts: {"a_fragment_present": 275, "a_orthology_panther": 57, "a_fragment": 34}. All rows are labelled exploratory, not significant. Null rows and per-input cascades are bulk: /tmp/pe_big/partB/specD34/O2_r1/.

## Caveats

- 3 null seeds (as the earlier specificity runs); empirical FDRs on null means of ~8-80 pairs are rough.
- Candidates files (`candidates.parquet`) were not read for writing; the cascade ran in scratch results roots with an emptied candidate table.
- O2 instability / species-tree lineage / PAE checks are `na` for specificity; O2 W members without an atlas site have no fingerprint.
- Design scope (D31): within-orthogroup sister pairs only.
