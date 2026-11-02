# Part B candidates — O1, snapshot 20261030T085136Z (development subset, preliminary)

Generated 2026-10-30 13:56 UTC by `scripts/partB/write_candidates.py`; prereg sha256 `52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83`; git `b2bb8f718761234093c47a9e8fb7597669362270+dirty`. Producer: phylogenetics agent; **not yet verified**. Gain/loss only (specificity needs footprint residues + fingerprints per member, not yet exported by the SOG stage).

## Pipeline

1. States: SOG agent, τ = 0.2 (primary) and 0.3 (sensitivity), primary A (pLDDT ≥ 90). D8: SOGs with footprint ≤ 3 columns excluded. **D11**: A requires the member footprint to span > 10 residues, else U (`scripts/partB/apply_d11.py`; A→U τ0.2: {'state': 556656, 'state_s70': 805246}, τ0.3: {'state': 1132315, 'state_s70': 1653454}).
2. Event-eligible SOGs (≥ 2 P and ≥ 2 A members; necessary for any candidate) → Sankoff on the O1 HOG gene tree, costs 2:1 (primary), 1:1, 3:1 → prereg gain/loss rules + FDR-calibrated cuts.
3. Rescue (k = 1, D9; hit = Jaccard ≥ 0.3, the rule k was calibrated with) of **every P and A member** of every SOG with a called candidate (A→W at ≥ 1/6 settings, P→W at 0/6), then events recomputed.

## Thresholds (re-derived on measured noise; event FDR 0.10)

Noise model (`results/partB/sim/measured/noise_primary_known_gt10.json`): per-species miss rate on known sites (> 10 residues, A_primary), rescue residual from error_rates.parquet (false-A after rescue k=1 = 0.054 for > 10 residues), **false-W on true absences 0.35** (control regions rescued), A-eligibility per species from the D11 states (pooled ≈ 1.4 % of non-P members are A → u_absent ≈ 0.98), false-P 0.02 (placeholder, not measured).

| model | cost | type | min_species | min_out_lineages | min_margin | recall | FDR |
|---|---|---|---|---|---|---|---|
| known-site miss (primary) | 1:1 | gain | 3 | 2 | 1.000 | 0.007 | 0.056 |
| known-site miss (primary) | 1:1 | loss | 2 | 2 | 2.000 | 0.000 | 0.000 |
| known-site miss (primary) | 2:1 | gain | 3 | 2 | 1.000 | 0.001 | 0.051 |
| known-site miss (primary) | 2:1 | loss | 2 | 3 | 2.000 | 0.000 | 0.000 |
| known-site miss (primary) | 3:1 | gain | 3 | 2 | 1.000 | 0.000 | 0.051 |
| known-site miss (primary) | 3:1 | loss | 2 | 3 | 2.000 | 0.000 | 0.000 |
| re-detection miss (sens.) | 1:1 | gain | 3 | 2 | 1.000 | 0.008 | 0.054 |
| re-detection miss (sens.) | 1:1 | loss | 3 | 4 | 1.000 | 0.000 | 0.000 |
| re-detection miss (sens.) | 2:1 | gain | 3 | 2 | 1.000 | 0.001 | 0.049 |
| re-detection miss (sens.) | 2:1 | loss | 2 | 3 | 2.000 | 0.000 | 0.000 |
| re-detection miss (sens.) | 3:1 | gain | 3 | 2 | 1.000 | 0.000 | 0.048 |
| re-detection miss (sens.) | 3:1 | loss | 2 | 3 | 2.000 | 0.000 | 0.000 |

Recall is per planted event on a random eligible branch of a real O1 tree; it is ≈ 0.1 % for gains and ≈ 0 for losses because true absences are almost never observable as A under the primary rule + D11 (u_absent ≈ 0.98) and 35 % of the observable ones are rescued to W. The loss cuts that meet the FDR target call essentially nothing in simulation: **real loss calls are not FDR-supported** and are listed for audit only.

## Funnel (real)

| tau | phase | cost | type | n_sogs | clade_size_ok | unambiguous | species | frac | outgroup | candidate | called | called_f1_single_event |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.2 | pre-rescue | 2:1 | gain | 9323 | 775 | 468 | 257 | 179 | 179 | 179 | 179 | 151 |
| 0.2 | pre-rescue | 2:1 | loss | 9323 | 6312 | 3125 | 3081 | 3060 | 2265 | 2265 | 63 | 32 |
| 0.2 | pre-rescue | 1:1 | gain | 9323 | 1468 | 835 | 559 | 503 | 503 | 503 | 503 | 379 |
| 0.2 | pre-rescue | 1:1 | loss | 9323 | 3725 | 943 | 931 | 860 | 860 | 860 | 9 | 6 |
| 0.2 | pre-rescue | 3:1 | gain | 9323 | 428 | 279 | 140 | 85 | 85 | 85 | 85 | 71 |
| 0.2 | pre-rescue | 3:1 | loss | 9323 | 8343 | 4874 | 4820 | 4807 | 3416 | 3416 | 79 | 35 |
| 0.2 | post-rescue | 2:1 | gain | 560 | 73 | 59 | 52 | 51 | 51 | 51 | 51 | 43 |
| 0.2 | post-rescue | 2:1 | loss | 560 | 341 | 175 | 174 | 170 | 106 | 106 | 18 | 9 |
| 0.2 | post-rescue | 1:1 | gain | 560 | 236 | 208 | 188 | 187 | 187 | 187 | 187 | 143 |
| 0.2 | post-rescue | 1:1 | loss | 560 | 169 | 77 | 77 | 64 | 49 | 49 | 5 | 4 |
| 0.2 | post-rescue | 3:1 | gain | 560 | 29 | 29 | 23 | 19 | 19 | 19 | 19 | 18 |
| 0.2 | post-rescue | 3:1 | loss | 560 | 430 | 271 | 271 | 270 | 143 | 143 | 19 | 10 |
| 0.3 | pre-rescue | 2:1 | gain | 13495 | 1571 | 894 | 500 | 366 | 366 | 366 | 366 | 312 |
| 0.3 | pre-rescue | 2:1 | loss | 13495 | 9020 | 4060 | 3988 | 3955 | 2787 | 2787 | 70 | 29 |
| 0.3 | pre-rescue | 1:1 | gain | 13495 | 2677 | 1406 | 913 | 821 | 821 | 821 | 821 | 643 |
| 0.3 | pre-rescue | 1:1 | loss | 13495 | 5059 | 1096 | 1081 | 985 | 985 | 985 | 13 | 10 |
| 0.3 | pre-rescue | 3:1 | gain | 13495 | 946 | 545 | 275 | 183 | 183 | 183 | 183 | 157 |
| 0.3 | pre-rescue | 3:1 | loss | 13495 | 12055 | 6462 | 6360 | 6345 | 4232 | 4232 | 108 | 35 |
| 0.3 | post-rescue | 2:1 | gain | 912 | 183 | 164 | 155 | 148 | 148 | 148 | 148 | 129 |
| 0.3 | post-rescue | 2:1 | loss | 912 | 582 | 287 | 284 | 279 | 187 | 187 | 24 | 14 |
| 0.3 | post-rescue | 1:1 | gain | 912 | 442 | 389 | 366 | 362 | 362 | 362 | 362 | 292 |
| 0.3 | post-rescue | 1:1 | loss | 912 | 313 | 151 | 150 | 126 | 104 | 104 | 8 | 7 |
| 0.3 | post-rescue | 3:1 | gain | 912 | 93 | 91 | 84 | 78 | 78 | 78 | 78 | 74 |
| 0.3 | post-rescue | 3:1 | loss | 912 | 779 | 453 | 452 | 451 | 269 | 269 | 34 | 15 |

Stage columns are cumulative: clade_size_ok = MPR event whose clade reaches the frozen species minimum; species/frac/outgroup = prereg gain/loss rules; candidate = all frozen rules; called = + calibrated cuts; called_f1 = + evidence filter 1 (single event spanning ≥ 2 species). Post-rescue rows cover only SOGs that had a called candidate before rescue.

Rescue τ0.2: 560 SOGs, 14460 proteins, A 7035→3660 (A→W 3375), P 8161→8030 (P→W 131), BURROW 3724 s.

SOGs dropped τ0.2: {'dropped_small_footprint_D8': 6595, 'dropped_lt2P_or_lt2A': 245073, 'og_not_in_orthology': 0}

Rescue τ0.3: 912 SOGs, 21429 proteins, A 14482→8676 (A→W 5806), P 10222→10222 (P→W 0), BURROW 1302 s.

SOGs dropped τ0.3: {'dropped_small_footprint_D8': 4531, 'dropped_lt2P_or_lt2A': 354673, 'og_not_in_orthology': 0}

## Real vs null inputs (pre-rescue D11 states, τ = 0.2, cost 2:1, same code and thresholds)

| input | seed | type | n_sogs | clade_size_ok | unambiguous | candidate | called | called_f1_single_event | note |
|---|---|---|---|---|---|---|---|---|---|
| real | – | gain | 9323 | 775.000 | 468.000 | 179.000 | 179.000 | 151.000 | – |
| real | – | loss | 9323 | 6312.000 | 3125.000 | 2265.000 | 63.000 | 32.000 | – |
| null_curveball | 0.000 | gain | 9323 | 454.000 | 323.000 | 32.000 | 32.000 | 17.000 | – |
| null_curveball | 0.000 | loss | 9323 | 5335.000 | 2724.000 | 2158.000 | 65.000 | 24.000 | – |
| null_pseudo_clade | 0.000 | gain | 9323 | 366.000 | 261.000 | 14.000 | 14.000 | 12.000 | – |
| null_pseudo_clade | 0.000 | loss | 9323 | 4344.000 | 2178.000 | 1747.000 | 23.000 | 5.000 | – |
| null_curveball | 1.000 | gain | 9323 | 434.000 | 298.000 | 27.000 | 27.000 | 20.000 | – |
| null_curveball | 1.000 | loss | 9323 | 5215.000 | 2734.000 | 2162.000 | 70.000 | 28.000 | – |
| null_pseudo_clade | 1.000 | gain | 9323 | 399.000 | 291.000 | 21.000 | 21.000 | 19.000 | – |
| null_pseudo_clade | 1.000 | loss | 9323 | 4177.000 | 2129.000 | 1722.000 | 25.000 | 9.000 | – |
| null_curveball | 2.000 | gain | 9323 | 427.000 | 288.000 | 29.000 | 29.000 | 17.000 | – |
| null_curveball | 2.000 | loss | 9323 | 5374.000 | 2866.000 | 2252.000 | 62.000 | 22.000 | – |
| null_pseudo_clade | 2.000 | gain | 9323 | 390.000 | 283.000 | 12.000 | 12.000 | 9.000 | – |
| null_pseudo_clade | 2.000 | loss | 9323 | 4228.000 | 2093.000 | 1680.000 | 23.000 | 13.000 | – |
| null_replicate | – | gain | 0 | – | – | – | – | – | not run: no replicate_states/ |
| null_replicate | – | loss | 0 | – | – | – | – | – | not run: no replicate_states/ |

### Real vs null at the same cut (D21)

| type | null | real_called | null_called_mean | null_called_max | real_f1 | null_f1_mean |
|---|---|---|---|---|---|---|
| gain | null_curveball | 179 | 29.333 | 32 | 151 | 18.000 |
| gain | null_pseudo_clade | 179 | 15.667 | 21 | 151 | 13.333 |
| loss | null_curveball | 63 | 65.667 | 70 | 32 | 24.667 |
| loss | null_pseudo_clade | 63 | 23.667 | 25 | 32 | 9.000 |

Each candidate row carries these counts (`null_*_called_mean_same_cut`, `real_called_same_cut_prerescue`). Pre-rescue comparison: the null inputs are permutations of the same pre-rescue D11 states.

## False-loss negative controls (losses overturned by rescue)

50 losses called before rescue (τ 0.2, 2:1) are not called after rescue. The five with the largest rescued fraction of clade A calls are the designated negative-control examples (`false_loss_controls.parquet`, `negative_control_example`). The other artifact checks were not run.

| negative_control_example | sog_id | level | n_A_before | n_A_after | n_W_after | frac_A_rescued | max_settings_hit | clade_species |
|---|---|---|---|---|---|---|---|---|
| True | HOG0035844:g7 | Eukaryota | 16 | 0 | 16 | 1.000 | 6 | UP000000437_7955;UP000000559_237561;UP000000589_10090;UP000001450_36329;UP000001631_447093;UP000001940_6239;UP000002059_502779;UP000002195_44689;UP000002311_559292;UP000002485_284812;UP000002494_10116;UP000005640_9606;UP000006672_6279;UP000008854_6183;UP000018087_1391915;UP000024404_6282;UP000030665_36087;UP000035681_6248;UP000053029_1442368;UP000078237_100816;UP000094526_86049;UP000270924_6293;UP000274756_318479 |
| True | HOG0035864:g3 | Mesangiospermae | 11 | 0 | 11 | 1.000 | 5 | UP000006548_3702;UP000007305_4577;UP000008827_3847;UP000059680_39947 |
| True | HOG0035379:g7 | Mesangiospermae | 9 | 0 | 9 | 1.000 | 6 | UP000006548_3702;UP000007305_4577;UP000008827_3847;UP000059680_39947 |
| True | HOG0036716:g3 | Rhabditida | 6 | 0 | 6 | 1.000 | 6 | UP000001940_6239;UP000006672_6279;UP000024404_6282;UP000035681_6248;UP000270924_6293;UP000274756_318479 |
| True | HOG0028153:g6 | rosids | 5 | 0 | 5 | 1.000 | 3 | UP000006548_3702;UP000008827_3847 |
| False | HOG0035431:g2 | Rhabditida | 5 | 0 | 5 | 1.000 | 4 | UP000001940_6239;UP000006672_6279;UP000024404_6282;UP000035681_6248;UP000270924_6293;UP000274756_318479 |
| False | HOG0033989:g0 | Enterobacteriaceae | 4 | 0 | 4 | 1.000 | 3 | UP000000625_83333;UP000001014_99287;UP000002716_300267;UP000007841_1125630 |
| False | HOG0035271:g3 | Enterobacteriaceae | 4 | 0 | 4 | 1.000 | 6 | UP000000625_83333;UP000001014_99287;UP000002716_300267;UP000007841_1125630 |

## Yield preflight (pi_freeze.yield_preflight)

Primary stratum (pLDDT ≥ 90 A, D11, τ 0.2, 2:1, post-rescue): **51 gains (FDR-calibrated) + 18 losses (not FDR-calibrated, D21) = 69 called** (52 also pass f1). Threshold 5 → secondary stratum not triggered.

## Primary candidates (top 25 by S_partial = persistence × min state confidence; evidence count pending filters 2/3)

| type | label | og_id | sog_id | level | clade_lca_species_tree | n_species_P | n_species_A | out_A_lineages | out_P_lineages | f1_single_event | score_S_partial | also_2to1_tau0.3 | also_1to1_tau0.2 | also_3to1_tau0.2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gain | FDR-calibrated candidate | HOG0028478 | HOG0028478:g6 | Bilateria | Bilateria | 10 | 0 | 4 | 0 | True | 10.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0035363 | HOG0035363:g6 | Opisthokonta | Opisthokonta | 9 | 0 | 3 | 0 | True | 9.000 | False | True | False |
| gain | FDR-calibrated candidate | HOG0035363 | HOG0035363:g0 | Opisthokonta | Opisthokonta | 13 | 0 | 3 | 0 | True | 6.500 | False | True | False |
| gain | FDR-calibrated candidate | HOG0026595 | HOG0026595:g6 | Ascomycota | Ascomycota | 6 | 0 | 3 | 0 | True | 6.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0029217 | HOG0029217:g0 | Rhabditida | Rhabditida | 6 | 0 | 5 | 0 | True | 6.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0036350 | HOG0036350:g6 | Rhabditida | Rhabditida | 6 | 0 | 6 | 0 | True | 6.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0025451 | HOG0025451:g0 | Ecdysozoa | Ecdysozoa | 5 | 0 | 3 | 0 | True | 5.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0025964 | HOG0025964:g2 | leotiomyceta | leotiomyceta | 5 | 1 | 3 | 0 | True | 5.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0030922 | HOG0030922:g6 | leotiomyceta | leotiomyceta | 5 | 0 | 5 | 0 | True | 5.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0026737 | HOG0026737:g5 | Spirurina | Spirurina | 4 | 0 | 6 | 0 | True | 4.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0029076 | HOG0029076:g3 | leotiomyceta | leotiomyceta | 4 | 0 | 3 | 0 | True | 4.000 | False | True | False |
| gain | FDR-calibrated candidate | HOG0029098 | HOG0029098:g4 | Rhabditida | Rhabditida | 4 | 0 | 5 | 0 | True | 4.000 | False | True | True |
| gain | FDR-calibrated candidate | HOG0029236 | HOG0029236:g20 | Rhabditida | Rhabditida | 4 | 0 | 5 | 1 | False | 4.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0029246 | HOG0029246:g3 | Eurotiomycetes | Eurotiomycetes | 4 | 0 | 8 | 0 | True | 4.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0029908 | HOG0029908:g16 | Bilateria | Bilateria | 4 | 0 | 5 | 0 | True | 4.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0030032 | HOG0030032:g5 | Rhabditida | Rhabditida | 4 | 0 | 5 | 0 | True | 4.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0030847 | HOG0030847:g4 | Rhabditida | Rhabditida | 4 | 0 | 8 | 0 | True | 4.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0031031 | HOG0031031:g3 | Euteleostomi | Euteleostomi | 4 | 0 | 4 | 0 | True | 4.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0033466 | HOG0033466:g4 | Enterobacteriaceae | Enterobacteriaceae | 4 | 0 | 3 | 0 | True | 4.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0034565 | HOG0034565:g13 | Euteleostomi | Euteleostomi | 4 | 0 | 5 | 0 | True | 4.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0034782 | HOG0034782:g3 | Mesangiospermae | Mesangiospermae | 4 | 0 | 29 | 9 | False | 4.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0034782 | HOG0034782:g3 | Mesangiospermae | Mesangiospermae | 4 | 0 | 29 | 9 | False | 4.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0034899 | HOG0034899:g8 | Pseudomonadati | Pseudomonadati | 4 | 0 | 3 | 0 | True | 4.000 | False | True | False |
| gain | FDR-calibrated candidate | HOG0035199 | HOG0035199:g6 | Bilateria | Bilateria | 4 | 0 | 4 | 0 | True | 4.000 | True | True | True |
| gain | FDR-calibrated candidate | HOG0035616 | HOG0035616:g5 | Enterobacteriaceae | Enterobacteriaceae | 4 | 0 | 4 | 0 | False | 4.000 | False | True | False |

Robustness of primary calls (fraction also called): {"also_1to1_tau0.2": 0.797, "also_1to1_tau0.3": 0.696, "also_3to1_tau0.2": 0.522, "also_3to1_tau0.3": 0.464, "also_2to1_tau0.2": 1.0, "also_2to1_tau0.3": 0.783}

## M3 tree consistency (D11 states, τ 0.2; `numbers.json` partA_descent.M3)

111685 multi-species SOGs: M3 = 0.932; checkerboard/curveball null 0.917 [0.917, 0.918] (1000 reps), p = 0.0010; noisy-vertical 0.848 ± 0.001. Swaps: 100 x min(#1, #0) (INTERP of prereg 100 x #1; see SIM.md).

## Caveats
- Development subset; O1 only (46 species). Deep (Eukaryota-level) events are flagged `deep` (orthology-uncertain).
- Artifact cascade (interface, TED/PAE, MSA-method flip, signal peptide already in states) not yet applied; required before any case is labelled.
- false-P is a placeholder; loss calls carry no FDR support.
