# Part B candidates — O1, snapshot 20261030T085136Z (development subset, preliminary)

Generated 2026-10-30 16:44 UTC by `scripts/partB/write_candidates.py`; prereg sha256 `52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83`; git `fe67f2c7270af9cc942c23de86049ca81ff53b28+dirty`. Producer: phylogenetics agent; **not yet verified**. Gain/loss only (specificity needs footprint residues + fingerprints per member, not yet exported by the SOG stage).

## Pipeline

1. States: SOG agent, τ = 0.2 (primary) and 0.3 (sensitivity), primary A (pLDDT ≥ 90). D8: SOGs with footprint ≤ 3 columns excluded. **D11**: A requires the member footprint to span > 10 residues, else U (applied by the SOG stage in release `O1_r3`; small_footprint / empty_footprint SOGs excluded).
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
| 0.2 | pre-rescue | 2:1 | gain | 9156 | 751 | 459 | 253 | 175 | 175 | 175 | 175 | 146 |
| 0.2 | pre-rescue | 2:1 | loss | 9156 | 6199 | 3060 | 3017 | 2995 | 2219 | 2219 | 61 | 31 |
| 0.2 | pre-rescue | 1:1 | gain | 9156 | 1440 | 821 | 550 | 494 | 494 | 494 | 494 | 371 |
| 0.2 | pre-rescue | 1:1 | loss | 9156 | 3673 | 938 | 925 | 854 | 854 | 854 | 9 | 6 |
| 0.2 | pre-rescue | 3:1 | gain | 9156 | 418 | 273 | 135 | 83 | 83 | 83 | 83 | 69 |
| 0.2 | pre-rescue | 3:1 | loss | 9156 | 8178 | 4766 | 4712 | 4699 | 3348 | 3348 | 76 | 33 |
| 0.2 | post-rescue | 2:1 | gain | 548 | 72 | 58 | 51 | 50 | 50 | 50 | 50 | 42 |
| 0.2 | post-rescue | 2:1 | loss | 548 | 338 | 178 | 177 | 173 | 109 | 109 | 17 | 8 |
| 0.2 | post-rescue | 1:1 | gain | 548 | 230 | 203 | 183 | 182 | 182 | 182 | 182 | 137 |
| 0.2 | post-rescue | 1:1 | loss | 548 | 165 | 75 | 75 | 62 | 48 | 48 | 5 | 4 |
| 0.2 | post-rescue | 3:1 | gain | 548 | 25 | 25 | 20 | 17 | 17 | 17 | 17 | 16 |
| 0.2 | post-rescue | 3:1 | loss | 548 | 430 | 273 | 273 | 272 | 146 | 146 | 18 | 9 |
| 0.3 | pre-rescue | 2:1 | gain | 13245 | 1536 | 879 | 494 | 360 | 360 | 360 | 360 | 305 |
| 0.3 | pre-rescue | 2:1 | loss | 13245 | 8860 | 3974 | 3905 | 3873 | 2730 | 2730 | 69 | 29 |
| 0.3 | pre-rescue | 1:1 | gain | 13245 | 2638 | 1387 | 894 | 803 | 803 | 803 | 803 | 632 |
| 0.3 | pre-rescue | 1:1 | loss | 13245 | 4954 | 1083 | 1067 | 971 | 971 | 971 | 13 | 10 |
| 0.3 | pre-rescue | 3:1 | gain | 13245 | 928 | 535 | 266 | 179 | 179 | 179 | 179 | 154 |
| 0.3 | pre-rescue | 3:1 | loss | 13245 | 11819 | 6330 | 6230 | 6215 | 4157 | 4157 | 107 | 35 |
| 0.3 | post-rescue | 2:1 | gain | 897 | 183 | 163 | 154 | 148 | 148 | 148 | 148 | 129 |
| 0.3 | post-rescue | 2:1 | loss | 897 | 570 | 282 | 279 | 274 | 185 | 185 | 24 | 14 |
| 0.3 | post-rescue | 1:1 | gain | 897 | 433 | 381 | 359 | 355 | 355 | 355 | 355 | 287 |
| 0.3 | post-rescue | 1:1 | loss | 897 | 309 | 146 | 145 | 121 | 99 | 99 | 8 | 7 |
| 0.3 | post-rescue | 3:1 | gain | 897 | 88 | 86 | 79 | 74 | 74 | 74 | 74 | 70 |
| 0.3 | post-rescue | 3:1 | loss | 897 | 774 | 452 | 451 | 450 | 269 | 269 | 34 | 15 |

Stage columns are cumulative: clade_size_ok = MPR event whose clade reaches the frozen species minimum; species/frac/outgroup = prereg gain/loss rules; candidate = all frozen rules; called = + calibrated cuts; called_f1 = + evidence filter 1 (single event spanning ≥ 2 species). Post-rescue rows cover only SOGs that had a called candidate before rescue.

Rescue τ0.2: 548 SOGs, 14261 proteins, A 6950→3612 (A→W 3338), P 8052→7921 (P→W 131), BURROW 33 s.

SOGs dropped τ0.2: {'dropped_small_footprint_D8': 6595, 'dropped_lt2P_or_lt2A': 245240, 'og_not_in_orthology': 0}

Rescue τ0.3: 897 SOGs, 21011 proteins, A 14153→8439 (A→W 5714), P 10069→10069 (P→W 0), BURROW 15 s.

SOGs dropped τ0.3: {'dropped_small_footprint_D8': 4531, 'dropped_lt2P_or_lt2A': 354923, 'og_not_in_orthology': 0}

## Real vs null inputs (pre-rescue D11 states, τ = 0.2, cost 2:1, same code and thresholds)

| input | seed | type | n_sogs | clade_size_ok | unambiguous | candidate | called | called_f1_single_event | note |
|---|---|---|---|---|---|---|---|---|---|
| real | – | gain | 9156 | 751.000 | 459.000 | 175.000 | 175.000 | 146.000 | – |
| real | – | loss | 9156 | 6199.000 | 3060.000 | 2219.000 | 61.000 | 31.000 | – |
| null_curveball | 0.000 | gain | 9156 | 397.000 | 285.000 | 28.000 | 28.000 | 16.000 | – |
| null_curveball | 0.000 | loss | 9156 | 5270.000 | 2698.000 | 2132.000 | 71.000 | 27.000 | – |
| null_pseudo_clade | 0.000 | gain | 9156 | 347.000 | 243.000 | 12.000 | 12.000 | 8.000 | – |
| null_pseudo_clade | 0.000 | loss | 9156 | 4237.000 | 2143.000 | 1739.000 | 29.000 | 6.000 | – |
| null_curveball | 1.000 | gain | 9156 | 421.000 | 300.000 | 39.000 | 39.000 | 22.000 | – |
| null_curveball | 1.000 | loss | 9156 | 5341.000 | 2742.000 | 2163.000 | 61.000 | 18.000 | – |
| null_pseudo_clade | 1.000 | gain | 9156 | 352.000 | 254.000 | 13.000 | 13.000 | 9.000 | – |
| null_pseudo_clade | 1.000 | loss | 9156 | 4230.000 | 2094.000 | 1638.000 | 20.000 | 9.000 | – |
| null_curveball | 2.000 | gain | 9156 | 427.000 | 302.000 | 35.000 | 35.000 | 22.000 | – |
| null_curveball | 2.000 | loss | 9156 | 5274.000 | 2720.000 | 2163.000 | 61.000 | 23.000 | – |
| null_pseudo_clade | 2.000 | gain | 9156 | 353.000 | 256.000 | 13.000 | 13.000 | 10.000 | – |
| null_pseudo_clade | 2.000 | loss | 9156 | 4233.000 | 2178.000 | 1744.000 | 24.000 | 12.000 | – |
| null_curveball | 3.000 | gain | 9156 | 453.000 | 310.000 | 39.000 | 39.000 | 25.000 | – |
| null_curveball | 3.000 | loss | 9156 | 5352.000 | 2788.000 | 2222.000 | 69.000 | 27.000 | – |
| null_pseudo_clade | 3.000 | gain | 9156 | 345.000 | 245.000 | 15.000 | 15.000 | 12.000 | – |
| null_pseudo_clade | 3.000 | loss | 9156 | 4278.000 | 2152.000 | 1705.000 | 25.000 | 12.000 | – |
| null_curveball | 4.000 | gain | 9156 | 404.000 | 274.000 | 35.000 | 35.000 | 21.000 | – |
| null_curveball | 4.000 | loss | 9156 | 5327.000 | 2747.000 | 2170.000 | 52.000 | 21.000 | – |
| null_pseudo_clade | 4.000 | gain | 9156 | 403.000 | 279.000 | 19.000 | 19.000 | 13.000 | – |
| null_pseudo_clade | 4.000 | loss | 9156 | 4348.000 | 2158.000 | 1728.000 | 22.000 | 9.000 | – |
| null_curveball | 5.000 | gain | 9156 | 414.000 | 300.000 | 34.000 | 34.000 | 17.000 | – |
| null_curveball | 5.000 | loss | 9156 | 5338.000 | 2722.000 | 2154.000 | 65.000 | 29.000 | – |
| null_pseudo_clade | 5.000 | gain | 9156 | 345.000 | 239.000 | 17.000 | 17.000 | 10.000 | – |
| null_pseudo_clade | 5.000 | loss | 9156 | 4290.000 | 2043.000 | 1637.000 | 25.000 | 9.000 | – |
| null_curveball | 6.000 | gain | 9156 | 433.000 | 334.000 | 38.000 | 38.000 | 26.000 | – |
| null_curveball | 6.000 | loss | 9156 | 5314.000 | 2732.000 | 2111.000 | 59.000 | 25.000 | – |
| null_pseudo_clade | 6.000 | gain | 9156 | 340.000 | 228.000 | 9.000 | 9.000 | 8.000 | – |
| null_pseudo_clade | 6.000 | loss | 9156 | 4249.000 | 2073.000 | 1659.000 | 26.000 | 9.000 | – |
| null_curveball | 7.000 | gain | 9156 | 458.000 | 313.000 | 47.000 | 47.000 | 34.000 | – |
| null_curveball | 7.000 | loss | 9156 | 5216.000 | 2757.000 | 2189.000 | 69.000 | 26.000 | – |
| null_pseudo_clade | 7.000 | gain | 9156 | 330.000 | 221.000 | 15.000 | 15.000 | 12.000 | – |
| null_pseudo_clade | 7.000 | loss | 9156 | 4337.000 | 2080.000 | 1658.000 | 21.000 | 7.000 | – |
| null_curveball | 8.000 | gain | 9156 | 396.000 | 270.000 | 32.000 | 32.000 | 22.000 | – |
| null_curveball | 8.000 | loss | 9156 | 5331.000 | 2781.000 | 2194.000 | 72.000 | 30.000 | – |
| null_pseudo_clade | 8.000 | gain | 9156 | 337.000 | 233.000 | 10.000 | 10.000 | 5.000 | – |
| null_pseudo_clade | 8.000 | loss | 9156 | 4258.000 | 2179.000 | 1701.000 | 26.000 | 7.000 | – |
| null_curveball | 9.000 | gain | 9156 | 419.000 | 290.000 | 40.000 | 40.000 | 22.000 | – |
| null_curveball | 9.000 | loss | 9156 | 5349.000 | 2766.000 | 2204.000 | 62.000 | 24.000 | – |
| null_pseudo_clade | 9.000 | gain | 9156 | 363.000 | 264.000 | 14.000 | 14.000 | 9.000 | – |
| null_pseudo_clade | 9.000 | loss | 9156 | 4318.000 | 2148.000 | 1686.000 | 33.000 | 13.000 | – |
| null_curveball | 10.000 | gain | 9156 | 409.000 | 282.000 | 27.000 | 27.000 | 15.000 | – |
| null_curveball | 10.000 | loss | 9156 | 5351.000 | 2836.000 | 2225.000 | 64.000 | 24.000 | – |
| null_pseudo_clade | 10.000 | gain | 9156 | 388.000 | 281.000 | 19.000 | 19.000 | 8.000 | – |
| null_pseudo_clade | 10.000 | loss | 9156 | 4316.000 | 2029.000 | 1638.000 | 23.000 | 10.000 | – |
| null_curveball | 11.000 | gain | 9156 | 413.000 | 285.000 | 27.000 | 27.000 | 20.000 | – |
| null_curveball | 11.000 | loss | 9156 | 5269.000 | 2751.000 | 2164.000 | 64.000 | 28.000 | – |
| null_pseudo_clade | 11.000 | gain | 9156 | 359.000 | 244.000 | 9.000 | 9.000 | 6.000 | – |
| null_pseudo_clade | 11.000 | loss | 9156 | 4118.000 | 2052.000 | 1620.000 | 17.000 | 6.000 | – |
| null_curveball | 12.000 | gain | 9156 | 397.000 | 273.000 | 25.000 | 25.000 | 13.000 | – |
| null_curveball | 12.000 | loss | 9156 | 5295.000 | 2784.000 | 2186.000 | 74.000 | 29.000 | – |
| null_pseudo_clade | 12.000 | gain | 9156 | 398.000 | 268.000 | 16.000 | 16.000 | 10.000 | – |
| null_pseudo_clade | 12.000 | loss | 9156 | 4145.000 | 2033.000 | 1604.000 | 25.000 | 9.000 | – |
| null_curveball | 13.000 | gain | 9156 | 427.000 | 307.000 | 43.000 | 43.000 | 27.000 | – |
| null_curveball | 13.000 | loss | 9156 | 5475.000 | 2879.000 | 2273.000 | 69.000 | 21.000 | – |
| null_pseudo_clade | 13.000 | gain | 9156 | 361.000 | 252.000 | 11.000 | 11.000 | 7.000 | – |
| null_pseudo_clade | 13.000 | loss | 9156 | 4253.000 | 2111.000 | 1683.000 | 27.000 | 5.000 | – |
| null_curveball | 14.000 | gain | 9156 | 427.000 | 296.000 | 42.000 | 42.000 | 26.000 | – |
| null_curveball | 14.000 | loss | 9156 | 5231.000 | 2709.000 | 2132.000 | 61.000 | 23.000 | – |
| null_pseudo_clade | 14.000 | gain | 9156 | 323.000 | 231.000 | 21.000 | 21.000 | 13.000 | – |
| null_pseudo_clade | 14.000 | loss | 9156 | 4218.000 | 2103.000 | 1664.000 | 28.000 | 10.000 | – |
| null_curveball | 15.000 | gain | 9156 | 425.000 | 306.000 | 31.000 | 31.000 | 19.000 | – |
| null_curveball | 15.000 | loss | 9156 | 5355.000 | 2836.000 | 2251.000 | 71.000 | 25.000 | – |
| null_pseudo_clade | 15.000 | gain | 9156 | 388.000 | 251.000 | 8.000 | 8.000 | 4.000 | – |
| null_pseudo_clade | 15.000 | loss | 9156 | 4247.000 | 2116.000 | 1698.000 | 23.000 | 4.000 | – |
| null_curveball | 16.000 | gain | 9156 | 432.000 | 306.000 | 40.000 | 40.000 | 22.000 | – |
| null_curveball | 16.000 | loss | 9156 | 5298.000 | 2778.000 | 2231.000 | 76.000 | 31.000 | – |
| null_pseudo_clade | 16.000 | gain | 9156 | 343.000 | 242.000 | 14.000 | 14.000 | 13.000 | – |
| null_pseudo_clade | 16.000 | loss | 9156 | 4240.000 | 2110.000 | 1681.000 | 24.000 | 9.000 | – |
| null_curveball | 17.000 | gain | 9156 | 404.000 | 290.000 | 27.000 | 27.000 | 16.000 | – |
| null_curveball | 17.000 | loss | 9156 | 5248.000 | 2706.000 | 2158.000 | 75.000 | 24.000 | – |
| null_pseudo_clade | 17.000 | gain | 9156 | 393.000 | 276.000 | 16.000 | 16.000 | 13.000 | – |
| null_pseudo_clade | 17.000 | loss | 9156 | 4309.000 | 2098.000 | 1692.000 | 25.000 | 14.000 | – |
| null_curveball | 18.000 | gain | 9156 | 419.000 | 288.000 | 32.000 | 32.000 | 23.000 | – |
| null_curveball | 18.000 | loss | 9156 | 5261.000 | 2784.000 | 2178.000 | 57.000 | 18.000 | – |
| null_pseudo_clade | 18.000 | gain | 9156 | 352.000 | 245.000 | 12.000 | 12.000 | 7.000 | – |
| null_pseudo_clade | 18.000 | loss | 9156 | 4260.000 | 2091.000 | 1645.000 | 23.000 | 5.000 | – |
| null_curveball | 19.000 | gain | 9156 | 425.000 | 312.000 | 36.000 | 36.000 | 27.000 | – |
| null_curveball | 19.000 | loss | 9156 | 5265.000 | 2737.000 | 2180.000 | 68.000 | 24.000 | – |
| null_pseudo_clade | 19.000 | gain | 9156 | 400.000 | 283.000 | 17.000 | 17.000 | 10.000 | – |
| null_pseudo_clade | 19.000 | loss | 9156 | 4229.000 | 2024.000 | 1607.000 | 30.000 | 11.000 | – |
| null_replicate | – | gain | 0 | – | – | – | – | – | not run: no replicate_states/ |
| null_replicate | – | loss | 0 | – | – | – | – | – | not run: no replicate_states/ |

### Real vs null at the same cut (D21, D24): empirical null FDR = mean null calls / real calls

| type | null | real_called | null_called_mean | null_called_max | real_f1 | null_f1_mean | n_null_seeds | empirical_null_FDR | simulated_FDR |
|---|---|---|---|---|---|---|---|---|---|
| gain | null_curveball | 175 | 34.850 | 47 | 146 | 21.750 | 20 | 0.199 | 0.051 |
| gain | null_pseudo_clade | 175 | 14.000 | 21 | 146 | 9.350 | 20 | 0.080 | 0.051 |
| loss | null_curveball | 61 | 66.000 | 76 | 31 | 24.850 | 20 | 1.000 | – |
| loss | null_pseudo_clade | 61 | 24.800 | 33 | 31 | 8.800 | 20 | 0.407 | – |

Each candidate row carries these counts (`null_*_called_mean_same_cut`, `real_called_same_cut_prerescue`). Pre-rescue comparison: the null inputs are permutations of the same pre-rescue D11 states.

## False-loss negative controls (D24)

**bad_model** (319 found; 5 marked): losses called on states before the quality filters (O1_r2 without D11; primary and pLDDT >= 70 strata) whose derived-clade A members are all poor models and whose call disappears once A_primary + D11 apply.

| negative_control_example | pre_stratum | sog_id | level | n_A_clade | frac_A_bad | min_A_footprint_plddt | bad_reasons | frac_A_rescued | clade_species |
|---|---|---|---|---|---|---|---|---|---|
| True | secondary_pLDDT70 | HOG0025405:g3 | Rhabditida | 4 | 1.000 | 72.750 | footprint_pLDDT<90 | – | UP000001940_6239;UP000006672_6279;UP000024404_6282;UP000035681_6248;UP000270924_6293;UP000274756_318479 |
| True | secondary_pLDDT70 | HOG0027801:g5 | Onchocercidae | 3 | 1.000 | 74.000 | footprint_pLDDT<90 | – | UP000006672_6279;UP000024404_6282;UP000270924_6293 |
| True | secondary_pLDDT70 | HOG0026077:g6 | Onchocercidae | 3 | 1.000 | 74.650 | footprint_pLDDT<90 | – | UP000006672_6279;UP000024404_6282;UP000270924_6293 |
| True | primary_pLDDT90 | HOG0025884:g3 | Euteleostomi | 8 | 1.000 | 90.250 | footprint<=10res | – | UP000000437_7955;UP000000589_10090;UP000002494_10116;UP000005640_9606 |
| True | primary_pLDDT90 | HOG0033916:g3 | Enterobacteriaceae | 3 | 1.000 | 90.250 | footprint<=10res | – | UP000000625_83333;UP000001014_99287;UP000002716_300267;UP000007841_1125630 |

**rescue_instability** (48 found; 5 marked): losses called before rescue on the release and not after rescue (high-quality models; pocket-finder instability).

| negative_control_example | pre_stratum | sog_id | level | n_A_clade | frac_A_bad | min_A_footprint_plddt | bad_reasons | frac_A_rescued | clade_species |
|---|---|---|---|---|---|---|---|---|---|
| True | O1_r3 primary (D11) | HOG0035844:g7 | Eukaryota | 16 | – | – | – | 1.000 | UP000000437_7955;UP000000559_237561;UP000000589_10090;UP000001450_36329;UP000001631_447093;UP000001940_6239;UP000002059_502779;UP000002195_44689;UP000002311_559292;UP000002485_284812;UP000002494_10116;UP000005640_9606;UP000006672_6279;UP000008854_6183;UP000018087_1391915;UP000024404_6282;UP000030665_36087;UP000035681_6248;UP000053029_1442368;UP000078237_100816;UP000094526_86049;UP000270924_6293;UP000274756_318479 |
| True | O1_r3 primary (D11) | HOG0035864:g3 | Mesangiospermae | 11 | – | – | – | 1.000 | UP000006548_3702;UP000007305_4577;UP000008827_3847;UP000059680_39947 |
| True | O1_r3 primary (D11) | HOG0035379:g7 | Mesangiospermae | 9 | – | – | – | 1.000 | UP000006548_3702;UP000007305_4577;UP000008827_3847;UP000059680_39947 |
| True | O1_r3 primary (D11) | HOG0036716:g3 | Rhabditida | 6 | – | – | – | 1.000 | UP000001940_6239;UP000006672_6279;UP000024404_6282;UP000035681_6248;UP000270924_6293;UP000274756_318479 |
| True | O1_r3 primary (D11) | HOG0028153:g6 | rosids | 5 | – | – | – | 1.000 | UP000006548_3702;UP000008827_3847 |

## Yield preflight (pi_freeze.yield_preflight)

Primary stratum (pLDDT ≥ 90 A, D11, τ 0.2, 2:1, post-rescue): **50 gains (FDR-calibrated) + 17 losses (not FDR-calibrated, D21) = 67 called** (50 also pass f1). Threshold 5 → secondary stratum not triggered.

## Primary candidates (top 25 by S_partial = persistence × min state confidence; evidence count pending filters 2/3)

| type | label | og_id | sog_id | level | clade_lca_species_tree | n_species_P | n_species_A | out_A_lineages | out_P_lineages | f1_single_event | score_S_partial | also_2to1_tau0.3 | also_1to1_tau0.2 | also_3to1_tau0.2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gain | FDR-calibrated candidate | HOG0028478 | HOG0028478:g6 | Bilateria | Bilateria | 10 | 0 | 3 | 0 | True | 10.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0035363 | HOG0035363:g6 | Opisthokonta | Opisthokonta | 9 | 0 | 3 | 0 | True | 9.000 | False | True | False |
| gain | FDR-calibrated candidate | HOG0035363 | HOG0035363:g0 | Opisthokonta | Opisthokonta | 13 | 0 | 3 | 0 | True | 6.500 | False | True | False |
| gain | FDR-calibrated candidate | HOG0026595 | HOG0026595:g6 | Ascomycota | Ascomycota | 6 | 0 | 3 | 0 | True | 6.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0029217 | HOG0029217:g0 | Rhabditida | Rhabditida | 6 | 0 | 5 | 0 | True | 6.000 | True | True | False |
| gain | FDR-calibrated candidate | HOG0036350 | HOG0036350:g6 | Rhabditida | Rhabditida | 6 | 0 | 4 | 0 | True | 6.000 | True | True | True |
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

Robustness of primary calls (fraction also called): {"also_1to1_tau0.2": 0.806, "also_1to1_tau0.3": 0.701, "also_3to1_tau0.2": 0.493, "also_3to1_tau0.3": 0.433, "also_2to1_tau0.2": 1.0, "also_2to1_tau0.3": 0.791}

## M3 tree consistency (D11 states, τ 0.2; `numbers.json` partA_descent.M3)

111685 multi-species SOGs: M3 = 0.933; checkerboard/curveball null 0.919 [0.918, 0.919] (1000 reps), p = 0.0010; noisy-vertical 0.942 ± 0.001. Swaps: 100 x min(#1, #0) (INTERP of prereg 100 x #1; see SIM.md).

### D22 swap-mixing check (null M3 at 50x / 100x / 200x min(#1, #0))

| swap_mult | replicates | swaps | null_mean | null_q025 | null_q975 | m3_real |
|---|---|---|---|---|---|---|
| 50.000 | 200.000 | 2245700.000 | 0.919 | 0.918 | 0.919 | 0.933 |
| 100.000 | 1000.000 | 4491400.000 | 0.919 | 0.918 | 0.919 | 0.933 |
| 200.000 | 200.000 | 8982800.000 | 0.919 | 0.918 | 0.919 | 0.933 |

A stable null mean across multipliers means the swap chain has mixed.

## Caveats
- Development subset; O1 only (46 species). Deep (Eukaryota-level) events are flagged `deep` (orthology-uncertain).
- Artifact cascade (interface, TED/PAE, MSA-method flip, signal peptide already in states) not yet applied; required before any case is labelled.
- false-P is a placeholder; loss calls carry no FDR support.

## Artifact cascade (run_artifacts.py; one column per check, nothing dropped)

Checks (fail -> pass_artifacts false): a_footprint_plddt, a_disorder, a_pae_ted, a_fragment, a_fragment_present, a_orthology_o2, a_orthology_panther, a_instability, a_species_tree_lineages, orthology_split. Flags (reported, not fails): flag_paralog_branch, flag_present_plddt_lt90. PAE was not fetched: the TED multi-domain rule turns such A calls into U. `na` = no data for that check (listed in `artifact_na`).

| type | n | pass_artifacts | pass_strict | paralog_branch | O2 paralog split (T/F/na) | present pLDDT<90 | fail counts |
|---|---|---|---|---|---|---|---|
| gain | 50 | 7 | 4 | 32 | 23/9/0 | 25 | {"a_footprint_plddt": 1, "a_disorder": 2, "a_pae_ted": 11, "a_fragment": 6, "a_fragment_present": 30, "a_orthology_o2": 27, "a_orthology_panther": 29, "a_species_tree_lineages": 5, "orthology_split": 41} |
| loss | 17 | 1 | 0 | 11 | 8/3/0 | 13 | {"a_pae_ted": 9, "a_fragment": 3, "a_fragment_present": 10, "a_orthology_o2": 5, "a_orthology_panther": 3, "a_species_tree_lineages": 6, "orthology_split": 15} |
| specificity | 5 | 0 | 0 | 3 | 0/0/0 | 2 | {"a_fragment": 1, "a_fragment_present": 3, "a_orthology_o2": 1, "a_orthology_panther": 3} |

## Specificity change (residue-first, prereg; empirical matched-depth null + BH)

| input | seed | n_pairs | residue_gate | residue_and_shift | bh_tested | called | called_f1_ge3_species |
|---|---|---|---|---|---|---|---|
| real | nan | 222489 | 10584 | 513 | 10584 | 0 | 0 |
| null_pseudo_clade | 0.0 | 222489 | 2230 | 88 | 2230 | 0 | 0 |
| null_pseudo_clade | 1.0 | 222489 | 2260 | 77 | 2260 | 0 | 0 |
| null_pseudo_clade | 2.0 | 222489 | 2243 | 71 | 2243 | 0 | 0 |

**No pair passes BH (negative result, 3 null seeds).** Pairs with >= 2 BLOSUM62<0 clade-consistent substitutions AND a fingerprint shift above the noise floor are listed as **exploratory, not significant** (`specificity/tau0.2/specificity_exploratory.parquet`); the top 5 (smallest q) are in `specificity_candidates.parquet` with artifact columns for the case builder's audit. No claims.

| sog_id | level | n_d | n_s | n_species_d | n_subst | substitutions | shift | p_empirical | q_bh | pass_artifacts | artifact_fails |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HOG0025600:g5 | Ecdysozoa | 3 | 2 | 3 | 3 | 9:G>D;10:G>I;11:F>A | 1.16 | 0.000128 | 0.289 | False | a_fragment_present;a_orthology_panther |
| HOG0025600:g5 | Ecdysozoa | 2 | 3 | 2 | 3 | 9:D>G;10:I>G;11:A>F | 1.16 | 0.000128 | 0.289 | False | a_orthology_panther |
| HOG0027899:g7 | Eukaryota | 10 | 1 | 10 | 4 | 0:G>R;2:G>L;3:L>A;7:P>A | 1 | 0.000397 | 0.289 | False | a_fragment;a_fragment_present |
| HOG0028186:g4 | Eukaryota | 8 | 1 | 4 | 2 | 5:F>E;7:N>A | 0.998 | 0.000397 | 0.289 | False | a_fragment_present |
| HOG0027961:g4 | Mesangiospermae | 4 | 7 | 4 | 2 | 3:R>I;5:W>R | 0.973 | 0.000223 | 0.289 | False | a_orthology_o2;a_orthology_panther |

Primary calls only (τ 0.2, 2:1); specificity = BH-called pairs (or the exploratory top 5). Per-candidate columns are in `candidates.parquet` / `specificity_candidates.parquet`.

