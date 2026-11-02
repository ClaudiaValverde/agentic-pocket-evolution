# results/sogs/20261030T085136Z/O2_r1_split: leakage-safe O2 split (PI decision D20)

States and SOGs of `../O2_r1` are unchanged; this is a split table next to it (`og_split_o2.parquet`,
`scripts/sog_split_o2.py`). An O2 root OG is `test` if its O2 hash half is test and no member belonged to a run O1 TUNE
family (O1_r2 inputs, >= 3 species); `tune` symmetrically (no O1-test member); otherwise `crossover`, excluded from
both tune and test reporting.

Eligible O2 OGs (>= 3 species, 35,093): tune 12,099 · test 12,220 · **crossover 10,774** (5,374 hash-test and 5,400
hash-tune). Over all 57,332 OGs: crossover 11,413. 5 control-list (forced-test) families are crossover.
Most annotated families are large and crossed: benchmark pairs (pos/neg) tune 4,399/1,196 · test 3,020/1,177 ·
crossover 30,727/6,961.

§8.1 benchmark re-reported on the clean halves (`tau_benchmark_o2_split.csv`; primary subset, CI = 200-replicate
bootstrap over OGs). τ is not re-chosen (frozen at 0.2, D8):

| τ | clean tune F1 [CI] | clean test F1 [CI] | crossover F1 |
|---|---|---|---|
| 0.2 | 0.723 [0.603, 0.809] | **0.695 [0.618, 0.775]** (R 0.633, P 0.771; 866 pos / 984 neg, 230 OGs) | 0.691 |
| 0.3 | 0.697 | 0.699 [0.638, 0.774] | 0.681 |
| 0.4 | 0.605 | 0.629 | 0.642 |
| 0.5 | 0.556 | 0.569 | 0.571 |
