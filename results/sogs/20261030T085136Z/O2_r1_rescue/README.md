# results/sogs/20261030T085136Z/O2_r1_rescue: BURROW rescue of O2_r1 A calls (τ = 0.2, k = 1)

Overlay on `../O2_r1` (its states are not modified). Script `scripts/rescue_sogs.py` (6 prereg settings via
`src/pocketevo/rescue.py`; hit = rescue pocket with residue Jaccard >= 0.3 against the member's footprint residues;
D9: k = 1, so A -> W iff any setting hits). CPU only (56 workers; GPUs untouched). Development subset, preliminary.

## Budget (`budget/`, 1,000 randomly sampled `rescue:pending` A calls, seed 0, 40 workers)
997 models; extraction 52 s (serial) + BURROW 73 s (early stop) -> about 0.125 s per model. All 442,783 pending A calls
(203,843 models) would take about 7 h at that rate, or about 4.4 h with parallel extraction. That is over the
90-minute limit, so only the Part B candidate pool was rescued, as the PI instructed: A calls in SOGs with P in
>= 3 species (100,967 rows, 80,418 models; projected 1.3–1.8 h). Sample result: 239/1,000 rescued.

## Run: candidate pool (`rescue_A.parquet`, `rescue_summary.json`)
Selection: `rescue:pending` A rows (model has atlas-filtered pockets) in SOGs with P in >= 3 species. SOGs flagged
small_footprint / empty_footprint would be skipped, but none of them have A calls (D11).
- Result: **W 50,308**, A 50,659 (`rescue:negative`) out of 100,967. No multi-fragment models were selected.
- By model tier (W / A): tier 1 21,511 / 21,538; tier 2 618 / 564; tier 3 9,172 / 9,782; tier L 19,007 / 18,775.
- Settings were run in order with early stop, so `n_settings_hit` is a lower bound and `hit_setting` is the first
  setting that hit. Rows first hit by each setting:
  | setting | rows hit / rows tested | models | wall |
  |---|---|---|---|
  | min_score_none | 36,209 / 100,967 | 80,418 | 780 s |
  | no_af_mode | 3,366 / 64,758 | 54,891 | 704 s |
  | seed_s1 | 1,306 / 61,392 | 52,330 | 593 s |
  | seed_s2 | 770 / 60,086 | 51,328 | 576 s |
  | preset_fast | 6,076 / 59,316 | 50,748 | 378 s |
  | preset_accurate | 2,581 / 53,240 | 46,012 | 1,185 s |
- Wall time 4,452 s (74 min): extraction 86 s, BURROW 4,216 s.
- Not rescued: the remaining 341,816 pending A calls outside the pool. They stay A with `rescue:pending` and must not
  be used as absences.
- Residual false-A after rescue at k = 1 (D9): about 0.13 overall and 0.054 for > 10-residue footprints. All A calls
  here satisfy the > 10-residue rule (D11).
