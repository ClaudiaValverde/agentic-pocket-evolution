# case_list_v2 (D25c/D28 co-fold pool), snapshot 20261030T085136Z

Headline: convincing_only (all three filters AND artifact cascade pass). Losses are not FDR-calibrated (D21).

## Empirical null FDR per category (mean null calls / real calls at the same cut, Part B)

| orthology:type | curveball | pseudo-clade | real called |
|---|---|---|---|
| O1:gain | 0.199 | 0.08 | 175 |
| O1:loss | 1.0 | 0.407 | 61 |
| O2:gain | 0.267 | 0.093 | 104 |
| O2:loss | 1.0 | 0.674 | 66 |

Category `duplication_gain_exploratory`: not preregistered; null is selection-conditioned; never in headline numbers.

## Selected

| case_id | orthology | category | selection | S | filter_label | f1 | f2 | f3 | f2+lit |
|---|---|---|---|---|---|---|---|---|---|
| control_falseloss_HOG0027591_g2_n66_O1 | O1 | control | negative_control | 0 | rejected | fail | fail | pending | nan |
| control_falseloss_HOG0028108_g15_n235_O1 | O1 | control | negative_control | 0 | rejected | fail | fail | pending | nan |
| control_falseloss_HOG0030328_g4_n103_O1 | O1 | control | negative_control | 0 | rejected | fail | fail | pending | nan |
| control_falseloss_HOG0034353_g3_n4659_O1 | O1 | control | negative_control | 0 | rejected | fail | fail | fail | nan |
| control_falseloss_HOG0035718_g5_n12_O1 | O1 | control | negative_control | 0 | rejected | fail | fail | pending | nan |
| dupgain_HOG0046679_g11_n33_O2 | O2 | duplication_gain_exploratory | exploratory_duplication | 0 | rejected | pass | fail | fail | fail |
| dupgain_HOG0054071_g9_n185_O2 | O2 | duplication_gain_exploratory | exploratory_duplication | 0 | rejected | pass | fail | fail | fail |
| gain_HOG0056135_g10_n70_O2 | O2 | gain | audit | 4 | supported | pass | fail | pass | fail |
| gain_HOG0034532_g4_n8_O2 | O2 | gain | top_k | 7 | supported | pass | fail | pass | fail |
| gain_HOG0035156_g10_n21_O1 | O1 | gain | top_k | 6 | supported | pass | fail | pass | fail |
| gain_HOG0025964_g2_n28_O1 | O1 | gain | top_k | 5 | supported | pass | fail | pass | fail |
| gain_HOG0026737_g5_n16_O1 | O1 | gain | top_k | 4 | supported | pass | fail | pass | fail |
| gain_HOG0035403_g5_n14_O1 | O1 | gain | top_k | 4 | supported | pass | fail | pass | fail |
| loss_HOG0027153_g11_n17_O1 | O1 | loss | top_k | 3 | supported | pass | fail | pass | fail |

## Counts

```
{
 "O1:control": {
  "n": 5,
  "f1": 0,
  "f2": 0,
  "f3": 0,
  "artifacts_pass": 0,
  "supported_or_convincing": 0,
  "convincing_headline": 0,
  "convincing_with_literature_post_hoc": 0,
  "eligible": 0,
  "duplicates": 0,
  "selected": 5
 },
 "O1:convergence": {
  "n": 3,
  "f1": 3,
  "f2": 0,
  "f3": 0,
  "artifacts_pass": 0,
  "supported_or_convincing": 0,
  "convincing_headline": 0,
  "convincing_with_literature_post_hoc": 0,
  "eligible": 0,
  "duplicates": 0,
  "selected": 0
 },
 "O1:gain": {
  "n": 50,
  "f1": 42,
  "f2": 1,
  "f3": 47,
  "artifacts_pass": 7,
  "supported_or_convincing": 39,
  "convincing_headline": 0,
  "convincing_with_literature_post_hoc": 0,
  "eligible": 5,
  "duplicates": 0,
  "selected": 4
 },
 "O1:loss": {
  "n": 17,
  "f1": 8,
  "f2": 2,
  "f3": 16,
  "artifacts_pass": 1,
  "supported_or_convincing": 8,
  "convincing_headline": 0,
  "convincing_with_literature_post_hoc": 0,
  "eligible": 1,
  "duplicates": 0,
  "selected": 1
 },
 "O1:specificity": {
  "n": 5,
  "f1": 2,
  "f2": 0,
  "f3": 0,
  "artifacts_pass": 0,
  "supported_or_convincing": 0,
  "convincing_headline": 0,
  "convincing_with_literature_post_hoc": 0,
  "eligible": 0,
  "duplicates": 0,
  "selected": 0
 },
 "O2:duplication_gain_exploratory": {
  "n": 2,
  "f1": 2,
  "f2": 0,
  "f3": 0,
  "artifacts_pass": 2,
  "supported_or_convincing": 0,
  "convincing_headline": 0,
  "convincing_with_literature_post_hoc": 0,
  "eligible": 0,
  "duplicates": 0,
  "selected": 2
 },
 "O2:gain": {
  "n": 104,
  "f1": 76,
  "f2": 0,
  "f3": 61,
  "artifacts_pass": 17,
  "supported_or_convincing": 47,
  "convincing_headline": 0,
  "convincing_with_literature_post_hoc": 0,
  "eligible": 8,
  "duplicates": 0,
  "selected": 2
 },
 "O2:loss": {
  "n": 66,
  "f1": 19,
  "f2": 3,
  "f3": 28,
  "artifacts_pass": 6,
  "supported_or_convincing": 7,
  "convincing_headline": 0,
  "convincing_with_literature_post_hoc": 0,
  "eligible": 0,
  "duplicates": 0,
  "selected": 0
 }
}
```
## Note (D30, added 2026-10-30 17:13 CET)
The M4 PANTHER/Pfam side-family check (`orthology_side_family`) was added after seeing results. It is **reported
only and does not enter selection**. The frozen run applied it as a gate, so it was re-run in report-only mode
(`--side-family-mode report`, now the default). That run gave the identical eligible sets (O1 6, O2 8) and the
identical selected cases: every supported/convincing case the check fails also fails `pass_artifacts` (Part B,
incl. the RBH `orthology_split` check). In the stage lists the check is still recorded as an artifact check.

## Note (17:19 CET): case files regenerated, selection unchanged
Bug fix (co-folding agent report): after dd2f3c1 the case.yaml files of the frozen list listed only P/W/A members, so
the controls lost the members that drove their false loss (now U). Every non-G member is now listed with `side`
(present/absent/other), state and reason. Controls also carry `state_pre_filter`, the state in the pre-D11 run that
called the false loss. The rerun gave the same selected case ids and selection labels (checked against 5e5a3ac).
Test: tests/test_case_files.py.
