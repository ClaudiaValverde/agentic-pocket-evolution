# Verification G4b: Part B on O2 / D28 cascade, validation track, cross-OG extension, REPORT v0

Verifier: independent agent (did not produce these results). Started 2026-10-30 23:07 CET, written 23:21 CET (`date`).
Runtime about 15 min wall (the longest step was an independent DIAMOND re-run of 10 RBH rows, 4 min). Scripts are in
`/tmp/verify_g4b/` (not committed). All commands were run from the repo root with `.venv/bin/python`. Files were read only.
The run was interrupted by a session restart at about 23:12; the O2 rescue check was re-run from scratch as `o2_rescue_v2.py`.

## Verdict summary

| # | claim | verdict |
|---|---|---|
| A1 | O2 empirical null FDR, gains: curveball 0.267, pseudo-clade 0.093 | PASS-with-notes |
| A2 | RBH logic (10 sampled supporters): best hit falls outside the HOG | FAIL as phrased; code is correct for what it tests |
| A3 | "617 of 782 outgroup A members in a different O2 root HOG" | PASS-with-notes (wrong label; the 79 % holds) |
| A4 | pass_artifacts O1 7/50, O2 17/104 | PASS (D28 ledger text not followed; see notes) |
| A5 | O2 rescue: A outside the rescue pool are U | PASS-with-notes (only "pending" A; 105,278 A kept as A) |
| B | validation table, 7 known events + 5 negative controls | PASS-with-notes for each row; FAIL for the headline reason in REPORT §0.5 and D31 |
| C1 | method commit 77db16c precedes all cross-OG outputs | PASS |
| C2 | 0/7 computed as specified | PASS |
| C3 | 610 vs 1 links; (b) 10 vs 33 | PASS-with-notes ("528 SOGs" counts one side only) |
| D | `check_report_numbers.py REPORT.md` | PASS (303 numbers, 0 unmatched) |
| D-spot | 10 claims checked for meaning | 5 PASS, 3 PASS-with-notes, 2 FAIL |

## A. Part B on O2 and the D28 cascade

### A1 O2 empirical null FDR: PASS-with-notes
`null_fdr.py` on `results/partB/nulls/20261030T085136Z/O2_r1/tau0.2/primary/funnel_real_vs_null.csv`:
- gain curveball: mean 27.8 null calls / 104 real = **0.2673** (median 0.255, range 17–45, 20 seeds)
- gain pseudo-clade: 9.65 / 104 = **0.0928** (range 1–14, 20 seeds)
- loss: 102.0 / 66 (curveball), 44.5 / 66 = 0.674 (pseudo-clade). numbers.json caps the curveball loss FDR at 1.0.

These match `partB.O2_funnel_and_null.empirical_null_FDR` in numbers.json exactly.

Notes:
1. Provenance. The sidecar and `manifest.json` record `seeds: [0,1,2]`, the prereg seeds, but argv has `--n-seeds 20` and
   the CSV holds seeds 0–19. `run_nulls.py:254,261` writes `seeds=seeds` rather than `null_seeds`. This is cosmetic, but
   a reader of the manifest would think 3 seeds were used.
2. The two FDRs are not comparable across orthologies. The O2 real and null runs use the A-rescued states
   (`states_resc_t0.2`), while the O1 0.199 is computed on the 175 pre-rescue calls. The REPORT §3.6 column header
   "called (pre-rescue)" is therefore wrong for O2: numbers.json itself calls 104 `called_post_rescue`.
3. `m3.json` in the same directory has 1 curveball replicate, with swaps_done 5.9 M against a target of 10.6 M. That is
   not used for the FDR, but M3 on O2 is not a valid result.

### A2 RBH logic (orthology_split): code is correct; the stated interpretation FAILS
How the code works (`rbh_check.py`): the forward best hit of the outgroup supporter q is searched over the
derived-clade species (q's own species excluded). The check passes only if that hit is a **derived-clade member** of the
same HOG tree and its reverse best hit in q's proteome is q. `run_artifacts.py` fails a call when fewer than 50 % of its
testable supporters pass.

Independent re-run: `rbh_rerun.py` takes 10 random O2 gain supporters (seed 7), extracts fresh queries from the
proteome FASTAs, and runs DIAMOND 2.2.8 `--sensitive --max-target-seqs 5` against each derived species. It reproduces
the stored forward best hit in **10/10** cases and the stored RBH = False in **10/10**. In 3 rows the reverse hit
differs (stored NaN, mine gives a hit), which does not change the verdict because the forward step already fails.

The question asked was whether the best hit really falls **outside the same HOG**. Mostly it does **not**. In **9 of 10**
sampled rows the best hit is in the same O2 HOG tree, just outside the derived clade. In 1 of 10 (F2U028) it is in
another root HOG.

Genome-wide check (`rbh_hog.py`, all 3,867 O2 rows): among 2,649 forward failures that have a hit, **2,335 (88 %) hit
the same O2 root HOG** and 314 (12 %) hit another root. In 1,408 of them the best hit is itself an outgroup supporter of
the same call.

So the O2 RBH check mostly detects within-HOG paralogy: the absent outgroups are orthologs of the sister copy, and the
gain sits on a paralog branch (consistent with 57/104 paralog-branch flags). It does not detect HOG misassignment.
REPORT §3.6 ("orthologs of another paralog copy") is the right reading. Any wording that says RBH shows the outgroups
are "in a different family/HOG" would be wrong.

The ledger is also out of date. D28 says the O2 check is a PANTHER-family split, but the code (commit 0cf2502,
"D28 revised") gates on RBH instead. As implemented, `panther_split` is `unknown` for all 104 O2 gains because it
requires every member to be annotated, and it does not enter `pass_artifacts`. **No D-entry records the switch to RBH**
(`grep -i "RBH\|reciprocal" DESIGN.md` returns nothing). This needs a ledger line.

### A3 617 / 782: PASS-with-notes
617 + 37 + 34 + 94 = 782, so the JSON is internally consistent. But `o1_o2_agreement.py` labels as "in another O2 root
HOG" what is actually `state2 == "not_in_O2_SOG_family"`: the protein has no row in the modal O2 **SOG** of the derived
P members. That is not the same as being in a different root HOG. The denominator also only covers O1 SOGs whose P
members map to a raw O2 SOG (inner merge).

Independent recount (`o1_root.py`): 50 O1 primary gains, their outgroup A members (states from the O1_r3 rescue, O1
trees), and O2 roots from `protein_hogs.parquet`. Result: **734 / 923 (79.5 %)** are in a different O2 root HOG than
the modal root of the derived P members, and 0 are missing from O2. 34 of 50 calls have more than 50 % in another root,
which matches `o2_root_split` fail = 34. The 79 % in D28 holds. The "617 of 782" pair should either be relabelled
("not in the derived clade's modal O2 SOG") or replaced by 734/923.

### A4 pass_artifacts counts: PASS
`pass_counts.py` recomputes `pass_artifacts = no CHECK == fail` from `candidates.parquet` (primary rows, no duplicate
keys):

| run | type | n | pass_artifacts | without orthology_split |
|---|---|---|---|---|
| O1 | gain | 50 | **7** | 9 |
| O1 | loss | 17 | 1 | — |
| O2 | gain | 104 | **17** | 34 |
| O2 | loss | 66 | 6 | — |

All match the stored values and `artifact_summary.json`. Strict counts are 4 and 13, and the paralog flags are 32 and
57, both as in the REPORT.

If D28 had been applied literally (O1: `o2_root_split` gate; O2: `panther_split` gate) instead of RBH, the counts would
be O1 7 and O2 34, because on O2 the PANTHER gate never fires. So the O2 count of 17 depends entirely on the
unrecorded RBH substitution (see A2).

### A5 O2 rescue, A outside the pool: PASS-with-notes
`o2_rescue_v2.py` compares every row of the raw `states_t0.2` against `states_resc_t0.2` across 1,024 shards (row
order checked). Every A row follows the rule in `prep_o2_states.py`:

| category | final state | rows |
|---|---|---|
| in overlay | A | 50,678 |
| in overlay | W | 50,309 |
| not in overlay, `n_filtered_pockets > 0` (pending) | U | 341,796 |
| not in overlay, 0 filtered pockets | A (kept) | 105,278 |

No non-A row changed. So A outside the pool become U **only if pending**. The 105,278 A with no atlas-filtered pockets
stay A without rescue, which is disclosed in CANDIDATES.md.

This matters for the results. Of the 2,570 outgroup-A supporters behind O2 gains, **594 (23 %) are such unrescued A**,
and the other 1,976 are overlay-confirmed A. On O1, rescue was run on every member of candidate SOGs. The O2 rescue is
therefore less complete than O1's for the A calls that support gains. P is not rescued on O2 (disclosed).

## B. Validation track (validation_table.md / .parquet, cases/val_*)

OG membership was checked independently: `ox_roots.py` parses both orthoXMLs and takes the top-level group ids, and
`ox_human.py` lists human members of a root group. Both agree with the `ogs` column for all 25 named proteins.

| item | O1 | O2 | verdict on stated reason |
|---|---|---|---|
| Steroid receptor (ER ancestral → AR/GR/MR/PR) | ER P03372, AR P10275, GR P04150 all in **HOG0025274**; MR P08235 HOG0022634; PR P06401 HOG0022342 | ER HOG0035169; **AR, GR, MR, PR all in HOG0035184** (sub-HOGs .35/.1/.24/.13) | O2 PASS. O1 PASS-with-notes: ER vs AR/GR is within one O1 OG, so "orthology split" is only part of the reason; the binding reason is "events but none called". **MR and GR are in different O1 root HOGs but the SAME O2 root HOG** |
| NR4A/Nurr1 loss | all 4 in HOG0025274 (same_og True) | DHR38 P41829 HOG0028532 vs NR4A HOG0035741 | PASS (O1: no A, not event-eligible; O2: split + not eligible) |
| PfLDH vs MDH | PfLDH Q76NM3 and PfMDH C6KT25 both HOG0035618; human MDH1 P40925 HOG0030579 | Q76NM3 and C6KT25 both HOG0019046; P40925 HOG0051690 | PASS-with-notes: the derived PfLDH and the ancestral PfMDH share an OG in both orthologies, so the binding reason is "SOG not event-eligible / no A", not orthology |
| IDH NADP→NAD | P08200 HOG0036048 vs P50213 HOG0036184 | HOG0054339 vs HOG0053286 | PASS |
| GID1 vs HSL | GID1 HOG0030582 vs HSL HOG0029749 | HOG0044281 vs HOG0042820 | PASS |
| D14 vs KAI2 | all HOG0034565 | all HOG0044130 | PASS (within-OG; specificity 0 / not run) |
| Lysozyme vs α-lactalbumin | all HOG0027604 | all HOG0047027 | PASS (within-OG) |
| Uricase (neg) | HOG0035545: 0 human members | HOG0054102: 0 human members | PASS |
| GULO (neg) | HOG0034483: 0 human members | HOG0052278: human LDHD Q86WU2 and D2HGDH Q8N465 (paralogs) | PASS-with-notes: on O2 the human rows are P/U (paralogs), not G. "Never A" holds, and there is no called loss in that OG. One called O2 *gain* in HOG0052278 includes human; it is irrelevant to the control |
| M. leprae truncation (neg) | A-rate leprae < tuberculosis; 0 called | — | PASS (numbers read from the table; not recomputed) |
| 3 detection-limited (proteases, CA, SOD) | reported, not counted | — | PASS (not counted as controls) |

**Headline-reason FAIL.** D31 and REPORT §0 point 5 say recall is 0/7 "because they lie between paralogous orthogroups".
By the table's own `same_og`, D14/KAI2 and lysozyme/α-lactalbumin are within one OG on both orthologies, and Nurr1 is
within one OG on O1. ER vs AR/GR and PfLDH vs PfMDH also share an OG on both. D31 cites "MR/GR" as a
between-orthogroup case, but on O2 MR and GR are in the same root HOG (HOG0035184). Orthology split is the binding
reason for at most IDH, GID1/HSL, and the O2 sides of steroid receptors and Nurr1. The others fail because SOGs are not
event-eligible (no A) or the specificity test has no derived-only clade.

Suggested wording: "recall 0/7: in 3–4 items the derived and ancestral proteins are in different orthogroups; in the
rest the SOG has no A members or no testable derived-only clade."

## C. Cross-OG extension

### C1 Method before outputs: PASS
- Commit order: `git log` shows 77db16c (17:31:22) → 0a9a242 (timestamp wording only) → 6f5bac5 code (17:37) →
  538bc2e → 461e082 validation outputs (17:46) → c79c4b3 → 57b7fd2 characterisation and CROSSOG.md (18:07).
- `git merge-base --is-ancestor 77db16c 461e082` → true.
- The earliest bulk output file, `/tmp/pe_big/partB_crossog/.../null_a_O1.parquet`, has mtime 17:39:21, after 77db16c.
- `git diff 77db16c HEAD -- docs/CROSSOG_METHOD.md` is a 1-line timestamp correction only.

Notes:
- `validation_recall.json` records git_commit fe67f2c, the REPORT commit made at 17:43, which descends from 77db16c.
  That is fine.
- The disclosed pre-run GR/MR alignment check looked at one validation pair before the run. The rule was not changed.

### C2 0/7 as specified: PASS
The code (`crossog.py` mode_validate) implements method §6:
- tested pairs = derived OG × ancestral OG
- specificity recovered = an (a) call on a link that contains the anchor SOG
- loss/gain recovered = a directional (b) call; the item counts as recovered if any expected type is

Raw outputs:

| run | tested links | call_a | min q_a | (b) rows | call_b | pre-rescue |
|---|---|---|---|---|---|---|
| O1 | 9 | 0 | 0.76 | 21 | 0 | 0 |
| O2 | 13 | 0 | 0.81 | 26 | 0 | 1 |

So 0/7 holds whatever the direction logic, and the lenient version is also 0. In-scope items = 5, matching CROSSOG.md.
The O2 `family_rule` is "None" for all items; per method §6 pairs are tested regardless, so it does not affect recall.

### C3 610 vs 1 links; (b) 10 vs 33: PASS-with-notes
Recomputed from the characterisation parquets:
- Links with J ≥ 0.2 and alntm ≥ 0.5: real **610**, null **1**. Tested: 606 vs 1.
- (a) calls 0 / 0. Real min p 0.0084, and 1.82 % of real links have p < 0.05.
- (b) recomputed from the counts (testable, alntm ≥ 0.5, coverage ≥ 0.8, P ≥ 3 species, no P/W, A ≥ 3 species):
  pre-rescue 20 / 35, post-rescue **10 / 33**, identical to the stored `call_b`.
- Unordered OG pairs with a (b) call: 8 / 27. Ordered pairs give 8 / 28. The text uses 27, the unordered count, which
  is fine.
- The null rows are the same projections as the real ones (same `sog_x` / `og_y`).

Note: "528 SOGs" in CROSSOG.md and the D32 text is the number of distinct `sog_x`. 1,051 distinct SOGs take part in at
least one link (548 distinct `sog_y`), so "link 528 of 14,284 eligible SOGs" understates linking.

## D. REPORT.md v0

`.venv/bin/python scripts/check_report_numbers.py REPORT.md` → `303 notable numbers, 6832 values in numbers.json,
0 unmatched; RESULT: PASS` (0.17 s). REPORT.md is unchanged since fe67f2c.

Spot-checks for meaning:

| # | claim (REPORT location) | verdict | evidence |
|---|---|---|---|
| 1 | §0.4 / §3.6: empirical curveball FDR 0.199 (O1), 0.267 (O2); pseudo-clade 0.080 / 0.093 | PASS | A1; case_list_v2.parquet.json carries the same values |
| 2 | §3.6 table: O2 gain 104 "called (pre-rescue)" | PASS-with-notes | O2 events and nulls run on the A-rescued states; numbers.json key `called_post_rescue` |
| 3 | §3.6: pass cascade 7 (4) / 1 (0) / 17 (13) / 6 (1); paralog flags 32/11/57/23 | PASS | A4; artifact_summary.json |
| 4 | §3.6: "2 of 50 re-called on O2 (0 of the 32 paralog-branch gains)" | PASS | o1_o2_gain_agreement.json, by_parent_kind duplication 0/32 |
| 5 | §3.6: "of 782 outgroup A members … 617 sit in a different O2 root group" | PASS-with-notes | A3: the metric is "not in modal O2 SOG"; the root-HOG recount is 734/923 (same 79 %) |
| 6 | §3.6: "Absent outgroups are mostly orthologs of another paralog copy" | PASS | A2: 88 % of failed forward hits are in the same root HOG, outside the derived clade (O2). The RBH numbers are still a TODO in the REPORT: suggest 0.197 overall, 0.165 gains, 0.260 losses, plus this 88 % |
| 7 | §3.6: specificity 222,489 pairs; 513 residue+shift; 0 BH; 3 null seeds | PASS | specificity/tau0.2/funnel_real_vs_null.csv |
| 8 | §1: "Thresholds were set to 10 % FDR on a simulated no-event null" | PASS | thresholds.json `event_fdr_target: 0.1`; the realised gain FDR is 0.051 at the prereg minima |
| 9 | §0.5 / §3.7: 0/7 "because they lie between paralogous orthogroups" | **FAIL** | B: false for D14/KAI2 and lysozyme (same OG in both), Nurr1 on O1, ER–AR/GR on O1, PfLDH–PfMDH in both |
| 10 | §3.7: "uricase and GULO are G (no gene) in human, never A" | **FAIL (minor)** | GULO on O2: human rows are P/U from paralogs LDHD/D2HGDH, not G. The "never A" part and the 5/5 pass hold |

## Requested fixes (priority)
1. Add a D-entry recording that D28's O2 PANTHER gate, and the O1 O2-root gate, were replaced by the DIAMOND RBH gate
   (commit 0cf2502). State that on O2 the PANTHER rule is never evaluable (104/104 unknown), and that O2 pass_artifacts
   would be 34, not 17, under D28 as written.
2. Reword the 0/7 explanation in REPORT §0.5, §3.7 and D31 (see B). Remove "MR/GR" as an example of a between-OG pair
   for O2.
3. Relabel or recompute the 617/782 figure (A3).
4. Change the REPORT §3.6 column header for O2 to post-A-rescue. Note that the O1 and O2 null FDRs sit on different
   rescue stages.
5. Disclose that 23 % of O2 gain outgroup-A supporters were never rescued (A5).
6. Minor: fix the null-seed provenance in `run_nulls.py` sidecars; fix "528 SOGs"; refine the GULO wording.
