# Final independent verification of REPORT.md

- Verifier: independent verifier agent (Claude Opus 5.5). I did not write REPORT.md or produce any result in it. Not a PI review.
- Date: Sat Oct 31 20:59:35 CET 2026 (`date`). Repo HEAD b83541b, branch `pocket-evolution`. Read-only on all project files.
- Inputs read: PROMPT.md, REPORT.md (471 lines), results/numbers.json, DESIGN.md §12 ledger (D7–D34, C-D7, C-D31),
  reviews/verify/*.md (verify_G4b, verify_D34_cofold and others), RESEARCH_TREE.md status and decision log, DATA.md,
  webapp/README.md, case_list_v2 sidecar, artifact_summary.json (O1_r3, O2_r1), o1_o2_gain_agreement.json,
  results/cofold/deep_dives/*, commit 21b5c9e (deep-dive criterion).
- The §5 "G5 second-snapshot comparison" item was ignored, as instructed.

## 1. Number check

`.venv/bin/python scripts/check_report_numbers.py REPORT.md` → `335 notable numbers, 7888 values in numbers.json, 0 unmatched`,
**RESULT: PASS** (exit 0).

Caveat: the checker matches values anywhere in numbers.json, including `partA_G5_20261031T020009Z`. I traced the Part A
headline values (0.819, 0.367, 0.935, 0.081, 0.021, 9,776,291, 2,627,453, 0.815, M1/M3 values) to their paths. All come from
`partA` / `partA_O2` / `partA_descent` with `snapshot_id 20261030T085136Z`, so they are snapshot-1 values (the 531f95f
restore holds). The co-folding numbers match `cofold.*` as exported after 464a0b9/7e1dc8b.

## 2. Meaning check (20 claims)

| # | line | claim | verdict |
|---|---|---|---|
| 1 | 16–18 | F_analog 0.819 [0.810, 0.826] O1 per family pair; 0.935 O2; per pair 0.367 on both | PASS (partA/partA_O2, snapshot 1) |
| 2 | 17, 130–131 | per pair is "granularity-robust" / "per pair does not [move]" | **OVERSTATED**: the test-half rows (l.174–175) give 0.285 (O1) and 0.136 (O2). Issue N2 |
| 3 | 19–20 | F_conv,func 0.081 [0.007, 0.178]; 0.021 [0, 0.062] | PASS (`partA.functional[0]`, `[3]`) |
| 4 | 20–22 | 3 of 2,627,453 directional passes; below toy false-pass 4.4e-5 | PASS |
| 5 | 23–25 | ortholog conjecture not supported on either orthology | PASS (O1 −0.0076 / −0.0110, O2 CI spans 0) |
| 6 | 26–28 | curveball FDR 0.199 O1, 0.267 O2 | PASS on values. Not comparable: O1 is pre-rescue (175 calls), O2 is post-rescue (verify_G4b A1 note 2). Issue N9 |
| 7 | 29–31 | specificity: 0 BH; population signal, FDR ≈ 0.15 (D34) | PASS on values. "several-fold" omits that, given the residue gate, the shift adds only 1.4–2.0× (verify_D34 A3). Issue N3 |
| 8 | 32–33 | gain failures "mainly on orthology (… another paralog subfamily of the same HOG, D33)" | **PARTLY WRONG for O1**: D33 restricts the same-HOG reading to O2. On O1 the dominant reason is a different O2 root HOG (734/923). Issue N1 |
| 9 | 34–37 | recall 0/7, mixed causes (C-D31) | PASS. Contradicted by §4 item 2 (B1) |
| 10 | 41–43 | 3 gate-passing cases, cognate ligand 12–34 Å away in "a different, top-ranked pocket" | PASS-with-notes: 12–34 Å holds for per-member medians and top models (samples 12–49 Å). The ligand pocket is rank 2–3 for MetE. Also, no cognate job "binds" by D23 (it is "placed"). Issue N4 |
| 11 | 43–45 | apo refolds: 3/4 testable gains, 1/4 testable false-loss controls | PASS (`cofold.apo_pool`). HOG0025964 is labelled "uninformative" in the table but counted as a testable non-reproduction in §0. Minor |
| 12 | 45–46, 327–331 | MR/GR resolved 3/3 seeds (+2.96/+2.68/+2.80) | PASS against the pre-committed criterion (21b5c9e). Caveats are missing (holo complexes are likely in AF3 training; all 4 cells pass D23). Issue N6 |
| 13 | 332–336 | LDH/MDH uninterpretable; independent rerun "reproduces the steroid numbers exactly and uses an LDHA holo gate that passes" | **INCOMPLETE**: in that rerun all 4 PfLDH/MDH test cells place the ligand 9.8–12.1 Å off-site (recall 0.17–0.33, d23_binds False). In the primary run LDH+PYR also fails D23. Issue N5 |
| 14 | 225–235 | Part B call counts and cascade 7(4)/1(0)/17(13)/6(1), paralog-branch 32/11/57/23 | PASS (artifact_summary.json) |
| 15 | 241–246 | 2 of 50 O1 gains re-called on O2, 0 of 32 paralog-branch; 79.5 % root-HOG recount | PASS |
| 16 | 247–249 | "Specificity: negative … O2 specificity is blocked: tier-3/L sites have no fingerprints" | **STALE/CONTRADICTS** D34 and §0.4. O2 was run (8eba86a; 743,349 pairs, 536 vs 66–84). Issue B2 |
| 17 | 252, 255, 120–121, 426–437 | dup-gain counts, false-loss pool sizes, Part B wall times, RBH fractions, robustness listed as TODO | **STALE**: all of these are in `numbers.json partB.*` (exported in 5e24c53). Issue B4 |
| 18 | 262–269 | per-category filter counts | PASS (case_list_v2 `counts`). "supported" is not defined in the text. Specificity row shows only the G4 top-5, not the D34 cascade candidates (200/55 O1, 230/50 O2) |
| 19 | 319–325 | false-loss controls: reappears 1 (0.60), not 3, 1 untestable | PASS on values. "so 'bad prediction' does not explain those absences" overstates. A fresh model from a related predictor does not prove a real loss; recall on true positives was only 0.44–1.00. Issue N7 |
| 20 | 160–162 | site-matching recall 0.591 / 0.633 | PASS-with-notes: this recall is conditional on both sites being detected. Unconditional recall is 0.31 (O2 test, `recall_unconditional`). Issue N8 |

## 3. Coverage against PROMPT.md

| requirement | status |
|---|---|
| Part A: three questions, with uncertainty | Yes (§0.1–0.3, §3.3–3.5; bootstrap CIs) |
| Part B: each category has cases or an explicit "none"; evidence per filter recorded | Yes for gain/loss/specificity/convergence (§3.7 table, all 0). Per-filter evidence for D34 specificity candidates is missing |
| False-loss negative control reported honestly | Mostly. §3.8 reports that the prediction is refuted for 3/5. §0.5 "negative controls produce no spurious losses" refers to the validation controls and is easy to misread. §4 does not carry the refutation (N7, N10) |
| Co-folding positive/negative controls, deep dives | Yes (calibration panel 24/72, gates, MR/GR, LDH/MDH). LDH rerun incomplete (N5) |
| Web app single command | Yes (l.419; webapp/README.md). Notes: the matrix view uses O1_r2 while results use O1_r3, and the tier-3 index lives in /tmp and must be rebuilt after a node restart (N11) |
| DATA.md with source/version/checksum/licence | Exists, but 4 rows have TBD version/checksum (AFDB v6 per-tar, PDB mirror date and checksum, AF3 cached MSAs, AF3 weights). REPORT never links DATA.md (N12) |
| Runtime with accuracy | Partial: table at l.105–118. Part B wall times are wrongly listed as TODO (B4) |
| Follow-ups, limitations, remaining ideas | Present (§4, §6). Limitations incomplete (B3) |

## 4. Required limitations

| topic | present in §4? |
|---|---|
| PI idle gap (no commits 18:07–22:56; PI idle 23:58 on 30 Oct → 20:39 on 31 Oct) and the missed 09:45 deadline | **MISSING** (not anywhere in REPORT) |
| Unblinded co-fold scoring | **MISSING from §4** (only §0 l.46 and §3.8 l.305) |
| Two-PI-session incident (reverted 48f713e/4d54c6b, false "PI review" attribution; O2 specificity run twice; deep-dive rerun from the other session) | **MISSING** (only an unexplained "parallel session's agent" at l.334) |
| G5 re-runnability failures (partB_O1, partB_O2 and the final Part A stage failed on snapshot-specific inputs; fix in b83541b) | **MISSING**. §7 l.410–412 presents the whole-pipeline rerun as working |
| Post hoc decisions | **Generic only** (l.375). Not named: D33 RBH rule chosen after O2 came back "unknown" (17 vs 34 O2 gains), D34 (b)/(c), D29 reframing, D25b literature f2, the duplication-gain track, cross-OG. Appendix A stops at D31 |

---

## Issues with proposed replacement text

### Blocking

**B1. l.342–343: stale cause for 0/7; contradicts C-D31, §0.5 and §3.7.**
Replace
`2. **Recall 0/7 on known events**, caused by the within-orthogroup scope (D31). The pipeline cannot see events between paralogous orthogroups, which is where most documented specificity shifts and losses happen.`
with
`2. **Recall 0/7 on known events**, with mixed causes (C-D31). MR/GR (O1), IDH and GID1/HSL lie between orthogroups, which the within-orthogroup pipeline cannot see. D14/KAI2, lysozyme/α-lactalbumin, Nurr1 and ER vs AR/GR (O1) and PfLDH/MDH are within one orthogroup but are missed for lack of confident absence calls, a pocket detected in both states, or a missing or non-significant specificity test. The cross-OG extension also recovers 0/7 (D32).`

**B2. Stale O2 specificity status; contradicts D34 and §0.4.**
- l.247–249: replace
`- **Specificity: negative.** 222,489 branch pairs tested on O1; 513 pass the residue gate and the shift; 0 pass BH, as in all 3 null seeds (…funnel_real_vs_null.csv). O2 specificity is blocked: tier-3/L sites have no fingerprints.`
with
`- **Specificity: 0 preregistered BH calls, which is uninformative by construction** (D34): the smallest attainable empirical p (3.5e-5 O1, 1.6e-5 O2) is above the single-call BH threshold (4.7e-6, 2.6e-6). O1: 222,489 branch pairs, 513 pass residue gate + shift vs 71–88 under the pseudo-clade null (3 seeds; D24 asks for ≥ 20). O2 (fingerprints read from atlas rows, 8eba86a): 743,349 pairs, 536 vs 66–84. Empirical FDR ≈ 0.15 at the funnel and 0.17–0.19 after the cascade (200 vs 30–44 on O1, 230 vs 32–45 on O2; 4 cascade checks are not evaluable on these pairs). Given the residue gate, the fingerprint shift adds only 1.4–2.0× (verify_D34_cofold A3). This is a population-level signal, post hoc, and no individual pair is supported.`
- l.349: replace `5. **Specificity: 0** BH-significant branch pairs; O2 specificity not run (no fingerprints for tier-3/L sites).` with `5. **Specificity: 0** BH-significant branch pairs on O1 and O2. BH cannot fire at the null's resolution (3 null seeds). The population-level excess (D34, post hoc) has no individually supported pair.`
- l.395: replace `3. **Fingerprints for tier-3/L sites**, so specificity can run on O2 with denser species sampling.` with `3. **≥ 20 pseudo-clade null seeds (D24) and a preregistered tail model for specificity**, so that individual specificity calls become attainable.`

**B3. §4 lacks the required process limitations.** Insert after l.376:
```
15. **Schedule and oversight.** No commits 18:07–22:56 on 30 Oct (co-folding stalled ~18:40–22:58), and the PI was idle
    23:58 (30 Oct) → 20:39 (31 Oct). The 09:45 delivery deadline was missed, and this report was finalised afterwards.
16. **Co-fold case scoring is unblinded.** A post-processing bug re-minted the blinding salt, and the state mapping was
    recovered by sequence before the final scoring pass. The folds ran blind. The deep dives stayed blinded.
17. **Two PI sessions ran in parallel.** Agent commits 48f713e/4d54c6b wrote a "PI review" with false attribution and were
    reverted. O2 specificity was run twice (8eba86a from the other session; the D34 rerun is bit-identical, which shows
    reproducibility, not independent replication). The LDHA-gated deep-dive rerun also comes from the other session.
18. **Re-runnability on a new snapshot failed at G5.** On snapshot 20261031T020009Z, partB_O1, partB_O2 and the final
    Part A stage failed on snapshot-specific inputs that were not carried over (rescue_k.json, the O1 events path,
    TimeTree pairwise). Fixes are in b83541b. Until a clean rerun passes, the §7 single-command rerun is not demonstrated.
19. **Post hoc choices.** On O2, the RBH orthology rule was chosen after PANTHER/O2-root came back "unknown" for all
    gains (D33: 17 pass instead of 34). Also post hoc: the D34 specificity population signal and tail fits, the D29
    per-pair "robust" framing, the D25b literature f2 line, duplication gains, and the cross-OG extension (D31/D32). None
    of these enters the convincing headline.
20. **Unrescued absences.** 105,278 O2 A calls with no filtered pockets were never rescued. They make up 23 % of the O2
    gain outgroup-A supporters (verify_G4b A5).
```

**B4. Stale status and stale TODOs.**
- l.3–5: replace `**Status: draft v0, written 2026-10-30 17:34 CET (`date`). Development subset, preliminary.**` and `Sections marked **PENDING (filled at G6)** are not written yet.` with `**Status: final, written <date output> (after the 09:45 deadline; see §4.15). Development subset, preliminary.**` and `Only the G5 second-snapshot comparison (§5) is pending.`
- l.120–121, l.252, l.255 and the TODO list (l.430–436): these values exist in `numbers.json partB` and must be cited.
  - Part B wall times: `partB.wall_times`.
  - RBH 0.095 O1 / 0.197 O2 and Pfam-arch match 0.731 / 0.835: `partB.rbh_orthology`.
  - Duplication gains: 8/32 O1, 18/55 O2; background ≥ 3: 0/2,597 and 4/5,276.
  - False-loss pools: bad-model 319, rescue-instability 48.
  - O1 gain robustness: 1:1 0.98, **3:1 0.34**, τ = 0.3 0.84.
  - O2 empirical null FDR: `partB.O2_funnel_and_null`.

  For example, l.252 `(counts: TODO numbers)` → `(8 of 32 on O1, 18 of 55 on O2; unselected background with ≥ 3 species: 0 of 2,597 and 4 of 5,276; numbers.json partB.duplication_gains)`. l.255 `(pool sizes: TODO numbers)` → `(pools on O1_r3: 319 bad-model, 48 rescue-instability; numbers.json partB.false_loss_controls)`. Add to §3.6: `- **Robustness of O1 gains:** 98 % are also called at cost 1:1 and 84 % at τ = 0.3, but only **34 % at 3:1** (numbers.json partB.O1_primary_gains).` The 3:1 fragility is material and is currently hidden behind a TODO.

**B5. Appendix A (l.439–471) stops at D31.** Append:
```
| D32 | Cross-OG extension 0/7, gate failed, no genome-wide run | none; negative result |
| D33 | O2 artifact orthology check = DIAMOND RBH (post hoc, after O2/PANTHER were 'unknown'); 17 vs 34 O2 gains under D28 as written; RBH failures = other paralog subfamily of the same HOG | fewer cases |
| C-D31 | 0/7 causes are mixed, not only between-OG | wording; 0/7 stands |
| D34 | Specificity: BH unattainable; population-level real-vs-null excess and tail-fit sensitivity reported post hoc | can only add; outside headline |
```

### Non-blocking

**N1. l.32–33.** Replace `mainly on orthology (comparisons against another paralog subfamily of the same HOG, D33) and present-side fragments.` with `mainly on orthology (on O1, outgroups that O2 places in a different root HOG, 79.5 %; on O2, comparisons against another paralog subfamily of the same HOG, D33) and present-side fragments.`

**N2. l.17–18 and l.130–131.** Replace `The **granularity-robust** per-pair estimate is` with `The per-pair estimate, which agrees between O1 and O2 on all pairs (post hoc framing, D29), is`. Replace `it moves with orthology granularity; per pair does not.` with `it moves with orthology granularity; per pair agrees on all pairs (0.367 vs 0.367) but not on the test halves (0.285 vs 0.136, different splits).`

**N3. l.30–31.** Replace `real sister pairs exceed the pseudo-clade null several-fold (empirical FDR ≈ 0.15; D34)` with `real sister pairs exceed the pseudo-clade null 6.5–7× (about 1.4–2× attributable to the fingerprint shift once the residue gate is passed; empirical FDR ≈ 0.15; D34, post hoc, 3 null seeds)`.

**N4. l.42–43.** Replace `the cognate ligand binds 12–34 Å away from the gained footprint (a different, top-ranked pocket)` with `AF3 places the cognate ligand 12–34 Å (per-member medians) from the gained footprint, in a different pocket (BURROW rank 1–3 in 17 of 18 models), and no cognate job passes D23`.

**N5. l.334–336.** After `reproduces the steroid numbers exactly`, add: `(same inputs, so reproducibility rather than replication); its LDHA gate passes, but in that rerun all four PfLDH/MDH cells place the ligand 9.8–12.1 Å from the site (recall ≤ 0.33, none passes D23), and in the primary run LDH+pyruvate also fails D23, so neither run supports an LDH/MDH resolution`.

**N6. l.330–331.** Replace `**AF3 resolves a known specificity difference**, but the effect (~3 pLDDT) is small next to the cognate-vs-decoy gap.` with `This meets the pre-committed "resolved" criterion, but the effect (~3 pLDDT; PAE interaction 0.03) is small. All four cells, including GR + aldosterone, pass D23, and the receptor–steroid complexes are likely in AF3's training data, so this is weak evidence of ligand-level discrimination.`

**N7. l.323–325.** Replace `so "bad prediction" does not explain those absences, and the control label is unsupported for them.` with `so a fresh model does not restore the pocket. Whether these are real losses discarded by the quality filter or errors shared across AF2/AF3 is unresolved, and the control label is unsupported for them.` Add to §4: `21. **False-loss controls:** for 3 of 5 bad-model controls, a fresh ligand-free AF3 model also lacks the pocket, so the A-quality filter may be discarding some real losses as well as artifacts.`

**N8. l.158.** Change the column header `recall / precision` to `recall (both sites detected) / precision`, and add below the table: `Unconditional recall (any benchmark pair) is 0.31 on the O2 test split (tau_benchmark_o2_split.csv recall_unconditional).`

**N9. l.27 and l.346.** After `0.267 on O2`, add `(O1 on 175 pre-rescue calls, O2 on rescued states; not directly comparable)`.

**N10. l.34 and l.37–38.** Fix the typo `events, . The causes` → `events. The causes`. Replace `The negative controls produce no spurious losses.` with `The validation negative controls (uricase, GULO, M. leprae truncations) produce no spurious losses; the separate false-loss co-fold controls are in item 6.`

**N11. l.102–103.** Add `The state matrix shows O1_r2. Part B results use O1_r3 (same 260,991 SOGs per funnel.json; state rules may differ). The tier-3 protein index lives in /tmp and must be rebuilt (python -m webapp build-index --tiers 3) after a node restart.`

**N12. DATA.md / REPORT.** REPORT never links DATA.md. Add to §1 `Data` (l.54): `Dataset sources, versions, checksums and licences: [DATA.md](DATA.md).` In DATA.md, 4 rows still have TBD (AFDB v6 checksum, PDB mirror date and checksum, AF3 MSA checksum, AF3 weights version and checksum). Fill these, or state that they are unknown for read-only mounts.

**N13. l.378.** The heading `## 5. Pending sections (filled at G6)` → `## 5. Exploratory extension and pending comparison`, because the cross-OG part is filled.
