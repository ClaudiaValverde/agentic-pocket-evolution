# Retrospective at G2 (fresh-process review; read PROMPT process rules, tree, LESSONS, DESIGN s0/s12, git log)

## What should we stop doing?
- Committing machine-regenerable bulk to git: ed79fdc adds ~264 `.parquet.json` sidecars + 10,563-line BULK_MANIFEST and 12 MB parquets in one commit; commit manifests/checksums only, keep bulk in /tmp/pe_big (D2) or an artifact dir.
- Committing parquets while home quota is ~50 GB (LESSONS 2): repo growth is an unforced quota risk.
- Writing thresholds as placeholders that later prove 8-13x wrong (false-A 0.05 vs 0.65, LESSONS 4); measure first, then simulate, or flag outputs "provisional" in sidecars automatically.
- Accepting prereg rules without a fallback: D9 (no k meets 5% target), D7 (thresholds must be re-derived), D8 (tau at grid edge, test-tied with 0.3) all became ad-hoc PI decisions.
- Setting tuning grids whose optimum sits on the edge (tau=0.2 is the grid edge); extend the grid before freezing.
- Accumulating deviations in the decision log (D7, D8, D9 "NEW conservative filter") rather than in DESIGN s12 as `check_prereg.py` requires; keep one deviation ledger.
- Duplicate/contradictory tree entries: M3 Part B engine appears twice ([x] and [~]); "Observation error rates" [~] while D9 already used its results; update status in place.
- Declaring the G2 slice "OK" on 5 hand-picked well-studied families (CDK2, ESR1/AR, LDHA, CA2) with one SOG each; this is a smoke test, not evidence of recall (test recall only 0.59).
- Starting parallel builders before the broadcast of resource rules (LESSONS 7: 3 agents hit quota, GPU 3 reserved mid-run); put rules in a shared brief template file.
- Running two design review rounds of ~2 h each for a 42-value config the reviewer could have drafted in round 1 (LESSONS 6); ask for the config at round 1.
- Leaving AF3 benchmark / co-folding controls [~] since M0 while they gate G4; time-box them.
- Status header still labeled "G1" after G2 reached (D8); refresh the header at every gate.
- Trusting same-agent summaries for pipeline numbers (F1 0.72/0.68, 66% missed sites): verifier checklist exists but no verification commit appears in the log yet.

## What worked well?
- Frozen snapshot manifest (D1) plus sidecar {snapshot, prereg sha, commit, seed}: all 22 commits are reproducible against 20261030T085136Z.
- Parallel low-regret work during review (LESSONS 3): orthology, annotations, sites landed by G1 and G2 hit at T+6 vs T+7 target.
- Error-rate measurement on known answers before trusting absences (v6->v4 re-detection 0.65, LIGYSIS 66% missed) reframed Part B honestly.
- Tune/test hash split by root HOG: tau reported with tune F1 0.72 vs test 0.68, tie with 0.3 disclosed rather than hidden.
- Conservative-direction deviations (D8 small_footprint, D9 k=1) all documented with direction of bias and "why".
- Controls reported as limits: AF3 places decoys in pockets (ipTM 0.72-0.86) so co-folding is demoted to relative evidence early.
- Detection-limited positives discovered at design time (catalytic triads not in BURROW pockets) and convergence repositioned to cofactor sites.
- Reviewer diversity (8 reports, round 1: PI/stats/data/prior work/structural/AI-failure) produced concrete estimand changes (F_analog vs convergence, tail-fitted null).
- Toy and 13.7M-SOG simulation validation of the Part B engine (null funnel exposed 2-species losses as noise).
- Small, frequent commits with tree updates after each decision; tree sufficient to resume.

## Proposed new LESSONS.md rules
1. Commit manifests and checksums, not regenerable bulk or per-file sidecars; one sidecar per directory. Evidence: ed79fdc added ~264 sidecars and 10.5k-line manifest.
2. Every tuned threshold needs the grid extended until the optimum is interior, or a pre-declared edge rule. Evidence: tau=0.2 at grid edge, test-tied with 0.3.
3. Every prereg threshold carries a `fallback:` key checked by `check_prereg.py`; no key, no freeze. Evidence: D9, D7.
4. Mark outputs computed with placeholder parameters `provisional=true` in the sidecar and block downstream consumption. Evidence: 5% placeholder vs 0.65 measured false-A.
5. A gate slice must include at least one negative (no-event, unrelated-family) case and a held-out family, not only famous ones. Evidence: G2 slice = 5 well-studied families, one SOG each.

## Existing rules to drop/merge
- Merge 2 and 7 into one "broadcast quota, storage and GPU rules via a shared brief template" rule (same root cause).
- Merge 3 and 6 into "parallelize design-independent work and ask reviewers for config drafts, not prose" (both save critical-path time).
- Drop 1 as a standalone rule (now embedded in DESIGN s0 and D1); keep only as an evidence line under the snapshot principle.
