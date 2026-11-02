# VERIFY_CHECKLIST — independent verifier
Rules: fresh context; committed artifacts only (code, manifests, outputs, config/prereg.yaml), never the
producer's summary. Record each item PASS / FAIL / N-A with the file+line or command used, in
`results/verify/<stage>_<run_id>.md`. Only the verifier marks a RESEARCH_TREE node `[x]` or `[!]`.

## 1. Provenance (all stages)
1. `sha256sum config/prereg.yaml` equals the value frozen in RESEARCH_TREE.md; else a dated deviation entry
   exists in DESIGN §12 with reason and before/after results.
2. Every output has a sidecar with `snapshot_id, prereg_sha256, git_commit, seed`; all equal across the stage.
3. `scripts/check_prereg.py` and `scripts/check_report_numbers.py` exit 0. `numbers.json` mixes no snapshots.
4. The stage manifest reports `n_in`, `n_out` and `n_dropped` by named reason, and the counts add up.

## 2. Re-run and hash
5. Re-run the stage on a declared small slice (one root HOG, one atlas unit, one toy world) from a clean
   checkout with `scripts/run_all.sh --snapshot <id> --slice toy`; diff outputs against committed ones.
6. Byte-identical for deterministic steps; for stochastic steps, identical under the prereg seeds.
7. Runtime recorded and within an order of magnitude of DESIGN §9.

## 3. Controls ran and passed
8. Toy worlds ran BEFORE real data: planted gains/losses (recall reported), no-event world (0 events),
   no-convergence world (delta null centred on 0), planted fragments/low-pLDDT (all U).
9. Null-input Part B funnel (curveball P/A, random pseudo-clades, same-protein fragments) is reported beside
   the real funnel, same code, same thresholds.
10. Negative controls: uricase/GULO give **G**, never A; M. leprae A-call rate reported before and after filters.
11. Co-folding: family holo positive control reproduces (ligand RMSD <= prereg), decoy and site-ablation
    controls behave. Failed gate => result labelled "uninterpretable", not dropped.
12. Detectability preflight present for every benchmark/control site; convergence recall is conditional on it.

## 4. No tuning leakage
13. tau, MSA method, rescue k and the absence pLDDT cut were chosen on the TUNE half only; headline numbers
    come from the TEST half; the full sensitivity grid is reported, not just the chosen cell.
14. Control-list families are in the test half; annotations used for tuning are not reused in evidence filter 2.

## 5. Spot checks by hand
15. Pick 10 random rows (seeded, record the seed) from the stage output; re-derive each by hand from inputs.
16. Pick 3 random A calls used by a candidate: confirm pLDDT, non-gap fraction, length ratio, SP/TM, interface
    and rescue all recorded and passing.
17. Resolve 5 cited identifiers (DOI/PMID/PDB/ChEMBL) against their API; a cached quote supports the claim.
18. Confirm every REPORT section carries the snapshot id and the "development subset, preliminary" statement,
    cases are never shown without their funnel denominator, and planned-but-not-completed items are listed.
