# EVIDENCE: M4 evidence filters and case selection

Code: `src/pocketevo/evidence.py` holds the rules and `scripts/cases/build_cases.py` does the I/O. Tests are in
`tests/test_evidence.py` (toy data with known answers). Outputs: `results/cases/<snapshot>/case_list_<v>[_<orth>].parquet`
with a `.json` sidecar (prereg sha256, git commit, seed, input sha256s, rules, counts, runtime), and
`cases/<case_id>/case.yaml` (schema case-v1, plus `evidence[].filter` and a pending `cofold` placeholder, see webapp/SCHEMA.md §4)
for the selected cases and the false-loss negative controls.

```bash
.venv/bin/python scripts/cases/build_cases.py --orthology O1 --version v1 --verify-ids
# O2: --orthology O2 --sog-release O2_r1 --candidates <O2 candidates> --events-base <O2 primary events dir>
#     --hog-trees results/partB/hog_trees/O2/hog_trees.parquet --false-loss <O2 controls> --tau-tuning <O2 tuning set>
```
The builder reads the candidate table again on every run, so it works unchanged when the table gains specificity rows or
artifact columns (`pass_artifacts`, `art_*`).

## Rules
Values from `config/prereg.yaml` (`part_b.evidence_filters`, `case_selection`) are read at run time: f1 needs at least 3
species, the BioLiP overlap must be at least 0.5, the experimental ECO codes are 0000269 and 0007744, pChEMBL must be at
least 6, and selection takes top-k = 5 per category plus 1 seeded audit case. The M4 constants were written before any
result was seen. They are the module constants in `evidence.py`.

- **f1 persistence** passes when all of these hold:
  - the change state is seen in at least 3 derived species;
  - at least 80 % of the informative derived members carry the change;
  - the event is single and tree-consistent (`f1_single_event`);
  - at least 2 outgroup species carry the ancestral state.

  The result also records whether the change forms a clean clade or is scattered. This uses the purity over all SOG
  members under the species-tree LCA.
- **f2 experimental** maps the footprint MSA columns to each member's residues. This is only done when the model
  sequence CRC64 equals the UniProt 2026_03 sequence. Site-level hits (an overlap of at least 0.5 of the annotated site)
  come from three sources:
  - BioLiP2 holo sites;
  - UniProt Binding or Active sites with experimental ECO codes only;
  - M-CSA *reference* residues. Homologue rows are reported but not counted.

  ChEMBL counts only at the protein level, so it never passes f2 (`about_the_change_required`). The two sides are
  compared:
  - A gain passes when the derived members that have the pocket carry a site-level hit and the outgroup members carry none.
  - A loss passes when the outgroup or ancestral members that have the pocket carry a hit and the derived members carry none.
  - A site-level hit on the absent side contradicts the call, and f2 then fails.
  - A missing annotation on the absent side is recorded as "weak".
  - A specificity change passes when both sides carry small-molecule ligand tokens at the footprint and the two token
    sets are disjoint.
  - A convergence pair passes when both sites share a site-level ligand token.

  Identifiers are checked against RCSB and UniProt (`--verify-ids`). An unresolvable ID (HTTP 404) fails f2. A network
  error leaves the ID unverified.
- **Tau-tuning disjointness** (`tune_test.annotation_disjointness`): an annotation used to tune τ is dropped from f2 for
  the SOG that contains its matched site. This is the literal prereg rule. Two other variants are reported in the case
  list but never used for labels: OG-wide exclusion (`f2_ogwide_exclusion`) and no exclusion (`f2_diag_no_exclusion`).
- **f3 computational** passes when at least one *independent* line passes:
  - **conservation:** at least 2 footprint columns are clade-specific, and the footprint rate is above the rate of the
    non-footprint background columns. A column is clade-specific when the derived majority residue is carried by at
    least 80 % of derived members and by less than 50 % of the outgroup.
  - **rescue_consistency:** at least 80 % of the present-side members are found in at least 3 of the 6 BURROW settings.
  - **source_consistency:** the site is re-detected in the v4 model (Jaccard ≥ 0.3) for at least 80 % of the
    present-side members that are in the v6/v4 re-detection sample.
  - **cofolding:** pending. The co-folding agent fills it under D23.

  For specificity, conservation is not independent of the call, so it is excluded there. For convergence, only
  source_consistency and cofolding count.
- **Labels:**
  - `filter_label` is convincing when all three filters pass, supported when f1 passes plus one of f2/f3, and rejected otherwise.
  - The final `label` is set to rejected when the artifact cascade fails. It is set to pending while `pass_artifacts`
    has not yet been written.
- **Score:** S = `score_S_partial` (persistence × min state confidence, from Part B) × evidence_count. The evidence count
  is the number of site-level kinds that pass f2 plus the number of passing independent f3 lines.
- **Selection:** top-5 by S among cases labelled `convincing`, plus 1 seeded random audit case from the rest of the
  convincing cases, per category. The bar is never relaxed, so a category with no convincing case gets no selected case.
  The false-loss controls with `negative_control_example` are always written as category `control`.
  Losses always carry `not FDR-calibrated (D21)`.

## Known weaknesses (open)
- rescue_consistency is close to tautological. Post-rescue P members already have at least 1 hit out of 6, because P at
  0/6 becomes W. It is weak computational support.
- Convergence sites from Part A have no SOG footprint in the partner OG, so conservation is not available for them.
- The artifact cascade comes from the Part B agent (`pass_artifacts`). Until it is written, every label stays pending.

## D25 changes (2026-10-30 15:20–15:35 CET)
- **(a) Bug fix in the disjointness rule.** `tau_tuning_annotations.parquet` lists *every* annotation of the tune-half
  OGs. The builder excluded the whole list, which over-excludes. It now excludes only annotations that entered a
  tune-half pair of `tau_benchmark_pairs.parquet` (`used_in_tuning`). That leaves 488 of the 827 listed rows for the
  candidate OGs. Regression test: `test_disjointness_excludes_only_annotations_used_in_tune_pairs`.
  - Q15257 (HOG0028478, tune half) appears in no pair, so it is no longer excluded.
  - P43004 CLR (HOG0035363, tune half) and the O43143/P53131 ADP rows (HOG0034468, tune half) are in tune pairs, so
    their exclusion stands.
- **Artifact check `orthology_side_family`** (DESIGN §5.4 "PANTHER disagreement"). It fails a gain or loss when the
  present-side and absent-side (A) members share no PANTHER family (subfamilies collapsed; Pfam if PANTHER is missing).
  - Why it was added: the literature review found unrelated families on the two sides (SLC1 transporters vs KAS,
    PTPA vs DEAD-box helicases, COX11 vs RFC).
  - It is post hoc and conservative (fewer cases).
  - It fails 28 of 51 gains and 6 of 18 losses on O1.
- **Manual inputs** in `cases/<id>/manual.yaml`:
  - `literature_evidence`: items with a PMID and a verbatim quote. Only items marked `about_change` and `review_ok`
    count, and only for `f2_with_literature` / `filter_label_with_literature` (D25b, post hoc).
  - `artifact_checks`: manual checks that the case adds to its own artifact checks.
- **`--pool supported`** (D25c, case_list_v2) applies the prereg top-5 + 1 audit rule to filter_label ∈ {convincing,
  supported} with artifacts ≠ fail. The headline stays convincing_only: convincing and artifacts not failed.

## Validation track and options (2026-10-30 evening)
- `scripts/cases/validation_track.py` (controls listed in `scripts/cases/validation_controls.yaml`) writes
  `results/cases/<snap>/validation_table.{parquet,md}` and `cases/val_<id>/case.yaml` (category validation_known).
  It runs the annotation-first track of DESIGN §5.3 on the O1/O2 outputs unchanged.
- `build_cases.py --snapshot <id>` rewrites every default path to that snapshot, and the defaults are the current
  O1_r3 inputs. The driver also takes it: `SNAPSHOT=<id> scripts/cases/run_case_list_v2.sh`.
  The side-family check is report-only by default (D30).
