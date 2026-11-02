# Part B candidates — O2 (O2_r1 + rescue overlay), snapshot 20261030T085136Z (development subset, preliminary)

Generated 2026-10-30 16:44 UTC by `scripts/partB/write_o2_report.py`; prereg sha256 `52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83`; git `fe67f2c7270af9cc942c23de86049ca81ff53b28+dirty`. Same code, cuts (O1 measured-noise thresholds) and nulls as O1_r3. Not yet verified.

## States

O2_r1 states τ 0.2 (D8 flags excluded, D11 applied by the SOG stage) + rescue overlays ['results/sogs/20261030T085136Z/O2_r1_rescue/rescue_A.parquet', 'results/sogs/20261030T085136Z/O2_r1_rescue2/rescue_A.parquet'] (k = 1): A in the rescue pool → W 50309 / A 50678; `rescue:pending` A outside the pool → U (341796); A without atlas-filtered pockets kept as A (105278, not re-described). P calls are not rescued on O2 (overlay is A-only).

## Funnel (real)

| input | cost | type | n_sogs | clade_size_ok | unambiguous | species | frac | outgroup | candidate | called | called_f1_single_event | chars_with_called |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| real | 2:1 | gain | 10185 | 654 | 399 | 153 | 104 | 104 | 104 | 104 | 80 | 104 |
| real | 2:1 | loss | 10185 | 7729 | 3695 | 3643 | 3562 | 2916 | 2916 | 66 | 21 | 62 |
| real | 1:1 | gain | 10185 | 1618 | 866 | 382 | 323 | 323 | 323 | 323 | 218 | 316 |
| real | 1:1 | loss | 10185 | 5094 | 1594 | 1576 | 1344 | 1274 | 1274 | 13 | 9 | 13 |
| real | 3:1 | gain | 10185 | 312 | 197 | 80 | 50 | 50 | 50 | 50 | 40 | 50 |
| real | 3:1 | loss | 10185 | 9272 | 5172 | 5110 | 5094 | 3978 | 3978 | 110 | 27 | 104 |

## Real vs null at the same cut (≥ 20 seeds; empirical null FDR = mean null calls / real calls)

| type | null | real_called | null_called_mean | null_called_max | n_null_seeds | empirical_null_FDR | simulated_FDR |
|---|---|---|---|---|---|---|---|
| gain | null_curveball | 104 | 27.800 | 45 | 20 | 0.267 | 0.051 |
| gain | null_pseudo_clade | 104 | 9.650 | 14 | 20 | 0.093 | 0.051 |
| loss | null_curveball | 66 | 102.000 | 122 | 20 | 1.000 | – |
| loss | null_pseudo_clade | 66 | 44.500 | 56 | 20 | 0.674 | – |

## Calls (2:1)

**104 gains (FDR-calibrated cut) + 66 losses (not FDR-calibrated, D21).** Yield preflight: ≥ 5 → secondary stratum not triggered.

| type | label | og_id | sog_id | level | clade_lca_species_tree | n_species_P | n_species_A | out_A_lineages | f1_single_event | score_S_partial |
|---|---|---|---|---|---|---|---|---|---|---|
| loss | not FDR-calibrated (audit) | HOG0047926 | HOG0047926:g7 | Opisthokonta | Opisthokonta | 1 | 38 | 5 | False | 38.000 |
| gain | FDR-calibrated candidate | HOG0054257 | HOG0054257:g2 | Opisthokonta | Opisthokonta | 34 | 2 | 6 | False | 34.000 |
| gain | FDR-calibrated candidate | HOG0049013 | HOG0049013:g2 | Eumetazoa | Eumetazoa | 25 | 0 | 5 | True | 25.000 |
| loss | not FDR-calibrated (audit) | HOG0052701 | HOG0052701:g14 | Eumetazoa | Eumetazoa | 5 | 28 | 5 | False | 14.000 |
| gain | FDR-calibrated candidate | HOG0056159 | HOG0056159:g12 | Fungi | Fungi | 14 | 0 | 7 | True | 14.000 |
| gain | FDR-calibrated candidate | HOG0054963 | HOG0054963:g20 | Eumetazoa | Eumetazoa | 14 | 0 | 22 | True | 14.000 |
| gain | FDR-calibrated candidate | HOG0055724 | HOG0055724:g2 | Bilateria | Bilateria | 21 | 2 | 3 | True | 10.500 |
| gain | FDR-calibrated candidate | HOG0047083 | HOG0047083:g6 | Eumetazoa | Eumetazoa | 10 | 0 | 5 | True | 10.000 |
| gain | FDR-calibrated candidate | HOG0008302 | HOG0008302:g0 | Pseudomonadati | Pseudomonadati | 10 | 0 | 5 | True | 10.000 |
| gain | FDR-calibrated candidate | HOG0054166 | HOG0054166:g12 | Embryophyta | Embryophyta | 8 | 0 | 3 | True | 8.000 |
| gain | FDR-calibrated candidate | HOG0054175 | HOG0054175:g6 | Fungi | Fungi | 16 | 6 | 4 | False | 8.000 |
| gain | FDR-calibrated candidate | HOG0057094 | HOG0057094:g6 | Thermoproteati | Thermoproteati | 8 | 0 | 3 | True | 8.000 |
| gain | FDR-calibrated candidate | HOG0056321 | HOG0056321:g9 | Embryophyta | Embryophyta | 8 | 0 | 4 | True | 8.000 |
| loss | not FDR-calibrated (audit) | HOG0051218 | HOG0051218:g9 | Mesangiospermae | Mesangiospermae | 1 | 7 | 4 | False | 7.000 |
| loss | not FDR-calibrated (audit) | HOG0056156 | HOG0056156:g15 | FCB group | FCB group | 0 | 7 | 13 | False | 7.000 |
| gain | FDR-calibrated candidate | HOG0050209 | HOG0050209:g8 | Rhabditida | Rhabditida | 7 | 0 | 6 | True | 7.000 |
| gain | FDR-calibrated candidate | HOG0053806 | HOG0053806:g5 | leotiomyceta | leotiomyceta | 7 | 1 | 3 | True | 7.000 |
| gain | FDR-calibrated candidate | HOG0034532 | HOG0034532:g4 | Amniota | Amniota | 7 | 1 | 3 | True | 7.000 |
| loss | not FDR-calibrated (audit) | HOG0043820 | HOG0043820:g10 | Mesangiospermae | Mesangiospermae | 0 | 7 | 4 | False | 7.000 |
| gain | FDR-calibrated candidate | HOG0056094 | HOG0056094:g11 | Pseudomonadati | Pseudomonadati | 7 | 0 | 4 | True | 7.000 |

## O1 → O2 agreement for gains (misassigned-orthology check)

50 O1_r3 primary gains: re-called on O2 2, any O2 MPR gain event 2, O1 derived members in one O2 SOG 50. O1 gains rest on outgroup absences; the share of those absent proteins that O2 places in a different root HOG than the derived clade measures orthology disagreement between O1 and O2.

Of the 782 outgroup A members behind the O1 gains, 617 are in a different O2 root HOG than the derived clade (O2 does not treat them as members of the same family).

```json
{
 "O1_state_vs_raw_O2_state_in_modal_O2_SOG": {
  "A": {
   "A": 37,
   "P": 34,
   "U": 94,
   "not_in_O2_SOG_family": 617
  },
  "P": {
   "A": 0,
   "P": 411,
   "U": 26,
   "not_in_O2_SOG_family": 17
  },
  "W": {
   "A": 8,
   "P": 9,
   "U": 15,
   "not_in_O2_SOG_family": 20
  }
 },
 "O1_A_that_are_U_on_O2_by_reason": {
  "footprint_residues<=10": 48,
  "partial_site_overlap": 25,
  "plddt<90": 9,
  "length_ratio": 8,
  "aligned<0.8": 3,
  "terminal_truncation": 1
 },
 "O1_A_members_in_another_O2_root_HOG": 617,
 "O1_A_members_total": 782
}
```

## Not done on O2
- Specificity: tier-3/L sites were built without fingerprints (`--no-fingerprints`).
- Artifact cascade and P rescue on O2 members; τ 0.3 sensitivity on O2.

## Artifact cascade (run_artifacts.py; one column per check, nothing dropped)

Checks (fail -> pass_artifacts false): a_footprint_plddt, a_disorder, a_pae_ted, a_fragment, a_fragment_present, a_orthology_o2, a_orthology_panther, a_instability, a_species_tree_lineages, orthology_split. Flags (reported, not fails): flag_paralog_branch, flag_present_plddt_lt90. PAE was not fetched: the TED multi-domain rule turns such A calls into U. `na` = no data for that check (listed in `artifact_na`).

| type | n | pass_artifacts | pass_strict | paralog_branch | O2 paralog split (T/F/na) | present pLDDT<90 | fail counts |
|---|---|---|---|---|---|---|---|
| gain | 104 | 17 | 13 | 57 | 57/0/0 | 63 | {"a_disorder": 2, "a_fragment": 7, "a_fragment_present": 45, "a_orthology_panther": 28, "a_species_tree_lineages": 13, "orthology_split": 80} |
| loss | 66 | 6 | 1 | 23 | 23/0/0 | 41 | {"a_disorder": 2, "a_pae_ted": 1, "a_fragment": 7, "a_fragment_present": 52, "a_orthology_panther": 4, "a_species_tree_lineages": 7, "orthology_split": 44} |

Primary calls only (τ 0.2, 2:1); specificity = BH-called pairs (or the exploratory top 5). Per-candidate columns are in `candidates.parquet` / `specificity_candidates.parquet`.

