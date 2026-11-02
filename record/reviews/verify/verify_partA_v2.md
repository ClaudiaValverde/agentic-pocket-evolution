# Independent verification: Part A v2 and Part A descent (M1/M2)

Verifier: independent agent (did not produce these results). Written 2026-10-30 14:33 CET (`date`).
Repo branch `pocket-evolution`, HEAD `92a11d7`. I edited no result, metric, control or test files.
Check scripts are in `/tmp/verify_partA/` (scratch, not committed). Each one reads only the stored intermediate tables and
does not call the producers' stage scripts:
- Part A bulk: `/tmp/pe_big/partA/20261030T085136Z/` (`cls_cos/pairs_cls.parquet`, `sample_rest.parquet`,
  `functional_site_labels.parquet`, `hidden_homology_*.parquet`, `partA.log`);
- descent bulk: `/tmp/pe_big/partA_descent/20261030T085136Z/O1/pairs.parquet`;
- raw tables: `data/annotations/*`, `data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz`,
  `data/timetree/*`.

**Total wall time of my checks: about 6 min of compute, 9 min including setup (17 min elapsed, 14:16 to 14:33).**

## Verdicts

| # | claim / check | verdict |
|---|---|---|
| 1 | Headline F (3 weightings), r_c per class, n_sig, q0 | **PASS** (exact match) |
| 2 | Hidden-homology ε and corrected F_analog | **PASS** |
| 3 | numbers.json equals the stage numbers files | **PASS-with-notes** (F_conv,func keys lost) |
| 4 | Homology classes on a 30-pair audit | **PASS** (30/30) |
| 5 | N class not contaminated by obvious homologs | **PASS-with-notes** |
| 6 | F_conv,func (n = 192, label subset, matched control) | **PASS-with-notes** |
| 7 | Directional test F_conv,dir (null, max-p, 3 passes) | **PASS-with-notes** |
| 8 | D17–D19 decided before the results they set | **PASS-with-notes** (D18 partly post hoc) |
| 9 | Tune/test split and τ sensitivity reporting | **FAIL (non-blocking)**: reporting gap |
| 10 | Descent M1/M2 numbers | **PASS** |
| 11 | Descent: ortholog vs paralog confounded by time? D15 symmetry | **PASS-with-notes** |
| 12 | Toy tests | **PASS in the producers' env; FAIL with the requested `uv run` command** (torch missing) |

No blocking issue for the claims as worded: F_analog is an upper bound, F_conv,func is modest, and F_conv,dir is about 0.
Two items should be fixed before a write-up (items 3 and 9), and the D18 wording should be corrected (item 8).

---

### 1. Headline estimands: PASS
`uv run --project . python /tmp/verify_partA/recompute_F.py` (6 s)

I recomputed the estimands from `pairs_cls.parquet` with my own code; nothing is imported from `pocketevo`:
- per-family-pair, per-site and per-pair weights;
- K_c and n_c (exact totals × SOG split; rest × sample class shares);
- q0 = 1.72 · 0.01 / N.

Results:
- n_sig = 9,776,291, and q0 = 5.3388e-9.
- All 12 F values match `numbers_partA.json` to ≤ 1e-18: F_anc, F_offsite, F_H4 and F_analog under per family pair,
  per site and per pair. Headline: F_analog = 0.8186 per family pair and 0.3671 per pair.
- All 10 r_adj match to 4 significant figures. Examples: H1s 1.0016e-2 and N 7.0038e-7. The ratios in the text
  check out: 14,300×, 79× and 182×.
- The totals partition N(N−1)/2 minus same-protein pairs exactly.
- In the primary set (E < 0.01), every search weight is 1. The secondary 1.2 % search sampling only affects the
  E < 1 sensitivity row.

Note (not an error): under the per-family-pair weighting, all the ancestral mass collapses into a few within-group cells.
If I restrict to pairs where both proteins are in HOG-derived groups (dropping singleton and proxy_sf pseudo-groups,
which carry 28 % of the N weight), the raw per-family-pair F_analog rises to 0.95
(`/tmp/verify_partA/famw.py`). The headline therefore depends strongly on the weighting; this is the preregistered
choice, but it should be stated next to the number.

### 2. Hidden homology: PASS
`/tmp/verify_partA/hh.py`
- All 200 sampled pairs are in the significant N set.
- max(qTM, tTM) ≥ 0.5 in 1/200 significant N pairs and 2/200 random N pairs; min ≥ 0.5 in 0 of each. Median max TM
  is 0.101.
- The corrections match: 0.8186 × 0.995 = 0.8145 and 0.3671 × 0.995 = 0.3652.
- ε for the significant set is not above ε for random pairs, so there is no sign of hidden homology enriched by
  significance.

### 3. numbers.json: PASS-with-notes
`/tmp/verify_partA/numbers_eq.py`
- `results/numbers.json["partA"]` is identical to `v2/numbers_partA.json`, and `["partA_descent"]` is identical to
  `descent_numbers.json`.
- **Bug:** in `scripts/partA/11_numbers.py:95`, `itertuples()._asdict()` renames the dotted columns. The F_conv,func
  values therefore sit under the anonymous keys `_13`, `_14`, `_16` and `_17` (point and CI, per family pair and per
  pair). The headline 0.0778 is present but unlabelled.
- CIs are stored as strings such as `"[np.float64(…)…]"`.
- Fix: use `to_dict("records")`.

### 4. Class audit on 30 pairs: PASS
`/tmp/verify_partA/audit30.py` (13 s)
- Sample: 3 significant pairs per class (H1s/x, H2s/x, H12s/x, H3j, H3s, H4, N), seed `random_state=7`.
- For each pair I recomputed the class independently:
  - my own OrthoXML parse and LCA (speciation or duplication);
  - Pfam clan or family, and TED CATH-H regions covering ≥ 50 % of the site residues, from the raw annotation tables;
  - PANTHER family from `uniprot_proteins`;
  - the MMseqs2/Foldseek rules applied to `rest_protein_pairs.parquet`.
- **30/30 agree.** The same-site s/x flag was taken from the table. Group membership for proxy-placed proteins was
  taken from `proteins.parquet`, so the proxy rule itself was not re-derived.

### 5. N-class contamination: PASS-with-notes
`/tmp/verify_partA/contam.py` (2 min) checked all 3,605,586 significant N pairs:

| shared label (raw tables) | sig N pairs | random rest-sample N |
|---|---|---|
| any PANTHER family | **0.0000** | 0.0000 |
| Pfam clan at protein level (anywhere in the protein) | 0.0165 | 0.0056 |
| Pfam family at protein level | 0.0003 | 0.0001 |
| TED CATH-H at protein level | 0.0081 | 0.0043 |
| Pfam clan or TED-H at the **site region** (20k sample) | **0** | — |

- There are no PANTHER or site-region homologs in N.
- 1.6 % share a Pfam clan elsewhere in the protein, about 3× the random rate. The prereg allows this ("Pfam/PANTHER
  absence NOT required"), but multi-domain proteins that share a clan in another domain are a plausible source of
  part of the N excess. A protein-level-clan-excluded sensitivity would cost little.

### 6. Functional convergence: PASS-with-notes
`uv run --project . --with rdkit python /tmp/verify_partA/func.py` (80 s; also `func2.py`)

**Reproduction.** With my own shared-class rule (common token, or ECFP4 Tanimoto ≥ 0.5 between non-metal CCDs) on
the stored site labels:
- test n = 192 (210 before the train exclusion);
- s_sig = 0.1562 (30/192); 192 distinct protein pairs.

**Label availability is the same for test and control.** Both sets require both sites to be annotated, and the
control draws only from annotated sites. "Shared" is defined only when both sides carry tokens.

**Independent control.** I drew a fresh control (seed 12345, 100k random annotated pairs, no MMseqs2/Foldseek
filter, which removes only about 0.2 %), matched on the same μ deciles:

| | stored | mine |
|---|---|---|
| s_ctl, unmatched | 0.0808 | 0.0799 |
| s_ctl, matched | 0.0612 | 0.0609 |
| δ | — | 0.0953 |

The value is stable to resampling.

**Covariate balance after μ matching** (per-covariate means, test vs matched control):

| covariate | test | matched control |
|---|---|---|
| log volume | 5.90 | 5.87 |
| pLDDT | 0.93 | 0.94 |
| aromatic fraction | 0.175 | 0.152 |
| polar fraction | 0.504 | 0.476 |
| SwissProt source | 0.40 | 0.42 |

Matching on the scalar μ leaves residual imbalance in polar and aromatic composition, and these covariates could
plausibly affect ligand class.

Notes:
- **Half of the "shared" calls are Tanimoto calls between nucleotides** (ADP↔AMP/ATP/SAM/CoA/UDP), and many of the
  rest are generic metals (Mg, Zn, Ca). Token-only rates are s_sig 0.083 vs s_ctl 0.029 matched (enrichment about
  2.9×), so the enrichment survives without Tanimoto. The excess is still "nucleotide/metal-binding site" sharing,
  not specific function.
- The matched control is lower than the unmatched one (0.061 vs 0.081), which raises δ. With the unmatched control,
  δ = 0.075 and the CI would very likely include 0.
- F_conv,func multiplies the per-family-pair F_analog by a per-**pair** δ, so the δ is not family-pair weighted.
  D18's wording ("per-family-pair excess") does not describe this.
- The CI lower bound of 0.004 rests on about 30 shared pairs. "Barely excludes 0" is the right wording.

### 7. Directional test: PASS-with-notes
**Arithmetic.**
- 3/2,627,453 = 1.14e-6 and 3/1,842,604 (N) = 1.63e-6, matching `directional_by_class.csv`.
- The 3 passes have p = 2.3e-15, 1.6e-13 and 9.5e-13 (normal-fit tail).

**Toy checks** (`/tmp/pe_big/envs/faiss/bin/python -I /tmp/verify_partA/dir_toy.py`, 20 s, plus `dir_toy2.py`). These
use the repo's own toy world, extended.

| check | result |
|---|---|
| `decide()` p equals the elementwise max of the single-null p's | holds exactly |
| no-convergence world, seeds 0–5: false passes with both nulls | **7 / 19,200 (3.6e-4)**: seeds 0, 1 (the ones the test uses) give 0; seeds 3 and 4 give 3 and 4 |
| false-pass rate, ortholog null alone | 2.0 % |
| false-pass rate, random null alone | 0.04 % |
| share of *selected* tests above q95 (nominal 5 %) | 11 % (random null), 18 % (ortholog null) |

- The nulls are **not conditioned on selection (Δ_sel > 0)**. In the toy, selection shifts mean Δ from 0.022 to 0.037,
  and the random null from 0.000 to 0.032 under its own selection. The q95 criterion is therefore anti-conservative.
- The real-data F_conv,dir (1.1e-6) is below the toy false-pass rate, so the 3 passes are consistent with false
  positives.
- Power was shown only for a strong planted effect (frac 0.5). At frac 0.1 and 0.2 the planted pairs never reach
  significance in the toy, so the power for weaker convergence is unknown.
- Recommendation: compute both nulls conditional on Δ_sel > 0.

**The 3 passing pairs** (`/tmp/verify_partA/cands.py`) all look like artifacts:

| pair | a side | b side |
|---|---|---|
| Q24298:s0 (DE-cadherin) ↔ E1JGJ8:s0 (release factor 2) | 11 residues | **6-residue** footprint, burrow 0.09 |
| P46723:s3 (ALAD) ↔ C0NED7:s2 (pseudouridine synthase) | **46 Å³**, burrow 0.04 | 835 Å³, pLDDT 76, 38 % of residues < 70 |
| Q22937:s1 (GPCR serpentine receptor) ↔ Q6R5P0:s1 (TLR11) | 53 Å³ TM-bundle crevice | 890 Å³, 50 residues, 42 % < 70 |

None has ligand labels. Tiny shallow pockets are paired with large low-pLDDT ones, and the volume mismatches are
8–17×. "No detectable directional convergence" is the correct reading, and the passes should not be presented as
candidates.

### 8. Commit order of D17–D19: PASS-with-notes
The ledger rows were written after the results; the code was mostly written before them. Times are CET.

| event | time |
|---|---|
| directional test code committed (`3f744ad`; includes `normal_p`, max over nulls, "beyond range" rule) | 13:03:42 |
| classify run at `bc706b0` | 13:16 |
| first functional run (`partA.log`) | 13:55:38 |
| first directional run | 13:59:19 |
| category-scaling code change (`3f11d4c`) | 13:59:58 |
| results committed (`1dda927`) | 14:13:42 |
| ledger rows D17–D19 committed (`92a11d7`) | 14:15:50 |

By deviation:
- **D17** (directional p-values): decided in code before any directional result; the code is unchanged since.
  Rerunning gave identical output.
- **D19** (τ for the split): τ = 0.2 was chosen on the tune half at 12:13 (`ed79fdc`), before the v2 code.
- **D18** (F_conv,func): the headline "any" row came from 13:03 code, before the result. But the **category scaling**
  (`frac_c`) was committed at 13:59:58, 4 minutes *after* the first functional output. That output had shown
  category F values above the total (metal 0.20, cofactor 0.56, against F_conv,func 0.078). The change is a sensible
  coherence fix, but the D18 statement "fixed before the result was seen" is false for the category rows. Amend the
  ledger.

### 9. Tune/test split and τ sensitivity: FAIL (reporting gap, non-blocking)
- prereg `sogs.tau_column_jaccard.reported` requires "test half + full grid for every headline metric", and D10
  promises "τ = 0.3 sensitivity reported for every headline metric".
- Neither PARTA.md v2 (r_H1s and the same-site split depend on τ) nor DESCENT.md (M1 and M2 depend on the τ = 0.2
  states) reports test-half or τ = 0.3 values.

I filled part of the descent gap (`/tmp/verify_partA/descent_halves.py`). Split key: sha1 of the smallest accession
in the root HOG (OrthoXML protId, which may differ slightly from the prereg's UniProt-accession key).

| half | M1 ortholog | M1 paralog | binned difference |
|---|---|---|---|
| test | 0.774 | 0.784 | −0.92 pp |
| tune | 0.773 | 0.777 | −0.45 pp |

The conclusions are unchanged. τ = 0.3 was not computed by me; `states_t0.3` exists in bulk.

### 10. Descent numbers: PASS
`/tmp/verify_partA/descent_check.py` (2 s)
- Cap pair counts: 880,234 ortholog and 373,273 paralog.
- M1 by variant from `pairs.parquet`, all matching DESCENT.md and the json:

  | variant | ortholog | paralog |
  |---|---|---|
  | qm90 | 0.7736 | 0.7803 |
  | lit | 0.8547 | 0.8632 |
  | s70 | 0.7581 | 0.7735 |
  | qm90nf | 0.8345 | 0.8380 |
  | qm90fp | 0.7768 | 0.7859 |

- The lit prereg-ceiling normalisation is 1.104, matching the json.
- The binned identity-standardised difference is −0.69 pp, and −1.12 pp for cross-species pairs; both match.

### 11. Descent confounding and D15 symmetry: PASS-with-notes
**Divergence time is strongly confounded with class at matched identity.**
- In every identity bin, orthologs come from species pairs that are 2–18× more divergent than the paralogs' species
  pairs. Example at identity 0.9–1.0: 678 vs 38 Mya.
- 30 % of paralog pairs are within one species (91 % in the 0.9–1.0 bin).
- Identity matching is the preregistered covariate and the standard design.

**The sign of the effect depends on what is controlled:**

| control | ortholog − paralog |
|---|---|
| identity only | −0.7 to −0.9 pp |
| identity × time (cross-species) | −4.6 pp (theirs: −3.95) |
| time only, < 100 Mya bin | **+10 pp** (0.868 vs 0.766) |
| time only, 1,200–2,000 Mya bin | −5 pp |

- Paralog duplication ages are unknown, so a time-matched comparison is ill-posed, as DESCENT.md says.
- "Not supported" is fair. "Orthologs slightly *less* conserved" should not be read as a real effect, and the
  time-only contrast should be shown.

**D15 symmetry.**
- In a 427k-row sample of the τ = 0.2 states, **100 % of A rows pass the qm90 P-criteria** (pLDDT, frac < 70,
  aligned fraction, length ratio), against 48 % of P rows. qm90 is therefore symmetric on every criterion that is
  stored for P.
- Truncation and SP/TM are not stored for P, as disclosed.
- The same mask is applied to ortholog and paralog rows.
- The descent states predate the D11 code (`91d5cbe`, 13:20; the O1 run was left untouched), so the asymmetry from
  D11's footprint > 10 rule does not arise here. It would arise on O1_r2/O2 states unless P is given the same
  > 10-residue requirement.

### 12. Toy tests
| command | result |
|---|---|
| `uv run --project . pytest tests/test_partA.py tests/test_descent.py -q` | **1 failed, 22 passed, 2 skipped** (6 s): `test_covmodel_recovers_size_effect` fails and 2 tests are skipped on `ModuleNotFoundError: torch`, because torch is not in the uv project env |
| `/tmp/pe_big/envs/faiss/bin/python -m pytest tests/test_partA.py tests/test_descent.py -q -p no:cacheprovider` | **25 passed** (13 s) |

Fix: add torch to the project env, or mark the tests `importorskip` consistently.

## Runtime
| check | wall time |
|---|---|
| recompute F | 7 s |
| contamination | 124 s |
| 30-pair audit | 13 s |
| functional (×2) | 160 s |
| directional toys | 25 s |
| candidates | 3 s |
| descent | 5 s |
| hidden homology | 3 s |
| pytest (×2) | 20 s |
| **total compute** | **≈ 6 min** (about 9 min including setup) |
