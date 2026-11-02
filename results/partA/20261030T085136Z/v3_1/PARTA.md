# Part A v3.1: both weightings side by side, and descent M1 with ledger D11 applied symmetrically. Snapshot 20261030T085136Z

**Status: development subset, preliminary.** Not yet verified. Written 2026-10-30 17:07 CET (`date`).
- This answers the PI's request after ledger D29: (a) show both weightings everywhere, (b) apply D11 symmetrically.
- v3 and v2.1 are not modified. Decomposition and functional numbers are read from `../v3/` (O2) and `../v2_1/` (O1).
  The descent outputs are new and live here.
- Code: `3966a91` (descent `qm90d11` + toy test) and `4f94aed` (`--variant`). Bulk descent tables:
  `/tmp/pe_big/partA_descent_v31/`.

**How to read the two weightings.**
- **Per pair** is the granularity-robust statement. It does not depend on how finely the orthology partitions
  families: O1 and O2 give identical values.
- **Per family pair** is the preregistered headline. It moves with the family partition (see v3: O2 packs H3j pairs
  into 3.6× fewer cells).

## 1. Decomposition (pooled tail, E < 0.01 in both directions; nominal expected false positives; 95 % bootstrap CIs)
| run | subset | F_anc, per pair | F_analog, per pair | F_anc, per family pair (prereg headline) | F_analog, per family pair (prereg headline) |
|---|---|---|---|---|---|
| O1 v2.1 | all | 0.615 [0.566, 0.661] | **0.367 [0.323, 0.416]** | 0.176 [0.169, 0.185] | 0.819 [0.810, 0.826] |
| O2 v3 | all | 0.616 [0.567, 0.658] | **0.367 [0.325, 0.416]** | 0.060 [0.057, 0.063] | 0.935 [0.932, 0.937] |
| O1 v2.1 | test half (hash split) | 0.695 [0.636, 0.737] | 0.285 [0.245, 0.340] | 0.175 [0.163, 0.189] | 0.819 [0.805, 0.830] |
| O2 v3 | test half (D20, crossover excluded) | 0.849 [0.802, 0.881] | 0.136 [0.104, 0.181] | 0.083 [0.077, 0.091] | 0.909 [0.900, 0.916] |
| O1 v2.1 | both pLDDT ≥ 90 | 0.820 [0.794, 0.840] | 0.158 [0.134, 0.181] | 0.401 [0.386, 0.418] | 0.588 [0.571, 0.604] |
| O2 v3 | both pLDDT ≥ 90 | 0.821 [0.796, 0.844] | 0.158 [0.135, 0.178] | 0.116 [0.109, 0.123] | 0.870 [0.862, 0.878] |

- The per-site weighting (F_analog 0.700 on both O1 and O2) is in the estimand CSVs.
- **The per-pair half values are not comparable with the "all" row.** Cross-half pairs and the excluded crossover
  families (63 % of sites on O2) drop out, and those are mostly N pairs.
- The per-class rates r_c are weighting-free. r_N is 7.0e-7 on both O1 and O2. r_H1s is 0.010 on O1 and 0.013 on O2.

## 2. Functional convergence F_conv,func (both weightings; pocketcompare training proteins excluded)
| run | subset | labels | n | F_conv,func, per pair | F_conv,func, per family pair (prereg headline) |
|---|---|---|---|---|---|
| O1 v2.1 | all | full | 192 | 0.035 [0.002, 0.078] | 0.081 [0.007, 0.178] |
| O2 v3 | all | full | 191 | 0.035 [0.005, 0.075] | 0.095 [0.021, 0.204] |
| O1 v2.1 | all | specific ligands (tokens, no metal) | 149 | 0.011 [0, 0.034] | 0.021 [0, 0.062] |
| O2 v3 | all | specific ligands (tokens, no metal) | 148 | 0.011 [0, 0.035] | 0.025 [0, 0.089] |
| O1 v2.1 | pLDDT ≥ 90 | full / specific | 100 / 79 | 0.011 / 0 | 0.044 / 0 |
| O2 v3 | pLDDT ≥ 90 | full / specific | 99 / 78 | 0.011 / 0 | 0.068 / 0 |
| O1 v2.1 | test half | full | 41 | 0.038 [0, 0.110] | 0.127 [0, 0.317] |
| O2 v3 | test half (D20) | full | **2** | not interpretable | not interpretable |

F_conv,dir stays ≈ 0 in both weightings (v2.1: 1.1e-6 per pair, 1.9e-6 per family pair; 3 passes, all flagged as
artifact-like).

## 3. Descent M1 with D11 applied symmetrically (`descent_{O1,O2}_{qm90,qm90d11}_summary.csv`)
- **Problem.** D11 (in the O1_r2 and O2_r1 states) makes an absence call A require a footprint that spans more than
  10 residues. P had no such rule, so small-footprint non-matches disappeared and M1 was inflated (v3 O2: 0.93).
- **Fix.** `qm90d11` (`src/pocketevo/descent.py`, documented in its docstring) applies the > 10-residue rule to every
  informative Y member, P or A, on both orthologies.
  - Footprint residues = SOG footprint columns × the member's non-gap fraction.
  - Toy test: `tests/test_descent.py::test_d11_asymmetry_inflates_m1_and_symmetric_rule_does_not`. The true rate is
    0.5; the asymmetric rule gives 0.67 in the X→Y direction and the symmetric rule gives exactly 0.5.
- **Run.** The pair tables were rebuilt with the unchanged builder plus the new variant: O1 at τ 0.2 with the verified
  O1 states, O2 via `17_descent_o2.py`.
- **Halves.** O1 test half = prereg hash on O1. O2 test half = D20 table, crossover excluded. Capped sampling; 200 OG
  bootstrap reps.

| variant | orthology | subset | M1 ortholog | M1 paralog | ortholog − paralog at matched identity | ortholog / paralog pairs |
|---|---|---|---|---|---|---|
| qm90 (asymmetric on O2) | O1 | all | 0.774 | 0.780 | −0.76 pp [−1.36, −0.19] | 708k / 294k |
| qm90 | O1 | **test** | 0.774 | 0.784 | −0.94 pp [−1.81, −0.21] | 358k / 151k |
| qm90 | O2 | all | 0.932 | 0.936 | −0.38 pp [−0.63, −0.14] | 1.33 M / 520k |
| qm90 | O2 | **test** | 0.928 | 0.931 | −0.23 pp [−0.61, +0.17] | 343k / 117k |
| **qm90d11 (symmetric)** | O1 | all | 0.889 | 0.897 | −0.90 pp [−1.33, −0.38] | 481k / 193k |
| **qm90d11** | O1 | **test** | **0.889** | **0.900** | −1.10 pp [−1.76, −0.48] | 245k / 101k |
| **qm90d11** | O2 | all | 0.890 | 0.892 | −0.16 pp [−0.49, +0.18] | 984k / 370k |
| **qm90d11** | O2 | **test** | **0.885** | **0.888** | −0.22 pp [−0.83, +0.51] | 260k / 87k |

Reading:
- **With D11 symmetric, O1 and O2 agree** (test half 0.889 vs 0.885 for orthologs), so the earlier O2 jump to 0.93 was
  the asymmetry. Both symmetric values sit above the O1 qm90 0.774 for a different reason: small-footprint SOGs, which
  have low retention, now leave the denominator on both sides. The symmetric rule therefore measures retention of
  sites with footprints larger than 10 residues.
- **The ortholog conjecture is not supported for presence** under every variant and on both orthologies. The
  difference is between −1.1 and +0.5 pp, never positive beyond its CI. This is not a claim that orthologs are less
  conserved.
- O2 has no TimeTree table for its non-tier-1 species, so the time-stratified contrasts are reported for O1 only (v2.1,
  `descent_2d_summary.csv`).

## Text corrections to v3 (verifier `reviews/verify/verify_G4.md`, c6c953c; numbers unchanged, edited 2026-10-30 17:28 CET)
1. **The family-pair shift is not caused only by tier-2 placement.** O2's root groups are coarser in general: F_analog
   per family pair also moves from 0.886 to 0.915 on tier-1-only pairs (verifier). v3's explanation ("O2 places tier-2
   proteins into fewer, larger groups") was too narrow. The general reason is that O2 has coarser groups for all
   proteins.
2. **v3's section "Pairs that move class because of real placement" should read "Pairs that move class between O1 and
   O2".** 39 % of the 1,529,928 class changes involve no tier-2 protein.
3. **N moves.** Report 165 pairs out of N and 190 pairs into N between O1 and O2, across all classes. This replaces v3's
   "63 out, 69 in", which counted H12↔N only.
4. **Test-half composition.** 2.7 % of the O2 D20 test-half sites belong to PANTHER-proxy proteins (proxy_hog) that
   inherit their assigned HOG's half.

## Deviations (PA ids)
- **PA28:** `qm90d11` is a new descent variant (symmetric D11). It is reported beside qm90; qm90 stays the prereg-literal
  primary. This changes the reported descent M1 level, not the conjecture conclusion. **For the ledger.**
- **PA29:** the PARTA.md presentation now shows both weightings. Per pair is labelled granularity-robust; per family
  pair stays the prereg headline (PI decision).

## Runtime
| step | wall |
|---|---|
| descent pair rebuild, O1 + O2 | 98 s |
| 4 descent 2-D runs | 103 s |
| tests | under 30 s |
| **total elapsed** | 17:02–17:08 |
