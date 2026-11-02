# Specificity change on O2 (O2_r1, run O2_r1_spec2): negative

Generated 2026-10-30 23:23 CET. Producer: Part B specificity agent; **not yet verified**. Script `scripts/partB/run_specificity.py --orthology O2 --release O2_r1 --run-id O2_r1_spec2 --tau 0.2 --workers 16 --n-null-seeds 3` (code commits 5fb62e4, 331928f, 8e8d456); prereg sha256 `52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83`; frozen thresholds unchanged (noise floor q95 = 0.3942 from 1798 BURROW-seed replicates, ≥2 BLOSUM62<0 clade-consistent substitutions, BH α 0.05). Runtime 312 s (CPU, 16 workers).

## Result

**0 sister pairs pass BH, in real data and in all 3 pseudo-clade null seeds.** Smallest q = 0.210 (smallest empirical p = 1.8e-05 over 19,102 tested pairs; 44 pairs have p < 0.001, 19 expected if uniform). This matches O1_r3 (0 pass BH; 513 residue+shift pairs vs 71–88 null).

| input | seed | n_pairs | residue_gate | residue_and_shift | bh_tested | called | called (≥3 species) |
|---|---|---|---|---|---|---|---|
| real |  | 743,349 | 19,462 | 536 | 19,102 | 0 | 0 |
| null_pseudo_clade | 0 | 743,349 | 5,596 | 84 | 5,461 | 0 | 0 |
| null_pseudo_clade | 1 | 743,349 | 5,295 | 66 | 5,199 | 0 | 0 |
| null_pseudo_clade | 2 | 743,349 | 5,493 | 80 | 5,394 | 0 | 0 |

Real vs null mean: residue gate 19,462 vs 5461 (3.6×); residue + shift 536 vs 77 (7.0×). The fraction of residue-gated pairs that also clear the shift floor is 2.8% real vs 1.4% null. The excess is mostly in the residue gate: real clades share substitutions and pseudo-clades do not, by construction. The shift test against the matched-depth empirical null gives no significant pair.

## Artifact cascade

No pair passed BH, so nothing is FDR-supported. As on O1, the 5 exploratory pairs with the smallest q went through `run_artifacts.py --orthology O2` (in a scratch results root, `/tmp/pe_big/partB/spec_O2_r1_spec2_cascade`, so the committed O2_r1 candidates files are untouched). They are **exploratory, not significant**, and no claims are made.

| sog_id | level | n_d | n_s | species_d | n_subst | substitutions | shift | p_emp | q_bh | pass_artifacts | strict | fails |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HOG0043018:g0 | Ecdysozoa | 3 | 3 | 3 | 3 | 13:Q>A;17:G>V;18:Q>Y | 1.166 | 1.8e-05 | 0.210 | False | False | a_fragment_present;a_orthology_panther |
| HOG0044827:g5 | Sar | 2 | 2 | 1 | 2 | 5:V>D;9:D>A | 1.131 | 4.7e-05 | 0.210 | False | False | a_fragment_present |
| HOG0044697:g11 | Rhabditida | 3 | 1 | 3 | 5 | 0:F>D;1:P>R;4:C>S;6:F>V;8:P>H | 1.072 | 3.1e-05 | 0.210 | False | False | a_fragment_present;a_orthology_panther |
| HOG0036068:g2 | Gnathostomata | 2 | 8 | 2 | 4 | 0:R>L;1:V>Y;7:T>L;8:G>K | 1.019 | 5.5e-05 | 0.210 | True | False |  |
| HOG0049595:g1 | Sordariomycetes | 2 | 8 | 2 | 2 | 6:G>T;7:C>G | 0.920 | 4.4e-05 | 0.210 | True | False |  |

2/5 pass the cascade, and 0/5 pass it strictly (quality flags). Instability, species-tree lineages and PAE are `na` for specificity on O2.

## Fingerprints (unblocking)

- A site's fingerprint is by definition the atlas row of its representative pocket (`scripts/build_sites.py`, `site_fp/`). Tier-3 sites had been built with `--no-fingerprints`, so `scripts/partB/extract_site_fp.py` read those atlas rows. It covered every P/W tier-3 member of O2_r1 (rescue-overlay states, τ 0.2): 516,652 sites from 143 atlas units, 0 missing, in 10 s. Atlas files matched the snapshot by size and mtime. The method record is pocketcompare 1.0.0, config sha256 `b3ea8fe963ccc180…`, identical to the snapshot's. Nothing was recomputed on GPU. Bulk output: `/tmp/pe_big/sites/20261030T085136Z/fingerprints_O2/tier3/`. This is a superset of the members of ≥3-species SOGs, because the empirical null pool also uses pairs from 2-species SOGs.
- Tier L (tier 4) was not blocked: `site_fp/tier4` already held fingerprints for all 1,837,819 tier-L sites with the same method record (CPU backend). Tiers 1 and 2 were also complete.
- Reproducibility: 200 seeded tier-1 O2 member sites were recomputed fresh with `pocketcompare describe --af-mode` on GPU 3 (CUDA_VISIBLE_DEVICES=3; 200 structures in 12.6 s, 0.063 s/structure). Cosine to the stored fingerprints: min 1.0000000, max|Δ| 6.1e-05 (float16 storage). The extracted tier-3 rows are bit-identical (max|Δ| 0) to the 2,329 overlapping sites in the earlier with-fingerprint tier-3 sample. Tier-3 models could not be re-described, because Foldcomp is not installed (docs/ENV.md).

## Caveats

- 50,309 W members (rescued, no atlas site: site_id null) have residues but no fingerprint, so they are left out of the shift statistic. O1_r3 had no such members.
- Fingerprints are keyed by (model tier, site_id): 88,260 site_ids collide between tier 3 and tier L. An earlier run (O2_r1_spec1, discarded and not committed) keyed by site_id only, and its funnel was identical.
- Model-source covariate (prereg regress_out) = 1 for AFDB v4 / Foldcomp models (tiers 3, L) and 0 for v6 (tiers 1, 2). It was constant on O1.
- 3 null seeds, as in O1_r3 (D24 asks for ≥ 20 for gain and loss FDR; with 0 calls in real data, more seeds would not change the conclusion).
- Design scope (D31): within-orthogroup sister pairs only. Most documented specificity changes are between paralogous orthogroups.

