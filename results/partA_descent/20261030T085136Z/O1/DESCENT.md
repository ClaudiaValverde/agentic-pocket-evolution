# DESCENT: do pockets descend by orthology? (Part A §4.1, M1 and M2)

Snapshot 20261030T085136Z · orthology O1 (Ancestra HOGs, 46 tier-1 proteomes) · SOG τ = 0.2 · prereg sha256
`52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83` · seed 0 (prereg meta.seeds[0]; bootstrap seed
rule `seed*1000+rep`, 200 reps over root OGs). Development subset: results are preliminary. Every number below comes
from `descent_numbers.json` (also in `results/numbers.json["partA_descent"]`), `m1_rates.csv`, `m2_cosine.csv` or
`controls.csv`. Producer: Part A descent agent. **An independent agent has not yet verified these results.**
M3 (tree consistency) and §4.2 (decomposition) are done by other agents and are not covered here.

## Headline
1. **Pocket positions are conserved at about 0.77 of the measurable rate.** Take a pocket site in protein X and an
   informative footprint in a homologous protein Y. The site is found at the same position in Y with probability
   M1 = **0.774 [0.769, 0.777]** for ortholog pairs and **0.780 [0.775, 0.787]** for paralog pairs. These figures use
   the primary, quality-matched definition (qm90), capped pair sampling and 95 % OG-bootstrap CIs.
   - Dividing by the same-protein re-detection ceiling, computed with the same observation rules, gives
     **0.85 [0.84, 0.85]** for both classes.
   - M1 decays with divergence. For orthologs it goes from 0.86 (< 100 Mya) to 0.73 (700–1,200 Mya), and with
     identity it goes from 0.89 (90–100 % identity) to 0.63 (10–20 %).
   - Below 10 % identity it falls to 0.14, which reflects alignment failure as much as biology.
2. **The data do not support a "pocket ortholog conjecture".** At matched global identity (identity spline plus class,
   binomial GLM), orthologs keep the pocket position slightly *less* often than paralogs:
   - marginal difference **−0.9 pp [−1.4, −0.3]**;
   - OR **0.95 [0.92, 0.98]**;
   - binned standardised difference −0.7 pp [−1.2, −0.1].

   The effect is tiny and its sign holds in every capped sensitivity (table below). The one exception is models with
   no atlas-filtered pockets (qm90nf), where it is −0.5 pp [−1.1, +0.2], i.e. null.
3. **Fingerprints of matched sites are far above the non-homologous background and decay with time.** The M2 median
   cosine for ortholog matches is 0.60 overall. The reference lines:

   | reference | median cosine |
   |---|---|
   | size-matched non-homologous background | 0.13 (p95 0.51) |
   | same protein, AFDB v4 vs v6 model | 0.84 |
   | same model, BURROW seed swap | 0.99 |

   Ortholog matches go from 0.77 (< 100 Mya; 90 % of the background-to-floor distance retained) to 0.56
   (> 1,200 Mya; 60 % retained). At matched identity, ortholog matches are slightly *more* similar than paralog
   matches: mean cosine difference **+0.015 [+0.012, +0.018]**. So the conjecture gets weak support from shape and
   chemistry, but none from position.

## Method (DESIGN v2 §4.1; prereg part_a.M1/M2; code `src/pocketevo/descent.py`)
- **Pairs.** All member pairs of the same O1 root OG that have at least one informative observation. The relation
  comes from the HOG tree's LCA node: speciation means ortholog, duplication means paralog
  (`homology.parse_orthoxml` + `lca_type`, as in Part A §4.2).
  - Fair sampling: per OG, ≤ 200 pairs (seeded per OG with `[seed, crc32(og_id), 1]`), each with weight 1. This is
    **"cap"**, the headline, so that big families cannot dominate.
  - **"Uncapped"**: ≤ 5,000 pairs per OG with Horvitz–Thompson weight n_candidates / n_sampled. This estimates the
    all-pairs rate.
  - Counts: 880,234 ortholog and 373,273 paralog cap pairs in 22,834 OGs; 22.7 M candidate pairs in total.
  - Pair counts are exact: sparse member × SOG products.
- **M1 observation.** For a pair (X, Y) and an SOG in which X is **P**:
  - The observation is informative if Y's footprint is informative, and a match if Y is P.
  - Both directions are counted, so a P/P SOG gives 2 of 2.
  - Y members that are U (partial overlap, low alignment, no model, …) leave the denominator.
- **Variants of "informative".**
  - `lit`: Y ∈ {P, A}, prereg-literal. A = A_primary (footprint pLDDT ≥ 90 …).
  - `qm90` (**primary**): Y = A, or Y = P with the A_primary quality criteria available for P rows. These are mean
    footprint pLDDT ≥ 90, frac < 70 ≤ 0.2, aligned ≥ 0.8 and length ratio 0.7–1.3. Truncation and SP/TM are not
    stored for P rows.
    - Why: in `lit`, a low-confidence Y counts when it is P but is never A, which inflates M1 by about 8 pp
      (0.855 vs 0.774). This is an interpretation, not a new metric: `lit` is reported alongside. It is the
      conservative direction.
  - `s70`: the prereg pLDDT ≥ 70 sensitivity stratum, quality-matched the same way.
  - `qm90nf`: qm90 restricted to Y models with no atlas-filtered pockets. This applies to P and A alike. A calls on
    the other models are `rescue:pending` (SOGS.md).
  - `qm90fp`: drops `small_footprint` SOGs (≤ 3 columns, decision D8).
- **Ceiling normalisation** (prereg: stratified by pLDDT band × model source). Each informative observation carries
  c(band of Y footprint, size of X site ≤ 10 / > 10 residues). The model source is fixed: Foldcomp v4 vs v6 atlas,
  τ = 0.2, atlas filters. This is the "different prediction of the same sequence" comparison, which is the relevant
  noise for two v6 models.
  - (a) **Prereg ceiling**: `error_rates.parquet` redetect rates, 0.56–0.87.
  - (b) **Rule-matched ceiling**, from `redetection_sites.parquet`: P = best Jaccard ≥ τ and A-like = no overlap. Partial
    overlaps (0 < J < τ) leave the denominator, exactly as `partial_site_overlap` members are U in the SOG states.
    This gives 0.71–0.97.
  - (a) is not on the M1 scale: M1/(a) exceeds 1 for close orthologs (up to 1.08 for qm90 and 1.23 for lit). (b) is
    the interpretable normalisation; both are reported.
- **Divergence time.** TimeTree 5 pairwise ages for the 1,035 tier-1 species pairs, fetched once
  (`data/timetree/…_pairwise.csv`, in DATA.md); 552 pairs have a TimeTree age. All 479 pairs whose MRCA is the root
  of cellular life (any Bacteria/Archaea vs Eukaryota/Archaea) get no TimeTree study. Fallbacks:
  - The **rank proxy is 4,200 Mya** (LUCA, Moody et al. 2024), flagged `rank_proxy`; named fallback counter
    `timetree_to_rank_proxy` = 479.
  - 3 pairs use a node median.
  - 1 pair is missing (Herpotrichiellaceae).
  - TimeTree ages are not ultrametric on the NCBI guide tree. For example, Bacteria (1,909) is younger than
    Pseudomonadati (2,849). They are used as returned.
- **Identity.** Identical / both-non-gap columns of the OG FAMSA MSA. This is the same definition as
  `sog.mean_pairwise_identity`.
- **Ortholog-conjecture test.**
  - Cells: OG × class × identity bin (width 0.02), with successes = matches and trials = informative observations.
  - Model: binomial GLM (IRLS) with logit = cubic B-spline(identity; knots 0.2/0.3/0.4/0.55/0.75) + class.
  - Effect: the marginal ortholog−paralog difference, standardised to the paralog identity distribution, plus the OR.
  - Also reported: a binned standardised difference (10 identity bins).
  - CIs: bootstrap over root OGs, 200 reps.
- **M2.** Cosine of the two sites' representative fingerprints (the Part A universe `fp.f16.npy`, tier 1) for every
  P/P observation of a cap pair: 1.55 M ortholog and 0.64 M paralog matches.
  - Noise floor: `redetection_sites.parquet` matched sites (J ≥ τ), v4 vs v6 and seeds.
  - Background: 10⁵ random tier-1 site pairs (seed 0). One site is drawn from the matched sites; the partner has the
    same n_res bin, a different root group and no shared protein-level label (Pfam clan / CATH SF / PANTHER;
    Part A `context.npz`).
  - "Retained" = (median − background median) / (floor median − background median).
- **Toy and unit tests**: `tests/test_descent.py`, 9 pass. They check:
  - exact counts on known two-member worlds;
  - that the quality matching removes low-pLDDT P;
  - the symmetric filtered-model rule;
  - site-pair extraction;
  - a synthetic pair set with known rates (0.6 vs 0.4 recovered to ±0.02). In a confounded world with no class effect,
    the naive difference is > 0.1 but the spline model gives ≈ 0;
  - the bootstrap seed rule;
  - identity;
  - the ceiling loaders.

## M1 results (cap, OG-bootstrap 95 % CI)
| variant | ortholog M1 | paralog M1 | ortholog M1/ceiling (rule-matched) | paralog M1/ceiling (rule-matched) |
|---|---|---|---|---|
| qm90 (primary) | 0.774 [0.769, 0.777] | 0.780 [0.775, 0.787] | 0.847 [0.843, 0.851] | 0.853 [0.847, 0.859] |
| lit (prereg-literal) | 0.855 [0.852, 0.858] | 0.863 [0.859, 0.868] | 0.952 | 0.960 |
| s70 (pLDDT ≥ 70 stratum) | 0.758 [0.754, 0.762] | 0.774 [0.768, 0.779] | 0.850 | 0.864 |
| qm90nf (no filtered pockets) | 0.835 [0.830, 0.839] | 0.838 [0.831, 0.844] | 0.919 | 0.920 |
| qm90fp (footprint > 3) | 0.777 [0.773, 0.781] | 0.786 [0.780, 0.792] | 0.850 | 0.858 |
| qm90, **uncapped** (HT-weighted) | 0.809 [0.795, 0.826] | 0.875 [0.809, 0.907] | 0.880 | 0.941 |

**By divergence time** (qm90, cap, orthologs; rule-matched normalisation in brackets):

| divergence (Mya) | M1 | M1/ceiling |
|---|---|---|
| 0.5–100 | 0.858 | 0.949 |
| 100–300 | 0.812 | 0.892 |
| 300–700 | 0.781 | 0.857 |
| 700–1,200 | 0.733 | 0.808 |
| 1,200–2,000 | 0.750 | 0.822 |
| > 2,000 (LUCA proxy) | 0.810 | 0.875 |

- The deepest bin rises. These are prokaryote–eukaryote orthologs, i.e. OGs that span domains, which are mostly core
  metabolic enzymes. That is a selection effect, not slower turnover.
- Paralogs within a species give 0.801. Their duplication age is unknown (species divergence 0 is only a lower
  bound), so paralog-vs-time curves are not interpretable as decay curves.

**By identity** (qm90, cap):

| class | 10–20 % | 20–30 % | 30–40 % | 50–60 % | 70–80 % | 90–100 % |
|---|---|---|---|---|---|---|
| orthologs | 0.63 | 0.74 | 0.78 | 0.82 | 0.86 | 0.89 |
| paralogs | 0.63 | 0.76 | 0.79 | 0.82 | 0.86 | 0.87 |

## Pocket ortholog conjecture (ortholog − paralog at matched identity)
| analysis | marginal Δ (pp) [95 % CI] | OR [95 % CI] | binned Δ (pp) |
|---|---|---|---|
| **qm90, cap (primary)** | **−0.88 [−1.43, −0.30]** | **0.948 [0.916, 0.982]** | −0.69 [−1.24, −0.11] |
| lit, cap | −0.90 [−1.28, −0.50] | 0.926 [0.896, 0.958] | −0.77 |
| s70, cap | −1.62 [−2.25, −0.98] | 0.909 [0.875, 0.944] | −1.38 |
| qm90nf, cap | −0.48 [−1.09, +0.16] | 0.964 [0.920, 1.012] | −0.37 [−0.96, +0.31] |
| qm90fp, cap | −1.07 [−1.63, −0.50] | 0.936 [0.905, 0.970] | −0.88 |
| qm90, cap, cross-species pairs only | −1.37 [−1.96, −0.82] | 0.922 [0.891, 0.954] | −1.12 |
| qm90, cap, cross-species, + time bins | −3.95 [−4.69, −3.35] | 0.798 [0.766, 0.826] | n/a |
| qm90, uncapped | −7.45 [−9.77, −2.35] | 0.559 [0.442, 0.846] | −8.02 |

Notes:
- Uncapped estimates are dominated by a few very large families: 32 M paralog vs 9.5 M ortholog weighted observations.
  Their bootstrap is unstable (lit/qm90nf uncapped OR CIs degenerate, see json). They are reported, not interpreted.
- Adding time to cross-species pairs makes the effect more negative. Paralogs of a given species divergence are older
  than that divergence, so this model over-adjusts and is shown only as a sensitivity.
- **Reading.** With identity held fixed, the orthology relation adds nothing to positional pocket conservation. The
  point estimate favours paralogs by about 1 pp, which is far below anything biologically meaningful and is within
  plausible residual confounding (e.g. within-species paralogs share proteome-level model properties). This mirrors the
  weak and contested evidence for the function-level ortholog conjecture (Nehrt 2011 vs Altenhoff 2012; LITERATURE.md).

## M2 results (cap; median cosine of matched sites, 95 % OG-bootstrap CI)
| stratum | ortholog | paralog |
|---|---|---|
| all | 0.601 [0.598, 0.604] (63 % > background p95) | 0.592 [0.588, 0.597] (62 %) |
| < 100 Mya | 0.773 (retained 0.90) | 0.605 |
| 100–300 Mya | 0.671 (0.76) | 0.593 |
| 300–700 Mya | 0.621 (0.69) | 0.555 |
| 700–1,200 Mya | 0.574 (0.62) | 0.557 |
| > 1,200 Mya | 0.56 (0.60) | 0.56 |
| identity 10–20 % → 90–100 % | 0.447 → 0.822 | 0.438 → 0.789 |

Reference lines:

| reference | median cosine | spread |
|---|---|---|
| same protein, v4 vs v6 | 0.843 | p05 0.441 |
| same model, BURROW seed swap | 0.986 | p05 0.60 |
| non-homologous size-matched background (n = 10⁵) | 0.130 | p95 0.509, p99 0.650 |

- Matched-identity difference (paralog weights) is +0.015 [+0.012, +0.018] for ortholog minus paralog.
- Within a 100 Myr window, ortholog matches retain almost all of the technical-floor similarity. Beyond about
  700 Mya the median sits at about 0.56, still 3.5× the background median.

## Positive controls (all pairs of the control OGs; `controls.csv`, `figures/controls_m1.png`)
| family (OG) | ortholog M1 qm90 (n pairs) | paralog M1 qm90 (n pairs) | M2 median cos ortho / para, known-pocket SOG |
|---|---|---|---|
| LDHA NADH/substrate (HOG0035618) | 0.966 (968) | 0.982 (306) | 0.61 / 0.54 (SOG g1) |
| CDK2 ATP site (HOG0027286) | 0.687 (612) | 0.662 (25) | 0.52 / 0.66 (SOG g0) |
| ESR1/AR ligand pocket (HOG0025274) | 0.586 (8,119) | 0.561 (9,727) | 0.56 / 0.53 (SOG g0) |

- M1 for the control OGs is computed over all of their SOGs, not only the known pocket.
- The known pockets form one multi-species SOG each (SOGS.md).
- LDHA behaves as a near-perfectly conserved pocket.
- The kinase and NR OGs are lower. They are large, multi-domain, low-identity OGs (median identity 0.52 and 0.23)
  with many secondary sites.

## Caveats
- **Circularity.** SOGs are built from site positions in the OG MSA. M1 therefore measures positional pocket
  conservation *as resolved by FAMSA*, and alignment error lowers it at low identity.
- **Partial overlaps are excluded.** Y members whose site overlaps the footprint below τ are U (SOGS.md choice a), so
  M1 is conditional on a clear yes/no. The rule-matched ceiling applies the same exclusion.
- **A calls are provisional.** 78 % of A rows are on models with atlas-filtered pockets, and rescue has not been
  applied to bulk A. In qm90nf, where every Y model has all pockets kept, M1 is higher (0.835) and the conjecture
  effect is null.
- **No interface rule.** The interface rule is not applied to M1 bulk rates (prereg; D6).
- **Tier 1 only, one model source.** Everything here is tier 1 (v6) and O1. Deep (Eukaryota-level) orthology is
  flagged uncertain in DESIGN §2. O2 and tier L reruns use the same commands with `--orthology` /
  `--orthology-version` / `--snapshot`.
- **Ceiling scale.** The ceilings come from a 2,000-protein sample (kingdom-pooled) and a slightly different matching
  rule: residue Jaccard to the best site, not footprint-column Jaccard.

## Reproduce (runtime on a shared 64-core node)
```bash
P=<repo>; cd $P
.venv/bin/python scripts/partA_descent/fetch_timetree.py                       # once, 455 s (network)
.venv/bin/python scripts/partA_descent/build_descent_pairs.py --workers 24     # 25 s -> /tmp/pe_big/partA_descent/<snap>/O1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/partA_descent/analyze_descent.py --workers 22  # 355 s
.venv/bin/python -m pytest -q tests/test_descent.py
```
Bulk (regenerable): `/tmp/pe_big/partA_descent/20261030T085136Z/O1/{pairs,m2_matches,m2_background}.parquet`
(+ sidecars). Figures: `figures/m1_vs_time.png`, `m1_vs_identity.png`, `m2_cosine.png`, `controls_m1.png`, each
with a `.json` sidecar.
