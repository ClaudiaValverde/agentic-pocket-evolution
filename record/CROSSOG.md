# Cross-orthogroup site matching — exploratory, not preregistered (D31)

Snapshot 20261030T085136Z, development subset, preliminary. Method frozen in `docs/CROSSOG_METHOD.md` (commit 77db16c,
17:31 CET), before any cross-OG output. Code `scripts/partB/crossog.py` (6f5bac5 and later), toy tests
`tests/test_crossog.py` (10 known-answer tests, all pass). Results: `results/partB_crossog/20261030T085136Z/`.
Commit order: method 77db16c → code 6f5bac5 → frozen validation 461e082 → post hoc diagnostic and characterisation.

## Headline
- **Validation recall: 0/7** (O1 0, O2 0) under the frozen method. 5 of the 7 items are in scope; D14/KAI2 and
  lysozyme/α-lactalbumin lie within one OG in both orthologies, and Nurr1 does too on O1.
- The gate of ≥ 3/7 failed, so **no genome-wide run was made** (method §6). As the method allows, I ran a
  characterisation on 2,000 random same-family O1 OG pairs, with real and shuffled site correspondence on the same pairs.
- Cross-OG **site linking works**: 610 linked SOG pairs with the real structural correspondence vs **1** with the shuffled one.
  The **tests on linked sites do not**:
  - (a) gives 0 calls in both. Its p-values are conservative: 1.8 % of links have p < 0.05.
  - (b) gives 10 calls with the real correspondence vs **33** shuffled. Absence at a random position is more common than absence at the matched position, so the (b) calls are at null level.

## Method (summary; full text in docs/CROSSOG_METHOD.md)
- **Eligible SOGs:** O1_r3 / O2_r1 SOGs at τ = 0.2, P in ≥ 3 species, not small or empty footprint.
- **Representatives:** for each SOG, the P member with the highest footprint pLDDT in the Foldseek DB. For each OG, the member with the largest summed pLDDT over its columns.
- **Same-family pairs:** CATH H-level TED labels held by ≥ 50 % of members. Pfam family or clan when an OG has no TED label.
- **Alignment and footprint mapping:**
  - Foldseek 10.941cd33 `structurealign` on exact representative pairs (a prefilter DB written with tsv2db, so no pair is dropped).
  - A pair is structurally supported if alntmscore ≥ 0.5.
  - The SOG footprint is mapped through the rep MSA row, then the CIGAR, then the target rep row into the target OG's MSA columns. Unmapped columns count in the union.
  - J = mean of the two directional column-Jaccards. SOGs are linked if J ≥ τ = 0.2.
  - Only each SOG's best partner is tested.
- **(a) Specificity divergence:**
  - D = median of 1 − cos over cross pairs of the 384-d site fingerprints (≤ 10 members per SOG, one per species).
  - Null: 200,000 within-SOG cross-species pairs at matched identity. The link's identity is the reps' Foldseek fident; null pairs use MSA identity. Bins are 0.1 wide and sparse bins merge downward.
  - Call: BH q ≤ 0.05, and D > both within-SOG medians.
- **(b) Presence/absence:**
  - The footprint is projected onto the other OG's rep. Coverage must be ≥ 80 % and alntm ≥ 0.5.
  - States come from `pocketevo.states.compute_states` (prereg A rules + D11).
  - Rescue k = 1 (6 BURROW settings) is applied to the A rows of pre-rescue passes.
  - Call: P in ≥ 3 species on one side; on the other side no P/W and A in ≥ 3 species.
- **Validation recovery:**
  - Tested pairs = (derived OG, ancestral OG) from `validation_controls.yaml`.
  - Specificity requires an (a) call on a link containing the validation table's anchor SOG.
  - Loss or gain requires a (b) call in the expected direction.

## Validation (results/partB_crossog/.../validation/, frozen method)
| item | O1 | O2 | why not recovered |
|---|---|---|---|
| steroid receptors (ER vs AR/GR/MR/PR) | in scope; 0 links tested | in scope; 0 links tested | footprints line up (J 0.46–0.48 on O1; 0.74 for the ER-OG ↔ MR/GR/AR/PR-OG LBP on O2), but full-chain alntmscore is 0.39–0.46 < 0.5: the long disordered N-terminal domains depress the whole-chain score |
| NR4A/Nurr1 loss | out of scope (one OG) | in scope; 0 (b) calls | the NHR-6 OG footprints projected into NR4A have either P members in NR4A (23) or no A at all (U by the quality rules); alntm 0.24–0.41 |
| PfLDH vs MDH | anchor linked J 0.54, p_a 0.76 | anchor linked J 0.33, p_a 0.53 | fingerprint distance within the matched-identity ortholog null |
| IDH NADP vs NAD | anchor linked J 0.69, p_a 0.54 | anchor linked J 0.94, p_a 0.39 | same |
| GID1 vs HSL | no link ≥ τ (best J 0.19); HSL OG over-merged (242 members) | anchor linked J 0.53, p_a 0.18 | same; (b) every projection has P/W on the HSL side or no A |
| D14 vs KAI2 | out of scope | out of scope | one OG |
| lysozyme vs α-lactalbumin | out of scope | out of scope | one OG |

Recall 0/7. Lenient specificity (any derived–ancestral link) is also 0. Funnels: `validation/O{1,2}_funnel.json`.

**Post hoc diagnostic** (`validation_posthoc_alntm0/`, `--alntm-min 0`, *not* a validation): with the structural
filter removed, recall is still **0/7**. The steroid links are then tested but give p_a 0.31–0.97, and Nurr1 still
has 0 (b) passes. The alntm rule is not the binding constraint; the tests are.

**Disclosure.** Before the validation run, a CIGAR-convention check aligned GR/MR (P04150/P08235) and E. coli IDH /
IDH3A (P08200/P50213). It showed alntmscore 0.42 for GR/MR. The frozen rule was not changed.

## Characterisation on 2,000 same-family O1 OG pairs (results/partB_crossog/.../characterization_sample2000/)
There are 1,172,645 same-family O1 pairs over 14,325 OGs (55 % TED-based). A seeded sample of 2,000 was run twice
through the identical pipeline: once with the real correspondence, once with every alignment's target residues
circularly shifted by a random offset in [L/4, 3L/4]. The (a) cut is the real run's BH threshold.

| step | real | shuffled null |
|---|---|---|
| OG pairs | 2,000 | 2,000 |
| SOG pairs considered / aligned / alntm ≥ 0.5 | 52,044 / 29,904 / 9,826 | same |
| linked (J ≥ 0.2) | **610** (528 SOGs, 335 OG pairs) | **1** |
| tested (best partner, fingerprints) | 606 | 1 |
| (a) calls | 0 | 0 |
| projections / testable | 20,056 / 1,974 | same |
| (b) pre-rescue / after rescue | 20 / **10** (8 OG pairs) | 35 / **33** (27 OG pairs) |

- **Linking is specific:** about 600× over shuffled. Pocket positions are conserved across paralogous OGs of the
  same family often enough to link 528 of 14,284 eligible SOGs in the sample.
- **Test (a) has no power and is conservative.** The real link p-values have a minimum of 0.008 and 1.8 % below
  0.05; BH gives q ≥ 0.88. Linked paralog sites are not more divergent in fingerprint space than within-OG pairs at
  the same identity: the global null is dominated by noisy sites (null median d 0.53–0.61 at identity < 0.3).
- **Test (b) is not specific:** null 33 > real 10, so the empirical FDR is ≥ 1. The 10 real calls are listed in
  `candidates.csv` as audit-level only (e.g. WASHC4/CYFIP1, MSTO1/TUBB8B, ESD/DPP4). They are not evidence of pocket loss or gain.

## Runtime (CPU, 40 workers for rescue, Foldseek 48 threads)
- Validation: O1 52 s, O2 108 s per pipeline pass (+ about 1 min of loading and the null build).
- Characterisation: 216 s real + 189 s null for 2,000 pairs (Foldseek 6 s for 36,000 rep pairs; link and state
  computation are serial Python, 75 + 61 s; rescue about 73 s).
- Extrapolated full O1 (1.17 M pairs): about 35 h serial. Parallelising the per-pair loop would bring it to roughly
  1 h on 40 cores. It was not run (gate failed).

## Limitations
- The **gate failed (0/7)**. Nothing here is a discovery. All outputs are exploratory and not preregistered (D31).
- **Whole-chain structural support** (alntmscore) penalises multi-domain proteins with disordered regions (NRs). A
  domain-restricted score would fix the nuclear receptor items. But the post hoc run shows recall stays 0 even without the filter.
- **The (a) null is global.** A site-specific statistic would be more sensitive. Example: D relative to the two
  SOGs' own within-OG spread (R = D / max(W_x, W_y)).
  - Hypothesis-generating only: the O2 ER ↔ 3-ketosteroid-receptor LBP link has R = 3.4, above the 98.5th
    percentile of R among sampled real links (1.5 % have R ≥ 3.4).
  - No other validation link has R > 1.5.
  - A v2 would need to be preregistered and validated on items not used here.
- **(b) has no detection-error null.** The shuffled null shows absence at a random position is common. Even after
  rescue, (b) calls cannot separate a real loss from a pocket-finder miss.
- **Coverage gaps:**
  - The interface rule, the TED two-domain rule and P-side rescue were not applied to the new A calls.
  - O2 reps are restricted to tier 1/2 models; tier 3 and tier L are not in the Foldseek DB.
  - Only one structural alignment is made per rep pair.
- **Orthology granularity decides scope** (D28): the merged O1 NR OG HOG0025274 holds ER, GR, AR and NR4A.
  Over-merged OGs (GID1 / HSL HOG0029749) dilute footprints and inflate the within-OG spread.

## Files
- `docs/CROSSOG_METHOD.md`: frozen method.
- `validation/`: frozen-method output (validation_recall.{csv,json}, links / tested / b parquets, funnels, run.log).
- `validation_posthoc_alntm0/`: post hoc diagnostic.
- `characterization_sample2000/`:
  - funnel.csv and O1_funnel.json (real all / real on the sample / null);
  - candidates.csv (audit-level (b) calls);
  - real and null links / tested / b parquets;
  - O1_og_pairs.parquet;
  - run.log.
- Bulk (regenerable): `/tmp/pe_big/partB_crossog/20261030T085136Z/` (Foldseek work dirs, rescue models, null_a_O{1,2}.parquet).
