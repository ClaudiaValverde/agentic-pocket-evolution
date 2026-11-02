# Site orthogroups (SOGs) and observation states: snapshot 20261030T085136Z, orthology O1

Development subset, preliminary (tier 1 = 46 complete proteomes; the atlas is still growing). prereg sha256
`52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83`. Built 2026-10-30.
Code: `src/pocketevo/{sog,states,afmodel}.py`, `scripts/build_sogs.py`, `scripts/sog_tau_benchmark.py`,
`scripts/sog_slice_report.py`. Tests: `tests/test_sog.py`, 22 pass (`pytest.log`).
Every output has a `<file>.json` sidecar with `snapshot_id, prereg_sha256, git_commit, seed` (seed 0 =
prereg seeds[0]), plus `code_sha256`.

## Method (DESIGN §3.2–3.3; all thresholds read from config/prereg.yaml, none hard-coded)
1. **Orthology.** The OrthoXML in the `--orthology` dir is parsed directly (`sog.read_orthoxml`; no Ancestra env
   needed). Root groups are the OGs, and `og_id` is the root group id. Root clade = TaxRange; for paralogGroup roots
   it is the LCA of the member species. OGs with ≥ 3 species are run (prereg `min_species_per_root_og`).
   - Tune/test: the key is the smallest member accession, sha1 % 2. OGs containing a control-list accession are
     forced into test.
   - O2: pass `--orthology data/orthology/ancestra_O2/<run> --orthology-version O2`. Nothing else changes, except
     the ID-gap mapping hook (open issue 3).
2. **Models.** AFDB v6 mmCIF is read straight from the proteome tars at their `results/webapp_index` offsets.
   - Parsed fields: the model sequence (`pdbx_seq_one_letter_code_can`) and per-residue pLDDT (`_ma_qa_metric_local`,
     equal to the B-factor).
   - pLDDT is cached as `floor(2·pLDDT)` in uint8. This is exact at the integer thresholds (70, 90); the mean has
     ≤ 0.25 error.
   - Fragments (215 human proteins) are combined in UniProt numbering (step 200), taking each position from its most
     central fragment.
   - Every combined sequence's CRC64 equals the model header's UniProt checksum: 552,810 / 552,810, 0 parse errors.
3. **Sequence check** (prereg `model_sequence_identity_check`). Model sequence md5 and length must equal the
   orthology-input proteome FASTA, and model CRC64 must equal the current UniProt (2026_03) `seq_crc64`. Otherwise the
   member is U (`model_sequence_mismatch`) in every SOG: 4,131 members.
4. **MSA and SOGs.** Per OG:
   - FAMSA 2.4.1-45c9b2b (default parameters) on the model sequences.
   - Each site's residues (UniProt numbering) are mapped to their MSA column set.
   - The site × column incidence product finds candidate pairs (the inverted index). Edges join sites of
     *different* members with column Jaccard ≥ τ.
   - Leiden (RBConfiguration = modularity, resolution 1.0, seed 0, iterate to convergence) gives the SOGs. Labels are
     canonicalised (largest first).
   - Footprint = columns present in ≥ 50 % of the SOG's member sites. Insertion columns are allowed. SOGs are built
     for every τ in the grid in one pass (`site_sog_all_tau.parquet`).
5. **States** for one τ (here **τ = 0.2**, chosen on the tune half: `TAU_BENCHMARK.md`, `tau_choice.json`):
   - **P**: the member has a site with Jaccard(site, footprint) ≥ τ, and its sequence is OK.
   - **A** (primary): every prereg A_primary criterion holds:
     - footprint non-gap fraction ≥ 0.8;
     - mean footprint pLDDT ≥ 90 and frac(< 70) ≤ 0.2;
     - length within 0.7–1.3 × the OG median;
     - no footprint column before the member's first or after its last residue (terminal truncation);
     - no UniProt SIGNAL/TRANSMEM overlap (`sptm_checked`);
     - model status ok.
   - `state_s70`: the same rule with the pLDDT ≥ 70 secondary stratum.
   - **U**: anything else. `reason` holds the first failing criterion, and `fail_bits` holds all of them.
   - **G**: a species of the OG's root clade with no member.
   - **W**: hook only (`states.apply_rescue`; see open issue 1).
6. **Interpretation choices** (not specified in prereg; please confirm):
   - (a) A member whose best site overlaps the footprint below τ, or that owns a Leiden-member site of the SOG
     failing the Jaccard rule, is U (`partial_site_overlap`), never A. A pocket is there, so absence is not observed.
   - (b) G is emitted only for species inside the OG's root clade. Species outside it have no row, which the web app
     draws as blank.
   - (c) Members without a parsed atlas model are U (`no_model`).
   - (d) P follows the prereg rule literally (Jaccard to the footprint). It is not the same as Leiden membership;
     `in_sog` records the latter.

## Funnel (τ = 0.2; `funnel.json`)
| step | count |
|---|---|
| root OGs (O1) | 36,736 |
| OGs with ≥ 3 species, run | 25,347 (tune 12,683 / test 12,664; 13 forced-test control OGs) |
| OGs with ≥ 1 atlas site / with a ≥ 2-species SOG | 23,798 / 21,094 |
| OGs, MSA identity < 0.25 (low-identity bin) | 1,633 |
| members in run OGs / with a parsed model / model-sequence mismatch | 417,970 / 382,896 / 4,131 |
| sites in run OGs (residues dropped outside the model: 0) | 1,431,296 |
| SOGs at τ = 0.2 / 0.3 / 0.4 / 0.5 | 260,991 / 372,699 / 486,879 / 613,994 |
| SOGs with ≥ 2 species at τ = 0.2 / 0.3 / 0.4 / 0.5 | 116,583 / 131,010 / 138,820 / 140,962 |
| SOGs at τ = 0.2: single-site / ≥ 3 species with sites / ≥ 3 P members | 136,542 / 79,331 / 88,759 |
| SOG × member rows (τ = 0.2) | 16,843,508 |

**States (τ = 0.2)**: P 1,078,376; A 748,649 (A at the ≥ 70 stratum: 1,101,363); U 12,831,364; G 2,185,119.

U by first reason:
| reason | rows |
|---|---|
| aligned < 0.8 | 6,589,509 |
| partial_site_overlap | 2,970,720 |
| no_model | 1,363,241 |
| length_ratio | 803,482 |
| plddt < 90 | 768,895 |
| model_sequence_mismatch | 179,913 |
| sp_tm_overlap | 84,510 |
| terminal_truncation | 71,025 |
| frac_lt70 > 0.2 | 69 |

**A is provisional.** 581,031 of the 748,649 A rows (78 %) are on models where the atlas filtered out pockets
(`n_pockets_all > n_pockets_kept`). The web export marks them `rescue:pending`. Re-description is required before
any of them is used (open issue 1).

Named fallback counters:
- famsa_to_mafft: 0
- model_seq_mismatch: 4,131
- msa_cache_miss: 0 in the final run; the first run computed all 25,336 MSAs
- PROPEP was not fetched by the annotation agent, so the SP/TM check covers SIGNAL and TRANSMEM only.

## G2 functional slice (5 control families; `slice_G2_sites.csv`, `slice_G2_ogs.csv`)
ESR1 and AR share one OG (the nuclear-receptor family, HOG0025274), so the five accessions give 4 OGs. The CA2 OG
(HOG0035671) falls in the **tune** half: CA2 is not on the prereg control list. The CDK2 OG is in test by hash.
The literature site residues used below are a functional check only.

| family | OG (half) | best atlas site: residues, coverage of known site | SOG: sites / species / footprint | P (species) / A / U / G |
|---|---|---|---|---|
| CDK2 ATP | HOG0027286 (test) | 68 res, 1.00 | 31 / 24 / 29 | 35 (28) / 0 / 6 / 1 |
| CA2 active site | HOG0035671 (tune) | 10 res, **0.23** (BURROW splits the cone) | 101 / 20 / 6 | 73 (17) / 5 / 123 / 20 |
| ESR1 + AR ligand pocket | HOG0025274 (test) | 36 res, 0.95 (ESR1); 26 res, 0.24 (AR) | 147 / 13 / 17 | 125 (13) / 8 / 104 / 0 |
| LDHA NADH/substrate | HOG0035618 (test) | 43 res, 0.53 | 55 / 24 / 26 | 50 (24) / 0 / 7 / 22 |

In all four OGs the known pocket forms one multi-species SOG. The ESR1 and AR ligand-binding pockets fall in the
same SOG, which is correct: they are positionally homologous across NR paralogs. Slice timings (first run, fresh
FAMSA): MSA 0.16–1.7 s per OG, graph + Leiden (4 τ) 0.5–1.0 s.

## Runtime (64-core node shared with other agents, load 60–200; 40 workers)
| stage | wall | CPU | notes |
|---|---|---|---|
| models (tier 1) | 80 s | 2,404 s | 555,846 mmCIF members, 48.8 GB read from the tars, cache 340 MB (`/tmp/pe_big/sogs/<snap>/models`) |
| inputs | 50 s | — | OrthoXML parse, split, DuckDB shard partition (1,024 shards) |
| sogs, first run | 1,159 s | 28,694 s | MSA 25,312 s summed over OGs (FAMSA uses 4–8 threads above 400 / 1,500 sequences); graph 1,864 s |
| sogs, re-run with cached MSAs | 448 s | 2,379 s | wall is set by the largest OG |
| states (per τ) | 42–53 s | 1.6–2.0 k s | 16.8 M (τ 0.2) to 99.7 M (τ 0.5) rows |
| export | 74 s | — | |
| τ benchmark | 554 s | 1 core | |

The largest OG (HOG0034448: 4,223 members, 18,373 sites, 8.8 M candidate edges) took 1,143 s by itself (FAMSA
707 s, graph 382 s), so it dominates the wall time. Scaling: shards are independent and the work is cached per
shard. Tier 3 needs a Foldcomp model reader (hook: `stage_models` refuses non-tar sources).

## Files
- Here (in git, 128 MB):
  - `web/`: web-app schema sogs-web-v1. See the `webapp/SCHEMA.md` §2 update; use `SOG_DIRS=<this dir>/web`.
  - `og_stats.parquet`: per OG, identity, timings, counts.
  - `sogs_all_tau.parquet`: SOG metadata for all τ.
  - `site_sog_all_tau.parquet`: site → SOG for all τ.
  - `funnel.json`, the τ benchmark files, the G2 slice files, `checksums.sha256`.
- Bulk (node-local, regenerable): `/tmp/pe_big/sogs/20261030T085136Z/` holds
  - `models/`: the pLDDT/sequence cache;
  - `O1/inputs/`, `O1/msa/`: MSAs, 106 MB;
  - `O1/sogs/`;
  - `O1/states_t{0.2,0.3,0.4,0.5}/`: full state tables with every diagnostic column.
  `BULK_MANIFEST.sha256` lists them.

```bash
export UV_CACHE_DIR=/tmp/pe_big/uvcache
A="--manifest data/snapshot/20261030T085136Z/manifest.parquet --orthology data/orthology/ancestra_tier1_v0/accurate \
   --orthology-version O1 --prereg config/prereg.yaml --out /tmp/pe_big/sogs --workers 40"
.venv/bin/python scripts/build_sogs.py $A --stage models,inputs,sogs
.venv/bin/python scripts/sog_tau_benchmark.py --sogs /tmp/pe_big/sogs/20261030T085136Z/O1 \
    --results results/sogs/20261030T085136Z/O1          # writes tau_choice.json
.venv/bin/python scripts/build_sogs.py $A --stage states,export      # tau from tau_choice.json
.venv/bin/python scripts/sog_slice_report.py --sogs /tmp/pe_big/sogs/20261030T085136Z/O1 --tau 0.2 \
    --out results/sogs/20261030T085136Z/O1
```
Determinism: two shards (413, 245) were re-run from scratch with fresh FAMSA and gave identical MSAs, SOGs and site
labels. The tests also check that labels do not depend on site order.

## Open issues (for the PI)
1. **W / rescue.**
   - The atlas stores only kept pockets, so W ("atlas-filtered pocket") needs re-description with `--min-score none`
     plus the other 5 rescue settings. That is the rescue agent's job.
   - `states.apply_rescue(states, rescue, k)` takes per-(sog_id, accession) `n_settings_hit/n_settings_run`.
     It is symmetric: A→W at ≥ k settings, P→W at < k.
   - It has not been applied yet, and 78 % of A rows are `rescue:pending`.
2. **The τ choice sits at the grid edge** (0.2). On the test half, τ 0.2 and 0.3 are indistinguishable.
3. **O2 ID gap.** Atlas accessions that are not the canonical gene (`ancestra_O2/idgap_mapping.parquet`) need a
   member → model accession map. Not implemented. O2 context species without tar models (Foldcomp v4) need a model
   reader.
4. **Interpretation choices 6(a)–(d) above.** (a) is conservative: 2.97 M U rows.
5. **Degenerate footprints.** 4,507 multi-species SOGs at τ = 0.2 have a footprint of ≤ 3 columns: Leiden chains
   heterogeneous small sites (e.g. CA2 g1 has a 6-column footprint). This could be flagged or tightened, but any fix
   needs a PI decision.
6. **Not applied here** (prereg: candidates and case members only): the interface rule, the TED/PAE rule, MSA-method
   flip demotion, and metal-only-site U (the hook `sog_flags` exists). FoldMason is not compared.
