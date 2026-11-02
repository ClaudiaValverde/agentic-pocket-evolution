# results/sogs/20261030T085136Z/O2_r1: SOGs + states on orthology O2 (+ tier L + atlas tier-3 context taxa)

Development subset, preliminary. prereg sha256 52f94618…9d83; τ = 0.2 frozen (D8; chosen on O1 tune half,
`../O1/tau_choice.json`). Rules: prereg states + D8 (`small_footprint` <= 3 columns, excluded from candidates),
verifier (`empty_footprint`, >= 2 sites with empty footprint, excluded; pLDDT lower bound for A, so 90.0–90.25 counts
as < 90), D11 (A needs the member footprint to span > 10 residues; `--a-min-footprint-residues 10`, U reason
`footprint_residues<=10`). Method and interpretation choices: `../O1/SOGS.md` (confirmed by the PI in D8).
Code: commit in every sidecar (`git_commit`). Web app: `SOG_DIRS=results/sogs/20261030T085136Z/O2_r1/web`.

## Inputs
- Orthology `data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz`: 157 species (46 tier 1 + 72 atlas-tier-3
  context + 39 tier-L context), 57,332 root OGs, 35,093 with >= 3 species (tune 17,499 / test 17,594; 14 forced test).
  The split key (smallest member accession) is recomputed on O2, so some families change half vs O1 (e.g. CDK2 OG: tune).
- Models (per member, lowest tier wins; the 57,222 tier-L accessions that are also in snapshot tier 3 are read once,
  from tier 3 = identical pockets per the tier-L manifest): tier 1 446,028 · tier 2 (Swiss-Prot) 8,027 · tier 3
  (Foldcomp v4) 217,095 · tier L (Foldcomp v4, `data/tierL/tierL_20261030T112000Z`) 449,932 members.
  Foldcomp pLDDT is lossy-discretised by Foldcomp and carries no UniProt checksum header (sequence check = proteome
  FASTA md5 only).
- ID gap (`data/orthology/ancestra_O2/idgap_mapping.parquet`): 30,976 members take the model of another atlas accession
  of the same gene (identical_seq 22,938 · gene_name+seq 3,208 · seq 2,864 · gene_name_only 888; one model per member and
  one member per model). Non-identical mappings fail the model-sequence check -> U. 47,227 members have no model.
- Sites: `/tmp/pe_big/sites/20261030T085136Z/sites/tier{1,2,3,4}` (tier 3: all 7,841 snapshot units, built
  `--no-fingerprints`; tier 4 = L from the tier-L manifest; site rule unchanged, J >= 0.3).

## Counts (τ = 0.2; `funnel.json`)
members 1,168,309 (with model 1,121,082; model-sequence mismatch 12,119) · sites 3,908,023 · OGs with sites 33,590,
low-identity (< 0.25) 5,331 · SOGs 463,014 (multi-species 192,317; >= 3 species 140,413; single-site 256,839) ·
small_footprint 15,751 (multi-species 12,076) · empty_footprint 3,619 · SOG x member rows 66,968,660.
States: P 2,780,379 · A 548,061 (s70: 864,577) · U 40,745,746 · G 22,894,474.
U by first reason: aligned<0.8 20.3 M · partial_site_overlap 9.11 M · no_model 3.32 M · plddt<90 2.62 M ·
length_ratio 2.55 M · footprint_residues<=10 (D11) 1.94 M · model_sequence_mismatch 0.59 M · terminal_truncation 0.22 M ·
sp_tm_overlap 0.08 M. A rows `rescue:pending` (atlas-filtered pockets exist, not yet re-described): 442,783 / 548,061.
P / A by model tier: 1: 1,129,600 / 231,511 · 2: 22,766 / 6,012 · 3: 556,868 / 107,046 · L: 1,071,145 / 203,492.
τ = 0.3 sensitivity (`sensitivity_tau0.3/`): 729,507 SOGs (229,123 multi-species), P 2,983,671, A 1,599,539.

## §8.1 benchmark on O2 (validation only, τ not re-chosen; `tau_benchmark.csv`)
primary subset F1, tune / test: τ 0.2 0.688 / 0.703 [0.668, 0.738]; 0.3 0.676 / 0.695; 0.4 0.603 / 0.669; 0.5 0.560 / 0.578.
71,365 annotated sites, 62 % detected. `tau_tuning_annotations.parquet` = annotations in O2 tune-half OGs.
G2 slice: `slice_G2_*.csv` (ESR1 and AR are separate Gnathostomata OGs in O2).

## Runtime (64-core node, load 180–260 from other jobs; 32 workers)
models (tiers 2/3/L new: 1.14 M entries read, 33.6 GB) 119 s wall / 2,745 CPU-s · inputs 66 s · sogs (fresh FAMSA)
2,882 s wall / 88,832 CPU-s (MSA 78,905 s summed, graph 6,974 s) · states 128 s (τ 0.2), 142 s (τ 0.3) · export 144 s ·
benchmark 1,265 s.
Bulk: `/tmp/pe_big/sogs/20261030T085136Z/O2_r1/` (`BULK_MANIFEST.sha256`).
