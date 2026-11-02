# results/sites/20261030T085136Z

The site parquet files (2.9 GB) live outside git, in `/tmp/pe_big/sites/20261030T085136Z/`. They are
node-local, so regenerate them with `scripts/build_sites.py`; the commands are in SITES.md.
Verify a copy with `cd /tmp/pe_big/sites/20261030T085136Z && sha256sum -c <P>/results/sites/20261030T085136Z/checksums.sha256`.

Files here:
- `SITES.md`: method, verified fragment rule, statistics, runtime.
- `site_summary_by_tier.csv`, `pockets_per_site_hist.csv`, `sites_per_protein_hist.csv`: statistics.
- `threshold_sensitivity.csv`: other Jaccard thresholds and within-model Jaccard pairs, 40 units per tier.
- `fragment_offsets.tsv`: `_struct_ref_seq` of all 3,251 fragmented tier-1 models, checked against the rule.
- `checksums.sha256`: sha256 of every output parquet, with paths relative to the output directory.
