# results/sogs/20261030T085136Z/O1_r3: O1_r2 + D11 (states and export only)

Same SOGs, MSAs and inputs as `../O1_r2` (bulk dirs `inputs`, `msa`, `sogs` under
`/tmp/pe_big/sogs/20261030T085136Z/O1_r3/` are symlinks to O1_r2). τ = 0.2 frozen (D8). Development subset, preliminary.
Adds PI decision D11: A needs the member's footprint to span > 10 residues (U reason `footprint_residues<=10`).
τ = 0.2: P 1,078,376 (unchanged) · A 188,247 (O1_r2: 735,962; 547,715 rows -> U by D11 alone) · A at the >= 70
stratum 296,108 · U 13,391,766 · G 2,185,119. small_footprint / empty_footprint flags as in O1_r2.
τ = 0.3 sensitivity (`sensitivity_tau0.3/`): A 462,355, P 1,163,181. G2 slice: `slice_G2_*.csv`.
Web app: `SOG_DIRS=results/sogs/20261030T085136Z/O1_r3/web`. Bulk: `BULK_MANIFEST.sha256`.
