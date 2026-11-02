# ColabFold MSA provenance

One `<label>_<seq_sha256[:12]>.provenance.json` per MSA fetched from the ColabFold
MMseqs2 server by `scripts/cofold/colabfold_msa.py`: server, mode, ticket id, UTC fetch
time, the returned database file names (the server exposes no version endpoint, so those
names are its DB version strings), sequence counts and a depth proxy. The same rows are
appended to `DATA.md`.

**The gzipped A3M payloads are not in the repo.** The home quota is tight, and each is
2-3 MB. They live at `/tmp/pe_big/msas/<seq_sha256>/unpaired.a3m.gz` (node-local, lost
with the job) together with the raw server tar. `MANIFEST.sha256` records their
checksums, so a refetch can be verified:

```bash
/tmp/af3_env/bin/python scripts/cofold/colabfold_msa.py --acc <ACC>   # refetch, cached
/tmp/af3_env/bin/python scripts/cofold/colabfold_msa.py --export-compact --record-data-md
( cd data/msas && sha256sum -c MANIFEST.sha256 )   # needs the payloads staged locally
```

An MSA is a deterministic-enough function of the sequence plus the server databases, and
the ticket id pins which database snapshot answered, so this is provenance, not data loss.
