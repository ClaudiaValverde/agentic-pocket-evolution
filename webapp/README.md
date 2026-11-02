# pocket-evolution research browser

A small FastAPI app with server-rendered HTML and 3Dmol.js (vendored, so it works offline). It reads index files
and the atlas directly. It never copies AFDB. Structures come from the AFDB v6 tars (tiers 1–2, read by member
offset) or from the Foldcomp DB (tier 3, read by entry offset, converted to mmCIF with gemmi). They are cached only
in an in-memory LRU.

## Run (from the repo root)
```bash
uv run --project . python -m webapp --port 8765
# then open http://127.0.0.1:8765/  (on Euler: ssh -L 8765:localhost:8765 <node>)
```
This is the only command needed. It starts in about 2 s from a fresh shell and prints the index tiers and the orthology
versions it found. Options: `--host`, `--port`, `--snapshot`, `--index-dirs a:b`. Dependencies are pinned in the optional extra
`webapp` of pyproject.toml (`webapp/requirements.txt`); they are already in `.venv`.

Environment variables (all optional; they override `webapp/config.py`):

| variable | default |
|---|---|
| `PE_WEBAPP_SNAPSHOT` | `20261030T085136Z` |
| `PE_WEBAPP_ATLAS` | `atlas_root` of `data/snapshot/<snapshot>/manifest.json` (read in place, never copied) |
| `PE_WEBAPP_AFDB_V6` | `/cluster/work/beltrao/shared/alphafold/v6` (tier 1-2 structures, read by tar member offset) |
| `PE_WEBAPP_FOLDCOMP` | `<atlas>/foldcomp` (tier 3) |
| `PE_WEBAPP_INDEX_DIRS` | `results/webapp_index/<snap>:/tmp/pe_big/webapp_index/<snap>` |
| `PE_WEBAPP_SOG_VERSIONS` | `O1_r2,O2_r1` (looked up as `results/sogs/<snap>/<v>/web`; absent ones are skipped) |
| `PE_WEBAPP_SOG_DIRS` | colon-separated `web/` dirs; replaces the version list and the toy fixture |
| `PE_WEBAPP_CASE_DIRS` | `cases:webapp/fixtures/cases` |

Orthology versions: the selector in the header switches between O1_r2 (25,347 OGs / 260,991 SOGs) and O2_r1 (35,093 OGs /
463,014 SOGs). The choice is kept in a cookie and can be forced with `?v=`. The TOY fixture is a third entry.

## Protein index (build once per snapshot; tier 3 must be rebuilt after a node restart)
One command regenerates all three tiers (about 2 min, tiers 1-2 -> `results/webapp_index/<snap>/` 29 MB, tier 3 ->
`/tmp/pe_big/webapp_index/<snap>/` 458 MB):
```bash
uv run --project . python -m webapp build-index            # = scripts/webapp_build_index.py --tiers 1,2,3
uv run --project . python -m webapp build-index --tiers 3  # tier 3 only (extra flags go to the script)
```
Measured 2026-10-30 (snapshot 20261030T085136Z, 48 threads): tier 1 = 555,846 structures, tier 2 = 436,886,
tier 3 = 22,962,643 (7,841 units). 0 structures lacked a tar member or Foldcomp entry. The index has 256 shards per
tier, keyed by a hash prefix of the accession, and each shard is sorted with 64k-row groups. One lookup reads one
row group per tier (about 10–30 ms when warm). The full tier 3 (214 M entries) would give about 0.8 M rows per
shard. The tier-3 index lives in node-local /tmp, so rebuild it after a job restart. Until then the app serves
tiers 1–2 only and the tier-3 test skips.

Accession -> orthogroup index (for the "Orthogroups and site orthogroups" table of the protein page; 40 MB, 15 s, sharded by
the same accession hash, `results/webapp_index/<snap>/og/<version>/bucket=XX.parquet`, regenerable, not committed):
```bash
uv run --project . python -m webapp build-og-index          # scripts/webapp_build_og_index.py; rerun after a new SOG export
```

## Pages and endpoints
| route | what |
|---|---|
| `/` | search, index status, first orthogroups of the selected version, case list |
| `/search?q=` | accession (also `AF-<acc>-F1-model_v6`), else gene name (exact, then prefix), else protein-name substring (from `data/annotations/uniprot_proteins.parquet`) |
| `/protein/{acc}[?tier=&fragment=&pocket=]` | atlas structures of any tier, pocket table, 3D viewer: click a pocket (or open with `?pocket=<pocket_id>`) to highlight its lining residues; the best pocket is preselected |
| `/api/protein/{acc}`, `/api/structure/{acc}` | JSON (structures + pockets), mmCIF text |
| `/og?q=&flags=&half=&sort=&desc=&page=` | orthogroups of the selected version, paginated (50/page) |
| `/sog?q=&event=&case=1&min_score=&sort=&page=` | all SOGs of the version (463k for O2_r1), filtered, sorted, paginated |
| `/og/{og_id}[?sp=&n=&ev=!none&sog=&full=1]` | SOG x species matrix on the species tree pruned to the species with rows (`full=1`: whole tree). 30 SOG columns per page; `sog=` focuses one SOG and lists every member (state, site, pocket link, Jaccard, pLDDT, reason) |
| `/api/{versions,ogs,sogs}`, `/api/og/{og_id}?sp=&n=&sog=` | paginated JSON |
| `/case`, `/case/{id}`, `/api/case/{id}`, `/case/{id}/file/{path}` | case view (see below) |
| `/api/status`, `/api/docs` | data sources; OpenAPI docs |

**Large data.** The app never loads a whole table. `orthogroups.parquet` (under 1 MB) and a column subset of `sogs.parquet`
(no footprint lists; 41 MB, 0.25 s once per version, kept in memory) serve the lists and filters. `sog_states` (about 40 M
rows) is read per request from the one shard `bucket=crc32(og_id)`, filtered to the OG and to the SOG ids of the page
(an OG can have 10^6 rows); results sit in a 64-entry LRU. The protein index reads one row group per tier; structures
are read by byte offset from the AFDB tar / Foldcomp DB and kept in a 64-entry in-memory LRU.

**Case view** (`cases/<case_id>/case.yaml`, schema `case-v1`, picked up on the next request): summary, proteins, evidence
filters with the evidence linked to each filter, artifact checks, evidence list (PDB/UniProt/ChEMBL/DOI/PMID/M-CSA
links, cached records), co-folding table, **ancestral | derived side-by-side 3D panels** (model selector per panel; ligand
as green sticks, residues within 4.5 Å of the ligand in cyan, case pocket residues in orange, pLDDT colouring), figures,
literature. Paths of `cofold[].file` are relative to the case dir (or declared repo-relative under `results/`).
See SCHEMA.md §3-4.

Input schemas for orthogroups and cases: **`webapp/SCHEMA.md`**. The toy fixtures in `webapp/fixtures/` (rebuilt by
`webapp/fixtures/make_fixtures.py`) are labelled TOY.

## Tests
`.venv/bin/python -m pytest -q tests/test_webapp.py -s` (about 12 s; real-data tests skip if `results/sogs/...`, `cases/` or the
tier-3 index are absent; the last test prints p50 latencies). Warm p50 on 2026-10-30 (TestClient, in-process): protein API 5 ms,
structure 4 ms (tiers 1 and 3), protein page 5 ms, SOG list page 5-23 ms, OG list 7 ms, OG matrix page of the
largest OG (582 SOGs, 1 M member rows) 65 ms. Cold: first protein call 3 s (gene table, warmed in a thread at start-up),
then about 40 ms; first OG matrix of a version about 0.3-0.5 s.

## Caveats
- Pocket residues use the model's author numbering (chain A). The table shows UniProt positions
  (fragment k: +200·(k−1)), and the viewer selects on model numbering.
- Tier-3 models are AFDB **v4** (Foldcomp; coordinates within about 0.1 Å after lossy compression). Tiers 1–2 are v6.
- The protein page's OG table needs the og index above (silently absent otherwise).
- Members with `model_tier` 4 (not in the atlas) are shown without a link.
