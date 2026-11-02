# Ancestra environment (installed 2026-10-30)

Prefix: `P/envs/ancestra` (P = <repo>). 3.0 GB on disk (pyham pulls ete3/Qt).
No conda/mamba/pixi on the node; a static micromamba 2.10.0 was put in `P/envs/bin/micromamba`.

**Home quota warning:** <home> has a ~50 GB quota shared by all agents (`dd` hit
"Disk quota exceeded" at ~47 GB during the first attempt). Per PI: use /tmp/pe_big (UV_CACHE_DIR=/tmp/pe_big/uvcache; Ancestra work in /tmp/pe_big/ancestra). Keep the micromamba package cache, pip cache,
uv cache, numba cache and all Ancestra run directories in /tmp (local NVMe, 1.7 TB free) and copy only final
outputs (hogs.orthoxml.gz, manifest.json) into P.

```bash
P=<repo>
mkdir -p $P/envs/bin
curl -sSL https://micro.mamba.pm/api/micromamba/linux-64/latest -o /tmp/mm.tar.bz2
tar -xjf /tmp/mm.tar.bz2 -C /tmp bin/micromamba && mv /tmp/bin/micromamba $P/envs/bin/
export MAMBA_ROOT_PREFIX=/tmp/pe_big/ancestra/mamba_root          # package cache OFF the home quota
$P/envs/bin/micromamba create -y -p $P/envs/ancestra -f $P/external/ancestra-mvp-ai/environment.yml   # 2m09s
cd $P/external/ancestra-mvp-ai
PIP_CACHE_DIR=/tmp/pe_big/ancestra/pipcache $P/envs/ancestra/bin/python -m pip install -e . --no-deps
```

Pinned versions (from `environment.yml`, all resolved exactly): python 3.11.18, numpy 1.26.4, scipy 1.12.0,
pandas 2.3.3, pyarrow 23.0.0, numba 0.68.0, python-igraph 1.0.0, leidenalg 0.12.0, pyham 1.1.11,
diamond 2.2.8, mmseqs2 18.8cc5c, famsa 2.4.1, FastTree 2.2.0 (double precision), foldseek 10.941cd33.

## Use

```bash
export PATH=$P/envs/ancestra/bin:$PATH NUMBA_CACHE_DIR=/tmp/pe_big/ancestra/numba TMPDIR=/tmp/pe_big/ancestra/tmp
ancestra --help                      # or: python -m ancestra run ...
```
No activation is needed; putting `envs/ancestra/bin` first on PATH is enough (the pipeline calls
diamond/mmseqs/famsa/FastTree from PATH).

## Checks run

```bash
cd $P/external/ancestra-mvp-ai
python -m pytest -q -m "not slow" -p no:cacheprovider   # 147 passed, 6 skipped, 1 deselected, 61 s
bash example/run_example.sh 4                          # 48 s; accurate P/R/F1 0.9976/0.9122/0.9530, fast 0.9979/0.9148/0.9545
```
(README states 148 passed / 5 skipped; the 6 skips here all need benchmark data — not investigated further.)
`pip install -e` writes `ancestra.egg-info` into the repo (git-ignored; `git status` clean).

## Extra data fetched (in /tmp, not in P)
- NCBI taxdump (`nodes.dmp`, `names.dmp`, `merged.dmp`) → /tmp/pe_big/ancestra/taxdump (the atlas copy in
  /cluster/scratch/gmagnusson/pocketcompare_afdb/ncbi has only merged.dmp + rankedlineage.dmp, which
  `ancestra.speciestree` cannot use).
- 46 UniProt current-release one-protein-per-gene reference proteomes for the tier-1 taxids →
  /tmp/pe_big/ancestra/tier1_fa (162 MB gz, 557,396 proteins), from
  `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/<Kingdom>/<UP>/<UP>_<taxid>.fasta.gz`
  (UP ids looked up via `https://rest.uniprot.org/proteomes/search?query=organism_id:<taxid>&format=json`).
