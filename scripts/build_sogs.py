#!/usr/bin/env python
"""Build site orthogroups (SOGs) and observation states (DESIGN.md §3.2-3.3; config/prereg.yaml).

Stages (run in order; each is resumable — finished work units are skipped when their input key is unchanged):
  models   per-residue sequence + pLDDT cache from AFDB v6 mmCIF (tier 1/2 tars via results/webapp_index),
           /tmp/pe_big/sogs/<snapshot>/models/<tar>/<chunk>.parquet   (orthology independent)
  inputs   parse orthology (OrthoXML), split tune/test, shard OGs (1024 shards, shard = crc32(og_id) & 0x3ff,
           webapp bucket = shard & 0xff), write per-shard model + site inputs
  sogs     per shard: FAMSA MSA per OG (cached), site -> column sets, column-Jaccard graph, Leiden for every tau
           in the prereg grid, footprints
  states   per shard: P/A/U(+G) per SOG x member for one tau (--tau, or tau_choice.json)
  export   webapp schema sogs-web-v1 + funnel + manifest into --results

    export UV_CACHE_DIR=/tmp/pe_big/uvcache
    .venv/bin/python scripts/build_sogs.py --manifest data/snapshot/20261030T085136Z/manifest.parquet \\
        --orthology data/orthology/ancestra_tier1_v0/accurate --orthology-version O1 \\
        --prereg config/prereg.yaml --out /tmp/pe_big/sogs --workers 32 --stage models,inputs,sogs
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import glob
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
import zlib
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "src"))
from pocketevo import afmodel, sog, states  # noqa: E402

AFDB_V6 = "/cluster/work/beltrao/shared/alphafold/v6"
FOLDCOMP_DIR = "/cluster/scratch/gmagnusson/pocketcompare_afdb/foldcomp"
ATLAS = "/cluster/scratch/gmagnusson/pocketcompare_afdb"
N_SHARDS = 1024
FAMSA_DEFAULT = str(P / "envs/ancestra/bin/famsa")
MAX_IDENTITY_PAIRS = 200


# ------------------------------------------------------------------------------------------------ provenance
CODE_FILES = ["src/pocketevo/sog.py", "src/pocketevo/states.py", "src/pocketevo/afmodel.py", "scripts/build_sogs.py",
              "scripts/sog_tau_benchmark.py", "scripts/sog_slice_report.py"]


def git_commit() -> str:
    try:
        c = subprocess.run(["git", "-C", str(P), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(P), "status", "--porcelain", "--", *CODE_FILES],
                               capture_output=True, text=True).stdout.strip()
        return c + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def code_sha256() -> str:
    h = hashlib.sha256()
    for f in ["src/pocketevo/sog.py", "src/pocketevo/states.py", "src/pocketevo/afmodel.py", "scripts/build_sogs.py"]:
        h.update((P / f).read_bytes())
    return h.hexdigest()


class Ctx:
    def __init__(self, a):
        self.a = a
        self.prereg = yaml.safe_load(open(a.prereg))
        self.params = sog.prereg_params(self.prereg)
        # PI decision D11 (DESIGN §12 deviation ledger; not a prereg key): A needs member footprint > N residues
        self.params["A_min_footprint_residues"] = a.a_min_footprint_residues
        self.prereg_sha = sog.prereg_sha256(a.prereg)
        self.snapshot = json.load(open(Path(a.manifest).with_suffix(".json")))["snapshot_id"]
        self.seed = int(self.params["seeds"][0])
        self.root = Path(a.out) / self.snapshot
        # every (re)build goes to its own run directory; committed runs are never overwritten
        self.run_id = a.run_id or a.orthology_version
        self.orth_root = self.root / self.run_id
        self.results = Path(a.results) if a.results else P / "results/sogs" / self.snapshot / self.run_id
        if (not a.allow_existing_results and self.results.exists()
                and subprocess.run(["git", "-C", str(P), "ls-files", "--error-unmatch", str(self.results)],
                                   capture_output=True).returncode == 0
                and any(st in a.stage for st in ("export",))):
            raise SystemExit(f"{self.results} is a committed run: use a new --run-id")
        self.git = git_commit()

    def sidecar(self, path, extra=None):
        d = {"snapshot_id": self.snapshot, "prereg_sha256": self.prereg_sha, "git_commit": self.git,
             "seed": self.seed, "code_version": self.git, "created_utc": pd.Timestamp.now("UTC").isoformat()}
        if extra:
            d.update(extra)
        with open(str(path) + ".json", "w") as fh:
            json.dump(d, fh, indent=1, default=str)


def write_parquet(df: pd.DataFrame | pa.Table, path, ctx: Ctx | None = None, extra=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    t = df if isinstance(df, pa.Table) else pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(t, tmp, compression="zstd")
    os.replace(tmp, path)
    if ctx is not None:
        ctx.sidecar(path, extra)


# ------------------------------------------------------------------------------------------------ stage: models
def _model_task(args):
    """One chunk of model files from one source file (AFDB v6 tar member = gzip mmCIF, or Foldcomp v4 entry)."""
    kind, src, tier, rows, out = args
    t0 = time.time()
    recs = []
    by_acc = {}
    if kind == "foldcomp":
        import foldcomp
    path = os.path.join(AFDB_V6, src) if kind == "tar_cif_gz" else os.path.join(FOLDCOMP_DIR, src)
    with open(path, "rb") as fh:
        for acc, frag, off, size in rows:
            fh.seek(off)
            blob = fh.read(size)
            try:
                if kind == "tar_cif_gz":
                    r = afmodel.parse_cif_gz(blob)
                    by_acc.setdefault(acc, {})[int(frag)] = (r["sequence"], r["plddt_q"], r["refdb"])
                else:
                    # Foldcomp keeps B-factors (= pLDDT) with lossy discretisation; no UniProt checksum header
                    d = foldcomp.get_data(blob)
                    seq = d["residues"]
                    b = np.asarray(d["b_factors"], dtype=np.float64)
                    if len(b) != len(seq):
                        raise ValueError(f"b_factors {len(b)} != residues {len(seq)}")
                    q = np.floor(np.clip(b, 0, 100) * afmodel.Q_SCALE).astype(np.uint8)
                    by_acc.setdefault(acc, {})[int(frag)] = (seq, q, {})
            except Exception as e:  # recorded, never silently dropped
                by_acc.setdefault(acc, {})[int(frag)] = None
                recs.append({"accession": acc, "error": f"F{frag}: {type(e).__name__}: {e}"[:300]})
    errs = {r["accession"]: r["error"] for r in recs}
    rows_out = []
    for acc, fr in by_acc.items():
        if any(v is None for v in fr.values()):
            rows_out.append({"accession": acc, "tier": tier, "model_source": kind, "n_res": 0, "sequence": None,
                             "seq_md5": None, "crc64": None, "crc64_ref": None, "crc_ok": False, "n_fragments": len(fr),
                             "fragment_conflict": False, "plddt_q": None,
                             "parse_error": errs.get(acc, "fragment parse error")})
            continue
        seq, q, conflict = afmodel.combine_fragments({k: (v[0], v[1]) for k, v in fr.items()})
        ref = fr[min(fr)][2]
        c = afmodel.crc64(seq)
        cref = ref.get("seq_db_sequence_checksum")
        rows_out.append({"accession": acc, "tier": tier, "model_source": kind, "n_res": len(seq), "sequence": seq,
                         "seq_md5": afmodel.seq_md5(seq), "crc64": c, "crc64_ref": cref,
                         "crc_ok": (c == cref) if cref is not None else True, "n_fragments": len(fr),
                         "fragment_conflict": bool(conflict), "plddt_q": q, "parse_error": None})
    schema = pa.schema([("accession", pa.string()), ("tier", pa.int8()), ("model_source", pa.string()),
                        ("n_res", pa.int32()), ("sequence", pa.string()),
                        ("seq_md5", pa.string()), ("crc64", pa.string()), ("crc64_ref", pa.string()),
                        ("crc_ok", pa.bool_()), ("n_fragments", pa.int16()), ("fragment_conflict", pa.bool_()),
                        ("plddt_q", pa.list_(pa.uint8())), ("parse_error", pa.string())])
    cols = {f.name: [r[f.name] for r in rows_out] for f in schema}
    cols["plddt_q"] = [None if v is None else v.tolist() for v in cols["plddt_q"]]
    t = pa.table(cols, schema=schema)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(t, out + ".tmp", compression="zstd")
    os.replace(out + ".tmp", out)
    return {"src": src, "out": out, "n": len(rows_out), "n_files": len(rows), "bytes": int(sum(r[3] for r in rows)),
            "errors": len(errs), "s": time.time() - t0}


IDX_COLS = ["accession", "fragment", "tier", "blob_kind", "blob_file", "blob_offset", "blob_size"]


def model_index(ctx: Ctx, tiers: list) -> pd.DataFrame:
    """accession/fragment -> model blob, for tiers 1/2 (AFDB v6 tars), 3 (Foldcomp v4, atlas snapshot) and
    4 = tier L (Foldcomp v4 entries listed by the observation agent)."""
    a = ctx.a
    frames = []
    for t in tiers:
        if t in (1, 2, 3):
            dirs = [P / "results/webapp_index" / ctx.snapshot, Path("/tmp/pe_big/webapp_index") / ctx.snapshot]
            files = [f for d in dirs for f in glob.glob(str(d / f"tier{t}" / "*.parquet"))]
            if not files:
                raise SystemExit(f"no model index for tier {t}")
            frames += [pd.read_parquet(f, columns=IDX_COLS) for f in files]
        elif t == 4:
            e = pd.read_parquet(a.tierl_entries, columns=["accession", "fragment", "off", "len"])
            frames.append(pd.DataFrame({"accession": e.accession, "fragment": e.fragment.astype("int16"), "tier": 4,
                                        "blob_kind": "foldcomp", "blob_file": "afdb_uniprot_v4",
                                        "blob_offset": e.off.astype("int64"), "blob_size": e.len.astype("int64")}))
    ix = pd.concat(frames, ignore_index=True)
    ix["tier"] = ix.tier.astype("int8")
    # one tier per accession (lowest tier wins, as in the atlas)
    best = ix.groupby("accession").tier.min()
    ix = ix[ix.tier.to_numpy() == ix.accession.map(best).to_numpy()]
    return ix.drop_duplicates(["accession", "fragment"])


def needed_model_accessions(ctx: Ctx) -> set:
    o = sog.read_orthoxml(os.path.join(ctx.a.orthology, "hogs.orthoxml.gz"))
    need = set(o.members.accession)
    if ctx.a.idmap:
        im = pd.read_parquet(ctx.a.idmap, columns=["atlas_accession", "hog_member_accession"])
        need |= set(im[im.hog_member_accession.isin(need)].atlas_accession)
    return need


def stage_models(ctx: Ctx):
    a = ctx.a
    tiers = [int(x) for x in a.model_tiers.split(",")]
    ix = model_index(ctx, tiers)
    # tiers 1-2 are cached completely (orthology independent); tier 3 / L only for the orthology members
    if any(t >= 3 for t in tiers):
        need = needed_model_accessions(ctx)
        ix = ix[(ix.tier <= 2) | ix.accession.isin(need)]
    out_root = ctx.root / "models"
    tasks = []
    chunk = a.model_chunk
    for (kind, src, tier), g in ix.groupby(["blob_kind", "blob_file", "tier"]):
        g = g.sort_values("blob_offset")
        accs = g.accession.drop_duplicates().tolist()   # all fragments of an accession stay in one chunk
        stem = src.replace(".tar", "") if kind == "tar_cif_gz" else f"foldcomp_{src}_tier{tier}"
        for k in range(0, len(accs), chunk):
            sub = g[g.accession.isin(set(accs[k:k + chunk]))]
            key = hashlib.sha1(("|".join(sub.accession + ":" + sub.fragment.astype(str))).encode()).hexdigest()[:12]
            out = str(out_root / stem / f"chunk{k // chunk:05d}_{key}.parquet")
            if os.path.exists(out):
                continue
            tasks.append((kind, src, int(tier),
                          list(sub[["accession", "fragment", "blob_offset", "blob_size"]].itertuples(index=False, name=None)), out))
    print(f"[models] {len(ix)} model files, {len(tasks)} chunks to do", flush=True)
    t0 = time.time()
    stats = []
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for i, r in enumerate(ex.map(_model_task, tasks, chunksize=1)):
            stats.append(r)
            if (i + 1) % 50 == 0:
                print(f"[models] {i + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    wall = time.time() - t0
    summ = {"stage": "models", "n_model_files": int(len(ix)), "n_chunks_done": len(tasks), "wall_s": wall,
            "cpu_s_sum": float(sum(s["s"] for s in stats)), "bytes_read": int(sum(s["bytes"] for s in stats)),
            "parse_errors": int(sum(s["errors"] for s in stats)), "tiers": tiers, "workers": a.workers,
            "by_tier": ix.groupby("tier").size().to_dict()}
    out_root.mkdir(parents=True, exist_ok=True)
    with open(out_root / f"run_models_{int(time.time())}.json", "w") as fh:
        json.dump(summ, fh, indent=1, default=int)
    ctx.sidecar(out_root / f"models_run_{int(time.time())}", summ)
    print(json.dumps(summ, default=int), flush=True)


# ------------------------------------------------------------------------------------------------ stage: inputs
def shard_of(og_id: str) -> int:
    return zlib.crc32(og_id.encode()) & (N_SHARDS - 1)


def read_proteome_md5(fasta_dir: str, accessions: set) -> pd.DataFrame:
    rows = []
    for f in sorted(glob.glob(os.path.join(fasta_dir, "*.fa*"))):
        name, buf = None, []
        with gzip.open(f, "rt") if f.endswith(".gz") else open(f) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name in accessions:
                        s = "".join(buf)
                        rows.append((name, len(s), hashlib.md5(s.encode()).hexdigest()))
                    h = line[1:].split()[0]
                    name, buf = (h.split("|")[1] if "|" in h else h), []
                else:
                    buf.append(line.strip())
            if name in accessions:
                s = "".join(buf)
                rows.append((name, len(s), hashlib.md5(s.encode()).hexdigest()))
    return pd.DataFrame(rows, columns=["accession", "src_len", "src_md5"]).drop_duplicates("accession")


def stage_inputs(ctx: Ctx):
    import duckdb
    a, p = ctx.a, ctx.params
    t0 = time.time()
    oxml = os.path.join(a.orthology, "hogs.orthoxml.gz")
    o = sog.read_orthoxml(oxml)
    m, og = o.members, o.ogs
    keys = m.groupby("og_id").accession.min().rename("split_key")
    og = og.merge(keys, left_on="og_id", right_index=True)
    og["half"] = og.split_key.map(sog.split_half)
    forced = set(p["forced_test_accessions"])
    forced_og = set(m[m.accession.isin(forced)].og_id)
    og["forced_test"] = og.og_id.isin(forced_og)
    og.loc[og.forced_test, "half"] = "test"
    og["shard"] = og.og_id.map(shard_of).astype("int32")
    og["bucket"] = (og.shard & 0xFF).astype("int32")
    og["eligible"] = og.n_species >= p["min_species"]
    if a.slice_accessions:
        sl = set(m[m.accession.isin(a.slice_accessions)].og_id)
        og["in_slice"] = og.og_id.isin(sl)
    else:
        og["in_slice"] = True
    m = m.merge(og[["og_id", "shard", "eligible", "in_slice"]], on="og_id")
    con = duckdb.connect()
    con.execute(f"SET threads={a.workers}")
    model_glob = str(ctx.root / "models" / "*" / "*.parquet")
    # model cache: one row per accession (lowest tier; tier-1 chunks written before the tier column existed = 1)
    cols = set(con.execute(f"DESCRIBE SELECT * FROM read_parquet('{model_glob}', union_by_name=true)").df().column_name)
    tier_expr = "coalesce(tier, 1)" if "tier" in cols else "1"
    con.execute(f"""CREATE TEMP VIEW cache AS SELECT * EXCLUDE (rn) FROM (
                      SELECT *, {tier_expr} AS model_tier,
                             row_number() OVER (PARTITION BY accession ORDER BY {tier_expr}) AS rn
                      FROM read_parquet('{model_glob}', union_by_name=true)) WHERE rn = 1""")
    cache_acc = con.execute("SELECT accession, model_tier FROM cache").df()
    tier_of = dict(zip(cache_acc.accession, cache_acc.model_tier))
    # member -> model accession: its own model, else the ID-gap map (atlas accession of the same gene)
    own = m.accession.isin(tier_of.keys())
    m["model_accession"] = np.where(own, m.accession, None)
    m["model_map_method"] = np.where(own, "own", None)
    n_idmap = 0
    if a.idmap:
        im = pd.read_parquet(a.idmap)
        prio = {"identical_seq": 0, "gene_name+seq": 1, "seq": 2, "gene_name_only": 3}
        im = im[im.method.isin(prio.keys()) & im.atlas_accession.isin(tier_of.keys())
                & im.hog_member_accession.isin(set(m.accession[~own]))].copy()
        im["prio"] = im.method.map(prio)
        im = im.sort_values(["hog_member_accession", "prio", "n_pockets_kept", "atlas_accession"],
                            ascending=[True, True, False, True]).drop_duplicates("hog_member_accession")
        # one model per member and one member per model (a model accession that is itself a member stays its own)
        im = im[~im.atlas_accession.isin(set(m.accession))].drop_duplicates("atlas_accession")
        mp = dict(zip(im.hog_member_accession, im.atlas_accession))
        mm_ = dict(zip(im.hog_member_accession, im.method))
        fill = ~own & m.accession.isin(mp.keys())
        m.loc[fill, "model_accession"] = m.loc[fill, "accession"].map(mp)
        m.loc[fill, "model_map_method"] = m.loc[fill, "accession"].map(mm_)
        n_idmap = int(fill.sum())
    m["model_tier"] = m.model_accession.map(tier_of).astype("Int8")
    # model status and kept-vs-all pocket counts (for the rescue flag), per model accession
    need = pd.DataFrame({"accession": sorted(set(m.model_accession.dropna()))})
    con.register("need", need)
    st_files = []
    tiers = [int(x) for x in a.model_tiers.split(",")]
    for t in tiers:
        if t in (1, 2):
            st_files += sorted(glob.glob(os.path.join(ATLAS, "structures", f"tier{t}", "*.parquet")))
        elif t == 3:
            man = pd.read_parquet(a.manifest)
            st_files += [os.path.join(ATLAS, x) for x in man[man.tier == 3].structures_path]
        elif t == 4:   # tier L: structure tables listed in the tier-L manifest (paths relative to --tierl-root)
            lm = pd.read_parquet(a.tierl_manifest)
            st_files += [os.path.join(a.tierl_root, x) for x in lm[lm.tier == 4].structures_path]
    fl = ", ".join(f"'{x}'" for x in st_files)
    st = con.execute(f"""SELECT s.accession AS model_accession,
                            CASE WHEN bool_and(s.status = 'ok') THEN 'ok' ELSE string_agg(DISTINCT s.status, ';') END AS status,
                            sum(s.n_pockets_all) AS n_pockets_all, sum(s.n_pockets_kept) AS n_pockets_kept
                         FROM read_parquet([{fl}], union_by_name=true) s JOIN need USING (accession)
                         GROUP BY s.accession""").df()
    m = m.merge(st, on="model_accession", how="left")
    # source sequence check (prereg: model seq == source seq by length + md5)
    src = read_proteome_md5(a.proteome_fasta, set(m.accession)) if a.proteome_fasta and os.path.isdir(a.proteome_fasta) \
        else pd.DataFrame(columns=["accession", "src_len", "src_md5"])
    m = m.merge(src, on="accession", how="left")
    m["src_kind"] = np.where(m.src_md5.notna(), "proteome_fasta", None)
    # current UniProt entry (annotation numbering): model CRC64 must equal UniProt seq_crc64 where annotated
    up = Path(a.annotations) / "uniprot_proteins.parquet"
    if up.exists():
        u = pd.read_parquet(up, columns=["accession", "seq_crc64", "length"]).rename(columns={"length": "uniprot_len"})
        m = m.merge(u.drop_duplicates("accession"), on="accession", how="left")
    else:
        m["seq_crc64"], m["uniprot_len"] = None, np.nan
    inp = ctx.orth_root / "inputs"
    inp.mkdir(parents=True, exist_ok=True)
    write_parquet(m, inp / "members.parquet", ctx, {"orthology": oxml})
    write_parquet(og, inp / "ogs.parquet", ctx, {"orthology": oxml})
    tree = o.newick(label=lambda x: str(sog.species_taxid(x)) if x in set(o.species) else x)
    (inp / "species_tree.nwk").write_text(tree + "\n")
    sp_rows = [(sog.species_taxid(s), s) for s in o.species]
    write_parquet(pd.DataFrame(sp_rows, columns=["taxid", "species"]), inp / "species.parquet", ctx)
    json.dump({"tree_children": o.tree_children, "root": o.root, "species": o.species},
              open(inp / "taxonomy.json", "w"))
    # per-shard model + site inputs (eligible OGs only)
    con.register("mem", m[m.eligible & m.in_slice & m.model_accession.notna()][["og_id", "accession", "model_accession", "shard"]])
    sites_glob = [g_ for t in tiers for g_ in [os.path.join(a.sites, "sites", f"tier{t}", "*.parquet")] if glob.glob(g_)]
    for sub in ["shard_models", "shard_sites"]:
        d = inp / sub
        if d.exists():
            import shutil
            shutil.rmtree(d)
    con.execute(f"""COPY (SELECT m.shard, m.og_id, m.accession, m.model_accession, c.model_tier, c.n_res, c.sequence,
                          c.seq_md5, c.crc64, c.crc_ok, c.fragment_conflict, c.plddt_q, c.parse_error
                   FROM mem m JOIN cache c ON c.accession = m.model_accession ORDER BY m.shard, m.og_id, m.accession)
                   TO '{inp / "shard_models"}' (FORMAT parquet, PARTITION_BY (shard), COMPRESSION zstd)""")
    sg = ", ".join(f"'{x}'" for x in sites_glob)
    # sites of the model accession, from the same tier as the model
    con.execute(f"""COPY (SELECT m.shard, m.og_id, s.site_id, m.accession, m.model_accession, s.site_index, s.residues,
                          s.n_res, s.rep_pocket_id, s.max_burrow_score
                   FROM mem m JOIN cache c ON c.accession = m.model_accession
                   JOIN read_parquet([{sg}], union_by_name=true) s ON s.accession = m.model_accession AND s.tier = c.model_tier
                   ORDER BY m.shard, m.og_id, m.accession, s.site_index)
                   TO '{inp / "shard_sites"}' (FORMAT parquet, PARTITION_BY (shard), COMPRESSION zstd)""")
    n_models = con.execute(f"SELECT count(*) FROM read_parquet('{inp}/shard_models/*/*.parquet')").fetchone()[0]
    n_sites = con.execute(f"SELECT count(*) FROM read_parquet('{inp}/shard_sites/*/*.parquet')").fetchone()[0]
    el = og[og.eligible & og.in_slice]
    summ = {"stage": "inputs", "orthology": oxml, "orthology_version": a.orthology_version,
            "n_ogs": int(len(og)), "n_ogs_eligible": int(og.eligible.sum()), "n_ogs_run": int(len(el)),
            "n_members_run": int(el.n_members.sum()), "n_members_with_model": int(n_models), "n_sites_run": int(n_sites),
            "n_tune": int((el.half == "tune").sum()), "n_test": int((el.half == "test").sum()),
            "n_forced_test_ogs": int(og.forced_test.sum()), "wall_s": time.time() - t0,
            "proteome_fasta": a.proteome_fasta, "n_members_src_seq": int(m.src_md5.notna().sum()),
            "n_members_model_via_idmap": n_idmap, "idmap": a.idmap,
            "members_by_model_tier": {str(k): int(v) for k, v in m[m.eligible].model_tier.value_counts().items()},
            "model_map_methods": {str(k): int(v) for k, v in m[m.eligible].model_map_method.value_counts(dropna=False).items()}}
    json.dump(summ, open(inp / "inputs_run.json", "w"), indent=1)
    ctx.sidecar(inp / "inputs_run", summ)
    print(json.dumps(summ), flush=True)


# ------------------------------------------------------------------------------------------------ stage: sogs
_G = {}


def _shard_inputs(ctx_args, shard):
    inp = Path(ctx_args["inp"])
    md = inp / "shard_models" / f"shard={shard}"
    sd = inp / "shard_sites" / f"shard={shard}"
    models = pd.read_parquet(md) if md.exists() else pd.DataFrame(columns=["og_id", "accession"])
    sites = pd.read_parquet(sd) if sd.exists() else pd.DataFrame(columns=["og_id", "accession", "site_id", "residues"])
    return models, sites


def _member_table(ctx_args, og_id, members_og, models_og):
    """Rows = members with a parsed model (MSA rows), plus flags."""
    mm = members_og.merge(models_og[["accession", "n_res", "sequence", "seq_md5", "crc64", "crc_ok", "fragment_conflict",
                                     "plddt_q", "parse_error"]], on="accession", how="left")
    mm["has_model"] = mm.sequence.notna()
    src_ok = (mm.src_md5.isna()) | ((mm.src_md5 == mm.seq_md5) & (mm.src_len == mm.n_res))
    up_ok = mm.seq_crc64.isna() | (mm.seq_crc64 == mm.crc64) if "seq_crc64" in mm else True
    mm["uniprot_seq_ok"] = up_ok
    mm["seq_ok"] = (mm.has_model & src_ok & up_ok & mm.crc_ok.fillna(False).astype(bool)
                    & ~mm.fragment_conflict.fillna(False).astype(bool))
    mm["status_ok"] = mm.status.fillna("") == "ok"
    return mm.sort_values("accession").reset_index(drop=True)


def _sogs_task(args):
    ctx_args, shard = args
    t0 = time.time()
    p = ctx_args["params"]
    out_dir = Path(ctx_args["orth_root"]) / "sogs"
    out = out_dir / f"shard={shard:04d}.parquet"
    key = ctx_args["key"]
    keyfile = Path(str(out) + ".key")
    if out.exists() and keyfile.exists() and keyfile.read_text() == key:
        return {"shard": shard, "skipped": True}
    if "members" not in _G:
        _G["members"] = pd.read_parquet(Path(ctx_args["inp"]) / "members.parquet")
        _G["ogs"] = pd.read_parquet(Path(ctx_args["inp"]) / "ogs.parquet")
    M, OG = _G["members"], _G["ogs"]
    ogs = OG[(OG.shard == shard) & OG.eligible & OG.in_slice].sort_values("og_id")
    models, sites = _shard_inputs(ctx_args, shard)
    msa_path = Path(ctx_args["orth_root"]) / "msa" / f"shard={shard:04d}.parquet"
    msa_cache = {}
    if msa_path.exists():
        mc = pd.read_parquet(msa_path)
        for og_id, g in mc.groupby("og_id"):
            msa_cache[og_id] = dict(zip(g.accession, g.aligned))
    msa_rows, site_rows, sog_rows, og_rows = [], [], [], []
    taus = p["tau_grid"]
    for og_id in ogs.og_id:
        t_og = time.time()
        mem = _member_table(ctx_args, og_id, M[M.og_id == og_id], models[models.og_id == og_id])
        mw = mem[mem.has_model]
        seqs = dict(zip(mw.accession, mw.sequence))
        rec = {"og_id": og_id, "n_members": len(mem), "n_with_model": len(mw), "n_seq_mismatch": int((mw.seq_ok == False).sum()),  # noqa: E712
               "n_species_with_model": int(mw.taxid.nunique())}
        if len(seqs) == 0:
            rec.update({"status": "no_models"})
            og_rows.append(rec)
            continue
        t1 = time.time()
        cached = msa_cache.get(og_id)
        if cached is not None and set(cached) == set(seqs) and all(cached[k].replace("-", "") == seqs[k] for k in seqs):
            aln = cached
            rec["msa_cached"] = True
        else:
            aln = sog.run_famsa(seqs, ctx_args["famsa"], threads=ctx_args["famsa_threads"](len(seqs)),
                                tmpdir=ctx_args["tmp"])
            rec["msa_cached"] = False
        rec["msa_s"] = time.time() - t1
        msa = sog.msa_from_aligned(aln, order=sorted(seqs))
        msa_rows += [(og_id, k, aln[k]) for k in msa.ids]
        rec["n_cols"] = msa.n_cols
        rec["identity"] = sog.mean_pairwise_identity(msa, MAX_IDENTITY_PAIRS, ctx_args["seed"])
        s = sites[(sites.og_id == og_id)]
        s = s[s.accession.isin(set(msa.ids))].sort_values(["accession", "site_index"])
        row_of = {k: i for i, k in enumerate(msa.ids)}
        sm = s.accession.map(row_of).to_numpy(np.int64)
        t2 = time.time()
        S, dropped = sog.site_columns(msa, sm, list(s.residues))
        res = sog.build_sogs(S, sm, taus, p["leiden_resolution"], ctx_args["seed"], p["footprint_member_fraction"])
        rec["graph_s"] = time.time() - t2
        rec["n_sites"] = int(S.shape[0])
        rec["n_site_res_dropped"] = int(dropped.sum())
        rec["n_pairs_ge_min_tau"] = int(len(res["_edges"][0]))
        taxid_of = dict(zip(mem.accession, mem.taxid))
        site_tax = s.accession.map(taxid_of).to_numpy()
        lab_all = {}
        for t in taus:
            r = res[t]
            lab = r["labels"]
            lab_all[t] = lab
            rec[f"n_sogs_t{t:g}"] = int(lab.max() + 1) if len(lab) else 0
            rec[f"n_edges_t{t:g}"] = r["n_edges"]
            F = r["footprint"]
            for k in range(F.shape[0]):
                idx = np.flatnonzero(lab == k)
                cols = F.indices[F.indptr[k]:F.indptr[k + 1]].astype(np.int32)
                sog_rows.append({"og_id": og_id, "tau": np.float32(t), "sog_index": np.int32(k),
                                 "sog_id": f"{og_id}:g{k}", "n_sites": np.int32(len(idx)),
                                 "n_members": np.int32(len(set(sm[idx]))), "n_species": np.int32(len(set(site_tax[idx]))),
                                 "footprint_columns": cols, "footprint_size": np.int32(len(cols)),
                                 "n_insertion_cols": np.int32((msa.nongap[:, cols].mean(axis=0) < 0.5).sum()) if len(cols) else 0})
        for i_site in range(S.shape[0]):
            row = {"og_id": og_id, "site_id": s.site_id.iloc[i_site], "accession": s.accession.iloc[i_site],
                   "n_cols": np.int32(S.indptr[i_site + 1] - S.indptr[i_site]), "n_res_dropped": np.int32(dropped[i_site])}
            for t in taus:
                row[f"sog_t{t:g}"] = np.int32(lab_all[t][i_site])
            site_rows.append(row)
        rec["status"] = "ok"
        rec["og_s"] = time.time() - t_og
        og_rows.append(rec)
    write_parquet(pd.DataFrame(msa_rows, columns=["og_id", "accession", "aligned"]), msa_path)
    write_parquet(pd.DataFrame(site_rows), out_dir / f"sites_shard={shard:04d}.parquet")
    write_parquet(pd.DataFrame(og_rows), out_dir / f"ogstats_shard={shard:04d}.parquet")
    write_parquet(pd.DataFrame(sog_rows), out)
    keyfile.write_text(key)
    return {"shard": shard, "skipped": False, "n_ogs": len(og_rows), "s": time.time() - t0}


def _famsa_threads(n):
    return 8 if n > 1500 else (4 if n > 400 else 1)


def _ctx_args(ctx: Ctx):
    return {"params": ctx.params, "inp": str(ctx.orth_root / "inputs"), "orth_root": str(ctx.orth_root),
            "famsa": ctx.a.famsa, "famsa_threads": _famsa_threads, "seed": ctx.seed,
            "tmp": str(ctx.root / "tmp"), "prereg_sha": ctx.prereg_sha, "snapshot": ctx.snapshot,
            "key": hashlib.sha256((code_sha256() + ctx.prereg_sha + str(os.path.getmtime(ctx.orth_root / "inputs" / "members.parquet"))).encode()).hexdigest()}


def _run_shards(ctx: Ctx, fn, name):
    og = pd.read_parquet(ctx.orth_root / "inputs" / "ogs.parquet")
    og = og[og.eligible & og.in_slice]
    # largest first (cost ~ n_members^1.5)
    cost = og.groupby("shard").n_members.apply(lambda x: float((x.astype(float) ** 1.5).sum())).sort_values(ascending=False)
    shards = [int(s) for s in cost.index]
    (ctx.root / "tmp").mkdir(parents=True, exist_ok=True)
    ca = _ctx_args(ctx)
    t0 = time.time()
    res = []
    with cf.ProcessPoolExecutor(ctx.a.workers) as ex:
        futs = [ex.submit(fn, (ca, s)) for s in shards]
        for i, f in enumerate(cf.as_completed(futs)):
            res.append(f.result())
            if (i + 1) % 64 == 0 or i + 1 == len(futs):
                print(f"[{name}] {i + 1}/{len(futs)} shards {time.time() - t0:.0f}s", flush=True)
    wall = time.time() - t0
    summ = {"stage": name, "n_shards": len(shards), "n_skipped": int(sum(r["skipped"] for r in res)),
            "wall_s": wall, "cpu_s_sum": float(sum(r.get("s", 0) for r in res)), "workers": ctx.a.workers}
    return summ


def stage_sogs(ctx: Ctx):
    summ = _run_shards(ctx, _sogs_task, "sogs")
    summ["famsa_version"] = sog.famsa_version(ctx.a.famsa)
    summ["tau_grid"] = ctx.params["tau_grid"]
    summ["leiden_resolution"] = ctx.params["leiden_resolution"]
    d = ctx.orth_root / "sogs"
    st = pd.concat([pd.read_parquet(f) for f in glob.glob(str(d / "ogstats_shard=*.parquet"))], ignore_index=True)
    summ["n_ogs"] = int(len(st))
    summ["msa_s_sum"] = float(st.msa_s.sum())
    summ["graph_s_sum"] = float(st.graph_s.sum())
    json.dump(summ, open(d / f"sogs_run_{int(time.time())}.json", "w"), indent=1)
    ctx.sidecar(d / "sogs_run", summ)
    print(json.dumps(summ), flush=True)


# ------------------------------------------------------------------------------------------------ stage: states
def _sptm_masks(ctx_args, accs, n_res):
    ann = _G.get("sptm")
    out = []
    for a_, n in zip(accs, n_res):
        if ann is None or a_ not in ann:
            out.append(None)
            continue
        sig, tm = ann[a_]
        segs = ";".join(x for x in [sig, tm] if isinstance(x, str) and x)
        out.append(states.segments_mask(segs, int(n)))
    return out


def _states_task(args):
    ctx_args, shard = args
    t0 = time.time()
    p, tau = ctx_args["params"], ctx_args["tau"]
    orth = Path(ctx_args["orth_root"])
    out = orth / f"states_t{tau:g}" / f"shard={shard:04d}.parquet"
    key = ctx_args["key"] + f"tau{tau}"
    keyfile = Path(str(out) + ".key")
    if out.exists() and keyfile.exists() and keyfile.read_text() == key:
        return {"shard": shard, "skipped": True}
    if "members" not in _G:
        _G["members"] = pd.read_parquet(Path(ctx_args["inp"]) / "members.parquet")
        _G["ogs"] = pd.read_parquet(Path(ctx_args["inp"]) / "ogs.parquet")
        tx = json.load(open(Path(ctx_args["inp"]) / "taxonomy.json"))
        _G["tax"] = sog.Orthology(members=None, ogs=None, tree_children=tx["tree_children"], species=tx["species"], root=tx["root"])
        ap = ctx_args.get("uniprot_proteins")
        if ap and os.path.exists(ap):
            u = pd.read_parquet(ap, columns=["accession", "signal", "transmem"])
            _G["sptm"] = {r.accession: (r.signal, r.transmem) for r in u.itertuples(index=False)}
    M, OG, TX = _G["members"], _G["ogs"], _G["tax"]
    ogs = OG[(OG.shard == shard) & OG.eligible & OG.in_slice].sort_values("og_id")
    models, sites = _shard_inputs(ctx_args, shard)
    msa_all = pd.read_parquet(orth / "msa" / f"shard={shard:04d}.parquet")
    sg = pd.read_parquet(orth / "sogs" / f"shard={shard:04d}.parquet")
    sg = sg[np.isclose(sg.tau, tau)]
    ss = pd.read_parquet(orth / "sogs" / f"sites_shard={shard:04d}.parquet")
    lab_col = f"sog_t{tau:g}"
    frames = []
    for og_id in ogs.og_id:
        mem = _member_table(ctx_args, og_id, M[M.og_id == og_id], models[models.og_id == og_id])
        g = msa_all[msa_all.og_id == og_id]
        if len(g) == 0:
            continue
        aln = dict(zip(g.accession, g.aligned))
        msa = sog.msa_from_aligned(aln, order=sorted(aln))
        mw = mem[mem.has_model].set_index("accession").loc[msa.ids].reset_index()
        mw["og_median_len"] = float(mw.n_res.median())
        mw["sptm"] = _sptm_masks(ctx_args, mw.accession, mw.n_res)
        s = sites[(sites.og_id == og_id)]
        s = s[s.accession.isin(set(msa.ids))].sort_values(["accession", "site_index"])
        row_of = {k: i for i, k in enumerate(msa.ids)}
        sm = s.accession.map(row_of).to_numpy(np.int64)
        S, _ = sog.site_columns(msa, sm, list(s.residues))
        so = ss[ss.og_id == og_id].set_index("site_id").loc[s.site_id]
        labels = so[lab_col].to_numpy(np.int64)
        gg = sg[sg.og_id == og_id].sort_values("sog_index")
        n_sogs = len(gg)
        if n_sogs == 0:   # OG without any site: no SOG, no state rows
            continue
        import scipy.sparse as sps
        rows = np.repeat(np.arange(n_sogs), gg.footprint_columns.map(len).to_numpy(np.int64))
        cols = np.concatenate(list(gg.footprint_columns)) if n_sogs else np.zeros(0, np.int32)
        F = sps.csr_matrix((np.ones(len(cols), np.int32), (rows, cols)), shape=(n_sogs, msa.n_cols))
        st = states.compute_states(msa, mw, F, S, sm, s.site_id.to_numpy(), labels, tau, p)
        if len(st) == 0:
            continue
        st.insert(0, "og_id", og_id)
        st.insert(1, "sog_id", og_id + ":g" + st.sog_index.astype(str))
        tax_of = dict(zip(mem.accession, mem.taxid))
        st["taxid"] = st.accession.map(tax_of).astype("int32")
        rep = dict(zip(s.site_id, s.rep_pocket_id))
        st["pocket_id"] = st.site_id.map(rep)
        npk = dict(zip(mem.accession, (mem.n_pockets_all.fillna(0) - mem.n_pockets_kept.fillna(0)).astype(int)))
        st["n_filtered_pockets"] = st.accession.map(npk).astype("int32")
        if "model_accession" in mem:
            st["model_accession"] = st.accession.map(dict(zip(mem.accession, mem.model_accession)))
            st["model_tier"] = st.accession.map(dict(zip(mem.accession, mem.model_tier))).astype("Int8")
        extra = []
        # members without a parsed model -> U(no_model), one row per SOG
        nm = mem[~mem.has_model]
        for r in nm.itertuples(index=False):
            reason = "no_model" if pd.isna(r.status) else ("model_parse_error" if isinstance(r.parse_error, str) else "no_model")
            extra.append(pd.DataFrame({"og_id": og_id, "sog_id": gg.sog_id.to_numpy(), "sog_index": gg.sog_index.to_numpy(),
                                       "accession": r.accession, "taxid": np.int32(r.taxid), "state": "U",
                                       "state_s70": "U", "reason": reason,
                                       "fail_bits": np.int16(0)}))
        # G: species of the OG's root clade (TaxRange) without any member
        level = ogs.loc[ogs.og_id == og_id, "level"].iloc[0]
        clade = set(TX.leaves(level)) if level in TX.tree_children or level in TX.species else set()
        have = set(mem.species)
        gsp = sorted(clade - have)
        for spn in gsp:
            extra.append(pd.DataFrame({"og_id": og_id, "sog_id": gg.sog_id.to_numpy(), "sog_index": gg.sog_index.to_numpy(),
                                       "accession": None, "taxid": np.int32(sog.species_taxid(spn)), "state": "G",
                                       "state_s70": "G", "reason": "", "fail_bits": np.int16(0)}))
        frames.append(pd.concat([st] + extra, ignore_index=True) if extra else st)
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    write_parquet(df, out)
    keyfile.write_text(key)
    return {"shard": shard, "skipped": False, "n_rows": len(df), "s": time.time() - t0}


def chosen_tau(ctx: Ctx) -> float:
    if ctx.a.tau is not None:
        return float(ctx.a.tau)
    f = ctx.results / "tau_choice.json"
    if f.exists():
        d = json.load(open(f))
        if d.get("prereg_sha256") != ctx.prereg_sha:
            raise SystemExit("tau_choice.json was made under a different prereg sha256")
        return float(d["tau"])
    raise SystemExit("no --tau and no tau_choice.json (run the benchmark first)")


def stage_states(ctx: Ctx):
    tau = chosen_tau(ctx)
    ca = _ctx_args(ctx)
    ca["tau"] = tau
    ca["uniprot_proteins"] = str(Path(ctx.a.annotations) / "uniprot_proteins.parquet")
    og = pd.read_parquet(ctx.orth_root / "inputs" / "ogs.parquet")
    og = og[og.eligible & og.in_slice]
    cost = og.groupby("shard").n_members.apply(lambda x: float((x.astype(float) ** 1.5).sum())).sort_values(ascending=False)
    t0 = time.time()
    res = []
    with cf.ProcessPoolExecutor(ctx.a.workers) as ex:
        futs = [ex.submit(_states_task, (ca, int(s))) for s in cost.index]
        for i, f in enumerate(cf.as_completed(futs)):
            res.append(f.result())
            if (i + 1) % 64 == 0 or i + 1 == len(futs):
                print(f"[states] {i + 1}/{len(futs)} {time.time() - t0:.0f}s", flush=True)
    summ = {"stage": "states", "tau": tau, "n_shards": len(res), "n_skipped": int(sum(r["skipped"] for r in res)),
            "wall_s": time.time() - t0, "cpu_s_sum": float(sum(r.get("s", 0) for r in res)),
            "n_rows": int(sum(r.get("n_rows", 0) for r in res)), "uniprot_proteins": ca["uniprot_proteins"],
            "sptm_available": os.path.exists(ca["uniprot_proteins"])}
    d = ctx.orth_root / f"states_t{tau:g}"
    json.dump(summ, open(d / f"states_run_{int(time.time())}.json", "w"), indent=1)
    ctx.sidecar(d / "states_run", summ)
    print(json.dumps(summ), flush=True)



# ------------------------------------------------------------------------------------------------ stage: export
def _ref_accession(mem: pd.DataFrame) -> str | None:
    mw = mem[mem.has_model] if "has_model" in mem else mem
    if len(mw) == 0:
        return None
    hs = mw[mw.taxid == 9606]
    return (hs if len(hs) else mw).accession.min()


def stage_export(ctx: Ctx):
    """Compact results (webapp schema sogs-web-v1 + funnel) for one tau into --results."""
    import scipy.sparse  # noqa: F401
    t0 = time.time()
    tau = chosen_tau(ctx)
    orth, res = ctx.orth_root, ctx.results
    res.mkdir(parents=True, exist_ok=True)
    inp = orth / "inputs"
    og = pd.read_parquet(inp / "ogs.parquet")
    M = pd.read_parquet(inp / "members.parquet")
    stats = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(orth / "sogs/ogstats_shard=*.parquet")))], ignore_index=True)
    sg_all = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(orth / "sogs/shard=*.parquet")))], ignore_index=True)
    ss = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(orth / "sogs/sites_shard=*.parquet")))], ignore_index=True)
    sdir = orth / f"states_t{tau:g}"
    sg = sg_all[np.isclose(sg_all.tau, tau)].copy()
    ann_p = Path(ctx.a.annotations) / "uniprot_proteins.parquet"
    gene = pd.read_parquet(ann_p, columns=["accession", "gene_primary"]).set_index("accession").gene_primary if ann_p.exists() else pd.Series(dtype=str)
    # ---- species
    sp = pd.read_parquet(inp / "species.parquet")
    rp_path = Path(ctx.a.orthology).parent / "ref_proteomes.tsv"   # O1 layout: taxid, UP id, kind, name, ...
    rp = pd.read_csv(rp_path, sep="\t", header=None) if rp_path.exists() else None
    names = {}
    if rp is not None:
        names = dict(zip(rp[0].astype(int), rp[3]))
    for cand in sorted(Path(ctx.a.orthology).parent.glob("*_proteomes.tsv")):   # O2 layout: taxid, species columns
        t_ = pd.read_csv(cand, sep="\t")
        if {"taxid", "species"} <= set(t_.columns):
            names.update({int(k): v for k, v in zip(t_.taxid, t_.species) if int(k) not in names})
    sp["name"] = [names.get(int(t), s_) for t, s_ in zip(sp.taxid, sp.species)]
    sp["code"] = None
    sp["clade"] = None
    web = res / "web"
    web.mkdir(parents=True, exist_ok=True)
    write_parquet(sp[["taxid", "name", "code", "clade"]].astype({"taxid": "int32"}), web / "species.parquet", ctx)
    (web / "species_tree.nwk").write_text((inp / "species_tree.nwk").read_text())
    # ---- per-shard states -> counts + webapp buckets
    counts, reason_rows, buckets = [], [], {}
    msa_ref = {}
    for f in sorted(glob.glob(str(sdir / "shard=*.parquet"))):
        st = pd.read_parquet(f)
        if len(st) == 0:
            continue
        shard = int(Path(f).stem.split("=")[1])
        c = st.groupby(["sog_id", "state"]).size().unstack(fill_value=0)
        c2 = st.groupby(["sog_id", "state_s70"]).size().unstack(fill_value=0).add_prefix("s70_")
        counts.append(c.join(c2, how="outer"))
        r = st.groupby(["state", "reason"]).size().rename("n").reset_index()
        reason_rows.append(r)
        web_cols = ["og_id", "sog_id", "taxid", "accession", "state", "site_id", "pocket_id", "jaccard",
                    "footprint_plddt", "reason", "state_s70", "fail_bits", "in_sog", "aligned_frac", "len_ratio",
                    "n_filtered_pockets", "model_accession", "model_tier"]
        w = st[[c_ for c_ in web_cols if c_ in st.columns]].copy()
        w["reason"] = np.where((w.state == "A") & (w.n_filtered_pockets.fillna(0) > 0), "rescue:pending", w.reason.fillna(""))
        buckets.setdefault(shard & 0xFF, []).append(w)
        # ref-numbering of footprints
        msa = pd.read_parquet(orth / "msa" / f"shard={shard:04d}.parquet")
        for og_id, g in msa.groupby("og_id"):
            msa_ref[og_id] = dict(zip(g.accession, g.aligned))
    for b, frames in (buckets.items() if not ctx.a.no_web_states else []):
        w = pd.concat(frames, ignore_index=True).sort_values(["og_id", "sog_id", "taxid", "accession"], na_position="last")
        w = w.astype({"taxid": "int32", "jaccard": "float32", "footprint_plddt": "float32"})
        write_parquet(w, web / "sog_states" / f"bucket={b:02x}.parquet", ctx, {"tau": tau})
    C = pd.concat(counts).fillna(0).astype(int) if counts else pd.DataFrame()
    for s_ in "PWAUG":
        if s_ not in C:
            C[s_] = 0
    # ---- orthogroups + sogs tables
    ran = og[og.eligible & og.in_slice].merge(stats, on=["og_id"], how="left", suffixes=("", "_st"))
    ref_of, label_of = {}, {}
    for og_id, mem in M[M.og_id.isin(set(ran.og_id))].groupby("og_id"):
        mem = mem.assign(has_model=mem.accession.isin(set(msa_ref.get(og_id, {}).keys())) if og_id in msa_ref else False)
        r = _ref_accession(mem)
        ref_of[og_id] = r
        label_of[og_id] = (gene.get(r) if r is not None and r in gene.index else None) or (r or og_id)
    flags = []
    for r in ran.itertuples():
        f = [r.half]
        if r.forced_test:
            f.append("control-family")
        if isinstance(r.identity, float) and not np.isnan(r.identity) and r.identity < ctx.params["low_identity"]:
            f.append("low-identity-msa")
        if r.status != "ok":
            f.append(str(r.status))
        flags.append(";".join(f))
    ran["flags"] = flags
    nsog = sg.groupby("og_id").size()
    ogt = pd.DataFrame({"og_id": ran.og_id, "label": ran.og_id.map(label_of), "n_species": ran.n_species.astype("int32"),
                        "n_members": ran.n_members.astype("int32"), "n_sogs": ran.og_id.map(nsog).fillna(0).astype("int32"),
                        "ref_accession": ran.og_id.map(ref_of), "flags": ran["flags"], "level": ran.level,
                        "half": ran.half, "identity": ran.identity.astype("float32")})
    write_parquet(ogt.sort_values("og_id"), web / "orthogroups.parquet", ctx)
    fp_ref = []
    for og_id, cols in zip(sg.og_id, sg.footprint_columns):
        r = ref_of.get(og_id)
        a = msa_ref.get(og_id, {}).get(r)
        if a is None:
            fp_ref.append(None)
            continue
        p2c = np.flatnonzero(np.frombuffer(a.encode(), np.uint8) != ord("-"))
        pos = np.searchsorted(p2c, cols)
        ok = (pos < len(p2c)) & (p2c[np.minimum(pos, len(p2c) - 1)] == cols)
        fp_ref.append(f"{r}:" + " ".join(str(int(x) + 1) for x in pos[ok]) if ok.any() else None)
    sg["footprint_ref"] = fp_ref
    sg = sg.join(C, on="sog_id")
    sgt = pd.DataFrame({"og_id": sg.og_id, "sog_id": sg.sog_id, "label": None,
                        "footprint_columns": sg.footprint_columns, "footprint_ref": sg.footprint_ref,
                        **{f"n_{s_}": sg[s_].fillna(0).astype("int32") for s_ in "PWAUG"},
                        "event": "none", "event_node": None, "score": np.nan, "case_id": None,
                        "n_sites": sg.n_sites, "n_species_sites": sg.n_species, "footprint_size": sg.footprint_size,
                        "small_footprint": sg.footprint_size <= ctx.a.small_footprint_max,
                        "empty_footprint": (sg.footprint_size == 0) & (sg.n_sites >= 2),
                        "n_insertion_cols": sg.n_insertion_cols,
                        **{f"n_{s_}_s70": sg.get(f"s70_{s_}", pd.Series(0, index=sg.index)).fillna(0).astype("int32") for s_ in "AU"}})
    write_parquet(sgt.sort_values(["og_id", "sog_id"]), web / "sogs.parquet", ctx)
    man = {"schema": "sogs-web-v1", "snapshot": ctx.snapshot, "orthology": ctx.a.orthology_version, "tau": tau,
           "created_utc": pd.Timestamp.now("UTC").isoformat(), "code_sha256": ctx.git, "toy": False,
           "prereg_sha256": ctx.prereg_sha, "git_commit": ctx.git, "seed": ctx.seed,
           "orthology_dir": str(ctx.a.orthology)}
    json.dump(man, open(web / "manifest.json", "w"), indent=1)
    # ---- funnel
    R = pd.concat(reason_rows).groupby(["state", "reason"]).n.sum().reset_index() if reason_rows else pd.DataFrame()
    el = og[og.eligible & og.in_slice]
    multi = sg[sg.n_species >= 2]
    fun = {
        "snapshot_id": ctx.snapshot, "orthology": ctx.a.orthology_version, "tau": tau,
        "ogs_total": int(len(og)), "ogs_ge_min_species": int(og.eligible.sum()), "ogs_run": int(len(el)),
        "ogs_tune": int((el.half == "tune").sum()), "ogs_test": int((el.half == "test").sum()),
        "ogs_with_sites": int((stats.n_sites.fillna(0) > 0).sum()),
        "ogs_with_multispecies_sog": int(multi.og_id.nunique()),
        "ogs_low_identity": int((stats.identity < ctx.params["low_identity"]).sum()),
        "members_run": int(el.n_members.sum()), "members_with_model": int(stats.n_with_model.sum()),
        "members_model_seq_mismatch": int(stats.n_seq_mismatch.sum()),
        "sites_in_run_ogs": int(stats.n_sites.fillna(0).sum()), "site_residues_dropped": int(stats.n_site_res_dropped.fillna(0).sum()),
        "sogs_by_tau": {f"{t:g}": int(np.isclose(sg_all.tau, t).sum()) for t in ctx.params["tau_grid"]},
        "multispecies_sogs_by_tau": {f"{t:g}": int((np.isclose(sg_all.tau, t) & (sg_all.n_species >= 2)).sum()) for t in ctx.params["tau_grid"]},
        "sogs": int(len(sg)), "sogs_singleton_site": int((sg.n_sites == 1).sum()), "sogs_multispecies": int(len(multi)),
        "sogs_ge3_species": int((sg.n_species >= 3).sum()),
        "small_footprint_max": ctx.a.small_footprint_max,
        "sogs_empty_footprint_multisite": int(((sg.footprint_size == 0) & (sg.n_sites >= 2)).sum()),
        "sogs_small_footprint": int((sg.footprint_size <= ctx.a.small_footprint_max).sum()),
        "sogs_small_footprint_multispecies": int(((sg.footprint_size <= ctx.a.small_footprint_max) & (sg.n_species >= 2)).sum()),
        "state_rows": int(R.n.sum()) if len(R) else 0,
        "states": R.groupby("state").n.sum().to_dict() if len(R) else {},
        "states_s70": {k[4:]: int(C[k].sum()) for k in C.columns if k.startswith("s70_")},
        "states_by_reason": {f"{a}|{b}": int(n) for a, b, n in R.itertuples(index=False)} if len(R) else {},
        "fallback_counters": {"famsa_to_mafft": 0, "msa_cache_miss": int((~stats.msa_cached.fillna(True).astype(bool)).sum()),
                              "model_seq_mismatch": int(stats.n_seq_mismatch.sum())},
        "sptm_annotations": ann_p.exists(), "wall_s": time.time() - t0,
    }
    fun = json.loads(json.dumps(fun, default=int))
    json.dump(fun, open(res / "funnel.json", "w"), indent=1)
    ctx.sidecar(res / "funnel.json", {"tau": tau})
    ogstats_out = stats.drop(columns=[c_ for c_ in stats.columns if c_.startswith("msa_cached")])
    write_parquet(ogstats_out, res / "og_stats.parquet", ctx)
    sg_all["small_footprint"] = sg_all.footprint_size <= ctx.a.small_footprint_max
    sg_all["empty_footprint"] = (sg_all.footprint_size == 0) & (sg_all.n_sites >= 2)
    write_parquet(sg_all.drop(columns=["footprint_columns"]), res / "sogs_all_tau.parquet", ctx,
                  {"small_footprint_max": ctx.a.small_footprint_max})
    write_parquet(ss, res / "site_sog_all_tau.parquet", ctx)
    # checksums of everything in results
    lines = []
    for f in sorted(res.rglob("*")):
        if f.is_file() and f.name != "checksums.sha256":
            lines.append(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(res)}")
    (res / "checksums.sha256").write_text("\n".join(lines) + "\n")
    print(json.dumps({k: v for k, v in fun.items() if not isinstance(v, dict)}), flush=True)


# ------------------------------------------------------------------------------------------------ main
STAGES = {"models": stage_models, "inputs": stage_inputs, "sogs": stage_sogs, "states": stage_states,
          "export": stage_export}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--orthology", required=True, help="dir with hogs.orthoxml.gz (O1 or O2)")
    ap.add_argument("--orthology-version", required=True, help="label, e.g. O1 / O2")
    ap.add_argument("--run-id", default=None, help="run directory name (default = orthology version), e.g. O1_r2")
    ap.add_argument("--allow-existing-results", action="store_true")
    ap.add_argument("--prereg", required=True)
    ap.add_argument("--out", default="/tmp/pe_big/sogs")
    ap.add_argument("--results", default=None, help="default results/sogs/<snapshot>/<orthology_version>")
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--stage", default="models,inputs,sogs")
    ap.add_argument("--sites", default=None, help="site dir (default /tmp/pe_big/sites/<snapshot>)")
    ap.add_argument("--model-tiers", default="1")
    ap.add_argument("--model-chunk", type=int, default=2000)
    ap.add_argument("--proteome-fasta", default="/tmp/pe_big/ancestra/tier1_fa",
                    help="orthology input proteomes (source sequences for the model-sequence check)")
    ap.add_argument("--annotations", default=str(P / "data/annotations"))
    ap.add_argument("--famsa", default=FAMSA_DEFAULT)
    ap.add_argument("--slice", default=None, help="comma-separated accessions: run only their OGs")
    ap.add_argument("--tau", type=float, default=None)
    ap.add_argument("--idmap", default=None, help="ID-gap table (atlas_accession, hog_member_accession, method, "
                    "n_pockets_kept): member -> model accession for members without an own atlas model")
    ap.add_argument("--small-footprint-max", type=int, default=3,
                    help="PI decision D8: SOGs with footprint <= this many columns get small_footprint=True "
                         "(excluded from event candidates downstream)")
    ap.add_argument("--a-min-footprint-residues", type=int, default=10,
                    help="PI decision D11: A requires the member's footprint to span > N residues (None-able via -1)")
    ap.add_argument("--no-web-states", action="store_true", help="export: skip web/sog_states (sensitivity runs)")
    ap.add_argument("--tierl-entries", default="/tmp/pe_big/tierL/fc_entries_tierL.parquet")
    ap.add_argument("--tierl-root", default="/tmp/pe_big/tierL", help="tier-L atlas-like root")
    ap.add_argument("--tierl-manifest", default=None, help="tier-L manifest.parquet (data/tierL/<id>/manifest.parquet)")
    a = ap.parse_args()
    a.slice_accessions = [x for x in a.slice.split(",")] if a.slice else None
    ctx = Ctx(a)
    if a.sites is None:
        a.sites = f"/tmp/pe_big/sites/{ctx.snapshot}"
    for st in a.stage.split(","):
        STAGES[st](ctx)


if __name__ == "__main__":
    main()
