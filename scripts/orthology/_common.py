"""Shared helpers for scripts/orthology/*.py (config, paths, downloads, atlas accession table)."""
from __future__ import annotations

import gzip
import hashlib
import os
import sys
import time
import urllib.request

import yaml

P = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(P, "src"))
ANCESTRA_BIN = os.path.join(P, "envs", "ancestra", "bin")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, file=sys.stderr, flush=True)


def load_config(path: str | None) -> dict:
    path = path or os.path.join(P, "config", "orthology.yaml")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    for k in ("snapshot_manifest", "tier1_proteomes", "out_dir"):
        if k in cfg and not os.path.isabs(cfg[k]):
            cfg[k] = os.path.join(P, cfg[k])
    cfg["_path"] = path
    return cfg


def snapshot_id(manifest: str) -> str:
    return os.path.basename(os.path.dirname(os.path.abspath(manifest)))


def fetch(url: str, dest: str, retries: int = 4, timeout: int = 300) -> str:
    """Download url to dest (cached: skipped if dest exists)."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r, open(tmp, "wb") as fh:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    fh.write(b)
            os.replace(tmp, dest)
            return dest
        except Exception as e:  # noqa: BLE001
            log(f"fetch failed ({i+1}/{retries}) {url}: {e}")
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"could not fetch {url}")


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def fasta_accessions(path: str) -> list[str]:
    """UniProt FASTA headers (>sp|ACC|NAME / >tr|ACC|NAME) -> accessions."""
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].split()[0]
                out.append(h.split("|")[1] if "|" in h else h)
    return out


def read_fasta(path: str):
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def acc_of_header(h: str) -> str:
    x = h.split()[0]
    return x.split("|")[1] if "|" in x else x


def atlas_accessions(cfg: dict, cache_dir: str | None = None) -> str:
    """Parquet (tier, accession, taxid, n_pockets_kept) of every structure row in the snapshot manifest.
    Cached in work_dir keyed by the snapshot id. Returns the parquet path."""
    import duckdb
    import pandas as pd
    sid = snapshot_id(cfg["snapshot_manifest"])
    cache_dir = cache_dir or cfg["work_dir"]
    os.makedirs(cache_dir, exist_ok=True)
    out = os.path.join(cache_dir, f"atlas_acc_{sid}.parquet")
    if os.path.exists(out):
        return out
    m = pd.read_parquet(cfg["snapshot_manifest"])
    files = [os.path.join(cfg["atlas_root"], p) for p in m.structures_path]
    con = duckdb.connect()
    con.execute(f"""COPY (SELECT tier::TINYINT tier, accession, taxid, n_pockets_kept
                         FROM read_parquet({files!r}, union_by_name=true)) TO '{out}' (FORMAT parquet)""")
    log(f"atlas accession table: {out}")
    return out
