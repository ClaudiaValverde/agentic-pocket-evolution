#!/usr/bin/env python
"""Build the web app's accession -> (orthology version, og_id, sog_ids) index from the SOG web exports.

Input : results/sogs/<snapshot>/<version>/web/sog_states/bucket=*.parquet (sogs-web-v1, sharded by og_id)
Output: results/webapp_index/<snapshot>/og/<version>/bucket=<BB>.parquet, BB = f"{crc32(accession) & 0xff:02x}"
        (the same hash as the protein index; sorted by accession). One row per (accession, og_id):
        accession, og_id, sog_k (list<int32>: SOG index k of sog_id `<og_id>:g<k>`), states (string, one letter
        P/W/A/U per entry of sog_k, same order). Members of a SOG with state G have no accession.
Usage : uv run --project . python -m webapp build-og-index [--snapshot S] [--sogs-root DIR] [--out DIR] [--versions O1_r2,O2_r1] [--threads 16]
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import time
import zlib
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[1]


def agg_shard(f):
    t = pq.read_table(f, columns=["og_id", "sog_id", "accession", "state"])
    t = t.filter(pc.is_valid(t["accession"]))
    if not t.num_rows:
        return None
    k = pc.cast(pc.utf8_slice_codeunits(t["sog_id"], 0), pa.string())
    k = pc.cast(pc.replace_substring_regex(k, r"^.*:g", ""), pa.int32())
    t = pa.table({"accession": t["accession"], "og_id": t["og_id"], "k": k, "state": t["state"]})
    g = t.group_by(["accession", "og_id"], use_threads=False).aggregate([("k", "list"), ("state", "list")])
    st = pa.array(["".join(x) for x in g["state_list"].to_pylist()])
    return pa.table({"accession": g["accession"], "og_id": g["og_id"], "sog_k": g["k_list"], "states": st})


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--versions", default="O1_r2,O2_r1")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--sogs-root", default=None, help="dir holding <snapshot>/<version>/web (default results/sogs)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    out = Path(a.out) if a.out else REPO / "results" / "webapp_index" / a.snapshot / "og"
    for v in a.versions.split(","):
        d = (Path(a.sogs_root) if a.sogs_root else REPO / "results" / "sogs") / a.snapshot / v / "web" / "sog_states"
        files = sorted(d.glob("bucket=*.parquet"))
        if not files:
            print(f"skip {v}: no {d}")
            continue
        t0 = time.time()
        with cf.ThreadPoolExecutor(a.threads) as ex:
            parts = [p for p in ex.map(agg_shard, files) if p is not None]
        # a given og_id lives in exactly one shard, so (accession, og_id) is unique across parts
        tab = pa.concat_tables(parts)
        bk = np.array([zlib.crc32(x.encode()) & 0xff for x in tab["accession"].to_pylist()], dtype=np.int16)
        o = out / v
        o.mkdir(parents=True, exist_ok=True)
        for old in o.glob("bucket=*.parquet"):
            old.unlink()
        for b in range(256):
            sub = tab.filter(pa.array(bk == b))
            if sub.num_rows:
                sub = sub.take(pc.sort_indices(sub, sort_keys=[("accession", "ascending"), ("og_id", "ascending")]))
                pq.write_table(sub, o / f"bucket={b:02x}.parquet", compression="zstd", row_group_size=65536)
        (o / "index.json").write_text(json.dumps({"version": v, "rows": tab.num_rows, "shards": 256,
                                                  "seconds": round(time.time() - t0, 1)}) + "\n")
        print(f"{v}: {tab.num_rows} (accession, og) rows in {time.time() - t0:.0f} s -> {o}")


if __name__ == "__main__":
    main()
