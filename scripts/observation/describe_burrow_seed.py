#!/usr/bin/env python
"""`pocketcompare describe --af-mode` with the pinned BURROW model swapped for another training seed
(burrow_all_default_s1 / _s2 from external/burrow/models). Everything else (merging, sph8, N1, heads) is the
frozen method, so the fingerprints measure "same site, different BURROW seed" (technical noise floor).
The output is NOT an atlas-method table: its schema metadata gets `pocketevo_burrow_model` = the model file and
sha256, and the method record is kept only to show the head/config are the frozen ones.

Runs in-process (--workers 1) so the patch is in effect; parallelise by running several processes on input chunks.

  python scripts/observation/describe_burrow_seed.py --seed s1 -o out.parquet model1.cif.gz model2.cif.gz ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

P = Path(__file__).resolve().parents[2]
MODELS = P / "external/burrow/models"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", choices=["s0", "s1", "s2"], required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("inputs", nargs="+")
    a = ap.parse_args()
    mp = MODELS / f"burrow_all_default_{a.seed}.json"
    import pocketcompare.method as M
    import pocketcompare.pockets as PK
    import pocketcompare.pockets_fast as PF
    M.default()                                  # build + verify the frozen method with the pinned s0 model first
    PK.model_path = lambda: mp                   # then swap the BURROW model used for pocket detection
    PF.model_path = lambda: mp
    PF.fast_model.cache_clear()
    from pocketcompare import cli
    sys.argv = ["pocketcompare", "describe", *a.inputs, "--af-mode", "--backend", "fast", "--workers", "1",
                "-o", a.out]
    rc = cli.main()
    import pyarrow.parquet as pq
    t = pq.read_table(a.out)
    md = dict(t.schema.metadata or {})
    md[b"pocketevo_burrow_model"] = json.dumps({"model": mp.name, "sha256": hashlib.sha256(mp.read_bytes()).hexdigest()}).encode()
    pq.write_table(t.replace_schema_metadata(md), a.out)
    return rc or 0


if __name__ == "__main__":
    sys.exit(main())
