#!/usr/bin/env python
"""Part A stage 1: site universe (tier-1/2 rep fingerprints) -> <work>/universe.parquet + fp.f16.npy."""
import json
import time

from _common import parser, resolve, logger, pairs

a = resolve(parser(__doc__).parse_args())
log = logger(a.work / "partA.log")
if (a.work / "universe.parquet").exists() and (a.work / "fp.f16.npy").exists():
    log("universe exists, skipping"); raise SystemExit(0)
t = time.time()
st = pairs.build_universe(a.sites, a.work, tiers=(1, 2))
st["wall_s"] = time.time() - t
log("universe", json.dumps(st, default=str))
pairs.write_sidecar(a.work / "universe.parquet", a.snapshot, None, {"stage": "universe", **st,
                    "sites_dir": str(a.sites)})
