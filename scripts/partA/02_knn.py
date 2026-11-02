#!/usr/bin/env python
"""Part A stage 2: exact inner-product kNN (FAISS GPU flat, k from prereg pairs.knn_k), resumable.
Run with the faiss venv and CUDA_VISIBLE_DEVICES set to the single GPU you may use."""
import json
import time

from _common import parser, resolve, logger, pairs

ap = parser(__doc__)
ap.add_argument("--block", type=int, default=32768)
ap.add_argument("--fp16-storage", action="store_true")
a = resolve(ap.parse_args())
log = logger(a.work / "partA.log")
k = int(a.prereg["pairs"]["knn_k"]["value"])
t = time.time()
st = pairs.knn_exact(a.work / "fp.f16.npy", a.work / "knn", k=k, block=a.block, fp16_storage=a.fp16_storage, log=log)
st["wall_s"] = time.time() - t
log("knn", json.dumps(st))
prev = a.work / "knn" / "I.i32.npy.json"
runs = json.loads(prev.read_text()).get("runs", []) if prev.exists() else []
pairs.write_sidecar(a.work / "knn" / "I.i32.npy", a.snapshot, None, {"stage": "knn", "runs": runs + [st]})
