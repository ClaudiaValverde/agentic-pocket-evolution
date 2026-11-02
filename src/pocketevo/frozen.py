"""Snapshot-invariant inputs (frozen decisions and species-level data), resolved independently of the running snapshot.

Some inputs were produced once, on the calibration snapshot (config/prereg.yaml meta.snapshot_id), and are frozen by
a ledger decision or are species-level: they must NOT be looked up under the snapshot being (re)run, otherwise every
new atlas snapshot fails on a missing file (G5) or, worse, silently skips the step. Resolve them here instead:

    from pocketevo import frozen
    p = frozen.path("rescue_k")            # results/observation/<calibration snapshot>/rescue_k.json
    prov = frozen.provenance("rescue_k")   # {name, path, calibration_snapshot, decision, sha256} for sidecars

Override for tests or a deliberate recalibration: $PE_FROZEN_SNAPSHOT (calibration snapshot id) or
$PE_FROZEN_<NAME> (full path for one input). Nothing here is re-tuned; this only says where the frozen value lives.
"""
from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PREREG = REPO / "config" / "prereg.yaml"

# name -> (path template relative to the repo, {cal} = calibration snapshot; why it is snapshot-invariant)
FROZEN = {
    "rescue_k": ("results/observation/{cal}/rescue_k.json",
                 "D9: rescue k = 1 (k_min_false_A), calibrated once on tiers 1-2 (TUNE half); not re-tuned per snapshot"),
    "error_rates": ("results/observation/{cal}/error_rates.parquet",
                    "observation-layer error rates (tiers 1-2 calibration), used by the descent pair builder"),
    "redetection_sites": ("results/observation/{cal}/redetection_sites.parquet",
                          "observation-layer redetection calibration (tiers 1-2)"),
    "timetree_O1": ("data/timetree/{cal}_O1_pairwise.csv",
                    "species-level TimeTree divergence times for the O1 species set (no atlas dependency)"),
}


def calibration_snapshot() -> str:
    """Snapshot on which the frozen inputs were produced (prereg meta.snapshot_id unless $PE_FROZEN_SNAPSHOT)."""
    if os.environ.get("PE_FROZEN_SNAPSHOT"):
        return os.environ["PE_FROZEN_SNAPSHOT"]
    # plain-text read: avoids a yaml dependency for this one key (Part A runs in a separate env)
    m = re.search(r'^\s*snapshot_id:\s*"?([0-9TZ]+)"?', PREREG.read_text(), re.M)
    if not m:
        raise RuntimeError(f"meta.snapshot_id not found in {PREREG}")
    return m.group(1)


def path(name: str) -> Path:
    if name not in FROZEN:
        raise KeyError(f"unknown frozen input {name!r} (known: {', '.join(FROZEN)})")
    env = os.environ.get(f"PE_FROZEN_{name.upper()}")
    if env:
        return Path(env)
    return REPO / FROZEN[name][0].format(cal=calibration_snapshot())


def provenance(name: str) -> dict:
    p = path(name)
    sha = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
    return {"name": name, "path": str(p.relative_to(REPO) if p.is_relative_to(REPO) else p),
            "calibration_snapshot": calibration_snapshot(), "decision": FROZEN[name][1], "sha256": sha}


def all_paths() -> dict:
    return {n: path(n) for n in FROZEN}
