#!/usr/bin/env python
"""O2 states for Part B: O2_r1 states (D8/D11 already applied) + the rescue overlay O2_r1_rescue/rescue_A.parquet.

- A in the overlay -> state_rescued (W if >= 1 of 6 settings hit, k = 1; else A, `rescue:negative`).
- A `rescue:pending` (n_filtered_pockets > 0) outside the overlay -> U (`rescue:pending_not_run`; PI instruction).
- A with no atlas-filtered pockets (not pending) stays A (not re-described on O2; counted).
P calls are not rescued on O2 (the overlay rescues A only; counted, stated in the report).
Writes /tmp/pe_big/partB/<snap>/O2_r1/states_resc_t0.2/shard=*.parquet (same columns) + a counts json.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg


def main():
    ap = argparse.ArgumentParser()
    recon.add_partb_args(ap, orthology="O2", sog_run="O2_r1")
    ap.add_argument("--overlays", nargs="*", default=None, help="rescue overlays (default <sog-run>_rescue, <sog-run>_rescue2)")
    a = ap.parse_args()
    R_ = a.sog_run
    t0 = time.time()
    src = Path(f"/tmp/pe_big/sogs/{SNAP}/{R_}/states_t0.2")
    out = Path(f"/tmp/pe_big/partB/{SNAP}/{R_}/states_resc_t0.2")
    out.mkdir(parents=True, exist_ok=True)
    if a.overlays:
        ovs = [Path(x) if Path(x).is_absolute() else ROOT / x for x in a.overlays]
    else:
        ovs = [ROOT / f"results/sogs/{SNAP}/{R_}_rescue/rescue_A.parquet"]
        r2 = ROOT / f"results/sogs/{SNAP}/{R_}_rescue2/rescue_A.parquet"     # targeted overlay, if present
        if r2.exists():
            ovs.append(r2)
    ov = pd.concat([pd.read_parquet(x, columns=["sog_id", "accession", "state_rescued", "reason_rescued"]) for x in ovs],
                   ignore_index=True).drop_duplicates(["sog_id", "accession"], keep="last")
    ovm = {(s, a): (st, r) for s, a, st, r in zip(ov["sog_id"], ov["accession"], ov["state_rescued"],
                                                   ov["reason_rescued"])}
    cnt = {"A_in": 0, "A_overlay_W": 0, "A_overlay_A": 0, "A_pending_to_U": 0, "A_not_pending_kept": 0}
    for f in sorted(glob.glob(str(src / "shard=*.parquet"))):
        d = pd.read_parquet(f)
        isA = d["state"] == "A"
        cnt["A_in"] += int(isA.sum())
        for i in d.index[isA]:
            k = (d.at[i, "sog_id"], d.at[i, "accession"])
            if k in ovm:
                st, r = ovm[k]
                d.at[i, "state"] = st
                d.at[i, "reason"] = r
                cnt["A_overlay_W" if st == "W" else "A_overlay_A"] += 1
            elif (d.at[i, "n_filtered_pockets"] or 0) > 0:
                d.at[i, "state"] = "U"
                d.at[i, "reason"] = "rescue:pending_not_run"
                cnt["A_pending_to_U"] += 1
            else:
                cnt["A_not_pending_kept"] += 1
        d.to_parquet(out / Path(f).name, index=False)
    cnt.update(overlays=[str(x.relative_to(ROOT)) for x in ovs], runtime_s=time.time() - t0, out=str(out), P_rescue="not applied on O2 (overlay rescues A only)")
    fo = ROOT / f"results/partB/{SNAP}/{R_}/o2_states_prep.json"
    fo.parent.mkdir(parents=True, exist_ok=True)
    fo.write_text(json.dumps(cnt, indent=1))
    recon.write_sidecar(fo, SNAP, 0, {"script": "scripts/partB/prep_o2_states.py"})
    print(json.dumps(cnt, indent=1))


if __name__ == "__main__":
    main()
