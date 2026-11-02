#!/usr/bin/env python
"""D34 (post hoc, exploratory): artifact cascade on ALL funnel-stage specificity pairs (residue gate + shift above the
noise floor), for real data and for every pseudo-clade null seed, so that an empirical null FDR after the cascade can
be computed (mean null passes / real passes, as D24 for gains).

Inputs: rows_<input>.parquet written by `run_specificity.py --rows-dir DIR`. Each input is cascaded by
`run_artifacts.py` in its own scratch results root (candidates.parquet emptied, so no gain/loss candidate file is
read for writing or modified; the committed candidates files are untouched). For null rows the cascade evaluates the
real proteins of the derived/sister clades (the pseudo-clade null permutes member data, not the tree).
Output: <rows-dir>/../cascade_<input>.parquet and, for the real input, the run's specificity_exploratory.parquet
(results/partB/<snap>/<run-id>/specificity/tau<tau>/) rewritten with D34 flags + cascade columns; summary json.

    python scripts/partB/spec_cascade_d34.py --orthology O2 --release O2_r1 --run-id O2_r1_specD34 \
        --rows-dir /tmp/pe_big/partB/specD34/O2_r1/rows
"""
from __future__ import annotations

import argparse
import glob
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orthology", choices=["O1", "O2"], required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--rows-dir", required=True)
    ap.add_argument("--tau", type=float, default=0.2)
    ap.add_argument("--snapshot", default=None)
    a = ap.parse_args()
    snap = a.snapshot or recon.default_snapshot()
    t0 = time.time()
    rows_dir = Path(a.rows_dir)
    work = rows_dir.parent
    src = ROOT / f"results/partB/{snap}/{a.release}"
    summary = {}
    casc_all = {}
    for f in sorted(rows_dir.glob("rows_*.parquet")):
        inp = f.stem.replace("rows_", "")
        df = pd.read_parquet(f)
        fun = df[df["residue_pass"] & df["shift_pass"]].copy()
        rr = work / f"casc_{inp}"
        if rr.exists():
            shutil.rmtree(rr)
        out = rr / snap / a.release
        (out / f"specificity/tau{a.tau}").mkdir(parents=True)
        pd.read_parquet(src / "candidates.parquet").iloc[:0].to_parquet(out / "candidates.parquet", index=False)
        if (src / "rbh_outgroups.parquet").exists():
            shutil.copy(src / "rbh_outgroups.parquet", out / "rbh_outgroups.parquet")
        (rr / "events").symlink_to(ROOT / "results/partB/events")
        fun.assign(tau=a.tau).to_parquet(out / f"specificity/tau{a.tau}/specificity_called.parquet", index=False)
        if len(fun):
            r = subprocess.run([sys.executable, str(ROOT / "scripts/partB/run_artifacts.py"), "--release", a.release,
                                "--orthology", a.orthology, "--snapshot", snap, "--tau", str(a.tau),
                                "--results-root", str(rr)], capture_output=True, text=True)
            if r.returncode:
                raise RuntimeError(f"{inp}: run_artifacts failed\n{r.stderr[-3000:]}")
            c = pd.read_parquet(out / "specificity_candidates.parquet")
        else:
            c = fun.assign(pass_artifacts=pd.Series(dtype=bool), pass_artifacts_strict=pd.Series(dtype=bool))
        c["input"] = inp
        c.to_parquet(work / f"cascade_{inp}.parquet", index=False)
        casc_all[inp] = c
        summary[inp] = dict(funnel_pairs=int(len(fun)), pass_artifacts=int(c["pass_artifacts"].sum()),
                            pass_artifacts_strict=int(c["pass_artifacts_strict"].sum()),
                            tail_normal_called=int(df["called_tail_normal"].sum()),
                            tail_gpd_called=int(df["called_tail_gpd"].sum()),
                            tail_normal_called_pass_artifacts=int((c["called_tail_normal"] & c["pass_artifacts"]).sum()),
                            tail_gpd_called_pass_artifacts=int((c["called_tail_gpd"] & c["pass_artifacts"]).sum()))
        print(f"[{time.time() - t0:.0f}s] {inp}: {summary[inp]}", flush=True)
    nulls = [v for k, v in summary.items() if k.startswith("null")]
    real = summary["real"]
    fdr = {}
    for k in ("funnel_pairs", "pass_artifacts", "pass_artifacts_strict", "tail_normal_called", "tail_gpd_called"):
        nm = sum(v[k] for v in nulls) / len(nulls)
        fdr[k] = dict(real=real[k], null_mean=nm, null_per_seed=[v[k] for v in nulls],
                      empirical_fdr=(min(1.0, nm / real[k]) if real[k] else None))
    # real exploratory table with flags (the run's specificity_exploratory.parquet, D34 version)
    c = casc_all["real"]
    exf = ROOT / f"results/partB/{snap}/{a.run_id}/specificity/tau{a.tau}/specificity_exploratory.parquet"
    c = c.sort_values(["q_bh", "shift"], ascending=[True, False]).copy()
    c["label"] = "exploratory, not significant (prereg BH: 0 calls); D34 post hoc"
    c["flag_bh_prereg_called"] = c["called"].astype(bool)
    c["flag_tail_normal_called"] = c["called_tail_normal"].astype(bool)
    c["flag_tail_gpd_called"] = c["called_tail_gpd"].astype(bool)
    c["flag_pass_artifacts"] = c["pass_artifacts"].astype(bool)
    c["flag_pass_artifacts_strict"] = c["pass_artifacts_strict"].astype(bool)
    c["flag_ge3_species_d"] = c["n_species_d"] >= 3
    c["flag_tail_and_cascade"] = c["flag_tail_normal_called"] & c["flag_pass_artifacts"]
    c.to_parquet(exf, index=False)
    recon.write_sidecar(exf, snap, 0, {"script": "scripts/partB/spec_cascade_d34.py", "argv": sys.argv[1:]})
    res = dict(release=a.release, run_id=a.run_id, orthology=a.orthology, per_input=summary,
               real_vs_null=fdr, runtime_s=time.time() - t0)
    (work / "cascade_summary.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(fdr, indent=1, default=float))


if __name__ == "__main__":
    main()
