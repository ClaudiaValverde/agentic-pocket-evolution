#!/usr/bin/env python3
"""Co-folding negative controls on the four positive-control systems.

For each of CDK2+ATP/Mg, CA2+Zn/AZM, PYGL+PLP, yADH1+NAD/Zn it builds and runs:

  <sys>_wt         wild type + cognate ligand(s)                  (reference)
  <sys>_gly        every known-site residue -> Gly  + cognate ligand(s)
  <sys>_phe        every known-site residue -> Phe  + cognate ligand(s)
  <sys>_decoy_cff  wild type + caffeine in place of the cognate organic ligand
  <sys>_apo        wild type, protein only (ligand-free refold, for BURROW re-description)

All runs are `--no-templates`, 3 seeds x 5 samples. WT and the ablated variants share
one MSA by construction (the cached AF3 jackhmmer MSA with the query row replaced), so
the WT-vs-ablation delta is purely sequence-driven and not an MSA-depth artifact.

Usage
  run_controls.py --build --run --summarise          # everything
  run_controls.py --summarise                        # re-score existing outputs
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import site_ablation as SA  # noqa: E402
from cofold_readouts import summarise_job  # noqa: E402

P = Path(__file__).resolve().parents[2]  # repo root
AF3_PY = os.environ.get("AF3_ENV", "/tmp/af3_env") + "/bin/python"
MK = P / "scripts/cofold/make_af3_input.py"
RUN = P / "scripts/cofold/run_af3.sh"
INP = Path("/tmp/pe_big/cofold/inputs/controls")
OUT_ROOT = Path("/tmp/pe_big/cofold/controls")
ABL = Path("/tmp/pe_big/ablation")
RES = P / "results/cofold/controls"
# GPU 3 is reserved for another agent (coordination decision, 2026-10-30): use 0-2 only.
# Systems sharing a GPU run sequentially, because run_af3.sh preallocates 95% of it.
GPU_OF = {"cdk2": 0, "ca2": 0, "pygl": 1, "yadh1": 2}
IONS = {"ZN", "MG", "CA", "NA", "K", "MN", "FE", "CU", "CO", "NI"}
DECOY = "CCD:CFF"          # caffeine: the common property-mismatched decoy


def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, text=True, capture_output=True, **kw)


def decoy_ligands(ligs: list[str]) -> list[str]:
    """Keep structural ions, replace the organic ligand(s) with one caffeine."""
    ions = [l for l in ligs if l.startswith("CCD:") and l.split(":")[1] in IONS]
    return ions + [DECOY]


def build() -> list[dict]:
    INP.mkdir(parents=True, exist_ok=True)
    jobs = []
    for s in sorted(SA.SYSTEMS):
        sj = json.loads((ABL / f"{s}_site.json").read_text())
        acc, ligs = sj["acc"], sj["ligands_for_cofold"]
        d = INP / s
        d.mkdir(parents=True, exist_ok=True)
        specs = [
            ("wt", sj["wt_seq"], ligs, False),
            ("gly", sj["gly_seq"], ligs, True),
            ("phe", sj["phe_seq"], ligs, True),
            ("decoy_cff", sj["wt_seq"], decoy_ligands(ligs), False),
            ("apo", sj["wt_seq"], [], False),
        ]
        for tag, seq, lg, variant in specs:
            name = f"{s}_{tag}"
            out = d / f"{name}.json"
            cmd = [AF3_PY, str(MK), "--acc", acc, "--seq", seq, "--no-uniprot",
                   "--no-templates", "--seeds", "1,2,3", "--name", name, "--out", str(out)]
            for l in lg:
                cmd += ["--ligand", l]
            if variant:
                cmd.append("--reuse-msa-for-variant")
            r = sh(cmd)
            if r.returncode != 0:
                print(f"[build] FAILED {name} rc={r.returncode}\n{r.stderr}", file=sys.stderr)
                continue
            # AF3's --input_dir globs *.json non-recursively, so the .meta.json
            # sidecars must not sit beside the job files.
            meta = out.with_suffix(".meta.json")
            if meta.exists():
                (d / "meta").mkdir(exist_ok=True)
                meta.rename(d / "meta" / meta.name)
            jobs.append({"system": s, "tag": tag, "name": name, "acc": acc,
                         "ligands": lg, "input": str(out), "gpu": GPU_OF[s],
                         "site": sj["site_uniprot_1based"]})
            print(f"[build] {name}: {len(seq)} aa, ligands={lg}, variant_msa={variant}")
    (INP / "jobs.json").write_text(json.dumps(jobs, indent=1))
    return jobs


def run(jobs: list[dict]) -> None:
    """One worker per GPU; systems assigned to the same GPU run one after another."""
    by_gpu: dict[int, list[str]] = {}
    for s, gpu in GPU_OF.items():
        if any((INP / s).glob("*.json")):
            by_gpu.setdefault(gpu, []).append(s)
    procs = {}
    for gpu, systems in sorted(by_gpu.items()):
        script = " && ".join(
            f'bash {RUN} -i {INP / s} -j controls_{s} -g {gpu}' for s in systems)
        logf = open(f"/tmp/pe_big/cofold/controls_g{gpu}.out", "w")
        procs[gpu] = (subprocess.Popen(["bash", "-c", script], stdout=logf,
                                       stderr=subprocess.STDOUT,
                                       env={**os.environ, "AF3_OUT_ROOT": str(OUT_ROOT)}),
                      logf, time.time(), systems)
        print(f"[run] GPU {gpu} <- {systems}")
    for gpu, (p, f, t0, systems) in procs.items():
        rc = p.wait()
        f.close()
        print(f"[run] GPU {gpu} {systems} finished rc={rc} in {time.time()-t0:.0f}s")


def summarise(jobs: list[dict]) -> None:
    RES.mkdir(parents=True, exist_ok=True)
    rows = []
    for j in jobs:
        # make_af3_input appends "_variantmsa" to the job name for the ablated variants
        base = OUT_ROOT / f"controls_{j['system']}"
        jd = next((p for p in sorted(base.glob(f"{j['name']}*")) if p.is_dir()),
                  base / j["name"])
        if not jd.exists():
            print(f"[summarise] missing {jd}", file=sys.stderr)
            continue
        if j["tag"] == "apo":
            continue  # no ligand: handled by the BURROW re-description step
        _, summ = summarise_job(jd, j["site"], pb_top=True)
        if summ is None:
            continue
        summ = {"system": j["system"], "tag": j["tag"],
                "ligands": "+".join(l.split(":")[1] for l in j["ligands"]), **summ}
        rows.append(summ)
        print(f"[summarise] {j['name']}: pLDDT={summ['top_lig_plddt']} "
              f"PAEmin={summ['top_prot_lig_pae_min']} d={summ['top_site_centroid_dist']}")
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(RES / "summary.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join("" if r.get(k) is None else str(r[k]) for k in keys) + "\n")
    # exchange rate: delta of each readout vs the system's own WT + cognate ligand run
    xr = []
    for s in sorted({r["system"] for r in rows}):
        wt = next((r for r in rows if r["system"] == s and r["tag"] == "wt"), None)
        if wt is None:
            continue
        for r in rows:
            if r["system"] != s or r["tag"] == "wt":
                continue
            xr.append({"system": s, "perturbation": r["tag"],
                       "d_lig_plddt": round(r["top_lig_plddt"] - wt["top_lig_plddt"], 1),
                       "d_pae_min": round(r["top_prot_lig_pae_min"] - wt["top_prot_lig_pae_min"], 2),
                       "d_centroid": round(r["top_site_centroid_dist"] - wt["top_site_centroid_dist"], 2),
                       "d_iptm": round(r["top_prot_lig_iptm"] - wt["top_prot_lig_iptm"], 3),
                       "d_site_recall": round(r["top_site_recall"] - wt["top_site_recall"], 2),
                       "wt_pb_valid": wt.get("pb_valid"), "pb_valid": r.get("pb_valid")})
    keys = list(xr[0].keys())
    with open(RES / "exchange_rate.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in xr:
            f.write("\t".join("" if r.get(k) is None else str(r[k]) for k in keys) + "\n")
    print(f"[summarise] wrote {RES}/summary.tsv and exchange_rate.tsv")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--summarise", action="store_true")
    a = ap.parse_args()
    jobs = build() if a.build else json.loads((INP / "jobs.json").read_text())
    if a.run:
        run(jobs)
    if a.summarise:
        summarise(jobs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
