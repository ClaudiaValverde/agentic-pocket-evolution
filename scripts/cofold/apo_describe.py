#!/usr/bin/env python3
"""Ligand-free refold control: does BURROW find the known site in the apo AF3 model?

Takes the `<sys>_apo` AF3 jobs produced by run_controls.py, writes the top-ranked
model (and optionally every sample) as plain mmCIF, runs
`pocketcompare describe <cif> --af-mode`, and scores every predicted pocket against
the known holo site (residues within 4.5 A of the crystal ligand, from site_ablation.py)
by residue recall / precision / Jaccard.

This separates "the AFDB model missed the pocket" from "the sequence does not encode
the pocket" (DESIGN.md Sec.6; review change 12a).

Usage
  apo_describe.py --out results/cofold/controls/apo_burrow.tsv
  apo_describe.py --all-samples           # all 15 samples per system, not just top
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import subprocess
import sys
from pathlib import Path

import zstandard

P = Path(__file__).resolve().parents[2]  # repo root
PC = P / ".venv/bin/pocketcompare"
OUT_ROOT = Path("/tmp/pe_big/cofold/controls")
ABL = Path("/tmp/pe_big/ablation")
WORK = Path("/tmp/pe_big/cofold/apo_describe")


def decompress(src: str, dst: Path) -> Path:
    raw = open(src, "rb").read()
    if src.endswith(".zst"):
        raw = zstandard.ZstdDecompressor().decompress(raw, max_output_size=2 ** 31)
    dst.write_bytes(raw)
    return dst


def models_for(system: str, all_samples: bool) -> list[tuple[str, str]]:
    jd = OUT_ROOT / f"controls_{system}" / f"{system}_apo"
    if all_samples:
        out = []
        for mp in sorted(glob.glob(f"{jd}/seed-*_sample-*/*_model.cif*")):
            m = re.search(r"seed-(\d+)_sample-(\d+)", mp)
            out.append((f"{system}_s{m.group(1)}x{m.group(2)}", mp))
        return out
    top = [f for f in glob.glob(f"{jd}/*_model.cif*") if "seed-" not in Path(f).parent.name]
    return [(f"{system}_top", top[0])] if top else []


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(P / "results/cofold/controls/apo_burrow.tsv"))
    ap.add_argument("--all-samples", action="store_true")
    ap.add_argument("--systems", default="ca2,cdk2,pygl,yadh1")
    a = ap.parse_args()
    import pandas as pd

    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in a.systems.split(","):
        sj = json.loads((ABL / f"{s}_site.json").read_text())
        site = {f"A:{p}" for p in sj["site_uniprot_1based"]}   # chain A, model numbering = UniProt
        for label, mp in models_for(s, a.all_samples):
            cif = decompress(mp, WORK / f"{label}.cif")
            pq = WORK / f"{label}.parquet"
            r = subprocess.run([str(PC), "describe", str(cif), "-o", str(pq), "--af-mode"],
                               capture_output=True, text=True)
            if r.returncode != 0:
                print(f"[apo_describe] {label}: describe failed\n{r.stderr[-800:]}",
                      file=sys.stderr)
                continue
            df = pd.read_parquet(pq)
            best = None
            for _, p in df.iterrows():
                res = set(str(p["residues"]).split())
                if not res:
                    continue
                inter = len(res & site)
                rec = inter / len(site)
                prec = inter / len(res)
                jac = inter / len(res | site)
                cand = {"recall": rec, "precision": prec, "jaccard": jac,
                        "pocket_id": p["pocket_id"], "burrow_rank": int(p["burrow_rank"]),
                        "burrow_score": round(float(p["burrow_score"]), 4),
                        "volume_A3": round(float(p["volume_A3"]), 1),
                        "plddt_mean": round(float(p.get("plddt_mean", float("nan"))), 1),
                        "n_res": len(res)}
                if best is None or cand["jaccard"] > best["jaccard"]:
                    best = cand
            rows.append({"system": s, "model": label, "n_pockets": len(df),
                         "n_site_res": len(site),
                         **({k: (round(v, 3) if isinstance(v, float) else v)
                             for k, v in best.items()} if best else {}),
                         "detected": bool(best and best["recall"] >= 0.5)})
            print(f"[apo_describe] {label}: {len(df)} pockets; best vs known site "
                  f"recall={best['recall']:.2f} jaccard={best['jaccard']:.2f} "
                  f"rank={best['burrow_rank']}" if best else f"[apo_describe] {label}: no pockets")
    if rows:
        keys = sorted({k for r in rows for k in r},
                      key=lambda k: ["system", "model", "n_pockets", "n_site_res", "recall",
                                     "precision", "jaccard", "detected"].index(k)
                      if k in ["system", "model", "n_pockets", "n_site_res", "recall",
                               "precision", "jaccard", "detected"] else 99)
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        with open(a.out, "w") as f:
            f.write("\t".join(keys) + "\n")
            for r in rows:
                f.write("\t".join("" if r.get(k) is None else str(r[k]) for k in keys) + "\n")
        print(f"[apo_describe] wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
