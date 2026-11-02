#!/usr/bin/env python3
"""Where does the cognate ligand bind, relative to the protein's own BURROW pockets?

For each gate-passing Tier 1 case, take every member's top-ranked COGNATE co-fold model:
  1. ligand site  = protein residues with a heavy atom within 4.5 A of the ligand;
  2. strip the ligand and run `pocketcompare describe --af-mode` on the protein alone;
  3. for every BURROW pocket compute recall/precision against (a) the ligand site and
     (b) the case footprint for that member;
  4. report the pocket that best covers the ligand (id, BURROW rank) and the pocket that
     best covers the footprint, and whether they are the same pocket.

This answers whether the gained/lost atlas pocket is a SECONDARY site relative to where
the cognate ligand actually binds. Reads blind_key.json only to know which jobs are
cognate and which member/footprint they belong to; the binding verdicts are already
fixed upstream, so this is a descriptive follow-up, not a re-scoring.

Usage: ligand_pocket_match.py --cases gain_HOG0035156_g10_n21_O1 ...
Writes results/cofold/cases/<case>/ligand_vs_burrow_pockets.tsv and a combined JSON.
"""
from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cofold_readouts import read_model, pick_ligand, heavy, is_aa  # noqa: E402

P = Path(__file__).resolve().parents[2]  # repo root
PC = P / ".venv/bin/pocketcompare"
BIG = Path("/tmp/pe_big/cofold/cases")
WORK = Path("/tmp/pe_big/cofold/ligand_pocket_match")
CUTOFF = 4.5
MIN_MATCH_RECALL = 0.3   # continuous recall, floor only to call a pocket "matched"


def top_model(job_dir: Path) -> str | None:
    f = [x for x in glob.glob(f"{job_dir}/*_model.cif*") if "seed-" not in Path(x).parent.name]
    return f[0] if f else None


def ligand_site_and_apo(model: str, out_cif: Path) -> tuple[list[int], int]:
    st = read_model(model)
    lc, lig = pick_ligand(st, None)
    if lig is None:
        return [], 0
    lpos = np.array([a.pos.tolist() for a in heavy(lig)])
    site = []
    keep = gemmi.Structure()
    keep.add_model(gemmi.Model("1"))
    for ch in st[0]:
        if not len(ch.get_polymer()):
            continue                                  # drop ligand / ion chains
        nc = gemmi.Chain(ch.name)
        for i, r in enumerate(ch):
            if not is_aa(r):
                continue
            nc.add_residue(r)
            xyz = np.array([a.pos.tolist() for a in heavy(r)])
            if len(xyz) and np.min(np.linalg.norm(xyz[:, None] - lpos[None], axis=-1)) < CUTOFF:
                site.append(r.seqid.num)
        keep[0].add_chain(nc)
    keep.setup_entities()
    doc = keep.make_mmcif_document()
    doc.write_file(str(out_cif))
    return sorted(site), len(lpos)


def score_pockets(pq: Path, ligand_site: list[int], footprint: list[int]) -> dict:
    import pandas as pd
    df = pd.read_parquet(pq)
    L = {f"A:{p}" for p in ligand_site}
    F = {f"A:{p}" for p in footprint}
    rows = []
    for _, p in df.iterrows():
        res = set(str(p["residues"]).split())
        if not res:
            continue
        rows.append({"pocket_id": p["pocket_id"], "burrow_rank": int(p["burrow_rank"]),
                     "burrow_score": round(float(p["burrow_score"]), 4), "n_res": len(res),
                     "lig_recall": round(len(res & L) / len(L), 3) if L else None,
                     "lig_precision": round(len(res & L) / len(res), 3) if L else None,
                     "fp_recall": round(len(res & F) / len(F), 3) if F else None})
    best_lig = max(rows, key=lambda r: (r["lig_recall"] or 0, r["lig_precision"] or 0),
                   default=None)
    best_fp = max(rows, key=lambda r: r["fp_recall"] or 0, default=None)
    return {"n_pockets": len(rows), "best_ligand_pocket": best_lig, "best_footprint_pocket": best_fp,
            # "Same pocket" needs real coverage of BOTH sites. A first version accepted any
            # overlap and reported same=True on a footprint recall of 0.083 (1 of 12
            # residues), which is not a match.
            "same_pocket": bool(best_lig and best_fp and best_lig["pocket_id"] == best_fp["pocket_id"]
                                and (best_lig["lig_recall"] or 0) >= MIN_MATCH_RECALL
                                and (best_fp["fp_recall"] or 0) >= MIN_MATCH_RECALL),
            "footprint_detected": bool(best_fp and (best_fp["fp_recall"] or 0) >= MIN_MATCH_RECALL),
            "ligand_site_fp_overlap": len(set(ligand_site) & set(footprint))}


def run_case(cid: str) -> list[dict]:
    wd = BIG / cid
    key = json.loads((wd / "blind_key.json").read_text())["map"]
    jobs = {j["blind_id"]: j for j in json.loads((wd / "jobs.json").read_text())}
    out = []
    (WORK / cid).mkdir(parents=True, exist_ok=True)
    for bid, k in sorted(key.items(), key=lambda x: (x[1]["state"], x[1]["acc"])):
        if k["ligand_kind"] != "cognate" or bid not in jobs:
            continue
        jd = next((p for p in (wd / "af3").glob(f"*/{jobs[bid]['name']}") if p.is_dir()), None)
        mp = top_model(jd) if jd else None
        if not mp:
            continue
        apo = WORK / cid / f"{k['acc']}_protein_only.cif"
        lsite, nlig = ligand_site_and_apo(mp, apo)
        pq = WORK / cid / f"{k['acc']}.parquet"
        r = subprocess.run([str(PC), "describe", str(apo), "-o", str(pq), "--af-mode"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"[match] {cid} {k['acc']}: describe failed {r.stderr[-200:]}", file=sys.stderr)
            continue
        s = score_pockets(pq, lsite, jobs[bid].get("site") or [])
        bl, bf = s["best_ligand_pocket"] or {}, s["best_footprint_pocket"] or {}
        row = {"case": cid, "acc": k["acc"], "state": k["state"], "n_lig_atoms": nlig,
               "ligand_site_n": len(lsite), "footprint_n": len(jobs[bid].get("site") or []),
               "ligand_site_fp_overlap": s["ligand_site_fp_overlap"], "n_pockets": s["n_pockets"],
               "ligand_pocket_id": bl.get("pocket_id"), "ligand_pocket_rank": bl.get("burrow_rank"),
               "ligand_pocket_recall": bl.get("lig_recall"),
               "ligand_pocket_precision": bl.get("lig_precision"),
               "footprint_pocket_id": bf.get("pocket_id"), "footprint_pocket_rank": bf.get("burrow_rank"),
               "footprint_pocket_recall": bf.get("fp_recall"),
               "footprint_detected": s["footprint_detected"], "same_pocket": s["same_pocket"]}
        out.append(row)
        print(f"[match] {cid[:24]:24s} {k['acc']:11s} {k['state']:9s} lig->{row['ligand_pocket_id']}"
              f"(rank {row['ligand_pocket_rank']}, rec {row['ligand_pocket_recall']}) "
              f"fp->{row['footprint_pocket_id']}(rank {row['footprint_pocket_rank']}, "
              f"rec {row['footprint_pocket_recall']}) same={row['same_pocket']}", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="+", required=True)
    a = ap.parse_args()
    import pandas as pd
    allrows = []
    for cid in a.cases:
        rows = run_case(cid)
        allrows += rows
        if rows:
            d = P / "results/cofold/cases" / cid
            d.mkdir(parents=True, exist_ok=True)
            tmp = d / "ligand_vs_burrow_pockets.tsv.tmp"
            pd.DataFrame(rows).to_csv(tmp, sep="\t", index=False)
            tmp.replace(d / "ligand_vs_burrow_pockets.tsv")
    tmp = P / "results/cofold/cases/ligand_vs_burrow_pockets.json.tmp"
    tmp.write_text(json.dumps(allrows, indent=1))
    tmp.replace(P / "results/cofold/cases/ligand_vs_burrow_pockets.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
