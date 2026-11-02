#!/usr/bin/env python3
"""Export representative co-fold models into cases/<id>/cofold/models/ for the web app.

Per case: for each TREE side (case.yaml proteins[].role: derived -> "derived",
outgroup -> "ancestral", which is what webapp/SCHEMA.md Sec.4 uses to fill the right and
left viewers), take the member whose COGNATE-ligand job has the best AF3 ranking score and
copy that job's top-ranked model (top seed) as gzipped mmCIF. The family positive control
is added with role "control" so the gate can be inspected. Metrics are read from the
already-computed blinded readouts, so nothing is re-scored here.

Using the tree role rather than the pocket state keeps a loss case the right way round:
there the derived lineage is the pocket-ABSENT side.

Appends one `cofold[]` entry per model (model_id, accession, ligand, seed, role, file,
metrics) after the existing summary entry; re-running replaces them (idempotent).
Refuses any file over --max-mb.

Usage: export_case_models.py --cases <id> ... [--workdir-root /tmp/pe_big/cofold/cases]
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import json
import sys
from pathlib import Path

import yaml
import zstandard

P = Path(__file__).resolve().parents[2]  # repo root


def read_tsv(path: Path) -> dict[str, dict]:
    with open(path) as f:
        return {r["blind_id"]: r for r in csv.DictReader(f, delimiter="\t")}


def top_model(job_dir: Path) -> Path | None:
    f = [Path(x) for x in glob.glob(f"{job_dir}/*_model.cif*") if "seed-" not in Path(x).parent.name]
    return f[0] if f else None


def ranking(job_dir: Path) -> tuple[float, int | None]:
    """(best ranking score, its seed) from the job's ranking CSV."""
    rc = glob.glob(f"{job_dir}/*_ranking_scores.csv")
    best, seed = -1.0, None
    if rc:
        for line in Path(rc[0]).read_text().splitlines()[1:]:
            if line.strip():
                s, _, r = line.split(",")
                if float(r) > best:
                    best, seed = float(r), int(s)
    return best, seed


def fnum(x):
    try:
        return round(float(x), 3)
    except (TypeError, ValueError):
        return None


def export(cid: str, root: Path, max_mb: float, case_dir: str | None = None,
           spec_path: str | None = None) -> list[dict]:
    """`cid` is the co-fold workdir name. `case_dir` is the cases/<dir> whose case.yaml
    gets the entries, when it differs (the deep dives are specified under cases/cofold/ but
    belong to the val_* case pages)."""
    wd = root / cid
    cdir = case_dir or cid
    key = json.loads((wd / "blind_key.json").read_text())["map"]
    jobs = {j["blind_id"]: j for j in json.loads((wd / "jobs.json").read_text())}
    readouts = read_tsv(wd / "blinded_readouts.tsv")
    doc = yaml.safe_load((P / "cases" / cdir / "case.yaml").read_text())
    tree_role = {p["accession"]: p.get("role") for p in doc.get("proteins") or []}
    if spec_path:                       # deep dives: roles come from the co-fold spec
        sp = yaml.safe_load((P / spec_path).read_text())
        tree_role = {m["accession"]: m["state"] for m in sp.get("members") or []}
    outd = P / "cases" / cdir / "cofold" / "models"
    outd.mkdir(parents=True, exist_ok=True)

    def job_dir(bid):
        return next((p for p in (wd / "af3").glob(f"*/{jobs[bid]['name']}") if p.is_dir()), None)

    # best cognate job per tree side
    best: dict[str, tuple] = {}
    for bid, k in key.items():
        if k["ligand_kind"] != "cognate" or bid not in jobs:
            continue
        jd = job_dir(bid)
        if jd is None:
            continue
        side = "derived" if tree_role.get(k["acc"]) == "derived" else "ancestral"
        score, seed = ranking(jd)
        if side not in best or score > best[side][0]:
            best[side] = (score, seed, bid, k, jd)
    picks = [(side, *v) for side, v in sorted(best.items())]
    # positive control
    pcj = next((j for j in jobs.values() if j["blind_id"] == "poscontrol"), None)
    if pcj:
        jd = next((p for p in (wd / "af3").glob(f"*/{pcj['name']}") if p.is_dir()), None)
        if jd:
            sc, sd = ranking(jd)
            pc = (doc.get("cofold") or [{}])[0].get("positive_control_gate") or {}
            picks.append(("control", sc, sd, "poscontrol",
                          {"acc": (pc.get("ref") or "poscontrol").split(":")[0], "state": "control",
                           "ligands": "+".join(l.split(":")[-1] for l in pcj["ligands"])}, jd))

    entries = []
    for side, score, seed, bid, k, jd in picks:
        mp = top_model(jd)
        if mp is None:
            continue
        raw = mp.read_bytes()
        if mp.suffix == ".zst":
            raw = zstandard.ZstdDecompressor().decompress(raw, max_output_size=2 ** 31)
        acc = k["acc"] if bid != "poscontrol" else "poscontrol"
        lig = (k.get("ligands") or "").replace("+", "-")
        dst = outd / f"{side}_{acc}_{lig}_top.cif.gz"
        tmp = dst.with_suffix(".gz.tmp")
        with gzip.open(tmp, "wb", compresslevel=9) as f:
            f.write(raw)
        mb = tmp.stat().st_size / 1e6
        if mb > max_mb:
            tmp.unlink()
            print(f"[export] {cid} {side}: {mb:.2f} MB > {max_mb} MB, skipped", file=sys.stderr)
            continue
        tmp.replace(dst)
        r = readouts.get(bid, {})
        entries.append({
            "model_id": f"{side}_{acc}",
            "accession": acc if bid != "poscontrol" else k["acc"],
            "ligand": k.get("ligands"),
            "seed": seed,
            "role": side,
            "file": str(dst.relative_to(P)),
            "ranking_score": fnum(score),
            "metrics": {"ligand_plddt": fnum(r.get("top_lig_plddt")),
                        "min_pae": fnum(r.get("top_prot_lig_pae_min")),
                        "centroid_dist_A": fnum(r.get("top_site_centroid_dist")),
                        "site_recall": fnum(r.get("top_site_recall")),
                        "posebusters_pass": {"True": True, "False": False}.get(r.get("pb_valid"))},
            "note": ("top-ranked sample of the best-ranked cognate job on this side; "
                     "centroid/recall are measured against the case footprint")
            if bid != "poscontrol" else
            "family holo positive control (gate); centroid/recall fields refer to the case "
            "footprint and are not meaningful for this protein",
        })
        print(f"[export] {cid} {side:9s} {acc:12s} {lig:10s} seed {seed} {mb:.2f} MB -> {dst.name}")

    # write back: keep the summary entry (no `file`), replace model entries
    cof = [c for c in (doc.get("cofold") or []) if not (isinstance(c, dict) and c.get("file"))]
    doc["cofold"] = cof + entries
    f = P / "cases" / cdir / "case.yaml"
    tmp = f.with_suffix(".yaml.tmp")
    tmp.write_text(yaml.safe_dump(doc, sort_keys=False, width=100))
    tmp.replace(f)
    return entries


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="+", required=True)
    ap.add_argument("--workdir-root", default="/tmp/pe_big/cofold/cases")
    ap.add_argument("--max-mb", type=float, default=2.0)
    ap.add_argument("--case-dir", help="cases/<dir> receiving the entries (single case)")
    ap.add_argument("--spec", help="co-fold spec whose members[].state gives the roles")
    a = ap.parse_args()
    n = 0
    for cid in a.cases:
        n += len(export(cid, Path(a.workdir_root), a.max_mb, a.case_dir, a.spec))
    print(f"[export] {n} models written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
