#!/usr/bin/env python3
"""Preregistered co-folding readouts (DESIGN.md Sec.6 / review change 10).

Per diffusion sample of an AF3 job:
  lig_plddt            mean pLDDT over the query ligand's heavy atoms
  prot_lig_pae_min     min protein-ligand PAE (summary_confidences chain_pair_pae_min)
  site_centroid_dist   |ligand centroid - centroid of the known site CA atoms|, computed
                       inside the predicted model, so it works for variants and needs no
                       holo superposition
  site_recall          fraction of known-site residues with a heavy atom within 5 A of
                       the predicted ligand
  prot_lig_iptm        chain_pair_iptm (recorded, never used alone: decoys reach 0.72-0.86)
  pb_valid             PoseBusters "dock" validity, top-ranked sample only (optional)

The site is given as 1-based positions in the folded sequence; because we always fold the
full UniProt (or variant) sequence, those indices map 1:1 onto model residues.

Importable (`summarise_job`) and runnable:
  cofold_readouts.py --job-dir /tmp/pe_big/cofold/controls/cdk2_wt --site-json /tmp/.../cdk2_site.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import statistics as stats
import sys
import subprocess
import tempfile
from pathlib import Path

import gemmi
import numpy as np
import zstandard

PB_PY = os.environ.get("PB_PY", "/tmp/pb_env/bin/python")

# SMILES for bond-order assignment when handing ligands to PoseBusters.
LIG_SMILES = {
    "ATP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
    "AZM": "CC(=O)Nc1nnc(s1)S(N)(=O)=O",
    "CFF": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "PLP": "Cc1ncc(COP(=O)(O)O)c(C=O)c1O",
    "NAD": "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(=O)([O-])OP(=O)(O)OC[C@H]2O[C@@H]"
           "(n3cnc4c(N)ncnc34)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O",
}
IONS = {"ZN", "MG", "CA", "NA", "K", "MN", "FE", "CU", "CO", "NI", "CL"}


def read_model(path: str) -> gemmi.Structure:
    raw = open(path, "rb").read()
    if path.endswith(".zst"):
        raw = zstandard.ZstdDecompressor().decompress(raw, max_output_size=2 ** 31)
    doc = gemmi.cif.read_string(raw.decode())
    st = gemmi.make_structure_from_block(doc.sole_block())
    st.setup_entities()
    return st


def heavy(res):
    return [a for a in res if a.element.name != "H"]


def is_aa(r) -> bool:
    t = gemmi.find_tabulated_residue(r.name)
    return bool(t and t.is_amino_acid())


def pick_ligand(st, want: str | None):
    """Return (chain_id, residue) for the query ligand: the named one, else the biggest
    non-ion heteroresidue.

    Ligand chains are identified by TOPOLOGY, not by residue name: AF3 emits each ligand as
    its own chain containing a single residue with no polymer. Name-based detection breaks
    whenever the ligand is itself a standard residue -- the MetE case folds free methionine
    (CCD MET), which is_aa() happily calls a protein residue, so every sample came back
    with no ligand at all.
    """
    cands = []
    for ch in st[0]:
        is_polymer = len(ch.get_polymer()) > 0
        for r in ch:
            if r.name == "HOH":
                continue
            if is_polymer and is_aa(r):
                continue
            cands.append((ch.name, r))
    if not cands:
        return None, None
    if want:
        for cid, r in cands:
            if r.name == want.upper():
                return cid, r
    non_ion = [(c, r) for c, r in cands if r.name not in IONS]
    pool = non_ion or cands
    return max(pool, key=lambda cr: len(heavy(cr[1])))


def posebusters_check(st, lig_chain: str, lig_res) -> dict:
    """PoseBusters 'dock' config on the predicted complex. Returns {} if unavailable."""
    if not Path(PB_PY).exists():
        return {}
    with tempfile.TemporaryDirectory() as td:
        prot = gemmi.Structure()
        prot.add_model(gemmi.Model("1"))
        for ch in st[0]:
            keep = [r for r in ch if is_aa(r)]
            if not keep:
                continue
            nc = gemmi.Chain(ch.name)
            for r in keep:
                nc.add_residue(r)
            prot[0].add_chain(nc)
        prot.setup_entities()
        pp = Path(td) / "prot.pdb"
        prot.write_pdb(str(pp))
        lp = Path(td) / "lig.pdb"
        lst = gemmi.Structure()
        lst.add_model(gemmi.Model("1"))
        nc = gemmi.Chain(lig_chain)
        nc.add_residue(lig_res)
        lst[0].add_chain(nc)
        lst.setup_entities()
        lst.write_pdb(str(lp))
        smi = LIG_SMILES.get(lig_res.name, "")
        code = r'''
import sys, json, warnings
warnings.filterwarnings("ignore")
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
prot, lig, smi, out = sys.argv[1:5]
m = Chem.MolFromPDBBlock(open(lig).read(), removeHs=False, sanitize=False)
if m is None:
    json.dump({"error": "ligand unreadable"}, open(out, "w")); sys.exit()
if smi:
    t = Chem.MolFromSmiles(smi)
    try:
        m = AllChem.AssignBondOrdersFromTemplate(t, m)
    except Exception as e:
        pass
try:
    Chem.SanitizeMol(m)
except Exception:
    pass
sdf = lig.replace(".pdb", ".sdf")
w = Chem.SDWriter(sdf); w.write(m); w.close()
from posebusters import PoseBusters
df = PoseBusters(config="dock").bust([sdf], None, prot, full_report=False)
rec = df.to_dict("records")[0]
rec = {k: (None if v is None else bool(v)) for k, v in rec.items() if isinstance(v, (bool, int))}
json.dump(rec, open(out, "w"))
'''
        sp = Path(td) / "pb.py"
        sp.write_text(code)
        op = Path(td) / "pb.json"
        try:
            subprocess.run([PB_PY, "-I", str(sp), str(pp), str(lp), smi, str(op)],
                           check=True, capture_output=True, timeout=300)
            return json.loads(op.read_text())
        except Exception as e:  # noqa: BLE001
            return {"error": str(e)[:200]}


def eval_sample(model_path: str, site_1based: list[int], want_lig: str | None,
                run_pb: bool = False) -> dict:
    st = read_model(model_path)
    lig_chain, lig = pick_ligand(st, want_lig)
    if lig is None:
        return {"lig_resname": None}
    latoms = heavy(lig)
    lpos = np.array([a.pos.tolist() for a in latoms])
    prot = [r for ch in st[0] for r in ch if is_aa(r)]
    site = [prot[i - 1] for i in site_1based if 0 < i <= len(prot)]
    ca = np.array([r["CA"][0].pos.tolist() for r in site if r.find_atom("CA", "*")])
    out = {
        "lig_resname": lig.name, "lig_chain": lig_chain, "n_lig_atoms": len(latoms),
        "lig_plddt": round(float(np.mean([a.b_iso for a in latoms])), 1),
        "site_centroid_dist": round(float(np.linalg.norm(lpos.mean(0) - ca.mean(0))), 2)
        if len(ca) else None,
    }
    rec = 0
    for r in site:
        xyz = np.array([a.pos.tolist() for a in heavy(r)])
        if len(xyz) and np.min(np.linalg.norm(xyz[:, None] - lpos[None], axis=-1)) < 5.0:
            rec += 1
    out["site_recall"] = round(rec / max(1, len(site)), 2)
    out["n_site"] = len(site)
    if run_pb:
        pb = posebusters_check(st, lig_chain, lig)
        checks = {k: v for k, v in pb.items() if isinstance(v, bool)}
        out["pb_valid"] = (all(checks.values()) if checks else None)
        out["pb_failed"] = ",".join(sorted(k for k, v in checks.items() if v is False))
        out["pb_n_checks"] = len(checks)
    return out


def summarise_job(job_dir: str | Path, site_1based: list[int], want_lig: str | None = None,
                  pb_top: bool = True) -> tuple[list[dict], dict | None]:
    """Return (per-sample rows, top-ranked summary row) for one AF3 job directory."""
    jd = Path(job_dir)
    rank = {}
    rc = glob.glob(f"{jd}/**/*_ranking_scores.csv", recursive=True)
    if rc:
        for line in Path(rc[0]).read_text().split("\n")[1:]:
            if line.strip():
                s, i, r = line.split(",")
                rank[(int(s), int(i))] = float(r)
    rows = []
    for mp in sorted(glob.glob(f"{jd}/**/seed-*_sample-*/*_model.cif*", recursive=True)):
        m = re.search(r"seed-(\d+)_sample-(\d+)", mp)
        sd, si = int(m.group(1)), int(m.group(2))
        scf = glob.glob(os.path.join(os.path.dirname(mp), "*summary_confidences.json"))
        sc = json.load(open(scf[0])) if scf else {}
        r = {"job": jd.name, "seed": sd, "sample": si, "ranking": rank.get((sd, si)),
             "path": mp}
        r.update(eval_sample(mp, site_1based, want_lig, run_pb=False))
        if sc and r.get("lig_chain"):
            ci = sc["chain_ids"].index(r["lig_chain"])
            r["iptm"] = sc.get("iptm")
            r["prot_lig_iptm"] = sc["chain_pair_iptm"][0][ci]
            r["prot_lig_pae_min"] = sc["chain_pair_pae_min"][0][ci]
        rows.append(r)
    n_all = len(rows)
    rows = [r for r in rows if r.get("site_recall") is not None]
    if n_all and not rows:
        print(f"[cofold_readouts] {jd.name}: no ligand chain in any of {n_all} samples",
              file=sys.stderr)
    if not rows:
        return [], None
    top = max(rows, key=lambda r: (r["ranking"] if r["ranking"] is not None else -1))
    summ = {"job": jd.name, "n_samples": len(rows), "n_seeds": len({r["seed"] for r in rows}),
            "top_ranking": round(top["ranking"], 4) if top["ranking"] is not None else None,
            "lig_resname": top["lig_resname"], "n_site": top.get("n_site")}
    for k in ("lig_plddt", "prot_lig_pae_min", "site_centroid_dist", "site_recall",
              "prot_lig_iptm", "iptm"):
        vals = [r[k] for r in rows if r.get(k) is not None]
        summ[f"top_{k}"] = top.get(k)
        summ[f"med_{k}"] = round(stats.median(vals), 3) if vals else None
        summ[f"sd_{k}"] = round(stats.pstdev(vals), 3) if len(vals) > 1 else None
    summ["frac_in_site"] = round(float(np.mean([r["site_recall"] >= 0.5 for r in rows])), 2)
    if pb_top:
        pbr = eval_sample(top["path"], site_1based, want_lig, run_pb=True)
        summ["pb_valid"] = pbr.get("pb_valid")
        summ["pb_failed"] = pbr.get("pb_failed")
    return rows, summ


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--job-dir", required=True)
    ap.add_argument("--site-json", help="site_ablation.py JSON (uses site_uniprot_1based)")
    ap.add_argument("--site", help="comma-separated 1-based positions")
    ap.add_argument("--ligand-resname")
    ap.add_argument("--no-posebusters", action="store_true")
    a = ap.parse_args()
    site = ([int(x) for x in a.site.split(",")] if a.site
            else json.loads(Path(a.site_json).read_text())["site_uniprot_1based"])
    rows, summ = summarise_job(a.job_dir, site, a.ligand_resname, not a.no_posebusters)
    print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
