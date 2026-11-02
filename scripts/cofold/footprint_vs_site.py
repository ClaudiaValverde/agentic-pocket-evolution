#!/usr/bin/env python3
"""How far is a case's BURROW footprint from the cognate ligand's crystal site?

Tier 1 co-folding returned "uninformative" for every case: the cognate ligand docked
essentially onto its crystal pose (positive-control RMSD 0.37-1.34 A) yet sat 12-34 A from
the footprint the case is about, with site recall 0.00. This script quantifies that
directly, so the mismatch is a measured number rather than an inference.

Method. Both residue sets are expressed in UniProt numbering on the positive control's
protein, which is itself a case member, so no cross-protein numbering is involved:
  * cognate site = BioLiP2 residues_uniprot for that accession and ligand;
  * case footprint = the case's present-side footprint, mapped onto that accession by
    pairwise alignment (identical to the member-site mapping in case_cofold);
  * distance = minimum C-alpha to C-alpha distance between the two sets, measured in the
    AF3 model of that protein that was folded for the gate (chain A, UniProt numbering).
Also reports the overlap, the centroid separation and each set's size.

Usage: footprint_vs_site.py --all
"""
from __future__ import annotations
import argparse, glob, json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_cofold as CC
import colabfold_msa as CF
from cofold_readouts import read_model, is_aa

P = Path(__file__).resolve().parents[2]  # repo root
import run_pool as RP


def ca_coords(model_path: str, positions: list[int]):
    st = read_model(model_path)
    prot = [r for ch in st[0] for r in ch if is_aa(r)]
    out = {}
    for p in positions:
        if 0 < p <= len(prot):
            a = prot[p - 1].find_atom("CA", "*")
            if a:
                out[p] = np.array(a.pos.tolist())
    return out


def analyse(cid: str, cfg: dict) -> dict | None:
    import pandas as pd
    pc = cfg["positive_control"]
    acc, lig = pc["acc"], pc["ligand"]
    bl = pd.read_parquet(P / "data/annotations/biolip_sites.parquet",
                         columns=["uniprot_acc", "pdb_id", "chain", "ligand",
                                  "residues_uniprot", "n_site_residues", "resolution"])
    rows = bl[(bl.uniprot_acc == acc) & (bl.ligand == lig) & (bl.pdb_id == pc["pdb"])]
    if rows.empty:
        rows = bl[(bl.uniprot_acc == acc) & (bl.ligand == lig)]
    if rows.empty:
        return None
    site = sorted({int(x) for x in rows.iloc[0].residues_uniprot})

    spec = CC.load_spec(P / "cases" / cid / "cofold/cofold_spec.yaml")
    ref = next((m for m in spec["members"] if m["state"] == "derived" and m.get("site")), None)
    if ref is None:
        return None
    pc_seq = CF.fetch_uniprot(acc)
    fp = (ref["site"] if ref["accession"] == acc
          else CC.map_positions(CF.fetch_uniprot(ref["accession"]), pc_seq, ref["site"]))

    wd = Path("/tmp/pe_big/cofold/cases") / cid
    mp = [f for f in glob.glob(f"{wd}/af3/*/{cid}_poscontrol/*_model.cif*")
          if "seed-" not in Path(f).parent.name]
    if not mp:
        return None
    ca = ca_coords(mp[0], sorted(set(site) | set(fp)))
    S = np.array([ca[p] for p in site if p in ca])
    F = np.array([ca[p] for p in fp if p in ca])
    if not len(S) or not len(F):
        return None
    D = np.linalg.norm(F[:, None] - S[None], axis=-1)
    return {
        "positive_control_acc": acc, "pdb": str(rows.iloc[0].pdb_id), "ligand": lig,
        "cognate_site_uniprot": site, "n_cognate_site": len(site),
        "case_footprint_on_pc": fp, "n_footprint": len(fp),
        "footprint_source_member": ref["accession"],
        "n_shared_residues": len(set(site) & set(fp)),
        "min_ca_distance_A": round(float(D.min()), 2),
        "median_min_ca_distance_A": round(float(np.median(D.min(axis=1))), 2),
        "centroid_separation_A": round(float(np.linalg.norm(F.mean(0) - S.mean(0))), 2),
        "measured_in": str(Path(mp[0]).relative_to(Path("/tmp/pe_big"))),
        "conclusion": None,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--case", nargs="*")
    a = ap.parse_args()
    cases = a.case or (sorted(RP.TIER1) if a.all else [])
    out = {}
    for cid in cases:
        r = analyse(cid, RP.TIER1[cid])
        if r is None:
            print(f"{cid}: could not measure", file=sys.stderr); continue
        far = r["min_ca_distance_A"] > 8.0 and r["n_shared_residues"] == 0
        r["conclusion"] = (
            "The case footprint and the cognate ligand site are DISJOINT: "
            f"{r['n_shared_residues']} shared residues, closest C-alpha pair "
            f"{r['min_ca_distance_A']} A apart, centroids {r['centroid_separation_A']} A apart. "
            "Co-folding with this ligand cannot test the case pocket, so the D23 result is "
            "uninformative about the gain/loss rather than evidence against it."
            if far else
            f"Footprint and cognate site are adjacent or overlapping "
            f"({r['n_shared_residues']} shared, min C-alpha {r['min_ca_distance_A']} A); "
            "co-folding with this ligand does probe the case pocket.")
        out[cid] = r
        print(f"{cid}: shared={r['n_shared_residues']} min_CA={r['min_ca_distance_A']} A "
              f"centroid={r['centroid_separation_A']} A "
              f"(footprint {r['n_footprint']} res vs site {r['n_cognate_site']} res, "
              f"{r['positive_control_acc']}/{r['pdb']}:{r['ligand']})")
    d = P / "results/cofold/cases"
    tmp = d / "footprint_vs_cognate_site.json.tmp"
    tmp.write_text(json.dumps(out, indent=1))
    tmp.replace(d / "footprint_vs_cognate_site.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
