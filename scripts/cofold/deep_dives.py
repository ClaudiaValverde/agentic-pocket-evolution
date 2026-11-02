#!/usr/bin/env python3
"""Known-answer deep dives: can co-folding see a documented specificity change at all?

Tier 1 could not test the pipeline's pockets, because those footprints turn out not to be
the families' ligand sites. These two systems remove that confound entirely: the pocket IS
the annotated ligand site, and the right answer is known from biochemistry. If co-folding
cannot separate these, it cannot be used as specificity evidence anywhere.

A. Steroid receptor LBPs (cases/val_steroid_receptor_LBP).
   MR (P08235) and GR (P04150) x aldosterone (AS4) and cortisol (HCY).
   Known answer: MR binds BOTH with high affinity; GR is selective for cortisol and binds
   aldosterone far more weakly. Expected ordering MR/AS4 ~ MR/HCY ~ GR/HCY > GR/AS4.
   Gate: MR + aldosterone against 2aa2 (1.95 A).

B. LDH vs MDH (cases/val_PfLDH_vs_MDH).
   P. falciparum LDH (Q76NM3), human MDH1 (P40925) and human LDHA (P00338, included so
   the panel has a protein with abundant holo structures to gate on) x pyruvate (PYR) and
   oxaloacetate (OAA), each with NADH (NAI).
   Known answer: LDH turns over pyruvate, MDH oxaloacetate. Expected: the diagonal
   (LDH/PYR, MDH/OAA) beats the off-diagonal.

Both are scored with the FROZEN D23 thresholds; the cognate/non-cognate contrast is the
readout, not a new threshold. Receptors over 500 aa are folded as the TED domain holding
the ligand site, so the LBD is folded rather than the full receptor.

Usage: deep_dives.py --build --run --score --gpus 0,1,2
"""
from __future__ import annotations
import argparse, glob, json, os, subprocess, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import colabfold_msa as CF
import run_apo_pool as AP
from cofold_readouts import summarise_job

P = Path(__file__).resolve().parents[2]  # repo root
RUN = P / "scripts/cofold/run_af3.sh"
BIG = Path("/tmp/pe_big/cofold/deep_dives")
RES = P / "results/cofold/deep_dives"

PANELS = {
    "steroid_receptor": {
        "proteins": {
            "MR_P08235": {"acc": "P08235", "site_ligand": "AS4", "role": "promiscuous"},
            "GR_P04150": {"acc": "P04150", "site_ligand": "HCY", "role": "cortisol-selective"},
        },
        "ligands": {"aldosterone": ["CCD:AS4"], "cortisol": ["CCD:HCY"]},
        "cognate": {"MR_P08235": ["aldosterone", "cortisol"], "GR_P04150": ["cortisol"]},
        "gate": {"protein": "MR_P08235", "ligand": "aldosterone",
                 "pdb": "2aa2", "chain": "A", "ccd": "AS4"},
        "max_aa": 500,
    },
    "ldh_mdh": {
        "proteins": {
            "PfLDH_Q76NM3": {"acc": "Q76NM3", "site_ligand": None, "role": "LDH (derived)",
                             "site_from": "P00338"},
            "MDH1_P40925": {"acc": "P40925", "site_ligand": "MLI", "role": "MDH (ancestral)"},
            "LDHA_P00338": {"acc": "P00338", "site_ligand": "NAI", "role": "LDH (gate)"},
        },
        "ligands": {"pyruvate": ["CCD:PYR", "CCD:NAI"],
                    "oxaloacetate": ["CCD:OAA", "CCD:NAI"]},
        "cognate": {"PfLDH_Q76NM3": ["pyruvate"], "LDHA_P00338": ["pyruvate"],
                    "MDH1_P40925": ["oxaloacetate"]},
        "gate": {"protein": "LDHA_P00338", "ligand": "pyruvate",
                 "pdb": "5w8j", "chain": "A", "ccd": "NAI"},
        "max_aa": 600,
    },
}


def log(m): print(f"[deep_dives] {m}", file=sys.stderr, flush=True)


def biolip_site(acc: str, ligand: str | None):
    import pandas as pd
    bl = pd.read_parquet(P / "data/annotations/biolip_sites.parquet",
                         columns=["uniprot_acc", "ligand", "residues_uniprot",
                                  "n_site_residues", "resolution", "pdb_id"])
    s = bl[bl.uniprot_acc == acc]
    if ligand:
        s = s[s.ligand == ligand]
    if s.empty:
        return [], None
    s = s.sort_values("n_site_residues", ascending=False)
    return sorted({int(x) for x in s.iloc[0].residues_uniprot}), str(s.iloc[0].pdb_id)


def build(panel: str) -> list[dict]:
    cfg = PANELS[panel]
    jobs = []
    for pname, pc in cfg["proteins"].items():
        acc = pc["acc"]
        seq = CF.fetch_uniprot(acc)
        site, pdb = biolip_site(acc, pc.get("site_ligand"))
        if not site and pc.get("site_from"):
            src = pc["site_from"]
            ssite, _ = biolip_site(src, PANELS[panel]["proteins"][f"LDHA_{src}"]["site_ligand"]
                                   if f"LDHA_{src}" in cfg["proteins"] else None)
            import case_cofold as CC
            site = CC.map_positions(CF.fetch_uniprot(src), seq, ssite)
            log(f"{pname}: site mapped from {src} ({len(site)} positions)")
        construct = "full_length"
        if len(seq) > cfg["max_aa"] and site:
            tc = AP.ted_construct(acc, site, seq)
            if tc:
                seq, site, construct = tc
                log(f"{pname}: {construct}")
        prov = CF.fetch_msa(f"{pname}_{construct}", seq, host=CF.HOST, mode="env",
                            cache=CF.CACHE, min_interval=20.0)
        for lname, ligs in cfg["ligands"].items():
            name = f"{panel}_{pname}_{lname}".lower()
            d = BIG / "inputs" / panel
            d.mkdir(parents=True, exist_ok=True)
            out = d / f"{name}.json"
            CF.write_af3_json(str(out), name, seq, CF.load_a3m(seq, CF.CACHE), ligs,
                              [1, 2, 3], prov)
            m = out.with_suffix(".meta.json")
            if m.exists():
                (d / "meta").mkdir(exist_ok=True); m.rename(d / "meta" / m.name)
            jobs.append({"panel": panel, "protein": pname, "acc": acc, "ligand": lname,
                         "ligand_ccds": ligs, "name": name, "length": len(seq),
                         "site": site, "construct": construct,
                         "cognate": lname in cfg["cognate"].get(pname, []),
                         "role": pc["role"], "site_pdb": pdb})
    return jobs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--gpus", default="0,1,2")
    a = ap.parse_args()
    BIG.mkdir(parents=True, exist_ok=True); RES.mkdir(parents=True, exist_ok=True)
    if a.build:
        jobs = []
        for panel in PANELS:
            jobs += build(panel)
        tmp = BIG / "jobs.json.tmp"; tmp.write_text(json.dumps(jobs, indent=1))
        tmp.replace(BIG / "jobs.json")
        log(f"built {len(jobs)} jobs")
    if a.run:
        gpus = [int(g) for g in a.gpus.split(",")]
        panels = list(PANELS)
        procs = []
        for i, panel in enumerate(panels):
            g = gpus[i % len(gpus)]
            lf = open(BIG / f"gpu{g}_{panel}.out", "w")
            procs.append((panel, subprocess.Popen(
                ["bash", str(RUN), "-i", str(BIG / "inputs" / panel), "-j", f"dd_{panel}",
                 "-g", str(g)], stdout=lf, stderr=subprocess.STDOUT,
                env={**os.environ, "AF3_OUT_ROOT": str(BIG / "af3")}), lf, time.time()))
            log(f"{panel} -> GPU {g}")
        for panel, p, f, t0 in procs:
            rc = p.wait(); f.close()
            log(f"{panel} done rc={rc} in {time.time()-t0:.0f}s")
    if a.score:
        import pandas as pd
        jobs = json.loads((BIG / "jobs.json").read_text())
        th = json.loads((P / "results/cofold/calibration/threshold.json").read_text())
        t = th["d23_composite"]["thresholds"]
        rows = []
        for j in jobs:
            jd = BIG / "af3" / f"dd_{j['panel']}" / j["name"]
            if not jd.exists():
                log(f"missing {jd}"); continue
            samp, summ = summarise_job(jd, j["site"], pb_top=True)
            if summ is None:
                continue
            per = {}
            for sd in sorted({s["seed"] for s in samp}):
                ss = [s for s in samp if s["seed"] == sd]
                per[sd] = {k: float(np.mean([s[k] for s in ss if s.get(k) is not None]))
                           for k in ("lig_plddt", "site_centroid_dist", "site_recall")}
            nb = sum(1 for v in per.values()
                     if v["lig_plddt"] >= t["t_plddt"]
                     and (v["site_centroid_dist"] <= t["t_centroid"]
                          or v["site_recall"] >= t["t_recall"]))
            rows.append({"panel": j["panel"], "protein": j["protein"], "role": j["role"],
                         "ligand": j["ligand"], "cognate": j["cognate"],
                         "construct": j["construct"], "n_site": len(j["site"]),
                         "lig_plddt": summ["med_lig_plddt"],
                         "pae_min": summ["med_prot_lig_pae_min"],
                         "centroid": summ["med_site_centroid_dist"],
                         "site_recall": summ["med_site_recall"],
                         "d23_binds": nb >= 2, "n_seeds_binding": nb,
                         "pb_valid": summ.get("pb_valid")})
            log(f"{j['panel']} {j['protein']:14s} {j['ligand']:13s} "
                f"cognate={j['cognate']} pLDDT={summ['med_lig_plddt']} "
                f"cen={summ['med_site_centroid_dist']} rec={summ['med_site_recall']} "
                f"binds={nb>=2}")
        df = pd.DataFrame(rows)
        tmp = RES / "deep_dives.tsv.tmp"; df.to_csv(tmp, sep="\t", index=False)
        tmp.replace(RES / "deep_dives.tsv")
        if len(df):
            print(df.to_string(index=False))
            summ = {"frozen_thresholds": t, "n_jobs": len(df),
                    "cognate_binds": int(df[df.cognate].d23_binds.sum()),
                    "n_cognate": int(df.cognate.sum()),
                    "noncognate_binds": int(df[~df.cognate].d23_binds.sum()),
                    "n_noncognate": int((~df.cognate).sum()),
                    "cognate_median_plddt": float(df[df.cognate].lig_plddt.median()),
                    "noncognate_median_plddt": float(df[~df.cognate].lig_plddt.median())}
            tmp = RES / "deep_dives_summary.json.tmp"
            tmp.write_text(json.dumps(summ, indent=1))
            tmp.replace(RES / "deep_dives_summary.json")
            print(json.dumps(summ, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
