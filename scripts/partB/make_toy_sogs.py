#!/usr/bin/env python
"""Toy SOG state table with known truth, in the real input layout (webapp/SCHEMA.md sogs-web-v1), for end-to-end
validation of scripts/partB/run_events.py and run_nulls.py before real SOG states exist.

Real O1 root-HOG trees (>= 6 species, <= 300 members, seeded sample); per OG 1-3 SOGs, each
vertical (70 %), planted gain (15 %) or planted loss (15 %) on an eligible branch; observation noise = the default
NoiseModel (simulate.py). In 15 % of the vertical SOGs a specificity change is planted at a sister pair:
3 clade-consistent footprint substitutions with BLOSUM62 < 0 + a fingerprint mean shift on the derived clade.
Also written: replicate observations of the same proteins (replicate_states/, for the same-protein null),
footprint residues, fingerprints with covariates, and replicate cosines (technical noise floor).
All ids start with TOY; manifest has toy: true.

    python scripts/partB/make_toy_sogs.py --out results/partB/sim/toy_sogs/TOY/O1
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import events, recon, simulate  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
CODE = recon.CODE_STATES


def bucket(key: str) -> str:
    return f"{zlib.crc32(key.encode()) & 0xff:02x}"


def neg_blosum_partner(a: str, rng) -> str:
    B = events.blosum62()
    opts = [b for b in AA if b != a and B[a][b] < 0]
    return str(rng.choice(opts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthoxml", default=str(ROOT / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz"))
    ap.add_argument("--guide", default=str(ROOT / "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk"))
    ap.add_argument("--n-ogs", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--fp-dim", type=int, default=32)
    ap.add_argument("--fp-shift", type=float, default=1.0)
    ap.add_argument("--out", default=str(ROOT / "results/partB/sim/toy_sogs/TOY/O1"))
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    if out.exists():
        shutil.rmtree(out)
    (out / "sog_states").mkdir(parents=True)
    (out / "replicate_states").mkdir(parents=True)
    trees = [t for t in recon.iter_orthoxml_trees(a.orthoxml, min_species=6) if t.n_leaves <= 300]
    pick = sorted(rng.choice(len(trees), size=min(a.n_ogs, len(trees)), replace=False).tolist())
    trees = [trees[i] for i in pick]
    st = recon.SpeciesTree.from_file(a.guide)
    all_species = st.tree.leaf_names()
    noise = simulate.NoiseModel()
    rows, reps, truth, sogs, ogs, resid, fps, repcos = [], [], [], [], [], [], [], []
    for t in trees:
        tix = events.TreeIndex(t)
        sp = t.leaf_species()
        names = t.leaf_names()
        elig = {k: simulate.eligible_branches(t, tix, k, min_species=3) for k in ("gain", "loss")}
        n_s = int(rng.integers(1, 4))
        og = t.tree_id
        ogs.append(dict(og_id=og, label=f"TOY {og}", n_species=len(set(sp)), n_members=t.n_leaves, n_sogs=n_s,
                        ref_accession=names[0], flags="toy"))
        for k in range(n_s):
            sog = f"{og}:g{k}"
            u = rng.random()
            scen = "vertical" if u < 0.7 else ("gain" if u < 0.85 else "loss")
            if scen != "vertical" and len(elig[scen]) == 0:
                scen = "vertical"
            tr, planted = simulate.simulate_truth(t, 1, scen, rng, elig.get(scen))
            rates = noise.rates(sp, rng)
            obs = simulate.observe(tr, rates, rng)[:, 0]
            obs2 = simulate.observe(tr, rates, rng)[:, 0]
            pl = int(planted[0])
            spec_node = -1
            spec_parent = -1
            # footprint residues + fingerprints
            L = 12
            anc = rng.choice(list(AA), size=L)
            R = np.tile(anc, (t.n_leaves, 1))
            mut = rng.random(R.shape) < 0.08
            R[mut] = rng.choice(list(AA), size=int(mut.sum()))
            base = np.abs(rng.normal(size=a.fp_dim))
            base /= np.linalg.norm(base)
            F = base[None, :] + rng.normal(scale=0.08, size=(t.n_leaves, a.fp_dim))
            vol = rng.normal(500, 80, size=t.n_leaves)
            plddt = rng.uniform(70, 100, size=t.n_leaves)
            src = rng.integers(0, 2, size=t.n_leaves)
            F += 0.0004 * (vol - 500)[:, None] * base[None, :]          # covariate effect to regress out
            if scen == "vertical" and rng.random() < 0.15:
                cand = [(v, c) for v, c, d, s in events.sister_pairs(t, obs, 2)
                        if len(d) >= 2]
                if cand:
                    v, c = cand[int(rng.integers(len(cand)))]
                    spec_node = c
                    spec_parent = v
                    lo, hi = t.leaf_lo[c], t.leaf_hi[c]
                    cols = rng.choice(L, size=3, replace=False)
                    for j in cols:
                        R[lo:hi, j] = neg_blosum_partner(anc[j], rng)
                        R[:lo, j] = anc[j]
                        R[hi:, j] = anc[j]
                    d = rng.normal(size=a.fp_dim)
                    d -= (d @ base) * base
                    d /= np.linalg.norm(d)
                    F[lo:hi] += a.fp_shift * d
            truth.append(dict(og_id=og, sog_id=sog, scenario=scen, planted_node=pl,
                              planted_name=t.node_name[pl] if pl >= 0 else "",
                              planted_clade=";".join(sorted(set(sp[i] for i in range(t.leaf_lo[pl], t.leaf_hi[pl]))))
                              if pl >= 0 else "",
                              specificity_node=spec_node, specificity_parent=spec_parent,
                              specificity_name=t.node_name[spec_node] if spec_node >= 0 else ""))
            cnt = {c: int((obs == recon.STATE_CODES[c]).sum()) for c in "PWAU"}
            present_sp = set(sp)
            g_sp = [s for s in all_species if s not in present_sp]
            sogs.append(dict(og_id=og, sog_id=sog, label=f"TOY {scen}", footprint_columns=list(range(L)),
                             footprint_ref=None, n_P=cnt["P"], n_W=cnt["W"], n_A=cnt["A"], n_U=cnt["U"],
                             n_G=len(g_sp), event="none", event_node=None, score=None, case_id=None))
            for i in range(t.n_leaves):
                tx = int(t.taxid[t.leaves[i]])
                st_ = CODE[obs[i]]
                for tbl, s_ in ((rows, st_), (reps, CODE[obs2[i]])):
                    tbl.append(dict(og_id=og, sog_id=sog, taxid=tx, accession=names[i], state=s_,
                                    site_id=f"{names[i]}:s0" if s_ in "PW" else None, pocket_id=None,
                                    jaccard=np.float32(0.5) if s_ in "PW" else None,
                                    footprint_plddt=np.float32(plddt[i] if s_ == "A" else plddt[i]),
                                    reason="" if s_ in "PA" else ("rescue:1/6" if s_ == "W" else "plddt<90")))
                resid.append(dict(sog_id=sog, accession=names[i], residues="".join(R[i])))
                if obs[i] in (recon.STATE_CODES["P"], recon.STATE_CODES["W"]):
                    fps.append(dict(sog_id=sog, accession=names[i], fingerprint=F[i].astype(np.float32).tolist(),
                                    site_volume_A3=float(vol[i]), footprint_plddt=float(plddt[i]),
                                    model_source=int(src[i])))
                    rep = F[i] + rng.normal(scale=0.08, size=a.fp_dim)
                    repcos.append(float(F[i] @ rep / np.linalg.norm(F[i]) / np.linalg.norm(rep)))
            for s_ in g_sp:
                rows.append(dict(og_id=og, sog_id=sog, taxid=recon.taxid_from_species_name(s_), accession=None,
                                 state="G", site_id=None, pocket_id=None, jaccard=None, footprint_plddt=None,
                                 reason=""))
    for name, tbl in (("sog_states", rows), ("replicate_states", reps)):
        df = pd.DataFrame(tbl)
        df["bucket"] = df["og_id"].map(bucket)
        for b, g in df.groupby("bucket"):
            g.drop(columns="bucket").sort_values(["og_id", "sog_id"]).to_parquet(
                out / name / f"bucket={b}.parquet", index=False)
    pd.DataFrame(sogs).to_parquet(out / "sogs.parquet", index=False)
    pd.DataFrame(ogs).to_parquet(out / "orthogroups.parquet", index=False)
    pd.DataFrame(truth).to_parquet(out / "truth.parquet", index=False)
    pd.DataFrame(resid).to_parquet(out / "footprint_residues.parquet", index=False)
    pd.DataFrame(fps).to_parquet(out / "fingerprints.parquet", index=False)
    pd.DataFrame({"cosine": repcos}).to_parquet(out / "replicate_cosines.parquet", index=False)
    shutil.copy(a.guide, out / "species_tree.nwk")
    pd.DataFrame(dict(taxid=[recon.taxid_from_species_name(s) for s in all_species], name=all_species)).to_parquet(
        out / "species.parquet", index=False)
    man = {"schema": "sogs-web-v1", "snapshot": "TOY", "orthology": "O1", "tau": None, "toy": True,
           "created_by": "scripts/partB/make_toy_sogs.py", "seed": a.seed, "n_ogs": len(trees),
           "n_sogs": len(sogs), "noise_model": {k: v for k, v in noise.__dict__.items() if k != "per_species"},
           "code_sha256": recon.code_sha256(recon, events, simulate, __file__)}
    (out / "manifest.json").write_text(json.dumps(man, indent=1))
    for f in ["truth.parquet", "sogs.parquet", "manifest.json", "fingerprints.parquet"]:
        recon.write_sidecar(out / f, "TOY", a.seed, {"script": "scripts/partB/make_toy_sogs.py"})
    tv = pd.DataFrame(truth)
    print(json.dumps({"n_ogs": len(trees), "n_sogs": len(sogs), "scenarios": tv["scenario"].value_counts().to_dict(),
                      "specificity_planted": int((tv["specificity_node"] >= 0).sum())}))


if __name__ == "__main__":
    main()
