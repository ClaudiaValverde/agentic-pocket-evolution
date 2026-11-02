#!/usr/bin/env python
"""False-A rate on known-present sites + rescue calibration (DESIGN §3.3 (ii), §12 round-2 item 7; prereg
states.rescue, states.measured_error_rates).

Known-present sites (experimental): LIGYSIS ligand sites transferred to AFDB v6 models (burrow af_set, UniProt
numbering), plus BioLiP2 (data/annotations/biolip_sites.parquet) and UniProt experimental binding sites
(data/annotations/uniprot_sites.parquet, ECO:0000269/ECO:0007744) when those files exist. Restricted to tier-1/2
proteins of the snapshot (single-fragment F1 models, so UniProt numbering = model numbering).

Per known site K (residue set):
  region pLDDT       mean and frac(<70) of K in the atlas v6 model (CA B-factors)
  A-eligible         primary: mean >= 90 and frac(<70) <= 0.2 (prereg states.A_primary); sensitivity: mean >= 70.
                     Same-protein evaluation, so aligned-fraction / length / truncation criteria are trivially met;
                     signal-peptide/TM overlap is applied only if data/annotations/uniprot_features.parquet exists.
  detected           an atlas site of the protein covers >= 50 % of K (prereg benchmarks.detectability_preflight)
  false-A            A-eligible and not detected. The state caller is blind to labels (it never reads them).
Rescue: the six settings (pocketevo.rescue) on all false-A cases (or a stratified sample of --max-cases) and on
true-negative controls: random surface patches of the same size in the same proteins that pass the same A
criteria and overlap no atlas site and no known site (seeded). k chosen by the prereg rule on the TUNE half.

  python scripts/observation/false_a.py --workers 6
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import yaml

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / "src"))
sys.path.insert(0, str(P / "scripts"))
from pocketevo import observation as O  # noqa: E402
from pocketevo import rescue as R  # noqa: E402
from pocketevo import sites as SI  # noqa: E402

ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")
AFDB_V6 = Path("/cluster/work/beltrao/shared/alphafold/v6")
SNAP = "20261030T085136Z"
AFSET = P / "external/burrow/benchmarks/af_set/sites.tsv.gz"
COFACTORS = {"HEM", "HEC", "HEA", "HEB", "SRM", "NAD", "NAP", "NDP", "NAI", "NAJ", "FAD", "FMN", "FDA", "SAM", "SAH",
             "COA", "ACO", "PLP", "PMP", "TPP", "TDP", "BTN", "B12", "CNC", "COB", "F43", "MGD", "PQQ", "TOP", "LPA",
             "GSH", "ATP", "ADP", "AMP", "ANP", "GTP", "GDP", "GNP", "UDP", "CTP", "UTP", "H4B", "THF", "FOL",
             "MTE", "BCL", "CLA", "CHL", "RET", "UQ1", "UQ2", "PQN", "MQ7", "SF4", "FES", "F3S"}
METAL_WORDS = ("metal", "ion")


# ------------------------------------------------------------------------------------------------ known sites
COFACTOR_WORDS = ("heme", "haem", "fad", "fmn", "nad", "coenzyme", "pyridoxal", "thiamine", "biotin", "cobalamin",
                  "s-adenosyl", "flavin", "lipoate", "molybdopterin", "pqq", "chlorophyll", "iron-sulfur", "[4fe-4s]",
                  "[2fe-2s]", "atp", "adp", "gtp", "gdp", "ubiquinone", "menaquinone", "tetrahydro")


def _dedupe(df: pd.DataFrame, jmax: float = 0.5) -> pd.DataFrame:
    """Per accession x ligand class: greedy (largest first) keep a site if residue Jaccard < jmax with every
    kept site (BioLiP lists the same site once per PDB entry / chain)."""
    keep = []
    for _, g in df.groupby(["accession", "ligand_class"], sort=False):
        kept = []
        for idx, res in sorted(zip(g.index, g.residues), key=lambda t: -len(t[1])):
            rs = set(res)
            if all(len(rs & k) / len(rs | k) < jmax for k in kept):
                kept.append(rs)
                keep.append(idx)
    return df.loc[sorted(keep)]


def _uniprot_class(lig: str) -> str:
    l = lig.lower()
    if re.search(r"\((\d)?\+\)$|\(\d\+\)", l) and len(l) < 16:
        return "metal"
    return "cofactor" if any(w in l for w in COFACTOR_WORDS) else "small_molecule"


def known_sites(tier12: set) -> tuple[pd.DataFrame, dict]:
    rows, cnt = [], {}
    a = pd.read_csv(AFSET, sep="\t", usecols=["structure_id", "site_id", "residues", "uniprot", "ion_prop"])
    for r in a.itertuples():
        rows.append(dict(source="LIGYSIS_af_set", accession=r.uniprot, known_id=f"L:{r.site_id}",
                         residues=sorted(O.parse_residues(r.residues)),
                         ligand_class="metal" if float(r.ion_prop) >= 0.999 else "small_molecule", ligand=None))
    L = pd.DataFrame(rows)
    parts = [L]
    bl = P / "data/annotations/biolip_sites.parquet"
    if bl.exists():
        b = pd.read_parquet(bl, columns=["accession", "pdb_id", "chain", "site_code", "ligand", "ligand_class",
                                         "residues_uniprot", "n_site_residues", "n_mapped"])
        cnt["biolip_rows"] = int(len(b))
        b = b[b.accession.isin(tier12)]
        cnt["biolip_rows_tier12"] = int(len(b))
        cnt["biolip_excluded_class"] = b.ligand_class.value_counts().to_dict()
        b = b[b.ligand_class.isin(["regular", "metal"]) & (b.n_mapped >= 0.8 * b.n_site_residues)]
        b = b.assign(residues=[sorted({int(x) for x in r}) for r in b.residues_uniprot],
                     ligand_class=np.where(b.ligand_class == "metal", "metal",
                                           np.where(b.ligand.isin(COFACTORS), "cofactor", "small_molecule")),
                     known_id="B:" + b.pdb_id + b.chain + ":" + b.site_code + ":" + b.ligand.astype(str),
                     source="BioLiP2")
        b = b[b.residues.map(len) >= 3].drop_duplicates("known_id")
        cnt["biolip_sites_before_dedupe"] = int(len(b))
        b = _dedupe(b.reset_index(drop=True))
        cnt["biolip_sites_after_dedupe"] = int(len(b))
        parts.append(b[["source", "accession", "known_id", "residues", "ligand_class", "ligand"]])
    us = P / "data/annotations/uniprot_site_features.parquet"
    if us.exists():
        u = pd.read_parquet(us, columns=["uniprot_acc", "type", "start", "end", "ligand", "eco_codes", "experimental"])
        u = u[(u.type == "Binding site") & u.uniprot_acc.isin(tier12) & u.ligand.notna()
              & u.eco_codes.astype(str).str.contains("ECO:0000269|ECO:0007744")]
        ur = []
        for (acc, lig), g in u.groupby(["uniprot_acc", "ligand"]):
            res = sorted({int(x) for s0, e0 in zip(g.start, g.end) if pd.notna(s0) and pd.notna(e0)
                          for x in range(int(s0), int(e0) + 1)})
            ur.append(dict(source="UniProt_ECO", accession=acc, known_id=f"U:{acc}:{lig}", residues=res,
                           ligand_class=_uniprot_class(str(lig)), ligand=lig))
        ur = pd.DataFrame(ur)
        cnt["uniprot_sites"] = int(len(ur))
        parts.append(ur)
    K = pd.concat(parts, ignore_index=True)
    K = K[K.accession.isin(tier12) & (K.residues.map(len) >= 3)]
    return K.reset_index(drop=True), cnt


def sp_tm_ranges(accs: set) -> dict:
    """accession -> set of residues in UniProt SIGNAL/TRANSMEM features (prereg overlap max 0)."""
    f = P / "data/annotations/uniprot_proteins.parquet"
    if not f.exists():
        return {}, {}
    p = pd.read_parquet(f, columns=["accession", "signal", "transmem", "seq_md5"])
    p = p[p.accession.isin(accs)]
    out, md5 = {}, {}
    for r in p.itertuples():
        rs = set()
        for txt in (r.signal, r.transmem):
            for m in re.finditer(r"(\d+)\s*-\s*(\d+)", str(txt or "")):
                rs |= set(range(int(m.group(1)), int(m.group(2)) + 1))
        out[r.accession] = rs
        md5[r.accession] = str(r.seq_md5 or "").lower()
    return out, md5


# ------------------------------------------------------------------------------------------------ helpers
def _model_info(path):
    try:
        seq, num, pl = O.read_model(path)
        import gemmi
        st = gemmi.read_structure(str(path))
        ch = st[0][0]
        xyz, burial_num = [], []
        for res in ch:
            a = res.find_atom("CB", "*") or res.find_atom("CA", "*")
            if a is None:
                continue
            xyz.append([a.pos.x, a.pos.y, a.pos.z])
            burial_num.append(res.seqid.num)
        X = np.asarray(xyz, float)
        D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
        nb = (D < 10.0).sum(1) - 1
        return dict(num=num, plddt=pl, seq=seq, xyz=X, nb=nb, err=None)
    except Exception as e:  # noqa: BLE001
        return dict(err=f"{type(e).__name__}: {e}")


def control_patch(info, size, forbid: set, rng, crit_mean=90.0, crit_frac=0.2, max_nb=18, tries=200):
    """Random surface patch of `size` residues (seed residue + nearest surface residues by CB distance) that
    passes the A pLDDT criteria and avoids `forbid`. None if not found."""
    num, pl, X, nb = info["num"], info["plddt"], info["xyz"], info["nb"]
    surf = np.flatnonzero((nb <= max_nb) & ~np.isin(num, list(forbid)))
    if len(surf) < size:
        return None
    for _ in range(tries):
        s = rng.choice(surf)
        d = np.sqrt(((X[surf] - X[s]) ** 2).sum(-1))
        pick = surf[np.argsort(d)[:size]]
        if d[np.argsort(d)[size - 1]] > 16.0:
            continue
        v = pl[pick]
        if v.mean() >= crit_mean and (v < 70).mean() <= crit_frac:
            return sorted(int(x) for x in num[pick])
    return None


def git_commit():
    r = subprocess.run(["git", "-C", str(P), "rev-parse", "HEAD"], capture_output=True, text=True)
    return r.stdout.strip() or None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/tmp/pe_big/obs/falseA")
    ap.add_argument("--out", default=str(P / f"results/observation/{SNAP}"))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--max-cases", type=int, default=1500)
    ap.add_argument("--tau", type=float, default=0.3, help="rescue hit rule: pocket-footprint residue Jaccard")
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--skip-rescue", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    T = {}
    W, out = Path(a.work), Path(a.out)
    W.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    pre = yaml.safe_load(open(P / "config/prereg.yaml"))
    seed = int(pre["meta"]["seeds"][0])
    Ap = pre["states"]["A_primary"]
    A_MEAN, A_FRAC = float(Ap["footprint_mean_plddt"]["value"]), float(Ap["footprint_frac_residues_below_70"]["value"])
    A_SENS = float(pre["states"]["A_sensitivity_stratum"]["footprint_mean_plddt"])
    OVL = float(pre["benchmarks"]["detectability_preflight"]["overlap_min"])
    TARGET = float(pre["states"]["rescue"]["false_A_target"]["value"])
    KGRID = [int(x) for x in pre["states"]["rescue"]["k_settings_to_rescue"]["grid"]]

    # ---- tier-1/2 structures (single-fragment, ok)
    man = pd.read_parquet(P / f"data/snapshot/{SNAP}/manifest.parquet")
    man = man[man.tier.isin([1, 2])]
    wl = pd.read_parquet(ATLAS / "worklist.parquet").set_index("unit_id")
    St = pd.concat([pd.read_parquet(ATLAS / p, columns=["unit_id", "tier", "structure_stem", "accession", "fragment",
                                                        "taxid", "organism", "n_res", "status", "n_pockets_kept"])
                    for p in man.structures_path], ignore_index=True)
    nfrag = St.groupby("accession").fragment.max()
    St = St[(St.status == "ok") & (St.fragment == 1) & St.accession.map(nfrag).eq(1)].drop_duplicates("accession")
    St["tar"] = St.unit_id.map(wl.source)
    K, kcnt = known_sites(set(St.accession))
    T["known_site_counts"] = kcnt
    K = K.merge(St[["accession", "unit_id", "tier", "structure_stem", "taxid", "organism", "tar"]], on="accession")
    T["known_sites_loaded_s"] = time.time() - t0
    # ---- models
    v6d = W / "v6"
    v6d.mkdir(exist_ok=True)
    from extract_afdb_member import extract
    pro = K.drop_duplicates("accession")
    for tar, g in pro.groupby("tar"):
        need = [f"{s}.cif.gz" for s in g.structure_stem if not (v6d / f"{s}.cif.gz").exists()]
        if need:
            extract(Path(tar) if "/" in str(tar) else AFDB_V6 / tar, need, v6d)
    paths = [str(v6d / f"{s}.cif.gz") for s in pro.structure_stem]
    with ProcessPoolExecutor(a.workers) as ex:
        infos = dict(zip(pro.accession, ex.map(_model_info, paths, chunksize=8)))
    T["models_s"] = time.time() - t0
    # ---- atlas sites of these proteins
    cols = SI.POCKET_COLUMNS
    aset = pa.array(sorted(pro.accession))
    parts = []
    for u in sorted(pro.unit_id.unique()):
        t = pq.read_table(ATLAS / man.set_index("unit_id").at[u, "pockets_path"], columns=cols)
        parts.append(t.filter(pc.is_in(t["accession"], aset)))
    sites, _, _ = SI.build_sites(pa.concat_tables(parts, promote_options="permissive"), with_fingerprints=False)
    sites = sites.to_pandas()
    asites = {k: [frozenset(x) for x in g.residues] for k, g in sites.groupby("accession")}
    # ---- per known site
    sptm, umd5 = sp_tm_ranges(set(pro.accession))
    rec = []
    for r in K.itertuples():
        inf = infos[r.accession]
        if inf.get("err"):
            rec.append(dict(known_id=r.known_id, skip="model_error"))
            continue
        ks = frozenset(r.residues)
        seq_ok = (r.source == "LIGYSIS_af_set" or r.accession not in umd5
                  or umd5[r.accession] == hashlib.md5(inf["seq"].encode()).hexdigest())
        if not seq_ok:
            rec.append(dict(known_id=r.known_id, skip="model_seq_ne_uniprot"))
            continue
        sp_ov = bool(ks & sptm.get(r.accession, set()))
        m = np.isin(inf["num"], list(ks))
        frac_found = m.sum() / len(ks)
        v = inf["plddt"][m]
        cov = max([O.coverage(ks, s) for s in asites.get(r.accession, [])], default=0.0)
        jac = max([O.jaccard(ks, s) for s in asites.get(r.accession, [])], default=0.0)
        rec.append(dict(known_id=r.known_id, skip=None if frac_found >= 0.8 else "residues_not_in_model",
                        sp_tm_overlap=sp_ov,
                        frac_found=frac_found, region_plddt=float(v.mean()) if len(v) else np.nan,
                        region_frac_lt70=float((v < 70).mean()) if len(v) else np.nan,
                        best_atlas_cov=cov, best_atlas_jaccard=jac, n_atlas_sites=len(asites.get(r.accession, []))))
    K = K.merge(pd.DataFrame(rec), on="known_id", how="left")
    K = K.drop_duplicates("known_id")
    K["detected"] = K.best_atlas_cov >= OVL
    K["sp_tm_overlap"] = K.sp_tm_overlap.fillna(False).astype(bool)
    K["A_primary"] = (K.region_plddt >= A_MEAN) & (K.region_frac_lt70 <= A_FRAC) & K.skip.isna() & ~K.sp_tm_overlap
    K["A_sens"] = (K.region_plddt >= A_SENS) & K.skip.isna() & ~K.sp_tm_overlap
    K["band"] = [O.plddt_band(x) for x in K.region_plddt]
    K["false_A"] = K.A_primary & ~K.detected
    K["size_bin"] = pd.cut(K.residues.map(len), [0, 5, 10, 20, 1000], labels=["3-5", "6-10", "11-20", ">20"]).astype(str)
    K["false_A_sens"] = K.A_sens & ~K.detected
    # tune/test split by root HOG (O1), control families forced to test
    kp = W / "hog_keys.parquet"
    if kp.exists():
        hk = pd.read_parquet(kp)
    else:
        d = O.root_hog_keys(P / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz")
        hk = pd.DataFrame({"accession": list(d), "key": list(d.values())})
        hk.to_parquet(kp, index=False)
    hkd = dict(zip(hk.accession, hk.key))
    forced = {hkd.get(x, x) for x in pre["tune_test"]["forced_test_half"]["root_hogs_containing"]}
    K["split_key"] = [hkd.get(x, x) for x in K.accession]
    K["half"] = ["test" if k in forced else O.split_half(k) for k in K.split_key]
    K["in_O1_hog"] = K.accession.isin(hkd)
    T["false_A_table_s"] = time.time() - t0
    # ---- rates
    def rates(df, flag="A_primary"):
        x = df[df[flag]]
        return dict(n_known=int(len(df)), n_A_eligible=int(len(x)), n_false_A=int((~x.detected).sum()),
                    false_A=float((~x.detected).mean()) if len(x) else None)
    strata = []
    for flag in ("A_primary", "A_sens"):
        for by in (None, "source", "band", "ligand_class", "size_bin", "tier", "half", "organism"):
            groups = [("all", K)] if by is None else list(K.groupby(by))
            for k, g in groups:
                strata.append(dict(criterion=flag, stratum=by or "all", value=str(k), **rates(g, flag)))
    Srt = pd.DataFrame(strata)
    Srt.to_csv(out / "false_a_strata.csv", index=False)
    # ---- rescue
    res = dict()
    if not a.skip_rescue:
        rng = np.random.default_rng(seed)
        cases = K[K.false_A | K.false_A_sens].copy()
        if len(cases) > a.max_cases:   # stratified by half x band x ligand class
            cases["st"] = cases.source + "|" + cases.half + "|" + cases.band + "|" + cases.ligand_class
            cases = cases.groupby("st", group_keys=False).apply(
                lambda g: g.sample(n=max(1, int(round(a.max_cases * len(g) / len(cases)))), random_state=seed))
        # controls: one per case, same protein, same size, A-primary pLDDT criteria, no overlap with atlas/known
        known_by_acc = K.groupby("accession").residues.apply(lambda s: set().union(*map(set, s)))
        ctrl = []
        for r in cases.itertuples():
            inf = infos[r.accession]
            forbid = set(known_by_acc.get(r.accession, set()))
            for s in asites.get(r.accession, []):
                forbid |= set(s)
            patch = control_patch(inf, len(r.residues), forbid, rng, A_MEAN if r.A_primary else A_SENS,
                                  A_FRAC if r.A_primary else 1.0)
            if patch is not None:
                ctrl.append(dict(known_id=f"CTRL:{r.known_id}", accession=r.accession, residues=patch,
                                 half=r.half, band=r.band, ligand_class="control", A_primary=r.A_primary,
                                 structure_stem=r.structure_stem))
        ctrl = pd.DataFrame(ctrl)
        models = sorted({str(v6d / f"{s}.cif.gz") for s in pd.concat([cases.structure_stem, ctrl.structure_stem])})
        t1 = time.time()
        outs = R.run_settings(models, W / "rescue", workers=a.workers)
        T["rescue_burrow_s"] = time.time() - t1
        T["rescue_n_models"] = len(models)
        hrows = []
        for kind, df in (("false_A_case", cases), ("control", ctrl)):
            for r in df.itertuples():
                for rule in ("jaccard", "coverage"):
                    h = R.hits_from_outputs(outs, r.structure_stem, r.residues, a.tau, rule, a.min_cov)
                    hrows.append(dict(kind=kind, rule=rule, known_id=r.known_id, accession=r.accession, half=r.half,
                                      band=r.band, ligand_class=r.ligand_class, A_primary=bool(r.A_primary),
                                      n_hits=R.n_hits(h), **{f"hit_{k}": v["hit"] for k, v in h.items()},
                                      **{f"jac_{k}": v["jaccard"] for k, v in h.items()}))
        H = pd.DataFrame(hrows)
        H.to_parquet(out / "rescue_hits.parquet", index=False)
        # false-A after rescue at k: falseA_k = falseA_0 * P(n_hits < k | false-A case) (case sample weights)
        table = []
        for rule in ("jaccard", "coverage"):
            for crit in ("A_primary", "A_sens"):
                fa_flag = "false_A" if crit == "A_primary" else "false_A_sens"
                for half in ("tune", "test", "all"):
                    Kh = K if half == "all" else K[K.half == half]
                    elig = Kh[Kh[crit]]
                    fa0 = float(Kh[fa_flag].sum() / max(1, len(elig)))
                    ids = set(Kh.loc[Kh[fa_flag], "known_id"])
                    hc = H[(H.kind == "false_A_case") & (H.rule == rule) & H.known_id.isin(ids)]
                    hn = H[(H.kind == "control") & (H.rule == rule) & H.known_id.str[5:].isin(ids)]
                    for k in KGRID:
                        still = float((hc.n_hits < k).mean()) if len(hc) else np.nan
                        table.append(dict(rule=rule, criterion=crit, half=half, k=k, n_A_eligible=int(len(elig)),
                                          false_A_before=fa0, n_cases_rescue_run=int(len(hc)),
                                          frac_cases_still_A=still, false_A_after=fa0 * still if len(hc) else np.nan,
                                          n_controls=int(len(hn)),
                                          control_rescued_frac=float((hn.n_hits >= k).mean()) if len(hn) else np.nan))
        RT = pd.DataFrame(table)
        RT.to_csv(out / "rescue_k_table.csv", index=False)
        tune = RT[(RT.rule == "jaccard") & (RT.criterion == "A_primary") & (RT.half == "tune")].sort_values("k")
        meet = tune[tune.false_A_after <= TARGET]
        k_lit = int(meet.k.min()) if len(meet) else None
        k_max = int(meet.k.max()) if len(meet) else None
        res = dict(rule_prereg="smallest k meeting false_A_target on the known-present set (TUNE half)",
                   false_A_target=TARGET, grid=KGRID, hit_rule=f"residue Jaccard >= {a.tau} (pocket vs footprint)",
                   k=k_lit, k_largest_meeting_target=k_max, target_met=bool(len(meet)),
                   k_min_false_A=int(tune.sort_values(["false_A_after", "k"]).k.iloc[0]),
                   if_target_not_met=("no k in the grid meets the target on the TUNE half: `k` is null (prereg "
                                      "gives no fallback; PI decision). k_min_false_A is the k with the lowest "
                                      "false-A after rescue (k=1, most conservative toward absence claims); "
                                      "residual false-A at that k must be propagated into the simulations."),
                   tune_table=tune[["k", "false_A_before", "frac_cases_still_A", "false_A_after",
                                    "control_rescued_frac", "n_cases_rescue_run", "n_controls"]].to_dict("records"),
                   test_table=RT[(RT.rule == "jaccard") & (RT.criterion == "A_primary") & (RT.half == "test")]
                   [["k", "false_A_before", "false_A_after", "control_rescued_frac"]].to_dict("records"),
                   note=("false-A rate decreases monotonically as k decreases, so the literal prereg rule picks the "
                         "smallest k that meets the target (k=1 whenever any k does). k_largest_meeting_target is "
                         "reported for the PI (stricter rescue, fewer true absences turned into W); the prereg value "
                         "is `k`."),
                   n_cases=int(len(cases)), n_controls=int(len(ctrl)),
                   n_controls_not_found=int(len(cases) - len(ctrl)))
    K.drop(columns=["tar"]).to_parquet(out / "false_a_sites.parquet", index=False)
    head = dict(snapshot_id=SNAP, prereg_sha256=hashlib.sha256((P / "config/prereg.yaml").read_bytes()).hexdigest(),
                git_commit=git_commit(), seed=seed, sources=K.source.value_counts().to_dict(),
                n_known=int(len(K)), n_proteins=int(K.accession.nunique()),
                n_skipped=K.skip.value_counts().to_dict(),
                A_primary=rates(K, "A_primary"), A_sens=rates(K, "A_sens"),
                by_band={b: rates(g) for b, g in K.groupby("band")},
                by_ligand_class={b: rates(g) for b, g in K.groupby("ligand_class")},
                by_size_bin={b: rates(g) for b, g in K.groupby("size_bin")},
                by_source={b: rates(g) for b, g in K.groupby("source")},
                false_A_primary_by_overlap_threshold={str(t): float((K[K.A_primary].best_atlas_cov < t).mean())
                                                      for t in (1e-9, 0.3, 0.5, 0.7)},
                by_tier={str(b): rates(g) for b, g in K.groupby("tier")},
                by_half={b: rates(g) for b, g in K.groupby("half")},
                n_in_O1_hog=int(K.in_O1_hog.sum()), timings=T, runtime_s=round(time.time() - t0, 1))
    (out / "false_a.json").write_text(json.dumps(head, indent=1, default=float))
    if res:
        res.update(snapshot_id=SNAP, prereg_sha256=head["prereg_sha256"], git_commit=head["git_commit"], seed=seed)
        (out / "rescue_k.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(head, indent=1, default=float))
    if res:
        print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
