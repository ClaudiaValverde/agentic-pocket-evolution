#!/usr/bin/env python3
"""Co-folding calibration panel: freeze the readout thresholds before any case is folded.

Positives are holo PDB complexes from BioLiP2 for proteins that already have a cached
AF3 MSA. Negatives are, per decision D23, BOTH property-matched decoy ligands AND the
site-ablated variant of the same protein with its cognate ligand, because a case compares
protein variants, not only ligands.

SELECTION RULE (deterministic; no hand-picking; dumped verbatim to selection_rule.json)
  R1  Candidates = UniProt accessions with a cached AF3 data JSON (human 25.05 or yeast
      25.06). The cached query sequence is the sequence folded.
  R2  BioLiP2 sites for those accessions, keeping rows with
        ligand_class in {regular, metal};  ligand is a 1-3 character CCD code present in
        the AF3 CCD;  resolution <= 2.5 A;  aln_coverage >= 0.95;  aln_identity >= 0.98;
        n_mapped == n_site_residues (every site residue mapped to UniProt);
        n_site_residues >= 5;  80 <= cached sequence length <= 480 (keeps the 512 token
        bucket, so per-job runtime stays near 2 min).
  R3  Drop the four accessions used as positive controls in docs/COFOLD_CONTROLS.md
      (P24941, P00918, P06737, P00330) so the panel is independent of them.
  R4  Ligand class from the CCD code: metal (single-metal codes), nucleotide, cofactor
      (fixed code lists below), else drug_like when 10 <= heavy atoms <= 60.
  R5  Collapse to one row per (accession, ligand): best resolution, ties by pdb_id asc.
  R6  Keep one row per accession: best resolution, ties by pdb_id asc.
  R7  Within each ligand class sort by (resolution asc, pdb_id asc, accession asc) and
      take the first --per-class rows.

DECOY RULE (prereg cofold.decoy_matching: MW, logP, HBD, HBA, rotatable bonds)
  Properties are computed with RDKit from the CCD SMILES. For each positive, the decoy
  pool is every ligand OF THE SAME LIGAND CLASS that survives R2/R4 anywhere in BioLiP2,
  minus every ligand BioLiP2 records as binding that accession in any structure. Within
  that pool the properties are z-scored and the two nearest neighbours of the cognate
  ligand are taken as decoys. Matching is within class on purpose: across classes the
  nearest neighbour of a Ca site is an organic molecule, which is a gross mismatch
  control rather than a decoy, and would flatter the readouts.

SITE
  BioLiP2 `residues_uniprot` (already UniProt-numbered and identity-checked). BioLiP's
  cutoff (sum of van der Waals radii + 0.5 A) is slightly tighter than the 4.5 A used in
  docs/COFOLD_CONTROLS.md; the panel is internally consistent, which is what calibration
  needs. The same positions define the Gly ablation.

MSA / templates
  Cached AF3 jackhmmer MSA, `--no-templates`, 3 seeds x 5 samples. Ablations reuse the
  cached MSA with the query row replaced, holding the MSA constant (same reasoning as
  COFOLD_CONTROLS.md Sec.1).

Usage
  calibration_panel.py --select            # write panel.tsv + selection_rule.json
  calibration_panel.py --build --run       # build inputs, fold on GPUs 0-2
  calibration_panel.py --score             # readouts, AUCs, threshold.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cofold_readouts import summarise_job  # noqa: E402

P = Path(__file__).resolve().parents[2]  # repo root
AF3_PY = os.environ.get("AF3_ENV", "/tmp/af3_env") + "/bin/python"
MK = P / "scripts/cofold/make_af3_input.py"
RUN = P / "scripts/cofold/run_af3.sh"
BIG = Path("/tmp/pe_big/cofold/calibration")
RES = P / "results/cofold/calibration"
CACHE_DIRS = ["/cluster/project/beltrao/shared/25.05_alphafold3_msas_human",
              "/cluster/project/beltrao/shared/25.06_alphafold3_msas_yeast"]
CONTROL_ACCS = {"P24941", "P00918", "P06737", "P00330"}

METALS = {"ZN", "MG", "CA", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO", "NA", "K",
          "CD", "HG", "3CO", "4MO", "MO", "W", "SR", "BA", "CS", "RB", "PB", "PT"}
NUCLEOTIDES = {"ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", "UTP", "UDP", "TTP",
               "ANP", "ACP", "AGS", "GNP", "GSP", "GCP", "APC", "ADX", "5GP", "ATR"}
COFACTORS = {"NAD", "NAI", "NAP", "NDP", "FAD", "FMN", "SAM", "SAH", "PLP", "HEM",
             "HEC", "COA", "ACO", "TPP", "TDP", "BTN", "B12", "MQ7", "UQ1", "F43",
             "THG", "THF", "MTE", "PQQ", "LPA", "GSH"}


def log(m: str) -> None:
    print(f"[calibration] {m}", file=sys.stderr, flush=True)


def cached_seq(acc: str) -> str | None:
    import gzip
    for d in CACHE_DIRS:
        p = Path(d) / f"{acc.lower()}_data.json.gz"
        if p.exists():
            with gzip.open(p, "rt") as f:
                dd = json.load(f)
            pr = [s["protein"] for s in dd["sequences"] if "protein" in s]
            return pr[0]["sequence"] if len(pr) == 1 else None
    return None


def cached_accessions() -> set[str]:
    out = set()
    for d in CACHE_DIRS:
        for p in Path(d).glob("*_data.json.gz"):
            out.add(p.name.split("_data.json.gz")[0].upper())
    return out


def ccd_props() -> dict:
    """CCD code -> (MW, logP, HBD, HBA, rotatable bonds, n_heavy); RDKit, from CCD SMILES."""
    from alphafold3.constants.chemical_components import Ccd
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    ccd = Ccd()
    return ccd, Chem, Descriptors, Lipinski, rdMolDescriptors


def props_for(code: str, ccd, Chem, Descriptors, Lipinski, rdMolDescriptors):
    try:
        d = ccd[code]
    except Exception:  # noqa: BLE001
        return None
    smi = dict(zip(d.get("_pdbx_chem_comp_descriptor.type", []),
                   d.get("_pdbx_chem_comp_descriptor.descriptor", []))).get("SMILES")
    if not smi:
        return None
    m = Chem.MolFromSmiles(smi, sanitize=True)
    if m is None:
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            return None
        try:
            Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL
                             ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:  # noqa: BLE001
            return None
    try:
        return dict(mw=float(Descriptors.MolWt(m)), logp=float(Descriptors.MolLogP(m)),
                    hbd=int(Lipinski.NumHDonors(m)), hba=int(Lipinski.NumHAcceptors(m)),
                    rotb=int(rdMolDescriptors.CalcNumRotatableBonds(m)),
                    n_heavy=int(m.GetNumHeavyAtoms()))
    except Exception:  # noqa: BLE001
        return None


def tanimoto(a: str, b: str, ccd, Chem, *_rest) -> float | None:
    """ECFP4 Tanimoto between two CCD codes, from their CCD SMILES."""
    from rdkit.Chem import rdFingerprintGenerator
    ms = []
    for code in (a, b):
        d = ccd[code]
        smi = dict(zip(d.get("_pdbx_chem_comp_descriptor.type", []),
                       d.get("_pdbx_chem_comp_descriptor.descriptor", []))).get("SMILES")
        m = Chem.MolFromSmiles(smi) if smi else None
        if m is None:
            return None
        ms.append(m)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    from rdkit import DataStructs
    return float(DataStructs.TanimotoSimilarity(gen.GetFingerprint(ms[0]),
                                                gen.GetFingerprint(ms[1])))


MAX_DECOY_TANIMOTO = 0.4


def classify(code: str, n_heavy: int) -> str | None:
    if code in METALS:
        return "metal"
    if code in NUCLEOTIDES:
        return "nucleotide"
    if code in COFACTORS:
        return "cofactor"
    if 10 <= n_heavy <= 60:
        return "drug_like"
    return None


def select(per_class: int) -> None:
    import pandas as pd
    import numpy as np
    RES.mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet(P / "data/annotations/biolip_sites.parquet")
    n0 = len(df)
    accs = cached_accessions()
    df = df[df.uniprot_acc.isin(accs) & ~df.uniprot_acc.isin(CONTROL_ACCS)]
    df = df[df.ligand_class.isin(["regular", "metal"])]
    df = df[df.ligand.str.len().le(3) & df.ligand.str.match(r"^[A-Z0-9]+$")]
    df = df[(df.resolution <= 2.5) & (df.aln_coverage >= 0.95) & (df.aln_identity >= 0.98)]
    df = df[(df.n_mapped == df.n_site_residues) & (df.n_site_residues >= 5)]
    log(f"R1-R3: {n0} -> {len(df)} BioLiP rows")

    # Length comes from uniprot_proteins.parquet, not from the cached AF3 JSONs: each of
    # those is ~30 MB gzipped and there are thousands of candidates. The cached sequence
    # itself is read only for the few proteins that survive selection.
    up = pd.read_parquet(P / "data/annotations/uniprot_proteins.parquet",
                         columns=["uniprot_acc", "length", "fragment"])
    up = up.drop_duplicates("uniprot_acc").set_index("uniprot_acc")
    df = df[df.uniprot_acc.isin(up.index)]
    df = df[~df.uniprot_acc.map(up.fragment).fillna(False).astype(bool)]
    df["seq_len"] = df.uniprot_acc.map(up.length).astype(int)
    df = df[(df.seq_len >= 80) & (df.seq_len <= 480)]
    log(f"R2 length filter: -> {len(df)} rows, {df.uniprot_acc.nunique()} accessions")

    ccdargs = ccd_props()
    pr = {c: props_for(c, *ccdargs) for c in df.ligand.unique()}
    df = df[df.ligand.map(lambda c: pr.get(c) is not None)]
    df["n_heavy"] = df.ligand.map(lambda c: pr[c]["n_heavy"])
    df["lig_class"] = [classify(c, h) for c, h in zip(df.ligand, df.n_heavy)]
    df = df[df.lig_class.notna()]
    log(f"R4: -> {len(df)} rows; classes {df.lig_class.value_counts().to_dict()}")

    df = df.sort_values(["resolution", "pdb_id", "uniprot_acc"])
    df = df.drop_duplicates(["uniprot_acc", "ligand"], keep="first")          # R5
    df = df.drop_duplicates(["uniprot_acc"], keep="first")                     # R6
    df = df.sort_values(["lig_class", "resolution", "pdb_id", "uniprot_acc"])
    panel = pd.concat([g.head(per_class) for _, g in df.groupby("lig_class", sort=True)],
                      ignore_index=True)                                       # R7
    log(f"R5-R7: panel of {len(panel)}; {panel.lig_class.value_counts().to_dict()}")

    # ---- decoys: 2 nearest ligands OF THE SAME LIGAND CLASS in z-scored 5-D property
    # space. Matching within class is what "property-matched" has to mean here: across
    # classes the nearest neighbour of a Ca site is an organic molecule, which is a
    # mismatch control, not a decoy. The pool is every ligand of that class surviving
    # R2/R4 anywhere in BioLiP (not just the panel), so it is large and real.
    keys = ["mw", "logp", "hbd", "hba", "rotb"]
    pools = {cls: sorted(g.ligand.unique()) for cls, g in df.groupby("lig_class")}
    zof, zpool = {}, {}
    for cls, pool_c in pools.items():
        M = np.array([[pr[c][k] for k in keys] for c in pool_c], float)
        Z = (M - M.mean(0)) / np.where(M.std(0) == 0, 1, M.std(0))
        zpool[cls] = pool_c
        zof[cls] = {c: Z[i] for i, c in enumerate(pool_c)}
    pool = sorted({c for v in pools.values() for c in v})
    binds = (pd.read_parquet(P / "data/annotations/biolip_sites.parquet",
                             columns=["uniprot_acc", "ligand"])
             .groupby("uniprot_acc").ligand.apply(set).to_dict())
    rows = []
    for _, r in panel.iterrows():
        own = binds.get(r.uniprot_acc, set())
        cls = r.lig_class
        cands = [(float(np.linalg.norm(zof[cls][c] - zof[cls][r.ligand])), c)
                 for c in zpool[cls] if c != r.ligand and c not in own]
        cands.sort()
        if len(cands) < 2:
            log(f"{r.uniprot_acc}: only {len(cands)} in-class decoys available")
            continue
        decoys = [c for _, c in cands[:2]]
        rows.append({
            "acc": r.uniprot_acc, "pdb_id": r.pdb_id, "chain": r.chain,
            "resolution": float(r.resolution), "ligand": r.ligand,
            "lig_class": r.lig_class, "lig_name": str(r.ligand_name)[:60],
            "n_heavy": int(r.n_heavy), "seq_len": int(r.seq_len),
            "n_site": int(r.n_site_residues),
            "site": ",".join(str(x) for x in r.residues_uniprot),
            "decoy1": decoys[0], "decoy2": decoys[1],
            "decoy1_dist": round(cands[0][0], 3), "decoy2_dist": round(cands[1][0], 3),
        })
    out = pd.DataFrame(rows).sort_values(["lig_class", "acc"])
    out.to_csv(RES / "panel.tsv", sep="\t", index=False)
    (RES / "selection_rule.json").write_text(json.dumps({
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "rule": __doc__.split("SELECTION RULE")[1].split("DECOY RULE")[0].strip(),
        "decoy_rule": __doc__.split("DECOY RULE")[1].split("SITE")[0].strip(),
        "per_class": per_class, "n_positives": len(out),
        "class_counts": out.lig_class.value_counts().to_dict(),
        "biolip_source": "data/annotations/biolip_sites.parquet (see DATA.md)",
        "control_accs_excluded": sorted(CONTROL_ACCS),
        "ligand_property_source": "AF3 CCD SMILES + RDKit",
    }, indent=1))
    (RES / "ligand_properties.tsv").write_text(
        "ligand\t" + "\t".join(keys + ["n_heavy"]) + "\n" +
        "".join(f"{c}\t" + "\t".join(str(round(pr[c][k], 3)) for k in keys)
                + f"\t{pr[c]['n_heavy']}\n" for c in pool))
    log(f"wrote {RES}/panel.tsv ({len(out)} positives, {len(pool)} distinct ligands)")


# ---------------------------------------------------------------- build / run
def panel_rows() -> list[dict]:
    import pandas as pd
    return pd.read_csv(RES / "panel.tsv", sep="\t").to_dict("records")


def build() -> list[dict]:
    inp = BIG / "inputs"
    jobs = []
    for r in panel_rows():
        acc = r["acc"]
        seq = cached_seq(acc)
        site = [int(x) for x in str(r["site"]).split(",")]
        # The site positions are UniProt-numbered; the folded sequence is the cached AF3
        # query sequence. Refuse the protein if they are not consistent.
        if seq is None or not site or max(site) > len(seq):
            log(f"SKIP {acc}: cached seq {None if seq is None else len(seq)} aa vs site "
                f"max {max(site) if site else None}")
            continue
        gly = list(seq)
        for p in site:
            gly[p - 1] = "G"
        gly = "".join(gly)
        variants = [
            ("cognate", seq, r["ligand"], False),
            ("decoy1", seq, r["decoy1"], False),
            ("decoy2", seq, r["decoy2"], False),
            ("ablation", gly, r["ligand"], True),
        ]
        for tag, s, lig, var in variants:
            name = f"{acc.lower()}_{tag}"
            d = inp / acc
            d.mkdir(parents=True, exist_ok=True)
            out = d / f"{name}.json"
            cmd = [AF3_PY, str(MK), "--acc", acc, "--seq", s, "--no-uniprot",
                   "--no-templates", "--seeds", "1,2,3", "--name", name,
                   "--ligand", f"CCD:{lig}", "--out", str(out)]
            if var:
                cmd.append("--reuse-msa-for-variant")
            res = subprocess.run(cmd, text=True, capture_output=True)
            if res.returncode != 0:
                log(f"BUILD FAILED {name} rc={res.returncode}: {res.stderr.strip()[-300:]}")
                continue
            meta = out.with_suffix(".meta.json")
            if meta.exists():
                (d / "meta").mkdir(exist_ok=True)
                meta.rename(d / "meta" / meta.name)
            jobs.append({"acc": acc, "tag": tag, "name": name, "ligand": lig,
                         "lig_class": r["lig_class"], "seq_len": r["seq_len"],
                         "site": site, "input": str(out)})
    (BIG / "jobs.json").write_text(json.dumps(jobs, indent=1))
    log(f"built {len(jobs)} jobs over {len({j['acc'] for j in jobs})} proteins")
    return jobs


def run(gpus: list[int]) -> None:
    jobs = json.loads((BIG / "jobs.json").read_text())
    accs = sorted({j["acc"] for j in jobs})
    load = {a: max(j["seq_len"] for j in jobs if j["acc"] == a) for a in accs}
    queues: dict[int, list[str]] = {g: [] for g in gpus}
    for a in sorted(accs, key=lambda a: -load[a]):           # longest-first balance
        g = min(queues, key=lambda g: sum(load[x] for x in queues[g]))
        queues[g].append(a)
    procs = []
    for g, qa in queues.items():
        if not qa:
            continue
        script = " ; ".join(f'bash {RUN} -i {BIG / "inputs" / a} -j calib_{a} -g {g}'
                            for a in qa)
        logf = open(BIG / f"gpu{g}.out", "w")
        procs.append((g, subprocess.Popen(["bash", "-c", script], stdout=logf,
                                          stderr=subprocess.STDOUT,
                                          env={**os.environ, "AF3_OUT_ROOT": str(BIG / "af3")}),
                      logf, time.time(), qa))
        log(f"GPU {g}: {len(qa)} proteins ({sum(load[x] for x in qa)} aa)")
    for g, p, f, t0, qa in procs:
        rc = p.wait()
        f.close()
        log(f"GPU {g} done rc={rc} in {time.time()-t0:.0f}s ({len(qa)} proteins)")


# ---------------------------------------------------------------- scoring
def collect() -> list[dict]:
    """Per-job, per-seed readouts."""
    import statistics as st
    jobs = json.loads((BIG / "jobs.json").read_text())
    rows = []
    for j in jobs:
        jd = BIG / "af3" / f"calib_{j['acc']}" / j["name"]
        if not jd.exists():
            cand = sorted((BIG / "af3" / f"calib_{j['acc']}").glob(j["name"] + "*"))
            jd = cand[0] if cand else jd
        if not jd.exists():
            log(f"missing {jd}")
            continue
        samp, summ = summarise_job(jd, j["site"], pb_top=True)
        if summ is None:
            continue
        per_seed = {}
        for sd in sorted({s["seed"] for s in samp}):
            ss = [s for s in samp if s["seed"] == sd]
            per_seed[sd] = {k: (st.mean([s[k] for s in ss if s.get(k) is not None])
                                if any(s.get(k) is not None for s in ss) else None)
                            for k in ("lig_plddt", "prot_lig_pae_min",
                                      "site_centroid_dist", "site_recall")}
        rows.append({"acc": j["acc"], "tag": j["tag"], "ligand": j["ligand"],
                     "lig_class": j["lig_class"], "seq_len": j["seq_len"],
                     "label": 1 if j["tag"] == "cognate" else 0,
                     "neg_kind": "" if j["tag"] == "cognate"
                     else ("ablation" if j["tag"] == "ablation" else "decoy"),
                     "lig_plddt": summ["med_lig_plddt"],
                     "pae_min": summ["med_prot_lig_pae_min"],
                     "centroid": summ["med_site_centroid_dist"],
                     "site_recall": summ["med_site_recall"],
                     "iptm": summ["med_prot_lig_iptm"],
                     "pb_valid": summ.get("pb_valid"),
                     "per_seed": per_seed})
    (BIG / "readouts.json").write_text(json.dumps(rows, indent=1))
    return rows


def auc_ci(scores, labels, n_boot=4000, seed=0):
    """AUC (higher score = positive) with a stratified bootstrap percentile CI."""
    import numpy as np
    s, y = np.asarray(scores, float), np.asarray(labels, int)
    ok = ~np.isnan(s)
    s, y = s[ok], y[ok]

    def _auc(s, y):
        pos, neg = s[y == 1], s[y == 0]
        if len(pos) == 0 or len(neg) == 0:
            return float("nan")
        order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
        ranks = np.empty(len(order), float)
        srt = np.concatenate([pos, neg])[order]
        i = 0
        while i < len(srt):                       # average ranks over ties
            j = i
            while j + 1 < len(srt) and srt[j + 1] == srt[i]:
                j += 1
            ranks[order[i:j + 1]] = (i + j) / 2.0 + 1
            i = j + 1
        return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))

    a = _auc(s, y)
    rng = np.random.default_rng(seed)
    ip, ineg = np.where(y == 1)[0], np.where(y == 0)[0]
    boots = []
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(ip, len(ip), replace=True),
                              rng.choice(ineg, len(ineg), replace=True)])
        boots.append(_auc(s[idx], y[idx]))
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return round(float(a), 4), round(float(lo), 4), round(float(hi), 4)


def d23_pass(row: dict, t_plddt: float, t_cen: float, t_rec: float, min_seeds: int = 2) -> bool:
    """D23 composite: ligand pLDDT AND (centroid OR site recall) must clear their
    thresholds in at least `min_seeds` of the 3 seeds."""
    n = 0
    for v in row["per_seed"].values():
        if v["lig_plddt"] is None:
            continue
        geom = ((v["site_centroid_dist"] is not None and v["site_centroid_dist"] <= t_cen)
                or (v["site_recall"] is not None and v["site_recall"] >= t_rec))
        if v["lig_plddt"] >= t_plddt and geom:
            n += 1
    return n >= min_seeds


def score() -> None:
    import numpy as np
    RES.mkdir(parents=True, exist_ok=True)
    rows = collect()
    y = [r["label"] for r in rows]
    singles = {}
    for key, sign in (("lig_plddt", +1), ("pae_min", -1), ("centroid", -1),
                      ("site_recall", +1), ("iptm", +1)):
        sc = [sign * (r[key] if r[key] is not None else float("nan")) for r in rows]
        singles[key] = dict(zip(("auc", "ci_lo", "ci_hi"), auc_ci(sc, y)))
        singles[key]["direction"] = "higher_is_positive" if sign > 0 else "lower_is_positive"
    # per ligand class -- essential, because the readouts do not behave the same way for
    # a lone metal ion as for an organic ligand (see the caveats block below)
    by_class = {}
    for cls in sorted({r["lig_class"] for r in rows}):
        sel = [i for i, r in enumerate(rows) if r["lig_class"] == cls]
        by_class[cls] = {"n_pos": sum(rows[i]["label"] == 1 for i in sel),
                         "n_neg": sum(rows[i]["label"] == 0 for i in sel)}
        for key, sign in (("lig_plddt", +1), ("pae_min", -1), ("centroid", -1),
                          ("site_recall", +1)):
            sc = [sign * (rows[i][key] if rows[i][key] is not None else float("nan"))
                  for i in sel]
            by_class[cls][key] = dict(zip(("auc", "ci_lo", "ci_hi"),
                                          auc_ci(sc, [rows[i]["label"] for i in sel])))

    # same, restricted to each negative kind
    by_neg = {}
    for kind in ("decoy", "ablation"):
        sel = [i for i, r in enumerate(rows) if r["label"] == 1 or r["neg_kind"] == kind]
        by_neg[kind] = {}
        for key, sign in (("lig_plddt", +1), ("pae_min", -1), ("centroid", -1),
                          ("site_recall", +1)):
            sc = [sign * (rows[i][key] if rows[i][key] is not None else float("nan"))
                  for i in sel]
            by_neg[kind][key] = dict(zip(("auc", "ci_lo", "ci_hi"),
                                         auc_ci(sc, [rows[i]["label"] for i in sel])))

    # ---- D23 composite: grid search, maximise Youden's J subject to specificity >= 0.9
    plddts = np.round(np.arange(55, 96, 1.0), 1)
    cens = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    recs = [0.5, 0.6, 0.7, 0.8, 0.9]
    best, grid = None, []
    for tp, tc, tr in itertools.product(plddts, cens, recs):
        calls = [d23_pass(r, tp, tc, tr) for r in rows]
        tp_ = sum(c and r["label"] == 1 for c, r in zip(calls, rows))
        fp_ = sum(c and r["label"] == 0 for c, r in zip(calls, rows))
        npos = sum(y); nneg = len(y) - npos
        sens, spec = tp_ / npos, 1 - fp_ / nneg
        rec = {"t_plddt": float(tp), "t_centroid": tc, "t_recall": tr,
               "sensitivity": round(sens, 4), "specificity": round(spec, 4),
               "youden_j": round(sens + spec - 1, 4)}
        grid.append(rec)
        if spec >= 0.90 and (best is None or rec["youden_j"] > best["youden_j"]
                             or (rec["youden_j"] == best["youden_j"]
                                 and rec["sensitivity"] > best["sensitivity"])):
            best = rec
    feasible = best is not None
    if not feasible:                       # report the best J overall, flagged
        best = max(grid, key=lambda r: r["youden_j"])

    def grid_search(subset: list[dict]) -> tuple[dict, bool]:
        b, npos = None, sum(r["label"] for r in subset)
        nneg = len(subset) - npos
        if npos == 0 or nneg == 0:
            return {}, False
        for tp, tc, tr in itertools.product(plddts, cens, recs):
            calls = [d23_pass(r, tp, tc, tr) for r in subset]
            tpr = sum(c and r["label"] == 1 for c, r in zip(calls, subset)) / npos
            spc = 1 - sum(c and r["label"] == 0 for c, r in zip(calls, subset)) / nneg
            rec = {"t_plddt": float(tp), "t_centroid": tc, "t_recall": tr,
                   "sensitivity": round(tpr, 4), "specificity": round(spc, 4),
                   "youden_j": round(tpr + spc - 1, 4)}
            if spc >= 0.90 and (b is None or rec["youden_j"] > b["youden_j"]
                                or (rec["youden_j"] == b["youden_j"]
                                    and rec["sensitivity"] > b["sensitivity"])):
                b = rec
        return (b, True) if b else (max(
            (r for r in [None] if False), default={}), False)

    # Pre-specified stratum, NOT a tuning knob: ligand class is fixed by the selection
    # rule before any folding, and the metal subset is known to invert ligand pLDDT.
    nonmetal = [r for r in rows if r["lig_class"] != "metal"]
    nm_best, nm_feasible = grid_search(nonmetal)
    metal = [r for r in rows if r["lig_class"] == "metal"]
    m_best, m_feasible = grid_search(metal)

    # composite ROC: sweep t_plddt with the geometry thresholds frozen at `best`
    roc = []
    for tp in np.round(np.arange(0, 101, 1.0), 1):
        calls = [d23_pass(r, tp, best["t_centroid"], best["t_recall"]) for r in rows]
        npos = sum(y); nneg = len(y) - npos
        roc.append({"t_plddt": float(tp),
                    "tpr": sum(c and r["label"] == 1 for c, r in zip(calls, rows)) / npos,
                    "fpr": sum(c and r["label"] == 0 for c, r in zip(calls, rows)) / nneg})
    roc_sorted = sorted(roc + [{"t_plddt": None, "tpr": 0.0, "fpr": 0.0},
                               {"t_plddt": None, "tpr": 1.0, "fpr": 1.0}],
                        key=lambda r: (r["fpr"], r["tpr"]))
    comp_auc = float(np.trapezoid([r["tpr"] for r in roc_sorted],
                                  [r["fpr"] for r in roc_sorted]))
    comp_scores = [sum(1 for v in r["per_seed"].values()
                       if v["lig_plddt"] is not None
                       and ((v["site_centroid_dist"] is not None
                             and v["site_centroid_dist"] <= best["t_centroid"])
                            or (v["site_recall"] is not None
                                and v["site_recall"] >= best["t_recall"]))) * 1000
                   + (r["lig_plddt"] or 0) for r in rows]
    comp_rank_auc = auc_ci(comp_scores, y)

    runtimes = []
    for lg in (BIG / "af3").glob("*/run_af3.log"):
        txt = lg.read_text()
        import re
        for m in re.finditer(r"Running model inference with seed \d+ took ([\d.]+) seconds", txt):
            runtimes.append(float(m.group(1)))
    walls = []
    for lg in (BIG / "af3").glob("*/run_af3.log"):
        import re
        m = re.search(r"# rc=\d+ wall_s=([\d.]+)", lg.read_text())
        if m:
            walls.append(float(m.group(1)))

    out = {
        "frozen_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "panel": {"n_positives": sum(y), "n_negatives": len(y) - sum(y),
                  "n_decoy_negatives": sum(r["neg_kind"] == "decoy" for r in rows),
                  "n_ablation_negatives": sum(r["neg_kind"] == "ablation" for r in rows),
                  "class_counts": {c: sum(r["lig_class"] == c and r["label"] == 1
                                          for r in rows)
                                   for c in sorted({r["lig_class"] for r in rows})},
                  "selection_rule": "results/cofold/calibration/selection_rule.json",
                  "job_level_statistic": "median over 3 seeds x 5 samples"},
        "single_readout_auc": singles,
        "single_readout_auc_by_ligand_class": by_class,
        "single_readout_auc_by_negative_kind": by_neg,
        "d23_composite": {
            "definition": "ligand pLDDT >= t_plddt AND (site_centroid_dist <= t_centroid "
                          "OR site_recall >= t_recall), in >= 2 of 3 seeds",
            "thresholds": {k: best[k] for k in ("t_plddt", "t_centroid", "t_recall")},
            "sensitivity": best["sensitivity"], "specificity": best["specificity"],
            "youden_j": best["youden_j"],
            "auc_threshold_sweep": round(comp_auc, 4),
            "auc_seedcount_rank": dict(zip(("auc", "ci_lo", "ci_hi"), comp_rank_auc)),
            "specificity_target_0.9_feasible": feasible,
        },
        "d23_composite_by_stratum": {
            "note": "ligand class is fixed by the selection rule before folding, so these "
                    "strata are pre-specified, not post-hoc tuning. The PRIMARY frozen "
                    "threshold is the pooled one above; these are reported so a metal-site "
                    "case is not scored with a rule calibrated on organic ligands.",
            "non_metal": {"n_pos": sum(r["label"] for r in nonmetal),
                          "n_neg": len(nonmetal) - sum(r["label"] for r in nonmetal),
                          "thresholds": nm_best,
                          "specificity_target_0.9_feasible": nm_feasible},
            "metal": {"n_pos": sum(r["label"] for r in metal),
                      "n_neg": len(metal) - sum(r["label"] for r in metal),
                      "thresholds": m_best,
                      "specificity_target_0.9_feasible": m_feasible},
        },
        "runtime": {"n_jobs": len(walls), "median_job_wall_s": round(float(np.median(walls)), 1)
                    if walls else None,
                    "median_seed_inference_s": round(float(np.median(runtimes)), 1)
                    if runtimes else None,
                    "total_gpu_s": round(float(np.sum(walls)), 1) if walls else None},
        "ipTM_note": "reported for completeness only; D23 removes ipTM as a readout",
        "posebusters_note": "reported only, never part of the call (D23). NOT VALID for "
                            "single-metal ligands: PoseBusters assumes van der Waals "
                            "contact, so a correctly coordinated metal at ~2.4 A always "
                            "fails minimum_distance_to_protein and "
                            "volume_overlap_with_protein. Treat pb_valid as missing for "
                            "lig_class == metal.",
        "caveats": [
            "Ligand pLDDT is unreliable for single-metal ligands: on the metal subset the "
            "cognate ion scores LOWER than the site-ablated variant, because AF3 gives a "
            "lone ion poor per-atom confidence wherever it puts it. Read the per-class "
            "AUCs, not the pooled number, and do not call a metal-site case on pLDDT.",
            "Geometry readouts (centroid distance, site recall) separate the metal class "
            "cleanly where pLDDT and PAE do not -- further support for the D23 "
            "concordance rule requiring pLDDT AND a geometry readout.",
            "Panel decoys are hard in-class near-misses, so these AUCs are a lower bound "
            "relative to the caffeine-style decoys in docs/COFOLD.md Sec.4.3.",
        ],
    }
    (RES / "threshold.json").write_text(json.dumps(out, indent=1))
    keys = ["acc", "tag", "neg_kind", "ligand", "lig_class", "label", "lig_plddt",
            "pae_min", "centroid", "site_recall", "iptm", "pb_valid"]
    with open(RES / "panel_readouts.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in sorted(rows, key=lambda r: (r["acc"], r["tag"])):
            f.write("\t".join("" if r.get(k) is None else str(r[k]) for k in keys) + "\n")
    with open(RES / "d23_grid.tsv", "w") as f:
        gk = list(grid[0].keys())
        f.write("\t".join(gk) + "\n")
        for g in grid:
            f.write("\t".join(str(g[k]) for k in gk) + "\n")
    print(json.dumps(out, indent=1))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--select", action="store_true")
    ap.add_argument("--per-class", type=int, default=6)
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--gpus", default="0,1,2")
    a = ap.parse_args()
    BIG.mkdir(parents=True, exist_ok=True)
    if a.select:
        select(a.per_class)
    if a.build:
        build()
    if a.run:
        run([int(g) for g in a.gpus.split(",")])
    if a.score:
        score()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
