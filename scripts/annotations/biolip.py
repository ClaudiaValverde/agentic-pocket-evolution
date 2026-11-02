#!/usr/bin/env python
"""BioLiP2 biologically relevant ligand-receptor sites, restricted to atlas accessions, in UniProt numbering.

Raw (to /tmp/pe_big/raw/biolip2/; untrusted, parsed as TSV only):
  https://zhanggroup.org/BioLiP/download/BioLiP.txt.gz   21-column annotation file (readme.txt)
  https://zhanggroup.org/BioLiP/data/ligand.tsv.gz       CCD ligand summary (formula, name, InChIKey, ChEMBL id)

BioLiP gives binding residues in PDB numbering (col 8) and renumbered from 1 on the observed receptor sequence
(col 9 indexes col 21). No UniProt numbering is provided, so each receptor sequence is aligned (semi-global,
Biopython PairwiseAligner, end gaps free) to the current UniProt sequence of its UniProt id (col 18), and the
renumbered residues are mapped through the alignment. A site residue is kept only if its aligned UniProt residue
has the same amino acid; per-site mapping stats are stored.

Output: data/annotations/biolip_sites.parquet
  accession, uniprot_acc, pdb_id, chain, resolution, site_code, ligand (CCD id or peptide/dna/rna),
  ligand_class (regular/metal/peptide/dna/rna), ligand_name, ligand_inchikey, ligand_chain, ligand_serial,
  residues_pdb (str), residues_uniprot (list<int>), n_site_residues, n_mapped, aln_identity, aln_coverage,
  catalytic_uniprot (list<int>), ec, affinity_manual, affinity_moad, affinity_pdbbind, affinity_bindingdb, pubmed
plus biolip_meta.json.

Usage: .venv/bin/python -I scripts/annotations/biolip.py --accessions data/annotations/accessions.parquet
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import re
import subprocess

import warnings

import numpy as np

warnings.filterwarnings("ignore", message=".*renamed.*")
import pandas as pd

URL = "https://zhanggroup.org/BioLiP/download/BioLiP.txt.gz"
URL_LIG = "https://zhanggroup.org/BioLiP/data/ligand.tsv.gz"
RAW = "/tmp/pe_big/raw/biolip2"
COLS = ["pdb_id", "chain", "resolution", "site_code", "ligand", "ligand_chain", "ligand_serial",
        "res_pdb", "res_renum", "cat_pdb", "cat_renum", "ec", "go", "affinity_manual", "affinity_moad",
        "affinity_pdbbind", "affinity_bindingdb", "uniprot", "pubmed", "ligand_auth_seq", "receptor_seq"]
RES_RE = re.compile(r"([A-Z])(-?\d+)")

_ALIGNER = None


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio import Align
        al = Align.PairwiseAligner()
        al.mode = "global"
        al.match_score, al.mismatch_score = 2, -1
        al.open_gap_score, al.extend_gap_score = -2, -0.1   # receptor seq has internal gaps (unobserved residues)
        al.end_open_gap_score, al.end_extend_gap_score = 0.0, 0.0  # semi-global: end gaps free
        _ALIGNER = al
    return _ALIGNER


def map_positions(args):
    """Return (dict receptor_pos(1-based) -> uniprot_pos(1-based), identity, coverage)."""
    rseq, useq = args
    if not isinstance(useq, str) or not rseq:
        return None, np.nan, np.nan
    i = useq.find(rseq)
    if i >= 0:  # exact substring: fast path
        return {k + 1: i + k + 1 for k in range(len(rseq))}, 1.0, 1.0
    al = _aligner()
    try:
        aln = al.align(useq, rseq)[0]
    except Exception:  # noqa: BLE001
        return None, np.nan, np.nan
    m, ident, n = {}, 0, 0
    for (ts, te), (qs, qe) in zip(*aln.aligned):
        for k in range(te - ts):
            m[qs + k + 1] = ts + k + 1
            n += 1
            ident += useq[ts + k] == rseq[qs + k]
    return m, ident / max(n, 1), n / len(rseq)


def parse_res(s):
    if not isinstance(s, str) or not s.strip():
        return []
    return [(aa, int(num)) for aa, num in RES_RE.findall(s)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", default="data/annotations/accessions.parquet")
    ap.add_argument("--proteins", default="data/annotations/uniprot_proteins.parquet")
    ap.add_argument("--sequences", default="/tmp/pe_big/work/annotations/uniprot_sequences.parquet")
    ap.add_argument("--outdir", default="data/annotations")
    ap.add_argument("--workers", type=int, default=32)
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    raw, raw_lig = os.path.join(RAW, "BioLiP.txt.gz"), os.path.join(RAW, "ligand.tsv.gz")
    for u, p in ((URL, raw), (URL_LIG, raw_lig)):
        if not os.path.exists(p):
            subprocess.run(["curl", "-sSf", "-D", p + ".headers", "-o", p, u], check=True)
    last_mod = next((l.split(":", 1)[1].strip() for l in open(raw + ".headers") if l.lower().startswith("last-modified")),
                    None) if os.path.exists(raw + ".headers") else None

    B = pd.read_csv(raw, sep="\t", header=None, names=COLS, dtype=str, quoting=3, keep_default_na=False)
    n_raw = len(B)
    P = pd.read_parquet(a.proteins, columns=["accession", "uniprot_acc", "secondary_accs"]).dropna(subset=["uniprot_acc"])
    amap = {}
    for ac, pa, sec in P.itertuples(index=False):
        amap.setdefault(ac, (ac, pa))
        amap.setdefault(pa, (ac, pa))
        for s in (sec or "").split(";"):
            if s:
                amap.setdefault(s, (ac, pa))
    B["uniprot"] = B.uniprot.str.split(",").str[0].str.strip()
    B = B[B.uniprot.isin(amap.keys())].copy()
    B["accession"] = B.uniprot.map(lambda u: amap[u][0])
    B["uniprot_acc"] = B.uniprot.map(lambda u: amap[u][1])

    S = pd.read_parquet(a.sequences).set_index("uniprot_acc").sequence
    pairs = B[["receptor_seq", "uniprot_acc"]].drop_duplicates().reset_index(drop=True)
    with mp.Pool(a.workers) as pool:
        res = pool.map(map_positions, [(r, S.get(u)) for r, u in pairs.itertuples(index=False)], chunksize=200)
    pairs["map"], pairs["aln_identity"], pairs["aln_coverage"] = zip(*res)
    B = B.merge(pairs, on=["receptor_seq", "uniprot_acc"], how="left")

    rows_u, rows_c, nsite = [], [], []
    for r in B[["uniprot_acc", "map", "res_renum", "cat_renum"]].itertuples(index=False):
        useq = S.get(r.uniprot_acc)
        def conv(s):
            out = []
            for aa, num in parse_res(s):
                u = r.map.get(num) if isinstance(r.map, dict) else None
                if u is not None and useq is not None and u <= len(useq) and useq[u - 1] == aa:
                    out.append(u)
            return out
        rr = parse_res(r.res_renum)
        nsite.append(len(rr))
        rows_u.append(conv(r.res_renum))
        rows_c.append(sorted(set(conv(r.cat_renum.replace(";", " ")))) if r.cat_renum else [])
    B["residues_uniprot"], B["catalytic_uniprot"], B["n_site_residues"] = rows_u, rows_c, nsite
    B["n_mapped"] = B.residues_uniprot.map(len)

    # ligand.tsv: CCD, formula, InChI, InChIKey, SMILES, name, ChEMBL, DrugBank, ZINC; names may contain tabs
    import gzip
    lrows = []
    with gzip.open(raw_lig, "rt", encoding="utf-8", errors="replace") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                f += [""] * (9 - len(f))
            lrows.append((f[0], f[1], f[3], "\t".join(f[5:-3]), f[-3]))
    L = pd.DataFrame(lrows, columns=["ligand", "ligand_formula", "ligand_inchikey", "ligand_name", "ligand_chembl"])
    B = B.merge(L[["ligand", "ligand_name", "ligand_inchikey", "ligand_formula", "ligand_chembl"]],
                on="ligand", how="left")

    def lclass(lig, formula):
        if lig == "peptide" or lig == "III":
            return "peptide"
        if lig in ("dna", "rna"):
            return lig
        f = formula.split() if isinstance(formula, str) else []
        if len(f) == 1 and re.fullmatch(r"[A-Z][a-z]?\d*", f[0]) and not re.match(r"^(C|O|N|S|P|H|Cl|F|Br|I)\d*$", f[0]):
            return "metal"
        return "regular"
    B["ligand_class"] = [lclass(l, f) for l, f in zip(B.ligand, B.ligand_formula)]
    B["resolution"] = pd.to_numeric(B.resolution, errors="coerce").replace(-1.0, np.nan)
    out = B[["accession", "uniprot_acc", "pdb_id", "chain", "resolution", "site_code", "ligand", "ligand_class",
             "ligand_name", "ligand_inchikey", "ligand_chembl", "ligand_chain", "ligand_serial", "res_pdb",
             "residues_uniprot", "n_site_residues", "n_mapped", "aln_identity", "aln_coverage",
             "catalytic_uniprot", "ec", "affinity_manual", "affinity_moad", "affinity_pdbbind",
             "affinity_bindingdb", "pubmed"]].rename(columns={"res_pdb": "residues_pdb"})
    out = out.replace({"": None})
    out["biolip_file_last_modified"] = last_mod
    out.to_parquet(os.path.join(a.outdir, "biolip_sites.parquet"), index=False, compression="zstd")
    meta = {"source": URL, "last_modified": last_mod, "sha256": {"BioLiP.txt.gz": sha256(raw),
                                                                 "ligand.tsv.gz": sha256(raw_lig)},
            "n_rows_raw": int(n_raw), "n_rows_atlas": int(len(out)),
            "n_accessions": int(out.accession.nunique()), "n_pdb_chains": int(out[["pdb_id", "chain"]].drop_duplicates().shape[0]),
            "frac_site_residues_mapped": float(out.n_mapped.sum() / max(out.n_site_residues.sum(), 1)),
            "frac_sites_fully_mapped": float((out.n_mapped == out.n_site_residues).mean()),
            "rows_by_ligand_class": out.ligand_class.value_counts().to_dict(),
            "n_accessions_regular_ligand": int(out[out.ligand_class == "regular"].accession.nunique())}
    json.dump(meta, open(os.path.join(a.outdir, "biolip_meta.json"), "w"), indent=1)
    print(meta)


if __name__ == "__main__":
    main()
