"""Observation-error helpers: site matching by residue Jaccard, method-record hashing, Foldcomp entry access,
model sequence / pLDDT extraction, tune/test split key.

All residue sets are UniProt positions (ints) of fragment-1 models ("A:123" -> 123).
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path

import numpy as np

# Fields of the pocketcompare method record that define the fingerprint. backend/device/storage are excluded:
# ref/fast/GPU are bit-identical (docs/tools_pocket.md §4) and storage is the fp16 copy.
METHOD_FIELDS = ("method", "config", "config_sha256", "dims", "metric", "pocketcompare_version", "af_mode", "unit",
                 "min_points", "site", "locate")
ATLAS_METHOD = {"method": "pocketcompare-final:A4b:n1_linear|raw", "config": "final_method.yaml",
                "config_sha256": "b3ea8fe963ccc180d912ddde40a01547626d082fde8ea481232d5eda66a36a0c", "dims": 384,
                "metric": "inner_product", "pocketcompare_version": "1.0.0", "af_mode": True, "unit": "merged",
                "min_points": 10, "site": None, "locate": "shell"}

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
             "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
             "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O"}


# ------------------------------------------------------------------------------------------------ method record
def method_record_hash(rec: dict | str | bytes) -> str:
    """sha256 of the canonical JSON of the fingerprint-defining fields of a pocketcompare method record."""
    if isinstance(rec, (bytes, str)):
        rec = json.loads(rec)
    core = {k: rec.get(k) for k in METHOD_FIELDS}
    return hashlib.sha256(json.dumps(core, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


ATLAS_METHOD_HASH = method_record_hash(ATLAS_METHOD)


def parquet_method_hash(path: str | Path) -> str:
    import pyarrow.parquet as pq
    md = pq.read_schema(path).metadata or {}
    if b"pocketcompare" not in md:
        raise ValueError(f"{path}: no pocketcompare method record in schema metadata")
    return method_record_hash(md[b"pocketcompare"])


# ------------------------------------------------------------------------------------------------ residues
def parse_residues(s, offset: int = 0) -> frozenset:
    """'A:53 A:79' -> frozenset({53, 79}) (+ offset). Accepts list/array of ints too."""
    if s is None:
        return frozenset()
    if isinstance(s, str):
        out = []
        for tok in s.split():
            r = tok.split(":")[-1]
            r = "".join(ch for ch in r if ch.isdigit() or ch == "-")   # drop insertion codes
            if r:
                out.append(int(r) + offset)
        return frozenset(out)
    return frozenset(int(x) + offset for x in s)


def jaccard(a, b) -> float:
    a, b = set(a), set(b)
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def coverage(known, other) -> float:
    """Fraction of `known` residues contained in `other`."""
    known = set(known)
    return len(known & set(other)) / len(known) if known else 0.0


def best_match(ref: list, query: list, metric: str = "jaccard"):
    """For each ref residue set, (index of best query set, best score); index -1 if no query overlaps.
    metric: 'jaccard' or 'coverage' (fraction of ref residues covered by the query set)."""
    out = []
    fn = jaccard if metric == "jaccard" else coverage
    for r in ref:
        best, bi = 0.0, -1
        for j, q in enumerate(query):
            if not (set(r) & set(q)):
                continue
            v = fn(r, q)
            if v > best:
                best, bi = v, j
        out.append((bi, best))
    return out


def match_sites(ref: list, query: list, tau: float):
    """Re-detection: each ref site is matched to its best-Jaccard query site; matched if Jaccard >= tau.
    Returns list of (matched: bool, query index or -1, jaccard)."""
    return [(bi >= 0 and v >= tau, bi if v >= tau else -1, v) for bi, v in best_match(ref, query)]


# ------------------------------------------------------------------------------------------------ Foldcomp
_FD: dict = {}


def foldcomp_read(db_path: str, off: int, length: int) -> tuple[str, str]:
    """Decompress one Foldcomp entry by (offset, length) -> (entry name, PDB text)."""
    import foldcomp
    fd = _FD.get(db_path)
    if fd is None:
        fd = _FD[db_path] = os.open(db_path, os.O_RDONLY)
    data = os.pread(fd, int(length), int(off))
    if len(data) != int(length):
        raise IOError(f"short read at {off}")
    name, pdb = foldcomp.decompress(data)
    return str(name), pdb


# ------------------------------------------------------------------------------------------------ models
def read_model(path: str | Path):
    """-> (sequence str, resnums int array, per-residue pLDDT array) from a PDB/mmCIF(.gz) AF model (chain A /
    first chain, CA B-factor = pLDDT)."""
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    ch = st[0][0]
    seq, num, pl = [], [], []
    for res in ch:
        ca = res.find_atom("CA", "*")
        if ca is None:
            continue
        seq.append(THREE2ONE.get(res.name, "X"))
        num.append(res.seqid.num)
        pl.append(ca.b_iso)
    return "".join(seq), np.asarray(num, int), np.asarray(pl, float)


def seq_md5(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()


def region_plddt(resnums: np.ndarray, plddt: np.ndarray, residues) -> dict:
    """Mean pLDDT and fraction < 70 over a residue set (UniProt numbering = model numbering for F1)."""
    m = np.isin(resnums, np.fromiter(residues, int))
    if not m.any():
        return {"mean": np.nan, "frac_lt70": np.nan, "n_found": 0}
    v = plddt[m]
    return {"mean": float(v.mean()), "frac_lt70": float((v < 70).mean()), "n_found": int(m.sum())}


def plddt_band(x: float) -> str:
    if x is None or not np.isfinite(x):
        return "na"
    return "<70" if x < 70 else ("70-90" if x < 90 else ">=90")


# ------------------------------------------------------------------------------------------------ tune/test split
def split_half(key: str) -> str:
    """prereg tune_test.rule: int(sha1(key)) % 2 == 0 -> tune else test."""
    return "tune" if int(hashlib.sha1(key.encode()).hexdigest(), 16) % 2 == 0 else "test"


def root_hog_keys(orthoxml_gz: str | Path) -> dict:
    """accession -> lexicographically smallest accession of its root HOG (O1 orthoXML). Accessions not in any
    HOG are absent (callers use the accession itself as its key: a singleton root HOG)."""
    import xml.etree.ElementTree as ET
    ns = "{http://orthoXML.org/2011/}"
    gene2acc: dict = {}
    keys: dict = {}
    depth = 0
    cur: list = []
    with gzip.open(orthoxml_gz, "rb") as fh:
        for ev, el in ET.iterparse(fh, events=("start", "end")):
            tag = el.tag.replace(ns, "")
            if ev == "start":
                if tag == "orthologGroup" or tag == "paralogGroup":
                    depth += 1
                continue
            if tag == "gene":
                gene2acc[el.get("id")] = el.get("protId")
                el.clear()
            elif tag == "geneRef":
                cur.append(gene2acc.get(el.get("id")))
            elif tag in ("orthologGroup", "paralogGroup"):
                depth -= 1
                if depth == 0:
                    accs = [a for a in cur if a]
                    if accs:
                        k = min(accs)
                        for a in accs:
                            keys[a] = k
                    cur = []
                el.clear()
    return keys
