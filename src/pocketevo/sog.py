"""Site orthogroups (SOGs): DESIGN.md §3.2, config/prereg.yaml `sogs`.

Per root orthogroup (OG) with >= min_species species:
  model sequences -> FAMSA MSA -> site -> MSA column set -> inverted index (sparse incidence product) ->
  column-Jaccard graph between sites of *different* members (edge iff J >= tau) -> Leiden communities = SOGs ->
  footprint = columns present in >= footprint_member_fraction of the SOG's member sites (insertion columns allowed).

Everything here is pure (no I/O except the FAMSA subprocess) so it can be unit tested on toy MSAs.
No threshold is hard-coded: callers pass values read from config/prereg.yaml (see `prereg_params`).
"""
from __future__ import annotations

import gzip
import hashlib
import io
import os
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp

SCHEMA_VERSION = "sogs-v1"
GAP = ord("-")


# ----------------------------------------------------------------------------------------------- prereg
def prereg_sha256(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def prereg_params(prereg: dict) -> dict:
    """Flatten the prereg values used by the SOG/state code. Raises KeyError if any is missing (never default)."""
    s, st, tt = prereg["sogs"], prereg["states"], prereg["tune_test"]
    a = st["A_primary"]
    lo, hi = a["protein_length_vs_og_median"]["value"]
    return {
        "min_species": int(s["min_species_per_root_og"]["value"]),
        "tau_default": float(s["tau_column_jaccard"]["value"]),
        "tau_grid": [float(x) for x in s["tau_column_jaccard"]["grid"]],
        "leiden_resolution": float(s["leiden_resolution"]["value"]),
        "footprint_member_fraction": float(s["footprint_member_fraction"]["value"]),
        "footprint_allows_insertion_columns": bool(s["footprint_allows_insertion_columns"]["value"]),
        "low_identity": _parse_lt(s["low_identity_bin"]["value"]),
        "identity_bins": [tuple(b) for b in s["identity_bins_reported"]["value"]],
        "A_aligned_frac": float(a["footprint_aligned_non_gap_frac"]["value"]),
        "A_plddt": float(a["footprint_mean_plddt"]["value"]),
        "A_frac_lt70": float(a["footprint_frac_residues_below_70"]["value"]),
        "A_len_lo": float(lo), "A_len_hi": float(hi),
        "A_sptm_overlap_max": float(a["signal_peptide_propeptide_tm_overlap_max"]["value"]),
        "A_model_status": a["model_status"]["value"],
        "A_plddt_secondary": float(st["A_sensitivity_stratum"]["footprint_mean_plddt"]),
        "rescue_k": int(st["rescue"]["k_settings_to_rescue"]["value"]),
        "rescue_n_settings": int(st["rescue"]["settings"]["n_settings"]),
        "split_rule": tt["rule"], "split_key": tt["key"],
        "forced_test_accessions": list(tt["forced_test_half"]["root_hogs_containing"]),
        "seeds": list(prereg["meta"]["seeds"]),
        "detect_overlap_min": float(prereg["benchmarks"]["detectability_preflight"]["overlap_min"]),
    }


def _parse_lt(text: str) -> float:
    import re
    m = re.search(r"<\s*([0-9.]+)", text)
    if not m:
        raise KeyError(f"cannot parse low-identity bin from {text!r}")
    return float(m.group(1))


# ----------------------------------------------------------------------------------------------- tune/test
def split_key(accessions) -> str:
    """prereg tune_test.key: lexicographically smallest UniProt accession among root-HOG members."""
    return min(accessions)


def split_half(key: str) -> str:
    """prereg tune_test.rule: int(sha1(key).hexdigest(), 16) % 2 == 0 -> tune, else test."""
    return "tune" if int(hashlib.sha1(key.encode()).hexdigest(), 16) % 2 == 0 else "test"


# ----------------------------------------------------------------------------------------------- orthology
@dataclass
class Orthology:
    members: pd.DataFrame          # og_id, accession, species, taxid
    ogs: pd.DataFrame              # og_id, level, n_members, n_species
    tree_children: dict            # taxon name -> list of child names (species are leaves)
    species: list                  # species names (orthoxml <species name>)
    root: str = ""
    leaves_under: dict = field(default_factory=dict)

    def leaves(self, level: str) -> list:
        if level not in self.leaves_under:
            out, stack = [], [level]
            while stack:
                x = stack.pop()
                ch = self.tree_children.get(x, [])
                if ch:
                    stack.extend(ch)
                else:
                    out.append(x)
            self.leaves_under[level] = sorted(out)
        return self.leaves_under[level]

    def parent_map(self) -> dict:
        if not hasattr(self, "_parent"):
            self._parent = {c: p for p, ch in self.tree_children.items() for c in ch}
        return self._parent

    def lca(self, names) -> str:
        par = self.parent_map()
        def path(x):
            out = [x]
            while x in par:
                x = par[x]
                out.append(x)
            return out
        names = list(dict.fromkeys(names))
        common = path(names[0])
        for nm in names[1:]:
            ps = set(path(nm))
            common = [x for x in common if x in ps]
        return common[0] if common else self.root

    def newick(self, label=lambda x: x) -> str:
        def rec(n):
            ch = self.tree_children.get(n, [])
            if not ch:
                return label(n)
            return "(" + ",".join(rec(c) for c in ch) + ")" + _nwk_name(n)
        return rec(self.root) + ";"


def _nwk_name(n: str) -> str:
    return "'" + n.replace("'", "''") + "'" if any(c in n for c in " ,;:()[]'") else n


def species_taxid(name: str) -> int:
    """Ancestra species name `UP000005640_9606` -> 9606 (qfo leaf naming)."""
    return int(name.rsplit("_", 1)[1])


def read_orthoxml(path: str) -> Orthology:
    """Parse an Ancestra OrthoXML (0.5). Root groups = direct children of <groups>; og_id = their `id`."""
    ns = "{http://orthoXML.org/2011/}"
    opener = gzip.open if path.endswith(".gz") else open
    gene_prot, gene_species, species = {}, {}, []
    rows, ogs = [], []
    tree_children, root = {}, ""
    with opener(path, "rb") as fh:
        tree = ET.parse(fh)
    r = tree.getroot()
    for spc in r.findall(f"{ns}species"):
        nm = spc.get("name")
        species.append(nm)
        for g in spc.iter(f"{ns}gene"):
            gene_prot[g.get("id")] = g.get("protId")
            gene_species[g.get("id")] = nm
    tax = r.find(f"{ns}taxonomy")
    if tax is not None:
        def rec(t):
            nm = t.get("name")
            kids = t.findall(f"{ns}taxon")
            tree_children[nm] = [k.get("name") for k in kids]
            for k in kids:
                rec(k)
        top = tax.findall(f"{ns}taxon")
        root = top[0].get("name")
        for t in top:
            rec(t)
    groups = r.find(f"{ns}groups")
    for g in groups:
        og_id = g.get("id")
        level = None
        for p in g.findall(f"{ns}property"):
            if p.get("name") == "TaxRange":
                level = p.get("value")
        accs = [gene_prot[x.get("id")] for x in g.iter(f"{ns}geneRef")]
        sps = [gene_species[x.get("id")] for x in g.iter(f"{ns}geneRef")]
        for a, s in zip(accs, sps):
            rows.append((og_id, a, s))
        ogs.append((og_id, level, g.tag.replace(ns, ""), len(accs), len(set(sps))))
    members = pd.DataFrame(rows, columns=["og_id", "accession", "species"])
    members["taxid"] = members.species.map(species_taxid).astype("int32")
    og = pd.DataFrame(ogs, columns=["og_id", "level", "root_kind", "n_members", "n_species"])
    o = Orthology(members=members, ogs=og, tree_children=tree_children, species=species, root=root)
    miss = og.level.isna()
    if miss.any():  # paralogGroup roots carry no TaxRange: use the LCA of the member species
        sp_of = members[members.og_id.isin(set(og.og_id[miss]))].groupby("og_id").species.apply(list)
        og.loc[miss, "level"] = og.loc[miss, "og_id"].map(lambda x: o.lca(sp_of[x]))
    return o


# ----------------------------------------------------------------------------------------------- MSA
def famsa_version(famsa: str) -> str:
    p = subprocess.run([famsa], capture_output=True, text=True)
    out = p.stdout + p.stderr
    for line in out.splitlines():
        if "version" in line:
            return line.strip()
    return "unknown"


def read_fasta(text: str) -> dict:
    out, name, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def run_famsa(seqs: dict, famsa: str, threads: int = 1, tmpdir: str | None = None) -> dict:
    """Align {id: seq}. Returns {id: aligned}. 1 sequence -> identity; 2+ -> FAMSA (default parameters)."""
    if len(seqs) == 1:
        return dict(seqs)
    with tempfile.TemporaryDirectory(dir=tmpdir) as d:
        fi, fo = os.path.join(d, "in.fa"), os.path.join(d, "out.fa")
        with open(fi, "w") as fh:
            for k in sorted(seqs):
                fh.write(f">{k}\n{seqs[k]}\n")
        p = subprocess.run([famsa, "-t", str(threads), fi, fo], capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"famsa failed: {p.stderr[-2000:]}")
        with open(fo) as fh:
            aln = read_fasta(fh.read())
    if set(aln) != set(seqs):
        raise RuntimeError("famsa output ids differ from input")
    for k, v in aln.items():
        if v.replace("-", "").upper() != seqs[k].upper():
            raise RuntimeError(f"famsa changed sequence {k}")
    return aln


@dataclass
class MSA:
    ids: list            # member accessions, row order
    mat: np.ndarray      # (n, L) uint8 ascii, '-' = gap
    pos2col: list        # per row: int32 array, UniProt position p (1-based) -> column = pos2col[p-1]

    @property
    def n_cols(self) -> int:
        return self.mat.shape[1]

    @property
    def nongap(self) -> np.ndarray:
        return self.mat != GAP


def msa_from_aligned(aln: dict, order=None) -> MSA:
    ids = list(order) if order is not None else sorted(aln)
    L = len(aln[ids[0]]) if ids else 0
    mat = np.frombuffer("".join(aln[i] for i in ids).encode(), dtype=np.uint8).reshape(len(ids), L).copy()
    p2c = [np.flatnonzero(row != GAP).astype(np.int32) for row in mat]
    return MSA(ids=ids, mat=mat, pos2col=p2c)


def mean_pairwise_identity(msa: MSA, max_pairs: int, seed: int) -> float:
    """Identity = identical / columns where both are non-gap; mean over up to max_pairs random pairs (seeded)."""
    n = len(msa.ids)
    if n < 2:
        return float("nan")
    rng = np.random.default_rng(seed)
    all_pairs = n * (n - 1) // 2
    if all_pairs <= max_pairs:
        ii, jj = np.triu_indices(n, 1)
    else:
        ii = rng.integers(0, n, max_pairs * 2)
        jj = rng.integers(0, n, max_pairs * 2)
        keep = ii != jj
        ii, jj = ii[keep][:max_pairs], jj[keep][:max_pairs]
    ng = msa.nongap
    ids = []
    for i, j in zip(ii, jj):
        both = ng[i] & ng[j]
        nb = both.sum()
        if nb:
            ids.append((msa.mat[i][both] == msa.mat[j][both]).sum() / nb)
    return float(np.mean(ids)) if ids else float("nan")


# ----------------------------------------------------------------------------------------------- sites -> columns
def site_columns(msa: MSA, site_member_row: np.ndarray, site_residues: list) -> tuple[sp.csr_matrix, np.ndarray]:
    """Binary CSR (n_sites x n_cols) of MSA columns per site. Residues outside the member sequence are dropped
    and counted (returned per site)."""
    indptr, indices, dropped = [0], [], np.zeros(len(site_residues), dtype=np.int32)
    for k, (r, res) in enumerate(zip(site_member_row, site_residues)):
        p2c = msa.pos2col[r]
        res = np.asarray(res, dtype=np.int64)
        ok = (res >= 1) & (res <= len(p2c))
        dropped[k] = int((~ok).sum())
        cols = np.unique(p2c[res[ok] - 1])
        indices.append(cols)
        indptr.append(indptr[-1] + len(cols))
    idx = np.concatenate(indices) if indices else np.zeros(0, dtype=np.int32)
    m = sp.csr_matrix((np.ones(len(idx), dtype=np.int32), idx, np.array(indptr)),
                      shape=(len(site_residues), msa.n_cols))
    return m, dropped


def jaccard_pairs(S: sp.csr_matrix, site_member: np.ndarray, tau_min: float):
    """Inverted-index candidate pairs (sharing >= 1 column) between sites of different members, Jaccard >= tau_min.
    Returns (i, j, jaccard) with i < j."""
    if S.shape[0] < 2:
        return (np.zeros(0, np.int64),) * 2 + (np.zeros(0, np.float64),)
    S = S.astype(np.int32)
    inter = (S @ S.T).tocoo()
    m = inter.row < inter.col
    i, j, x = inter.row[m].astype(np.int64), inter.col[m].astype(np.int64), inter.data[m].astype(np.float64)
    diff = site_member[i] != site_member[j]
    i, j, x = i[diff], j[diff], x[diff]
    size = np.asarray(S.sum(axis=1)).ravel().astype(np.float64)
    jac = x / (size[i] + size[j] - x)
    keep = jac >= tau_min - 1e-12
    o = np.lexsort((j[keep], i[keep]))
    return i[keep][o], j[keep][o], jac[keep][o]


def leiden_labels(n: int, i: np.ndarray, j: np.ndarray, w: np.ndarray, resolution: float, seed: int) -> np.ndarray:
    """Leiden (RBConfiguration = modularity with resolution) on the weighted site graph. Deterministic given seed
    and node order. Labels are renumbered by (size desc, min node index asc). Isolated nodes = singleton SOGs."""
    import igraph as ig
    import leidenalg as la
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    if len(i) == 0:
        return np.arange(n, dtype=np.int64)
    g = ig.Graph(n=n, edges=list(zip(i.tolist(), j.tolist())), directed=False)
    g.es["weight"] = w.tolist()
    part = la.find_partition(g, la.RBConfigurationVertexPartition, weights="weight",
                             resolution_parameter=resolution, seed=seed, n_iterations=-1)
    lab = np.asarray(part.membership, dtype=np.int64)
    return canonical_labels(lab)


def canonical_labels(lab: np.ndarray) -> np.ndarray:
    """Renumber communities: larger first, ties by smallest node index (stable, input-order independent)."""
    u, inv, cnt = np.unique(lab, return_inverse=True, return_counts=True)
    first = np.full(len(u), np.iinfo(np.int64).max)
    np.minimum.at(first, inv, np.arange(len(lab)))
    order = np.lexsort((first, -cnt))
    rank = np.empty(len(u), dtype=np.int64)
    rank[order] = np.arange(len(u))
    return rank[inv]


def footprints(S: sp.csr_matrix, labels: np.ndarray, frac: float) -> sp.csr_matrix:
    """(n_sogs x n_cols) binary: columns present in >= frac of the SOG's member sites."""
    n_sogs = int(labels.max()) + 1 if len(labels) else 0
    M = sp.csr_matrix((np.ones(len(labels)), (labels, np.arange(len(labels)))), shape=(n_sogs, len(labels)))
    counts = (M @ S).tocsr()
    n_sites = np.asarray(M.sum(axis=1)).ravel()
    counts = counts.tocoo()
    keep = counts.data >= frac * n_sites[counts.row] - 1e-9
    F = sp.csr_matrix((np.ones(keep.sum(), dtype=np.int32), (counts.row[keep], counts.col[keep])),
                      shape=(n_sogs, S.shape[1]))
    F.sort_indices()
    return F


def build_sogs(S: sp.csr_matrix, site_member: np.ndarray, taus: list, resolution: float, seed: int,
               footprint_frac: float) -> dict:
    """For each tau: labels (per site) and footprint matrix. Edges computed once at min(taus)."""
    i, j, jac = jaccard_pairs(S, site_member, min(taus))
    out = {}
    for t in taus:
        m = jac >= t - 1e-12
        lab = leiden_labels(S.shape[0], i[m], j[m], jac[m], resolution, seed)
        out[t] = {"labels": lab, "footprint": footprints(S, lab, footprint_frac), "n_edges": int(m.sum())}
    out["_edges"] = (i, j, jac)
    return out


def footprint_jaccard(F: sp.csr_matrix, S: sp.csr_matrix) -> sp.coo_matrix:
    """Column Jaccard between every footprint (rows of F) and every site (rows of S) sharing >= 1 column."""
    inter = (F.astype(np.int32) @ S.astype(np.int32).T).tocoo()
    fs = np.asarray(F.sum(axis=1)).ravel().astype(np.float64)
    ss = np.asarray(S.sum(axis=1)).ravel().astype(np.float64)
    x = inter.data.astype(np.float64)
    jac = x / (fs[inter.row] + ss[inter.col] - x)
    return sp.coo_matrix((jac, (inter.row, inter.col)), shape=inter.shape)
