"""Homology classes for Part A site pairs (DESIGN §4.2; prereg pairs.homology_classes; round-2 amendments 1, 3).

Classes, first match wins, for an unordered pair of sites (i, j):
  H0   same protein (excluded from every estimand)
  H1   ortholog: same root group, both proteins are leaves of the orthology run and their LCA is a speciation node
  H2   paralog: same root group, LCA is a duplication node
  H12  same root group but at least one protein is placed by proxy (tier 2 / ID gap; no ortholog/paralog call)
  H3j  sequence-homologous by joins: the two sites share a region label (Pfam clan at the site region, TED CATH
       superfamily at the site region) or the proteins share a PANTHER family
  rest everything else -> resolved only by searches (MMseqs2 -> H3s, Foldseek -> H4, negative -> N, else unknown)

Root groups (a *partition* of proteins, used for H1/H2/H12, the bootstrap unit and the family-pair weighting):
  * protein in the orthology run            -> its root HOG
  * otherwise, PANTHER subfamily proxy      -> the root HOG holding most orthology-run proteins of the same PANTHER
                                               subfamily (flag 'proxy_hog'); if none, a pseudo-group per PANTHER
                                               subfamily ('proxy_sf'); placement tables (O2) override this
  * otherwise                               -> singleton group of the protein ('singleton')

Labels are integer codes; a site carries up to K_SITE region labels + the PANTHER family of its protein; a protein
carries up to K_PROT protein-level labels (all clans / CATH SFs / PANTHER family) used to exclude homologs from
null backgrounds. Overflow is counted, never silent.
"""
from __future__ import annotations

import gzip
import itertools
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass

import numpy as np
import pandas as pd

H0, H1, H2, H12, H3J, REST = 0, 1, 2, 3, 4, 5
CLASS_NAMES = {H0: "H0", H1: "H1", H2: "H2", H12: "H12", H3J: "H3j", REST: "rest"}
K_SITE = 6
K_PROT = 16


# ----------------------------------------------------------------------------------------------- orthoxml
@dataclass
class HogForest:
    gene_acc: np.ndarray      # (G,) accession of leaf g
    gene_parent: np.ndarray   # (G,) node index of the leaf's parent group (-1: not in any group)
    gene_pos: np.ndarray      # (G,) DFS leaf position
    node_parent: np.ndarray   # (M,) parent node (-1 for roots)
    node_type: np.ndarray     # (M,) 1 orthologGroup (speciation) / 2 paralogGroup (duplication)
    node_lo: np.ndarray       # (M,) min leaf position under node
    node_hi: np.ndarray       # (M,) max leaf position under node
    node_root: np.ndarray     # (M,) root node index
    root_names: list          # root node id strings, indexed by root order
    node_root_rank: np.ndarray  # (M,) root rank (0..R-1)


def parse_orthoxml(path: str) -> HogForest:
    """Stream-parse an OrthoXML file (Ancestra output) into flat arrays."""
    ns = "{http://orthoXML.org/2011/}"
    gid2acc = {}
    node_parent, node_type, node_root, node_lo, node_hi = [], [], [], [], []
    root_names, node_root_rank = [], []
    gene_acc, gene_parent, gene_pos = [], [], []
    stack: list[int] = []
    in_groups = False
    pos = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        for ev, el in ET.iterparse(fh, events=("start", "end")):
            tag = el.tag.replace(ns, "")
            if ev == "start":
                if tag == "groups":
                    in_groups = True
                elif in_groups and tag in ("orthologGroup", "paralogGroup"):
                    idx = len(node_parent)
                    par = stack[-1] if stack else -1
                    node_parent.append(par)
                    node_type.append(1 if tag == "orthologGroup" else 2)
                    if par < 0:
                        node_root.append(idx)
                        node_root_rank.append(len(root_names))
                        root_names.append(el.get("id") or f"root{len(root_names)}")
                    else:
                        node_root.append(node_root[par])
                        node_root_rank.append(node_root_rank[par])
                    node_lo.append(pos)
                    node_hi.append(-1)
                    stack.append(idx)
                elif in_groups and tag == "geneRef":
                    gene_acc.append(gid2acc[el.get("id")])
                    gene_parent.append(stack[-1] if stack else -1)
                    gene_pos.append(pos)
                    pos += 1
            else:
                if tag == "gene" and not in_groups:
                    gid2acc[el.get("id")] = el.get("protId")
                    el.clear()
                elif in_groups and tag in ("orthologGroup", "paralogGroup"):
                    idx = stack.pop()
                    node_hi[idx] = pos - 1
                    el.clear()
                elif tag == "species":
                    el.clear()
    return HogForest(np.array(gene_acc, dtype=object), np.array(gene_parent, np.int64), np.array(gene_pos, np.int64),
                     np.array(node_parent, np.int64), np.array(node_type, np.int8), np.array(node_lo, np.int64),
                     np.array(node_hi, np.int64), np.array(node_root, np.int64), root_names,
                     np.array(node_root_rank, np.int64))


def lca_type(f: HogForest, leaf_a: np.ndarray, leaf_b: np.ndarray) -> np.ndarray:
    """Vectorised LCA type for gene-index pairs in the same root: 1 speciation, 2 duplication, 0 if a == b."""
    node = f.gene_parent[leaf_a].copy()
    pb = f.gene_pos[leaf_b]
    out = np.zeros(len(leaf_a), np.int8)
    todo = leaf_a != leaf_b
    for _ in range(10_000):
        if not todo.any():
            break
        ok = todo & (f.node_lo[node] <= pb) & (pb <= f.node_hi[node])
        out[ok] = f.node_type[node[ok]]
        todo &= ~ok
        node[todo] = f.node_parent[node[todo]]
        if (node[todo] < 0).any():
            raise ValueError("pair not in the same root")
    return out


def exact_ortho_para_counts(f: HogForest, gene_weight: np.ndarray) -> pd.DataFrame:
    """Per root: Σ over unordered leaf pairs of w_a w_b split by LCA type (exact; one pass over the forest).

    At node v with children c_1..c_m (subgroups and direct leaves), pairs between different children have LCA v:
    Σ_{a<b} W_a W_b = ((Σ W)^2 − Σ W^2) / 2.
    """
    M = len(f.node_parent)
    W = np.zeros(M)            # subtree weight
    W2 = np.zeros(M)           # Σ over children of child weight^2
    has_leaf = f.gene_parent >= 0
    np.add.at(W, f.gene_parent[has_leaf], gene_weight[has_leaf])
    np.add.at(W2, f.gene_parent[has_leaf], gene_weight[has_leaf] ** 2)
    # children appear after parents in preorder -> process in reverse order to accumulate
    for v in range(M - 1, -1, -1):
        p = f.node_parent[v]
        if p >= 0:
            W[p] += W[v]
            W2[p] += W[v] ** 2
    pairs_at = (W ** 2 - W2) / 2
    df = pd.DataFrame({"root": f.node_root_rank, "type": f.node_type, "pairs": pairs_at})
    out = df.pivot_table(index="root", columns="type", values="pairs", aggfunc="sum", fill_value=0.0)
    out = out.rename(columns={1: "ortho", 2: "para"})
    for c in ("ortho", "para"):
        if c not in out:
            out[c] = 0.0
    return out[["ortho", "para"]]


# ----------------------------------------------------------------------------------------------- labels
def parse_chopping(chop: str) -> list[tuple[int, int]]:
    segs = []
    for s in str(chop).split("_"):
        if "-" in s:
            a, b = s.split("-")[:2]
            segs.append((int(a), int(b)))
    return segs


def region_labels_for_sites(site_acc: np.ndarray, site_res: list, regions: pd.DataFrame, min_overlap: float,
                            prefix: str) -> tuple[list[list[str]], int]:
    """For each site, labels of regions (accession, start, end, label) that cover >= min_overlap of its residues.

    `regions` may hold several rows per (accession, label) for multi-segment domains; overlap is computed per
    region id (column `rid`), so segments of one domain add up.
    """
    out: list[list[str]] = [[] for _ in range(len(site_acc))]
    if regions is None or len(regions) == 0:
        return out, 0
    by_acc = {a: g for a, g in regions.groupby("accession", sort=False)}
    n_hit = 0
    for k, (acc, res) in enumerate(zip(site_acc, site_res)):
        g = by_acc.get(acc)
        if g is None or len(res) == 0:
            continue
        res = np.asarray(res)
        cov: dict = {}
        for rid, lab, s, e in zip(g.rid.values, g.label.values, g.start.values, g.end.values):
            c = int(((res >= s) & (res <= e)).sum())
            if c:
                cov[(rid, lab)] = cov.get((rid, lab), 0) + c
        labs = sorted({lab for (rid, lab), c in cov.items() if c / len(res) >= min_overlap})
        if labs:
            out[k] = [prefix + x for x in labs]
            n_hit += 1
    return out, n_hit


def encode_label_slots(lists: list[list[str]], vocab: dict, K: int) -> tuple[np.ndarray, int]:
    """List-of-label-lists -> (n, K) int32 codes padded with -1; returns overflow count (labels dropped)."""
    arr = np.full((len(lists), K), -1, np.int32)
    over = 0
    for i, labs in enumerate(lists):
        codes = sorted({vocab.setdefault(x, len(vocab)) for x in labs})
        if len(codes) > K:
            over += len(codes) - K
            codes = codes[:K]
        arr[i, :len(codes)] = codes
    return arr, over


def share_any(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Row-wise: do padded label slot rows A[r] and B[r] share any code >= 0?"""
    out = np.zeros(len(A), bool)
    for a in range(A.shape[1]):
        col = A[:, a][:, None]
        out |= ((col == B) & (col >= 0)).any(axis=1)
    return out


# ----------------------------------------------------------------------------------------------- classification
@dataclass
class PairContext:
    """Everything needed to classify site pairs; arrays indexed by site (si) or protein (pi)."""
    site_prot: np.ndarray       # (N,) protein index
    site_lab: np.ndarray        # (N, K_SITE+1) region labels + PANTHER family (join labels)
    prot_group: np.ndarray      # (P,) root-group index (partition)
    prot_gene: np.ndarray       # (P,) gene index in the HogForest, -1 if proxy-placed / absent
    forest: HogForest | None


def classify_pairs(ctx: PairContext, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    pi, pj = ctx.site_prot[i], ctx.site_prot[j]
    out = np.full(len(i), REST, np.int8)
    same_p = pi == pj
    same_g = (ctx.prot_group[pi] == ctx.prot_group[pj]) & ~same_p
    out[same_p] = H0
    ga, gb = ctx.prot_gene[pi], ctx.prot_gene[pj]
    both = same_g & (ga >= 0) & (gb >= 0)
    out[same_g & ~both] = H12
    if both.any():
        t = lca_type(ctx.forest, ga[both], gb[both])
        o = np.where(t == 1, H1, np.where(t == 2, H2, H12)).astype(np.int8)
        out[both] = o
    rest = out == REST
    if rest.any():
        sh = share_any(ctx.site_lab[i[rest]], ctx.site_lab[j[rest]])
        tmp = out[rest]
        tmp[sh] = H3J
        out[rest] = tmp
    return out


# ----------------------------------------------------------------------------------------------- exact counts
def _site_weight_groups(site_mask, site_prot, prot_group, site_weight=None):
    """Sites (or summed site weights) per protein and per group among masked sites."""
    sp = site_prot[site_mask]
    sw = None if site_weight is None else np.asarray(site_weight, float)[site_mask]
    s_p = np.bincount(sp, weights=sw, minlength=len(prot_group)).astype(np.float64)
    s_g = np.bincount(prot_group, weights=s_p, minlength=prot_group.max() + 1)
    return s_p, s_g


def exact_class_totals(ctx: PairContext, site_mask: np.ndarray | None = None, group_weight: np.ndarray | None = None,
                       subsets: "LabelSubsets | None" = None, site_weight: np.ndarray | None = None) -> dict:
    """Exact pair totals for H1, H2, H12, H3j and the residual 'rest' (cross-group, no shared join label).

    Bootstrap weights: group g gets multiplicity w_g; a cross-group pair (g, h) gets w_g w_h, a within-group pair
    gets w_g (a cluster drawn twice contributes its internal pairs twice; pairs between the two copies of the same
    cluster are not real pairs). Both are unbiased (E w = 1, E w_g w_h = 1).
    """
    N = len(ctx.site_prot)
    if site_mask is None:
        site_mask = np.ones(N, bool)
    P = len(ctx.prot_group)
    G = int(ctx.prot_group.max()) + 1
    w = np.ones(G) if group_weight is None else np.asarray(group_weight, float)
    s_p, s_g = _site_weight_groups(site_mask, ctx.site_prot, ctx.prot_group, site_weight)
    # all different-protein pairs
    tot_w = (s_g * w).sum()
    cross = (tot_w ** 2 - ((s_g * w) ** 2).sum()) / 2
    within_diffprot_g = (s_g ** 2 - np.bincount(ctx.prot_group, weights=s_p ** 2, minlength=G)) / 2
    within = (w * within_diffprot_g).sum()
    # H1/H2 from the forest; H12 = within-group remainder
    if ctx.forest is not None:
        f = ctx.forest
        gw = np.zeros(len(f.gene_acc))
        has = ctx.prot_gene >= 0
        gw[ctx.prot_gene[has]] = s_p[has]
        op = exact_ortho_para_counts(f, gw)
        # map forest root rank -> group index of its proteins (groups for forest roots are root ranks)
        root_w = w[op.index.values]
        h1 = float((root_w * op["ortho"].values).sum())
        h2 = float((root_w * op["para"].values).sum())
    else:
        h1 = h2 = 0.0
    h12 = within - h1 - h2
    # H3j: cross-group pairs sharing >= 1 join label, inclusion-exclusion over label subsets
    if subsets is None:
        subsets = LabelSubsets.build(ctx, site_mask, site_weight=site_weight)
    h3 = subsets.cross_group_union(w)
    return {"total": cross + within, "within_group": within, "cross_group": cross, "H1": h1, "H2": h2, "H12": h12,
            "H3j": h3, "rest": cross - h3}


@dataclass
class LabelSubsets:
    """Sparse (label-subset S × group g) counts m_{S,g} of sites whose join-label set contains S."""
    sign: np.ndarray     # (nS,) +1/-1 for |S| odd/even
    S_idx: np.ndarray    # (nnz,) subset index
    g_idx: np.ndarray    # (nnz,) group index
    m: np.ndarray        # (nnz,) counts
    n_subsets: int

    @staticmethod
    def build(ctx: PairContext, site_mask: np.ndarray, max_labels: int = K_SITE + 1,
              site_weight: np.ndarray | None = None) -> "LabelSubsets":
        """Pairwise weights w_i w_j are supported through site_weight (counts become weight sums)."""
        lab = ctx.site_lab[site_mask]
        grp = ctx.prot_group[ctx.site_prot[site_mask]]
        swt = np.ones(int(site_mask.sum())) if site_weight is None else np.asarray(site_weight, float)[site_mask]
        nlab = (lab >= 0).sum(1)
        cnt: dict = {}
        keys = {}
        # group sites by (label tuple, group) first: far fewer distinct rows than sites
        df = pd.DataFrame(lab)
        df["g"] = grp
        df["_w"] = swt
        df = df[nlab > 0]
        agg = df.groupby([c for c in df.columns if c != "_w"], sort=False)["_w"].sum().reset_index(name="n")
        labcols = list(range(lab.shape[1]))
        for row in agg.itertuples(index=False):
            labs = tuple(int(x) for x in row[:len(labcols)] if x >= 0)
            g, n = int(row[len(labcols)]), float(row[len(labcols) + 1])
            for r in range(1, len(labs) + 1):
                for S in itertools.combinations(labs, r):
                    k = keys.setdefault(S, len(keys))
                    cnt[(k, g)] = cnt.get((k, g), 0.0) + n
        if cnt:
            arr = np.array([(k, g, n) for (k, g), n in cnt.items()], dtype=np.float64)
        else:
            arr = np.zeros((0, 3), np.float64)
        sizes = np.zeros(len(keys), np.int64)
        for S, k in keys.items():
            sizes[k] = len(S)
        sign = np.where(sizes % 2 == 1, 1.0, -1.0)
        return LabelSubsets(sign, arr[:, 0].astype(np.int64), arr[:, 1].astype(np.int64), arr[:, 2], len(keys))

    def cross_group_union(self, w: np.ndarray) -> float:
        if self.n_subsets == 0:
            return 0.0
        ww = w[self.g_idx]
        a = np.bincount(self.S_idx, weights=ww * self.m, minlength=self.n_subsets)
        b = np.bincount(self.S_idx, weights=(ww * self.m) ** 2, minlength=self.n_subsets)
        return float((self.sign * (a ** 2 - b) / 2).sum())


# ----------------------------------------------------------------------------------------------- pair searches
def run_pair_search(pairs: "pd.DataFrame", db: str, lookup: dict, workdir, tool: str, binary: str,
                    threads: int = 32, evalue: float = 10.0, extra=(), fields=None, log=print) -> "pd.DataFrame":
    """Align exactly the given (qname, tname) pairs with MMseqs2 `align` / Foldseek `structurealign`.

    A prefilter-type DB holding only these pairs is written with `tsv2db`, so no prefilter can miss a pair.
    `lookup` maps entry name -> DB key. Pairs whose name is not in the DB are returned with status 'missing'.
    Returns one row per input pair: qname, tname, status ('hit' | 'no_alignment' | 'missing') + requested fields.
    """
    import os
    import subprocess
    from pathlib import Path
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    qk = pairs.qname.map(lookup)
    tk = pairs.tname.map(lookup)
    ok = qk.notna() & tk.notna()
    sub = pd.DataFrame({"q": qk[ok].astype(np.int64), "t": tk[ok].astype(np.int64)}).drop_duplicates()
    sub = sub.sort_values(["q", "t"])
    tsv = wd / "pairs.tsv"
    with open(tsv, "w") as fh:
        for q, t in zip(sub.q.values, sub.t.values):
            fh.write(f"{q}\t{t}\t0\t0\n")
    for p in wd.glob("pref*"):
        p.unlink()
    for p in wd.glob("aln*"):
        p.unlink()
    fields = fields or (["query", "target", "evalue", "bits", "qstart", "qend", "tstart", "tend", "cigar"]
                        if tool == "mmseqs" else ["query", "target", "evalue", "bits", "alntmscore", "qtmscore",
                                                   "ttmscore", "prob"])
    run = lambda *c: subprocess.run([binary, *map(str, c)], check=True, capture_output=True, text=True)
    run("tsv2db", tsv, wd / "pref", "--output-dbtype", 7)
    if tool == "mmseqs":
        run("align", db, db, wd / "pref", wd / "aln", "-e", evalue, "-a", "--threads", threads, *extra)
    else:
        run("structurealign", db, db, wd / "pref", wd / "aln", "-e", evalue, "-a", "--threads", threads, *extra)
    out_tsv = wd / "aln.tsv"
    run("convertalis", db, db, wd / "aln", out_tsv, "--format-output", ",".join(fields), "--threads", threads)
    res = pd.read_csv(out_tsv, sep="\t", header=None, names=fields)
    res = res.sort_values("evalue").drop_duplicates(["query", "target"])
    base = pd.DataFrame({"qname": pairs.qname.values, "tname": pairs.tname.values})
    merged = base.merge(res.rename(columns={"query": "qname", "target": "tname"}), on=["qname", "tname"], how="left")
    merged["status"] = np.where(~ok.values, "missing", np.where(merged.evalue.notna(), "hit", "no_alignment"))
    log(f"{tool}: {len(pairs)} pairs, {int(ok.sum())} searchable, {int((merged.status == 'hit').sum())} aligned")
    return merged


def cigar_map(cigar: str, qstart: int, tstart: int) -> tuple[np.ndarray, np.ndarray]:
    """Aligned (match) position pairs, 1-based, from an MMseqs2 CIGAR (M/I/D; I = insertion in query)."""
    qs, ts = [], []
    q, t = qstart, tstart
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num) if num else 1
        num = ""
        if ch == "M":
            qs.append(np.arange(q, q + n)); ts.append(np.arange(t, t + n))
            q += n; t += n
        elif ch == "I":
            q += n
        elif ch == "D":
            t += n
    if not qs:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(qs), np.concatenate(ts)


def mapped_jaccard(res_q: np.ndarray, res_t: np.ndarray, qpos: np.ndarray, tpos: np.ndarray) -> float:
    """Jaccard of site q's residues mapped through the alignment onto t vs site t's residues (unmapped q residues
    stay in the union)."""
    if len(res_q) == 0 or len(res_t) == 0:
        return 0.0
    m = dict(zip(qpos.tolist(), tpos.tolist()))
    mapped = {m[r] for r in res_q.tolist() if r in m}
    n_unmapped = sum(1 for r in res_q.tolist() if r not in m)
    st = set(res_t.tolist())
    inter = len(mapped & st)
    union = len(mapped | st) + n_unmapped
    return inter / union if union else 0.0


# ----------------------------------------------------------------------------------------------- tune / test halves
def group_halves(prot_acc: np.ndarray, prot_group: np.ndarray, forest: "HogForest | None",
                 forced_test_accessions=()) -> np.ndarray:
    """Tune/test half per root group (prereg tune_test): key = lexicographically smallest UniProt accession among
    the root-HOG members (all forest genes of the root for HOG groups; the group's proteins otherwise);
    sha1(key) % 2 == 0 -> tune (0), else test (1); groups holding a forced-test accession -> test.
    Returns (G,) int8 array: 0 tune, 1 test."""
    import hashlib
    G = int(prot_group.max()) + 1
    key = np.full(G, None, dtype=object)
    if forest is not None:
        has = forest.gene_parent >= 0
        roots = forest.node_root_rank[forest.gene_parent[has]]
        df = pd.DataFrame({"g": roots, "acc": forest.gene_acc[has].astype(str)})
        mn = df.groupby("g").acc.min()
        key[mn.index.values] = mn.values
    pf = pd.DataFrame({"g": prot_group, "acc": prot_acc.astype(str)}).groupby("g").acc.min()
    miss = np.array([key[g] is None for g in pf.index.values])
    key[pf.index.values[miss]] = pf.values[miss]
    half = np.array([0 if k is None else (1 if int(hashlib.sha1(k.encode()).hexdigest(), 16) % 2 else 0)
                     for k in key], np.int8)
    forced = set(forced_test_accessions)
    if forced:
        fg = np.unique(prot_group[np.isin(prot_acc.astype(str), list(forced))])
        half[fg] = 1
        if forest is not None:
            gi = np.isin(forest.gene_acc.astype(str), list(forced)) & (forest.gene_parent >= 0)
            half[np.unique(forest.node_root_rank[forest.gene_parent[gi]])] = 1
    return half


def group_halves_from_table(n_groups: int, forest: "HogForest", table: pd.DataFrame,
                            col: str = "split_o2") -> np.ndarray:
    """Halves from a precomputed split table (og_id -> 'tune' | 'test' | 'crossover'; ledger D20).
    Root groups of the forest take the table's label; crossover, unlisted and non-HOG groups (proxy / singleton)
    get -1 and are excluded from half reporting. Returns (G,) int8: 0 tune, 1 test, -1 excluded."""
    lab = dict(zip(table.og_id.astype(str), table[col].astype(str)))
    half = np.full(n_groups, -1, np.int8)
    for r, name in enumerate(forest.root_names):
        if r >= n_groups:
            break
        v = lab.get(name)
        half[r] = 1 if v == "test" else (0 if v == "tune" else -1)
    return half
