"""Part B reconstruction engine (DESIGN.md §5.1): HOG gene trees, species tree, batched Sankoff parsimony.

Gene trees come from Ancestra OrthoXML (docs/tools_ancestra.md §2): every root group (orthologGroup or paralogGroup
at the top of <groups>) becomes one tree. Internal nodes: orthologGroup = speciation (kind SPEC, level = TaxRange),
paralogGroup = duplication (kind DUP, level = taxon of `taxonId`, i.e. the level just below the duplication).
Leaves = proteins (UniProt accessions) with their species. Nodes are multifurcating (NCBI polytomies, in-paralog
groups) and carry no branch lengths.

Sankoff: binary character "pocket present" (state 1) / "absent" (state 0). Cost matrix C[parent, child]:
C[0,1] = gain cost, C[1,0] = loss cost (prereg part_b.reconstruction: 2:1 primary; 1:1, 3:1 sensitivity).
Leaf costs from observation states (DESIGN §3.3):
  P -> present (absent forbidden); A -> absent (present forbidden); U -> missing (both free);
  W -> "present-ambiguous" (prereg W_handling): default mode "present" = W constrains the reconstruction like P
       (W is never absence), but W never counts as P evidence in candidate thresholds (events.py) and has state
       confidence 0.5. Modes "missing" and "soft:<x>" (absent costs x) exist for sensitivity.
  G (gene absent) is not a leaf (no member protein) and never enters the tree.
All characters that share one tree (all SOGs of a root HOG, or simulation replicates) are solved in one batch with
numba kernels: down pass, up (outside) pass, a single most-parsimonious reconstruction (MPR) with "no change"
tie-breaking, and per-node optimal costs conditional on each state, which identify nodes whose state is the same in
every MPR (unambiguous) and the cost margin of each event.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
import re
import subprocess
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numba import njit

SCHEMA_VERSION = "partb-recon-v1"
LEAF, SPEC, DUP = 0, 1, 2
KIND_NAME = {LEAF: "leaf", SPEC: "speciation", DUP: "duplication"}
BIG = 1.0e9
STATE_CODES = {"P": 0, "W": 1, "A": 2, "U": 3, "G": 4}
CODE_STATES = "PWAUG"

PROJECT = Path(__file__).resolve().parents[2]
PREREG = PROJECT / "config" / "prereg.yaml"


# --------------------------------------------------------------------------------------------- provenance
def file_sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def prereg_sha256() -> str:
    return file_sha256(PREREG)


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", str(PROJECT), "rev-parse", "HEAD"], capture_output=True, text=True,
                             timeout=10)
        commit = out.stdout.strip() or "unknown"
        dirty = subprocess.run(["git", "-C", str(PROJECT), "status", "--porcelain", "--", "src", "scripts"],
                               capture_output=True, text=True, timeout=30).stdout.strip()
        return commit + ("+dirty" if dirty else "")
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def code_sha256(*modules) -> str:
    h = hashlib.sha256()
    for m in modules or (__file__,):
        h.update(Path(m if isinstance(m, (str, Path)) else m.__file__).read_bytes())
    return h.hexdigest()


def write_sidecar(output: str | Path, snapshot_id: str, seed, extra: dict | None = None) -> Path:
    """Write `<output>.json` with the prereg-required fields {snapshot_id, prereg_sha256, git_commit, seed}."""
    output = Path(output)
    rec = {"output": output.name, "snapshot_id": snapshot_id, "prereg_sha256": prereg_sha256(),
           "git_commit": git_commit(), "seed": seed, "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if output.exists() and output.is_file():
        rec["sha256"] = file_sha256(output)
    rec.update(extra or {})
    side = output.with_name(output.name + ".json")
    side.write_text(json.dumps(rec, indent=1, default=str))
    return side


def load_prereg() -> dict:
    import yaml  # noqa: PLC0415  (pyyaml is in the venv; imported lazily)
    return yaml.safe_load(PREREG.read_text())


def costs_from_prereg(pr: dict | None = None) -> dict:
    pr = pr or load_prereg()
    r = pr["part_b"]["reconstruction"]
    grid = [tuple(x) for x in r["cost_sensitivity_grid"]]
    return {"primary": (float(r["gain_cost"]), float(r["loss_cost"])), "grid": grid}


# --------------------------------------------------------------------------------------------- trees
@dataclass
class Tree:
    """Rooted multifurcating tree in array form. Node 0 is the root; nodes are in preorder."""
    tree_id: str
    parent: np.ndarray            # int32, -1 for root
    kind: np.ndarray              # int8 LEAF/SPEC/DUP
    level: list                   # taxonomic level name per node (leaf: species name)
    node_name: list               # group id (internal) / accession (leaf)
    species: list                 # species name per node (leaf), "" internal
    taxid: np.ndarray             # int64 per node (leaf taxid, -1 internal)
    child_ptr: np.ndarray = field(init=False)
    child_idx: np.ndarray = field(init=False)
    post: np.ndarray = field(init=False)
    leaves: np.ndarray = field(init=False)          # leaf node indices in preorder (= left-to-right)
    leaf_lo: np.ndarray = field(init=False)         # per node: first leaf rank under node
    leaf_hi: np.ndarray = field(init=False)         # per node: one past last leaf rank under node
    leaf_rank: np.ndarray = field(init=False)       # per node: rank among leaves (-1 internal)

    def __post_init__(self):
        n = len(self.parent)
        self.parent = np.asarray(self.parent, dtype=np.int32)
        self.kind = np.asarray(self.kind, dtype=np.int8)
        self.taxid = np.asarray(self.taxid, dtype=np.int64)
        if n > 1 and np.any(self.parent[1:] >= np.arange(1, n)):
            raise ValueError("nodes must be in preorder (parent index < child index)")
        counts = np.bincount(self.parent[1:], minlength=n) if n > 1 else np.zeros(n, dtype=np.int64)
        self.child_ptr = np.zeros(n + 1, dtype=np.int32)
        self.child_ptr[1:] = np.cumsum(counts)
        order = np.argsort(self.parent[1:], kind="stable") + 1
        self.child_idx = order.astype(np.int32)
        self.post = _postorder(self.child_ptr, self.child_idx, n)
        is_leaf = counts == 0
        self.leaves = np.flatnonzero(is_leaf).astype(np.int32)
        self.leaf_rank = np.full(n, -1, dtype=np.int32)
        self.leaf_rank[self.leaves] = np.arange(len(self.leaves), dtype=np.int32)
        lo = np.full(n, np.iinfo(np.int32).max, dtype=np.int64)
        hi = np.full(n, -1, dtype=np.int64)
        lo[self.leaves] = self.leaf_rank[self.leaves]
        hi[self.leaves] = self.leaf_rank[self.leaves] + 1
        for v in self.post:
            p = self.parent[v]
            if p >= 0:
                lo[p] = min(lo[p], lo[v])
                hi[p] = max(hi[p], hi[v])
        self.leaf_lo, self.leaf_hi = lo.astype(np.int32), hi.astype(np.int32)
        if np.any(self.kind[self.leaves] != LEAF):
            bad = [self.node_name[i] for i in self.leaves if self.kind[i] != LEAF][:3]
            raise ValueError(f"{self.tree_id}: childless internal nodes {bad}")

    @property
    def n_nodes(self) -> int:
        return len(self.parent)

    @property
    def n_leaves(self) -> int:
        return len(self.leaves)

    def children(self, v: int) -> np.ndarray:
        return self.child_idx[self.child_ptr[v]:self.child_ptr[v + 1]]

    def leaf_names(self) -> list:
        return [self.node_name[i] for i in self.leaves]

    def leaf_species(self) -> list:
        return [self.species[i] for i in self.leaves]

    def leaf_taxids(self) -> np.ndarray:
        return self.taxid[self.leaves]

    def ancestors(self, v: int) -> list:
        out = []
        p = self.parent[v]
        while p >= 0:
            out.append(int(p))
            p = self.parent[p]
        return out

    def to_newick(self, v: int = 0) -> str:
        ch = self.children(v)
        lab = self.node_name[v] if self.kind[v] == LEAF else f"{KIND_NAME[self.kind[v]][0].upper()}_{self.level[v]}"
        lab = re.sub(r"[ (),:;']", "_", str(lab))
        if len(ch) == 0:
            return lab
        return "(" + ",".join(self.to_newick(int(c)) for c in ch) + ")" + lab


@njit(cache=True)
def _postorder(child_ptr, child_idx, n):
    out = np.empty(n, dtype=np.int32)
    stack = np.empty(n, dtype=np.int32)
    seen = np.zeros(n, dtype=np.bool_)
    sp = 0
    k = 0
    if n == 0:
        return out
    stack[0] = 0
    sp = 1
    while sp > 0:
        v = stack[sp - 1]
        if not seen[v]:
            seen[v] = True
            for j in range(child_ptr[v + 1] - 1, child_ptr[v] - 1, -1):
                stack[sp] = child_idx[j]
                sp += 1
        else:
            sp -= 1
            out[k] = v
            k += 1
    return out


def tree_from_nested(tree_id: str, nested, species_of: dict | None = None, taxid_of: dict | None = None) -> Tree:
    """Build a Tree from nested tuples (kind, level, name, [children]) / leaf strings (accession).

    Leaf = str (accession) or ("leaf", species, accession). Internal = ("S"|"D", level, name, [children]).
    `species_of` maps accession -> species name; `taxid_of` maps species name -> taxid.
    """
    parent, kind, level, name, species, taxid = [], [], [], [], [], []

    def rec(x, p):
        i = len(parent)
        if isinstance(x, str) or (isinstance(x, tuple) and x[0] == "leaf"):
            acc = x if isinstance(x, str) else x[2]
            sp = (species_of or {}).get(acc, x[1] if isinstance(x, tuple) else "")
            parent.append(p); kind.append(LEAF); level.append(sp); name.append(acc); species.append(sp)
            taxid.append(int((taxid_of or {}).get(sp, -1)))
            return
        k, lev, nm, ch = x
        parent.append(p); kind.append(SPEC if k in ("S", "speciation") else DUP); level.append(lev); name.append(nm)
        species.append(""); taxid.append(-1)
        for c in ch:
            rec(c, i)

    rec(nested, -1)
    return Tree(tree_id, np.array(parent), np.array(kind), level, name, species, np.array(taxid))


def taxid_from_species_name(name: str) -> int:
    """Ancestra species names are `UP000005640_9606` (qfo leaf naming) -> 9606; plain integers pass through."""
    m = re.search(r"(\d+)$", str(name))
    return int(m.group(1)) if m else -1


def iter_orthoxml_trees(path: str | Path, min_species: int = 1, keep=None):
    """Stream root groups of an Ancestra OrthoXML file as Tree objects (one per root HOG).

    `keep` (optional set of root group ids) restricts parsing output. Genes not in any group (singletons) are
    skipped. Yields Tree with tree_id = root group id (e.g. HOG0000001).
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    ns = "{http://orthoXML.org/2011/}"
    gene_acc, gene_sp, taxon_name = {}, {}, {}
    cur_species = None
    depth = 0
    in_groups = False
    with opener(path, "rb") as fh:
        ctx = ET.iterparse(fh, events=("start", "end"))
        for ev, el in ctx:
            tag = el.tag.replace(ns, "")
            if ev == "start":
                depth += 1
                if tag == "species":
                    cur_species = el.get("name")
                elif tag == "groups":
                    in_groups = True
                    groups_depth = depth
                continue
            # end
            if tag == "gene":
                gene_acc[el.get("id")] = el.get("protId")
                gene_sp[el.get("id")] = cur_species
                el.clear()
            elif tag == "taxon":
                taxon_name[el.get("id")] = el.get("name")
            elif in_groups and tag in ("orthologGroup", "paralogGroup") and depth == groups_depth + 1:
                gid = el.get("id")
                if keep is None or gid in keep:
                    t = _tree_from_element(gid, el, ns, gene_acc, gene_sp, taxon_name)
                    if len(set(t.species[i] for i in t.leaves)) >= min_species:
                        yield t
                el.clear()
            elif tag == "species":
                el.clear()
            depth -= 1


def _tree_from_element(gid, el, ns, gene_acc, gene_sp, taxon_name) -> Tree:
    parent, kind, level, name, species, taxid = [], [], [], [], [], []
    stack = [(el, -1)]
    while stack:
        e, p = stack.pop()
        tag = e.tag.replace(ns, "")
        i = len(parent)
        if tag == "geneRef":
            g = e.get("id")
            sp = gene_sp.get(g, "")
            parent.append(p); kind.append(LEAF); level.append(sp); name.append(gene_acc.get(g, g))
            species.append(sp); taxid.append(taxid_from_species_name(sp))
            continue
        if tag == "orthologGroup":
            lev = None
            for c in e:
                if c.tag.replace(ns, "") == "property" and c.get("name") == "TaxRange":
                    lev = c.get("value")
            lev = lev or taxon_name.get(e.get("taxonId"), "")
            k = SPEC
        else:
            lev = taxon_name.get(e.get("taxonId"), "")
            k = DUP
        parent.append(p); kind.append(k); level.append(lev); name.append(e.get("id") or f"{gid}.n{i}")
        species.append(""); taxid.append(-1)
        kids = [c for c in e if c.tag.replace(ns, "") in ("orthologGroup", "paralogGroup", "geneRef")]
        for c in reversed(kids):     # stack -> preorder in document order
            stack.append((c, i))
    # stack-based DFS already produced preorder with parent < child
    return Tree(gid, np.array(parent), np.array(kind), level, name, species, np.array(taxid))


# --------------------------------------------------------------------------------------------- species tree
def parse_newick(s: str):
    """Minimal Newick parser -> (parent list, names list) in preorder. Handles quoted labels, lengths ignored."""
    s = s.strip().rstrip(";")
    pos = 0
    parent, names = [], []

    def label():
        nonlocal pos
        if pos < len(s) and s[pos] == "'":
            j = s.index("'", pos + 1)
            lab = s[pos + 1:j]
            pos = j + 1
        else:
            m = re.match(r"[^,():;]*", s[pos:])
            lab = m.group(0)
            pos += len(lab)
        if pos < len(s) and s[pos] == ":":
            m = re.match(r":[^,();]*", s[pos:])
            pos += len(m.group(0))
        return lab.strip()

    def node(p):
        nonlocal pos
        i = len(parent)
        parent.append(p); names.append("")
        if s[pos] == "(":
            pos += 1
            while True:
                node(i)
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"bad newick at {pos}")
        names[i] = label()

    node(-1)
    return parent, names


class SpeciesTree:
    """NCBI guide tree (polytomies kept). Leaves keyed by taxid; internal nodes by level name."""

    def __init__(self, newick: str):
        parent, names = parse_newick(newick)
        n = len(parent)
        kind = np.array([SPEC] * n)
        has_child = np.zeros(n, bool)
        for p in parent[1:]:
            has_child[p] = True
        kind[~has_child] = LEAF
        taxid = np.array([taxid_from_species_name(nm) if not has_child[i] else -1 for i, nm in enumerate(names)])
        self.tree = Tree("species_tree", np.array(parent), kind, names, names,
                         [nm if not has_child[i] else "" for i, nm in enumerate(names)], taxid)
        t = self.tree
        self.names = names
        self.leaf_of_taxid = {int(t.taxid[v]): int(v) for v in t.leaves}
        self.leaf_of_species = {names[v]: int(v) for v in t.leaves}
        self.node_of_level = {nm: i for i, nm in enumerate(names) if has_child[i]}
        self.depth = np.zeros(n, dtype=np.int32)
        for v in range(1, n):
            self.depth[v] = self.depth[parent[v]] + 1
        self.n_species_below = (t.leaf_hi - t.leaf_lo).astype(np.int32)

    @classmethod
    def from_file(cls, path):
        return cls(Path(path).read_text())

    def lca(self, nodes) -> int:
        nodes = list(nodes)
        if not nodes:
            return -1
        t = self.tree
        lo = min(int(t.leaf_lo[v]) for v in nodes)
        hi = max(int(t.leaf_hi[v]) for v in nodes)
        v = nodes[0]
        while not (t.leaf_lo[v] <= lo and t.leaf_hi[v] >= hi):
            v = t.parent[v]
        return int(v)

    def lca_of_taxids(self, taxids) -> int:
        return self.lca([self.leaf_of_taxid[int(x)] for x in taxids if int(x) in self.leaf_of_taxid])

    def level_depth(self, level: str) -> int:
        if level in self.node_of_level:
            return int(self.depth[self.node_of_level[level]])
        if level in self.leaf_of_species:
            return int(self.depth[self.leaf_of_species[level]])
        return -1


# --------------------------------------------------------------------------------------------- Sankoff
def cost_matrix(gain: float, loss: float) -> np.ndarray:
    return np.array([[0.0, gain], [loss, 0.0]])


def leaf_costs(codes: np.ndarray, w_mode: str = "present") -> np.ndarray:
    """codes: int array (n_leaves, R) of STATE_CODES -> cost array (n_leaves, R, 2) [absent, present]."""
    codes = np.asarray(codes)
    lc = np.zeros(codes.shape + (2,), dtype=np.float64)
    P, W, A = STATE_CODES["P"], STATE_CODES["W"], STATE_CODES["A"]
    lc[codes == P, 0] = BIG
    lc[codes == A, 1] = BIG
    if w_mode == "present":
        lc[codes == W, 0] = BIG
    elif w_mode.startswith("soft:"):
        lc[codes == W, 0] = float(w_mode.split(":", 1)[1])
    elif w_mode != "missing":
        raise ValueError(f"unknown w_mode {w_mode}")
    if np.any(codes == STATE_CODES["G"]):
        raise ValueError("state G is not a leaf state (gene absent proteins are not tree leaves)")
    return lc


@njit(cache=True)
def _sankoff_kernel(post, parent, child_ptr, child_idx, leaf_rank, lc, C):
    n = post.shape[0]
    R = lc.shape[1]
    down = np.zeros((n, R, 2))
    for k in range(n):
        v = post[k]
        if child_ptr[v + 1] == child_ptr[v]:
            lr = leaf_rank[v]
            for r in range(R):
                down[v, r, 0] = lc[lr, r, 0]
                down[v, r, 1] = lc[lr, r, 1]
        else:
            for j in range(child_ptr[v], child_ptr[v + 1]):
                c = child_idx[j]
                for r in range(R):
                    for s in range(2):
                        a = down[c, r, 0] + C[s, 0]
                        b = down[c, r, 1] + C[s, 1]
                        down[v, r, s] += a if a < b else b
    up = np.zeros((n, R, 2))
    for k in range(n - 1, -1, -1):
        v = post[k]
        for j in range(child_ptr[v], child_ptr[v + 1]):
            c = child_idx[j]
            for r in range(R):
                for t in range(2):
                    best = 1e300
                    for s in range(2):
                        a = down[c, r, 0] + C[s, 0]
                        b = down[c, r, 1] + C[s, 1]
                        m = a if a < b else b
                        val = up[v, r, s] + down[v, r, s] - m + C[s, t]
                        if val < best:
                            best = val
                    up[c, r, t] = best
    # single MPR, preorder (node index order is preorder), ties -> no change; root tie -> present
    mpr = np.zeros((n, R), dtype=np.int8)
    for r in range(R):
        mpr[0, r] = 1 if down[0, r, 1] <= down[0, r, 0] else 0
    for v in range(1, n):
        p = parent[v]
        for r in range(R):
            s = mpr[p, r]
            c_same = down[v, r, s] + C[s, s]
            c_diff = down[v, r, 1 - s] + C[s, 1 - s]
            mpr[v, r] = s if c_same <= c_diff else 1 - s
    return down, up, mpr


@dataclass
class SankoffResult:
    tree: Tree
    down: np.ndarray        # (n_nodes, R, 2)
    up: np.ndarray
    mpr: np.ndarray         # (n_nodes, R) int8
    opt: np.ndarray         # (R,)
    total: np.ndarray       # (n_nodes, R, 2) optimal cost given node state
    C: np.ndarray

    @property
    def unambiguous(self) -> np.ndarray:
        """(n_nodes, R) bool: node state identical in every MPR."""
        tol = 1e-9
        alt = np.take_along_axis(self.total, (1 - self.mpr)[..., None].astype(np.int64), axis=2)[..., 0]
        return alt > self.opt[None, :] + tol

    def margin(self) -> np.ndarray:
        """(n_nodes, R) extra cost of the best reconstruction that flips this node's MPR state."""
        alt = np.take_along_axis(self.total, (1 - self.mpr)[..., None].astype(np.int64), axis=2)[..., 0]
        return np.minimum(alt - self.opt[None, :], BIG)

    def p_state_ambiguous(self) -> np.ndarray:
        return ~self.unambiguous


def sankoff(tree: Tree, codes: np.ndarray, gain: float = 2.0, loss: float = 1.0, w_mode: str = "present"):
    """Batched Sankoff. codes (n_leaves, R) STATE_CODES in leaf order (tree.leaves)."""
    codes = np.asarray(codes)
    if codes.ndim == 1:
        codes = codes[:, None]
    if codes.shape[0] != tree.n_leaves:
        raise ValueError(f"codes rows {codes.shape[0]} != n_leaves {tree.n_leaves}")
    lc = leaf_costs(codes, w_mode)
    C = cost_matrix(gain, loss)
    down, up, mpr = _sankoff_kernel(tree.post, tree.parent, tree.child_ptr, tree.child_idx, tree.leaf_rank, lc, C)
    total = down + up
    opt = np.minimum(down[0, :, 0], down[0, :, 1])
    return SankoffResult(tree, down, up, mpr, opt, total, C)


def encode_states(states) -> np.ndarray:
    """Iterable of 'P'/'W'/'A'/'U' strings -> int8 codes; anything unknown (incl. None) -> U."""
    m = STATE_CODES
    return np.array([m.get(s, m["U"]) if s != "G" else m["U"] for s in states], dtype=np.int8)


# --------------------------------------------------------------------------------------------- Dollo (M3)
@njit(cache=True)
def _dollo_kernel(post, parent, child_ptr, leaf_rank, X):
    """X (n_chars, n_leaves) int8: 1 present, 0 absent, <0 missing. Returns (n_events, lca) per char.

    Dollo: one gain at the LCA of present leaves (not counted when the LCA is the root, i.e. ancestral), then the
    minimal number of losses = number of maximal P-free subtrees below the LCA containing >= 1 absent leaf.
    Characters without any present leaf get n_events = -1 (undefined, excluded from M3).
    """
    m, nl = X.shape
    n = post.shape[0]
    ev = np.zeros(m, dtype=np.int32)
    lcas = np.full(m, -1, dtype=np.int32)
    cP = np.zeros(n, dtype=np.int32)
    cA = np.zeros(n, dtype=np.int32)
    for i in range(m):
        cP[:] = 0
        cA[:] = 0
        totP = 0
        for k in range(n):
            v = post[k]
            if child_ptr[v + 1] == child_ptr[v]:
                x = X[i, leaf_rank[v]]
                if x == 1:
                    cP[v] = 1
                    totP += 1
                elif x == 0:
                    cA[v] = 1
            p = parent[v]
            if p >= 0:
                cP[p] += cP[v]
                cA[p] += cA[v]
        if totP == 0:
            ev[i] = -1
            continue
        lca = -1
        for k in range(n):
            v = post[k]
            if cP[v] == totP:
                lca = v
                break
        lcas[i] = lca
        # nodes in subtree(lca): preorder indices lca..; check membership via ancestry walk
        losses = 0
        for v in range(lca + 1, n):
            p = parent[v]
            if cP[v] == 0 and cA[v] > 0 and cP[p] > 0:
                # v must be inside subtree(lca)
                u = p
                inside = False
                while u >= 0:
                    if u == lca:
                        inside = True
                        break
                    u = parent[u]
                if inside:
                    losses += 1
        ev[i] = losses + (1 if lca != 0 else 0)
    return ev, lcas


def dollo_events(tree: Tree, X: np.ndarray):
    X = np.ascontiguousarray(np.asarray(X, dtype=np.int8))
    if X.ndim == 1:
        X = X[None, :]
    return _dollo_kernel(tree.post, tree.parent, tree.child_ptr, tree.leaf_rank, X)


# --------------------------------------------------------------------------------------------- tree tables
def trees_to_table(trees):
    """Flatten Trees into column lists (for parquet): one row per node."""
    cols = {k: [] for k in ["og_id", "node", "parent", "kind", "level", "name", "species", "taxid"]}
    for t in trees:
        n = t.n_nodes
        cols["og_id"].extend([t.tree_id] * n)
        cols["node"].extend(range(n))
        cols["parent"].extend(t.parent.tolist())
        cols["kind"].extend(t.kind.tolist())
        cols["level"].extend(t.level)
        cols["name"].extend(t.node_name)
        cols["species"].extend(t.species)
        cols["taxid"].extend(t.taxid.tolist())
    return cols


def trees_from_table(table) -> dict:
    """Inverse of trees_to_table (pyarrow Table or pandas DataFrame) -> {og_id: Tree}."""
    import pandas as pd  # noqa: PLC0415
    df = table.to_pandas() if hasattr(table, "to_pandas") else pd.DataFrame(table)
    out = {}
    for og, g in df.groupby("og_id", sort=False):
        g = g.sort_values("node")
        out[og] = Tree(og, g["parent"].to_numpy(), g["kind"].to_numpy(), g["level"].tolist(), g["name"].tolist(),
                       g["species"].tolist(), g["taxid"].to_numpy())
    return out


def n_workers(default: int = 16) -> int:
    return max(1, min(default, int(os.environ.get("PE_WORKERS", default))))


# --------------------------------------------------------------------------------------------- CLI helpers
def default_snapshot() -> str:
    """Snapshot id: $PE_SNAPSHOT, else config/prereg.yaml meta.snapshot_id (no hard-coded data ids in scripts)."""
    return os.environ.get("PE_SNAPSHOT") or str(load_prereg()["meta"]["snapshot_id"])


def cli_value(flag: str, default=None, argv=None):
    """Peek a `--flag value` / `--flag=value` from argv without consuming it (module-level constants of drivers)."""
    argv = sys.argv[1:] if argv is None else argv
    for i, x in enumerate(argv):
        if x == flag and i + 1 < len(argv):
            return argv[i + 1]
        if x.startswith(flag + "="):
            return x.split("=", 1)[1]
    return default


def add_partb_args(ap, orthology: str = "O1", sog_run: str | None = None, skip: tuple = ()):
    """Common Part B entry-point options (defaults = current values, so results do not change):
    --snapshot (default_snapshot()), --orthology {O1,O2}, --sog-run (SOG release, e.g. O1_r3), --run-id (label)."""
    if "snapshot" not in skip:
        ap.add_argument("--snapshot", default=default_snapshot())
    if "orthology" not in skip:
        ap.add_argument("--orthology", default=orthology, choices=["O1", "O2"])
    if "sog_run" not in skip:
        ap.add_argument("--sog-run", default=sog_run, help="SOG release id (results/sogs/<snapshot>/<sog-run>)")
    if "run_id" not in skip:
        ap.add_argument("--run-id", default=None, help="label recorded in sidecars (default: sog-run)")
    return ap
