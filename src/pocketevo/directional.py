"""Cross-fitted directional convergence test (DESIGN v2 §4.2; prereg pairs.directional_convergence).

For a significant pair (a in family F, b in family G) where b belongs to a multi-species SOG S of G:
  * derived clade D = smallest species-tree clade containing b's species that holds >= MIN_DERIVED other members
    of S (from >= 2 species) while leaving >= MIN_OUT members of S outside it; outside members = sister lineages /
    outgroups (ancestral-state proxy)
  * derived members other than b are split (seeded) into a selection half (with b) and a test half
  * selection: delta_sel = mean z(a, selection) - mean z(a, outside) > 0
  * statistic:  delta = mean z(a, test) - mean z(a, outside)       (b's own value never enters)
Calibration: the same statistic for (i) significant same-SOG ortholog pairs (vertical descent) and (ii) random
non-homologous partners of the same b sites; p = empirical upper-tail p against the null; pass iff
delta > null q95 AND BH(q = 0.05) across tested pairs. Conservative p = max over the two nulls.
"""
from __future__ import annotations

import numpy as np

MIN_DERIVED = 4
MIN_DERIVED_SPECIES = 2
MIN_OUT = 2
MAX_SET = 64


# ----------------------------------------------------------------------------------------------- species tree
def parse_newick(s: str):
    """-> (parent array, names list, leaf indices); internal nodes get their label or ''."""
    s = s.strip().rstrip(";")
    parent, names = [], []
    stack, cur, i, tok = [], None, 0, ""

    def new(par):
        parent.append(par); names.append("")
        return len(parent) - 1

    root = new(-1)
    cur = root
    pending_name_for = None
    while i < len(s):
        ch = s[i]
        if ch == "(":
            n = new(cur); stack.append(cur); cur = n
            i += 1
        elif ch == ",":
            if tok:
                names[cur] = tok.split(":")[0]; tok = ""
            cur = new(stack[-1])
            i += 1
        elif ch == ")":
            if tok:
                names[cur] = tok.split(":")[0]; tok = ""
            cur = stack.pop()
            i += 1
        else:
            tok += ch
            i += 1
    if tok:
        names[cur] = tok.split(":")[0]
    parent = np.array(parent)
    is_leaf = np.ones(len(parent), bool)
    is_leaf[parent[parent >= 0]] = False
    return parent, names, np.flatnonzero(is_leaf)


class SpeciesTree:
    def __init__(self, newick: str, leaf_key=lambda name: name):
        self.parent, self.names, leaves = parse_newick(newick)
        self.leaf_of = {leaf_key(self.names[l]): l for l in leaves}
        n = len(self.parent)
        self.anc = []
        for v in range(n):
            path, u = [], v
            while u >= 0:
                path.append(u)
                u = self.parent[u]
            self.anc.append(path)
        self.leafset = [set() for _ in range(n)]
        for l in leaves:
            for u in self.anc[l]:
                self.leafset[u].add(l)


def derived_split(tree: SpeciesTree, b_leaf: int, mem_leaf: np.ndarray, rng) -> tuple | None:
    """Members (by species leaf) -> (derived idx, outside idx) for the smallest qualifying clade above b."""
    for v in tree.anc[b_leaf]:
        ls = tree.leafset[v]
        inside = np.array([l in ls for l in mem_leaf], bool)
        if inside.sum() >= MIN_DERIVED and len(set(mem_leaf[inside].tolist())) >= MIN_DERIVED_SPECIES:
            out = np.flatnonzero(~inside)
            if len(out) < MIN_OUT:
                return None
            return np.flatnonzero(inside), out
    return None


def halves(derived: np.ndarray, rng) -> tuple[np.ndarray, np.ndarray]:
    perm = rng.permutation(derived)
    k = len(perm) // 2
    return perm[k:], perm[:k]          # selection (larger half), test


def cap(idx: np.ndarray, rng, m=MAX_SET):
    return idx if len(idx) <= m else rng.choice(idx, m, replace=False)


# ----------------------------------------------------------------------------------------------- decision
def empirical_p(x: np.ndarray, null: np.ndarray) -> np.ndarray:
    null = np.sort(null[np.isfinite(null)])
    ge = len(null) - np.searchsorted(null, x, side="left")
    return (ge + 1) / (len(null) + 1)


def bh(p: np.ndarray, q: float = 0.05) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return np.zeros(0, bool)
    o = np.argsort(p)
    thr = q * np.arange(1, n + 1) / n
    ok = p[o] <= thr
    k = np.max(np.flatnonzero(ok)) + 1 if ok.any() else 0
    out = np.zeros(n, bool)
    out[o[:k]] = True
    return out


def normal_p(x: np.ndarray, null: np.ndarray) -> np.ndarray:
    """Upper-tail p under a normal fitted to the null (delta is a difference of means -> near-normal); needed
    because empirical p cannot go below 1/(n_null + 1), which makes BH over many tests powerless."""
    from scipy.stats import norm
    null = null[np.isfinite(null)]
    return norm.sf((x - null.mean()) / null.std())


def decide(delta, delta_sel, nulls: list[np.ndarray], q=0.05):
    """pass = selected & delta > every null's q95 & BH on the conservative p (max over nulls of the larger of the
    empirical and normal-fit p in the bulk; normal-fit p beyond the null's range)."""
    ps = []
    for n in nulls:
        pe, pn = empirical_p(delta, n), normal_p(delta, n)
        beyond = delta > np.nanmax(n)
        ps.append(np.where(beyond, pn, np.maximum(pe, pn)))
    p = np.max(np.stack(ps), 0)
    q95 = max(np.nanquantile(n, 0.95) for n in nulls)
    sel = delta_sel > 0
    passed = np.zeros(len(delta), bool)
    idx = np.flatnonzero(sel & np.isfinite(delta))
    if len(idx):
        passed[idx] = bh(p[idx], q) & (delta[idx] > q95)
    return passed, p, q95


# ----------------------------------------------------------------------------------------------- toy world
def _balanced(depth):
    def rec(d, k):
        if d == 0: return f"s{k}", k+1
        l, k = rec(d-1, k); r, k = rec(d-1, k); return f"({l},{r})", k
    return rec(depth, 0)[0] + ";"

def toy_world(seed, planted=False, n_fam=20, d=32, n_a=2000, sb=0.15, noise=0.15, frac=0.5):
    rng = np.random.default_rng(seed)
    tree = SpeciesTree(_balanced(5)); P = tree.parent
    leaves = np.array([tree.leaf_of[f"s{k}"] for k in range(32)])
    nrm = lambda M: M / np.linalg.norm(M, axis=-1, keepdims=True)
    A = nrm(rng.normal(size=(n_a, d)))
    fams = []
    for f in range(n_fam):
        X = np.zeros((len(P), d)); X[0] = rng.normal(size=d)
        for v in range(1, len(P)): X[v] = X[P[v]] + rng.normal(size=d) * sb
        G = X[leaves] + rng.normal(size=(32, d)) * noise
        if planted and f == 0:
            G[:8] = (1 - frac) * nrm(G[:8]) + frac * A[0]
        fams.append(nrm(G))
    def delta(av, G, b, excl=None):
        mem = np.array([m for m in range(32) if m != b and m != excl])
        sp = derived_split(tree, leaves[b], leaves[mem], rng)
        if sp is None: return np.nan, np.nan
        der, outs = mem[sp[0]], mem[sp[1]]
        sel, tst = halves(der, rng)
        z = lambda idx: (G[idx] @ av).mean()
        return z(tst) - z(outs), z(np.r_[sel, b]) - z(outs)
    allS = np.concatenate([A @ G.T for G in fams], 1)
    thr = np.quantile(allS, 0.995)
    D, tag = [], []
    for f, G in enumerate(fams):
        ia, ib = np.nonzero(A @ G.T > thr)
        for x, y in zip(ia, ib):
            D.append(delta(A[x], G, y)); tag.append(f == 0 and x == 0)
    D = np.array(D); tag = np.array(tag)
    N1, N2 = [], []
    for G in fams:
        gg = G @ G.T
        for x in range(32):
            for y in range(32):
                if x != y and gg[x, y] > np.median(gg): N1.append(delta(G[x], G, y, excl=x)[0])
        for y in range(32):
            for x in rng.integers(0, n_a, 30): N2.append(delta(A[x], G, y)[0])
    passed, p, q95 = decide(D[:, 0], D[:, 1], [np.array(N1), np.array(N2)])
    return dict(n=len(D), F_dir=passed.mean(), false_pass=passed[~tag].mean(), power=passed[tag].mean() if tag.any() else np.nan,
                n_planted=int(tag.sum()), q95=q95, meanD=np.nanmean(D[:, 0]), nN1=len(N1), nN2=len(N2))



