"""Part B simulations and null inputs (DESIGN.md §5.1, §5.4, §4.1 M3; prereg part_b.event_fdr, part_a.M3).

Truth models on a real gene tree (leaves = member proteins):
  vertical  : SOG present at the root, never changes (no-event world, prereg "simulated no-event null")
  gain      : absent at the root, one gain on a planted branch c (clade below c present)
  loss      : present at the root, one loss on a planted branch c (clade below c absent)
Observation model per leaf (per-species rates; defaults below until measured rates from
results/observation/... replace them via NoiseModel.from_json):
  true present -> U (u_present) | W (w_present) | A (false_A) | P
  true absent  -> U (u_absent)  | P (false_P)   | W (w_absent) | A
Rates are applied in that order (U first, then the rest among the non-U draws).
Default false_A = prereg states.rescue.false_A_target (0.05, post-rescue target); the other defaults are
placeholders (flagged in SIM.md) and are swept by a noise multiplier.

Null inputs for the whole Part B pipeline (prereg part_b.negative_controls.null_inputs):
  curveball  : SOG x species P/A matrix permuted by a masked curveball (per-SOG and per-species P totals and the
               informative/missing mask preserved), mapped back to members;
  pseudo_clade: member states shuffled among the leaves of the same SOG tree (destroys clade structure);
  fragments  : the same proteins observed twice (replicate observations) as two sister copies of the SOG tree
               -> any gain/loss separating the copies is false (replicate_pair_tree).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict

import numpy as np
from numba import njit

from . import recon
from .recon import STATE_CODES, Tree

P_, W_, A_, U_ = STATE_CODES["P"], STATE_CODES["W"], STATE_CODES["A"], STATE_CODES["U"]


@dataclass
class NoiseModel:
    false_A: float = 0.05
    false_P: float = 0.02
    w_present: float = 0.05
    u_present: float = 0.15
    u_absent: float = 0.50
    w_absent: float = 0.0                               # true absent rescued to W (rescue on true negatives)
    per_species: dict = field(default_factory=dict)     # species -> {rate: value} overrides
    heterogeneity: float = 0.0                          # >0: per-species Beta-distributed rates, concentration
    source: str = "defaults (prereg false_A_target; other rates placeholders)"

    def scaled(self, lam: float) -> "NoiseModel":
        d = asdict(self)
        d["false_A"] = min(1.0, self.false_A * lam)
        d["false_P"] = min(1.0, self.false_P * lam)
        d["per_species"] = {s: {k: (min(1.0, v * lam) if k in ("false_A", "false_P") else v) for k, v in r.items()}
                            for s, r in self.per_species.items()}
        d["source"] = f"{self.source}; false_A/false_P x {lam}"
        return NoiseModel(**d)

    @classmethod
    def from_json(cls, path) -> "NoiseModel":
        """Measured rates (observation agent): {"global": {...}, "per_species": {species: {...}}, "source": str}."""
        d = json.loads(open(path).read())
        g = d.get("global", {})
        return cls(**{k: g[k] for k in ("false_A", "false_P", "w_present", "u_present", "u_absent", "w_absent")
                      if k in g},
                   per_species=d.get("per_species", {}), source=d.get("source", str(path)))

    def rates(self, species: list, rng: np.random.Generator | None = None) -> dict:
        """Per-leaf rate arrays for the leaf species list."""
        keys = ("false_A", "false_P", "w_present", "u_present", "u_absent", "w_absent")
        uniq = sorted(set(species))
        table = {}
        for s in uniq:
            r = {k: getattr(self, k) for k in keys}
            r.update(self.per_species.get(s, {}))
            if self.heterogeneity > 0 and rng is not None:
                c = self.heterogeneity
                r = {k: (float(rng.beta(max(v * c, 1e-3), max((1 - v) * c, 1e-3))) if 0 < v < 1 else v)
                     for k, v in r.items()}
            table[s] = r
        return {k: np.array([table[s][k] for s in species]) for k in keys}


def observe(truth: np.ndarray, rates: dict, rng: np.random.Generator) -> np.ndarray:
    """truth (n_leaves, R) {0,1} -> observed state codes (n_leaves, R)."""
    n, R = truth.shape
    col = lambda k: rates[k][:, None]  # noqa: E731
    u = rng.random((n, R))
    x = rng.random((n, R))
    out = np.full((n, R), P_, dtype=np.int8)
    pres = truth == 1
    # true present
    isU = pres & (u < col("u_present"))
    isW = pres & ~isU & (x < col("w_present"))
    isA = pres & ~isU & ~isW & (x < col("w_present") + col("false_A"))
    out[isU] = U_
    out[isW] = W_
    out[isA] = A_
    # true absent
    ab = ~pres
    aU = ab & (u < col("u_absent"))
    aP = ab & ~aU & (x < col("false_P"))
    aW = ab & ~aU & ~aP & (x < col("false_P") + col("w_absent"))
    out[ab] = A_
    out[aU] = U_
    out[aP] = P_
    out[aW] = W_
    return out


def eligible_branches(tree: Tree, tix, typ: str, min_species: int = 3, min_out_lineages: int = 2,
                      min_out_species: int = 2) -> np.ndarray:
    """Non-root nodes whose clade has >= min_species species, with >= min_out_lineages outgroup lineages and
    >= min_out_species species outside (the event is detectable in principle under the prereg rules)."""
    all_sp = len(tix.sp_names)
    out = []
    for c in range(1, tree.n_nodes):
        if tix.n_species_below[c] < min_species:
            continue
        lin = tix.lineages(c)
        if len(lin) < min_out_lineages:
            continue
        lo, hi = tree.leaf_lo[c], tree.leaf_hi[c]
        outside = np.concatenate([tix.leaf_sp[:lo], tix.leaf_sp[hi:]])
        if len(np.unique(outside)) < min_out_species:
            continue
        out.append(c)
    _ = all_sp
    return np.array(out, dtype=np.int32)


def simulate_truth(tree: Tree, R: int, scenario: str, rng: np.random.Generator, branches: np.ndarray | None = None):
    """Return (truth (n_leaves, R) int8, planted node per replicate (-1 for vertical))."""
    n = tree.n_leaves
    if scenario == "vertical":
        return np.ones((n, R), dtype=np.int8), np.full(R, -1, dtype=np.int32)
    if branches is None or len(branches) == 0:
        raise ValueError("no eligible branch to plant an event")
    planted = rng.choice(branches, size=R)
    base, inside = (0, 1) if scenario == "gain" else (1, 0)
    truth = np.full((n, R), base, dtype=np.int8)
    for r, c in enumerate(planted):
        truth[tree.leaf_lo[c]:tree.leaf_hi[c], r] = inside
    return truth, planted.astype(np.int32)


# --------------------------------------------------------------------------------------------- curveball
@njit(cache=True)
def _curveball(X, n_swaps, seed):
    """Masked curveball (Strona et al. 2014) on rows of X (1 present, 0 absent, <0 fixed missing).
    Trades between random row pairs among columns informative in both rows where they differ.
    Preserves row sums, column sums of 1s and the missing mask. Stops when the cumulative number of exchanged
    pairs (min(|A|,|B|) per trade; = equivalent checkerboard swaps) reaches n_swaps, or after a trade cap."""
    np.random.seed(seed)
    m, n = X.shape
    Y = X.copy()
    done = 0
    trades = 0
    cap = 50 * n_swaps + 1000
    a_idx = np.empty(n, dtype=np.int64)
    b_idx = np.empty(n, dtype=np.int64)
    pool = np.empty(n, dtype=np.int64)
    while done < n_swaps and trades < cap:
        trades += 1
        i = np.random.randint(m)
        j = np.random.randint(m)
        if i == j:
            continue
        na = 0
        nb = 0
        for k in range(n):
            yi = Y[i, k]
            yj = Y[j, k]
            if yi == 1 and yj == 0:
                a_idx[na] = k
                na += 1
            elif yi == 0 and yj == 1:
                b_idx[nb] = k
                nb += 1
        if na == 0 or nb == 0:
            continue
        tot = na + nb
        for k in range(na):
            pool[k] = a_idx[k]
        for k in range(nb):
            pool[na + k] = b_idx[k]
        # Fisher-Yates shuffle of pool[:tot]
        for k in range(tot - 1, 0, -1):
            r = np.random.randint(k + 1)
            tmp = pool[k]
            pool[k] = pool[r]
            pool[r] = tmp
        for k in range(tot):
            c = pool[k]
            if k < na:
                Y[i, c] = 1
                Y[j, c] = 0
            else:
                Y[i, c] = 0
                Y[j, c] = 1
        done += na if na < nb else nb
    return Y, done, trades


def curveball(X: np.ndarray, n_swaps: int | None = None, seed: int = 0, swaps_per_one: int = 100):
    """prereg part_a.M3 nulls.curveball: n_swaps default = 100 x number of 1s."""
    X = np.ascontiguousarray(np.asarray(X, dtype=np.int8))
    if n_swaps is None:
        n_swaps = int(swaps_per_one * int((X == 1).sum()))
    Y, done, trades = _curveball(X, np.int64(n_swaps), np.int64(seed))
    return Y, int(done), int(trades)


# --------------------------------------------------------------------------------------------- SOG matrices
def species_matrix(sog_leaf_codes: list, sog_leaf_species: list, species_order: list, w_present: bool = True):
    """Collapse member states to SOG x species (strongest state P>W>A>U). Returns int8 matrix:
    1 present (P, and W if w_present), 0 absent (A), -1 uninformative (U or W if not w_present), -2 no member (G)."""
    col = {s: i for i, s in enumerate(species_order)}
    M = np.full((len(sog_leaf_codes), len(species_order)), -2, dtype=np.int8)
    rank = {P_: 3, W_: 2, A_: 1, U_: 0}
    for i, (codes, sp) in enumerate(zip(sog_leaf_codes, sog_leaf_species)):
        best = {}
        for c, s in zip(codes, sp):
            if s not in col:
                continue
            if s not in best or rank[int(c)] > rank[best[s]]:
                best[s] = int(c)
        for s, c in best.items():
            if c == P_ or (c == W_ and w_present):
                M[i, col[s]] = 1
            elif c == A_:
                M[i, col[s]] = 0
            else:
                M[i, col[s]] = -1
    return M


def remap_members(sog_leaf_codes: list, sog_leaf_species: list, species_order: list, M_orig: np.ndarray,
                  M_perm: np.ndarray) -> list:
    """Map a permuted species matrix back to members: where a species cell flipped 1->0 all informative P/W
    members become A; 0->1, all A members become P. Unchanged cells keep member states."""
    col = {s: i for i, s in enumerate(species_order)}
    out = []
    for i, (codes, sp) in enumerate(zip(sog_leaf_codes, sog_leaf_species)):
        c = np.array(codes, dtype=np.int8).copy()
        for k, s in enumerate(sp):
            j = col.get(s)
            if j is None or M_orig[i, j] == M_perm[i, j]:
                continue
            if M_perm[i, j] == 0 and c[k] in (P_, W_):
                c[k] = A_
            elif M_perm[i, j] == 1 and c[k] == A_:
                c[k] = P_
        out.append(c)
    return out


def pseudo_clades(codes: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle member states among the leaves of the tree (per column)."""
    codes = np.array(codes, copy=True)
    if codes.ndim == 1:
        return rng.permutation(codes)
    for r in range(codes.shape[1]):
        codes[:, r] = rng.permutation(codes[:, r])
    return codes


def replicate_pair_tree(tree: Tree) -> Tree:
    """Same-protein replicate null: root (level 'replicate') with two copies of the SOG tree, the second holding
    a replicate observation of the same proteins (other fragment window / model source / seed; leaf names
    '<acc>|rep'). Leaf order = original leaves then replicate leaves. Any gain/loss called on the stem of either
    copy (nodes 1 and 1 + tree.n_nodes) separates a protein set from itself and is false by construction."""
    n = tree.n_nodes
    parent = np.concatenate([[-1], tree.parent + 1, tree.parent + 1 + n]).astype(np.int32)
    parent[1] = 0
    parent[1 + n] = 0
    kind = np.concatenate([[recon.DUP], tree.kind, tree.kind])
    level = ["replicate"] + list(tree.level) * 2
    name = ["replicate_root"] + list(tree.node_name) + [x + "|rep" for x in tree.node_name]
    species = [""] + list(tree.species) * 2
    taxid = np.concatenate([[-1], tree.taxid, tree.taxid])
    return Tree(tree.tree_id + "|reppairs", parent, kind, level, name, species, taxid)


@njit(cache=True)
def _checkerboard(X, n_swaps, max_attempts, seed):
    """Masked checkerboard swaps. Cells of the minority value v (0 or 1) are sampled in pairs: (i,b), (j,a) both v
    with (i,a), (j,b) both 1-v (informative) -> swap the 2x2. Preserves row/column sums and the missing (<0) mask;
    O(1) per attempt with high acceptance when informative cells are mostly of one value."""
    np.random.seed(seed)
    Y = X.copy()
    m, n = Y.shape
    n1 = 0
    n0 = 0
    for i in range(m):
        for c in range(n):
            if Y[i, c] == 1:
                n1 += 1
            elif Y[i, c] == 0:
                n0 += 1
    v = 0 if n0 <= n1 else 1
    k = n0 if v == 0 else n1
    ri = np.empty(k, dtype=np.int64)
    ci = np.empty(k, dtype=np.int64)
    k = 0
    for i in range(m):
        for c in range(n):
            if Y[i, c] == v:
                ri[k] = i
                ci[k] = c
                k += 1
    done = 0
    att = 0
    if k < 2:
        return Y, done, att
    w = 1 - v
    while done < n_swaps and att < max_attempts:
        att += 1
        p = np.random.randint(k)
        q = np.random.randint(k)
        i, b, j, a = ri[p], ci[p], ri[q], ci[q]
        if i == j or a == b:
            continue
        if Y[i, a] != w or Y[j, b] != w:
            continue
        Y[i, b] = w
        Y[i, a] = v
        Y[j, a] = w
        Y[j, b] = v
        ci[p] = a
        ci[q] = b
        done += 1
    return Y, done, att


def checkerboard(X: np.ndarray, n_swaps: int | None = None, seed: int = 0, swaps_per_one: int = 100,
                 attempts_per_swap: int = 200):
    """prereg part_a.M3 curveball null, swap form: n_swaps = 100 x number of 1s successful checkerboard swaps
    (attempt cap attempts_per_swap x n_swaps; achieved count returned)."""
    X = np.ascontiguousarray(np.asarray(X, dtype=np.int8))
    if n_swaps is None:
        n_swaps = int(swaps_per_one * int((X == 1).sum()))
    Y, done, att = _checkerboard(X, np.int64(n_swaps), np.int64(attempts_per_swap * n_swaps), np.int64(seed))
    return Y, int(done), int(att)
