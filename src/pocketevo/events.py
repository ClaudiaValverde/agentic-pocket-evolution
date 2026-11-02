"""Part B event extraction, candidate calling, specificity test, FDR utilities and funnels (DESIGN.md §5.1–5.4).

Definitions (prereg part_b; interpretation choices are marked INTERP and listed in results/partB/sim/SIM.md):
- An *event* is a branch parent->child of the gene tree whose MPR states differ (0->1 gain, 1->0 loss).
  It is *unambiguous* if both endpoint states are identical in every most-parsimonious reconstruction (MPR).
  Its *margin* = extra parsimony cost of the best reconstruction that flips either endpoint (support).
- Clade C/D = all leaves (member proteins) below the child node. Informative members = P, W, A.
- Outgroup *lineages* (INTERP): for an outgroup member m, its lineage is the child subtree of LCA(m, clade node)
  that contains m (sister lineage at the first ancestor, then successive outgroups). Distinct lineages = distinct
  such subtrees.
- Gain candidate: unambiguous 0->1 on C's stem; P in >= 3 species of C; (P+W)/informative(C) >= 0.80;
  >= 2 A members outside C from >= 2 distinct outgroup lineages.  W never counts toward the species threshold
  (INTERP of "present-ambiguous": W supports presence but is not P evidence).
- Loss candidate: unambiguous 1->0 on D's stem; P in >= 2 outgroup lineages; A in >= 2 species of D;
  A/informative(D) >= 0.80 (W in D counts against the loss).
- Calibrated thresholds (event FDR 10%): additional stricter cuts on (min species in derived clade, min outgroup
  lineages supporting the ancestral state, min margin);
  chosen in scripts/partB/run_sim.py, stored in results/partB/sim/thresholds.json. The frozen prereg minima are
  never relaxed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict

import numpy as np

from . import recon
from .recon import STATE_CODES, Tree

P_, W_, A_, U_ = STATE_CODES["P"], STATE_CODES["W"], STATE_CODES["A"], STATE_CODES["U"]
DEEP_LEVELS = {"cellular organisms", "Eukaryota", "Bacteria", "Archaea"}


@dataclass(frozen=True)
class Rules:
    gain_min_species: int = 3
    gain_min_frac_present: float = 0.80
    gain_min_absent_members_out: int = 2
    gain_min_out_lineages: int = 2
    loss_min_out_lineages_present: int = 2
    loss_min_species_absent: int = 2
    loss_min_frac_absent: float = 0.80
    # calibrated (FDR) extra cuts; defaults = no extra cut
    gain_min_species_cal: int = 3
    loss_min_species_cal: int = 2
    gain_min_margin: float = 0.0
    loss_min_margin: float = 0.0
    gain_min_out_lineages_cal: int = 2
    loss_min_out_lineages_cal: int = 2
    require_unambiguous: bool = True

    @classmethod
    def from_prereg(cls, pr: dict | None = None, calibrated: dict | None = None) -> "Rules":
        pr = pr or recon.load_prereg()
        g, l = pr["part_b"]["gain"], pr["part_b"]["loss"]
        kw = dict(gain_min_species=g["min_species_present_in_clade"]["value"],
                  gain_min_frac_present=g["min_informative_fraction_present"]["value"],
                  gain_min_absent_members_out=g["min_absent_members_outgroup"]["value"],
                  gain_min_out_lineages=g["min_outgroup_lineages"]["value"],
                  loss_min_out_lineages_present=l["min_outgroup_lineages_present"]["value"],
                  loss_min_species_absent=l["min_species_absent_in_derived_clade"]["value"],
                  loss_min_frac_absent=l["min_informative_fraction_absent"]["value"])
        kw["gain_min_species_cal"] = kw["gain_min_species"]
        kw["loss_min_species_cal"] = kw["loss_min_species_absent"]
        kw["gain_min_out_lineages_cal"] = kw["gain_min_out_lineages"]
        kw["loss_min_out_lineages_cal"] = kw["loss_min_out_lineages_present"]
        if calibrated:
            for k, v in calibrated.items():
                if k in cls.__dataclass_fields__:
                    kw[k] = v
        r = cls(**kw)
        r.check_not_relaxed()
        return r

    def check_not_relaxed(self):
        if (self.gain_min_species_cal < self.gain_min_species or self.loss_min_species_cal < self.loss_min_species_absent
                or self.gain_min_out_lineages_cal < self.gain_min_out_lineages
                or self.loss_min_out_lineages_cal < self.loss_min_out_lineages_present):
            raise ValueError("calibrated thresholds may not relax the frozen prereg minima")

    def with_cal(self, **kw) -> "Rules":
        d = asdict(self)
        d.update(kw)
        r = Rules(**d)
        r.check_not_relaxed()
        return r


class TreeIndex:
    """Static per-tree helpers: species codes per leaf, outgroup lineage intervals per node (cached)."""

    def __init__(self, tree: Tree, species_tree: recon.SpeciesTree | None = None):
        self.t = tree
        sp = tree.leaf_species()
        self.sp_names = sorted(set(sp))
        code = {s: i for i, s in enumerate(self.sp_names)}
        self.leaf_sp = np.array([code[s] for s in sp], dtype=np.int32)
        self.st = species_tree
        # number of distinct species below each node (static)
        n = tree.n_nodes
        self.n_species_below = np.zeros(n, dtype=np.int32)
        for v in range(n):
            lo, hi = tree.leaf_lo[v], tree.leaf_hi[v]
            self.n_species_below[v] = len(np.unique(self.leaf_sp[lo:hi]))
        self._lin = {}

    def lineages(self, c: int):
        """Outgroup lineages of node c: list of (ancestor, lineage_root, lo, hi)."""
        if c in self._lin:
            return self._lin[c]
        t = self.t
        out = []
        on_path = c
        for a in t.ancestors(c):
            for x in t.children(a):
                if x != on_path:
                    out.append((a, int(x), int(t.leaf_lo[x]), int(t.leaf_hi[x])))
            on_path = a
        self._lin[c] = out
        return out


def _prefix(mask: np.ndarray) -> np.ndarray:
    """mask (n_leaves, R) -> prefix sums (n_leaves+1, R)."""
    out = np.zeros((mask.shape[0] + 1, mask.shape[1]), dtype=np.int32)
    np.cumsum(mask, axis=0, out=out[1:])
    return out


def extract_events(tree: Tree, codes: np.ndarray, res: recon.SankoffResult, rules: Rules,
                   tix: TreeIndex | None = None, confidence: np.ndarray | None = None,
                   char_ids=None, prefilter: bool = True) -> list:
    """All MPR events of each character (column of codes) with clade statistics and rule flags.

    prefilter=True skips events whose clade cannot reach the frozen species minimum (counted as 'n_small').
    Returns list of dicts; the special key 'n_events_char' gives the MPR event count of that character.
    confidence: optional (n_leaves, R) state confidence (prereg candidate_score); default P/A 1.0, W 0.5.
    """
    tix = tix or TreeIndex(tree)
    t = tree
    codes = np.asarray(codes)
    if codes.ndim == 1:
        codes = codes[:, None]
    R = codes.shape[1]
    mpr = res.mpr
    par = t.parent
    child_nodes = np.arange(1, t.n_nodes)
    diff = mpr[child_nodes] != mpr[par[child_nodes]]          # (n-1, R)
    n_ev_char = diff.sum(axis=0)
    # events whose clade spans >= 2 species (evidence filter f1 "single tree-consistent event")
    n_clade_ev_char = (diff & (tix.n_species_below[1:] >= 2)[:, None]).sum(axis=0)
    if not diff.any():
        return []
    unamb = res.unambiguous
    margin = res.margin()
    isP, isW, isA = codes == P_, codes == W_, codes == A_
    isInf = isP | isW | isA
    cP, cW, cA, cI = _prefix(isP), _prefix(isW), _prefix(isA), _prefix(isInf)
    if confidence is None:
        confidence = np.where(isW, 0.5, 1.0)
    out = []
    ev_i, ev_r = np.nonzero(diff)
    for i, r in zip(ev_i.tolist(), ev_r.tolist()):
        c = i + 1
        p = int(par[c])
        typ = "gain" if mpr[c, r] == 1 else "loss"
        nsp = int(tix.n_species_below[c])
        if prefilter and nsp < (rules.gain_min_species if typ == "gain" else rules.loss_min_species_absent):
            continue
        lo, hi = int(t.leaf_lo[c]), int(t.leaf_hi[c])
        sp_slice = tix.leaf_sp[lo:hi]
        nP = int(cP[hi, r] - cP[lo, r]); nW = int(cW[hi, r] - cW[lo, r]); nA = int(cA[hi, r] - cA[lo, r])
        nI = nP + nW + nA
        spP = np.unique(sp_slice[isP[lo:hi, r]])
        spA = np.unique(sp_slice[isA[lo:hi, r]])
        lin = tix.lineages(c)
        outA_members = outP_members = 0
        linA = linP = 0
        outA_sp, outP_sp = set(), set()
        for (_a, x, llo, lhi) in lin:
            a_ = int(cA[lhi, r] - cA[llo, r]); p_ = int(cP[lhi, r] - cP[llo, r])
            outA_members += a_; outP_members += p_
            if a_ > 0:
                linA += 1
                outA_sp.update(tix.leaf_sp[llo:lhi][isA[llo:lhi, r]].tolist())
            if p_ > 0:
                linP += 1
                outP_sp.update(tix.leaf_sp[llo:lhi][isP[llo:lhi, r]].tolist())
        un = bool(unamb[c, r] and unamb[p, r])
        mg = float(min(margin[c, r], margin[p, r])) if p >= 0 else float(margin[c, r])
        conf_clade = confidence[lo:hi, r][isInf[lo:hi, r]]
        rec = dict(char=(char_ids[r] if char_ids is not None else r), r=int(r), node=c, parent=p, type=typ,
                   node_name=t.node_name[c], level=t.level[c], parent_level=t.level[p],
                   parent_kind=recon.KIND_NAME[int(t.kind[p])], node_kind=recon.KIND_NAME[int(t.kind[c])],
                   unambiguous=un, margin=mg, parent_is_root=bool(p == 0), n_events_char=int(n_ev_char[r]),
                   n_clade_events_char=int(n_clade_ev_char[r]),
                   n_members=hi - lo, n_species=nsp, n_inf=nI, n_P=nP, n_W=nW, n_A=nA,
                   n_species_P=len(spP), n_species_A=len(spA),
                   frac_present=(nP + nW) / nI if nI else 0.0, frac_absent=nA / nI if nI else 0.0,
                   out_A_members=outA_members, out_A_lineages=linA, out_P_members=outP_members,
                   out_P_lineages=linP, n_out_lineages=len(lin),
                   clade_species=";".join(tix.sp_names[s] for s in np.unique(sp_slice)),
                   out_A_species=";".join(tix.sp_names[s] for s in sorted(outA_sp)),
                   out_P_species=";".join(tix.sp_names[s] for s in sorted(outP_sp)),
                   min_conf_clade=float(conf_clade.min()) if conf_clade.size else 0.0,
                   deep=bool(t.level[c] in DEEP_LEVELS or t.level[p] in DEEP_LEVELS))
        if tix.st is not None:
            taxids = [recon.taxid_from_species_name(tix.sp_names[s]) for s in np.unique(sp_slice)]
            v = tix.st.lca_of_taxids(taxids)
            rec["clade_lca_species_tree"] = tix.st.names[v] if v >= 0 else ""
        out.append(rec)
    return out


def apply_rules(ev: dict, rules: Rules) -> dict:
    """Add rule-pass flags and 'candidate' (frozen rules) / 'called' (frozen + calibrated) to one event."""
    f = {}
    f["unambiguous"] = ev["unambiguous"] or not rules.require_unambiguous
    if ev["type"] == "gain":
        f["species"] = ev["n_species_P"] >= rules.gain_min_species
        f["frac"] = ev["frac_present"] >= rules.gain_min_frac_present - 1e-12
        f["outgroup"] = (ev["out_A_members"] >= rules.gain_min_absent_members_out and
                         ev["out_A_lineages"] >= rules.gain_min_out_lineages)
        f["cal_species"] = ev["n_species_P"] >= rules.gain_min_species_cal
        f["cal_margin"] = ev["margin"] >= rules.gain_min_margin - 1e-12
        f["cal_lineages"] = ev["out_A_lineages"] >= rules.gain_min_out_lineages_cal
    else:
        f["species"] = ev["n_species_A"] >= rules.loss_min_species_absent
        f["frac"] = ev["frac_absent"] >= rules.loss_min_frac_absent - 1e-12
        f["outgroup"] = ev["out_P_lineages"] >= rules.loss_min_out_lineages_present
        f["cal_species"] = ev["n_species_A"] >= rules.loss_min_species_cal
        f["cal_margin"] = ev["margin"] >= rules.loss_min_margin - 1e-12
        f["cal_lineages"] = ev["out_P_lineages"] >= rules.loss_min_out_lineages_cal
    ev = dict(ev)
    for k, v in f.items():
        ev["pass_" + k] = bool(v)
    ev["candidate"] = bool(f["unambiguous"] and f["species"] and f["frac"] and f["outgroup"])
    ev["called"] = bool(ev["candidate"] and f["cal_species"] and f["cal_margin"] and f["cal_lineages"])
    ev["f1_single_event"] = bool(ev["n_clade_events_char"] == 1)
    ev["called_f1"] = bool(ev["called"] and ev["f1_single_event"])
    n_pers = ev["n_species_P"] if ev["type"] == "gain" else ev["n_species_A"]
    ev["score_S_partial"] = float(n_pers * ev["min_conf_clade"])   # x evidence_count (filters 2/3) later
    return ev


FUNNEL_STAGES = ["mpr_event", "clade_size_ok", "unambiguous", "species", "frac", "outgroup", "candidate",
                 "called"]


def funnel(events: list, n_chars: int) -> dict:
    """Counts per stage (cumulative AND) per event type; events must have been through apply_rules."""
    out = {"n_chars": int(n_chars)}
    for typ in ("gain", "loss"):
        e = [x for x in events if x["type"] == typ]
        row = {"clade_size_ok": len(e)}
        keep = e
        for k in ("unambiguous", "species", "frac", "outgroup"):
            keep = [x for x in keep if x["pass_" + k]]
            row[k] = len(keep)
        row["candidate"] = sum(x["candidate"] for x in e)
        row["called"] = sum(x["called"] for x in e)
        row["called_f1_single_event"] = sum(x["called_f1"] for x in e)
        row["chars_with_called"] = len({x["char"] for x in e if x["called"]})
        out[typ] = row
    return out


def yield_preflight(events: list, min_candidates: int = 5) -> dict:
    """prereg pi_freeze.yield_preflight: < 5 gain+loss candidates in the primary stratum -> secondary stratum."""
    n = sum(1 for e in events if e.get("called"))
    return {"n_gain_loss_called": n, "threshold": min_candidates, "run_secondary_stratum": n < min_candidates,
            "rule": "pi_freeze.yield_preflight (pLDDT>=70 stratum = 'secondary stratum')"}


# --------------------------------------------------------------------------------------------- statistics
def bh(p: np.ndarray, alpha: float = 0.05):
    """Benjamini-Hochberg: returns (q-values, reject mask)."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([]), np.array([], dtype=bool)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    qq = np.empty(n)
    qq[o] = np.minimum(q, 1.0)
    return qq, qq <= alpha


def mixture_fdr(f0: np.ndarray, f1: np.ndarray, t1: np.ndarray, pi: float) -> np.ndarray:
    """FDR of a world where a fraction pi of SOGs carry a planted event.

    f0 = false calls per null SOG, f1 = false calls per planted SOG, t1 = true calls per planted SOG.
    """
    false = (1 - pi) * np.asarray(f0) + pi * np.asarray(f1)
    tot = false + pi * np.asarray(t1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(tot > 0, false / tot, 0.0)


def empirical_fdr(null_calls_per_sog: float, n_real_sogs: int, n_real_calls: int) -> float:
    """Real-data FDR estimate from the simulated no-event null: E[false calls] / observed calls (capped at 1)."""
    if n_real_calls == 0:
        return 0.0
    return float(min(1.0, null_calls_per_sog * n_real_sogs / n_real_calls))


# --------------------------------------------------------------------------------------------- specificity
_BLOSUM = None


def blosum62():
    global _BLOSUM
    if _BLOSUM is None:
        from Bio.Align import substitution_matrices  # noqa: PLC0415
        _BLOSUM = substitution_matrices.load("BLOSUM62")
    return _BLOSUM


def clade_consistent_substitutions(res_d: list, res_s: list, max_blosum: float = -1e-9, min_members: int = 2):
    """Footprint columns where all non-gap D members share residue a, all non-gap sister members share b != a,
    and BLOSUM62(a, b) < 0. res_d / res_s: lists of equal-length strings (one char per footprint column, '-' gap).
    Returns list of (column_index, a, b, score)."""
    if not res_d or not res_s:
        return []
    L = len(res_d[0])
    B = blosum62()
    out = []
    for j in range(L):
        a = {s[j] for s in res_d if s[j] not in "-.X"}
        b = {s[j] for s in res_s if s[j] not in "-.X"}
        nd = sum(1 for s in res_d if s[j] not in "-.X")
        ns = sum(1 for s in res_s if s[j] not in "-.X")
        if len(a) == 1 and len(b) == 1 and nd >= min_members and ns >= 1:
            aa, bb = next(iter(a)), next(iter(b))
            if aa != bb:
                sc = float(B[aa][bb])
                if sc <= max_blosum:
                    out.append((j, aa, bb, sc))
    return out


def residualize(fp: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Regress fingerprints (n, d) on covariates (n, k) (+intercept) and return residuals + global mean
    (keeps the cosine geometry; prereg regress_out: volume, footprint pLDDT, model source)."""
    fp = np.asarray(fp, dtype=float)
    if cov is None or cov.size == 0:
        return fp
    X = np.column_stack([np.ones(len(fp)), np.asarray(cov, dtype=float)])
    beta, *_ = np.linalg.lstsq(X, fp, rcond=None)
    return fp - X @ beta + fp.mean(axis=0, keepdims=True)


def shift_stat(fp_d: np.ndarray, fp_s: np.ndarray) -> float:
    """1 - cosine(mean fingerprint of D, mean fingerprint of sister)."""
    a, b = np.asarray(fp_d).mean(axis=0), np.asarray(fp_s).mean(axis=0)
    den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(1.0 - (a @ b) / den) if den > 0 else 0.0


def noise_floor(replicate_cosines: np.ndarray, quantile: float = 0.95) -> float:
    """Technical noise floor in shift units: quantile of (1 - same-site replicate cosine)."""
    return float(np.quantile(1.0 - np.asarray(replicate_cosines, dtype=float), quantile))


def sister_pairs(tree: Tree, codes: np.ndarray, min_members: int = 2, min_frac_present: float = 0.8):
    """Candidate (D, sister) pairs at every internal node where the SOG is present throughout (INTERP: in each of
    D and S >= min_frac_present of informative members are P/W, mirroring the gain rule; A members are left out
    of the comparison). Binary node -> one pair (D = smaller child); polytomy -> each child vs the rest.
    Yields (node, D child, D leaf ranks, S leaf ranks)."""
    t = tree
    codes = np.asarray(codes)
    pres = np.isin(codes, (P_, W_))
    inf = pres | (codes == A_)
    for v in range(t.n_nodes):
        ch = t.children(v)
        if len(ch) < 2:
            continue
        groups, ok = [], []
        for c in ch:
            g = np.arange(t.leaf_lo[c], t.leaf_hi[c])
            ni = int(inf[g].sum())
            ok.append(ni == 0 or pres[g].sum() / ni >= min_frac_present - 1e-12)
            groups.append(g[pres[g]])
        if len(ch) == 2:
            pairs = [(0, 1)] if len(groups[0]) <= len(groups[1]) else [(1, 0)]
        else:
            pairs = [(i, -1) for i in range(len(ch))]
        for i, j in pairs:
            d = groups[i]
            if j >= 0:
                s, s_ok = groups[j], ok[j]
            else:
                s = np.concatenate([g for k, g in enumerate(groups) if k != i])
                rest = np.concatenate([np.arange(t.leaf_lo[c], t.leaf_hi[c]) for k, c in enumerate(ch) if k != i])
                ni = int(inf[rest].sum())
                s_ok = ni == 0 or pres[rest].sum() / ni >= min_frac_present - 1e-12
            if ok[i] and s_ok and len(d) >= min_members and len(s) >= 1:
                yield v, int(ch[i]), d, s


def specificity_tests(tree: Tree, codes: np.ndarray, residues: dict | None, fingerprints: dict | None,
                      covariates: dict | None, floor: float, min_subst: int = 2, min_members: int = 2,
                      species_tree: recon.SpeciesTree | None = None, sog_id: str = ""):
    """Residue-first specificity statistics for every sister pair of one SOG.

    residues: accession -> aligned footprint residue string; fingerprints: accession -> vector;
    covariates: accession -> vector (volume, pLDDT, model_source code). Returns list of dict rows with the shift
    statistic (NaN if fingerprints missing); p-values are assigned across SOGs by `empirical_null_pvalues`.
    """
    names = tree.leaf_names()
    rows = []
    for v, c, d, s in sister_pairs(tree, codes, min_members):
        acc_d, acc_s = [names[i] for i in d], [names[i] for i in s]
        sub = []
        if residues:
            rd = [residues[a] for a in acc_d if a in residues]
            rs = [residues[a] for a in acc_s if a in residues]
            sub = clade_consistent_substitutions(rd, rs, min_members=min_members)
        stat = float("nan")
        if fingerprints:
            have = [a for a in acc_d + acc_s if a in fingerprints]
            if len(set(have) & set(acc_d)) >= min_members and set(have) & set(acc_s):
                F = np.array([fingerprints[a] for a in have], dtype=float)
                if covariates:
                    F = residualize(F, np.array([covariates[a] for a in have], dtype=float))
                isd = np.array([a in set(acc_d) for a in have])
                stat = shift_stat(F[isd], F[~isd])
        lev = tree.level[v]
        depth = species_tree.level_depth(lev) if species_tree is not None else -1
        rows.append(dict(sog_id=sog_id, node=int(v), d_node=int(c), level=lev, depth=int(depth),
                         node_kind=recon.KIND_NAME[int(tree.kind[v])], n_d=len(d), n_s=len(s),
                         n_species_d=len({tree.species[tree.leaves[i]] for i in d}),
                         n_subst=len(sub), substitutions=";".join(f"{j}:{a}>{b}" for j, a, b, _ in sub),
                         residue_pass=len(sub) >= min_subst, shift=stat,
                         shift_pass=bool(stat > floor) if not math.isnan(stat) else False))
    return rows


def empirical_null_pvalues(rows: list, min_null: int = 20) -> np.ndarray:
    """p-value of each residue-gated row's shift against the empirical null (prereg pi_freeze.specificity_null):
    shifts of the *other* sister-clade pairs (those without the residue-first signal, n_subst < min) of matched
    depth (species-tree depth of the node level) and matched derived-clade size (log2 bin). The match is widened
    to +-1, +-2, any when fewer than `min_null` comparators exist. Non-gated rows get NaN (not tested).
    Implemented with per-(depth, size) sorted pools (searchsorted), so it scales to millions of pairs."""
    shift = np.array([r["shift"] for r in rows], dtype=float)
    depth = np.array([r["depth"] for r in rows], dtype=int)
    size = np.floor(np.log2(np.maximum(1, np.array([r["n_d"] for r in rows])))).astype(int)
    gate = np.array([r["residue_pass"] for r in rows], dtype=bool)
    p = np.full(len(rows), np.nan)
    have = ~np.isnan(shift)
    pool = have & ~gate
    bins = {}
    for d_, s_, x in zip(depth[pool], size[pool], shift[pool]):
        bins.setdefault((int(d_), int(s_)), []).append(x)
    bins = {k: np.sort(np.array(v)) for k, v in bins.items()}
    allpool = np.sort(shift[pool])
    for i in np.flatnonzero(have & gate):
        x = shift[i]
        for w in (0, 1, 2, 99):
            if w == 99:
                arrs = [allpool]
            else:
                arrs = [v for (d_, s_), v in bins.items() if abs(d_ - depth[i]) <= w and abs(s_ - size[i]) <= w]
            nn = sum(len(v) for v in arrs)
            if nn >= min_null or w == 99:
                break
        ge = sum(len(v) - int(np.searchsorted(v, x, side="left")) for v in arrs)
        p[i] = (1 + ge) / (1 + nn)
    return p


def _gpd_pwm(y: np.ndarray):
    """GPD (shape xi, scale sc) by probability-weighted moments (Hosking & Wallis 1987), as nulls.gpd_pwm_rows."""
    y = np.sort(np.asarray(y, dtype=float))
    m = len(y)
    a0 = y.mean()
    pk = (np.arange(1, m + 1) - 0.35) / m
    a1 = np.mean((1 - pk) * y)
    den = a0 - 2 * a1
    return 2 - a0 / den, 2 * a0 * a1 / den


def specificity_tail_pvalues(rows: list, min_null: int = 20, q_u: float = 0.95, min_exc: int = 20) -> dict:
    """D34 sensitivity (post hoc): tail-fitted p-values for the residue-gated rows, with the SAME matched comparator
    pool as `empirical_null_pvalues` (depth x log2 clade-size bins, widened +-1, +-2, any to >= min_null).
    p_normal (D17 rule, src/pocketevo/directional.decide): normal fitted to the pool; beyond the pool's maximum
    the normal-fit p, else max(empirical, normal). p_gpd: GPD fitted by PWM to the exceedances over the pool's
    q95 (>= min_exc exceedances, else the normal-fit p); beyond the pool maximum the GPD p, else
    max(empirical, GPD) above u and empirical below. Also returns the pool size (min attainable empirical p is
    1/(1 + n_pool)) and the pool maximum."""
    from scipy.stats import norm
    shift = np.array([r["shift"] for r in rows], dtype=float)
    depth = np.array([r["depth"] for r in rows], dtype=int)
    size = np.floor(np.log2(np.maximum(1, np.array([r["n_d"] for r in rows])))).astype(int)
    gate = np.array([r["residue_pass"] for r in rows], dtype=bool)
    have = ~np.isnan(shift)
    pool = have & ~gate
    bins = {}
    for d_, s_, x in zip(depth[pool], size[pool], shift[pool]):
        bins.setdefault((int(d_), int(s_)), []).append(x)
    bins = {k: np.sort(np.array(v)) for k, v in bins.items()}
    allpool = np.sort(shift[pool])
    n = len(rows)
    out = {k: np.full(n, np.nan) for k in ("p_emp", "p_normal", "p_gpd", "n_pool", "pool_max", "gpd_xi")}
    cache = {}
    for i in np.flatnonzero(have & gate):
        key = (int(depth[i]), int(size[i]))
        if key not in cache:
            for w in (0, 1, 2, 99):
                arrs = [allpool] if w == 99 else [v for (d_, s_), v in bins.items()
                                                   if abs(d_ - key[0]) <= w and abs(s_ - key[1]) <= w]
                if sum(len(v) for v in arrs) >= min_null or w == 99:
                    break
            P = np.sort(np.concatenate(arrs)) if arrs else np.zeros(0)
            u = np.quantile(P, q_u) if len(P) else np.nan
            y = P[P > u] - u if len(P) else np.zeros(0)
            xi, sc = _gpd_pwm(y) if len(y) >= min_exc else (np.nan, np.nan)
            cache[key] = (P, P.mean() if len(P) else np.nan, P.std() if len(P) else np.nan, u, xi, sc,
                          len(y) / max(1, len(P)))
        P, mu, sd, u, xi, sc, zeta = cache[key]
        x = shift[i]
        nn = len(P)
        pe = (1 + nn - int(np.searchsorted(P, x, side="left"))) / (1 + nn)
        pn = float(norm.sf((x - mu) / sd)) if sd > 0 else (0.0 if x > mu else 1.0)
        beyond = nn > 0 and x > P[-1]
        if np.isfinite(xi) and sc > 0:
            yy = max(x - u, 0.0)
            if abs(xi) < 1e-6:
                pg = zeta * np.exp(-yy / sc)
            else:
                t = 1 + xi * yy / sc
                pg = zeta * t ** (-1 / xi) if t > 0 else 0.0
            pgpd = pg if beyond else (max(pe, pg) if x > u else pe)
        else:
            pgpd = pn if beyond else max(pe, pn)
        out["p_emp"][i] = pe
        out["p_normal"][i] = pn if beyond else max(pe, pn)
        out["p_gpd"][i] = pgpd
        out["n_pool"][i] = nn
        out["pool_max"][i] = P[-1] if nn else np.nan
        out["gpd_xi"][i] = xi
    return out


def call_specificity(rows: list, alpha: float = 0.05, min_null: int = 20):
    """Tested = residue-first gate passed (>= 2 clade-consistent BLOSUM62<0 substitutions) and a shift statistic
    exists; empirical-null p; BH across all tested SOG x branch pairs; call = BH-significant & shift above the
    technical noise floor. Returns (rows with p/q/called, number tested)."""
    p = empirical_null_pvalues(rows, min_null)
    tested = ~np.isnan(p)
    q = np.full(len(rows), np.nan)
    if tested.any():
        qq, _ = bh(p[tested], alpha)
        q[tested] = qq
    out = []
    for r, pi, qi in zip(rows, p, q):
        r = dict(r, tested=bool(not math.isnan(pi)), p_empirical=float(pi), q_bh=float(qi))
        r["called"] = bool(r["residue_pass"] and r["shift_pass"] and not math.isnan(qi) and qi <= alpha)
        out.append(r)
    return out, int(tested.sum())


def permutation_pvalue(F: np.ndarray, isd: np.ndarray, n: int, rng: np.random.Generator) -> float:
    """Secondary null (prereg specificity_change.permutation): permute D/S labels among the pair's members."""
    obs = shift_stat(F[isd], F[~isd])
    k = 0
    for _ in range(n):
        perm = rng.permutation(isd)
        if shift_stat(F[perm], F[~perm]) >= obs:
            k += 1
    return (1 + k) / (1 + n)


# --------------------------------------------------------------------------------------------- SOG tables
def load_sog_states(sog_dir, states_glob: str | None = None, state_col: str = "state") -> "pd.DataFrame":
    """Read `<sog_dir>/sog_states/bucket=*.parquet` (webapp/SCHEMA.md sogs-web-v1) or, with `states_glob`, the bulk
    state shards (e.g. /tmp/pe_big/sogs/<snap>/O1/states_t0.2/shard=*.parquet). `state_col` selects the stratum
    (`state` primary pLDDT>=90 A; `state_s70` secondary) and is returned as column `state`."""
    import glob  # noqa: PLC0415
    import pyarrow.dataset as pads  # noqa: PLC0415
    files = sorted(glob.glob(states_glob)) if states_glob else sorted(glob.glob(f"{sog_dir}/sog_states/bucket=*.parquet"))
    if not files:
        raise FileNotFoundError(f"no state shards ({states_glob or sog_dir})")
    cols = ["og_id", "sog_id", "taxid", "accession", state_col, "footprint_plddt", "reason"]
    dset = pads.dataset(files, format="parquet")
    have = [c for c in cols if c in dset.schema.names]
    df = dset.to_table(columns=have).to_pandas()
    for c in cols:
        if c not in df.columns:
            df[c] = None
    if state_col != "state":
        df = df.drop(columns=[c for c in ["state"] if c in df.columns]).rename(columns={state_col: "state"})
    return df[["og_id", "sog_id", "taxid", "accession", "state", "footprint_plddt", "reason"]]


def eligible_sogs(states, min_P: int = 2, min_A: int = 2):
    """SOG ids that can possibly yield a gain or loss candidate under the frozen rules (gain: >= 3 P species and
    >= 2 A members; loss: >= 2 A and P in >= 2 lineages) -> necessary condition n_P >= min_P and n_A >= min_A."""
    c = states.groupby("sog_id")["state"].agg(nP=lambda x: int((x == "P").sum()), nA=lambda x: int((x == "A").sum()))
    return set(c.index[(c["nP"] >= min_P) & (c["nA"] >= min_A)])


def member_confidence(state: str, footprint_plddt, reason) -> float:
    """prereg candidate_score.state_confidence: P 1.0, W 0.5, A 1.0 if footprint pLDDT >= 90 and rescue negative
    in 6/6 settings else 0.5; U 0."""
    import re as _re  # noqa: PLC0415
    if state == "P":
        return 1.0
    if state == "W":
        return 0.5
    if state == "A":
        try:
            ok = float(footprint_plddt) >= 90
        except (TypeError, ValueError):
            ok = False
        if reason and _re.search(r"rescue:[1-9]", str(reason)):
            ok = False
        return 1.0 if ok else 0.5
    return 0.0


def assemble(trees: dict, states, counters: dict | None = None):
    """Group SOG member states by root OG and map them onto the OG's gene tree.

    Returns list of (og_id, sog_ids, codes (n_leaves, n_sogs) int8, confidence (n_leaves, n_sogs)).
    Members of the tree without a row in a SOG -> U (counted 'member_without_row'); rows whose accession is
    not a leaf of the tree -> dropped ('accession_not_in_tree'); OG not in the orthology -> ('og_not_in_orthology');
    G rows are not leaves ('G_rows', informational)."""
    cnt = counters if counters is not None else {}
    for k in ("og_not_in_orthology", "accession_not_in_tree", "member_without_row", "G_rows", "taxid_mismatch",
              "sogs_in", "rows_in"):
        cnt.setdefault(k, 0)
    out = []
    cnt["rows_in"] += len(states)
    for og, g in states.groupby("og_id", sort=True):
        sog_ids = sorted(g["sog_id"].unique())
        cnt["sogs_in"] += len(sog_ids)
        t = trees.get(og)
        if t is None:
            cnt["og_not_in_orthology"] += len(sog_ids)
            continue
        names = t.leaf_names()
        pos = {a: i for i, a in enumerate(names)}
        txs = t.leaf_taxids()
        col = {s: j for j, s in enumerate(sog_ids)}
        codes = np.full((t.n_leaves, len(sog_ids)), U_, dtype=np.int8)
        conf = np.zeros((t.n_leaves, len(sog_ids)))
        seen = np.zeros((t.n_leaves, len(sog_ids)), dtype=bool)
        for sog, acc, st, tx, pl, rs in zip(g["sog_id"], g["accession"], g["state"], g["taxid"],
                                            g["footprint_plddt"], g["reason"]):
            if st == "G" or acc is None:
                cnt["G_rows"] += 1
                continue
            i = pos.get(acc)
            if i is None:
                cnt["accession_not_in_tree"] += 1
                continue
            if tx is not None and int(tx) != int(txs[i]) and txs[i] >= 0:
                cnt["taxid_mismatch"] += 1
            j = col[sog]
            codes[i, j] = STATE_CODES.get(st, U_)
            conf[i, j] = member_confidence(st, pl, rs)
            seen[i, j] = True
        cnt["member_without_row"] += int((~seen).sum())
        out.append((og, sog_ids, codes, conf))
    return out


def call_og(t: Tree, sog_ids: list, codes: np.ndarray, conf: np.ndarray, rules: Rules, costs=(2.0, 1.0),
            w_mode: str = "present", species_tree=None, tix: TreeIndex | None = None, keep_states: bool = False):
    """Reconstruct all SOGs of one OG and return (events with rule flags, node-state rows or None)."""
    res = recon.sankoff(t, codes, costs[0], costs[1], w_mode)
    tix = tix or TreeIndex(t, species_tree)
    ev = [apply_rules(e, rules) for e in extract_events(t, codes, res, rules, tix, confidence=conf,
                                                        char_ids=sog_ids)]
    for e in ev:
        e["og_id"] = t.tree_id
        e["gain_cost"], e["loss_cost"] = costs
    states = None
    if keep_states:
        un = res.unambiguous
        states = []
        for j, s in enumerate(sog_ids):
            for v in range(t.n_nodes):
                states.append((t.tree_id, s, v, t.node_name[v], recon.KIND_NAME[int(t.kind[v])], t.level[v],
                               int(t.parent[v]), int(res.mpr[v, j]), bool(un[v, j]),
                               float(min(res.total[v, j, 0], recon.BIG)), float(min(res.total[v, j, 1], recon.BIG))))
    return ev, states, res


_RUN = {}


def _run_one(args):
    og, sog_ids, codes, conf = args
    t = _RUN["trees"][og]
    out, states = [], []
    for costs, rules in _RUN["settings"]:
        keep = _RUN["keep_states"] and costs == _RUN["settings"][0][0]
        ev, stt, _ = call_og(t, sog_ids, codes, conf, rules, costs, _RUN["w_mode"], _RUN["species_tree"],
                             keep_states=keep)
        out.extend(ev)
        if stt:
            states.extend(stt)
    return out, states


def run_assembled(assembled: list, trees: dict, settings: list, species_tree=None, w_mode: str = "present",
                  workers: int = 1, keep_states: bool = False):
    """Reconstruct + call every SOG of `assembled` under each (costs, Rules) in `settings` (first = primary).
    Returns (events list, node-state rows of the primary setting)."""
    _RUN.update(trees=trees, settings=settings, species_tree=species_tree, w_mode=w_mode, keep_states=keep_states)
    ev, st = [], []
    if workers <= 1:
        for x in assembled:
            e, s = _run_one(x)
            ev.extend(e); st.extend(s)
    else:
        from multiprocessing import get_context  # noqa: PLC0415
        with get_context("fork").Pool(workers) as pool:
            for e, s in pool.imap(_run_one, assembled, chunksize=8):
                ev.extend(e); st.extend(s)
    return ev, st


NODE_STATE_COLUMNS = ["og_id", "sog_id", "node", "node_name", "kind", "level", "parent", "state_mpr",
                      "unambiguous", "cost_absent", "cost_present"]


def rules_by_cost(thresholds: dict | None, cost_grid, pr: dict | None = None) -> list:
    """[(costs, Rules)] with calibrated cuts from results/partB/sim/thresholds.json (chosen[<g>:<l>])."""
    out = []
    for g, l in cost_grid:
        cal = {}
        ch = (thresholds or {}).get("chosen", {}).get(f"{int(g)}:{int(l)}", {})
        for typ in ("gain", "loss"):
            c = ch.get(typ)
            if c:
                cal[f"{typ}_min_species_cal"] = c["min_species"]
                cal[f"{typ}_min_margin"] = c["min_margin"]
                cal[f"{typ}_min_out_lineages_cal"] = c.get("min_out_lineages", 2)
        out.append(((float(g), float(l)), Rules.from_prereg(pr, cal)))
    return out
