"""Tests for the Part B engine (src/pocketevo/recon.py, events.py, simulate.py).

Hand-computed Sankoff results (gain cost g, loss cost l; C[parent, child]); trees written as nested tuples.
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pocketevo import events, recon, simulate  # noqa: E402
from pocketevo.recon import encode_states, sankoff, tree_from_nested  # noqa: E402


def T(nested, species=None):
    species = species or {}
    return tree_from_nested("toy", nested, species_of=species)


def leaf_states(tree, states: dict):
    return encode_states([states[a] for a in tree.leaf_names()])


# ------------------------------------------------------------------------------------------- hand-computed
def test_simple_loss():
    t = T(("S", "root", "r", [("S", "AB", "ab", ["a", "b"]), "c"]))
    res = sankoff(t, leaf_states(t, {"a": "P", "b": "P", "c": "A"}), 2, 1)
    assert res.opt[0] == 1.0                        # root present + loss on c (gain alternative costs 2)
    assert res.mpr[:, 0].tolist() == [1, 1, 1, 1, 0]
    assert res.unambiguous[:, 0].all()
    # conditional optimal costs at the root: present 1, absent 2 (gain on ab stem)
    assert res.total[0, 0].tolist() == [2.0, 1.0]


def test_cost_ratio_changes_root():
    # polytomy: 2 P and 3 A children of the root
    t = T(("S", "root", "r", ["a", "b", "c", "d", "e"]))
    x = leaf_states(t, dict(a="P", b="P", c="A", d="A", e="A"))
    r21 = sankoff(t, x, 2, 1)
    assert r21.opt[0] == 3.0 and r21.mpr[0, 0] == 1     # 3 losses (3) < 2 gains (4)
    r11 = sankoff(t, x, 1, 1)
    assert r11.opt[0] == 2.0 and r11.mpr[0, 0] == 0     # 2 gains (2) < 3 losses (3)
    r31 = sankoff(t, x, 3, 1)
    assert r31.opt[0] == 3.0 and r31.mpr[0, 0] == 1


def test_tie_is_ambiguous():
    # (a:P, (b:A, c:A)) with 1:1 -> root absent + gain at a (1) == root present + loss on bc (1)
    t = T(("S", "root", "r", ["a", ("S", "BC", "bc", ["b", "c"])]))
    x = leaf_states(t, dict(a="P", b="A", c="A"))
    r = sankoff(t, x, 1, 1)
    assert r.opt[0] == 1.0
    assert not r.unambiguous[0, 0]
    r2 = sankoff(t, x, 2, 1)                       # 2:1 breaks the tie: loss on bc
    assert r2.opt[0] == 1.0 and r2.unambiguous[:, 0].all() and r2.mpr[2, 0] == 0


def test_duplication_loss_on_one_copy():
    # D(S(h1:P, m1:P), S(h2:A, m2:A)): one loss on the second paralog copy
    t = T(("D", "Euteleostomi", "d", [("S", "Murinae", "s1", ["h1", "m1"]), ("S", "Murinae", "s2", ["h2", "m2"])]))
    r = sankoff(t, leaf_states(t, dict(h1="P", m1="P", h2="A", m2="A")), 2, 1)
    assert r.opt[0] == 1.0
    assert r.mpr[:, 0].tolist() == [1, 1, 1, 1, 0, 0, 0]
    ev = events.extract_events(t, leaf_states(t, dict(h1="P", m1="P", h2="A", m2="A")), r, events.Rules(),
                               prefilter=False)
    assert len(ev) == 1 and ev[0]["type"] == "loss" and ev[0]["node_name"] == "s2"
    assert ev[0]["parent_kind"] == "duplication" and ev[0]["unambiguous"]


def test_missing_and_weak():
    t = T(("S", "root", "r", ["a", "b", "c"]))
    r = sankoff(t, leaf_states(t, dict(a="P", b="U", c="P")), 2, 1)
    assert r.opt[0] == 0.0 and (r.mpr[:, 0] == 1).all()
    # W in "present" mode: (a:W, b:A, c:A): gain at a (2) == losses at b, c (2) -> ambiguous root
    x = leaf_states(t, dict(a="W", b="A", c="A"))
    r = sankoff(t, x, 2, 1, w_mode="present")
    assert r.opt[0] == 2.0 and not r.unambiguous[0, 0]
    r = sankoff(t, x, 2, 1, w_mode="missing")
    assert r.opt[0] == 0.0
    with pytest.raises(ValueError):
        recon.leaf_costs(np.array([[recon.STATE_CODES["G"]]]))


def _brute(tree, codes, g, l):
    """Enumerate all internal assignments: min total cost and per-node conditional minima."""
    C = recon.cost_matrix(g, l)
    internal = [v for v in range(tree.n_nodes) if tree.kind[v] != recon.LEAF]
    lc = recon.leaf_costs(codes[:, None])[:, 0, :]
    best = np.full((tree.n_nodes, 2), np.inf)
    opt = np.inf
    for assign in itertools.product([0, 1], repeat=len(internal)):
        st = np.zeros(tree.n_nodes, int)
        st[internal] = assign
        base = 0.0
        for lv in tree.leaves:                            # leaves free given costs: choose best state
            pass
        # leaves: enumerate is too large; leaves choose best state given parent
        cost = 0.0
        leaf_state = {}
        for v in range(1, tree.n_nodes):
            p = tree.parent[v]
            if tree.kind[v] == recon.LEAF:
                k = tree.leaf_rank[v]
                opts = [lc[k, s] + C[st[p], s] for s in (0, 1)]
                s = int(np.argmin(opts))
                leaf_state[v] = s
                cost += opts[s]
            else:
                cost += C[st[p], st[v]]
        _ = base
        opt = min(opt, cost)
        for v in internal:
            best[v, st[v]] = min(best[v, st[v]], cost)
    return opt, best, internal


def test_against_brute_force():
    rng = np.random.default_rng(0)
    for trial in range(40):
        # random tree with <= 7 internal nodes
        nested = _random_nested(rng, depth=3)
        t = T(nested)
        codes = rng.choice([0, 1, 2, 3], size=t.n_leaves).astype(np.int8)
        for g, l in [(2, 1), (1, 1), (3, 1)]:
            r = sankoff(t, codes, g, l)
            opt, best, internal = _brute(t, codes, g, l)
            assert np.isclose(r.opt[0], opt)
            for v in internal:
                assert np.allclose(np.minimum(r.total[v, 0], 1e8), np.minimum(best[v], 1e8)), (trial, v)


def _random_nested(rng, depth, counter=[0]):
    k = int(rng.integers(2, 4))
    ch = []
    for _ in range(k):
        if depth > 0 and rng.random() < 0.5:
            ch.append(_random_nested(rng, depth - 1))
        else:
            counter[0] += 1
            ch.append(f"x{counter[0]}")
    counter[0] += 1
    return ("S" if rng.random() < 0.7 else "D", "lev", f"n{counter[0]}", ch)


def test_orthoxml_roundtrip(tmp_path):
    xml = """<?xml version="1.0"?><orthoXML xmlns="http://orthoXML.org/2011/" version="0.5">
<species name="UP1_9606"><database name="d"><genes><gene id="1" protId="H1"/><gene id="2" protId="H2"/></genes></database></species>
<species name="UP2_10090"><database name="d"><genes><gene id="3" protId="M1"/></genes></database></species>
<species name="UP3_7227"><database name="d"><genes><gene id="4" protId="F1"/></genes></database></species>
<taxonomy><taxon id="1" name="Bilateria"><taxon id="2" name="Euarchontoglires"><taxon id="3" name="UP1_9606"/><taxon id="4" name="UP2_10090"/></taxon><taxon id="5" name="UP3_7227"/></taxon></taxonomy>
<groups><orthologGroup id="HOG1" taxonId="1"><property name="TaxRange" value="Bilateria"/>
<paralogGroup id="HOG1.1" taxonId="2"><orthologGroup id="HOG1.1a" taxonId="2"><property name="TaxRange" value="Euarchontoglires"/><geneRef id="1"/><geneRef id="3"/></orthologGroup><geneRef id="2"/></paralogGroup>
<geneRef id="4"/></orthologGroup></groups></orthoXML>"""
    p = tmp_path / "h.orthoxml"
    p.write_text(xml)
    trees = list(recon.iter_orthoxml_trees(p))
    assert len(trees) == 1
    t = trees[0]
    assert t.tree_id == "HOG1" and t.leaf_names() == ["H1", "M1", "H2", "F1"]
    assert [recon.KIND_NAME[k] for k in t.kind] == ["speciation", "duplication", "speciation", "leaf", "leaf",
                                                     "leaf", "leaf"]
    assert t.level[1] == "Euarchontoglires" and t.leaf_taxids().tolist() == [9606, 10090, 9606, 7227]


def test_species_tree():
    st = recon.SpeciesTree("((UP1_9606,UP2_10090)Euarch,UP3_7227,UP4_6239)'Bila teria';")
    assert st.names[0] == "Bila teria"
    assert st.tree.n_leaves == 4 and len(st.tree.children(0)) == 3      # polytomy kept
    v = st.lca_of_taxids([9606, 10090])
    assert st.names[v] == "Euarch" and st.level_depth("Euarch") == 1


# ------------------------------------------------------------------------------------------- planted events
def _balanced(n_species=12):
    sp = [f"UP{i}_{1000 + i}" for i in range(n_species)]

    def build(lo, hi, k=[0]):
        if hi - lo == 1:
            return ("leaf", sp[lo], f"P{lo}")
        mid = (lo + hi) // 2
        k[0] += 1
        return ("S", f"L{lo}_{hi}", f"n{lo}_{hi}", [build(lo, mid), build(mid, hi)])

    return tree_from_nested("bal", build(0, n_species)), sp


@pytest.mark.parametrize("scenario", ["gain", "loss"])
def test_planted_single_event_recovered(scenario):
    # 32 species: a gain is uniquely parsimonious at 2:1 only if the alternative needs >= 3 losses (clade nested
    # >= 3 splits deep); top-level "gains" are reconstructed as losses of the sister (no outgroup) -> not called
    t, _ = _balanced(32)
    tix = events.TreeIndex(t)
    rules = events.Rules()
    br = simulate.eligible_branches(t, tix, scenario, min_species=3)
    assert len(br) > 0
    n_rec = 0
    for c in br:
        truth, _ = simulate.simulate_truth(t, 1, scenario, np.random.default_rng(0), branches=np.array([c]))
        codes = np.where(truth == 1, recon.STATE_CODES["P"], recon.STATE_CODES["A"]).astype(np.int8)
        r = sankoff(t, codes, 2, 1)
        ev = [events.apply_rules(e, rules) for e in events.extract_events(t, codes, r, rules, tix)]
        called = [e for e in ev if e["called"]]
        assert all(e["node"] == c and e["type"] == scenario for e in called), (scenario, c, called)
        if r.unambiguous[:, 0].all():          # unique MPR -> the planted event must be called
            assert len(called) == 1, (scenario, c, ev)
            n_rec += 1
        else:                                  # tie (e.g. nested gain = 2 losses at 2:1) -> nothing called
            assert scenario == "gain" and not called
        if scenario == "gain" and t.parent[c] == 0:
            assert not called
    assert n_rec >= len(br) // 2


def test_gain_rules_reject_small_or_unsupported():
    t, _ = _balanced(8)
    tix = events.TreeIndex(t)
    rules = events.Rules()
    # gain in a 2-species clade (n4_6) -> fails the 3-species rule
    c = t.node_name.index("n4_6")
    truth, _ = simulate.simulate_truth(t, 1, "gain", np.random.default_rng(0), branches=np.array([c]))
    codes = np.where(truth == 1, 0, 2).astype(np.int8)
    r = sankoff(t, codes, 2, 1)
    ev = [events.apply_rules(e, rules) for e in events.extract_events(t, codes, r, rules, tix, prefilter=False)]
    assert not any(e["called"] for e in ev)
    # gain in a 4-species clade but all outgroup members U -> no absent outgroup -> no event at all
    c = t.node_name.index("n4_8")
    truth, _ = simulate.simulate_truth(t, 1, "gain", np.random.default_rng(0), branches=np.array([c]))
    codes = np.where(truth == 1, 0, 3).astype(np.int8)
    r = sankoff(t, codes, 2, 1)
    ev = [events.apply_rules(e, rules) for e in events.extract_events(t, codes, r, rules, tix)]
    assert not any(e["called"] for e in ev)


def test_no_event_world_false_calls_bounded():
    """Vertical inheritance + default noise: called events per SOG stay small; zero noise -> zero events."""
    t, sp = _balanced(16)
    tix = events.TreeIndex(t)
    rules = events.Rules()
    rng = np.random.default_rng(1)
    truth, _ = simulate.simulate_truth(t, 2000, "vertical", rng)
    clean = np.where(truth == 1, 0, 2).astype(np.int8)
    r = sankoff(t, clean, 2, 1)
    assert r.opt.max() == 0 and not events.extract_events(t, clean, r, rules, tix)
    nm = simulate.NoiseModel()
    codes = simulate.observe(truth, nm.rates(t.leaf_species()), rng)
    r = sankoff(t, codes, 2, 1)
    ev = [events.apply_rules(e, rules) for e in events.extract_events(t, codes, r, rules, tix)]
    per_sog = sum(e["called"] for e in ev) / 2000
    # a loss call needs >= 2 false-A species in one clade with >= 80 % of informative members A: with
    # false_A = 0.05 the expected rate is << 0.05 per SOG (pairs of sister leaves: ~ 8 * 0.05^2 * ...)
    assert per_sog < 0.02, per_sog


# ------------------------------------------------------------------------------------------- M3 / nulls
def test_dollo_and_curveball():
    t, _ = _balanced(8)
    X = np.array([[1] * 8,                        # 0 events
                  [0, 0, 0, 0, 1, 1, 1, 1],       # one gain (LCA n4_8)
                  [1, 1, 0, 0, 1, 1, 1, 1],       # one loss at n2_4
                  [1, 0, 1, 0, 1, 0, 1, 0],       # many
                  [-1] * 8], dtype=np.int8)
    ev, lca = recon.dollo_events(t, X)
    assert ev.tolist() == [0, 1, 1, 4, -1]
    rng = np.random.default_rng(0)
    M = rng.integers(0, 2, size=(40, 10)).astype(np.int8)
    M[rng.random(M.shape) < 0.2] = -1
    Y, done, _ = simulate.curveball(M, seed=3)
    assert done >= 100 * (M == 1).sum()
    assert (Y.sum(axis=1) == M.sum(axis=1)).all() and ((Y == 1).sum(0) == (M == 1).sum(0)).all()
    assert ((Y < 0) == (M < 0)).all() and (Y != M).any()


def test_replicate_pair_null_no_events():
    t, _ = _balanced(8)
    tt = simulate.replicate_pair_tree(t)
    assert tt.n_leaves == 2 * t.n_leaves and tt.n_nodes == 2 * t.n_nodes + 1
    codes = np.zeros(tt.n_leaves, np.int8)        # same proteins, same states -> nothing
    r = sankoff(tt, codes, 2, 1)
    assert r.opt[0] == 0


def test_specificity_toy():
    t, _ = _balanced(8)
    codes = np.zeros(t.n_leaves, np.int8)
    names = t.leaf_names()
    # clade n0_4 has 'W' 'D' at two footprint columns where the sister has 'F' 'K' (BLOSUM62 W/F=1 >= 0, D/K=-1)
    res = {a: ("WDE" if i < 4 else "FKE") for i, a in enumerate(names)}
    res2 = {a: ("GDE" if i < 4 else "FKE") for i, a in enumerate(names)}   # G/F=-3, D/K=-1 -> 2 substitutions
    rng = np.random.default_rng(0)
    fp = {a: rng.normal(size=16) + (3.0 if i < 4 else 0.0) * np.eye(16)[0] for i, a in enumerate(names)}
    rows = events.specificity_tests(t, codes, res, fp, None, floor=0.05)
    root = [r for r in rows if r["node"] == 0][0]
    assert root["n_subst"] == 1 and not root["residue_pass"]
    rows2 = events.specificity_tests(t, codes, res2, fp, None, floor=0.05)
    root2 = [r for r in rows2 if r["node"] == 0][0]
    assert root2["n_subst"] == 2 and root2["residue_pass"] and root2["shift_pass"]
    q, rej = events.bh(np.array([0.01, 0.04, 0.03, 0.5]), 0.05)
    assert rej.tolist() == [True, False, False, False]
    assert np.allclose(q, [0.04, 0.0533333, 0.0533333, 0.5], atol=1e-6)


def test_rules_from_prereg_frozen_values():
    r = events.Rules.from_prereg()
    assert (r.gain_min_species, r.gain_min_frac_present, r.gain_min_absent_members_out, r.gain_min_out_lineages) \
        == (3, 0.8, 2, 2)
    assert (r.loss_min_out_lineages_present, r.loss_min_species_absent, r.loss_min_frac_absent) == (2, 2, 0.8)
    with pytest.raises(ValueError):
        r.with_cal(gain_min_species_cal=2)
    c = recon.costs_from_prereg()
    assert c["primary"] == (2.0, 1.0) and (1, 1) in c["grid"] and (3, 1) in c["grid"]


def test_checkerboard_preserves_margins_and_mask():
    rng = np.random.default_rng(1)
    M = (rng.random((300, 12)) < 0.85).astype(np.int8)
    M[rng.random(M.shape) < 0.5] = -1
    Y, done, _ = simulate.checkerboard(M, n_swaps=5000, seed=2)
    assert done == 5000
    assert (Y.sum(axis=1) == M.sum(axis=1)).all() and ((Y == 1).sum(0) == (M == 1).sum(0)).all()
    assert ((Y < 0) == (M < 0)).all() and (Y != M).any()


def test_false_w_on_true_absences():
    nm = simulate.NoiseModel(u_absent=0.0, false_P=0.0, w_absent=0.35)
    truth = np.zeros((2000, 5), np.int8)
    obs = simulate.observe(truth, nm.rates(["s"] * 2000), np.random.default_rng(0))
    frac_w = (obs == recon.STATE_CODES["W"]).mean()
    assert abs(frac_w - 0.35) < 0.02 and not (obs == recon.STATE_CODES["P"]).any()
