"""Toy worlds with known answers for Part A (DESIGN §8.3: synthetic pair sets with known class composition).

Run:  /tmp/pe_big/envs/faiss/bin/python -m pytest -q tests/test_partA.py
(needs numpy/pandas/scipy; the GPD and CovModel block-form tests use torch on CPU and are skipped with a message
when torch is absent, e.g. under `uv run --project . pytest` - torch lives in /tmp/pe_big/envs/faiss, docs/ENV.md.)
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pocketevo import estimands as est  # noqa: E402
from pocketevo import homology as hom  # noqa: E402
from pocketevo import nulls  # noqa: E402

ORTHOXML = """<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.5">
 <species name="A" taxonId="1"><database name="x" version="1"><genes>
  <gene id="1" protId="A1"/><gene id="2" protId="A2"/><gene id="7" protId="A3"/></genes></database></species>
 <species name="B" taxonId="2"><database name="x" version="1"><genes>
  <gene id="3" protId="B1"/><gene id="4" protId="B2"/><gene id="8" protId="B3"/></genes></database></species>
 <species name="C" taxonId="3"><database name="x" version="1"><genes>
  <gene id="5" protId="C1"/><gene id="6" protId="C2"/><gene id="9" protId="C9"/></genes></database></species>
 <groups>
  <orthologGroup id="R1">
   <paralogGroup id="R1.1">
     <orthologGroup id="R1.1a"><geneRef id="1"/><geneRef id="3"/></orthologGroup>
     <orthologGroup id="R1.1b"><geneRef id="2"/><geneRef id="4"/></orthologGroup>
   </paralogGroup>
   <geneRef id="5"/>
  </orthologGroup>
  <paralogGroup id="R2"><geneRef id="6"/><orthologGroup id="R2a"><geneRef id="7"/><geneRef id="8"/></orthologGroup>
  </paralogGroup>
 </groups>
</orthoXML>
"""


@pytest.fixture
def forest(tmp_path):
    p = tmp_path / "h.orthoxml"
    p.write_text(ORTHOXML)
    return hom.parse_orthoxml(str(p))


def test_orthoxml_lca_and_exact_counts(forest):
    f = forest
    acc = list(f.gene_acc)
    assert acc == ["A1", "B1", "A2", "B2", "C1", "C2", "A3", "B3"]   # C9 is in no group
    g = {a: k for k, a in enumerate(acc)}
    # A1-B1 ortholog (LCA R1.1a), A1-A2 paralog (R1.1), A1-C1 ortholog (R1), C2-A3 paralog (R2), A3-B3 ortholog
    pairs = [("A1", "B1", 1), ("A1", "A2", 2), ("A1", "B2", 2), ("A1", "C1", 1), ("C2", "A3", 2), ("A3", "B3", 1)]
    a = np.array([g[x] for x, _, _ in pairs]); b = np.array([g[y] for _, y, _ in pairs])
    assert list(hom.lca_type(f, a, b)) == [t for *_, t in pairs]
    assert list(hom.lca_type(f, b, a)) == [t for *_, t in pairs]
    w = np.array([2, 1, 3, 1, 2, 1, 1, 4], float)
    ex = hom.exact_ortho_para_counts(f, w)
    # brute force
    bf = {0: [0.0, 0.0], 1: [0.0, 0.0]}
    for x, y in itertools.combinations(range(len(acc)), 2):
        rx, ry = f.node_root_rank[f.gene_parent[x]], f.node_root_rank[f.gene_parent[y]]
        if rx != ry:
            continue
        t = hom.lca_type(f, np.array([x]), np.array([y]))[0]
        bf[rx][t - 1] += w[x] * w[y]
    for r in (0, 1):
        assert ex.loc[r, "ortho"] == pytest.approx(bf[r][0])
        assert ex.loc[r, "para"] == pytest.approx(bf[r][1])


def random_world(rng, n_prot=120, n_groups=25, n_lab=15):
    """Proteins in groups (some groups proxy-only), 0-4 sites each, random region labels."""
    prot_group = rng.integers(0, n_groups, n_prot)
    prot_gene = np.full(n_prot, -1)
    site_prot = np.repeat(np.arange(n_prot), rng.integers(0, 5, n_prot))
    N = len(site_prot)
    lab = np.full((N, hom.K_SITE + 1), -1, np.int32)
    for s in range(N):
        k = rng.integers(0, 4)
        codes = np.sort(rng.choice(n_lab, k, replace=False))
        lab[s, :k] = codes
    return hom.PairContext(site_prot, lab, prot_group, prot_gene, None)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_exact_totals_match_bruteforce(seed):
    rng = np.random.default_rng(seed)
    ctx = random_world(rng)
    N = len(ctx.site_prot)
    i, j = np.triu_indices(N, 1)
    cls = hom.classify_pairs(ctx, i, j)
    gw = rng.poisson(1.0, ctx.prot_group.max() + 1).astype(float)
    gi_, gj_ = ctx.prot_group[ctx.site_prot[i]], ctx.prot_group[ctx.site_prot[j]]
    pw = est.pair_group_weight(gw, gi_, gj_)
    for w, pwt in ((None, np.ones(len(i))), (gw, pw)):
        tot = hom.exact_class_totals(ctx, group_weight=w)
        assert tot["H12"] == pytest.approx(pwt[cls == hom.H12].sum())   # no forest: all within-group -> H12
        assert tot["H3j"] == pytest.approx(pwt[cls == hom.H3J].sum())
        assert tot["rest"] == pytest.approx(pwt[cls == hom.REST].sum())
        assert tot["total"] == pytest.approx(pwt[cls != hom.H0].sum())


def test_site_mask_totals():
    rng = np.random.default_rng(5)
    ctx = random_world(rng)
    N = len(ctx.site_prot)
    mask = rng.random(N) < 0.6
    idx = np.flatnonzero(mask)
    a, b = np.triu_indices(len(idx), 1)
    i, j = idx[a], idx[b]
    cls = hom.classify_pairs(ctx, i, j)
    tot = hom.exact_class_totals(ctx, site_mask=mask)
    assert tot["H3j"] == (cls == hom.H3J).sum()
    assert tot["rest"] == (cls == hom.REST).sum()


def test_region_labels():
    regions = pd.DataFrame({"accession": ["X", "X", "X"], "rid": [0, 0, 1], "label": ["a", "a", "b"],
                            "start": [1, 50, 20], "end": [10, 60, 30]})
    labs, n = hom.region_labels_for_sites(np.array(["X", "X", "Y"]),
                                          [np.array([2, 3, 55, 25]), np.array([21, 22, 23]), np.array([1])],
                                          regions, 0.5, "P:")
    assert labs == [["P:a"], ["P:b"], []] and n == 2


def _gpd_case(kind, rng, Q, B, M):
    if kind == "gamma":        # exponential-type tail, true xi = 0
        shape = rng.uniform(2, 8, Q)[:, None]
        return rng.gamma(shape, 1.0, (Q, B)), rng.gamma(shape, 1.0, (Q, M))
    if kind == "beta":         # bounded tail, true xi < 0 (like a cosine)
        a = rng.uniform(2, 5, Q)[:, None]
        return rng.beta(a, 8, (Q, B)), rng.beta(a, 8, (Q, M))
    raise ValueError(kind)


@pytest.mark.parametrize("kind", ["gamma", "beta"])
def test_gpd_tail_calibration(kind):
    """Per-query GPD on 10^4 null draws, extrapolated 5,000x below the threshold quantile (1e-5 vs 0.05).

    Known answer: the per-query PWM shape is noisy and its negative errors make extreme tails too light
    (observed >> expected; documented, this is why 'perq' is not primary). The pooled-shape tail is calibrated
    within a factor 2 for the exponential-type world, and the xi>=0 variant is never anti-conservative.
    """
    torch = pytest.importorskip("torch", reason="torch not installed in this env (GPU path optional; docs/ENV.md)")
    rng = np.random.default_rng(0)
    Q, B, M = 200, 10_000, 200_000
    Zb, fresh = _gpd_case(kind, rng, Q, B, M)
    fit = nulls.gpd_pwm_rows(torch.tensor(Zb), torch.ones(Q, B, dtype=torch.bool))
    xi_pool = float(torch.median(fit["xi"]))
    fresh = torch.tensor(fresh)
    res = {}
    for var in ("perq", "pooled", "xi0"):
        sc, xi = nulls.tail_variant(fit, var, xi_pool)
        p = nulls.gpd_sf(fresh, fit["u"][:, None], sc[:, None], xi[:, None], fit["zeta"][:, None])
        res[var] = {t: float((p < t).sum()) / (Q * M * t) for t in (1e-3, 1e-4, 1e-5)}
    for t in (1e-3, 1e-4):
        assert 0.5 < res["pooled"][t] < 2.0, (kind, res)
    assert res["xi0"][1e-5] < 1.6, (kind, res)
    if kind == "gamma":
        assert 0.5 < res["pooled"][1e-5] < 2.0, res
        assert res["perq"][1e-5] > res["pooled"][1e-5]


def test_covmodel_recovers_size_effect():
    rng = np.random.default_rng(1)
    n = 200_000
    Fi, Fj = rng.normal(size=(n, 2)), rng.normal(size=(n, 2))
    s = 0.1 + 0.05 * (Fi[:, 0] + Fj[:, 0]) + 0.03 * Fi[:, 0] * Fj[:, 0] + rng.normal(size=n) * 0.1 * np.exp(
        0.2 * (Fi[:, 1] + Fj[:, 1]))
    m = nulls.CovModel.fit(Fi, Fj, s, mean=np.zeros(2), sd=np.ones(2))
    c, a, b, W = nulls.unpack(m.beta_mu, 2)
    assert c == pytest.approx(0.1, abs=0.005) and a[0] == pytest.approx(0.05, abs=0.003)
    assert W[0, 0] == pytest.approx(0.03, abs=0.003)
    mu, ls = m.mu_logsd_pairs(Fi, Fj)
    z = m.z(s, mu, ls)
    assert abs(np.std(z) - 1) < 0.05


def test_covmodel_block_form_equals_pair_form():
    """GPU/torch path (optional: the production null runs in /tmp/pe_big/envs/faiss with torch; see docs/ENV.md)."""
    torch = pytest.importorskip("torch", reason="torch not installed in this env (GPU path optional; docs/ENV.md)")
    rng = np.random.default_rng(1)
    n = 20_000
    Fi, Fj = rng.normal(size=(n, 2)), rng.normal(size=(n, 2))
    s = 0.1 + 0.05 * (Fi[:, 0] + Fj[:, 0]) + rng.normal(size=n) * 0.1
    m = nulls.CovModel.fit(Fi, Fj, s, mean=np.zeros(2), sd=np.ones(2))
    fq, fb = torch.tensor(Fi[:5]), torch.tensor(Fj[:7])
    mub, lsb = m.mu_logsd_block(fq, fb)
    mup, lsp = m.mu_logsd_pairs(np.repeat(Fi[:5], 7, 0), np.tile(Fj[:7], (5, 1)))
    assert np.allclose(mub.numpy().ravel(), mup) and np.allclose(lsb.numpy().ravel(), lsp)


def toy_pair_world(seed=0, q0=2e-4, r_true=None, n_planted=300):
    """Known world: exact class sizes, per-class significance rates, null FP rate q0, planted analog pairs."""
    rng = np.random.default_rng(seed)
    r_true = r_true or {"H1s": 0.6, "H1x": 0.01, "H2s": 0.3, "H2x": 0.005, "H12s": 0.4, "H12x": 0.01, "H3j": 0.05,
                        "H3s": 0.02, "H4": 0.004, "N": 0.0}
    n_groups = 400
    n_true = {"H1s": 20_000, "H1x": 180_000, "H2s": 10_000, "H2x": 90_000, "H12s": 5_000, "H12x": 45_000,
              "H3j": 300_000, "H3s": 50_000, "H4": 150_000, "N": 3_000_000}
    rows, all_groups = [], {}
    for c, n in n_true.items():
        within = c[:2] in ("H1", "H2") or c.startswith("H12")
        gi = rng.integers(0, n_groups, n)
        gj = gi if within else rng.integers(0, n_groups, n)
        all_groups[c] = (gi, gj)
        sig = rng.random(n) < r_true[c] + q0            # true signal + null false positives
        if c == "N":
            planted = rng.choice(n, n_planted, replace=False)
            sig[planted] = True
        k = np.flatnonzero(sig)
        rows.append(pd.DataFrame({"i": rng.integers(0, 10**7, len(k)), "j": rng.integers(0, 10**7, len(k)),
                                  "cls": c, "gi": gi[k], "gj": gj[k]}))
        if c == "N":
            n_sig_N = len(k)
    sig = pd.concat(rows, ignore_index=True)
    totals = {"H1": n_true["H1s"] + n_true["H1x"], "H2": n_true["H2s"] + n_true["H2x"],
              "H12": n_true["H12s"] + n_true["H12x"], "H3j": n_true["H3j"],
              "rest": n_true["H3s"] + n_true["H4"] + n_true["N"]}
    split = {"H1": 0.1, "H2": 0.1, "H12": 0.1}
    # uniform sample of 'rest' pairs with oracle search labels
    M = 100_000
    rest = ["H3s", "H4", "N"]
    rest_cls = np.concatenate([[c] * n_true[c] for c in rest])
    rest_g = [np.concatenate([all_groups[c][k] for c in rest]) for k in (0, 1)]
    pick = rng.integers(0, len(rest_cls), M)
    sample = pd.DataFrame({"cls": rest_cls[pick], "gi": rest_g[0][pick], "gj": rest_g[1][pick]})

    def totals_fn(w):
        t = {}
        for key, cs in (("H1", ["H1s", "H1x"]), ("H2", ["H2s", "H2x"]), ("H12", ["H12s", "H12x"]),
                        ("H3j", ["H3j"]), ("rest", rest)):
            t[key] = float(sum(est.pair_group_weight(w, *all_groups[c]).sum() for c in cs))
        return t
    return sig, totals, split, sample, q0, r_true, n_true, n_planted, n_groups, totals_fn


def test_estimand_recovery_toy_world():
    sig, totals, split, sample, q0, r_true, n_true, n_planted, G, totals_fn = toy_pair_world()
    e = est.estimands(sig, totals, split, sample, q0)
    pp = e["per_pair"]
    for c in ("H1s", "H2s", "H3j", "H4"):
        assert pp["r_adj"][c] == pytest.approx(r_true[c], rel=0.08), c
    # N: true rate 0 + planted pairs; FP subtraction must recover planted / n_N
    assert pp["r_adj"]["N"] == pytest.approx(n_planted / n_true["N"], rel=0.5, abs=3e-5)
    true_K = {c: r_true[c] * n_true[c] for c in r_true}
    true_K["N"] = n_planted
    tot = sum(true_K.values())
    true_Fanc = sum(true_K[c] for c in est.ANC) / tot
    true_Fan = true_K["N"] / tot
    assert pp["F_adj"]["F_anc"] == pytest.approx(true_Fanc, rel=0.05)
    assert pp["F_adj"]["F_analog"] == pytest.approx(true_Fan, rel=0.5)
    # without FP subtraction F_analog is inflated by q0 * n_N (the review's point 2)
    assert pp["F_raw"]["F_analog"] > 2 * true_Fan
    # bootstrap CI covers the truth for F_anc
    reps = est.bootstrap(sig, totals_fn, split, sample, q0, G, 40, 0, ("per_pair",))
    lo, hi = np.quantile(reps["per_pair.F_adj.F_anc"], [0.025, 0.975])
    assert lo - 0.01 < true_Fanc < hi + 0.01


def test_no_convergence_world_gives_analog_near_zero_after_subtraction():
    sig, totals, split, sample, q0, r_true, n_true, _, G, _ = toy_pair_world(seed=3, n_planted=0)
    e = est.estimands(sig, totals, split, sample, q0)
    exp_raw = q0 * n_true["N"] / (sum(r_true[c] * n_true[c] for c in r_true) + q0 * sum(n_true.values()))
    assert e["per_pair"]["F_raw"]["F_analog"] == pytest.approx(exp_raw, rel=0.15)  # raw > 0 from null FPs
    assert e["per_pair"]["F_adj"]["F_analog"] < 0.005


def test_family_pair_weighting_counts_cells():
    sig = pd.DataFrame({"i": [0, 1, 2, 3], "j": [4, 5, 6, 7], "cls": ["H1s", "H1s", "N", "N"],
                        "gi": [0, 0, 1, 2], "gj": [0, 0, 3, 4]})
    w = est.pair_weights(sig, "per_family_pair")
    assert list(w) == [0.5, 0.5, 1.0, 1.0]
    ws = est.pair_weights(sig, "per_site")
    assert list(ws) == [1.0, 1.0, 1.0, 1.0]


# ----------------------------------------------------------------------------------------------- directional test
from pocketevo import directional as dr  # noqa: E402


_dir_world = dr.toy_world


def test_newick_and_derived_split():
    t = dr.SpeciesTree("((a,b)x,(c,(d,e)y)z)r;")
    assert t.names[0] == "r" and set(t.leaf_of) == {"a", "b", "c", "d", "e"}
    mem = np.array([t.leaf_of[k] for k in "abcdde"])
    old_min = dr.MIN_DERIVED
    dr.MIN_DERIVED = 2
    try:
        der, out = dr.derived_split(t, t.leaf_of["d"], mem, np.random.default_rng(0))
    finally:
        dr.MIN_DERIVED = old_min
    assert set(der) == {3, 4, 5} and set(out) == {0, 1, 2}      # clade y=(d,e): 3 members, 2 species


@pytest.mark.parametrize("seed", [0, 1])
def test_directional_no_convergence_world_reads_zero(seed):
    r = _dir_world(seed, n_fam=10)
    assert r["F_dir"] <= 0.002 and r["false_pass"] == 0.0


def test_directional_planted_convergence_detected():
    r = _dir_world(0, planted=True, n_fam=10)
    assert r["power"] >= 0.75 and r["false_pass"] <= 0.002


def test_weighted_exact_totals_match_bruteforce():
    """site weights (redundancy reweighting sensitivity): totals = sum over pairs of w_i w_j."""
    rng = np.random.default_rng(11)
    ctx = random_world(rng)
    N = len(ctx.site_prot)
    sw = rng.uniform(0.1, 1.0, N)
    i, j = np.triu_indices(N, 1)
    cls = hom.classify_pairs(ctx, i, j)
    pw = sw[i] * sw[j]
    tot = hom.exact_class_totals(ctx, site_weight=sw)
    assert tot["H3j"] == pytest.approx(pw[cls == hom.H3J].sum())
    assert tot["rest"] == pytest.approx(pw[cls == hom.REST].sum())
    assert tot["H12"] == pytest.approx(pw[cls == hom.H12].sum())
