"""Toy checks for pocketevo.descent (Part A M1/M2) with known answers."""
import numpy as np
import pandas as pd
import pytest

from pocketevo import descent

CRIT = {"plddt": 90, "plddt_s70": 70, "frac_lt70": 0.2, "aligned": 0.8, "len_lo": 0.7, "len_hi": 1.3}
CEIL = descent.Ceiling({(">=90", "<=10"): 0.5, (">=90", ">10"): 0.8, ("70-90", "<=10"): 0.4, ("70-90", ">10"): 0.6},
                       "toy", 0.2)


def rows_from(states: dict, plddt: dict | None = None, nfilt: dict | None = None, big=None):
    """states: {(member, sog): 'P'|'A'|'U'}."""
    recs = []
    for (m, s), st in states.items():
        pl = (plddt or {}).get((m, s), 95.0)
        recs.append({"member": m, "sog": s, "state": st, "state_s70": st if st != "U" or pl < 70 else "A",
                     "footprint_plddt": pl, "frac_lt70": 0.0, "aligned_frac": 1.0, "len_ratio": 1.0,
                     "n_filtered_pockets": (nfilt or {}).get(m, 0), "big_fp": True if big is None else big(s)})
    return pd.DataFrame(recs)


def build(states, n_mem, n_sog, small=None, **kw):
    r = rows_from(states, **kw)
    masks = descent.informative_masks(r, CRIT)
    sm = np.zeros(len(r), bool) if small is None else small(r)
    si = np.arange(len(r))
    M = descent.og_matrices(r.member.to_numpy(), r.sog.to_numpy(), n_mem, n_sog, masks, sm, si)
    return r, M, descent.pair_counts(M, CEIL)


def test_exact_counts_two_members():
    # X (0) is P in 10 SOGs; Y (1) is P in 7 of them and A in 3 -> X->Y 7/10, Y->X 7/7, both directions 14/17
    st = {}
    for s in range(10):
        st[(0, s)] = "P"
        st[(1, s)] = "P" if s < 7 else "A"
    r, M, C = build(st, 2, 10)
    assert C["lit"]["n_inf"][0, 1] == 17 and C["lit"]["n_match"][0, 1] == 14
    assert C["lit"]["n_inf"][0, 1] == C["lit"]["n_inf"][1, 0]
    # ceiling: all large? small=None -> all sites '>10', band >=90 -> 0.8 per informative observation
    assert np.isclose(C["lit"]["ceil_sum"][0, 1], 17 * 0.8)
    assert C["lit"]["n_inf"][0, 0] == 0


def test_quality_matching_removes_low_plddt_P():
    # Y's P calls in SOG 0..4 have pLDDT 80 -> qm90 drops them; lit keeps them; A rows always count
    st, pl = {}, {}
    for s in range(10):
        st[(0, s)] = "P"
        st[(1, s)] = "P" if s < 5 else "A"
        if s < 5:
            pl[(1, s)] = 80.0
    r, M, C = build(st, 2, 10, plddt=pl)
    # lit: X->Y 5/10 informative all, matches 5; Y->X: 5/5 (X high plddt) -> (5+5)/(10+5)
    assert C["lit"]["n_match"][0, 1] / C["lit"]["n_inf"][0, 1] == pytest.approx(10 / 15)
    # qm90: X->Y only the 5 A are informative -> 0/5; Y->X 5/5 -> 5/10
    assert C["qm90"]["n_match"][0, 1] == 5 and C["qm90"]["n_inf"][0, 1] == 10
    # s70: low-pLDDT P passes the 70 stratum -> same as lit here
    assert C["s70"]["n_match"][0, 1] / C["s70"]["n_inf"][0, 1] == pytest.approx(10 / 15)


def test_filtered_model_variant_symmetric():
    st = {(0, 0): "P", (1, 0): "A", (0, 1): "P", (1, 1): "P"}
    _, _, C = build(st, 2, 2, nfilt={1: 3})
    # member 1 has filtered pockets -> never informative as Y in qm90nf; X=1 -> Y=0 still counts
    assert C["qm90nf"]["n_inf"][0, 1] == 1 and C["qm90nf"]["n_match"][0, 1] == 1


def test_matched_site_pairs_ids():
    st = {(0, 0): "P", (1, 0): "P", (2, 0): "A", (0, 1): "P", (2, 1): "P"}
    r, M, C = build(st, 3, 2)
    mm = descent.matched_site_pairs(M, np.array([0, 0, 1]), np.array([1, 2, 2]))
    si = {(m, s): k for k, (m, s) in enumerate(zip(r.member, r.sog))}
    got = {(int(k), int(s)): (int(a), int(b)) for k, s, a, b in mm[["k", "sog", "si_a", "si_b"]].itertuples(index=False)}
    assert got == {(0, 0): (si[(0, 0)], si[(1, 0)]), (1, 1): (si[(0, 1)], si[(2, 1)])}


def test_synthetic_known_match_rates_and_effect():
    """Synthetic pair set: P(match) = 0.6 for orthologs and 0.4 for paralogs at every identity -> marginal 0.2."""
    rng = np.random.default_rng(0)
    n = 4000
    ident = rng.uniform(0.2, 0.9, n)
    rel = np.where(rng.random(n) < 0.5, descent.REL_ORTHO, descent.REL_PARA)
    trials = rng.integers(5, 30, n)
    p = np.where(rel == descent.REL_ORTHO, 0.6, 0.4)
    succ = rng.binomial(trials, p)
    cells = pd.DataFrame({"identity": ident, "rel": rel, "succ": succ, "trials": trials})
    eff = descent.ortholog_effect(cells, np.array([0.3, 0.5, 0.7]))
    assert eff["marginal_diff"] == pytest.approx(0.2, abs=0.02)
    assert eff["log_or"] == pytest.approx(np.log(1.5) - np.log(2 / 3), abs=0.08)
    b = descent.binned_effect(cells, np.linspace(0, 1, 11))
    assert b["std_diff"] == pytest.approx(0.2, abs=0.02)
    # identity confounding: paralogs at low identity with lower rates but no class effect -> effect ~ 0
    rel2 = np.where(ident > 0.55, descent.REL_ORTHO, descent.REL_PARA)
    rel2[rng.random(n) < 0.2] = descent.REL_PARA
    p2 = descent.expit(-2 + 4 * ident)
    cells2 = pd.DataFrame({"identity": ident, "rel": rel2, "succ": rng.binomial(trials, p2), "trials": trials})
    raw = (cells2.succ[rel2 == 1].sum() / cells2.trials[rel2 == 1].sum()
           - cells2.succ[rel2 == 2].sum() / cells2.trials[rel2 == 2].sum())
    assert raw > 0.1                    # naive comparison is confounded by identity
    eff2 = descent.ortholog_effect(cells2, np.array([0.3, 0.5, 0.7]))
    assert abs(eff2["marginal_diff"]) < 0.03


def test_bootstrap_weights_and_sampling():
    codes = np.repeat(np.arange(50), 3)
    w = descent.bootstrap_weights(codes, 50, rep=7, seed=0)
    assert w.sum() == 150 and np.allclose(w.reshape(50, 3).std(1), 0)
    w2 = descent.bootstrap_weights(codes, 50, rep=7, seed=0)
    assert np.array_equal(w, w2)
    k = descent.sample_pairs(np.arange(1000), np.arange(1000), 200, np.random.default_rng(1))
    assert len(k) == 200 and len(set(k)) == 200
    assert len(descent.sample_pairs(np.arange(10), np.arange(10), 200, np.random.default_rng(1))) == 10


def test_identity():
    seqs = ["ACDE-F", "ACDQ-F", "------"]
    mat = np.frombuffer("".join(seqs).encode(), np.uint8).reshape(3, 6)
    ident, nb = descent.pair_identity(mat, np.array([0, 0]), np.array([1, 2]))
    assert nb[0] == 5 and ident[0] == pytest.approx(4 / 5)
    assert nb[1] == 0 and np.isnan(ident[1])


def test_ceiling_loader():
    er = []
    for band in (">=90", "70-90"):
        for size in ("<=10", ">10"):
            er.append({"metric": "redetect", "model_source": "x", "tau": 0.2, "kingdom": "all", "plddt_band": band,
                       "site_size_class": size, "rate": 0.5})
    c = descent.load_ceiling(pd.DataFrame(er), "x", 0.2)
    assert c.value(">=90", ">10") == 0.5
    with pytest.raises(ValueError):
        descent.load_ceiling(pd.DataFrame(er), "y", 0.2)


def test_strata_sum_to_informative_and_rule_matched_ceiling():
    st = {(0, s): "P" for s in range(6)}
    st.update({(1, s): ("P" if s < 4 else "A") for s in range(6)})
    r, M, C = build(st, 2, 6, small=lambda r: r.sog.to_numpy() < 2)
    tot = sum(C["lit"]["strata"][k] for k in descent.STRATA)
    assert np.array_equal(tot, C["lit"]["n_inf"])
    assert np.isclose(descent.ceil_from_strata(C["lit"]["strata"], CEIL)[0, 1], C["lit"]["ceil_sum"][0, 1])
    red = pd.DataFrame({"comparison": "x", "n_res": [5, 5, 5, 20, 20, 20, 20] * 2,
                        "band": [">=90"] * 7 + ["70-90"] * 7,
                        "best_jaccard": [1.0, 0.0, 0.1, 0.5, 0.5, 0.5, 0.0] * 2})
    c = descent.rule_matched_ceiling(red, "x", 0.2)
    assert c.value(">=90", "<=10") == pytest.approx(0.5)      # partial (0.1) leaves the denominator
    assert c.value(">=90", ">10") == pytest.approx(0.75)


def test_d11_asymmetry_inflates_m1_and_symmetric_rule_does_not():
    """Known answer: Y keeps the pocket in exactly half of X's SOGs, independent of footprint size.

    40 SOGs where X is P; SOGs 0..19 have small footprints (5 residues), 20..39 large (20 residues). Y is P in the
    even-numbered SOGs and absent (A) in the odd ones, so the true retention rate is 0.5 in both size classes.
    D11 (O1_r2 / O2_r1 states) turns A on <= 10-residue footprints into U, but P has no such rule: qm90 then counts
    10 small-footprint P matches with no small-footprint non-matches -> inflated M1 (X->Y 20/30). qm90d11 applies
    the > 10-residue rule to P too -> 10/20 = 0.5.
    """
    st, recs = {}, []
    for s in range(40):
        small_fp = s < 20
        y = "P" if s % 2 == 0 else ("U" if small_fp else "A")      # D11: small-footprint A -> U
        for m, state in ((0, "P"), (1, y)):
            recs.append({"member": m, "sog": s, "state": state, "state_s70": state, "footprint_plddt": 95.0,
                         "frac_lt70": 0.0, "aligned_frac": 1.0, "len_ratio": 1.0, "n_filtered_pockets": 0,
                         "big_fp": True, "fp_res": 5 if small_fp else 20})
    r = pd.DataFrame(recs)
    masks = descent.informative_masks(r, CRIT)
    M = descent.og_matrices(r.member.to_numpy(), r.sog.to_numpy(), 2, 40, masks, np.zeros(len(r), bool),
                            np.arange(len(r)))
    C = descent.pair_counts(M, CEIL)
    asym = C["qm90"]["n_match"][0, 1] / C["qm90"]["n_inf"][0, 1]
    sym = C["qm90d11"]["n_match"][0, 1] / C["qm90d11"]["n_inf"][0, 1]
    # pair totals count both directions; X->Y part alone: asym 20/30, sym 10/20; Y->X adds 20/20 resp. 10/10
    assert asym > 0.7          # (20 + 20) / (30 + 20) = 0.8: inflated above the true 0.5 + Y->X P/P
    assert sym == pytest.approx((10 + 10) / (20 + 10))
    assert asym > sym
    # X->Y direction only (Y informative given X present): exactly the true rate under the symmetric rule
    xy_sym = (M["P"] @ M["IP_qm90d11"].T)[0, 1] / (M["P"] @ M["I_qm90d11"].T)[0, 1]
    xy_asym = (M["P"] @ M["IP_qm90"].T)[0, 1] / (M["P"] @ M["I_qm90"].T)[0, 1]
    assert xy_sym == pytest.approx(0.5) and xy_asym == pytest.approx(20 / 30)
