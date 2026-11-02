"""Tests for SOG building (src/pocketevo/sog.py), states (states.py) and the model parser (afmodel.py).

Toy world: 4 members (4 species) on a 30-column MSA.
  site1  conserved site at columns 2..6 in A, B, C; in D shifted by one column (3..7)
  site2  insertion site at columns 20..24, present only in B and C (A, D have gaps there)
  site3  A-only site at columns 12..15 (B, C aligned there, D truncated before column 11)
"""
import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "src"))
from pocketevo import afmodel, sog, states  # noqa: E402

PARAMS = sog.prereg_params(yaml.safe_load(open(P / "config/prereg.yaml")))
L = 30
INS = set(range(20, 25))


def _aligned(gaps):
    s, aa = [], "ACDEFGHIKLMNPQRSTVWY"
    for c in range(L):
        s.append("-" if c in gaps else aa[c % 20])
    return "".join(s)


def toy_msa():
    aln = {
        "A": _aligned(INS),
        "B": _aligned(set()),
        "C": _aligned(set()),
        "D": _aligned(INS | set(range(11, L))),   # D ends at column 10 (C-terminal truncation)
    }
    return sog.msa_from_aligned(aln, order=["A", "B", "C", "D"])


def cols_to_res(msa, row, cols):
    """UniProt positions (1-based) of member `row` at MSA columns `cols`."""
    p2c = msa.pos2col[row]
    return [int(np.flatnonzero(p2c == c)[0]) + 1 for c in cols]


def toy_sites(msa):
    sites = [  # (site_id, member row, columns)
        ("A:s0", 0, range(2, 7)), ("B:s0", 1, range(2, 7)), ("C:s0", 2, range(2, 7)), ("D:s0", 3, range(3, 8)),
        ("B:s1", 1, range(20, 25)), ("C:s1", 2, range(20, 25)),
        ("A:s1", 0, range(12, 16)),
    ]
    ids = np.array([s[0] for s in sites])
    sm = np.array([s[1] for s in sites])
    res = [cols_to_res(msa, r, c) for _, r, c in sites]
    return ids, sm, res


def build(tau=0.5, seed=0):
    msa = toy_msa()
    ids, sm, res = toy_sites(msa)
    S, dropped = sog.site_columns(msa, sm, res)
    out = sog.build_sogs(S, sm, [tau], PARAMS["leiden_resolution"], seed, PARAMS["footprint_member_fraction"])
    return msa, ids, sm, S, out[tau]


def groups(ids, labels):
    g = {}
    for i, l in zip(ids, labels):
        g.setdefault(int(l), set()).add(i)
    return sorted(map(frozenset, g.values()), key=lambda x: sorted(x))


# ----------------------------------------------------------------------------------------------- columns / sogs
def test_site_columns_roundtrip():
    msa = toy_msa()
    ids, sm, res = toy_sites(msa)
    S, dropped = sog.site_columns(msa, sm, res)
    assert dropped.sum() == 0
    assert list(S[4].indices) == list(range(20, 25))      # insertion site in B maps to insertion columns
    assert list(S[3].indices) == list(range(3, 8))        # D shifted site


def test_site_residues_beyond_sequence_dropped():
    msa = toy_msa()
    S, dropped = sog.site_columns(msa, np.array([3]), [[1, 2, 500]])
    assert dropped.tolist() == [1] and S[0].nnz == 2


def test_sog_correspondences():
    msa, ids, sm, S, r = build(tau=0.5)
    g = groups(ids, r["labels"])
    assert frozenset({"A:s0", "B:s0", "C:s0", "D:s0"}) in g      # shifted site joins (J = 4/6 >= 0.5)
    assert frozenset({"B:s1", "C:s1"}) in g                       # insertion-column site
    assert frozenset({"A:s1"}) in g                               # singleton
    assert len(g) == 3


def test_shifted_site_drops_out_at_high_tau():
    msa = toy_msa()
    ids, sm, res = toy_sites(msa)
    S, _ = sog.site_columns(msa, sm, res)
    # J(D, others) = 4/6 = 0.667 < 0.7
    out = sog.build_sogs(S, sm, [0.7], 1.0, 0, 0.5)
    g = groups(ids, out[0.7]["labels"])
    assert frozenset({"D:s0"}) in g


def test_footprint_insertion_columns_allowed():
    msa, ids, sm, S, r = build()
    lab = r["labels"]
    k = int(lab[list(ids).index("B:s1")])
    F = r["footprint"]
    cols = set(F.indices[F.indptr[k]:F.indptr[k + 1]].tolist())
    assert cols == INS
    assert (msa.nongap[:, sorted(cols)].mean(axis=0) == 0.5).all()   # half the members have gaps there


def test_footprint_member_fraction():
    msa, ids, sm, S, r = build()
    k = int(r["labels"][0])
    F = r["footprint"]
    cols = set(F.indices[F.indptr[k]:F.indptr[k + 1]].tolist())
    # columns 2..7: 3..6 in all 4 sites, 2 in 3/4, 7 in 1/4 -> 7 excluded at frac 0.5
    assert cols == {2, 3, 4, 5, 6}


def test_no_edges_within_member():
    msa = toy_msa()
    # two identical sites in the same member must not be linked
    S, _ = sog.site_columns(msa, np.array([1, 1]), [[3, 4, 5], [3, 4, 5]])
    i, j, w = sog.jaccard_pairs(S, np.array([1, 1]), 0.2)
    assert len(i) == 0


def test_determinism_and_order_independence():
    a = build(seed=0)[4]["labels"]
    b = build(seed=0)[4]["labels"]
    assert np.array_equal(a, b)
    msa = toy_msa()
    ids, sm, res = toy_sites(msa)
    perm = np.random.default_rng(3).permutation(len(ids))
    S, _ = sog.site_columns(msa, sm[perm], [res[i] for i in perm])
    out = sog.build_sogs(S, sm[perm], [0.5], 1.0, 0, 0.5)
    assert groups(ids[perm], out[0.5]["labels"]) == groups(ids, a)


def test_canonical_labels():
    lab = sog.canonical_labels(np.array([5, 5, 2, 9, 9, 9]))
    assert lab.tolist() == [1, 1, 2, 0, 0, 0]


# ----------------------------------------------------------------------------------------------- states
def member_table(msa, plddt_by_member, seq_ok=None, sptm=None):
    n_res = [len(p) for p in msa.pos2col]
    rows = []
    for i, a in enumerate(msa.ids):
        q = np.full(n_res[i], afmodel.q_threshold(plddt_by_member[a]), dtype=np.uint8)
        rows.append({"accession": a, "n_res": n_res[i], "plddt_q": q,
                     "seq_ok": True if seq_ok is None else seq_ok[a], "status_ok": True,
                     "sptm": None if sptm is None else sptm.get(a)})
    df = pd.DataFrame(rows)
    df["og_median_len"] = float(np.median(n_res))
    return df


def run_states(plddt, tau=0.5, **kw):
    msa, ids, sm, S, r = build(tau=tau)
    mem = member_table(msa, plddt, **kw)
    st = states.compute_states(msa, mem, r["footprint"], S, sm, ids, r["labels"], tau, PARAMS)
    lab = r["labels"]
    k_site1, k_ins, k_a = int(lab[0]), int(lab[list(ids).index("B:s1")]), int(lab[list(ids).index("A:s1")])
    get = lambda k, a: st[(st.sog_index == k) & (st.accession == a)].iloc[0]  # noqa: E731
    return st, get, (k_site1, k_ins, k_a)


def test_states_basic():
    st, get, (k1, kins, ka) = run_states({"A": 95, "B": 95, "C": 80, "D": 95})
    for a in "ABCD":
        assert get(k1, a).state == "P"
    assert get(k1, "D").site_id == "D:s0" and get(k1, "D").in_sog
    # insertion SOG: A and D have no residues at the footprint -> U(aligned)
    assert get(kins, "B").state == "P"
    assert get(kins, "A").state == "U" and get(kins, "A").reason.startswith("aligned<")
    # A-only SOG
    assert get(ka, "A").state == "P"
    b = get(ka, "B")
    assert b.state == "A" and b.reason == "" and b.fail_bits == 0
    c = get(ka, "C")                                   # pLDDT 80: U in primary, A in the >=70 stratum
    assert c.state == "U" and c.reason.startswith("plddt<") and c.state_s70 == "A"
    d = get(ka, "D")                                   # truncated before the footprint
    assert d.state == "U" and d.fail_bits & states.BIT["terminal_truncation"]


def test_states_frac_below_70():
    msa, ids, sm, S, r = build()
    mem = member_table(msa, {"A": 95, "B": 95, "C": 95, "D": 95})
    # B: mean high but 2/4 footprint residues < 70 at the A-only footprint (cols 12..15)
    q = mem.at[1, "plddt_q"].copy()
    pos = cols_to_res(msa, 1, [12, 13])
    q[np.array(pos) - 1] = afmodel.q_threshold(60)
    q[np.array(cols_to_res(msa, 1, [14, 15])) - 1] = 200
    mem.at[1, "plddt_q"] = q
    st = states.compute_states(msa, mem, r["footprint"], S, sm, ids, r["labels"], 0.5, PARAMS)
    ka = int(r["labels"][list(ids).index("A:s1")])
    row = st[(st.sog_index == ka) & (st.accession == "B")].iloc[0]
    assert row.state == "U" and row.fail_bits & states.BIT["frac_lt70>0.2"]
    assert abs(row.footprint_plddt - 80.125) < 0.3


def test_states_seq_mismatch_forces_U():
    st, get, (k1, kins, ka) = run_states({"A": 95, "B": 95, "C": 95, "D": 95},
                                         seq_ok={"A": True, "B": False, "C": True, "D": True})
    assert get(k1, "B").state == "U" and get(k1, "B").reason == "model_sequence_mismatch"
    assert get(ka, "B").state == "U"


def test_states_sptm_overlap():
    msa = toy_msa()
    nB = len(msa.pos2col[1])
    mask = np.zeros(nB, dtype=bool)
    mask[np.array(cols_to_res(msa, 1, [13])) - 1] = True
    st, get, (k1, kins, ka) = run_states({"A": 95, "B": 95, "C": 95, "D": 95}, sptm={"B": mask})
    assert get(ka, "B").state == "U" and get(ka, "B").reason == "sp_tm_overlap"
    assert get(ka, "C").state == "A" and not get(ka, "C").sptm_checked


def test_states_partial_overlap_is_not_absence():
    msa, ids, sm, S, r = build()
    # give C a small site touching the A-only footprint (1 of 4 columns): J = 1/4... make it 1/7
    ids2 = np.append(ids, "C:s2")
    sm2 = np.append(sm, 2)
    res2 = toy_sites(msa)[2] + [cols_to_res(msa, 2, [15, 16, 17, 18])]
    S2, _ = sog.site_columns(msa, sm2, res2)
    out = sog.build_sogs(S2, sm2, [0.5], 1.0, 0, 0.5)[0.5]
    mem = member_table(msa, {"A": 95, "B": 95, "C": 95, "D": 95})
    st = states.compute_states(msa, mem, out["footprint"], S2, sm2, ids2, out["labels"], 0.5, PARAMS)
    ka = int(out["labels"][list(ids2).index("A:s1")])
    row = st[(st.sog_index == ka) & (st.accession == "C")].iloc[0]
    assert row.state == "U" and row.reason == "partial_site_overlap"


def test_states_length_ratio():
    msa, ids, sm, S, r = build()
    mem = member_table(msa, {"A": 95, "B": 95, "C": 95, "D": 95})
    mem.loc[mem.accession == "B", "n_res"] = int(mem.og_median_len.iloc[0] * 1.5)
    st = states.compute_states(msa, mem, r["footprint"], S, sm, ids, r["labels"], 0.5, PARAMS)
    ka = int(r["labels"][list(ids).index("A:s1")])
    assert st[(st.sog_index == ka) & (st.accession == "B")].iloc[0].reason == "length_ratio"


def test_rescue_hook():
    st, get, (k1, kins, ka) = run_states({"A": 95, "B": 95, "C": 95, "D": 95})
    st["sog_id"] = "OG:g" + st.sog_index.astype(str)
    resc = pd.DataFrame({"sog_id": [f"OG:g{ka}", f"OG:g{k1}"], "accession": ["B", "A"],
                         "n_settings_hit": [PARAMS["rescue_k"], PARAMS["rescue_k"] - 1],
                         "n_settings_run": [PARAMS["rescue_n_settings"]] * 2})
    out = states.apply_rescue(st, resc, PARAMS["rescue_k"])
    assert out[(out.sog_id == f"OG:g{ka}") & (out.accession == "B")].state.iloc[0] == "W"
    assert out[(out.sog_id == f"OG:g{k1}") & (out.accession == "A")].state.iloc[0] == "W"
    assert out[(out.sog_id == f"OG:g{ka}") & (out.accession == "C")].state.iloc[0] == "A"


def test_states_deterministic():
    a = run_states({"A": 95, "B": 95, "C": 80, "D": 95})[0]
    b = run_states({"A": 95, "B": 95, "C": 80, "D": 95})[0]
    pd.testing.assert_frame_equal(a, b)


# ----------------------------------------------------------------------------------------------- split / prereg
def test_split_rule():
    import hashlib
    for key in ["P24941", "A0A024R1R8", "Q9Y6K9"]:
        exp = "tune" if int(hashlib.sha1(key.encode()).hexdigest(), 16) % 2 == 0 else "test"
        assert sog.split_half(key) == exp
    assert sog.split_key(["Q1", "A2", "P3"]) == "A2"


def test_prereg_missing_key_raises():
    pr = yaml.safe_load(open(P / "config/prereg.yaml"))
    bad = copy.deepcopy(pr)
    del bad["sogs"]["tau_column_jaccard"]
    with pytest.raises(KeyError):
        sog.prereg_params(bad)


# ----------------------------------------------------------------------------------------------- afmodel
CIF = """data_x
#
_entity_poly.pdbx_seq_one_letter_code_can MKV
#
loop_
_ma_qa_metric_local.label_asym_id
_ma_qa_metric_local.label_comp_id
_ma_qa_metric_local.label_seq_id
_ma_qa_metric_local.metric_id
_ma_qa_metric_local.metric_value
_ma_qa_metric_local.model_id
_ma_qa_metric_local.ordinal_id
A MET 1 2 69.99 1 1
A LYS 2 2 70.00 1 2
A VAL 3 2 90.50 1 3
#
_ma_target_ref_db_details.seq_db_sequence_checksum     ABC
"""


def test_parse_cif_inline_sequence_and_q():
    r = afmodel.parse_cif_text(CIF.encode())
    assert r["sequence"] == "MKV"
    q = r["plddt_q"]
    assert q.tolist() == [139, 140, 181]
    assert (q < afmodel.q_threshold(70)).tolist() == [True, False, False]
    assert abs(afmodel.plddt_from_q(q)[2] - 90.75) < 1e-6


def test_combine_fragments():
    s = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(1700))
    q = np.arange(1700) % 200
    f1 = (s[:1400], q[:1400].astype(np.uint8))
    f2 = (s[200:1600], q[200:1600].astype(np.uint8))
    f3 = (s[400:1700], q[400:1700].astype(np.uint8))
    seq, qq, conflict = afmodel.combine_fragments({1: f1, 2: f2, 3: f3})
    assert seq == s and not conflict and np.array_equal(qq, q.astype(np.uint8))


def test_crc64_known():
    # UniProt P24941 (CDK2) CRC64 checked against its AFDB model header in development
    assert afmodel.crc64("") == "0000000000000000"


def test_plddt_quantisation_is_conservative_for_A():
    # true pLDDT 89.9 everywhere: stored q = 179 (89.5); the midpoint estimate (89.75) is reported, but the
    # A threshold uses the lower bound, so it must be U(plddt<90) - and 90.0 exactly must still pass
    msa, ids, sm, S, r = build()
    for val, exp in [(89.9, "U"), (90.0, "A")]:
        mem = member_table(msa, {"A": 95, "B": 95, "C": 95, "D": 95})
        q = np.full(mem.at[1, "n_res"], int(np.floor(val * 2)), dtype=np.uint8)
        mem.at[1, "plddt_q"] = q
        st = states.compute_states(msa, mem, r["footprint"], S, sm, ids, r["labels"], 0.5, PARAMS)
        ka = int(r["labels"][list(ids).index("A:s1")])
        assert st[(st.sog_index == ka) & (st.accession == "B")].iloc[0].state == exp


def test_D11_min_footprint_residues():
    # A-only footprint has 4 columns: with the D11 rule (> 10 residues) B can no longer be A
    msa, ids, sm, S, r = build()
    mem = member_table(msa, {"A": 95, "B": 95, "C": 95, "D": 95})
    prm = dict(PARAMS, A_min_footprint_residues=10)
    st = states.compute_states(msa, mem, r["footprint"], S, sm, ids, r["labels"], 0.5, prm)
    ka = int(r["labels"][list(ids).index("A:s1")])
    row = st[(st.sog_index == ka) & (st.accession == "B")].iloc[0]
    assert row.state == "U" and row.fail_bits & states.BIT["footprint_residues<=10"]
    assert row.reason == "footprint_residues<=10"
    prm3 = dict(PARAMS, A_min_footprint_residues=3)
    st = states.compute_states(msa, mem, r["footprint"], S, sm, ids, r["labels"], 0.5, prm3)
    assert st[(st.sog_index == ka) & (st.accession == "B")].iloc[0].state == "A"
