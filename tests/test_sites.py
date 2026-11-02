"""Tests for pocketevo.sites and scripts/build_sites.py (toy data with known answers + real-atlas checks).

Run:  .venv/bin/python -m pytest tests/test_sites.py -v
Real-atlas tests are skipped when the atlas / AFDB tars are not mounted.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pytest

from pocketevo import sites as S

P = Path(__file__).resolve().parents[1]
ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")
HUMAN_TAR = Path("/cluster/work/beltrao/shared/alphafold/v6/UP000005640_9606_HUMAN_v6.tar")
SNAP = P / "data" / "snapshot" / "20261030T085136Z"
needs_atlas = pytest.mark.skipif(not ATLAS.exists(), reason="atlas not mounted")


def toy(rows):
    """rows: (accession, fragment, pocket_suffix, score, residue list)"""
    n = len(rows)
    rng = np.random.default_rng(0)
    fp = rng.standard_normal((n, 384)).astype(np.float16)
    return pa.table({
        "unit_id": pa.array([0] * n, pa.int32()), "tier": pa.array([1] * n, pa.int8()),
        "accession": [r[0] for r in rows], "fragment": pa.array([r[1] for r in rows], pa.int16()),
        "model_version": pa.array([6] * n, pa.int8()), "taxid": pa.array([9606] * n, pa.int32()),
        "organism": ["Homo sapiens"] * n,
        "pocket_id": [f"AF-{r[0]}-F{r[1]}-model_v6:m{r[2]}" for r in rows],
        "burrow_score": [float(r[3]) for r in rows], "n_points": [10] * n, "volume_A3": [10.0] * n,
        "site_volume_A3": [5.0] * n, "residues": [" ".join(f"A:{x}" for x in r[4]) for r in rows],
        "plddt_mean": [90.0 - i for i in range(n)], "plddt_min": [70.0 - i for i in range(n)],
        "plddt_frac_lt70": [0.01 * i for i in range(n)],
        "fingerprint": pa.FixedSizeListArray.from_arrays(pa.array(fp.ravel(), pa.float16()), 384),
    })


def site_of(ps: pa.Table) -> dict:
    return dict(zip(ps["pocket_id"].to_pylist(), ps["site_id"].to_pylist()))


# ------------------------------------------------------------------------------------- parsing / numbering
def test_parse_residues():
    row, num = S.parse_residues(pa.array(["A:53 A:79 A:82", "", None, "A:1"]))
    assert row.tolist() == [0, 0, 0, 3] and num.tolist() == [53, 79, 82, 1]
    with pytest.raises(ValueError):
        S.parse_residues(pa.array(["A:53A"]))
    with pytest.raises(ValueError):
        S.parse_residues(pa.array(["53"]))


def test_fragment_offset_rule():
    assert S.fragment_offset(1) == 0 and S.fragment_offset(2) == 200 and S.fragment_offset(9) == 1600
    assert S.fragment_offset(np.array([1, 3])).tolist() == [0, 400]


# ------------------------------------------------------------------------------------- merging (known answers)
def r(a, b):
    return list(range(a, b))


def test_chain_of_overlaps_merges_transitively():
    # A={1..10}, B={6..15}: J=5/15=0.333; B,C={11..20}: 0.333; A,C: 0 -> one site via B
    t = toy([("X", 1, 0, 0.9, r(1, 11)), ("X", 1, 1, 0.5, r(6, 16)), ("X", 1, 2, 0.2, r(11, 21)),
             ("X", 1, 3, 0.95, r(100, 110))])
    sites, fps, ps = S.build_sites(t)
    m = site_of(ps)
    assert m["AF-X-F1-model_v6:m0"] == m["AF-X-F1-model_v6:m1"] == m["AF-X-F1-model_v6:m2"]
    assert m["AF-X-F1-model_v6:m3"] != m["AF-X-F1-model_v6:m0"]
    # numbering by max score: m3 (0.95) -> s0, chain (0.9) -> s1
    assert m["AF-X-F1-model_v6:m3"] == "X:s0" and m["AF-X-F1-model_v6:m0"] == "X:s1"
    s1 = sites.filter(pc.equal(sites["site_id"], "X:s1")).to_pylist()[0]
    assert s1["n_pockets"] == 3 and s1["residues"] == r(1, 21) and s1["n_res"] == 20
    assert s1["member_pocket_ids"] == [f"AF-X-F1-model_v6:m{i}" for i in (0, 1, 2)]
    assert s1["rep_pocket_id"] == "AF-X-F1-model_v6:m0" and s1["max_burrow_score"] == pytest.approx(0.9)
    assert s1["mean_burrow_score"] == pytest.approx((0.9 + 0.5 + 0.2) / 3)
    assert s1["plddt_mean_min"] == pytest.approx(88.0) and s1["plddt_min_min"] == pytest.approx(68.0)
    assert s1["volume_sum"] == pytest.approx(30.0) and s1["res_min"] == 1 and s1["res_max"] == 20
    # rep fingerprint is the rep pocket's
    fp_s1 = fps.filter(pc.equal(fps["site_id"], "X:s1"))["fingerprint"].to_pylist()[0]
    assert fp_s1 == t["fingerprint"].to_pylist()[0]


@pytest.mark.parametrize("inter,only_a,only_b,expect", [
    (3, 4, 3, True),    # 3/10 = 0.300 exactly -> merge (>=)
    (4, 5, 4, True),    # 4/13 = 0.3077 just above
    (2, 3, 2, False),   # 2/7  = 0.2857 just below
    (29, 34, 34, False),  # 29/97 = 0.29897 just below
    (0, 5, 5, False),
])
def test_threshold_boundary(inter, only_a, only_b, expect):
    a = r(1, 1 + inter + only_a)
    b = r(1 + only_a, 1 + only_a + inter + only_b)
    assert len(set(a) & set(b)) == inter
    J = inter / len(set(a) | set(b))
    t = toy([("X", 1, 0, 0.5, a), ("X", 1, 1, 0.4, b)])
    sites, _, _ = S.build_sites(t, jaccard=0.3)
    assert (sites.num_rows == 1) == expect, J


def test_threshold_is_a_parameter():
    t = toy([("X", 1, 0, 0.5, r(1, 6)), ("X", 1, 1, 0.4, r(4, 11))])   # J = 2/10 = 0.2
    assert S.build_sites(t, jaccard=0.3)[0].num_rows == 2
    assert S.build_sites(t, jaccard=0.2)[0].num_rows == 1


def test_no_merge_across_accessions():
    t = toy([("X", 1, 0, 0.5, r(1, 11)), ("Y", 1, 0, 0.4, r(1, 11))])
    sites, _, ps = S.build_sites(t)
    assert sites.num_rows == 2 and set(sites["site_id"].to_pylist()) == {"X:s0", "Y:s0"}


def test_cross_fragment_merge_after_offset():
    # F2 local 1..10 = UniProt 201..210 == F1 201..210 (duplicate in the overlap window)
    t = toy([("L", 1, 0, 0.3, r(201, 211)), ("L", 2, 0, 0.6, r(1, 11)), ("L", 2, 1, 0.2, r(201, 211))])
    sites, _, ps = S.build_sites(t)
    m = site_of(ps)
    assert m["AF-L-F1-model_v6:m0"] == m["AF-L-F2-model_v6:m0"] == "L:s0"
    s0 = sites.to_pylist()[0]
    assert s0["residues"] == r(201, 211) and s0["fragments"] == [1, 2] and s0["n_fragments"] == 2
    assert s0["rep_pocket_id"] == "AF-L-F2-model_v6:m0" and s0["rep_fragment"] == 2
    # F2 local 201..210 = UniProt 401..410: separate site
    assert m["AF-L-F2-model_v6:m1"] == "L:s1"
    assert sites.to_pylist()[1]["residues"] == r(401, 411)


def test_toy_determinism_under_row_shuffle():
    rows = [("X", 1, i, s, r(st, st + 10)) for i, (s, st) in
            enumerate([(0.9, 1), (0.5, 6), (0.5, 50), (0.2, 11), (0.7, 200), (0.7, 300)])] + \
           [("Y", 2, 0, 0.1, r(1, 5))]
    t = toy(rows)
    a = S.build_sites(t)
    perm = np.random.default_rng(1).permutation(t.num_rows)
    b = S.build_sites(t.take(pa.array(perm)))
    assert a[0].equals(b[0]) and a[1].equals(b[1])
    # tie on score 0.7 broken by rep pocket_id: m4 before m5
    m = site_of(a[2])
    assert int(m["AF-X-F1-model_v6:m4"].split(":s")[1]) < int(m["AF-X-F1-model_v6:m5"].split(":s")[1])


def test_empty_input():
    t = toy([("X", 1, 0, 0.5, r(1, 3))]).slice(0, 0)
    sites, fps, ps = S.build_sites(t)
    assert sites.num_rows == 0 and sites.schema.equals(S.SITES_SCHEMA) and fps.schema.equals(S.FP_SCHEMA)


# ------------------------------------------------------------------------------------- real data
def _check_invariants(t: pa.Table, sites: pa.Table, ps: pa.Table, jaccard=0.3):
    pids = t["pocket_id"].to_pylist()
    members = [p for lst in sites["member_pocket_ids"].to_pylist() for p in lst]
    assert sorted(members) == sorted(pids) and len(set(members)) == len(members)  # exactly one site each
    assert int(pc.sum(sites["n_pockets"]).as_py()) == t.num_rows
    assert len(sites["site_id"].unique()) == sites.num_rows
    df = t.select(["pocket_id", "accession", "fragment", "residues", "burrow_score"]).to_pandas()
    df["site_id"] = df.pocket_id.map(site_of(ps))
    assert (df.site_id.str.split(":s").str[0] == df.accession).all()
    srec = sites.to_pandas().set_index("site_id")
    df["pos"] = [frozenset(int(x.split(":")[1]) + 200 * (f - 1) for x in res.split()) for res, f in
                 zip(df.residues, df.fragment)]
    for sid, g in df.groupby("site_id"):
        s = srec.loc[sid]
        assert sorted(set().union(*g.pos)) == list(s.residues)
        assert g.burrow_score.max() == pytest.approx(s.max_burrow_score)
        assert s.rep_pocket_id == g.sort_values(["burrow_score", "pocket_id"], ascending=[False, True]).pocket_id.iloc[0]
    # sites of one accession are numbered by descending max score
    for acc, g in srec.reset_index().groupby("accession"):
        g = g.sort_values("site_index")
        assert list(g.site_index) == list(range(len(g)))
        assert np.all(np.diff(g.max_burrow_score.to_numpy()) <= 0)
    # sites of the same accession: no pair of member pockets across two sites has J >= threshold
    for acc, g in df.groupby("accession"):
        if g.site_id.nunique() < 2 or len(g) > 400:
            continue
        P_ = g.pos.tolist()
        sid = g.site_id.tolist()
        for i in range(len(P_)):
            for j in range(i + 1, len(P_)):
                if sid[i] != sid[j]:
                    u = len(P_[i] | P_[j])
                    assert u == 0 or len(P_[i] & P_[j]) / u < jaccard


@needs_atlas
def test_real_unit_invariants_and_determinism():
    path = ATLAS / "pockets/tier2/000200.parquet"
    t = pq.read_table(path, columns=S.POCKET_COLUMNS + ["fingerprint"])
    a = S.build_sites(t)
    b = S.build_sites(t)
    assert a[0].equals(b[0]) and a[1].equals(b[1]) and a[2].equals(b[2])
    _check_invariants(t, a[0], a[2])
    assert a[1]["site_id"].equals(a[0]["site_id"])


@needs_atlas
def test_real_multifragment_protein_P04275():
    """VWF (P04275, 2813 aa, 9 AFDB v6 fragments): offsets from the mmCIF == rule; duplicates merge."""
    if not HUMAN_TAR.exists():
        pytest.skip("AFDB v6 human tar not mounted")
    idx = pd.read_parquet(ATLAS / "tar_index" / f"{HUMAN_TAR.name}.parquet").set_index("name")
    with open(HUMAN_TAR, "rb") as fh:
        for k in (1, 2, 5, 9):
            n = f"AF-P04275-F{k}-model_v6.cif.gz"
            fh.seek(int(idx.at[n, "offset_data"]))
            d = S.struct_ref_seq(fh.read(int(idx.at[n, "size"])))
            assert d["pdbx_db_accession"] == "P04275"
            assert int(d["db_align_beg"]) - int(d["pdbx_auth_seq_align_beg"]) == S.fragment_offset(k)
            assert int(d["pdbx_auth_seq_align_beg"]) == 1
    end9 = int(d["db_align_end"])
    assert end9 == 2813
    t = pq.read_table(ATLAS / "pockets/tier1", columns=S.POCKET_COLUMNS + ["fingerprint"],
                      filters=pc.field("accession") == "P04275")
    assert set(t["fragment"].to_pylist()) >= {1, 2, 9}
    sites, fps, ps = S.build_sites(t)
    _check_invariants(t, sites, ps)
    assert sites.num_rows < t.num_rows                    # overlapping windows -> duplicates merged
    assert max(sites["n_fragments"].to_pylist()) >= 2
    assert max(sites["res_max"].to_pylist()) <= end9      # all positions inside the UniProt sequence


@needs_atlas
def test_fragment_offsets_table_all_ok():
    f = P / "results/sites/20261030T085136Z/fragment_offsets.tsv"
    if not f.exists():
        pytest.skip("run scripts/verify_fragment_offsets.py first")
    df = pd.read_csv(f, sep="\t")
    assert len(df) > 3000 and df.ok_offset.all() and df.ok_window.all()


@needs_atlas
def test_build_script_resume_and_determinism(tmp_path):
    if not (SNAP / "manifest.parquet").exists():
        pytest.skip("snapshot missing")
    # units 1 and 2 share the split accession P98161 -> exercises the cross-unit path
    man = pd.read_parquet(SNAP / "manifest.parquet")
    sub = SNAP / "manifest.parquet"
    cmd = [sys.executable, str(P / "scripts/build_sites.py"), "--manifest", str(sub), "--tiers", "1",
           "--workers", "4", "--sample", "0"]
    # restrict to a small manifest copy: units 1, 2, 7
    small = tmp_path / "snap"
    small.mkdir()
    m = man[man.unit_id.isin([1, 2, 7])]
    m.to_parquet(small / "manifest.parquet")
    import hashlib
    mj = json.loads((SNAP / "manifest.json").read_text())
    mj["manifest_sha256"] = hashlib.sha256((small / "manifest.parquet").read_bytes()).hexdigest()
    (small / "manifest.json").write_text(json.dumps(mj))
    cmd[3] = str(small / "manifest.parquet")
    outs = []
    for k in range(2):
        out = tmp_path / f"o{k}"
        subprocess.run(cmd + ["--out", str(out)], check=True, capture_output=True)
        outs.append(out)
    res = subprocess.run(cmd + ["--out", str(outs[0])], check=True, capture_output=True, text=True)
    summ = json.loads((outs[0] / "run_tiers1.json").read_text())
    assert summ["n_skipped"] == 3 and summ["n_processed"] == 0
    tabs = []
    for out in outs:
        s = pa.concat_tables([pq.read_table(f) for f in sorted((out / "sites/tier1").glob("*.parquet"))])
        f = pa.concat_tables([pq.read_table(f) for f in sorted((out / "site_fp/tier1").glob("*.parquet"))])
        tabs.append((s.replace_schema_metadata(None), f.replace_schema_metadata(None)))
    assert tabs[0][0].equals(tabs[1][0]) and tabs[0][1].equals(tabs[1][1])
    s = tabs[0][0]
    # P98161 sites only in unit 1 (owner), none in unit 2; all its pockets accounted for
    p = s.filter(pc.equal(s["accession"], "P98161"))
    assert set(p["unit_id"].to_pylist()) == {1}
    u2 = pq.read_table(outs[0] / "sites/tier1/000002.parquet", columns=["accession"])
    assert "P98161" not in set(u2["accession"].to_pylist())
    n_atlas = pq.read_table(ATLAS / "pockets/tier1", columns=["pocket_id"],
                            filters=pc.field("accession") == "P98161").num_rows
    assert int(pc.sum(p["n_pockets"]).as_py()) == n_atlas
    # global: every pocket of the 3 units in exactly one site
    pk = pa.concat_tables([pq.read_table(ATLAS / r_, columns=["pocket_id"]) for r_ in m.pockets_path])
    mem = pc.list_flatten(s["member_pocket_ids"])
    assert len(mem) == pk.num_rows and pc.count_distinct(mem).as_py() == pk.num_rows
