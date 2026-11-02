"""Tests for the observation-error tooling (src/pocketevo/observation.py, src/pocketevo/rescue.py).

  .venv/bin/python -m pytest -q tests/test_observation.py
Toy sites with known answers; one Foldcomp entry round trip (skipped if the atlas DB is not mounted);
method-record hash against an atlas file.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "src"))
from pocketevo import observation as O  # noqa: E402
from pocketevo import rescue as R  # noqa: E402

ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")
FC_DB = ATLAS / "foldcomp/afdb_uniprot_v4"
FC_T12 = Path("/tmp/pe_big/obs/fc_t12_entries.parquet")


# ------------------------------------------------------------------ matching on toy sites
def test_parse_residues():
    assert O.parse_residues("A:53 A:79 A:82") == frozenset({53, 79, 82})
    assert O.parse_residues("A:10A A:11") == frozenset({10, 11})        # insertion code dropped
    assert O.parse_residues("A:1 A:2", offset=200) == frozenset({201, 202})
    assert O.parse_residues([3, 4]) == frozenset({3, 4})
    assert O.parse_residues(None) == frozenset()


def test_jaccard_coverage():
    assert O.jaccard({1, 2, 3}, {2, 3, 4}) == pytest.approx(0.5)
    assert O.jaccard(set(), set()) == 0.0
    assert O.coverage({1, 2, 3, 4}, {3, 4, 5, 6, 7}) == pytest.approx(0.5)
    assert O.coverage(set(), {1}) == 0.0


def test_match_sites_known_answer():
    ref = [{1, 2, 3, 4}, {10, 11, 12}, {20, 21}]
    qry = [{1, 2, 3, 5}, {11, 12, 13, 14, 15, 16}, {30, 31}]
    m = O.match_sites(ref, qry, tau=0.3)
    assert m[0][:2] == (True, 0) and m[0][2] == pytest.approx(3 / 5)
    assert m[1][0] is False and m[1][1] == -1 and m[1][2] == pytest.approx(2 / 7)   # below tau
    assert m[2] == (False, -1, 0.0)                                                 # no overlap at all
    m2 = O.match_sites(ref, qry, tau=0.25)
    assert m2[1][:2] == (True, 1)                                                   # tau is a parameter
    # best match is the max-Jaccard one, not the first overlapping one
    assert O.best_match([{1, 2, 3}], [{1, 9, 8, 7}, {1, 2, 3, 4}])[0][0] == 1
    # empty query
    assert O.match_sites([{1}], [], 0.3) == [(False, -1, 0.0)]


def test_rescue_hit_rules_and_k():
    pockets = [(0.5, frozenset({1, 2, 3, 4, 5, 6, 7, 8})), (0.1, frozenset({50, 51}))]
    fp = {1, 2, 3}
    j = R.setting_hit(pockets, fp, tau=0.3, rule="jaccard")
    c = R.setting_hit(pockets, fp, tau=0.3, rule="coverage", min_cov=0.5)
    assert j["jaccard"] == pytest.approx(3 / 8) and j["hit"]
    assert c["coverage"] == 1.0 and c["hit"]
    assert not R.setting_hit(pockets, {100, 101}, tau=0.3)["hit"]
    assert not R.setting_hit(pockets, {1, 20, 21, 22, 23, 24, 25, 26, 27}, tau=0.3)["hit"]  # Jaccard 1/16
    hits = {k: {"hit": i < 2} for i, k in enumerate(R.SETTING_NAMES)}
    assert R.n_hits(hits) == 2 and R.rescued(hits, 2) and not R.rescued(hits, 3)
    assert len(R.SETTINGS) == 6       # prereg states.rescue.settings.n_settings


def test_read_burrow_csv(tmp_path):
    f = tmp_path / "x.csv"
    f.write_text("rank,score,center_x,center_y,center_z,residues\n1,0.4,0,0,0,A:1 A:2 A:3\n2,0.01,0,0,0,\n")
    p = R.read_burrow_csv(f)
    assert p[0] == (0.4, frozenset({1, 2, 3})) and p[1] == (0.01, frozenset())


def test_split_half_deterministic():
    assert O.split_half("P04035") == O.split_half("P04035")
    halves = [O.split_half(f"Q{i:05d}") for i in range(2000)]
    assert 0.45 < halves.count("tune") / 2000 < 0.55


def test_region_plddt_and_band():
    num = np.array([1, 2, 3, 4])
    pl = np.array([95.0, 92.0, 60.0, 99.0])
    r = O.region_plddt(num, pl, {1, 2, 3})
    assert r["mean"] == pytest.approx(82.333, abs=1e-3) and r["frac_lt70"] == pytest.approx(1 / 3)
    assert O.plddt_band(69.9) == "<70" and O.plddt_band(70) == "70-90" and O.plddt_band(90) == ">=90"


# ------------------------------------------------------------------ method record hash
def test_method_hash_invariances():
    h = O.method_record_hash(O.ATLAS_METHOD)
    assert h == O.ATLAS_METHOD_HASH
    assert O.method_record_hash({**O.ATLAS_METHOD, "device": "cuda:3", "backend": "ref", "storage": "float16"}) == h
    assert O.method_record_hash({**O.ATLAS_METHOD, "af_mode": False}) != h
    assert O.method_record_hash({**O.ATLAS_METHOD, "config_sha256": "0" * 64}) != h
    assert O.method_record_hash(json.dumps(O.ATLAS_METHOD)) == h


@pytest.mark.skipif(not (ATLAS / "pockets/tier1/000000.parquet").exists(), reason="atlas not mounted")
def test_method_hash_atlas_file():
    assert O.parquet_method_hash(ATLAS / "pockets/tier1/000000.parquet") == O.ATLAS_METHOD_HASH


@pytest.mark.skipif(not Path("/tmp/pe_big/tierL/pockets/tierL").exists(), reason="tier L not built")
def test_method_hash_tierL_files():
    fs = sorted(Path("/tmp/pe_big/tierL/pockets/tierL").glob("[0-9]*.parquet"))[:5]
    if not fs:
        pytest.skip("no tier-L files yet")
    assert {O.parquet_method_hash(f) for f in fs} == {O.ATLAS_METHOD_HASH}


# ------------------------------------------------------------------ Foldcomp round trip
@pytest.mark.skipif(not FC_DB.exists() or not FC_T12.exists(), reason="Foldcomp DB / entry list not available")
def test_foldcomp_roundtrip_one_entry(tmp_path):
    import foldcomp
    import pandas as pd
    e = pd.read_parquet(FC_T12)
    r = e[e.acc == "P00918"].iloc[0]          # human carbonic anhydrase II (tier 1)
    name, pdb = O.foldcomp_read(str(FC_DB), int(r.off), int(r.len))
    assert name.startswith("AF-P00918-F1-model_v4")
    f = tmp_path / "a.pdb"
    f.write_text(pdb)
    seq, num, pl = O.read_model(f)
    assert len(seq) == 260 and seq.startswith("MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDT")
    assert num[0] == 1 and 0 <= pl.min() and pl.max() <= 100
    # compress the decompressed model again -> decompress -> same sequence, coordinates within 0.1 A
    fcz = foldcomp.compress(name, pdb)
    name2, pdb2 = foldcomp.decompress(fcz)
    f2 = tmp_path / "b.pdb"
    f2.write_text(pdb2)
    seq2, num2, pl2 = O.read_model(f2)
    assert seq2 == seq and np.allclose(pl2, pl, atol=0.2)       # pLDDT is quantised again on re-compression
    import gemmi
    a = np.array([[at.pos.x, at.pos.y, at.pos.z] for res in gemmi.read_structure(str(f))[0][0] for at in res])
    b = np.array([[at.pos.x, at.pos.y, at.pos.z] for res in gemmi.read_structure(str(f2))[0][0] for at in res])
    assert a.shape == b.shape and np.sqrt(((a - b) ** 2).sum(1)).max() < 0.5
