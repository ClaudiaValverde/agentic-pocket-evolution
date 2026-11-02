"""Sanity checks for the atlas snapshot (known answers + toy positive/negative controls).

Run:
  uv run --with pytest --with pyarrow --with duckdb --with pandas pytest tests/test_atlas_snapshot.py -v

The snapshot under test is $ATLAS_SNAPSHOT (a directory containing manifest.parquet) or, by default, the
most recent directory in data/snapshot/. Checks that touch pocket rows use random unit samples with a
fixed seed, so they stay cheap on the complete atlas.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import random
import subprocess
import sys

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 20261030
SNAP_SCRIPT = os.path.join(P, "scripts", "atlas_snapshot.py")


def _snapshot_dir() -> str:
    d = os.environ.get("ATLAS_SNAPSHOT")
    if d:
        return d
    cands = sorted(glob.glob(os.path.join(P, "data", "snapshot", "*", "manifest.parquet")))
    if not cands:
        pytest.skip("no snapshot found; run scripts/atlas_snapshot.py first")
    return os.path.dirname(cands[-1])


@pytest.fixture(scope="module")
def snap():
    d = _snapshot_dir()
    summary = json.load(open(os.path.join(d, "manifest.json")))
    man = pq.read_table(os.path.join(d, "manifest.parquet")).to_pandas()
    return dict(dir=d, summary=summary, man=man, atlas=summary["atlas_root"])


@pytest.fixture(scope="module")
def con():
    c = duckdb.connect()
    c.execute("SET threads=32")
    return c


def _files(snap, kind, tiers=None, n=None, seed=SEED):
    m = snap["man"]
    if tiers is not None:
        m = m[m.tier.isin(tiers)]
    paths = [os.path.join(snap["atlas"], p) for p in m[f"{kind}_path"]]
    if n is not None and len(paths) > n:
        paths = random.Random(seed).sample(paths, n)
    return paths


# ---------------------------------------------------------------- manifest integrity
def test_manifest_sha256(snap):
    h = hashlib.sha256(open(os.path.join(snap["dir"], "manifest.parquet"), "rb").read()).hexdigest()
    assert h == snap["summary"]["manifest_sha256"]


def test_manifest_units_unique_and_complete(snap):
    m = snap["man"]
    assert m.unit_id.is_unique
    assert len(m) == snap["summary"]["n_units"]
    assert m.rows_match_json.all(), "parquet row counts disagree with done json counts"
    # tiers 1 and 2 are complete per the guide
    pt = snap["summary"]["per_tier"]
    for t in ("1", "2"):
        assert pt[t]["units"] == pt[t]["worklist_units"]
    assert pt["1"]["structures"] == 555_846  # POCKETOME §2
    assert pt["2"]["structures"] == 436_886


def test_manifest_files_unchanged(snap):
    """Frozen files must still exist with the recorded size (atlas files are written once)."""
    m = snap["man"].sample(n=min(500, len(snap["man"])), random_state=SEED)
    for r in m.itertuples():
        for kind in ("pockets", "structures"):
            p = os.path.join(snap["atlas"], getattr(r, f"{kind}_path"))
            assert os.path.getsize(p) == getattr(r, f"{kind}_size"), p
            assert pq.ParquetFile(p).metadata.num_rows == getattr(r, f"{kind}_rows"), p
        assert os.path.exists(os.path.join(snap["atlas"], r.done_json))


def test_single_method_record(snap):
    assert len(snap["summary"]["method_config_sha256"]) == 1


# ---------------------------------------------------------------- known answers
def test_human_tier1_counts(snap, con):
    files = _files(snap, "structures", tiers=[1])
    r = con.execute("""SELECT count(*), count(DISTINCT accession), count(DISTINCT taxid), min(taxid), max(taxid)
                       FROM read_parquet(?, union_by_name=true) WHERE organism='Homo sapiens'""", [files]).fetchone()
    n_struct, n_acc, n_tax, tmin, tmax = r
    assert 20_000 <= n_acc <= 23_000, n_acc           # human reference proteome ~20.5k entries
    assert n_acc <= n_struct <= 26_000, n_struct      # plus fragments of long proteins
    assert n_tax == 1 and tmin == tmax == 9606


def test_human_rows_taxid_9606(snap, con):
    """All rows of the HUMAN worklist units carry taxid 9606 and organism Homo sapiens, in both tables."""
    wl = pq.read_table(os.path.join(snap["atlas"], "worklist.parquet"), columns=["unit_id", "code"]).to_pandas()
    human_units = set(wl.unit_id[wl.code == "HUMAN"])
    assert human_units
    m = snap["man"][snap["man"].unit_id.isin(human_units)]
    for kind in ("structures", "pockets"):
        files = [os.path.join(snap["atlas"], p) for p in m[f"{kind}_path"]]
        bad = con.execute("""SELECT count(*) FROM read_parquet(?, union_by_name=true)
                             WHERE taxid<>9606 OR organism<>'Homo sapiens'""", [files]).fetchone()[0]
        assert bad == 0, kind


def test_tier_model_version_source(snap, con):
    """tier / model_version / source agree with POCKETOME §2-3.

    Known discrepancy with the guide (observed 2026-10-30): tier-1 *structure* tables carry the per-proteome tar
    name (e.g. 'UP000005640_9606_HUMAN_v6.tar') in `source`, not 'afdb_v6_tar'; tier-1 pocket tables have no
    `source` column at all. Both are accepted here as 'v6 tar' and must be normalised downstream.
    """
    for tier, mv, ok_src in [(1, 6, "(coalesce(source,'afdb_v6_tar') = 'afdb_v6_tar' OR regexp_matches(source, '^UP[0-9]+_[0-9]+_[A-Z0-9]+_v6\\.tar$'))"),
                             (2, 6, "coalesce(source,'afdb_swissprot_v6_tar') = 'afdb_swissprot_v6_tar'"),
                             (3, 4, "coalesce(source,'') = 'foldcomp:afdb_uniprot_v4'")]:
        for kind in ("structures", "pockets"):
            files = _files(snap, kind, tiers=[tier], n=300)
            if kind == "pockets" and tier == 3:   # empty pocket files lack the column
                files = [f for f in files if pq.ParquetFile(f).metadata.num_rows > 0]
            has_src = any("source" in pq.read_schema(f).names for f in files)
            src_expr = "source" if has_src else "NULL::VARCHAR AS source"
            if tier == 3:
                assert has_src, "tier-3 files must carry `source`"
            r = con.execute(f"""SELECT count(*) FILTER (WHERE tier<>{tier}), count(*) FILTER (WHERE model_version<>{mv}),
                       count(*) FILTER (WHERE NOT {ok_src})
                       FROM (SELECT tier, model_version, {src_expr} FROM read_parquet(?, union_by_name=true))""",
                              [files]).fetchone()
            assert r == (0, 0, 0), (tier, kind, r)


# ---------------------------------------------------------------- pocket table checks
def test_pocket_id_unique_tiers12(snap, con):
    files = _files(snap, "pockets", tiers=[1, 2])
    n, nd = con.execute("SELECT count(*), count(DISTINCT pocket_id) FROM read_parquet(?, union_by_name=true)",
                        [files]).fetchone()
    assert n == nd and n > 3_000_000


def test_pocket_id_unique_and_wellformed_sample(snap, con):
    files = _files(snap, "pockets", n=400)
    r = con.execute("""SELECT count(*), count(DISTINCT pocket_id),
          count(*) FILTER (WHERE NOT starts_with(pocket_id,
              'AF-' || accession || '-F' || fragment || '-model_v' || model_version || ':m'))
        FROM read_parquet(?, union_by_name=true)""", [files]).fetchone()
    assert r[0] == r[1], "duplicate pocket_id in sample"
    assert r[2] == 0, "pocket_id inconsistent with accession/fragment/model_version"


def test_pockets_pass_filters(snap, con):
    files = _files(snap, "pockets", n=400)
    r = con.execute("""SELECT min(burrow_score), min(plddt_mean), min(n_points)
                       FROM read_parquet(?, union_by_name=true)""", [files]).fetchone()
    assert r[0] >= 0.03 - 1e-9 and r[1] >= 70.0 - 1e-9 and r[2] >= 10


def test_fingerprints_unit_norm_10k(snap):
    rng = np.random.default_rng(SEED)
    files = _files(snap, "pockets", n=600)
    vecs, got, per = [], 0, 40
    for f in files:
        if got >= 10_000:
            break
        pf = pq.ParquetFile(f)
        if pf.metadata.num_rows == 0:
            continue
        t = pf.read_row_group(int(rng.integers(pf.num_row_groups)), columns=["fingerprint"])
        idx = rng.choice(t.num_rows, size=min(per, t.num_rows), replace=False)
        col = t.column("fingerprint").take(pa.array(idx)).combine_chunks()
        vecs.append(np.asarray(col.flatten().to_numpy(zero_copy_only=False), dtype=np.float32).reshape(-1, 384))
        got += len(vecs[-1])
    V = np.concatenate(vecs)[:10_000]
    assert len(V) == 10_000
    norms = np.linalg.norm(V, axis=1)
    assert np.all(np.isfinite(V))
    assert np.abs(norms - 1).max() < 5e-3, norms.min()   # float16 storage
    # negative control: vectors are not degenerate
    assert np.median(V[:1000] @ V[1000:2000].T) < 0.99


def test_structure_pocket_counts_consistent(snap, con):
    m = snap["man"].sample(n=min(200, len(snap["man"])), random_state=SEED)
    for r in m.itertuples():
        s = pq.read_table(os.path.join(snap["atlas"], r.structures_path), columns=["n_pockets_kept", "n_pockets_all"])
        kept = int(pa.compute.sum(s["n_pockets_kept"]).as_py() or 0)
        allp = int(pa.compute.sum(s["n_pockets_all"]).as_py() or 0)
        assert kept == r.pockets_rows, r.unit_id
        assert kept <= allp


def test_accession_in_one_tier_only(snap, con):
    f1 = _files(snap, "structures", tiers=[1])
    f2 = _files(snap, "structures", tiers=[2])
    f3 = _files(snap, "structures", tiers=[3], n=300)
    con.execute("CREATE OR REPLACE TEMP TABLE a12 AS SELECT DISTINCT accession, tier FROM read_parquet(?, union_by_name=true)",
                [f1 + f2])
    assert con.execute("SELECT count(*) FROM (SELECT accession FROM a12 GROUP BY 1 HAVING count(*)>1)").fetchone()[0] == 0
    n = con.execute("""SELECT count(*) FROM read_parquet(?, union_by_name=true) s JOIN a12 USING (accession)""",
                    [f3]).fetchone()[0]
    assert n == 0


# ---------------------------------------------------------------- toy atlas: positive/negative controls
def _toy_unit(root, unit, tier, n_struct, n_pockets, write_done=True, corrupt=None):
    os.makedirs(os.path.join(root, "state", "done"), exist_ok=True)
    for kind, n in (("pockets", n_pockets), ("structures", n_struct)):
        d = os.path.join(root, kind, f"tier{tier}")
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, f"{unit:06d}.parquet")
        if corrupt == kind:
            open(p, "wb").write(b"not a parquet file")
        else:
            pq.write_table(pa.table({"unit_id": [unit] * n, "tier": [tier] * n}), p)
    if write_done:
        json.dump({"unit_id": unit, "tier": tier, "structures": n_struct, "structures_ok": n_struct,
                   "structures_error": 0, "pockets_all": n_pockets + 1, "pockets_kept": n_pockets,
                   "time": "2026-10-30T00:00:00"},
                  open(os.path.join(root, "state", "done", f"{unit}.json"), "w"))


def test_snapshot_toy_atlas(tmp_path):
    root = tmp_path / "atlas"
    _toy_unit(str(root), 0, 1, 5, 7)                                   # good
    _toy_unit(str(root), 300, 3, 3, 0)                                 # good, empty pocket file
    _toy_unit(str(root), 301, 3, 3, 4, write_done=False)               # in progress: no done json
    _toy_unit(str(root), 302, 3, 3, 4, corrupt="pockets")              # unreadable parquet
    _toy_unit(str(root), 303, 3, 3, 4)
    os.remove(root / "structures" / "tier3" / "000303.parquet")        # missing structures file
    pq.write_table(pa.table({"tier": [1, 3, 3, 3, 3]}), root / "worklist.parquet")
    out = tmp_path / "snap"
    subprocess.run([sys.executable, SNAP_SCRIPT, "--atlas-root", str(root), "--out", str(out),
                    "--snapshot-id", "toy", "--threads", "2"], check=True, capture_output=True)
    man = pq.read_table(out / "toy" / "manifest.parquet").to_pandas()
    exc = pq.read_table(out / "toy" / "excluded.parquet").to_pandas()
    summ = json.load(open(out / "toy" / "manifest.json"))
    assert sorted(man.unit_id) == [0, 300]
    assert man.set_index("unit_id").loc[0, "pockets_rows"] == 7
    assert man.set_index("unit_id").loc[0, "structures_rows"] == 5
    assert man.rows_match_json.all()
    reasons = dict(zip(exc.unit_id, exc.excluded_reason))
    assert reasons[302].startswith("unreadable_pockets_parquet")
    assert reasons[303] == "missing_structures_parquet"
    assert reasons[301].startswith("no_done_json")
    assert summ["per_tier"]["3"]["worklist_units"] == 4
    assert summ["manifest_sha256"] == hashlib.sha256(open(out / "toy" / "manifest.parquet", "rb").read()).hexdigest()
