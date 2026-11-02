"""Web app tests (FastAPI TestClient). Needs the atlas mount and the protein index
(scripts/webapp_build_index.py); tier-3 tests skip if the node-local tier-3 index is absent."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import gemmi
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from fastapi.testclient import TestClient  # noqa: E402

from webapp.app import create_app  # noqa: E402
from webapp.atlas import bucket_of  # noqa: E402
from webapp.config import Config  # noqa: E402

TIER3_ACC = "A0A633RIM1"  # first entry of tier-3 unit 272 (Foldcomp afdb_uniprot_v4)


@pytest.fixture(scope="module")
def client():
    cfg = Config()
    if 1 not in create_app(cfg).state.atlas.tier_dirs():
        pytest.skip("protein index missing: run scripts/webapp_build_index.py --tiers 1,2")
    return TestClient(create_app(cfg))


def _parse_cif(txt: str):
    doc = gemmi.cif.read_string(txt)
    st = gemmi.make_structure_from_block(doc.sole_block())
    return st


def test_bucket_matches_builder():
    spec = importlib.util.spec_from_file_location("wbi", REPO / "scripts" / "webapp_build_index.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    for a in ("P24941", "A0A633RIM1", "Q9XYZ1"):
        assert f"{m.bucket_of(a):02x}" == bucket_of(a)


def test_protein_cdk2(client):
    r = client.get("/api/protein/P24941")
    assert r.status_code == 200
    d = r.json()
    assert [s["tier"] for s in d["structures"]] == [1]
    assert d["structures"][0]["blob_file"] == "UP000005640_9606_HUMAN_v6.tar"
    assert len(d["pockets"]) >= 1
    p = d["pockets"][0]
    assert p["pocket_id"].startswith("AF-P24941-F1-model_v6:m")
    assert p["burrow_score"] >= 0.03 and p["plddt_mean"] >= 70  # atlas filters
    assert all(x.startswith("A:") for x in p["residues"].split())
    html = client.get("/protein/P24941").text
    assert "Pockets (" in html and 'class="pocket"' in html and "3Dmol-min.js" in html


def test_protein_lowercase_and_missing(client):
    assert client.get("/api/protein/p24941").status_code == 200
    assert client.get("/api/protein/NOTANACC").status_code == 404
    assert client.get("/search?q=p24941", follow_redirects=False).status_code == 303


def test_fragment_numbering(client):
    """VWF (P04275) has 9 fragments; residues of fragment k map to UniProt pos = resnum + 200(k-1)."""
    d = client.get("/api/protein/P04275").json()
    assert sorted({s["fragment"] for s in d["structures"]}) == list(range(1, 10))
    p9 = [p for p in d["pockets"] if p["fragment"] == 9]
    if p9:
        r = int(p9[0]["residues"].split()[0].split(":")[1])
        assert int(p9[0]["uniprot_positions"].split()[0]) == r + 1600


def test_structure_tier1_mmcif(client):
    r = client.get("/api/structure/P24941")
    assert r.status_code == 200 and r.headers["x-atlas-tier"] == "1"
    st = _parse_cif(r.text)
    ch = st[0]["A"]
    assert len(ch) == 298  # CDK2 length
    assert ch[0].name == "MET"
    # pLDDT in B-factor, within [0, 100]
    assert 0 <= ch[0][0].b_iso <= 100
    # pocket residues of the top pocket exist in the model
    p = client.get("/api/protein/P24941").json()["pockets"][0]
    nums = {res.seqid.num for res in ch}
    assert all(int(x.split(":")[1]) in nums for x in p["residues"].split())


def test_structure_tier3_mmcif(client):
    if 3 not in client.app.state.atlas.tier_dirs():
        pytest.skip("tier-3 index absent (/tmp/pe_big/webapp_index): run scripts/webapp_build_index.py --tiers 3")
    r = client.get(f"/api/structure/{TIER3_ACC}")
    assert r.status_code == 200 and r.headers["x-atlas-tier"] == "3"
    assert r.headers["x-atlas-source"] == "foldcomp:afdb_uniprot_v4"
    st = _parse_cif(r.text)
    ch = st[0]["A"]
    s = client.get(f"/api/protein/{TIER3_ACC}").json()["structures"][0]
    assert len(ch) == s["n_res"]
    assert ch[0][0].b_iso > 0


def test_toy_orthogroup(client):
    r = client.get("/og/TOYOG1")
    assert r.status_code == 200
    t = r.text
    assert "TOY fixture" in t and "<svg" in t and "Homo sapiens" in t and "Escherichia coli" in t
    assert "#1b7837" in t and "#d73027" in t  # P and A cells drawn
    d = client.get("/api/og/TOYOG1").json()
    states = {(x["sog_id"], x["taxid"], x["state"]) for x in d["states"]}
    assert ("TOYOG1:g1", 7227, "A") in states and ("TOYOG1:g0", 83333, "G") in states
    assert client.get("/og/NOPE").status_code == 404
    assert "TOYOG1" in client.get("/og?v=TOY").text


def test_toy_case(client):
    r = client.get("/case/TOY_gain_01")
    assert r.status_code == 200
    t = r.text
    for s in ("TOY fixture case", "persistence", "Co-folding models", "toy_tree.svg", "/og/TOYOG1",
              'data-role="ancestral"', 'data-role="derived"', "TOY-EVIDENCE-1", "toy_m2.cif"):
        assert s in t
    assert client.get("/case/TOY_gain_01/file/figures/toy_tree.svg").status_code == 200
    m = client.get("/case/TOY_gain_01/file/cofold/toy_m1.cif")
    assert m.status_code == 200 and len(_parse_cif(m.text)[0]) >= 1
    # path traversal refused
    assert client.get("/case/TOY_gain_01/file/..%2F..%2Fsogs%2Fmanifest.json").status_code == 404
    assert client.get("/case/..%2Fx").status_code == 404
    assert "TOY_gain_01" in client.get("/case").text


def test_home(client):
    r = client.get("/")
    assert r.status_code == 200 and "tier 1" in r.text


def test_gene_search_with_table(tmp_path):
    """Gene-name search reads data/annotations/{proteins,uniprot_proteins}.parquet when present; toy table here."""
    import pandas as pd
    f = tmp_path / "proteins.parquet"
    pd.DataFrame({"accession": ["P24941"], "gene_primary": ["CDK2"],
                  "protein_name": ["Cyclin-dependent kinase 2"]}).to_parquet(f)
    c = TestClient(create_app(Config(annotations=[f])))
    r = c.get("/search?q=cdk2", follow_redirects=False)
    assert r.status_code == 303 and r.headers["location"] == "/protein/P24941"
    if not any((REPO / "results/webapp_index").rglob("*.parquet")):
        pytest.skip("atlas protein index not built (bulk data, not in the public repo)")
    assert "CDK2" in c.get("/protein/P24941").text


def test_case_undeclared_repo_file_refused(client):
    # only files declared in case.yaml may be served from outside the case dir
    assert client.get("/case/TOY_gain_01/file/results/cofold/x.cif").status_code == 404
    assert client.get("/case/TOY_gain_01/file/pyproject.toml").status_code == 404


def test_search_variants(client):
    r = client.get("/search?q=AF-P24941-F1-model_v6", follow_redirects=False)
    assert r.status_code == 303 and r.headers["location"] == "/protein/P24941"
    if client.app.state.atlas.gene_table() is None:
        pytest.skip("no gene table")
    j = client.get("/api/search?q=CDK2").json()
    assert j["by"] in ("gene", "gene prefix") and "P24941" in j["accessions"]
    assert client.get("/search?q=cyclin-dependent kinase 2").status_code == 200


def test_protein_pocket_param_and_ids(client):
    pid = client.get("/api/protein/P24941").json()["pockets"][0]["pocket_id"]
    t = client.get(f"/protein/P24941?pocket={pid}").text
    assert f'data-pid="{pid}"' in t and "/api/structure/P24941?tier=1&fragment=1" in t


# ---------------------------------------------------------------- real data (skip when absent)
SNAP = "20261030T085136Z"
REAL = {v: REPO / "results" / "sogs" / SNAP / v / "web" for v in ("O1_r2", "O2_r1")}


def _real(v):
    if not (REAL[v] / "sogs.parquet").exists():
        pytest.skip(f"{REAL[v]} absent")


def test_real_versions_selector(client):
    _real("O1_r2")
    ids = [v["id"] for v in client.get("/api/versions").json()]
    assert ids[0] == "O1_r2" and "TOY" in ids
    assert 'id="vsel"' in client.get("/og").text
    # selector state: ?v= wins over the cookie; unknown version falls back to the default
    assert "O1_r2" in client.get("/og?v=nonsense").text


def test_real_sog_matrix_paginated(client):
    _real("O2_r1")
    j = client.get("/api/sogs?v=O2_r1&size=7&page=2").json()
    assert j["total"] == 463014 and len(j["rows"]) == 7 and "footprint_ref" not in j["rows"][0]
    top = client.get("/api/ogs?v=O2_r1&sort=n_sogs&desc=true&size=1").json()["rows"][0]
    og = top["og_id"]
    assert top["n_sogs"] > 100
    d = client.get(f"/api/og/{og}?v=O2_r1&n=5").json()
    assert len(d["sogs"]) == 5 and d["n_sogs_total"] == top["n_sogs"]
    assert {r["sog_id"] for r in d["states"]} <= {s["sog_id"] for s in d["sogs"]}   # only the page's rows were read
    h = client.get(f"/og/{og}?v=O2_r1").text
    assert "page 1/" in h and h.count('class="cell s') > 100 and len(h) < 4_000_000
    sog = d["sogs"][0]["sog_id"]
    h2 = client.get(f"/og/{og}?v=O2_r1&sog={sog}").text
    assert f"Members of {sog}" in h2
    assert client.get("/sog?v=O2_r1&page=9000&sort=n_P").status_code == 200


def test_real_version_o1_matrix_and_fallback(client):
    _real("O1_r2")
    og = client.get("/api/ogs?v=O1_r2&size=1").json()["rows"][0]["og_id"]
    assert client.get(f"/og/{og}").status_code == 200   # default version
    assert client.get(f"/og/{og}?v=nonsense").status_code == 200


def test_real_member_links_resolve_to_atlas(client):
    """A P/A member of a real OG that is in atlas tiers 1-3 opens in the protein view with its pocket selected."""
    _real("O2_r1")
    og = client.get("/api/ogs?v=O2_r1&sort=n_species&desc=true&size=1").json()["rows"][0]["og_id"]
    d = client.get(f"/api/og/{og}?v=O2_r1&n=50").json()
    row = next((r for r in d["states"] if r["state"] == "P" and r.get("pocket_id") and r.get("model_tier") in (1, 2)), None)
    if row is None:
        pytest.skip("no tier 1-2 present member on the first SOG page")
    p = client.get(f"/api/protein/{row['accession']}").json()
    assert row["pocket_id"] in {x["pocket_id"] for x in p["pockets"]}


def test_real_cases_render(client):
    """Whatever real cases exist under cases/ (case-v1) must render, and their co-folding files must resolve."""
    import yaml
    ys = [y for y in sorted((REPO / "cases").glob("*/case.yaml"))]
    if not ys:
        pytest.skip("no cases/<id>/case.yaml yet")
    for y in ys:
        cid = y.parent.name
        r = client.get(f"/case/{cid}")
        assert r.status_code == 200, cid
        c = yaml.safe_load(y.read_text()) or {}
        for m in c.get("cofold") or []:
            if m.get("file"):
                assert client.get(f"/case/{cid}/file/{m['file']}").status_code == 200, (cid, m["file"])


def test_every_real_case_renders_page_and_api():
    """Render every cases/*/case.yaml (page + JSON API); needs no atlas index. New fields (cofold verdict, blinding,
    evidence[].filter, post_hoc label, calibration banner) must show; declared files must resolve."""
    import yaml
    ys = sorted((REPO / "cases").glob("*/case.yaml")) if (REPO / "cases").is_dir() else []
    if not ys:
        pytest.skip("cases/ absent or empty")
    cl = TestClient(create_app(Config()))
    assert cl.get("/case").status_code == 200
    v = cl.get("/validation")
    assert v.status_code == 200 and "Validation (known events)" in v.text
    assert cl.get("/api/validation").status_code == 200
    for y in ys:
        cid = y.parent.name
        c = yaml.safe_load(y.read_text()) or {}
        r = cl.get(f"/case/{cid}")
        assert r.status_code == 200, cid
        assert cl.get(f"/api/case/{cid}").status_code == 200, cid
        for m in c.get("cofold") or []:
            if not isinstance(m, dict):
                continue
            if m.get("status") and m["status"] != "pending":
                assert m["status"] in r.text, cid
            if isinstance(m.get("blinding"), dict):
                assert m["blinding"]["status"] in r.text, cid
            if m.get("file") and ((REPO / m["file"]).exists() or (y.parent / m["file"]).exists()):
                assert cl.get(f"/case/{cid}/file/{m['file']}").status_code == 200, (cid, m["file"])
        if c.get("filter_label_with_literature_post_hoc"):
            assert "post hoc" in r.text, cid
        if "fdr_supported=False" in str(c.get("notes")):
            assert "NOT FDR-calibrated" in r.text, cid
        for ev in c.get("evidence") or []:
            if ev.get("file"):
                assert cl.get(f"/case/{cid}/file/{ev['file']}").status_code == 200, (cid, ev["file"])


def test_endpoint_latency_p50(client):
    """Time a few endpoints (warm, 7 repeats) and print p50; loose bounds only. Run with -s to see the table."""
    import statistics
    import time
    urls = ["/api/protein/P24941", "/api/structure/P24941", "/protein/P24941", "/api/search?q=P24941", "/og?v=TOY"]
    if (REAL["O2_r1"] / "sogs.parquet").exists():
        top = client.get("/api/ogs?v=O2_r1&sort=n_sogs&desc=true&size=1").json()["rows"][0]["og_id"]
        urls += ["/api/sogs?v=O2_r1&page=100", "/api/sogs?v=O2_r1&sort=n_P&desc=true", "/og?v=O2_r1",
                 f"/og/{top}?v=O2_r1", f"/og/{top}?v=O2_r1&sp=3"]
    if 3 in client.app.state.atlas.tier_dirs():
        urls += [f"/api/protein/{TIER3_ACC}", f"/api/structure/{TIER3_ACC}"]
    out = {}
    for u in urls:
        client.get(u)  # warm
        ts = []
        for _ in range(7):
            t = time.perf_counter()
            assert client.get(u).status_code == 200, u
            ts.append(time.perf_counter() - t)
        out[u] = statistics.median(ts) * 1000
    for u, ms in out.items():
        print(f"p50 {ms:8.1f} ms  {u}")
    assert max(out.values()) < 3000


def test_real_protein_lists_orthogroups(client):
    """accession -> OG/SOG index (python -m webapp build-og-index): CDK2 page lists its OGs and links to the matrix."""
    idx = REPO / "results" / "webapp_index" / SNAP / "og"
    if not idx.is_dir() or not (REAL["O1_r2"] / "sogs.parquet").exists():
        pytest.skip("og index absent: run python -m webapp build-og-index")
    d = client.get("/api/protein/P24941").json()
    assert d["orthogroups"] and all(o["sogs"] for o in d["orthogroups"])
    o = d["orthogroups"][0]
    h = client.get("/protein/P24941").text
    assert "Orthogroups and site orthogroups" in h and f'/og/{o["og_id"]}?v={o["version"]}' in h


def test_build_og_index_options(tmp_path):
    """build-og-index with --snapshot/--sogs-root/--out on the toy fixture (any snapshot id, any input root)."""
    import subprocess
    import pyarrow.parquet as pq
    web = tmp_path / "sogs" / "SNAPX" / "OT" / "web"
    web.parent.mkdir(parents=True)
    web.symlink_to(REPO / "webapp" / "fixtures" / "sogs")
    out = tmp_path / "out"
    r = subprocess.run([sys.executable, "-m", "webapp", "build-og-index", "--snapshot", "SNAPX", "--sogs-root",
                        str(tmp_path / "sogs"), "--versions", "OT", "--out", str(out)], cwd=REPO, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = [x for f in (out / "OT").glob("bucket=*.parquet") for x in pq.read_table(f).to_pylist()]
    assert rows and {"accession", "og_id", "sog_k", "states"} <= set(rows[0])
