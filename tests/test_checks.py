import importlib.util, json, sys
from pathlib import Path
import pytest

S = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(S))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, S / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


cp, crn, cprov = _load("check_prereg"), _load("check_report_numbers"), _load("check_provenance")
import provlib as L

DESIGN = "# D\n### Deviation ledger\n| id | date |\n|---|---|\n| D7 | x |\n| C1 | y |\n"


def mk(tmp, prereg_text="a: 1\n", sidecar_hash=None, doc="no refs\n"):
    (tmp / "config").mkdir()
    (tmp / "config/prereg.yaml").write_text(prereg_text)
    h = L.sha256_file(tmp / "config/prereg.yaml")
    (tmp / "record").mkdir()
    (tmp / "record/DESIGN.md").write_text(DESIGN)
    (tmp / "results/s").mkdir(parents=True)
    (tmp / "cases").mkdir()
    (tmp / "results/s/x.parquet.json").write_text(json.dumps({"prereg_sha256": sidecar_hash or h, "seed": 0}))
    (tmp / "results/s/DOC.md").write_text(doc)
    return h


def test_good_prereg(tmp_path):
    h = mk(tmp_path)
    r = cp.check(tmp_path, frozen=h)
    assert r["prereg_ok"] and not r["fail"]


def test_bad_prereg_hash_and_sidecar(tmp_path):
    h = mk(tmp_path, sidecar_hash="0" * 64)
    assert cp.check(tmp_path, frozen="f" * 64)["prereg_ok"] is False
    r = cp.check(tmp_path, frozen=h)
    assert r["sidecars"]["differ"] and r["fail"]


def test_missing_sidecar_hash(tmp_path):
    h = mk(tmp_path)
    (tmp_path / "results/s/y.csv.json").write_text(json.dumps({"seed": 1}))
    r = cp.check(tmp_path, frozen=h)
    assert r["sidecars"]["missing"] == ["results/s/y.csv.json"]


def test_prereg_variant_field_names():
    assert L.prereg_of({"meta": {"prereg_sha256": "SHA256:" + "A" * 64}}) == "a" * 64
    assert L.commit_of({"code_commit": "abc"}) == "abc"


def test_ledger_ids(tmp_path):
    h = mk(tmp_path, doc="see D7 and D99 and `D5` and\n```\nD6\n```\n")
    r = cp.check(tmp_path, frozen=h)
    assert set(r["unknown_ids"]) == {"D99"}


def test_ledger_ok(tmp_path):
    h = mk(tmp_path, doc="Deviation D7 applies.\n")
    assert not cp.check(tmp_path, frozen=h)["unknown_ids"]


def _report(tmp, text, numbers):
    (tmp / "results").mkdir(exist_ok=True)
    (tmp / "results/numbers.json").write_text(json.dumps(numbers))
    (tmp / "REPORT.md").write_text(text)
    return crn.check(tmp / "REPORT.md", tmp)


def test_numbers_match_with_rounding(tmp_path):
    n, un, _ = _report(tmp_path, "Recall was 0.65 and F1 72% with 3,221,702 sites.\n",
                       {"a": {"recall": 0.6529, "f1": 0.7179, "n": 3221702}})
    assert n == 3 and un == []


def test_numbers_unmatched(tmp_path):
    n, un, _ = _report(tmp_path, "Recall was 0.75.\n\nTable 99.5 here.\n", {"recall": 0.6529})
    assert [u[1] for u in un] == ["0.75", "99.5"] and un[0][0] == 1


def test_numbers_skips(tmp_path):
    txt = ("## 4.2 Heading\nDate 2026-10-30 10:47, snapshot 20261030T085136Z, doi 10.1038/s41586-021-03819-2, "
           "PDB 1abc, §5.4, D12, AF3, in 2024, 7 sites.\n```\nx = 12345.678\n```\nuse `0.123` inline\n")
    n, un, _ = _report(tmp_path, txt, {})
    assert un == []


def test_linked_sidecar_and_allowlist(tmp_path):
    (tmp_path / "results/x").mkdir(parents=True)
    (tmp_path / "results/x/t.csv").write_text("a,b\n1,0.4321\n")
    (tmp_path / "config").mkdir()
    (tmp_path / "config/report_numbers_allowlist.txt").write_text("94.7  # Smith 2020\n")
    n, un, _ = _report(tmp_path, "Value 0.43 ([t](results/x/t.csv)) and literature 94.7%.\n\nElsewhere 0.43.\n", {})
    assert not any(u[1] == "0.43" and u[0] == 1 for u in un)
    # line 1 matched via link, 94.7 via allowlist; line 3 is outside the +-3 window? it is within, so also matched
    assert all(u[1] != "94.7%" for u in un)


def test_manifest_verify_and_missing(tmp_path):
    d = tmp_path / "results/s"
    d.mkdir(parents=True)
    (d / "a.txt").write_text("hi")
    (d / "b.txt").write_text("bye")
    ha, hb = L.sha256_file(d / "a.txt"), L.sha256_file(d / "b.txt")
    (d / "BULK_MANIFEST.sha256").write_text(f"{ha}  a.txt\n{'0'*64}  b.txt\n{hb}  gone.txt\n")
    v = cprov.verify_manifest(d / "BULK_MANIFEST.sha256", tmp_path / "results")
    assert v["ok"] == 1 and v["bad"] == ["b.txt"] and v["missing"] == ["gone.txt"]
