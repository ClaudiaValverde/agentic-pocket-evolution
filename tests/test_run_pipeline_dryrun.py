"""G5 re-runnability: `run_pipeline.py --dry-run` on a snapshot id that has never existed must find every input.

Failure mode this guards (G5, snapshot 20261031T020009Z): entry points looked up frozen / species-level inputs under
the RUNNING snapshot (results/observation/<snap>/rescue_k.json, data/timetree/<snap>_O1_pairwise.csv), and the O2
lane read O1 events of a run that only snapshot 1 had. Also: no G5 stage may write results/numbers.json["partA"].
"""
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PY = str(REPO / ".venv/bin/python")
FAKE = "29991231T235959Z"


def dry(*extra, env=None):
    return subprocess.run([PY, "scripts/run_pipeline.py", "--snapshot", FAKE, "--run-id", "drytest", "--dry-run",
                           *extra], cwd=REPO, capture_output=True, text=True, timeout=600,
                          env={**os.environ, **(env or {})})


@pytest.fixture(scope="module")
def out():
    if not Path("/tmp/pe_big").exists():
        pytest.skip("bulk root /tmp/pe_big not present (cluster-only test)")
    r = dry()
    return r


def test_dry_run_finds_every_input(out):
    assert out.returncode == 0, out.stdout[-3000:] + out.stderr[-3000:]
    assert "# input check: all inputs found" in out.stdout
    assert "MISSING INPUT" not in out.stdout
    assert "BLOCKED" not in out.stdout, "a stage is blocked on a new snapshot"


def test_input_check_detects_missing_frozen_input():
    if not Path("/tmp/pe_big").exists():
        pytest.skip("cluster-only")
    r = dry("--stages", "partB_O1", env={"PE_FROZEN_RESCUE_K": "/nonexistent/rescue_k.json"})
    assert r.returncode != 0
    assert "MISSING INPUT [partB_O1] /nonexistent/rescue_k.json" in r.stdout


def test_frozen_inputs_do_not_follow_the_running_snapshot():
    sys.path.insert(0, str(REPO / "src"))
    from pocketevo import frozen
    a = {k: str(v) for k, v in frozen.all_paths().items()}
    os.environ["PE_SNAPSHOT"] = FAKE
    try:
        b = {k: str(v) for k, v in frozen.all_paths().items()}
    finally:
        del os.environ["PE_SNAPSHOT"]
    assert a == b and not any(FAKE in v for v in b.values())
    if any(v.endswith(".parquet") and not Path(v).exists() for v in a.values()):
        pytest.skip("frozen parquet inputs are bulk data, not in the public repo")
    assert all(Path(v).exists() for v in a.values()), a


def test_lanes_and_numbers_are_snapshot_safe(out):
    lines = out.stdout.splitlines()
    o2 = [x for x in lines if x.startswith("[partB_O2]")]
    assert o2 and all("agreement" not in x for x in o2), "the O2 lane must not read O1 events"
    agr = [x for x in lines if x.startswith("[partB_O1O2]")]
    assert agr and all("--o1-run O1_drytest" in x for x in agr), "agreement must use this snapshot's O1 run"
    num = [x for x in lines if "scripts/partA/11_numbers.py" in x]
    assert num and all("--numbers-json none" in x for x in num), "G5 must not write results/numbers.json[partA]"
    d2 = [x for x in lines if "scripts/partA/16_descent_2d.py" in x]
    assert d2 and all("--timetree" in x and "--datasets" in x for x in d2)
    o1 = [x for x in lines if x.startswith("[partB_O1]")]
    assert o1 and all("rbh,artifacts" in x for x in o1), "the artifact cascade needs the RBH check (D28) first"
    cases = [x for x in lines if x.startswith("[cases]")]
    assert cases and all(f"--cases-dir /tmp/pe_big/cases_files/{FAKE}/drytest" in x for x in cases)


# entry points the driver runs with internal defaults: no code line may key a frozen input on the running snapshot
PIPELINE_SCRIPTS = ["scripts/partB/run_rescue.py", "scripts/partB/run_m3.py", "scripts/partA/16_descent_2d.py",
                    "scripts/partB/run_events.py", "scripts/partB/run_nulls.py", "scripts/partB/write_candidates.py",
                    "scripts/partB/run_artifacts.py", "scripts/partB/prep_o2_states.py",
                    "scripts/partB/o1_o2_agreement.py", "scripts/partB/write_o2_report.py"]
BAD = re.compile(r"(results/observation/\{(SNAP|snap|a\.snapshot)\}|data/timetree/\{(SNAP|snap|a\.snapshot)\})")


@pytest.mark.parametrize("script", PIPELINE_SCRIPTS)
def test_no_snapshot_keyed_frozen_inputs(script):
    code = [ln for ln in (REPO / script).read_text().splitlines() if not ln.lstrip().startswith("#")]
    hits = [ln for ln in code if BAD.search(ln)]
    assert not hits, f"{script}: frozen input keyed on the running snapshot: {hits}"


def test_partA_numbers_guard_present():
    s = (REPO / "scripts/partA/11_numbers.py").read_text()
    assert "--numbers-json" in s and "refusing to overwrite results/numbers.json['partA']" in s
