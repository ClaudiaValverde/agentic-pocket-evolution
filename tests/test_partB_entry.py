"""The Part B entry point with default options reproduces the committed O1_r3 candidates by content.

Runs scripts/partB/run_partB.sh (stages post, candidates, artifacts) into a temporary results root, starting from
the committed rescued states, rescue hits and null funnel, and compares the gain/loss candidate table
(identity, rule outcome, FDR label, artifact verdicts) with results/partB/<snap>/O1_r3/candidates.parquet.
Needs the node-local SOG bulk (/tmp/pe_big/sogs/<snap>/O1_r3) for the artifact checks; skipped without it.
Runtime ~2-4 min.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.default_snapshot()
RUN = "O1_r3"
COMMITTED = ROOT / f"results/partB/{SNAP}/{RUN}/candidates.parquet"
KEY = ["primary", "tau", "gain_cost", "loss_cost", "type", "og_id", "sog_id", "node"]
VALS = ["n_species_P", "n_species_A", "out_A_lineages", "out_P_lineages", "fdr_supported", "pass_artifacts",
        "artifact_fails", "flag_paralog_branch", "a_species_tree_lineages", "a_fragment_present", "orthology_split"]


@pytest.mark.skipif(not Path(f"/tmp/pe_big/sogs/{SNAP}/{RUN}/msa").exists() or not COMMITTED.exists(),
                    reason="SOG bulk or committed candidates missing")
def test_entry_point_reproduces_o1_r3_candidates(tmp_path):
    rr = tmp_path / "partB"
    src = ROOT / "results/partB"
    for t in ("0.2", "0.3"):
        d = f"events/{SNAP}/{RUN}/tau{t}/primary"
        (rr / d / "rescue").mkdir(parents=True)
        for f in ("funnel.csv", "manifest.json", "events.parquet"):
            shutil.copy(src / d / f, rr / d / f)
        for f in ("rescue_hits.parquet", "rescue_summary.json"):
            shutil.copy(src / d / "rescue" / f, rr / d / "rescue" / f)
        shutil.copytree(src / d / "rescue/sog_states", rr / d / "rescue/sog_states")
    for d in (f"nulls/{SNAP}/{RUN}/tau0.2/primary", f"m3/{SNAP}/{RUN}/tau0.2"):
        shutil.copytree(src / d, rr / d)
    (rr / f"{SNAP}/{RUN}").mkdir(parents=True)
    for f in ("false_loss_controls.parquet", "rbh_outgroups.parquet"):     # inputs from earlier stages
        if (src / f"{SNAP}/{RUN}/{f}").exists():
            shutil.copy(src / f"{SNAP}/{RUN}/{f}", rr / f"{SNAP}/{RUN}/{f}")
    if (src / f"{SNAP}/{RUN}/specificity").exists():
        shutil.copytree(src / f"{SNAP}/{RUN}/specificity", rr / f"{SNAP}/{RUN}/specificity")
    r = subprocess.run([str(ROOT / "scripts/partB/run_partB.sh"), "--orthology", "O1", "--stages",
                        "post,candidates,artifacts", "--results-root", str(rr), "--workers", "4"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    new = pd.read_parquet(rr / f"{SNAP}/{RUN}/candidates.parquet")
    old = pd.read_parquet(COMMITTED)
    a = new[KEY + VALS].sort_values(KEY).reset_index(drop=True)
    b = old[KEY + VALS].sort_values(KEY).reset_index(drop=True)
    assert len(a) == len(b)
    pd.testing.assert_frame_equal(a, b, check_dtype=False)
