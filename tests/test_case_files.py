"""Committed case files (results/cases/<snap>/case_list_v2.parquet -> cases/<id>/case.yaml): every gain/loss/control/
specificity case lists members on both sides of the change, absent side with its state (bug report 2026-10-30)."""
from pathlib import Path

import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
LIST = ROOT / "results/cases/20261030T085136Z/case_list_v2.parquet"
SIDED = {"gain", "loss", "control", "duplication_gain_exploratory"}


def _cases():
    if not LIST.exists():
        return []
    d = pd.read_parquet(LIST)
    return [(r.case_id, r.category) for r in d[d.selected].itertuples() if r.category in SIDED]


@pytest.mark.parametrize("case_id,category", _cases())
def test_case_lists_both_sides(case_id, category):
    c = yaml.safe_load((ROOT / "cases" / case_id / "case.yaml").read_text())
    prot = c.get("proteins") or []
    present = [p for p in prot if p.get("side") == "present"]
    absent = [p for p in prot if p.get("side") == "absent"]
    assert present, f"{case_id}: no present-side member"
    assert absent, f"{case_id}: no absent-side member"
    assert all(p.get("state") for p in absent)
    # the absent side carries the absence: A, W (after rescue) or U, never P
    assert all(p["state"] in ("A", "W", "U") for p in absent)
    if category == "control":
        # false-loss controls: the members that drove the original loss call are listed with their pre-filter state
        assert any(str(p.get("state_pre_filter", "")).startswith("A") for p in absent), case_id
