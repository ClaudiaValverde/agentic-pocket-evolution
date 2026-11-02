"""Case records, schema `case-v1` (webapp/SCHEMA.md §3)."""
from __future__ import annotations

import re
from pathlib import Path

import yaml

from .config import REPO

CASE_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
REQUIRED = ("schema", "case_id", "title", "category", "label")
KNOWN = set(REQUIRED) | {"track", "toy", "summary", "og_id", "sog_id", "score", "proteins", "filters",
                         "artifact_checks", "evidence", "cofold", "figures", "literature", "notes", "headline",
                         "filter_label_with_literature_post_hoc", "provenance"}


def fdr_status(c: dict) -> dict:
    """Calibration status of a case, derived from category/notes ("fdr_supported=True/False", "selection: ...").
    Returns {"kind": fdr|not_fdr|validation|unknown, "text": str, "selection": str|None, "exploratory": bool}."""
    notes = str(c.get("notes") or "")
    m = re.search(r"fdr_supported=(True|False)", notes)
    sel = re.search(r"selection: ([A-Za-z_]+)", notes)
    sel = sel.group(1) if sel else None
    expl = "exploratory" in str(c.get("category") or "").lower() or "exploratory" in (sel or "").lower()
    if str(c.get("category") or "").startswith("validation"):
        return {"kind": "validation", "text": "known-answer validation, not a discovery and not FDR-calibrated",
                "selection": sel, "exploratory": False}
    if m and m.group(1) == "True" and not expl:
        return {"kind": "fdr", "text": "FDR-calibrated candidate", "selection": sel, "exploratory": False}
    if m or expl:
        return {"kind": "not_fdr", "text": "NOT FDR-calibrated" + (" - exploratory, not preregistered" if expl else ""),
                "selection": sel, "exploratory": expl}
    return {"kind": "unknown", "text": "", "selection": sel, "exploratory": False}


class Cases:
    def __init__(self, dirs):
        self.dirs = [Path(d) for d in dirs]

    def _dir(self, case_id: str) -> Path | None:
        if not CASE_ID_RE.match(case_id):
            return None
        for root in self.dirs:
            d = root / case_id
            if (d / "case.yaml").exists():
                return d
        return None

    def list(self) -> list[dict]:
        out, seen = [], set()
        for root in self.dirs:
            if not root.is_dir():
                continue
            for y in sorted(root.glob("*/case.yaml")):
                cid = y.parent.name
                if cid in seen or not CASE_ID_RE.match(cid):
                    continue
                seen.add(cid)
                try:
                    c = yaml.safe_load(y.read_text()) or {}
                except yaml.YAMLError as e:
                    c = {"title": f"(unparseable case.yaml: {e})"}
                out.append({"case_id": cid, "title": c.get("title"), "category": c.get("category"),
                            "label": c.get("label"), "toy": bool(c.get("toy")),
                            "fdr": fdr_status(c) if isinstance(c, dict) else {"kind": "unknown", "text": ""}, "dir": str(y.parent)})
        return out

    def get(self, case_id: str) -> dict | None:
        d = self._dir(case_id)
        if d is None:
            return None
        try:
            c = yaml.safe_load((d / "case.yaml").read_text()) or {}
        except yaml.YAMLError as ex:
            c = {"title": f"(unparseable case.yaml: {ex})"}
        if not isinstance(c, dict):
            c = {"title": "(case.yaml is not a mapping)"}
        c["_missing_required"] = [k for k in REQUIRED if k not in c]
        c["_extra"] = {k: v for k, v in c.items() if k not in KNOWN and not k.startswith("_")}
        c["_dir"] = str(d)
        c["_fdr"] = fdr_status(c)
        return c

    def file(self, case_id: str, rel: str) -> Path | None:
        """Resolve a file of a case. Normally inside the case dir (paths that leave it are refused). Exception
        (SCHEMA.md §3): a `file` declared in case.yaml under cofold/figures/evidence may be repo-relative and live
        under results/ or cases/ (e.g. `results/cofold/<run>/model.cif`); only declared files are served."""
        d = self._dir(case_id)
        if d is None:
            return None
        p = (d / rel).resolve()
        if d.resolve() in p.parents and p.is_file():
            return p
        try:
            c = yaml.safe_load((d / "case.yaml").read_text()) or {}
        except yaml.YAMLError:
            return None
        declared = {x.get("file") for k in ("cofold", "figures", "evidence") for x in (c.get(k) or [])
                    if isinstance(x, dict)}
        if rel not in declared:
            return None
        # repo-relative: this repo, or the repo that holds the case dir (<root>/cases/<id>, e.g. a full workspace)
        for root in dict.fromkeys((REPO, d.parent.parent)):
            p = (root / rel).resolve()
            if any(r.resolve() in p.parents for r in (root / "results", root / "cases")) and p.is_file():
                return p
        return None
