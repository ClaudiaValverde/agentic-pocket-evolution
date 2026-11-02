"""Paths and settings of the research browser. Every path can be overridden by an environment variable
(PE_WEBAPP_*) or by the command-line flags of `python -m webapp`."""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WEBAPP = Path(__file__).resolve().parent
DEFAULT_SNAPSHOT = "20261030T085136Z"
# orthology versions looked for under results/sogs/<snapshot>/<v>/web (first existing = default; absent ones are skipped)
DEFAULT_ORTHOLOGY_VERSIONS = ("O1_r2", "O2_r1")


def _env(name, default):
    return os.environ.get(f"PE_WEBAPP_{name}", default)


@dataclass
class Config:
    snapshot: str = field(default_factory=lambda: _env("SNAPSHOT", DEFAULT_SNAPSHOT))
    atlas_root: Path = None          # from the snapshot manifest
    afdb_v6_dir: Path = field(default_factory=lambda: Path(_env("AFDB_V6", "/cluster/work/beltrao/shared/alphafold/v6")))
    foldcomp_dir: Path = None        # default <atlas_root>/foldcomp
    index_dirs: list = None          # dirs holding tier<T>/bucket=XX.parquet
    annotations: list = None         # candidate gene-name tables
    sog_dirs: list = None            # orthogroup data dirs (results/sogs/<snapshot>, fixtures)
    case_dirs: list = None           # case roots (cases/, fixtures/cases)
    structure_cache: int = 64        # in-memory LRU entries (mmCIF text)

    def __post_init__(self):
        snap = REPO / "data" / "snapshot" / self.snapshot
        if self.atlas_root is None:
            root = _env("ATLAS", None)
            if root is None and (snap / "manifest.json").exists():
                root = json.loads((snap / "manifest.json").read_text())["atlas_root"]
            self.atlas_root = Path(root or "/cluster/scratch/gmagnusson/pocketcompare_afdb")
        if self.foldcomp_dir is None:
            self.foldcomp_dir = Path(_env("FOLDCOMP", str(self.atlas_root / "foldcomp")))
        if self.index_dirs is None:
            env = _env("INDEX_DIRS", None)
            self.index_dirs = [Path(p) for p in env.split(":")] if env else [
                REPO / "results" / "webapp_index" / self.snapshot,
                Path("/tmp/pe_big/webapp_index") / self.snapshot]
        if self.annotations is None:
            self.annotations = [REPO / "data" / "annotations" / "proteins.parquet",
                                REPO / "data" / "annotations" / "uniprot_proteins.parquet"]
        if self.sog_dirs is None:
            env = _env("SOG_DIRS", None)
            self.sog_dirs = [Path(p) for p in env.split(":")] if env else [
                REPO / "results" / "sogs" / self.snapshot / v / "web"
                for v in _env("SOG_VERSIONS", ",".join(DEFAULT_ORTHOLOGY_VERSIONS)).split(",") if v] + [
                WEBAPP / "fixtures" / "sogs"]
        if self.case_dirs is None:
            env = _env("CASE_DIRS", None)
            self.case_dirs = [Path(p) for p in env.split(":")] if env else [
                REPO / "cases", WEBAPP / "fixtures" / "cases"]
