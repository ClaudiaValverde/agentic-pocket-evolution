"""Rescue of absence (A) calls by alternative BURROW settings (DESIGN §3.3; prereg states.rescue).

Six settings (prereg states.rescue.settings, in order): "--min-score none", "no AF mode", "seed s1", "seed s2",
"preset fast", "preset accurate". Each is ONE change from the base run, which mirrors the atlas detector
(`burrow predict --af-mode --output-rule region`: AF mode, region output rule, pinned default model s0):

    min_score_none   --af-mode --min-score none     (full NMS'd list, no score floor; burrow refuses
                                                    --output-rule region together with --min-score)
    no_af_mode       base without --af-mode         (no pLDDT pruning / soft occlusion)
    seed_s1          base + --model burrow_all_default_s1.json
    seed_s2          base + --model burrow_all_default_s2.json
    preset_fast      base + --preset fast
    preset_accurate  base + --preset accurate

No pLDDT or score filter is applied to rescue pockets (the point is to find anything there). A setting "hits"
a footprint if one of its pockets has residue Jaccard >= tau with the footprint (rule "jaccard", the P-call rule)
or covers >= min_cov of the footprint residues (rule "coverage", sensitivity). A is turned into W if the number
of hitting settings >= k (k from results/observation/<snapshot>/rescue_k.json).

Usage
    from pocketevo.rescue import rescue
    hits = rescue("AF-P04035-F1-model_v6.cif.gz", footprint={123, 124, 127}, workdir="/tmp/x")
    # -> {"min_score_none": {"hit": True, "jaccard": 0.41, "coverage": 0.9, "n_pockets": 23}, ...}

Batch use (many models, one `burrow predict` per setting): run_settings(models, outdir, workers) then
hits_from_outputs(outdir, stem, footprint).
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

P = Path(__file__).resolve().parents[2]
MODELS = P / "external/burrow/models"
BURROW = str(Path(sys.executable).with_name("burrow"))
BASE = ["--af-mode", "--output-rule", "region"]
SETTINGS: dict[str, list[str]] = {
    "min_score_none": ["--af-mode", "--min-score", "none"],   # burrow forbids region + min-score; full list
                                                              # (<= 50 NMS'd pockets) is a superset of region
    "no_af_mode": ["--output-rule", "region"],
    "seed_s1": BASE + ["--model", str(MODELS / "burrow_all_default_s1.json")],
    "seed_s2": BASE + ["--model", str(MODELS / "burrow_all_default_s2.json")],
    "preset_fast": BASE + ["--preset", "fast"],
    "preset_accurate": BASE + ["--preset", "accurate"],
}
SETTING_NAMES = list(SETTINGS)


def model_stem(path) -> str:
    n = Path(path).name
    for suf in (".gz", ".cif", ".pdb", ".ent"):
        n = n.removesuffix(suf)
    return n


def run_settings(models: list, outdir, workers: int = 4, settings: dict | None = None, backend: str = "auto",
                 env: dict | None = None) -> dict:
    """Run `burrow predict` once per setting over all models -> {setting: output dir}. Existing outputs (one CSV
    per model stem) are reused, so the call is resumable."""
    settings = settings or SETTINGS
    outdir = Path(outdir)
    e = dict(os.environ if env is None else env)
    for k in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        e.setdefault(k, "1")
    out = {}
    for name, args in settings.items():
        d = outdir / name
        d.mkdir(parents=True, exist_ok=True)
        todo = [str(m) for m in models if not (d / f"{model_stem(m)}.csv").exists()]
        if todo:
            # long argument lists: burrow accepts many inputs; chunk to stay under ARG_MAX
            for i in range(0, len(todo), 500):
                cmd = [BURROW, "predict", *todo[i:i + 500], "--out", str(d), "--workers", str(workers),
                       "--backend", backend, "--no-points", *args]
                r = subprocess.run(cmd, env=e, capture_output=True, text=True)
                if r.returncode != 0:
                    raise RuntimeError(f"burrow predict ({name}) failed rc={r.returncode}: {r.stderr[-1500:]}")
        out[name] = d
    return out


def read_burrow_csv(path, offset: int = 0) -> list[tuple[float, frozenset]]:
    """[(score, residue set)] from a burrow predict CSV (rank,score,center_x,y,z,residues)."""
    out = []
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    with open(p) as fh:
        for row in csv.DictReader(fh):
            res = frozenset(int("".join(c for c in tok.split(":")[-1] if c.isdigit() or c == "-")) + offset
                            for tok in (row.get("residues") or "").split() if tok)
            out.append((float(row["score"]), res))
    return out


def setting_hit(pockets: list, footprint, tau: float = 0.3, rule: str = "jaccard", min_cov: float = 0.5) -> dict:
    fp = frozenset(int(x) for x in footprint)
    bj = bc = 0.0
    hit = False
    for _, res in pockets:
        inter = len(fp & res)
        if not inter:
            continue
        j = inter / len(fp | res)
        c = inter / len(fp)
        bj, bc = max(bj, j), max(bc, c)
        if (rule == "jaccard" and j >= tau) or (rule == "coverage" and c >= min_cov) or \
                (rule == "either" and (j >= tau or c >= min_cov)):
            hit = True
    return {"hit": hit, "jaccard": bj, "coverage": bc, "n_pockets": len(pockets)}


def hits_from_outputs(outdirs: dict, stem: str, footprint, tau: float = 0.3, rule: str = "jaccard",
                      min_cov: float = 0.5, offset: int = 0) -> dict:
    """Per-setting hit records for one model stem from run_settings() output dirs."""
    return {name: setting_hit(read_burrow_csv(Path(d) / f"{stem}.csv", offset), footprint, tau, rule, min_cov)
            for name, d in outdirs.items()}


def n_hits(hits: dict) -> int:
    return int(sum(bool(h["hit"]) for h in hits.values()))


def rescued(hits: dict, k: int) -> bool:
    """A -> W if >= k of the settings hit the footprint."""
    return n_hits(hits) >= int(k)


def rescue(model, footprint, workdir=None, tau: float = 0.3, rule: str = "jaccard", min_cov: float = 0.5,
           workers: int = 1, offset: int = 0) -> dict:
    """Single model convenience: run the six settings on `model` and return per-setting hits for `footprint`
    (UniProt residue numbers; `offset` is added to model residue numbers, 0 for fragment-1 models)."""
    tmp = None
    if workdir is None:
        tmp = tempfile.TemporaryDirectory(prefix="rescue_")
        workdir = tmp.name
    try:
        outs = run_settings([model], workdir, workers=workers)
        return hits_from_outputs(outs, model_stem(model), footprint, tau, rule, min_cov, offset)
    finally:
        if tmp is not None:
            tmp.cleanup()
