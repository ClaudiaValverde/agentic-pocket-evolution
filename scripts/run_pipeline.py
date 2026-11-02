#!/usr/bin/env python
"""One-command rerun of the whole pipeline on an atlas snapshot (G5): snapshot -> sites -> web index ->
SOGs O1 + O2 -> rescue pool -> Part A -> Part B -> cases -> OG web index.

Every stage calls the owning agent's documented entry point unchanged (commands taken from their READMEs, docs and
commit messages; see STAGES below). Every stage takes the snapshot id and the run id, and refuses to write into an
output that is tracked by git (a committed run) or that already exists (unless --resume).
Stages whose entry point cannot take a new snapshot / run id are listed as BLOCKED with the reason and are not
patched here (see `blockers()`; printed by --dry-run).

    .venv/bin/python scripts/run_pipeline.py --snapshot new --run-id g5 --dry-run           # print commands
    .venv/bin/python scripts/run_pipeline.py --snapshot 20261030T085136Z --run-id pt1 \\
        --families P24941,P00918,P03372,P00338,P10275 --stages sites,sogs_O1,sogs_O2,rescue_O2,partB_O1_events \\
        --bulk-root /tmp/pe_big/pipeline_test --results-root /tmp/pe_big/pipeline_test/results
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "src"))
from pocketevo import frozen  # noqa: E402  (snapshot-invariant inputs: frozen decisions, species-level data)
PY = str(P / ".venv/bin/python")
FAISS_PY = "/tmp/pe_big/envs/faiss/bin/python"           # Part A env (results/partA/.../PARTA.md)
ATLAS = "/cluster/scratch/gmagnusson/pocketcompare_afdb"
CUR_SNAP = "20261030T085136Z"                              # snapshot hard-coded in Part B / cases scripts
TIERL_MANIFEST = "data/tierL/tierL_20261030T112000Z/manifest.parquet"
O1_ORTH = "data/orthology/ancestra_tier1_v0/accurate"
O2_ORTH = "data/orthology/ancestra_O2/accurate"
TAU = "0.2"                                                # D8: frozen
STAGE_ORDER = ["snapshot", "sites", "webindex", "sogs_O1", "sogs_O2", "bench_O1", "bench_O2", "split_O2", "rescue_O2", "partA",
               "partB_O1_events", "partB_O1", "partB_O2", "partB_O1O2", "cases", "og_index"]
# Dependencies between stages (inputs a stage reads from another stage's outputs). partB_O1O2 (O1 -> O2 gain
# agreement + O2 report) needs BOTH Part B lanes of the same snapshot; cases needs Part B O1 and Part A.
DEPENDS = {"webindex": ["sites"], "sogs_O1": ["sites"], "sogs_O2": ["sites"], "bench_O1": ["sogs_O1"],
           "bench_O2": ["sogs_O2"], "split_O2": ["sogs_O1", "sogs_O2"], "rescue_O2": ["sogs_O2"],
           "partA": ["sites", "sogs_O1"], "partB_O1_events": ["sogs_O1"], "partB_O1": ["sogs_O1", "webindex"],
           "partB_O2": ["rescue_O2"], "partB_O1O2": ["partB_O1", "partB_O2"], "cases": ["partB_O1", "partA"],
           "og_index": ["sogs_O1", "sogs_O2"]}
OUTPUT_FLAGS = {"--out", "--out12", "--out3", "--results", "--work", "--cases-dir", "--log"}


def git_tracked(path: Path) -> bool:
    r = subprocess.run(["git", "-C", str(P), "ls-files", "--", str(path)], capture_output=True, text=True)
    return bool(r.stdout.strip())


def help_text(cmd: list) -> str:
    try:
        r = subprocess.run(cmd + ["--help"], cwd=P, capture_output=True, text=True, timeout=120)
        return r.stdout + r.stderr
    except Exception:
        return ""


PARTB_SH = P / "scripts/partB/run_partB.sh"     # announced by the Part B agent (D27): --snapshot / --orthology


def partb_driver_ok() -> tuple[bool, str]:
    """run_partB.sh is used only if it exists and its --help lists the options used below."""
    if not PARTB_SH.exists():
        return False, "scripts/partB/run_partB.sh not present yet"
    h = help_text(["bash", str(PARTB_SH)])
    miss = [o for o in ("--snapshot", "--orthology", "--sog-run", "--run-id", "--stages") if o not in h]
    return (not miss), ("run_partB.sh --help lacks " + ", ".join(miss)) if miss else "ok"


def blockers(a) -> dict:
    """Stages that cannot be rerun unchanged on another snapshot / run id (not patched; reported to the PI)."""
    b = {}
    pb_ok, pb_why = partb_driver_ok()
    cases_snap = "--snapshot" in help_text([PY, "scripts/cases/build_cases.py"])
    if a.snapshot != CUR_SNAP and not pb_ok:
        b["partB_O1"] = (pb_why + "; scripts/partB/run_o1_pipeline.sh hard-codes S=20261030T085136Z; run_artifacts.py, "
                         "false_loss_controls.py, run_specificity.py, run_m3.py, write_candidates.py, run_rescue.py, "
                         "apply_d11.py, make_noise_json.py define SNAP = '20261030T085136Z' (no --snapshot option)")
    if a.snapshot != CUR_SNAP and not cases_snap:
        b["cases"] = ("scripts/cases/build_cases.py: SNAP constant in all input defaults (no --snapshot); works only "
                      "if every input path is passed explicitly")
    if not pb_ok:
        b["partB_O2"] = "no O2 driver for Part B yet (" + pb_why + ")"
    if a.results_root != str(P / "results"):
        b.setdefault("partB_O1", "Part B scripts write to fixed results/partB/... paths (no output-root option); "
                                 "only run_events.py (--out) can be redirected (stage partB_O1_events)")
        b.setdefault("cases", "build_cases.py writes cases/<id>/ in the repo unless --cases-dir/--out are given")
        if "--input-root" not in help_text([PY, "scripts/webapp_build_og_index.py"]):
            b["og_index"] = "webapp_build_og_index.py reads inputs from results/sogs/<snap>/<version>/web only (no input root)"
    return b


def stages(a) -> dict:
    """name -> list of (command, [outputs that must not be committed/existing])."""
    S, R = a.snapshot, a.run_id
    res, bulk = Path(a.results_root), Path(a.bulk_root)
    man = f"{a.snapshot_root}/{S}/manifest.parquet"
    sites = bulk / "sites" / S
    o1, o2 = f"O1_{R}", f"O2_{R}"
    sl = ["--slice", a.families] if a.families else []
    common = ["--manifest", man, "--prereg", "config/prereg.yaml", "--out", str(bulk / "sogs"), "--workers", str(a.workers),
              "--sites", str(sites)]
    o1a = [PY, "scripts/build_sogs.py", *common, "--orthology", O1_ORTH, "--orthology-version", "O1", "--run-id", o1]
    o2a = [PY, "scripts/build_sogs.py", *common, "--orthology", O2_ORTH, "--orthology-version", "O2", "--run-id", o2,
           "--model-tiers", "1,2,3,4", "--idmap", "data/orthology/ancestra_O2/idgap_mapping.parquet",
           "--proteome-fasta", "/tmp/pe_big/ancestra/o2/input_O2", "--tierl-manifest", TIERL_MANIFEST]
    r1, r2 = res / "sogs" / S / o1, res / "sogs" / S / o2
    # MSA cache seeding (MSAs depend only on member sequences; build_sogs re-checks every cached MSA against them)
    seed1 = [] if not a.msa_seed else [(["bash", "-c", f"mkdir -p {bulk}/sogs/{S}/{o1} && cp -rn "
                                          f"{a.msa_seed}/{CUR_SNAP}/O1_r2/msa {bulk}/sogs/{S}/{o1}/"], [])]
    seed2 = [] if not a.msa_seed else [(["bash", "-c", f"mkdir -p {bulk}/sogs/{S}/{o2} && cp -rn "
                                          f"{a.msa_seed}/{CUR_SNAP}/O2_r1/msa {bulk}/sogs/{S}/{o2}/"], [])]
    st = {}
    st["snapshot"] = [([PY, "scripts/atlas_snapshot.py", "--atlas-root", ATLAS, "--out", a.snapshot_root,
                        "--threads", "32"], [])]
    st["sites"] = [([PY, "scripts/build_sites.py", "--manifest", man, "--tiers", "1,2", "--out", str(sites),
                     "--workers", str(a.workers)] + (["--no-fingerprints"] if a.families else []), []),
                   ([PY, "scripts/build_sites.py", "--manifest", man, "--tiers", "3", "--out", str(sites),
                     "--workers", str(a.workers), "--no-fingerprints"], []),
                   ([PY, "scripts/build_sites.py", "--manifest", TIERL_MANIFEST, "--tiers", "4", "--out", str(sites),
                     "--workers", str(a.workers)], [])]
    st["webindex"] = [([PY, "-I", "scripts/webapp_build_index.py", "--snapshot", S, "--tiers", "1,2,3",
                        "--out12", str(res / "webapp_index" / S), "--out3", str(bulk / "webapp_index" / S)],
                       [res / "webapp_index" / S / "tier1"])]
    st["sogs_O1"] = seed1 + [(o1a + ["--stage", "models,inputs,sogs", *sl], [bulk / "sogs" / S / o1 / "sogs"]),
                     (o1a + ["--stage", "states,export", "--tau", TAU, "--results", str(r1)], [r1]),
                     (o1a + ["--stage", "states", "--tau", "0.3"], []),
                     (o1a + ["--stage", "export", "--tau", "0.3", "--no-web-states", "--results",
                             str(r1 / "sensitivity_tau0.3")], []),
                     ([PY, "scripts/sog_slice_report.py", "--sogs", str(bulk / "sogs" / S / o1), "--tau", TAU,
                       "--out", str(r1)], [])]
    st["sogs_O2"] = seed2 + [(o2a + ["--stage", "models,inputs,sogs", *sl], [bulk / "sogs" / S / o2 / "sogs"]),
                     (o2a + ["--stage", "states,export", "--tau", TAU, "--results", str(r2)], [r2]),
                     (o2a + ["--stage", "states", "--tau", "0.3"], []),
                     (o2a + ["--stage", "export", "--tau", "0.3", "--no-web-states", "--results",
                             str(r2 / "sensitivity_tau0.3")], [])]
    # §8.1 benchmark as validation (tau stays frozen at 0.2: --no-choice); O2 pairs feed the D20 split report
    for k_, bdir, rdir in (("bench_O1", bulk / "sogs" / S / o1, r1), ("bench_O2", bulk / "sogs" / S / o2, r2)):
        st[k_] = [([PY, "scripts/sog_tau_benchmark.py", "--sogs", str(bdir), "--annotations", "data/annotations",
                    "--prereg", "config/prereg.yaml", "--results", str(rdir), "--no-choice"],
                   [rdir / "tau_benchmark.csv"])]
    st["split_O2"] = [([PY, "scripts/sog_split_o2.py", "--o1", str(bulk / "sogs" / S / o1), "--o2",
                        str(bulk / "sogs" / S / o2), "--o2-results", str(r2), "--out", str(res / "sogs" / S / f"{o2}_split")],
                       [res / "sogs" / S / f"{o2}_split"])]
    st["rescue_O2"] = [([PY, "scripts/rescue_sogs.py", "--run", str(bulk / "sogs" / S / o2), "--tau", TAU, "--web",
                         str(r2 / "web"), "--out", str(res / "sogs" / S / f"{o2}_rescue"), "--min-p-species", "3",
                         "--workers", str(a.workers), "--work", a.rescue_work or str(bulk / "rescue" / S / o2)],
                        [res / "sogs" / S / f"{o2}_rescue"])]
    pa_work, pa_res = bulk / "partA" / S / R, res / "partA" / S / R
    descent_out = bulk / f"partA_descent_{R}"
    pa = []
    for s in ["01_universe", "02_knn", "03_context", "04_null", "04b_tail_audit", "05_sig", "06_classify",
              "07_estimands", "08_topk_baselines", "09_ligands", "10_figures", "12_functional",
              "13_directional", "14_hidden_homology", "15_v21_tables", "16_descent_2d", "11_numbers"]:
        # 11_numbers last: it collects functional / directional / hidden-homology / descent-2d outputs of 12-16
        extra = []
        if s == "06_classify":
            extra = ["--sogs", str(r1 / "site_sog_all_tau.parquet"), "--tau", TAU]
        if s == "11_numbers":   # never results/numbers.json["partA"] (= the reported snapshot): per-run file only
            extra = ["--numbers-json", "none"]
        if s == "16_descent_2d":
            # descent pair table of THIS snapshot's O1 SOGs (built below) + the frozen species-level TimeTree table
            pa.append(([PY, "scripts/partA_descent/build_descent_pairs.py", "--snapshot", S, "--orthology-version", "O1",
                        "--orthology", O1_ORTH, "--sogs-bulk", str(bulk / "sogs" / S / o1), "--tau", TAU,
                        "--sogs-web", str(r1 / "web" / "sogs.parquet"), "--universe", str(pa_work / "universe.parquet"),
                        "--fingerprints", str(pa_work / "fp.f16.npy"), "--error-rates", str(frozen.path("error_rates")),
                        "--workers", str(min(a.workers, 24)), "--out", str(descent_out)], [descent_out / S / "O1"]))
            extra = ["--datasets", f"O1_tau0.2={descent_out / S / 'O1' / 'pairs.parquet'}",
                     "--timetree", str(frozen.path("timetree_O1")),
                     "--species-tree", str(r1 / "web" / "species_tree.nwk")]
        pa.append(([FAISS_PY, f"scripts/partA/{s}.py", "--snapshot", S, "--sites", str(sites), "--work", str(pa_work),
                    "--results", str(pa_res), *extra], [pa_res] if s == "01_universe" else []))
    st["partA"] = pa
    ev = res / "partB" / "events" / S / o1 / f"tau{TAU}" / "primary"
    th = "results/partB/sim/measured/d11_primary_known_gt10/thresholds.json"
    st["partB_O1_events"] = [([PY, "scripts/partB/run_events.py", "--sogs", str(r1 / "web"), "--orthology", "O1",
                               "--snapshot", S, "--states-glob", str(bulk / "sogs" / S / o1 / f"states_t{TAU}" / "shard=*.parquet"),
                               "--state-col", "state", "--sog-meta", str(r1 / "sogs_all_tau.parquet"), "--tau", TAU,
                               "--thresholds", th, "--prefilter", "--stratum", "primary", "--workers", "16",
                               "--out", str(ev)], [ev])]
    if partb_driver_ok()[0]:
        # O1: pred11/controls are skipped on a new snapshot (they need a pre-D11 SOG release, i.e. an O1_r2
        # equivalent, which the G5 run does not build); the current snapshot keeps the defaults.
        # rbh (D28 orthology_split) must run before artifacts: without it orthology_split = unknown and the cascade
        # passes candidates the reported run rejected (G5: 9 vs 7 gains)
        o1_st = "all" if S == CUR_SNAP else "pre,rescue,post,nulls,candidates,rbh,artifacts"
        pbw = str(min(a.workers, 32))
        st["partB_O1"] = [(["bash", str(PARTB_SH), "--snapshot", S, "--orthology", "O1", "--sog-run", o1,
                            "--run-id", o1, "--stages", o1_st, "--workers", pbw],
                           [res / "partB" / "events" / S / o1, res / "partB" / S / o1 / "candidates.parquet"])]
        # O2 lane: everything that reads only O2 inputs (report is written without the agreement section) ...
        st["partB_O2"] = [(["bash", str(PARTB_SH), "--snapshot", S, "--orthology", "O2", "--sog-run", o2,
                            "--run-id", o2, "--stages", "prep,events,nulls,report,rbh,artifacts", "--o1-run", o1,
                            "--workers", pbw], [res / "partB" / "events" / S / o2])]
        # ... and the O1 -> O2 agreement on the O1 run of the SAME snapshot, after both lanes. report re-writes
        # candidates.parquet, so the artifact cascade (which adds its columns to it) is re-run too (rbh is unchanged)
        st["partB_O1O2"] = [(["bash", str(PARTB_SH), "--snapshot", S, "--orthology", "O2", "--sog-run", o2,
                              "--run-id", o2, "--stages", "agreement,report,artifacts", "--o1-run", o1, "--workers", pbw],
                             [res / "partB" / S / o1 / "o1_o2_gain_agreement.json"])]
    else:   # legacy O1 chain (current snapshot only; see blockers())
        st["partB_O2"] = []
        st["partB_O1O2"] = []
        st["partB_O1"] = [(["bash", "scripts/partB/run_o1_pipeline.sh", o1], [res / "partB" / "events" / S / o1])] + [
        ([PY, f"scripts/partB/{s}.py", "--release", o1], []) for s in
        ["run_artifacts", "false_loss_controls", "run_specificity"]] + [
        ([PY, "scripts/partB/run_m3.py", "--release", o1, "--states-dir", str(bulk / "sogs" / S / o1 / f"states_t{TAU}")], []),
        ([PY, "scripts/partB/write_candidates.py", "--release", o1, "--out", str(res / "partB" / S / o1)],
         [res / "partB" / S / o1])]
    # cases: every input explicit. False-loss controls need a pre-D11 SOG release (not built on a new snapshot):
    # passed as '' (explicitly off) instead of a path that silently does not exist. Case dirs go to the bulk root,
    # never into the repo's cases/ (manual.yaml inputs are only read from there).
    fl = res / "partB" / S / o1 / "false_loss_controls.parquet"
    new = S != CUR_SNAP
    st["cases"] = [([PY, "scripts/cases/build_cases.py", "--orthology", "O1", "--sog-release", o1,
                     "--sog-bulk", str(bulk / "sogs" / S / o1), "--candidates", str(res / "partB" / S / o1 / "candidates.parquet"),
                     "--false-loss", "" if new else str(fl), "--pre-filter-states", "" if new else
                     f"results/sogs/{S}/O1_r2/web/sog_states/bucket=*.parquet",
                     "--redetection", str(frozen.path("redetection_sites")),
                     "--convergence", str(pa_res / "convergence_candidates.parquet"),
                     "--cases-dir", str(bulk / "cases_files" / S / R), "--manual-dir", "cases",
                     "--events-base", str(res / "partB" / "events" / S / o1 / f"tau{TAU}" / "primary"),
                     "--tau-tuning", f"results/sogs/{CUR_SNAP}/O1/tau_tuning_annotations.parquet",
                     "--tau-pairs", f"results/sogs/{CUR_SNAP}/O1/tau_benchmark_pairs.parquet",
                     "--sites", str(sites / "sites"), "--out", str(res / "cases" / S / R), "--version", R,
                     *(["--snapshot", S] if "--snapshot" in help_text([PY, "scripts/cases/build_cases.py"]) else [])],
                    [res / "cases" / S / R])]
    og_extra = (["--input-root", str(res)] if "--input-root" in help_text([PY, "scripts/webapp_build_og_index.py"]) else [])
    st["og_index"] = [([PY, "scripts/webapp_build_og_index.py", "--snapshot", S, "--versions", f"{o1},{o2}", *og_extra],
                       [res / "webapp_index" / S / "og" / o1])]
    return st


def hidden_inputs(a) -> dict:
    """Inputs the entry points resolve INTERNALLY (not visible on their command line). Every one must be either
    snapshot-invariant (pocketevo.frozen, repo data) or produced by an earlier stage of this run; a path keyed on
    the running snapshot that only the first snapshot ever had is exactly the G5 failure mode."""
    th = P / "results/partB/sim/measured/d11_primary_known_gt10/thresholds.json"
    o1x, o2x = P / O1_ORTH / "hogs.orthoxml.gz", P / O2_ORTH / "hogs.orthoxml.gz"
    return {
        "partA": [Path(FAISS_PY), Path("/tmp/pe_big/partA/foldseek/afdb_v6.fasta"),
                  Path("/tmp/pe_big/work/annotations/uniprot_sequences.parquet"),
                  P / "data/annotations/uniprot_proteins.parquet", frozen.path("timetree_O1"),
                  frozen.path("error_rates")],
        "partB_O1_events": [th, o1x],
        "partB_O1": [frozen.path("rescue_k"), th, o1x, P / "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk",
                     P / "data/annotations/uniprot_proteins.parquet", P / "data/annotations/ted_domains.parquet",
                     P / "data/orthology/ancestra_O2/tables/protein_hogs.parquet"],
        "partB_O2": [th, o2x, P / "data/orthology/ancestra_O2/o2_guide_tree.nwk"],
        "partB_O1O2": [o1x, o2x, th],
        "cases": [P / "results/partB/hog_trees/O1/hog_trees.parquet", frozen.path("redetection_sites")],
    }


def _is_produced(tok: str, a) -> bool:
    """A path token is an upstream output of this run if it names the run's snapshot or lives under its run dirs."""
    return a.snapshot in tok or f"_{a.run_id}" in tok or f"/{a.run_id}" in tok


def check_inputs(a, st, want) -> list:
    """Dry-run input check: every path argument that is not an output and not produced upstream, plus every hidden
    input, must exist now. Returns [(stage, path)] of missing inputs."""
    import glob as _glob
    miss = []
    hid = hidden_inputs(a)
    for name in [x for x in STAGE_ORDER if x in want and x in st]:
        paths = list(hid.get(name, []))
        for cmd, _ in st[name]:
            if cmd[:2] == ["bash", "-c"]:      # MSA seed copy: check its source dir
                paths += [Path(t) for t in cmd[2].split() if t.startswith("/") and t.endswith("/msa") and a.snapshot not in t]
                continue
            prev = ""
            for t in cmd[1:]:
                vals = [t.split("=", 1)[1]] if ("=" in t and t.split("=", 1)[1].startswith("/")) else [t]
                for v in vals:
                    if prev not in OUTPUT_FLAGS and ("/" in v) and not v.startswith("-") and not _is_produced(v, a) \
                            and not v.startswith("http"):
                        paths.append(Path(v) if Path(v).is_absolute() else P / v)
                prev = t
        for p_ in paths:
            if not (_glob.glob(str(p_)) if any(c in str(p_) for c in "*?[") else p_.exists()):
                miss.append((name, str(p_)))
    for n, p_ in frozen.all_paths().items():
        if not p_.exists():
            miss.append(("frozen:" + n, str(p_)))
    return miss


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True, help="snapshot id, or 'new' to take one first")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--stages", default=",".join(STAGE_ORDER))
    ap.add_argument("--families", default=None, help="accessions; run only their OGs (SOG stages: --slice)")
    ap.add_argument("--bulk-root", default="/tmp/pe_big")
    ap.add_argument("--results-root", default=str(P / "results"))
    ap.add_argument("--snapshot-root", default=str(P / "data/snapshot"))
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--msa-seed", default=None, help="bulk sogs root holding O1_r2/O2_r1 MSA caches "
                    "(e.g. /tmp/pe_big/sogs); saves the FAMSA time (~19 + 48 min)")
    ap.add_argument("--resume", action="store_true", help="allow existing (uncommitted) outputs")
    ap.add_argument("--rescue-work", default=None, help="BURROW work dir to reuse (outputs are reused per model stem)")
    ap.add_argument("--log", default=None)
    ap.add_argument("--substeps", default=None, help="comma list of substrings: within the chosen stages run only "
                    "the commands containing one of them (e.g. --stages partA --substeps 04b_,06_,07_)")
    a = ap.parse_args()
    want = a.stages.split(",")
    if a.snapshot == "new":
        if a.dry_run:
            a.snapshot = "<NEW_SNAPSHOT_ID>"
        else:
            r = subprocess.run(stages(argparse.Namespace(**{**vars(a), "snapshot": "x"}))["snapshot"][0][0], cwd=P,
                               capture_output=True, text=True)
            if r.returncode:
                raise SystemExit(r.stderr[-2000:])
            a.snapshot = sorted(os.listdir(a.snapshot_root))[-1]
            want = [w for w in want if w != "snapshot"]
    st, bl = stages(a), blockers(a)
    if a.substeps:
        subs = a.substeps.split(",")
        st = {k: [c for c in v if any(x in " ".join(c[0]) for x in subs)] for k, v in st.items()}
    log = Path(a.log or Path(a.bulk_root) / "pipeline_logs" / f"{a.snapshot}_{a.run_id}.jsonl")
    if not a.dry_run:
        log.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, UV_CACHE_DIR="/tmp/pe_big/uvcache")
    for name in [s for s in STAGE_ORDER if s in want]:
        if name in bl:
            print(f"# [{name}] BLOCKED: {bl[name]}")
            if not a.dry_run:
                with open(log, "a") as fh:
                    fh.write(json.dumps({"stage": name, "rc": None, "blocked": bl[name], "wall_s": 0}) + "\n")
            continue
        for cmd, outs in st[name]:
            for o in outs:
                if git_tracked(Path(o)):
                    raise SystemExit(f"[{name}] refuses: {o} is a committed run (choose another --run-id)")
                if Path(o).exists() and not (a.resume or a.dry_run):
                    raise SystemExit(f"[{name}] refuses: {o} exists (use --resume or another --run-id)")
            e = dict(env)
            if name == "partA":
                e["CUDA_VISIBLE_DEVICES"] = "3"
            elif name == "rescue_O2":
                e["CUDA_VISIBLE_DEVICES"] = ""
            line = (f"CUDA_VISIBLE_DEVICES={e['CUDA_VISIBLE_DEVICES']} " if "CUDA_VISIBLE_DEVICES" in e and name in
                    ("partA", "rescue_O2") else "") + shlex.join(cmd)
            print(f"[{name}] {line}", flush=True)
            if a.dry_run:
                continue
            t0 = time.time()
            r = subprocess.run(cmd, cwd=P, env=e, capture_output=True, text=True)
            rec = {"stage": name, "cmd": cmd, "rc": r.returncode, "wall_s": round(time.time() - t0, 1),
                   "stderr_tail": r.stderr[-1500:] if r.returncode else "", "stdout_tail": r.stdout[-600:]}
            with open(log, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            print(f"    rc={r.returncode} wall={rec['wall_s']}s", flush=True)
            if r.returncode:
                print(r.stderr[-1500:])
                raise SystemExit(f"[{name}] failed")
    if a.dry_run:
        miss = check_inputs(a, st, [w for w in want if w not in bl])
        for name, p_ in miss:
            print(f"# MISSING INPUT [{name}] {p_}")
        print(f"# input check: {'all inputs found' if not miss else f'{len(miss)} missing'} "
              f"(snapshot {a.snapshot}, run {a.run_id})")
        if miss:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
