#!/usr/bin/env python3
"""Drive the frozen case pool through co-folding (Tier 1) and write results into case.yaml.

Tier 1 cases have a cognate ligand in BioLiP2/UniProt and get the full ligand co-folding
treatment: 3 present-side + 3 absent-side members x (cognate + 1 in-class property-matched
decoy with a Tanimoto floor) x 3 seeds, no templates, blinded, plus the family's holo PDB
positive-control gate. Scored with the FROZEN threshold in
results/cofold/calibration/threshold.json -- never recomputed here.

Writes cases/<id>/cofold/cofold_spec.yaml (the resolved spec, auditable) and updates the
`cofold:` section of cases/<id>/case.yaml with status, readouts and the interpretation
rules that follow from the calibration.

Usage: run_pool.py --tier1 --gpus 0,1,2 [--only CASE_ID]
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
from pathlib import Path
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_cofold as CC

P = Path(__file__).resolve().parents[2]  # repo root
AF3_PY = os.environ.get("AF3_ENV", "/tmp/af3_env") + "/bin/python"
THRESH = P / "results/cofold/calibration/threshold.json"

# case_id -> cognate ligand CCD codes, decoy CCD codes, family holo positive control.
# Ligands come from BioLiP2 rows of the case's own members (see the commit message);
# positive controls are themselves case members, so the gate is family-matched.
TIER1 = {
    "gain_HOG0035403_g5_n14_O1": {
        "ligands": ["CCD:OAA"], "decoys": [["CCD:OKG"]],
        "family": "citrate synthase (EC 2.3.3.16)",
        "positive_control": {"acc": "P0ABH7", "pdb": "4jag", "chain": "A",
                             "ligand": "OAA", "seqnum": None},
    },
    "gain_HOG0035156_g10_n21_O1": {
        "ligands": ["CCD:MET"], "decoys": [["CCD:9X6"]],
        "family": "MetE methyltransferase (EC 2.1.1.14)",
        "positive_control": {"acc": "O50008", "pdb": "1u1j", "chain": "A",
                             "ligand": "MET", "seqnum": None},
    },
    "dupgain_HOG0054071_g9_n185_O2": {
        "ligands": ["CCD:ADP", "CCD:MG"], "decoys": [["CCD:TTP", "CCD:MG"]],
        "family": "choline/ethanolamine kinase (EC 2.7.1.32)",
        "positive_control": {"acc": "P35790", "pdb": "3g15", "chain": "A",
                             "ligand": "ADP", "seqnum": None},
    },
    "loss_HOG0027153_g11_n17_O1": {
        "ligands": ["CCD:APC", "CCD:MG"], "decoys": [["CCD:TTP", "CCD:MG"]],
        "family": "DNA-directed RNA polymerase (EC 2.7.7.6)",
        "positive_control": {"acc": "O00411", "pdb": "8u8v", "chain": "E",
                             "ligand": "APC", "seqnum": None},
    },
}

BLINDING_TIER1 = {
    "status": "UNBLINDED",
    "prereg_requirement": "cofold.blinding: derived/ancestral labels hashed in file names; "
                          "the agent scoring models sees no state label",
    "what_actually_happened":
        "The folds themselves were run blind: inputs carried HMAC-hashed names and the "
        "state labels existed only in blind_key.json. But build() minted a fresh salt on "
        "every call, which destroyed the original key, and recovering it required matching "
        "each AF3 output back to its member by sequence identity. That reconstruction, and "
        "the debugging around it, happened BEFORE the final scoring pass, so the state "
        "labels were visible to me while the readouts were being produced.",
    "honest_claim": "These four cases are reported as UNBLINDED. The readout code is "
                    "deterministic and never receives a state label, so the measurements "
                    "themselves are label-independent, but the preregistered blinding "
                    "claim cannot be made for them.",
    "remedy": "The salt is now persisted, so cases folded from here on can be scored "
              "without ever reconstructing the mapping, and will be labelled blinded.",
}

INTERPRETATION = [
    "A composite 'binds' call is informative: specificity 0.96 on the frozen panel.",
    "A composite 'does not bind' call is WEAK: sensitivity 0.42, so co-folding misses "
    "most true complexes. A gain or loss is 'co-fold-consistent' only when the present "
    "side binds AND the absent side does not, and the absent-side half stays weak.",
    "Metal-ion ligands are uninterpretable by co-folding (panel AUC 0.58 for ligand "
    "pLDDT, 0.56 for min PAE, i.e. at chance).",
]


def log(m): print(f"[run_pool] {m}", file=sys.stderr, flush=True)


def spec_path(cid: str) -> Path:
    d = P / "cases" / cid / "cofold"
    d.mkdir(parents=True, exist_ok=True)
    return d / "cofold_spec.yaml"


def write_spec(cid: str, cfg: dict) -> Path:
    doc = yaml.safe_load((P / "cases" / cid / "case.yaml").read_text())
    doc["cofold_ligands"] = cfg["ligands"]
    doc["cofold_decoys"] = cfg["decoys"]
    doc["positive_control"] = cfg["positive_control"]
    doc["cofold_family"] = cfg["family"]
    sp = spec_path(cid)
    sp.write_text(yaml.safe_dump(doc, sort_keys=False, width=100))
    return sp


def update_case_yaml(cid: str, dec: dict, table: Path, cfg: dict) -> None:
    f = P / "cases" / cid / "case.yaml"
    doc = yaml.safe_load(f.read_text())
    gate = dec.get("positive_control_gate", {})
    # MERGE, never replace: other producers add to this section (model entries carrying a
    # `file`, the ligand-vs-pocket analysis on the summary entry). A rescore that replaced
    # the whole list silently deleted them once already.
    prev = [c for c in (doc.get("cofold") or []) if isinstance(c, dict)]
    prev_summary = next((c for c in prev if not c.get("file")), {})
    prev_models = [c for c in prev if c.get("file")]
    doc["cofold"] = [{
        "model_id": cid,
        "status": dec.get("verdict"),
        "role": "case",
        "family": cfg["family"],
        "ligands": cfg["ligands"], "decoys": cfg["decoys"],
        "msa": "ColabFold MMseqs2 for all members (prereg cofold.msa)",
        "templates": "none (prereg cofold.templates.primary)",
        "seeds": 3, "samples_per_seed": 5,
        "threshold_file": "results/cofold/calibration/threshold.json",
        "threshold_sha256": (P / "results/cofold/calibration/threshold.json.sha256"
                             ).read_text().split()[0]
        if (P / "results/cofold/calibration/threshold.json.sha256").exists() else None,
        "positive_control_gate": gate,
        "verdict_basis": dec.get("verdict_basis"),
        "legacy_prereg_rule": dec.get("legacy_prereg_rule"),
        "readouts_table": str(table.relative_to(P)) if table.exists() else None,
        "interpretation": INTERPRETATION,
        "blinding": BLINDING_TIER1,
        "d23_composite": dec.get("d23_composite"),
        "note": "scored against the frozen threshold, which was fixed before any case was "
                "folded; see `blinding` for why these four cases are labelled UNBLINDED",
    }]
    keep = ("cognate_site_vs_atlas_pockets", "footprint_vs_cognate_site")
    for k in keep:
        if k in prev_summary and k not in doc["cofold"][0]:
            doc["cofold"][0][k] = prev_summary[k]
    doc["cofold"] += prev_models
    tmp = f.with_suffix(".yaml.tmp")
    tmp.write_text(yaml.safe_dump(doc, sort_keys=False, width=100))
    tmp.replace(f)          # atomic: never leaves a half-written case.yaml
    log(f"{cid}: case.yaml cofold section updated -> {dec.get('verdict')}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier1", action="store_true")
    ap.add_argument("--gpus", default="0,1,2")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--skip-run", action="store_true")
    ap.add_argument("--score-only", action="store_true",
                    help="score existing outputs; do not rebuild inputs or keys")
    a = ap.parse_args()
    if not THRESH.exists():
        raise SystemExit("frozen threshold missing; refusing to score cases")
    cases = {k: v for k, v in TIER1.items() if not a.only or k in a.only}
    prereg = CC.load_prereg()
    for cid, cfg in cases.items():
        log(f"=== {cid} ({cfg['family']})")
        sp = write_spec(cid, cfg)
        spec = CC.load_spec(sp)
        wd = Path("/tmp/pe_big/cofold/cases") / cid
        wd.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        try:
            if a.score_only:
                # Never call build() when scoring existing outputs: it rewrites jobs.json
                # and blind_key.json, which orphans the already-folded AF3 directories.
                if not (wd / "jobs.json").exists():
                    raise SystemExit(f"{wd}: no jobs.json, nothing to score")
            else:
                CC.build(spec, "colabfold", wd)
                if not a.skip_run:
                    CC.run(wd, [int(g) for g in a.gpus.split(",")])
            CC.collect(wd, spec)
            res = P / "results/cofold/cases" / cid
            dec = CC.unblind(wd, spec, prereg, res)
            dec["blinding"] = BLINDING_TIER1          # all four Tier 1 cases: see BLINDING_TIER1
            tmpd = res / "decision.json.tmp"
            tmpd.write_text(json.dumps(dec, indent=1)); tmpd.replace(res / "decision.json")
            update_case_yaml(cid, dec, res / "case_table.tsv", cfg)
            log(f"{cid}: {dec['verdict']} in {time.time()-t0:.0f}s")
        except SystemExit as e:
            log(f"{cid}: ABORTED {e}")
        except Exception as e:  # noqa: BLE001
            import traceback; traceback.print_exc()
            log(f"{cid}: FAILED {e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
