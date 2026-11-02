#!/usr/bin/env python
"""Rescue of SOG A calls (prereg states.rescue; D9: k = 1; hit = rescue pocket with residue Jaccard >= 0.3 against
the member's footprint residues, the rule k was calibrated with). Uses src/pocketevo/rescue.py (6 BURROW settings).

Selection: A rows of a SOG run with `rescue:pending` (model has atlas-filtered pockets, n_pockets_all > kept), or all
A rows with --all-A. Optional pool restriction --min-p-species N (SOGs with P in >= N species = Part B candidate pool).
SOGs flagged small_footprint / empty_footprint are skipped (excluded from candidates, D8 / verifier).

k = 1 makes the call "A -> W iff any setting hits", so settings run sequentially (prereg order) and a protein leaves
the queue once all its selected rows are hit (early stop; identical W/A result for k = 1; `n_settings_hit` is then a
lower bound and `n_settings_run` is the number actually run for that row).
Models: tier 1/2 AFDB v6 mmCIF from the tars, tier 3 / L Foldcomp v4 decompressed to PDB (B-factor = pLDDT).
Multi-fragment models are not rescued: A -> U (`rescue:not_run_multifragment`), as in scripts/partB/run_rescue.py.

    .venv/bin/python scripts/rescue_sogs.py --run /tmp/pe_big/sogs/20261030T085136Z/O2_r1 --tau 0.2 \\
        --web results/sogs/20261030T085136Z/O2_r1/web --out results/sogs/20261030T085136Z/O2_r1_rescue \\
        [--sample 1000 --seed 0]   # budgeting
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")   # CPU only (GPUs reserved for co-folding)

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "src"))
sys.path.insert(0, str(P / "scripts"))
from pocketevo import rescue as R  # noqa: E402
from build_sogs import git_commit, model_index  # noqa: E402

AFDB_V6 = Path("/cluster/work/beltrao/shared/alphafold/v6")
FOLDCOMP_DIR = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb/foldcomp")


def member_footprints(rows: pd.DataFrame, run: Path, tau: float) -> dict:
    """{(sog_id, accession): set of member (= model) residue positions at the footprint columns}."""
    out = {}
    for shard, g in rows.groupby("shard"):
        sg = pd.read_parquet(run / "sogs" / f"shard={shard:04d}.parquet", columns=["og_id", "sog_id", "tau", "footprint_columns"])
        sg = sg[np.isclose(sg.tau, tau) & sg.sog_id.isin(set(g.sog_id))].set_index("sog_id")
        msa = pd.read_parquet(run / "msa" / f"shard={shard:04d}.parquet")
        msa = msa[msa.og_id.isin(set(g.og_id))].set_index(["og_id", "accession"]).aligned
        for og, sog, acc in zip(g.og_id, g.sog_id, g.accession):
            a = np.frombuffer(msa[(og, acc)].encode(), np.uint8) != ord("-")
            pos = np.cumsum(a)
            c = np.asarray(sg.at[sog, "footprint_columns"], dtype=int)
            out[(sog, acc)] = set(pos[c][a[c]].tolist())
    return out


def _extract_chunk(args):
    import foldcomp
    kind, src, rows, outdir = args
    path = AFDB_V6 / src if kind == "tar_cif_gz" else FOLDCOMP_DIR / src
    out = {}
    with open(path, "rb") as fh:
        for acc, off, size in rows:
            ext, ver = (".cif.gz", "v6") if kind == "tar_cif_gz" else (".pdb", "v4")
            p = Path(outdir) / f"AF-{acc}-F1-model_{ver}{ext}"
            if not p.exists():
                fh.seek(int(off))
                blob = fh.read(int(size))
                if kind == "tar_cif_gz":
                    p.write_bytes(blob)
                else:
                    _, pdb = foldcomp.decompress(blob)
                    p.write_text(pdb)
            out[acc] = str(p)
    return out


def extract(models: pd.DataFrame, outdir: Path, workers: int = 16) -> dict:
    """models: model_accession, blob_kind, blob_file, blob_offset, blob_size (fragment 1 only) -> {model_acc: path}.
    Parallel over offset-sorted chunks of each source file."""
    import concurrent.futures as cf
    outdir.mkdir(parents=True, exist_ok=True)
    tasks = []
    for (kind, src), g in models.groupby(["blob_kind", "blob_file"]):
        g = g.sort_values("blob_offset")
        rows = list(zip(g.model_accession, g.blob_offset, g.blob_size))
        for k in range(0, len(rows), 500):
            tasks.append((kind, src, rows[k:k + 500], str(outdir)))
    paths = {}
    with cf.ProcessPoolExecutor(workers) as ex:
        for r in ex.map(_extract_chunk, tasks):
            paths.update({k: Path(v) for k, v in r.items()})
    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="bulk SOG run dir (/tmp/pe_big/sogs/<snap>/<run_id>)")
    ap.add_argument("--web", required=True, help="results web dir of the run (sogs.parquet flags)")
    ap.add_argument("--tau", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=1, help="D9: k = 1")
    ap.add_argument("--hit-tau", type=float, default=0.3)
    ap.add_argument("--all-A", action="store_true")
    ap.add_argument("--min-p-species", type=int, default=0)
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=40)
    ap.add_argument("--work", default=None)
    ap.add_argument("--sogs-file", default=None, help="restrict to SOG ids listed in this file (one per line)")
    ap.add_argument("--skip-overlay", default=None, help="skip (sog_id, accession) rows already in this rescue_A.parquet")
    a = ap.parse_args()
    if a.k != 1:
        raise SystemExit("early stopping is only valid for k = 1")
    t0 = time.time()
    run = Path(a.run)
    snap = run.parent.name
    st = ds.dataset(sorted(glob.glob(str(run / f"states_t{a.tau:g}" / "shard=*.parquet"))), format="parquet").to_table(
        columns=["og_id", "sog_id", "accession", "taxid", "state", "reason", "n_filtered_pockets", "model_accession",
                 "model_tier", "footprint_plddt"], filter=ds.field("state") == "A").to_pandas()
    sg = pd.read_parquet(Path(a.web) / "sogs.parquet", columns=["sog_id", "n_P", "small_footprint", "empty_footprint"])
    nP_sp = None
    n_all = len(st)
    sel = st if a.all_A else st[st.n_filtered_pockets > 0]
    n_pending = len(sel)
    bad = set(sg.sog_id[sg.small_footprint | sg.empty_footprint])
    sel = sel[~sel.sog_id.isin(bad)]
    n_flag_skipped = n_pending - len(sel)
    if a.min_p_species:
        stP = ds.dataset(sorted(glob.glob(str(run / f"states_t{a.tau:g}" / "shard=*.parquet"))), format="parquet").to_table(
            columns=["sog_id", "taxid"], filter=ds.field("state") == "P").to_pandas()
        nP_sp = stP.groupby("sog_id").taxid.nunique()
        sel = sel[sel.sog_id.map(nP_sp).fillna(0) >= a.min_p_species]
    if a.sogs_file:
        keep = set(Path(a.sogs_file).read_text().split())
        sel = sel[sel.sog_id.isin(keep)]
    if a.skip_overlay:
        prev = pd.read_parquet(a.skip_overlay, columns=["sog_id", "accession"])
        prev = set(zip(prev.sog_id, prev.accession))
        sel = sel[[k not in prev for k in zip(sel.sog_id, sel.accession)]]
    if a.sample:
        sel = sel.sample(n=min(a.sample, len(sel)), random_state=a.seed)
    sel = sel.reset_index(drop=True)
    ogs = pd.read_parquet(run / "inputs/ogs.parquet", columns=["og_id", "shard"])
    sel = sel.merge(ogs, on="og_id")
    print(f"[{time.time() - t0:.0f}s] A rows {n_all}, pending {n_pending}, flag-skipped {n_flag_skipped}, "
          f"selected {len(sel)} rows / {sel.model_accession.nunique()} models", flush=True)
    fps = member_footprints(sel, run, a.tau)
    # models
    class _C:  # minimal ctx for model_index
        pass
    ctx = _C()
    ctx.snapshot = snap
    ctx.a = argparse.Namespace(tierl_entries="/tmp/pe_big/tierL/fc_entries_tierL.parquet")
    ix = model_index(ctx, [1, 2, 3, 4])
    ix = ix[ix.accession.isin(set(sel.model_accession))]
    nfrag = ix.groupby("accession").fragment.max()
    multi = set(nfrag.index[nfrag > 1])
    one = ix[(ix.fragment == 1) & ~ix.accession.isin(multi)].rename(columns={"accession": "model_accession"})
    work = Path(a.work or f"/tmp/pe_big/rescue/{snap}/{run.name}")
    t1 = time.time()
    paths = extract(one, work / "models")
    t_extract = time.time() - t1
    print(f"[{time.time() - t0:.0f}s] extracted {len(paths)} models ({t_extract:.0f}s), multi-fragment {len(multi)}",
          flush=True)
    # sequential settings with early stop (k = 1)
    sel["key"] = list(zip(sel.sog_id, sel.accession))
    todo = sel[sel.model_accession.isin(paths.keys())].copy()
    hits = {k: 0 for k in todo.key}
    nrun = {k: 0 for k in todo.key}
    hit_by = {k: None for k in todo.key}
    per_setting = []
    open_keys = set(todo.key)
    for name in R.SETTING_NAMES:
        g = todo[todo.key.isin(open_keys)]
        models = sorted(set(g.model_accession))
        if not models:
            break
        t2 = time.time()
        outs = R.run_settings([str(paths[m]) for m in models], work / "burrow", workers=a.workers,
                              settings={name: R.SETTINGS[name]})
        dt = time.time() - t2
        n_hit = 0
        for key, macc in zip(g.key, g.model_accession):
            h = R.setting_hit(R.read_burrow_csv(Path(outs[name]) / f"{R.model_stem(paths[macc])}.csv"),
                              fps[key], a.hit_tau, "jaccard")
            nrun[key] += 1
            if h["hit"]:
                hits[key] += 1
                hit_by[key] = hit_by[key] or name
                open_keys.discard(key)
                n_hit += 1
        per_setting.append({"setting": name, "n_models": len(models), "n_rows": len(g), "n_hit": n_hit, "wall_s": dt,
                            "s_per_model": dt / len(models)})
        print(f"[{time.time() - t0:.0f}s] {name}: {len(models)} models, {n_hit}/{len(g)} rows hit, {dt:.0f}s", flush=True)
    sel["n_settings_hit"] = sel.key.map(hits)
    sel["n_settings_run"] = sel.key.map(nrun)
    sel["hit_setting"] = sel.key.map(hit_by)
    sel["state_rescued"] = np.where(sel.n_settings_hit.fillna(0) >= a.k, "W", "A")
    sel["reason_rescued"] = np.where(sel.state_rescued == "W", "rescue:" + sel.hit_setting.fillna(""), "rescue:negative")
    nr = sel.n_settings_run.isna()
    sel.loc[nr, "state_rescued"] = "U"
    sel.loc[nr, "reason_rescued"] = np.where(sel.loc[nr, "model_accession"].isin(multi), "rescue:not_run_multifragment",
                                             "rescue:not_run")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"_sample{a.sample}_seed{a.seed}" if a.sample else ""
    res = sel.drop(columns=["key"])
    res.to_parquet(out / f"rescue_A{tag}.parquet", index=False)
    wall = time.time() - t0
    burrow_s = float(sum(r["wall_s"] for r in per_setting))
    summ = {"snapshot_id": snap, "run": str(run), "tau": a.tau, "k": a.k, "hit_tau": a.hit_tau, "all_A": a.all_A,
            "min_p_species": a.min_p_species, "sample": a.sample, "seed": a.seed, "workers": a.workers,
            "n_A_rows": n_all, "n_pending": n_pending, "n_flag_skipped": n_flag_skipped, "n_selected": len(sel),
            "n_models": int(sel.model_accession.nunique()), "n_multifragment": len(multi),
            "per_setting": per_setting, "extract_s": t_extract, "burrow_s": burrow_s, "wall_s": wall,
            "result": res.state_rescued.value_counts().to_dict(),
            "result_by_tier": {f"{k[0]}|{k[1]}": int(v) for k, v in res.groupby(["model_tier", "state_rescued"]).size().items()},
            "git_commit": git_commit(), "prereg_sha256": json.load(open(run / "inputs/members.parquet.json"))["prereg_sha256"]}
    json.dump(summ, open(out / f"rescue_summary{tag}.json", "w"), indent=1, default=str)
    json.dump({k: summ[k] for k in ("snapshot_id", "prereg_sha256", "git_commit")} | {"seed": a.seed, "code_version": summ["git_commit"]},
              open(out / f"rescue_A{tag}.parquet.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != "per_setting"}, indent=1, default=str))


if __name__ == "__main__":
    main()
