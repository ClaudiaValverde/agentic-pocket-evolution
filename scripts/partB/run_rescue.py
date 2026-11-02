#!/usr/bin/env python
"""Rescue (prereg states.rescue; D9 k = 1) for every member of SOGs that carry a surviving gain/loss candidate,
then re-run the event calling on the rescued states (DESIGN §3.3, §5.2, §5.4).

Symmetric rule (states.apply_rescue): A -> W if >= k of the 6 settings hit the footprint; P -> W if < k settings
hit. Hit = a rescue pocket with residue Jaccard >= --hit-tau (0.3 = the rule k was calibrated with,
results/observation/<snap>/rescue_k.json) against the member's footprint residues (MSA footprint columns mapped to
UniProt positions of that member). Models: AFDB v6 tar members via results/webapp_index (blob offsets).
Multi-fragment models are not rescued: their A calls become U (`rescue:not_run_multifragment`; conservative), P
unchanged (counted).

    python scripts/partB/run_rescue.py --tau 0.2 --events results/partB/events/<snap>/O1/tau0.2/primary
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import frozen, recon, rescue as R, states as ST  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg
AFDB_V6 = Path("/cluster/work/beltrao/shared/alphafold/v6")


def member_footprints(sogs: set, tau: float, bulk: Path, accs_by_sog: dict) -> dict:
    """{(sog_id, accession): set of UniProt positions at the footprint columns (non-gap)}."""
    fp = ds.dataset(sorted(glob.glob(str(bulk / "sogs/shard=*.parquet"))), format="parquet").to_table(
        columns=["og_id", "sog_id", "tau", "footprint_columns"]).to_pandas()
    fp = fp[np.isclose(fp["tau"], tau) & fp["sog_id"].isin(sogs)]
    ogs = set(fp["og_id"])
    msa = ds.dataset(sorted(glob.glob(str(bulk / "msa/shard=*.parquet"))), format="parquet").to_table().to_pandas()
    msa = msa[msa["og_id"].isin(ogs)]
    pos = {}
    for og, acc, aln in zip(msa["og_id"], msa["accession"], msa["aligned"]):
        a = np.frombuffer(aln.encode(), dtype=np.uint8) != ord("-")
        p = np.cumsum(a)               # UniProt position at each column (if non-gap)
        pos[(og, acc)] = (a, p)
    out = {}
    for og, sog, cols in zip(fp["og_id"], fp["sog_id"], fp["footprint_columns"]):
        cols = np.asarray(cols, dtype=int)
        for acc in accs_by_sog.get(sog, ()):
            if (og, acc) not in pos:
                continue
            a, p = pos[(og, acc)]
            c = cols[cols < len(a)]
            out[(sog, acc)] = set(p[c][a[c]].tolist())
    return out


def extract_models(accs: set, outdir: Path) -> tuple[dict, dict]:
    """{accession: model path} for single-fragment tier-1/2 models; counters."""
    outdir.mkdir(parents=True, exist_ok=True)
    idx = pd.concat([pd.read_parquet(f, columns=["accession", "fragment", "blob_file", "blob_offset", "blob_size"])
                     for t in ("tier1", "tier2") for f in sorted(glob.glob(
                         str(ROOT / f"results/webapp_index/{SNAP}/{t}/bucket=*.parquet")))], ignore_index=True)
    idx = idx[idx["accession"].isin(accs)]
    nfrag = idx.groupby("accession")["fragment"].max()
    multi = set(nfrag.index[nfrag > 1])
    one = idx[~idx["accession"].isin(multi) & (idx["fragment"] == 1)].drop_duplicates("accession")
    paths, cnt = {}, {"multi_fragment": len(multi), "no_model": len(accs - set(idx["accession"]))}
    for tar, g in one.groupby("blob_file"):
        fh = None
        for acc, off, size in zip(g["accession"], g["blob_offset"], g["blob_size"]):
            p = outdir / f"AF-{acc}-F1-model_v6.cif.gz"
            if not p.exists():
                if fh is None:
                    fh = open(AFDB_V6 / tar, "rb")
                fh.seek(int(off))
                p.write_bytes(fh.read(int(size)))
            paths[acc] = p
        if fh is not None:
            fh.close()
    cnt["extracted"] = len(paths)
    return paths, cnt, multi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tau", type=float, required=True)
    ap.add_argument("--events", required=True, help="run_events.py output dir (events.parquet)")
    ap.add_argument("--state-col", default="state")
    ap.add_argument("--include-frozen-loss", action="store_true",
                    help="also rescue SOGs whose only candidates are frozen-rule (uncalibrated) losses")
    ap.add_argument("--hit-tau", type=float, default=0.3)
    ap.add_argument("--k", type=int, default=None, help="default: frozen D9 value, results/observation/"
                    "<calibration snapshot>/rescue_k.json k_min_false_A (= 1; snapshot-invariant, pocketevo.frozen)")
    ap.add_argument("--work", default=f"/tmp/pe_big/partB/{SNAP}/rescue")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--out", default=None)
    ap.add_argument("--bulk", default=f"/tmp/pe_big/sogs/{SNAP}/O1", help="SOG bulk dir (msa/, sogs/, states_t*)")
    ap.add_argument("--states-dir", default=None, help="state shards dir (default bulk states_t<tau>; D11: "
                                                       "/tmp/pe_big/partB/<snap>/O1/states_d11_t<tau>)")
    recon.add_partb_args(ap)
    a = ap.parse_args()
    t0 = time.time()
    rk = json.loads(frozen.path("rescue_k").read_text())   # D9 frozen decision: not keyed on the running snapshot
    k = a.k or int(rk["k"] if rk.get("k") else rk["k_min_false_A"])
    k_prov = {"source": "--k" if a.k else "frozen", **({} if a.k else frozen.provenance("rescue_k"))}
    bulk = Path(a.bulk)
    ev = pd.read_parquet(Path(a.events) / "events.parquet")
    sel = ev["called"] | (ev["candidate"] & (ev["type"] == "loss") & a.include_frozen_loss)
    sogs = set(ev.loc[sel, "char"])
    sdir = Path(a.states_dir) if a.states_dir else bulk / f"states_t{a.tau}"
    st = ds.dataset(sorted(glob.glob(str(sdir / "shard=*.parquet"))), format="parquet").to_table(
        columns=["og_id", "sog_id", "taxid", "accession", a.state_col, "footprint_plddt", "reason"]).to_pandas()
    st = st[st["sog_id"].isin(sogs)].rename(columns={a.state_col: "state"})
    todo = st[st["state"].isin(["P", "A"])]
    accs_by_sog = todo.groupby("sog_id")["accession"].apply(set).to_dict()
    print(f"[{time.time() - t0:.0f}s] {len(sogs)} SOGs, {len(todo)} P/A rows, {todo['accession'].nunique()} proteins",
          flush=True)
    fps = member_footprints(sogs, a.tau, bulk, accs_by_sog)
    paths, mcnt, multi = extract_models(set(todo["accession"]), Path(a.work) / "models")
    print(f"[{time.time() - t0:.0f}s] models {mcnt}", flush=True)
    t1 = time.time()
    outs = R.run_settings([str(p) for p in sorted(paths.values())], Path(a.work) / "burrow", workers=a.workers)
    t_burrow = time.time() - t1
    print(f"[{time.time() - t0:.0f}s] burrow done ({t_burrow:.0f}s)", flush=True)
    rows, nfail = [], {"no_footprint": 0, "burrow_output_missing": 0}
    for sog, acc in zip(todo["sog_id"], todo["accession"]):
        if acc not in paths:
            continue
        f = fps.get((sog, acc))
        if not f:
            nfail["no_footprint"] += 1
            continue
        try:
            h = R.hits_from_outputs(outs, R.model_stem(paths[acc]), f, a.hit_tau, "jaccard")
        except FileNotFoundError:
            nfail["burrow_output_missing"] += 1
            continue
        rows.append(dict(sog_id=sog, accession=acc, n_settings_hit=R.n_hits(h), n_settings_run=len(h),
                         **{f"hit_{s}": v["hit"] for s, v in h.items()}))
    rs = pd.DataFrame(rows)
    out = Path(a.out) if a.out else Path(a.events) / "rescue"
    out.mkdir(parents=True, exist_ok=True)
    resc = ST.apply_rescue(st, rs, k)
    # members that could not be rescued: A -> U (conservative), P unchanged
    nr = (resc["state"] == "A") & ~resc["rescue_applied"]
    resc.loc[nr, "state"] = "U"
    resc.loc[nr, "reason"] = np.where(resc.loc[nr, "accession"].isin(multi), "rescue:not_run_multifragment",
                                      "rescue:not_run")
    summ = dict(k=k, k_provenance=k_prov, hit_tau=a.hit_tau, tau=a.tau, n_sogs=len(sogs), n_PA_rows=len(todo),
                n_proteins=int(todo["accession"].nunique()), models=mcnt, failures=nfail, burrow_s=t_burrow,
                n_rescue_rows=len(rs),
                A_to_W=int(((st["state"].values == "A") & (resc["state"].values == "W")).sum()),
                P_to_W=int(((st["state"].values == "P") & (resc["state"].values == "W")).sum()),
                A_to_U_not_rescued=int(nr.sum()),
                A_before=int((st["state"] == "A").sum()), A_after=int((resc["state"] == "A").sum()),
                P_before=int((st["state"] == "P").sum()), P_after=int((resc["state"] == "P").sum()),
                runtime_s=time.time() - t0)
    files = []
    f = out / "rescue_hits.parquet"; rs.to_parquet(f, index=False); files.append(f)
    sd = out / "sog_states"
    sd.mkdir(exist_ok=True)
    f = sd / "bucket=all.parquet"
    resc.drop(columns=["n_settings_hit", "n_settings_run"]).to_parquet(f, index=False); files.append(f)
    f = out / "rescue_summary.json"; f.write_text(json.dumps(summ, indent=1, default=int)); files.append(f)
    for f in files:
        recon.write_sidecar(f, SNAP, 0, {"script": "scripts/partB/run_rescue.py", "argv": sys.argv[1:]})
    print(json.dumps(summ, indent=1, default=int))


if __name__ == "__main__":
    main()
