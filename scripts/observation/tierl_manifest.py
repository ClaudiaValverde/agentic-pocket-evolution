#!/usr/bin/env python
"""Freeze a snapshot-format manifest of the tier-L (local pockets) run, verify it, and write
data/tierL/<id>/{manifest.parquet, manifest.json, excluded.parquet, tier3_overlap.parquet, verify_fp.json}.

The manifest has exactly the columns of data/snapshot/<id>/manifest.parquet (scripts/atlas_snapshot.py
scan_unit) plus tier_label = 'L', proteome, taxid. `tier` stays the int code 4 used inside the tables (int8 in the
atlas schema and in src/pocketevo/sites.py); paths are relative to atlas_root = the tier-L root
(pockets/tierL/<unit>.parquet).

Checks:
  * every worklist unit has a done marker + both parquet files, rows match the done json (else excluded)
  * method record hash of every pocket file == the atlas method hash (pocketevo.observation.ATLAS_METHOD_HASH)
  * 20 random structures (seeded) re-described from the Foldcomp entry with `pocketcompare describe --af-mode
    --backend ref` and the atlas filters -> same pocket ids, residues and fingerprints (cos >= 0.99999)
  * tier3_overlap.parquet: tier-L accessions that are also in the main snapshot's tier-3 files (downstream should
    read them from ONE source; both are the same Foldcomp v4 entry + method, so they must be identical: checked on
    the overlap of one unit)

  python scripts/observation/tierl_manifest.py --root /tmp/pe_big/tierL --id tierL_<UTC>
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / "src"))
sys.path.insert(0, str(P / "scripts"))
from pocketevo import observation as O  # noqa: E402
from atlas_snapshot import scan_unit  # noqa: E402

ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")
FC_DB = str(ATLAS / "foldcomp/afdb_uniprot_v4")
SNAP = "20261030T085136Z"


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def verify_fingerprints(root: Path, man: pd.DataFrame, n: int, seed: int, work: Path) -> dict:
    rng = np.random.default_rng(seed)
    S = pd.concat([pd.read_parquet(root / p, columns=["structure_stem", "accession", "status", "n_pockets_kept"])
                   .assign(unit_id=u) for u, p in zip(man.unit_id, man.structures_path)], ignore_index=True)
    S = S[(S.status == "ok") & (S.n_pockets_kept > 0)]
    pick = S.iloc[rng.choice(len(S), size=min(n, len(S)), replace=False)]
    ent = pd.read_parquet(root / "fc_entries_tierL.parquet", columns=["accession", "fragment", "version", "off", "len"])
    ent["stem"] = "AF-" + ent.accession + "-F" + ent.fragment.astype(str) + "-model_v" + ent.version.astype(str)
    ent = ent.set_index("stem")
    work.mkdir(parents=True, exist_ok=True)
    paths = []
    for s in pick.structure_stem:
        name, pdb = O.foldcomp_read(FC_DB, int(ent.at[s, "off"]), int(ent.at[s, "len"]))
        f = work / f"{s}.pdb"
        f.write_text(pdb)
        paths.append(str(f))
    out = work / "verify.parquet"
    env = dict(os.environ, OMP_NUM_THREADS="1", NUMBA_NUM_THREADS="1")
    subprocess.run([str(P / ".venv/bin/pocketcompare"), "describe", *paths, "--af-mode", "--backend", "ref",
                    "--workers", "4", "-o", str(out)], check=True, env=env, capture_output=True)
    fresh = pq.read_table(out).to_pandas()
    fresh = fresh[fresh.ok & (fresh.burrow_score >= 0.03) & (fresh.plddt_mean.fillna(-1) >= 70)]
    stored = []
    for u, g in pick.groupby("unit_id"):
        t = pq.read_table(root / man.set_index("unit_id").at[u, "pockets_path"],
                          columns=["pocket_id", "accession", "residues", "fingerprint"]).to_pandas()
        stored.append(t[t.accession.isin(g.accession)])
    stored = pd.concat(stored, ignore_index=True)
    a, b = set(fresh.pocket_id), set(stored.pocket_id)
    fr = fresh.set_index("pocket_id")
    st = stored.set_index("pocket_id")
    cos, maxd, res_eq = [], [], 0
    for pid in sorted(a & b):
        x = np.asarray(fr.at[pid, "fingerprint"], np.float32)
        y = np.asarray(st.at[pid, "fingerprint"], np.float32)
        cos.append(float(x @ y / np.linalg.norm(x) / np.linalg.norm(y)))
        maxd.append(float(np.abs(x - y).max()))
        res_eq += fr.at[pid, "residues"] == st.at[pid, "residues"]
    return dict(n_structures=int(len(pick)), n_fresh_kept=len(a), n_stored=len(b), only_fresh=sorted(a - b),
                only_stored=sorted(b - a), n_common=len(a & b), min_cos=min(cos) if cos else None,
                max_abs_diff=max(maxd) if maxd else None, residues_equal=int(res_eq),
                ok=bool(a == b and cos and min(cos) >= 0.99999 and res_eq == len(a & b)),
                accessions=sorted(pick.accession), backend="ref (atlas was written with fast; identical by design)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/tmp/pe_big/tierL")
    ap.add_argument("--id", default=None)
    ap.add_argument("--n-verify", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    root = Path(a.root).resolve()
    sid = a.id or "tierL_" + dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outd = P / "data/tierL" / sid
    outd.mkdir(parents=True, exist_ok=False)
    wl = pd.read_parquet(root / "worklist.parquet")
    recs, exc = [], []
    for r in wl.itertuples():
        dp = root / "state/done" / f"{r.unit_id}.json"
        if not dp.exists():
            exc.append(dict(unit_id=r.unit_id, tier=int(r.tier), excluded_reason="not_done"))
            continue
        rec = scan_unit(str(root), int(r.unit_id), str(dp))
        if rec.get("excluded_reason"):
            exc.append(dict(unit_id=r.unit_id, tier=int(r.tier), excluded_reason=rec["excluded_reason"]))
            continue
        for k in ("pockets", "structures"):
            rec[f"{k}_path"] = rec[f"{k}_path"].replace("/tier4/", "/tierL/")
            assert (root / rec[f"{k}_path"]).exists()
        rec.update(tier_label="L", proteome=r.proteome, taxid=int(r.taxid), organism=r.organism)
        recs.append(rec)
    man = pd.DataFrame(recs).drop(columns=["excluded_reason"])
    snap_cols = list(pd.read_parquet(P / f"data/snapshot/{SNAP}/manifest.parquet").columns)
    man = man[snap_cols + ["tier_label", "proteome", "taxid", "organism"]]
    hashes = {u: O.parquet_method_hash(root / p) for u, p in zip(man.unit_id, man.pockets_path)}
    man["method_record_sha256"] = man.unit_id.map(hashes)
    bad_hash = int((man.method_record_sha256 != O.ATLAS_METHOD_HASH).sum())
    man.to_parquet(outd / "manifest.parquet", index=False)
    pd.DataFrame(exc, columns=["unit_id", "tier", "excluded_reason"]).to_parquet(outd / "excluded.parquet", index=False)
    # tier-3 overlap (same accession also in a tier-3 unit of the main snapshot)
    acc = pd.concat([pd.read_parquet(root / p, columns=["accession"]) for p in man.structures_path]).accession
    smain = pd.read_parquet(P / f"data/snapshot/{SNAP}/manifest.parquet")
    s3 = smain[smain.tier == 3]
    accset = pa.array(sorted(set(acc)))
    import pyarrow.compute as pc
    ov = []
    for u, p in zip(s3.unit_id, s3.structures_path):
        t = pq.read_table(ATLAS / p, columns=["accession", "taxid"])
        t = t.filter(pc.is_in(t["accession"], accset))
        if t.num_rows:
            ov.append(t.to_pandas().assign(tier3_unit_id=u))
    ov = pd.concat(ov, ignore_index=True) if ov else pd.DataFrame(columns=["accession", "taxid", "tier3_unit_id"])
    ov.to_parquet(outd / "tier3_overlap.parquet", index=False)
    # identity check on the overlap: one tier-3 unit's pockets vs tier-L pockets for the same accessions
    ident = None
    if len(ov):
        u3 = int(ov.tier3_unit_id.value_counts().index[0])
        accs = set(ov[ov.tier3_unit_id == u3].accession)
        t3 = pq.read_table(ATLAS / s3.set_index("unit_id").at[u3, "pockets_path"],
                           columns=["pocket_id", "accession", "residues", "fingerprint"]).to_pandas()
        t3 = t3[t3.accession.isin(accs)].set_index("pocket_id")
        tl = pd.concat([pq.read_table(root / p, columns=["pocket_id", "accession", "residues", "fingerprint"]).to_pandas()
                        for p in man.pockets_path])
        tl = tl[tl.accession.isin(accs)].set_index("pocket_id")
        common = sorted(set(t3.index) & set(tl.index))
        d = [float(np.abs(np.asarray(t3.at[i, "fingerprint"], np.float32) - np.asarray(tl.at[i, "fingerprint"], np.float32)).max()) for i in common]
        ident = dict(tier3_unit=u3, n_accessions=len(accs), n_t3=len(t3), n_tierL=len(tl), n_common=len(common),
                     max_abs_diff=max(d) if d else None,
                     residues_equal=int(sum(t3.at[i, "residues"] == tl.at[i, "residues"] for i in common)))
    ver = verify_fingerprints(root, man, a.n_verify, a.seed, root / "verify")
    (outd / "verify_fp.json").write_text(json.dumps(ver, indent=1))
    done = [json.loads((root / "state/done" / f"{u}.json").read_text()) for u in man.unit_id]
    t_first = min(d["t"] - d["wall_s"] for d in done)
    t_last = max(d["t"] for d in done)
    summ = dict(snapshot_id=sid, created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), atlas_root=str(root),
                tier_code=4, tier_label="L", manifest="manifest.parquet",
                manifest_sha256=sha256(outd / "manifest.parquet"),
                worklist_sha256=sha256(root / "worklist.parquet"),
                fc_entries_sha256=sha256(root / "fc_entries_tierL.parquet"),
                n_units=int(len(man)), n_excluded=len(exc),
                excluded_reasons=pd.Series([e["excluded_reason"] for e in exc]).value_counts().to_dict(),
                method_config_sha256=sorted(man.method_config_sha256.unique()),
                method_record_sha256=sorted(man.method_record_sha256.unique()), n_bad_method_hash=bad_hash,
                per_tier={"L": dict(units=int(len(man)), structures=int(man.structures_rows.sum()),
                                    pockets=int(man.pockets_rows.sum()),
                                    structures_ok=int(man.json_structures_ok.sum()),
                                    structures_error=int(man.json_structures_error.sum()),
                                    pockets_all=int(man.json_pockets_all.sum()),
                                    pockets_bytes=int(man.pockets_size.sum()),
                                    rows_mismatch=int((~man.rows_match_json).sum()),
                                    worklist_units=int(len(wl)), proteomes=int(man.proteome.nunique()))},
                run_wall_s=round(t_last - t_first, 1),
                run_structures_per_s=round(float(man.structures_rows.sum()) / max(1.0, t_last - t_first), 1),
                tier3_overlap_accessions=int(ov.accession.nunique()), tier3_overlap_identity=ident,
                verify_fingerprints=dict(ok=ver["ok"], n_structures=ver["n_structures"], n_common=ver["n_common"],
                                         min_cos=ver["min_cos"]),
                script=str(Path(__file__).resolve()), argv=sys.argv, runtime_s=round(time.time() - t0, 1))
    (outd / "manifest.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
