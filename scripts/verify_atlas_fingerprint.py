#!/usr/bin/env python
"""Compare a fresh `pocketcompare describe --af-mode` table with the stored atlas rows of the same model.

Applies the atlas worker filters (ok, burrow_score >= 0.03, plddt_mean >= 70; see pocketcompare/afdb/README.md
'Pocket filters') to the fresh table, joins on pocket_id, and reports per-pocket cosine, residue identity and
max |delta|. Exit code 1 if any matched cosine < 0.9999 or the kept pocket sets differ.

    python scripts/verify_atlas_fingerprint.py fresh.parquet --tier 1
"""
import argparse
import sys
from pathlib import Path

import duckdb
import numpy as np

from pocketcompare.fpindex import read_table

ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fresh", help="parquet written by `pocketcompare describe`")
    ap.add_argument("--tier", type=int, default=1)
    ap.add_argument("--atlas", default=str(ATLAS))
    ap.add_argument("--min-score", type=float, default=0.03)
    ap.add_argument("--min-plddt", type=float, default=70.0)
    a = ap.parse_args()

    meta, V, rec = read_table(a.fresh)
    meta = meta.reset_index(drop=True)
    keep = meta.ok & (meta.burrow_score >= a.min_score) & (meta.plddt_mean >= a.min_plddt)
    stems = sorted({p.split(":")[0] for p in meta.pocket_id})
    accs = [s.split("-")[1] for s in stems]
    glob = f"{a.atlas}/pockets/tier{a.tier}/*.parquet"
    files = duckdb.sql(
        f"select distinct filename from read_parquet('{glob}', union_by_name=true, filename=true) "
        f"where accession in ({','.join(repr(x) for x in accs)})").fetchall()
    rows = []
    for (f,) in files:
        m2, V2, rec2 = read_table(f)
        sel = m2.pocket_id.str.split(":").str[0].isin(stems).to_numpy()
        rows.append((m2[sel].reset_index(drop=True), V2[sel]))
    if not rows:
        print("no atlas rows found"); return 1
    am = rows[0][0] if len(rows) == 1 else __import__("pandas").concat([r[0] for r in rows], ignore_index=True)
    AV = np.concatenate([r[1] for r in rows])
    print(f"method fresh={rec.get('method')} atlas={rec2.get('method')} cfg_sha equal={rec.get('config_sha256') == rec2.get('config_sha256')}")
    fresh_ids = set(meta.pocket_id[keep]); atlas_ids = set(am.pocket_id)
    print(f"fresh pockets all={len(meta)} kept_by_filters={len(fresh_ids)}; atlas rows={len(atlas_ids)}; "
          f"only_fresh={sorted(fresh_ids - atlas_ids)} only_atlas={sorted(atlas_ids - fresh_ids)}")
    bad = fresh_ids != atlas_ids
    fi = {p: i for i, p in enumerate(meta.pocket_id)}
    for j, p in enumerate(am.pocket_id):
        if p not in fi:
            continue
        i = fi[p]
        v, w = V[i].astype(np.float64), AV[j].astype(np.float64)
        cos = float(v @ w / (np.linalg.norm(v) * np.linalg.norm(w)))
        same_res = meta.residues[i] == am.residues[j]
        print(f"{p:32s} score={am.burrow_score[j]:.3f} cos={cos:.7f} max|d|={np.abs(v - w).max():.2e} residues_equal={same_res}")
        bad |= (cos < 0.9999) or not same_res
    return int(bad)


if __name__ == "__main__":
    sys.exit(main())
