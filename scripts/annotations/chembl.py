#!/usr/bin/env python
"""ChEMBL bioactivity evidence per atlas accession: number of distinct compounds with pChEMBL >= 6.

Raw (EBI FTP ChEMBLdb/latest; to /tmp/pe_big/raw/chembl/): chembl_<v>_sqlite.tar.gz (sha256 checked against
checksums.txt), extracted to /tmp/pe_big/work/chembl/ (~30 GB). The full SQLite is used instead of the REST API
because distinct-compound counts need all activity rows (millions), which is slow and incomplete via paging.

Counting rules (pChEMBL is only defined for '=' relations of standardised IC50/EC50/Ki/Kd/... values):
  activities.pchembl_value >= 6, data_validity_comment is null or 'Manually validated', potential_duplicate = 0.
  * n_cpd_pchembl6_single  : target_type = 'SINGLE PROTEIN' and assay confidence_score = 9 (direct single protein)
  * n_cpd_pchembl6_any     : any target containing the protein as component (complexes, families, ...)
Output: data/annotations/chembl_targets.parquet (accession x ChEMBL target) and chembl_by_accession.parquet
(one row per atlas accession with any ChEMBL target) + chembl_meta.json.

Usage: .venv/bin/python -I scripts/annotations/chembl.py --accessions data/annotations/accessions.parquet
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import sqlite3
import subprocess
import tarfile

import pandas as pd

FTP = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
RAW = "/tmp/pe_big/raw/chembl"
WORK = "/tmp/pe_big/work/chembl"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", default="data/annotations/accessions.parquet")
    ap.add_argument("--proteins", default="data/annotations/uniprot_proteins.parquet")
    ap.add_argument("--outdir", default="data/annotations")
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    ck = os.path.join(RAW, "checksums.txt")
    if not os.path.exists(ck):
        subprocess.run(["curl", "-sSf", "-o", ck, FTP + "checksums.txt"], check=True)
    sums = {l.split()[1]: l.split()[0] for l in open(ck) if l.strip() and not l.startswith("SHA256")}
    tgz_name = next(f for f in sums if f.endswith("_sqlite.tar.gz"))
    version = tgz_name.split("_sqlite")[0]
    tgz = os.path.join(RAW, tgz_name)
    if not os.path.exists(tgz):
        subprocess.run(["curl", "-sSf", "-o", tgz, FTP + tgz_name], check=True)
    digest = sha256(tgz)
    assert digest == sums[tgz_name], (digest, sums[tgz_name])
    dbs = glob.glob(os.path.join(WORK, "**", "*.db"), recursive=True)
    if not dbs:
        with tarfile.open(tgz) as tf:
            members = [m for m in tf.getmembers() if m.isfile() and m.name.endswith(".db")
                       and not os.path.isabs(m.name) and ".." not in m.name.split("/")]
            tf.extractall(WORK, members=members, filter="data")
        dbs = glob.glob(os.path.join(WORK, "**", "*.db"), recursive=True)
    db = dbs[0]
    tmpd = os.path.join(WORK, "sqlite_tmp")   # /var/tmp is a 64 MB tmpfs in the container
    os.makedirs(tmpd, exist_ok=True)
    os.environ["SQLITE_TMPDIR"] = tmpd
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.execute(f"PRAGMA temp_store_directory = '{tmpd}'")
    con.execute("PRAGMA cache_size = -4000000")
    q = """
    select cs.accession as uniprot_acc, td.chembl_id as target_chembl_id, td.pref_name as target_name,
           td.target_type, td.tax_id as target_tax_id,
           count(distinct case when ass.confidence_score = 9 then act.molregno end) as n_cpd_pchembl6_conf9,
           count(distinct act.molregno) as n_cpd_pchembl6,
           count(*) as n_act_pchembl6,
           max(act.pchembl_value) as max_pchembl
    from component_sequences cs
    join target_components tc on tc.component_id = cs.component_id
    join target_dictionary td on td.tid = tc.tid
    join assays ass on ass.tid = td.tid
    join activities act on act.assay_id = ass.assay_id
    where act.pchembl_value >= 6
      and (act.data_validity_comment is null or act.data_validity_comment = 'Manually validated')
      and act.potential_duplicate = 0
      and cs.accession is not null
    group by cs.accession, td.chembl_id
    """
    T = pd.read_sql_query(q, con)
    # targets with no potent compound still matter (target exists) -> list all targets of our proteins
    allt = pd.read_sql_query("""select cs.accession as uniprot_acc, td.chembl_id as target_chembl_id,
                                td.pref_name as target_name, td.target_type, td.tax_id as target_tax_id
                                from component_sequences cs join target_components tc on tc.component_id=cs.component_id
                                join target_dictionary td on td.tid=tc.tid where cs.accession is not null""", con)
    T = allt.merge(T, on=["uniprot_acc", "target_chembl_id", "target_name", "target_type", "target_tax_id"],
                   how="left").fillna({"n_cpd_pchembl6_conf9": 0, "n_cpd_pchembl6": 0, "n_act_pchembl6": 0})
    # per-accession distinct compounds (single protein, conf 9) — computed in SQL to avoid double counting
    per = pd.read_sql_query("""
    select cs.accession as uniprot_acc,
           count(distinct case when td.target_type='SINGLE PROTEIN' and ass.confidence_score=9 then act.molregno end)
               as n_cpd_pchembl6_single,
           count(distinct act.molregno) as n_cpd_pchembl6_any
    from component_sequences cs
    join target_components tc on tc.component_id = cs.component_id
    join target_dictionary td on td.tid = tc.tid
    join assays ass on ass.tid = td.tid
    join activities act on act.assay_id = ass.assay_id
    where act.pchembl_value >= 6
      and (act.data_validity_comment is null or act.data_validity_comment = 'Manually validated')
      and act.potential_duplicate = 0 and cs.accession is not null
    group by cs.accession""", con)

    acc = pd.read_parquet(a.accessions, columns=["accession"]).drop_duplicates()
    amap = {x: x for x in acc.accession}
    if os.path.exists(a.proteins):
        P = pd.read_parquet(a.proteins, columns=["accession", "uniprot_acc", "secondary_accs"]).dropna(subset=["uniprot_acc"])
        for ac, pa, sec in P.itertuples(index=False):
            amap.setdefault(pa, ac)
            for s in (sec or "").split(";"):
                if s:
                    amap.setdefault(s, ac)
    for D in (T, per):
        D["accession"] = D.uniprot_acc.map(amap)
    T = T[T.accession.notna()].copy()
    per = per[per.accession.notna()].copy()
    for c in ("n_cpd_pchembl6_conf9", "n_cpd_pchembl6", "n_act_pchembl6"):
        T[c] = T[c].astype("int64")
    T["chembl_version"] = version
    T.to_parquet(os.path.join(a.outdir, "chembl_targets.parquet"), index=False, compression="zstd")
    tgt = T.groupby("accession").agg(n_chembl_targets=("target_chembl_id", "nunique"),
                                     single_protein_targets=("target_chembl_id",
                                         lambda s: ";".join(sorted(s[T.loc[s.index, "target_type"] == "SINGLE PROTEIN"]))),
                                     max_pchembl=("max_pchembl", "max")).reset_index()
    by = tgt.merge(per.drop(columns=["uniprot_acc"]).groupby("accession").max().reset_index(), on="accession", how="left")
    by = by.fillna({"n_cpd_pchembl6_single": 0, "n_cpd_pchembl6_any": 0})
    by["n_cpd_pchembl6_single"] = by.n_cpd_pchembl6_single.astype("int64")
    by["n_cpd_pchembl6_any"] = by.n_cpd_pchembl6_any.astype("int64")
    by["chembl_version"] = version
    by.to_parquet(os.path.join(a.outdir, "chembl_by_accession.parquet"), index=False, compression="zstd")
    meta = {"source": FTP + tgz_name, "chembl_version": version, "sha256": digest, "license": "CC BY-SA 3.0",
            "n_accessions_with_target": int(len(by)),
            "n_accessions_pchembl6_single_ge1": int((by.n_cpd_pchembl6_single >= 1).sum()),
            "n_accessions_pchembl6_any_ge1": int((by.n_cpd_pchembl6_any >= 1).sum()),
            "n_accessions_total": int(len(acc))}
    json.dump(meta, open(os.path.join(a.outdir, "chembl_meta.json"), "w"), indent=1)
    print(meta)


if __name__ == "__main__":
    main()
