#!/usr/bin/env python
"""M-CSA catalytic residues (UniProt numbering) + Riziotis et al. 2025 convergent-enzyme cases.

Downloads (to /tmp/pe_big/raw/mcsa/ and /tmp/pe_big/raw/riziotis2025/; untrusted, parsed as JSON/XML only):
  api/entries/?format=json (paginated)      entries: name, EC, reference UniProt, CATH, mechanisms
  api/residues/?format=json                  reference catalytic residues (PDB + UniProt numbering)
  api/homologues_residues.json               catalytic residues transferred to homologues (UniProt numbering)
  Europe PMC full text XML of PMC11796326    Riziotis et al. FEBS J 2025 (doi:10.1111/febs.17332) Tables 1-3

Outputs (data/annotations/):
  mcsa_entries.parquet        one row per M-CSA entry
  mcsa_residues.parquet       one row per (mcsa_id, uniprot_acc, resid); is_reference, in_atlas flag
  mcsa_convergent_cases.parquet  Riziotis 2025 analogue cases (34 cases; one row per member enzyme)
  mcsa_meta.json

Usage: .venv/bin/python -I scripts/annotations/mcsa.py --accessions data/annotations/accessions.parquet
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time

import pandas as pd
import requests

API = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/"
RAW = "/tmp/pe_big/raw/mcsa"
RAW_RIZ = "/tmp/pe_big/raw/riziotis2025"
EPMC_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11796326/fullTextXML"
CATH_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def download(url, path, force=False):
    if os.path.exists(path) and not force:
        return path
    for i in range(6):
        try:
            with requests.get(url, stream=True, timeout=600) as r:
                r.raise_for_status()
                with open(path + ".part", "wb") as fh:
                    for b in r.iter_content(1 << 20):
                        fh.write(b)
            os.replace(path + ".part", path)
            return path
        except Exception as e:  # noqa: BLE001
            print("retry", url, e, file=sys.stderr)
            time.sleep(2 ** i)
    raise RuntimeError(url)


def fetch_entries(raw):
    path = os.path.join(raw, "entries_all.json")
    if os.path.exists(path):
        return json.load(open(path))
    out, url = [], API + "entries/?format=json"
    while url:
        r = requests.get(url, timeout=300)
        r.raise_for_status()
        d = r.json()
        out += d["results"]
        url = d.get("next")
    json.dump(out, open(path, "w"))
    return out


def parse_entries(entries):
    rows = []
    for e in entries:
        ecs = e.get("all_ecs") or []
        rx = e.get("reaction") or {}
        mechs = rx.get("mechanisms") or []
        cath = sorted({rc.get("domain_cath_id") for r in e.get("residues", []) for rc in r.get("residue_chains", [])
                       if rc.get("domain_cath_id")})
        pdbs = sorted({rc.get("pdb_id") for r in e.get("residues", []) for rc in r.get("residue_chains", [])
                       if rc.get("is_reference")})
        seqs = [s.get("uniprot_id") for s in (e.get("protein") or {}).get("sequences", [])]
        rows.append({
            "mcsa_id": e["mcsa_id"],
            "enzyme_name": e.get("enzyme_name"),
            "reference_uniprot": e.get("reference_uniprot_id"),
            "is_reference_uniprot": e.get("is_reference_uniprot_id"),
            "uniprot_ids": ";".join(filter(None, seqs)),
            "ec": rx.get("ec"),
            "all_ecs": ";".join(ecs) if isinstance(ecs, list) else str(ecs),
            "reference_pdb": ";".join(pdbs),
            "cath": ";".join(cath),
            "n_residues": len(e.get("residues", [])),
            "n_mechanisms": len(mechs),
            "mechanism_rating": max([m.get("rating") or 0 for m in mechs], default=None),
            "is_detailed_mechanism": any(m.get("is_detailed") for m in mechs),
            "cofactors": ";".join(sorted({f"{c.get('name')}|CHEBI:{c.get('chebi_id')}" for c in rx.get("compounds", [])
                                          if c.get("type") == "cofactor"})),
        })
    return pd.DataFrame(rows)


def parse_residues(res, is_homologue_file):
    rows = []
    for r in res:
        chains = r.get("residue_chains") or []
        ref_chain = next((c for c in chains if c.get("is_reference")), chains[0] if chains else {})
        for s in r.get("residue_sequences") or []:
            if s.get("uniprot_id") is None or s.get("resid") is None:
                continue
            rows.append({
                "mcsa_id": r["mcsa_id"],
                "uniprot_acc": s["uniprot_id"],
                "resid": int(s["resid"]),
                "code": s.get("code"),
                "is_reference": bool(s.get("is_reference")),
                "source": "homologue" if (is_homologue_file and not s.get("is_reference")) else "reference",
                "roles_summary": r.get("roles_summary"),
                "function_location": r.get("function_location_abv") or None,
                "ptm": r.get("ptm") or None,
                "ref_pdb": ref_chain.get("pdb_id"),
                "ref_chain": ref_chain.get("chain_name"),
                "ref_auth_resid": ref_chain.get("auth_resid"),
                "ref_cath": ref_chain.get("domain_cath_id"),
            })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ Riziotis 2025 tables
def _cells(row):
    return [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", c)).strip()
            for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S)]


def parse_riziotis(xml_path):
    """Tables 1-3: one case starts with a row whose 1st cell is the enzyme name and 2nd the EC; each member
    row then has CATH | M-CSA id | PDB | residue set | ... Pairwise rows without CATH (extra superpositions
    in 3-member cases) are skipped. Returns one row per member enzyme."""
    x = open(xml_path, encoding="utf-8").read()
    cats = {"febs17332-tbl-0001": "structural+mechanistic", "febs17332-tbl-0002": "mechanistic-only",
            "febs17332-tbl-0003": "reaction-only"}
    out, case_no = [], 0
    for m in re.finditer(r'<table-wrap id="([^"]*)".*?</table-wrap>', x, re.S):
        tid = m.group(1)
        body = re.search(r"<tbody>(.*?)</tbody>", m.group(0), re.S)
        rows = re.findall(r"<tr>(.*?)</tr>", body.group(1) if body else m.group(0), re.S)
        cur = None
        for r in rows:
            c = _cells(r)
            if not c:
                continue
            if len(c) > 2 and re.match(r"^\d+\.\d+\.\d+\.[\d-]+", c[1]) and CATH_RE.match(c[2]):
                case_no += 1
                cur = {"case_id": f"R{case_no:02d}", "category": cats.get(tid, tid), "table": tid,
                       "enzyme": c[0], "ec": c[1].split(" ")[0]}
                mem = c[2:]
            elif CATH_RE.match(c[0]) and cur is not None:
                mem = c
            else:
                continue
            if len(mem) < 3 or not mem[1].isdigit():
                continue
            out.append({**cur, "cath": mem[0], "mcsa_id": int(mem[1]), "pdb": mem[2].lower(),
                        "residue_types": mem[3] if len(mem) > 3 else None,
                        "superposed_residues": mem[4] if len(mem) > 4 and re.search(r"[A-Z][a-z]{2}\d", mem[4]) else None})
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", default="data/annotations/accessions.parquet")
    ap.add_argument("--proteins", default="data/annotations/uniprot_proteins.parquet",
                    help="optional; used to map secondary accessions")
    ap.add_argument("--outdir", default="data/annotations")
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(RAW_RIZ, exist_ok=True)

    entries = fetch_entries(RAW)
    p_res = download(API + "residues/?format=json", os.path.join(RAW, "residues.json"))
    p_hom = download(API + "homologues_residues.json", os.path.join(RAW, "homologues_residues.json"))
    p_xml = download(EPMC_XML, os.path.join(RAW_RIZ, "PMC11796326_fullText.xml"))

    E = parse_entries(entries)
    R = pd.concat([parse_residues(json.load(open(p_res)), False),
                   parse_residues(json.load(open(p_hom)), True)], ignore_index=True)
    R = (R.sort_values(["mcsa_id", "uniprot_acc", "resid", "is_reference"], ascending=[True, True, True, False])
           .drop_duplicates(["mcsa_id", "uniprot_acc", "resid"]))

    acc = pd.read_parquet(a.accessions)
    wanted = set(acc.accession)
    # map UniProt accessions used by M-CSA onto atlas accessions (primary or via secondary accessions)
    m = {x: x for x in wanted}
    if os.path.exists(a.proteins):
        P = pd.read_parquet(a.proteins, columns=["accession", "uniprot_acc", "secondary_accs"])
        for ac, pa, sec in P.dropna(subset=["uniprot_acc"]).itertuples(index=False):
            m.setdefault(pa, ac)
            for s in (sec or "").split(";"):
                if s:
                    m.setdefault(s, ac)
    R["accession"] = R.uniprot_acc.map(m)
    R["in_atlas"] = R.accession.notna()
    R.to_parquet(os.path.join(a.outdir, "mcsa_residues.parquet"), index=False, compression="zstd")
    E["reference_in_atlas"] = E.reference_uniprot.fillna("").map(
        lambda v: any(x.strip() in m for x in v.split(",") if x.strip()))
    E.to_parquet(os.path.join(a.outdir, "mcsa_entries.parquet"), index=False, compression="zstd")

    C = parse_riziotis(p_xml)
    C = C.merge(E[["mcsa_id", "enzyme_name", "reference_uniprot", "reference_in_atlas"]], on="mcsa_id", how="left")
    # atlas accessions carrying M-CSA residues for this entry (reference or homologue)
    hits = R[R.in_atlas].groupby("mcsa_id").accession.nunique().rename("n_atlas_accessions_with_mcsa_residues")
    C = C.merge(hits, on="mcsa_id", how="left").fillna({"n_atlas_accessions_with_mcsa_residues": 0})
    C.to_parquet(os.path.join(a.outdir, "mcsa_convergent_cases.parquet"), index=False, compression="zstd")

    meta = {
        "source": API, "downloaded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_entries": int(len(E)), "n_residue_rows": int(len(R)),
        "n_residue_rows_in_atlas": int(R.in_atlas.sum()),
        "n_atlas_accessions_with_residues": int(R.loc[R.in_atlas, "accession"].nunique()),
        "n_atlas_accessions_reference": int(R.loc[R.in_atlas & R.is_reference, "accession"].nunique()),
        "riziotis_cases": int(C.case_id.nunique()), "riziotis_members": int(len(C)),
        "riziotis_cases_by_category": C.drop_duplicates("case_id").category.value_counts().to_dict(),
        "sha256": {os.path.basename(p): sha256(p) for p in
                   [os.path.join(RAW, "entries_all.json"), p_res, p_hom, p_xml]},
    }
    json.dump(meta, open(os.path.join(a.outdir, "mcsa_meta.json"), "w"), indent=1)
    print({k: v for k, v in meta.items() if k != "sha256"})


if __name__ == "__main__":
    main()
