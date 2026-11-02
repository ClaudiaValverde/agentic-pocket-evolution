#!/usr/bin/env python
"""Write the toy fixtures used by the web app's orthogroup and case views (schemas: webapp/SCHEMA.md).

Everything here is synthetic (ids start with TOY); it exists only so that the views and tests work end-to-end
before the pipeline writes real results. Deterministic: rerunning produces identical files.

  .venv/bin/python -I webapp/fixtures/make_fixtures.py
"""
from __future__ import annotations

import json
import math
import zlib
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
SOG = HERE / "sogs"
CASE = HERE / "cases" / "TOY_gain_01"

# toy species (real NCBI taxids so the tree looks familiar; states are synthetic)
SPECIES = [(9606, "Homo sapiens", "HUMAN", "Mammalia"), (10090, "Mus musculus", "MOUSE", "Mammalia"),
           (7955, "Danio rerio", "DANRE", "Actinopterygii"), (7227, "Drosophila melanogaster", "DROME", "Insecta"),
           (6239, "Caenorhabditis elegans", "CAEEL", "Nematoda"), (559292, "Saccharomyces cerevisiae", "YEAST", "Fungi"),
           (83333, "Escherichia coli K-12", "ECOLI", "Bacteria")]
TREE = "((((((9606,10090)Mammalia,7955)Vertebrata,(7227,6239)Ecdysozoa)Bilateria,559292)Opisthokonta)Eukaryota,83333)root;"


def bucket(key: str) -> str:
    return f"{zlib.crc32(key.encode()) & 0xff:02x}"


def write(t: pa.Table, p: Path):
    p.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(t, p, compression="zstd")


def sogs():
    (SOG / "sog_states").mkdir(parents=True, exist_ok=True)
    for f in (SOG / "sog_states").glob("*.parquet"):
        f.unlink()
    (SOG / "manifest.json").write_text(json.dumps({
        "schema": "sogs-web-v1", "snapshot": "TOY", "orthology": "TOY", "tau": 0.3,
        "created_utc": "1970-01-01T00:00:00Z", "code_sha256": None, "toy": True}, indent=1) + "\n")
    (SOG / "species_tree.nwk").write_text(TREE + "\n")
    write(pa.table({"taxid": pa.array([s[0] for s in SPECIES], pa.int32()), "name": [s[1] for s in SPECIES],
                    "code": [s[2] for s in SPECIES], "clade": [s[3] for s in SPECIES]}), SOG / "species.parquet")
    ogs = [("TOYOG1", "toy kinase family", "TOYHS1", ""), ("TOYOG2", "toy dehydrogenase family", "TOYHS2",
                                                            "orthology-uncertain")]
    # per OG: list of SOGs; per SOG: states per species (paralogs as a list)
    design = {
        "TOYOG1": [("ATP site", "none", None, {9606: "P", 10090: "P", 7955: "P", 7227: "P", 6239: "PW",
                                                559292: "P", 83333: "G"}),
                   ("allosteric site (toy gain in Vertebrata)", "gain", "Vertebrata",
                    {9606: "P", 10090: "P", 7955: "P", 7227: "A", 6239: "A", 559292: "U", 83333: "G"})],
        "TOYOG2": [("cofactor site", "none", None, {9606: "P", 10090: "P", 7955: "W", 7227: "P", 6239: "U",
                                                     559292: "P", 83333: "P"}),
                   ("lid pocket (toy loss in Mammalia)", "loss", "Mammalia",
                    {9606: "A", 10090: "A", 7955: "P", 7227: "P", 6239: "P", 559292: "U", 83333: "A"})]}
    og_rows, sog_rows, st_rows = [], [], []
    for og_id, label, ref, flags in ogs:
        members = set()
        for k, (slabel, event, node, states) in enumerate(design[og_id]):
            sog_id = f"{og_id}:g{k}"
            cnt = {s: 0 for s in "PWAUG"}
            for taxid, sts in states.items():
                for j, s in enumerate(sts):
                    cnt[s] += 1
                    code = next(x[2] for x in SPECIES if x[0] == taxid)
                    acc = None if s == "G" else f"TOY{og_id[-1]}{code[:3]}{j}"
                    if acc:
                        members.add(acc)
                    st_rows.append({"og_id": og_id, "sog_id": sog_id, "taxid": taxid, "accession": acc, "state": s,
                                    "site_id": f"{acc}:s{k}" if s in "PW" else None,
                                    "pocket_id": f"AF-{acc}-F1-model_v6:m{k}" if s in "PW" else None,
                                    "jaccard": {"P": 0.62, "W": 0.41}.get(s),
                                    "footprint_plddt": None if s == "G" else {"U": 61.5}.get(s, 93.0),
                                    "reason": {"U": "plddt<90", "W": "rescue:2/6", "A": ""}.get(s, "")})
            sog_rows.append({"og_id": og_id, "sog_id": sog_id, "label": slabel, "footprint_columns": [10 + k, 11 + k, 40],
                             "footprint_ref": f"{ref}:{12 + k} {13 + k} {44}", **{f"n_{s}": cnt[s] for s in "PWAUG"},
                             "event": event, "event_node": node, "score": 12.5 if event != "none" else None,
                             "case_id": "TOY_gain_01" if (og_id, k) == ("TOYOG1", 1) else None})
        og_rows.append({"og_id": og_id, "label": label, "n_species": len(SPECIES), "n_members": len(members),
                        "n_sogs": len(design[og_id]), "ref_accession": ref, "flags": flags})
    write(pa.Table.from_pylist(og_rows, schema=pa.schema([
        ("og_id", pa.string()), ("label", pa.string()), ("n_species", pa.int32()), ("n_members", pa.int32()),
        ("n_sogs", pa.int32()), ("ref_accession", pa.string()), ("flags", pa.string())])), SOG / "orthogroups.parquet")
    write(pa.Table.from_pylist(sog_rows, schema=pa.schema([
        ("og_id", pa.string()), ("sog_id", pa.string()), ("label", pa.string()),
        ("footprint_columns", pa.list_(pa.int32())), ("footprint_ref", pa.string()),
        *[(f"n_{s}", pa.int32()) for s in "PWAUG"], ("event", pa.string()), ("event_node", pa.string()),
        ("score", pa.float64()), ("case_id", pa.string())])), SOG / "sogs.parquet")
    st_schema = pa.schema([("og_id", pa.string()), ("sog_id", pa.string()), ("taxid", pa.int32()),
                           ("accession", pa.string()), ("state", pa.string()), ("site_id", pa.string()),
                           ("pocket_id", pa.string()), ("jaccard", pa.float32()), ("footprint_plddt", pa.float32()),
                           ("reason", pa.string())])
    by_b = {}
    for r in st_rows:
        by_b.setdefault(bucket(r["og_id"]), []).append(r)
    for b, rows in by_b.items():
        rows.sort(key=lambda r: (r["og_id"], r["sog_id"], r["taxid"]))
        write(pa.Table.from_pylist(rows, schema=st_schema), SOG / "sog_states" / f"bucket={b}.parquet")


def toy_model_cif(ligand_dz: float = 0.0) -> str:
    """A 12-residue poly-Ala helix plus one toy ligand atom set, as mmCIF (pLDDT in B_iso)."""
    import gemmi
    lines = []
    n = 0
    for i in range(12):
        ang = math.radians(100 * i)
        for name, el, r, dz in (("N", "N", 1.55, -0.6), ("CA", "C", 2.3, 0.0), ("C", "C", 1.65, 0.7),
                                ("O", "O", 1.9, 1.6), ("CB", "C", 3.3, 0.3)):
            n += 1
            x, y, z = r * math.cos(ang), r * math.sin(ang), 1.5 * i + dz
            lines.append(f"ATOM  {n:5d}  {name:<3s} ALA A{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{80 + i:6.2f}"
                         f"           {el}")
    for j, (x, y, z) in enumerate(((0.0, 0.0, 6.0 + ligand_dz), (0.5, 0.8, 7.2 + ligand_dz), (-0.6, 0.7, 8.4 + ligand_dz), (0.1, -0.4, 9.5 + ligand_dz))):
        n += 1
        lines.append(f"HETATM{n:5d}  C{j + 1:<2d} LIG B   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 85.00           C")
    st = gemmi.read_pdb_string("\n".join(lines) + "\nEND\n")
    st.name = "TOY_cofold_m1"
    st.setup_entities()
    return "# TOY co-folding model (synthetic helix + ligand), webapp fixture\n" + st.make_mmcif_document().as_string()


CASE_YAML = """\
schema: case-v1
case_id: TOY_gain_01
title: "TOY: gain of an allosteric pocket in vertebrates"
category: gain
track: discovery
label: supported
toy: true
summary: >
  Synthetic fixture. The allosteric toy SOG is present in three vertebrates and absent-observed in fly and worm;
  parsimony places one gain on the Vertebrata stem. Nothing in this case is a real result.
og_id: TOYOG1
sog_id: "TOYOG1:g1"
score: 12.5
proteins:
  - {accession: TOY1HUM0, taxid: 9606, species: Homo sapiens, role: derived, state: P, site_id: "TOY1HUM0:s1", residues: "A:13 A:14 A:44"}
  - {accession: TOY1DAN0, taxid: 7955, species: Danio rerio, role: derived, state: P, site_id: "TOY1DAN0:s1", residues: "A:13 A:14 A:44"}
  - {accession: TOY1DRO0, taxid: 7227, species: Drosophila melanogaster, role: sister, state: A, residues: ""}
  - {accession: P24941, taxid: 9606, species: Homo sapiens, role: control, state: P, residues: "", note: "real atlas protein, link test only"}
filters:
  - {name: persistence, status: pass, detail: "P in 3/3 vertebrates (toy)"}
  - {name: experimental, status: pending, detail: "no evidence (toy)"}
  - {name: computational, status: pass, detail: "toy co-folding: ligand pLDDT 85 vs decoy 55"}
artifact_checks:
  - {name: footprint_plddt, status: pass, detail: "A members footprint pLDDT 93"}
  - {name: truncation, status: pass, detail: "lengths 0.9-1.1x OG median"}
  - {name: rescue, status: pending, detail: ""}
evidence:
  - {type: other, id: TOY-EVIDENCE-1, about_change: false, detail: "placeholder record", file: evidence/toy_record.json, filter: persistence}
cofold:
  - {model_id: m1, accession: TOY1HUM0, ligand: LIG, seed: 1, role: derived, file: cofold/toy_m1.cif,
     metrics: {ligand_plddt: 85.0, min_pae: 1.4, centroid_dist_A: 1.2, posebusters_pass: true}}
  - {model_id: m2, accession: TOY1DRO0, ligand: LIG, seed: 1, role: ancestral, file: cofold/toy_m2.cif,
     metrics: {ligand_plddt: 52.0, min_pae: 9.8, centroid_dist_A: 14.0, posebusters_pass: false}}
figures:
  - {file: figures/toy_tree.svg, caption: "Toy SOG states on the species tree"}
literature:
  - {doi: "10.0000/toy.0001", title: "Placeholder reference (toy)", note: "not a real DOI"}
notes: "Generated by webapp/fixtures/make_fixtures.py."
"""

TOY_SVG = """<svg xmlns="http://www.w3.org/2000/svg" width="320" height="140" viewBox="0 0 320 140">
<rect width="320" height="140" fill="#fff"/>
<g stroke="#333" fill="none"><path d="M10 70 H60 V30 H120 M60 70 V110 H120 M120 30 V15 H180 M120 30 V45 H180"/></g>
<g font-family="sans-serif" font-size="11"><text x="186" y="19">human</text><text x="186" y="49">zebrafish</text>
<text x="126" y="114">fly</text></g>
<circle cx="250" cy="15" r="6" fill="#1b7837"/><circle cx="250" cy="45" r="6" fill="#1b7837"/>
<circle cx="190" cy="110" r="6" fill="#d73027"/><text x="10" y="135" font-size="10" font-family="sans-serif">TOY figure</text>
</svg>
"""


def case():
    for d in ("figures", "cofold", "evidence"):
        (CASE / d).mkdir(parents=True, exist_ok=True)
    (CASE / "case.yaml").write_text(CASE_YAML)
    (CASE / "figures" / "toy_tree.svg").write_text(TOY_SVG)
    (CASE / "cofold" / "toy_m1.cif").write_text(toy_model_cif())
    (CASE / "cofold" / "toy_m2.cif").write_text(toy_model_cif(14.0))
    (CASE / "evidence" / "toy_record.json").write_text(json.dumps({"toy": True, "record": "placeholder"}) + "\n")


if __name__ == "__main__":
    sogs()
    case()
    print("fixtures written to", HERE)
