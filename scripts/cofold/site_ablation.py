#!/usr/bin/env python3
"""Site-ablation variants: known binding-site residues -> Gly (or Phe).

The site is defined on the *holo PDB* structure: every polymer residue with a
heavy atom within --cutoff (default 4.5 A) of the crystal ligand. Those residues
are mapped onto the UniProt sequence by pairwise alignment of the PDB SEQRES-
observed chain to the UniProt sequence, then substituted.

Gly ablation (Masters 2025 style) removes the side chains that line the pocket;
the Phe variant instead fills the cavity with bulk. Both keep the sequence length,
so the cached AF3 MSA can be reused with the query row replaced, which holds the
MSA *constant* and makes the comparison purely sequence-driven (the hardest form
of the control). Use --msa colabfold for an independent MSA instead.

Usage
  site_ablation.py --system cdk2 --mode gly --out-fasta /tmp/x.fa --out-sites /tmp/x.json
  site_ablation.py --all --out-dir /tmp/pe_big/ablation
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import gemmi
import numpy as np

PDB = "/cluster/scratch/jjaenes/pdb"

# system -> (uniprot acc, pdb id, chain, ligand resname, ligand seqnum or None, cofactor ccds)
SYSTEMS = {
    "cdk2":  ("P24941", "1hck", "A", "ATP", 400, ["CCD:ATP", "CCD:MG"]),
    "ca2":   ("P00918", "3hs4", "A", "AZM", None, ["CCD:ZN", "CCD:AZM"]),
    "pygl":  ("P06737", "1fa9", "A", "PLP", 999, ["CCD:PLP"]),
    "yadh1": ("P00330", "4w6z", "A", "8ID", 403, ["CCD:NAD", "CCD:ZN", "CCD:ZN"]),
}


def heavy(res):
    return [a for a in res if a.element.name != "H"]


def uniprot_seq(acc: str) -> str:
    import urllib.request
    url = f"https://rest.uniprot.org/uniprotkb/{acc.upper()}.fasta"
    with urllib.request.urlopen(url, timeout=30) as r:
        txt = r.read().decode()
    return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">")).upper()


def site_residues(pdb: str, chain: str, lig: str, num: int | None, cutoff: float):
    """Return (polymer residues, ligand residue, indices of site residues in polymer)."""
    st = gemmi.read_structure(f"{PDB}/{pdb[1:3]}/{pdb}.cif.gz")
    st.remove_hydrogens()
    ch = st[0][chain]
    cands = [r for r in ch if r.name == lig and (num is None or r.seqid.num == num)]
    if not cands:
        raise SystemExit(f"{pdb}:{chain} has no {lig} {num}")
    if num is None:
        zn = [r for r in ch if r.name == "ZN"]
        if zn:
            zp = zn[0][0].pos
            cands.sort(key=lambda r: min(a.pos.dist(zp) for a in r))
    ligres = cands[0]
    lpos = np.array([a.pos.tolist() for a in heavy(ligres)])
    poly = [r for r in ch if r.het_flag != "H" and gemmi.find_tabulated_residue(r.name)
            and gemmi.find_tabulated_residue(r.name).is_amino_acid()]
    site = []
    for i, r in enumerate(poly):
        xyz = np.array([a.pos.tolist() for a in heavy(r)])
        if np.min(np.linalg.norm(xyz[:, None] - lpos[None], axis=-1)) < cutoff:
            site.append(i)
    return poly, ligres, site


def map_to_uniprot(poly, useq: str) -> dict[int, int]:
    """polymer index -> 0-based index in the UniProt sequence."""
    s1 = gemmi.one_letter_code([r.name for r in poly])
    aln = gemmi.align_string_sequences(list(s1), list(useq), [])
    m, i, j = {}, 0, 0
    for a, b in zip(aln.add_gaps(s1, 1), aln.add_gaps(useq, 2)):
        if a != "-" and b != "-":
            if a.upper() == b.upper():
                m[i] = j
        i += a != "-"
        j += b != "-"
    return m


def build(system: str, cutoff: float, seq: str | None = None) -> dict:
    acc, pdb, chain, lig, num, ligands = SYSTEMS[system]
    useq = seq or uniprot_seq(acc)
    poly, ligres, site = site_residues(pdb, chain, lig, num, cutoff)
    m = map_to_uniprot(poly, useq)
    mapped, unmapped = [], []
    for i in site:
        (mapped if i in m else unmapped).append(i)
    pos = sorted(m[i] for i in mapped)          # 0-based UniProt positions
    wt = "".join(useq[p] for p in pos)
    out = {"system": system, "acc": acc, "pdb": pdb, "chain": chain, "ligand": lig,
           "ligand_seqnum": num, "cutoff_A": cutoff, "ligands_for_cofold": ligands,
           "length": len(useq), "n_site_pdb": len(site), "n_site_mapped": len(pos),
           "n_site_unmapped": len(unmapped),
           "site_uniprot_1based": [p + 1 for p in pos],
           "site_wt_residues": [f"{useq[p]}{p+1}" for p in pos],
           "wt_seq": useq}
    for mode, aa in (("gly", "G"), ("phe", "F")):
        s = list(useq)
        n = 0
        for p in pos:
            if s[p] != aa:
                s[p] = aa
                n += 1
        out[f"{mode}_seq"] = "".join(s)
        out[f"n_sub_{mode}"] = n
    out["site_wt_string"] = wt
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", choices=sorted(SYSTEMS))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--out-dir", default="/tmp/pe_big/ablation")
    a = ap.parse_args()
    systems = sorted(SYSTEMS) if a.all else ([a.system] if a.system else [])
    if not systems:
        ap.error("give --system or --all")
    od = Path(a.out_dir)
    od.mkdir(parents=True, exist_ok=True)
    for s in systems:
        d = build(s, a.cutoff)
        (od / f"{s}_site.json").write_text(json.dumps(d, indent=1))
        for mode in ("wt", "gly", "phe"):
            (od / f"{s}_{mode}.fa").write_text(f">{s}_{mode}\n{d[mode + '_seq']}\n")
        print(f"{s} {d['acc']} {d['pdb']}:{d['ligand']}  len={d['length']}  "
              f"site={d['n_site_mapped']}/{d['n_site_pdb']} mapped "
              f"(unmapped {d['n_site_unmapped']})  nsub gly={d['n_sub_gly']} phe={d['n_sub_phe']}")
        print("   " + " ".join(d["site_wt_residues"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
