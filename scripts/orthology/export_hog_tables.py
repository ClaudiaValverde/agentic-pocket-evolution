#!/usr/bin/env python
"""Export compact tables from an Ancestra run (hogs.orthoxml[.gz]) for downstream use.

Writes into <out>/tables/:
  hogs.parquet, protein_hogs.parquet   (schema: src/pocketevo/orthology.py)
  levels.tsv                           guide-tree node name -> NCBI taxid, rank, parent, #leaves
  guide_tree_taxids.nwk                guide tree with every node (leaves and internal) labelled by NCBI taxid
  stats.json                           counts (root HOGs, groups, proteins in groups, per species fraction)
Run with the ancestra env (needs ancestra.speciestree for the taxid labels):
  PYTHONPATH=src envs/ancestra/bin/python scripts/orthology/export_hog_tables.py \
      --orthoxml data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz --tree data/orthology/ancestra_O2/o2_guide_tree.nwk \
      --taxdump /tmp/pe_big/ancestra/taxdump --out data/orthology/ancestra_O2
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

P = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(P, "src"))
from pocketevo.orthology import parse_orthoxml, write_tables  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthoxml", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--taxdump", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from ancestra import newick, speciestree
    t0 = time.time()
    guide = newick.read_newick(a.tree) if hasattr(newick, "read_newick") else newick.parse(open(a.tree).read())
    leaves = guide.leaf_names()
    leaf_taxid = {l: int(l.rsplit("_", 1)[1]) for l in leaves}
    tax = speciestree.load_taxdump(a.taxdump)
    nt = speciestree.ncbi_tree(leaf_taxid, tax=tax)
    # the guide tree was built by the same function: check identical node names
    gnames = {n.name for n in guide.preorder()}
    level_taxid, rows = {}, []
    for n in nt.preorder():
        t = int(n.data.get("taxid", -1))
        level_taxid[n.name] = t
        rows.append((n.name, t, n.data.get("rank", ""), n.parent.name if n.parent is not None else "",
                     len(n.leaf_names()), n.name in gnames))
    missing = gnames - set(level_taxid)
    if missing:
        raise SystemExit(f"guide-tree nodes not reproduced from taxdump: {sorted(missing)[:10]}")
    tdir = os.path.join(a.out, "tables")
    os.makedirs(tdir, exist_ok=True)
    with open(os.path.join(tdir, "levels.tsv"), "w") as fh:
        fh.write("level\ttaxid\trank\tparent\tn_leaves\tin_guide_tree\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    def nwk(n):
        lab = str(level_taxid[n.name])
        return ("(" + ",".join(nwk(c) for c in n.children) + ")" + lab) if n.children else lab
    with open(os.path.join(tdir, "guide_tree_taxids.nwk"), "w") as fh:
        fh.write(nwk(guide) + ";\n")

    px = parse_orthoxml(a.orthoxml)
    h, p = write_tables(px, tdir, level_taxid)
    unk = sorted(set(h.level[h.level_taxid < 0]))
    p["in_group"] = p.root >= 0
    per_sp = p.groupby("species").in_group.mean()
    st = dict(orthoxml=a.orthoxml, n_proteins=int(len(p)), n_in_groups=int(p.in_group.sum()),
              n_root_hogs=int((h.parent == -1).sum()), n_groups=int(len(h)),
              n_ortholog_groups=int((h.kind == "S").sum()), n_paralog_groups=int((h.kind == "D").sum()),
              max_root_size=int(h[h.parent == -1].n_genes.max()),
              frac_in_groups_per_species=dict(median=float(per_sp.median()), min=float(per_sp.min()),
                                              lowest=per_sp.nsmallest(8).round(3).to_dict()),
              levels_without_taxid=unk[:20], seconds=round(time.time() - t0, 1))
    json.dump(st, open(os.path.join(tdir, "stats.json"), "w"), indent=1)
    print(json.dumps(st, indent=1))


if __name__ == "__main__":
    main()
