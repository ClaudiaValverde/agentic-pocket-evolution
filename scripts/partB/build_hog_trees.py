#!/usr/bin/env python
"""Per-root-HOG gene trees from an Ancestra OrthoXML (DESIGN §5.1; docs/tools_ancestra.md §2-3).

Nested orthologGroup (speciation, level = TaxRange) / paralogGroup (duplication, level = taxon just below the
duplication) -> one multifurcating tree per root group; leaves = proteins. Checks every leaf species against the
NCBI guide tree and every node level against guide-tree level names. Writes
  <out>/hog_trees.parquet   one row per node: og_id, node (preorder), parent, kind (0 leaf/1 spec/2 dup), level,
                            name (group id / accession), species, taxid
  <out>/summary.json        counts (trees, nodes, duplications, polytomies, species / level mismatches)
(+ sidecars). Reload with recon.trees_from_table(pyarrow.parquet.read_table(...)).

    python scripts/partB/build_hog_trees.py --orthoxml data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz \
        --guide data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk --out results/partB/hog_trees/O1
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthoxml", default=str(ROOT / "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz"))
    ap.add_argument("--guide", default=str(ROOT / "data/orthology/ancestra_tier1_v0/tier1_guide_tree.nwk"))
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--out", default=str(ROOT / "results/partB/hog_trees/O1"))
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    st = recon.SpeciesTree.from_file(a.guide)
    guide_species = set(st.tree.leaf_names())
    guide_levels = set(st.node_of_level)
    trees, c = [], Counter()
    sp_bad, lev_bad = Counter(), Counter()
    nch_hist = Counter()
    for t in recon.iter_orthoxml_trees(a.orthoxml):
        trees.append(t)
        c["trees"] += 1
        c["nodes"] += t.n_nodes
        c["leaves"] += t.n_leaves
        c["dup_nodes"] += int((t.kind == recon.DUP).sum())
        c["spec_nodes"] += int((t.kind == recon.SPEC).sum())
        nsp = len(set(t.leaf_species()))
        c["trees_ge3_species"] += nsp >= 3
        for v in range(t.n_nodes):
            k = len(t.children(v))
            if k:
                nch_hist[min(k, 10)] += 1
            if t.kind[v] == recon.LEAF:
                if t.species[v] not in guide_species:
                    sp_bad[t.species[v]] += 1
            elif t.level[v] not in guide_levels and t.level[v] not in guide_species:
                lev_bad[t.level[v]] += 1
        # duplication nodes at the root
        c["root_is_dup"] += int(t.kind[0] == recon.DUP)
    cols = recon.trees_to_table(trees)
    tbl = pa.table({"og_id": pa.array(cols["og_id"], pa.string()), "node": pa.array(cols["node"], pa.int32()),
                    "parent": pa.array(cols["parent"], pa.int32()), "kind": pa.array(cols["kind"], pa.int8()),
                    "level": pa.array(cols["level"], pa.string()), "name": pa.array(cols["name"], pa.string()),
                    "species": pa.array(cols["species"], pa.string()), "taxid": pa.array(cols["taxid"], pa.int64())})
    f = out / "hog_trees.parquet"
    pq.write_table(tbl, f, compression="zstd")
    summ = dict(c, children_per_internal_node={str(k): v for k, v in sorted(nch_hist.items())},
                frac_internal_polytomy=float(sum(v for k, v in nch_hist.items() if k > 2) / max(1, sum(nch_hist.values()))),
                leaf_species_not_in_guide=dict(sp_bad), levels_not_in_guide=dict(lev_bad.most_common(20)),
                n_level_mismatch_nodes=int(sum(lev_bad.values())),
                orthoxml=a.orthoxml, orthoxml_sha256=recon.file_sha256(a.orthoxml), guide=a.guide,
                guide_sha256=recon.file_sha256(a.guide), n_guide_species=len(guide_species),
                runtime_s=time.time() - t0)
    s = out / "summary.json"
    s.write_text(json.dumps(summ, indent=1))
    for x in (f, s):
        recon.write_sidecar(x, a.snapshot, None, {"script": "scripts/partB/build_hog_trees.py"})
    print(json.dumps({k: v for k, v in summ.items() if k not in ("levels_not_in_guide",)}, indent=1))


if __name__ == "__main__":
    main()
