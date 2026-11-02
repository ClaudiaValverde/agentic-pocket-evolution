"""NCBI taxonomy helpers (taxdump nodes.dmp / names.dmp / merged.dmp). No external dependencies.

Taxonomy(dir).lineage(t) -> [1, ..., t] (root first), after merged.dmp remapping.
"""
from __future__ import annotations

import os


class Taxonomy:
    def __init__(self, taxdump_dir: str, load_names: bool = True):
        self.dir = taxdump_dir
        self.parent: dict[int, int] = {}
        self.rank: dict[int, str] = {}
        with open(os.path.join(taxdump_dir, "nodes.dmp")) as fh:
            for line in fh:
                f = line.split("\t|\t", 3)
                t = int(f[0])
                self.parent[t] = int(f[1])
                self.rank[t] = f[2]
        self.merged: dict[int, int] = {}
        mp = os.path.join(taxdump_dir, "merged.dmp")
        if os.path.exists(mp):
            with open(mp) as fh:
                for line in fh:
                    f = line.split("\t|\t")
                    self.merged[int(f[0])] = int(f[1].split("\t")[0])
        self.name: dict[int, str] = {}
        if load_names:
            with open(os.path.join(taxdump_dir, "names.dmp")) as fh:
                for line in fh:
                    if "scientific name" not in line:
                        continue
                    f = line.split("\t|\t")
                    if f[3].startswith("scientific name"):
                        self.name[int(f[0])] = f[1]
        self._lin: dict[int, list[int]] = {}

    def resolve(self, t: int) -> int:
        return self.merged.get(int(t), int(t))

    def lineage(self, t: int) -> list[int]:
        t = self.resolve(t)
        if t in self._lin:
            return self._lin[t]
        if t not in self.parent:
            return []
        path = [t]
        while path[-1] != 1 and self.parent[path[-1]] != path[-1]:
            path.append(self.parent[path[-1]])
        lin = path[::-1]
        self._lin[t] = lin
        return lin

    def at_rank(self, t: int, rank: str):
        for x in self.lineage(t):
            if self.rank.get(x) == rank:
                return x
        return None

    def name_at_rank(self, t: int, rank: str) -> str:
        x = self.at_rank(t, rank)
        return self.name.get(x, "") if x is not None else ""

    def is_descendant(self, t: int, anc: int) -> bool:
        return anc in set(self.lineage(t))

    def lineage_names(self, t: int, skip_root: bool = True) -> str:
        lin = self.lineage(t)
        if skip_root:
            lin = [x for x in lin if x not in (1, 131567)]
        return ";".join(self.name.get(x, str(x)) for x in lin)


def taxdump_info(taxdump_dir: str) -> dict:
    """mtime + md5 of the taxdump files (for DATA.md / provenance)."""
    import hashlib
    import time
    out = {}
    for f in ("nodes.dmp", "names.dmp", "merged.dmp", "taxdump.tar.gz"):
        p = os.path.join(taxdump_dir, f)
        if os.path.exists(p):
            h = hashlib.md5()
            with open(p, "rb") as fh:
                for b in iter(lambda: fh.read(1 << 22), b""):
                    h.update(b)
            out[f] = {"md5": h.hexdigest(), "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(os.path.getmtime(p)))}
    return out
