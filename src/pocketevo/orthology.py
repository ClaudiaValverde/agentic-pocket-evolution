"""Compact HOG tables from an Ancestra OrthoXML, and ortholog/paralog queries on them.

Tables (written by scripts/orthology/export_hog_tables.py; any OrthoXML 0.5 with nested orthologGroup /
paralogGroup works):
  hogs.parquet          one row per group (HOG): hog (int32, preorder id), hog_name (OrthoXML id), parent (int32,
                        -1 for a root HOG), root (int32), kind ('S' = orthologGroup/speciation, 'D' =
                        paralogGroup/duplication), level (taxon name of the group: TaxRange for S; for D the
                        lineage in which the duplication happened, i.e. Ancestra's taxonId = level just below the
                        duplication), level_taxid (NCBI taxid, -1 if unknown), depth (root = 0), n_genes
  protein_hogs.parquet  one row per input protein: protein, species (OrthoXML species name), root (-1 if the
                        protein is in no group = singleton), path (list<int32>: HOG ids root -> deepest group)

Relation of two proteins a, b (QfO flattening, as Ancestra's hogs_to_pairs): let g = deepest group containing
both. g is an orthologGroup -> 'ortholog'; a paralogGroup -> 'paralog'; no common group -> 'unrelated' (different
root HOGs or singletons); proteins absent from the tables -> 'unknown'. The returned level is g's level.
"""
from __future__ import annotations

import gzip
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass

NS = "{http://orthoXML.org/2011/}"


@dataclass
class ParsedOrthoXML:
    hogs: list           # rows (hog, hog_name, parent, root, kind, level, depth, n_genes)
    proteins: list       # rows (protein, species, root, path)
    taxa: dict           # orthoxml taxon id -> name
    taxon_parent: dict   # orthoxml taxon id -> parent id


def parse_orthoxml(path: str) -> ParsedOrthoXML:
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rb") as fh:
        root = ET.parse(fh).getroot()
    gene_prot, gene_sp = {}, {}
    for sp in root.findall(f"{NS}species"):
        name = sp.get("name")
        for g in sp.iter(f"{NS}gene"):
            gene_prot[g.get("id")] = g.get("protId")
            gene_sp[g.get("id")] = name
    taxa, tparent = {}, {}
    tx = root.find(f"{NS}taxonomy")
    if tx is not None:
        def walk_t(el, par):
            for c in el.findall(f"{NS}taxon"):
                taxa[c.get("id")] = c.get("name")
                tparent[c.get("id")] = par
                walk_t(c, c.get("id"))
        walk_t(tx, None)
    hogs, path_of = [], {}

    def level_of(el):
        for p in el.findall(f"{NS}property"):
            if p.get("name") == "TaxRange":
                return p.get("value")
        t = el.get("taxonId")
        return taxa.get(t, t) if t is not None else ""

    def walk(el, parent, rootid, depth, path):
        hid = len(hogs)
        kind = "S" if el.tag == f"{NS}orthologGroup" else "D"
        row = [hid, el.get("id") or "", parent, hid if parent < 0 else rootid, kind, level_of(el), depth, 0]
        hogs.append(row)
        rootid = row[3]
        path = path + [hid]
        n = 0
        for c in el:
            if c.tag == f"{NS}geneRef":
                path_of[c.get("id")] = path
                n += 1
            elif c.tag in (f"{NS}orthologGroup", f"{NS}paralogGroup"):
                n += walk(c, hid, rootid, depth + 1, path)
        row[7] = n
        return n

    groups = root.find(f"{NS}groups")
    if groups is not None:
        for el in groups:
            if el.tag in (f"{NS}orthologGroup", f"{NS}paralogGroup"):
                walk(el, -1, -1, 0, [])
    proteins = []
    for gid, prot in gene_prot.items():
        p = path_of.get(gid)
        proteins.append((prot, gene_sp[gid], p[0] if p else -1, p or []))
    return ParsedOrthoXML(hogs, proteins, taxa, tparent)


def to_frames(px: ParsedOrthoXML, level_taxid: dict | None = None):
    import pandas as pd
    h = pd.DataFrame(px.hogs, columns=["hog", "hog_name", "parent", "root", "kind", "level", "depth", "n_genes"])
    h["level_taxid"] = h.level.map(level_taxid or {}).fillna(-1).astype("int64")
    for c in ("hog", "parent", "root", "depth", "n_genes"):
        h[c] = h[c].astype("int32")
    p = pd.DataFrame(px.proteins, columns=["protein", "species", "root", "path"])
    p["root"] = p.root.astype("int32")
    return h, p


def write_tables(px: ParsedOrthoXML, out_dir: str, level_taxid: dict | None = None):
    import pyarrow as pa
    import pyarrow.parquet as pq
    h, p = to_frames(px, level_taxid)
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(h, preserve_index=False), os.path.join(out_dir, "hogs.parquet"),
                   compression="zstd")
    tbl = pa.table({"protein": p.protein.values, "species": p.species.values,
                    "root": pa.array(p.root.values, pa.int32()),
                    "path": pa.array(p.path.tolist(), pa.list_(pa.int32()))})
    pq.write_table(tbl, os.path.join(out_dir, "protein_hogs.parquet"), compression="zstd")
    return h, p


class HogIndex:
    """In-memory index over hogs.parquet + protein_hogs.parquet (or a ParsedOrthoXML).

    alias: optional {other_accession: hog_member_accession} (e.g. the ID-gap mapping or placements), applied
    before lookup."""

    def __init__(self, hogs, proteins, alias: dict | None = None):
        self.kind = dict(zip(hogs.hog, hogs.kind))
        self.level = dict(zip(hogs.hog, hogs.level))
        self.level_taxid = dict(zip(hogs.hog, hogs.level_taxid)) if "level_taxid" in hogs else {}
        self.hog_name = dict(zip(hogs.hog, hogs.hog_name))
        self.path = {a: list(p) for a, p in zip(proteins.protein, proteins.path)}
        self.species = dict(zip(proteins.protein, proteins.species))
        self.alias = alias or {}

    @classmethod
    def load(cls, table_dir: str, alias: dict | None = None):
        import pandas as pd
        h = pd.read_parquet(os.path.join(table_dir, "hogs.parquet"))
        p = pd.read_parquet(os.path.join(table_dir, "protein_hogs.parquet"))
        return cls(h, p, alias)

    @classmethod
    def from_orthoxml(cls, path: str, level_taxid: dict | None = None):
        h, p = to_frames(parse_orthoxml(path), level_taxid)
        return cls(h, p)

    def _resolve(self, a):
        return self.alias.get(a, a)

    def root(self, a):
        p = self.path.get(self._resolve(a))
        return p[0] if p else None

    def relation(self, a: str, b: str) -> dict:
        """-> {'relation': ortholog|paralog|unrelated|same|unknown, 'level', 'level_taxid', 'mrca_hog'}"""
        ra, rb = self._resolve(a), self._resolve(b)
        out = {"relation": "unknown", "level": None, "level_taxid": None, "mrca_hog": None}
        if ra not in self.path or rb not in self.path:
            return out
        if ra == rb:
            out["relation"] = "same"
            return out
        pa_, pb_ = self.path[ra], self.path[rb]
        g = None
        for x, y in zip(pa_, pb_):
            if x != y:
                break
            g = x
        if g is None:
            out["relation"] = "unrelated"
            return out
        out.update(relation="ortholog" if self.kind[g] == "S" else "paralog", level=self.level[g],
                   level_taxid=self.level_taxid.get(g), mrca_hog=self.hog_name.get(g) or int(g))
        return out

    def hogs_of(self, a: str) -> list[tuple]:
        """Nested HOGs of a protein, root -> deepest: [(hog, hog_name, kind, level, level_taxid)]."""
        return [(g, self.hog_name[g], self.kind[g], self.level[g], self.level_taxid.get(g))
                for g in self.path.get(self._resolve(a), [])]


def relation(a: str, b: str, table_dir: str | None = None, _cache: dict = {}) -> dict:  # noqa: B006
    """Convenience: relation of two proteins using the O2 tables (cached index).
    table_dir defaults to <repo>/data/orthology/ancestra_O2/accurate (or $POCKETEVO_HOG_TABLES)."""
    if table_dir is None:
        here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        table_dir = os.environ.get("POCKETEVO_HOG_TABLES",
                                   os.path.join(here, "data", "orthology", "ancestra_O2", "tables"))
    if table_dir not in _cache:
        alias = {}
        mp = os.path.join(os.path.dirname(table_dir), "idgap_mapping.parquet")
        if os.path.exists(mp):
            import pandas as pd
            m = pd.read_parquet(mp, columns=["atlas_accession", "hog_member_accession", "method"])
            m = m[m.method.isin(["identical_seq", "merged_into_member", "gene_name+seq", "seq"])]
            alias = dict(zip(m.atlas_accession, m.hog_member_accession))
        _cache[table_dir] = HogIndex.load(table_dir, alias)
    return _cache[table_dir].relation(a, b)
