"""Read-only access to the pocket atlas: accession index, pocket rows, model structures.

Nothing is written to disk. Structures are read straight from the AFDB v6 tars (tiers 1-2, by member offset)
or the Foldcomp DB (tier 3, by entry offset) and kept only in an in-memory LRU.
"""
from __future__ import annotations

import gzip
import re
import threading
import zlib
from collections import OrderedDict
from functools import lru_cache
from pathlib import Path

import pyarrow.parquet as pq

from .config import Config

N_BUCKETS = 256
_AFID = re.compile(r"^AF-([A-Za-z0-9]+)-F\d+(?:-model_v\d+)?$")
TIERS = (1, 2, 3)
POCKET_COLS = ["pocket_id", "accession", "fragment", "tier", "unit_id", "model_version", "merged_index",
               "burrow_rank", "burrow_score", "n_points", "volume_A3", "site_volume_A3", "center_x", "center_y",
               "center_z", "residues", "plddt_mean", "plddt_min", "plddt_frac_lt70"]


def bucket_of(acc: str) -> str:
    """Must match scripts/webapp_build_index.py."""
    return f"{zlib.crc32(acc.encode()) & (N_BUCKETS - 1):02x}"


def fragment_offset(fragment: int) -> int:
    """UniProt position = author resnum + 200*(k-1) (verified, results/sites/<snap>/SITES.md)."""
    return 200 * (int(fragment) - 1)


class LRU:
    def __init__(self, n):
        self.n, self.d, self.lock = n, OrderedDict(), threading.Lock()

    def get(self, k):
        with self.lock:
            if k in self.d:
                self.d.move_to_end(k)
                return self.d[k]
        return None

    def put(self, k, v):
        with self.lock:
            self.d[k] = v
            self.d.move_to_end(k)
            while len(self.d) > self.n:
                self.d.popitem(last=False)


class Atlas:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.struct_cache = LRU(cfg.structure_cache)
        self._genes = None
        self._genes_lock = threading.Lock()

    # ---------- index ----------
    def tier_dirs(self):
        """{tier: dir} for every tier whose index exists (first index dir wins)."""
        out = {}
        for d in self.cfg.index_dirs:
            for t in TIERS:
                p = Path(d) / f"tier{t}"
                if t not in out and p.is_dir() and any(p.glob("bucket=*.parquet")):
                    out[t] = p
        return out

    def lookup(self, acc: str) -> list[dict]:
        """All atlas structures of an accession (any tier, any fragment), lowest tier first."""
        acc = acc.strip().upper()
        rows = []
        for t, d in sorted(self.tier_dirs().items()):
            f = d / f"bucket={bucket_of(acc)}.parquet"
            if f.exists():
                rows += _read_bucket(str(f), acc)
        return sorted(rows, key=lambda r: (r["tier"], r["fragment"]))

    def status(self):
        return {"snapshot": self.cfg.snapshot, "atlas_root": str(self.cfg.atlas_root),
                "index_tiers": {t: str(p) for t, p in self.tier_dirs().items()},
                "gene_table": str(self.gene_table()) if self.gene_table() else None}

    # ---------- genes ----------
    def gene_table(self):
        for p in self.cfg.annotations:
            if Path(p).exists():
                return Path(p)
        return None

    def genes(self):
        """{GENE_UPPER: [accession,...]} and {accession: (gene, protein_name)}, loaded lazily once."""
        with self._genes_lock:
            if self._genes is None:
                g2a, a2g = {}, {}
                p = self.gene_table()
                if p is not None:
                    names = pq.read_schema(p).names
                    gcol = next((c for c in ("gene", "gene_primary", "gene_name") if c in names), None)
                    ncol = next((c for c in ("protein_name", "name") if c in names), None)
                    if gcol and "accession" in names:
                        df = pq.read_table(p, columns=["accession", gcol] + ([ncol] if ncol else [])).to_pandas()
                        for r in df.itertuples(index=False):
                            gene = getattr(r, gcol)
                            pname = getattr(r, ncol) if ncol else None
                            a2g[r.accession] = (gene, pname)
                            if isinstance(gene, str) and gene:
                                g2a.setdefault(gene.upper(), []).append(r.accession)
                self._genes = (g2a, a2g)
            return self._genes

    def search(self, q: str, limit: int = 100) -> dict:
        """Accession (any tier; also `AF-<acc>-F<k>-model_v<n>`), else gene name (exact, then prefix), else
        protein-name substring. `hits` = [{accession, gene, protein_name, match}] for the gene/name cases."""
        q = q.strip()
        if not q:
            return {"query": q, "accessions": [], "by": None, "hits": []}
        m = _AFID.match(q)
        qa = (m.group(1) if m else q).upper()
        if self.lookup(qa):
            return {"query": q, "accessions": [qa], "by": "accession", "hits": []}
        g2a, a2g = self.genes()
        gene_search = self.gene_table() is not None
        ql = q.upper()
        accs = sorted(set(g2a.get(ql, [])))
        by = "gene" if accs else None
        hits = [(a, "gene") for a in accs]
        if not accs and gene_search and len(q) >= 2:
            pre = [(a, "gene prefix") for g, aa in g2a.items() if g.startswith(ql) for a in aa]
            hits = sorted(set(pre))[:limit]
            if not hits and len(q) >= 3:
                qlow = q.lower()
                hits = sorted({(a, "protein name") for a, (_, pn) in a2g.items()
                               if isinstance(pn, str) and qlow in pn.lower()})[:limit]
            by = hits[0][1] if hits else None
        return {"query": q, "accessions": [a for a, _ in hits], "by": by, "gene_search": gene_search,
                "hits": [{"accession": a, "gene": a2g.get(a, (None, None))[0],
                          "protein_name": a2g.get(a, (None, None))[1], "match": how} for a, how in hits]}

    # ---------- pockets ----------
    def pockets(self, acc: str) -> list[dict]:
        acc = acc.strip().upper()
        out = []
        for s in self.lookup(acc):
            f = self.cfg.atlas_root / "pockets" / f"tier{s['tier']}" / f"{s['unit_id']:06d}.parquet"
            out += [r for r in _read_pockets(str(f), acc) if r["fragment"] == s["fragment"]]
        for r in out:
            off = fragment_offset(r["fragment"])
            r["uniprot_positions"] = " ".join(str(int(x.split(":")[1]) + off) for x in r["residues"].split())
        return sorted(out, key=lambda r: (r["tier"], r["fragment"], -r["burrow_score"]))

    # ---------- structures ----------
    def structure(self, acc: str, tier: int | None = None, fragment: int | None = None):
        """Return (mmCIF text, structure record). Raises KeyError if not in the atlas index."""
        rows = self.lookup(acc)
        if tier is not None:
            rows = [r for r in rows if r["tier"] == int(tier)]
        if fragment is not None:
            rows = [r for r in rows if r["fragment"] == int(fragment)]
        if not rows:
            raise KeyError(acc)
        s = rows[0]
        key = (s["blob_kind"], s["blob_file"], s["blob_offset"], s["blob_size"])
        txt = self.struct_cache.get(key)
        if txt is None:
            txt = self._read_blob(s)
            self.struct_cache.put(key, txt)
        return txt, s

    def _read_blob(self, s) -> str:
        if s["blob_offset"] is None:
            raise KeyError(f"no source offset for {s['accession']}")
        if s["blob_kind"] == "tar_cif_gz":
            path = self.cfg.afdb_v6_dir / s["blob_file"]
            with open(path, "rb") as f:
                f.seek(int(s["blob_offset"]))
                raw = f.read(int(s["blob_size"]))
            return gzip.decompress(raw).decode()
        if s["blob_kind"] == "foldcomp":
            import foldcomp
            import gemmi
            path = self.cfg.foldcomp_dir / s["blob_file"]
            with open(path, "rb") as f:
                f.seek(int(s["blob_offset"]))
                raw = f.read(int(s["blob_size"]))
            name, pdb = foldcomp.decompress(raw)
            st = gemmi.read_pdb_string(pdb)
            st.name = name.replace(".cif.gz", "").replace(".pdb", "")
            st.setup_entities()
            doc = st.make_mmcif_document()
            return (f"# Converted by pocket-evolution webapp from Foldcomp {s['blob_file']} entry {name} "
                    f"(gemmi {gemmi.__version__}); B_iso_or_equiv = pLDDT\n" + doc.as_string())
        raise ValueError(s["blob_kind"])


@lru_cache(maxsize=4096)
def _read_bucket_cached(path: str, acc: str):
    # shards are sorted by accession with 64k-row groups: the filter prunes row groups by min/max stats
    return tuple(pq.read_table(path, filters=[("accession", "=", acc)]).to_pylist())


def _read_bucket(path: str, acc: str) -> list[dict]:
    return [dict(r) for r in _read_bucket_cached(path, acc)]


@lru_cache(maxsize=256)
def _read_pockets_cached(path: str, acc: str):
    names = pq.read_schema(path).names
    cols = [c for c in POCKET_COLS if c in names]
    return tuple(pq.read_table(path, columns=cols, filters=[("accession", "=", acc)]).to_pylist())


def _read_pockets(path: str, acc: str) -> list[dict]:
    return [dict(r) for r in _read_pockets_cached(path, acc)]
