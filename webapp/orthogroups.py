"""Orthogroup (SOG x species) data access, schema `sogs-web-v1` (webapp/SCHEMA.md §2)."""
from __future__ import annotations

import json
import threading
import zlib
from collections import OrderedDict
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

STATE_ORDER = "PWAUG"  # strongest first
STATE_COLORS = {"P": "#1b7837", "W": "#a6dba0", "A": "#d73027", "U": "#bdbdbd", "G": "#ffffff"}
STATE_NAMES = {"P": "present", "W": "weak-present", "A": "absent-observed", "U": "uninformative",
               "G": "gene absent"}


def bucket(key: str) -> str:
    return f"{zlib.crc32(key.encode()) & 0xff:02x}"


# ---------------- Newick ----------------
def parse_newick(s: str) -> dict:
    """Minimal Newick parser -> nested {"name", "children"}; branch lengths are dropped."""
    s = s.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        n = {"name": "", "children": []}
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                n["children"].append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"bad newick at {pos}")
        start = pos
        while pos < len(s) and s[pos] not in ",()":
            pos += 1
        lab = s[start:pos].split(":")[0].strip().strip("'\"")
        n["name"] = lab
        return n

    t = node()
    return t


def prune(t: dict, keep: set) -> dict | None:
    if not t["children"]:
        return t if t["name"] in keep else None
    ch = [c for c in (prune(c, keep) for c in t["children"]) if c is not None]
    if not ch:
        return None
    if len(ch) == 1 and not t["name"]:
        return ch[0]
    return {"name": t["name"], "children": ch}


def leaves(t: dict) -> list[str]:
    return [t["name"]] if not t["children"] else [x for c in t["children"] for x in leaves(c)]


def layout(t: dict, dx=14, dy=18):
    """Cladogram coordinates: leaves at max depth. Returns (segments, leaf_y {name: y}, width, internal labels)."""
    def depth(n):
        return 0 if not n["children"] else 1 + max(depth(c) for c in n["children"])
    D = depth(t)
    segs, leaf_y, labels = [], {}, []
    counter = [0]

    def rec(n, d):
        if not n["children"]:
            y = counter[0] * dy + dy / 2
            counter[0] += 1
            leaf_y[n["name"]] = y
            return D * dx, y
        pts = [rec(c, d + 1) for c in n["children"]]
        x = d * dx
        ys = [p[1] for p in pts]
        for cx, cy in pts:
            segs.append((x, cy, cx, cy))
        segs.append((x, min(ys), x, max(ys)))
        y = (min(ys) + max(ys)) / 2
        if n["name"]:
            labels.append((x, y, n["name"]))
        return x, y

    rec(t, 0)
    return segs, leaf_y, D * dx, labels


SOG_COMPACT = ["og_id", "sog_id", "label", "n_P", "n_W", "n_A", "n_U", "n_G", "event", "event_node", "score",
               "case_id", "footprint_size", "n_species_sites", "footprint_ref"]
STATE_COLS = ["og_id", "sog_id", "taxid", "accession", "state", "site_id", "pocket_id", "jaccard",
              "footprint_plddt", "reason", "state_s70", "model_accession", "model_tier"]


def _cols(path, wanted):
    have = set(pq.read_schema(path).names)
    return [c for c in wanted if c in have]


class Version:
    """One orthology version (a `.../web/` directory). Everything is read lazily:
    - species/tree/orthogroups.parquet (small) on first use,
    - sogs.parquet only as a column subset (no footprint_columns), on the first SOG query,
    - sog_states: only the one shard `bucket=<crc32(og_id)>` of the requested OG, filtered to that OG (LRU)."""

    def __init__(self, vid, d, manifest):
        self.id, self.dir, self.manifest = vid, Path(d), manifest
        self.toy = bool(manifest.get("toy"))
        self.sig = Version.signature(d)
        self._lock = threading.RLock()
        self._c = {}
        self._states = OrderedDict()

    @staticmethod
    def signature(d):
        d = Path(d)
        return tuple((d / f).stat().st_mtime_ns if (d / f).exists() else 0
                     for f in ("manifest.json", "orthogroups.parquet", "sogs.parquet"))

    def _once(self, key, fn):
        with self._lock:
            if key not in self._c:
                self._c[key] = fn()
            return self._c[key]

    def species(self):
        return self._once("species", lambda: {r["taxid"]: r for r in
                                              pq.read_table(self.dir / "species.parquet").to_pylist()})

    def tree(self):
        return self._once("tree", lambda: parse_newick((self.dir / "species_tree.nwk").read_text()))

    def ogs(self) -> pa.Table:
        return self._once("ogs", lambda: pq.read_table(self.dir / "orthogroups.parquet"))

    def sogs(self) -> pa.Table:
        p = self.dir / "sogs.parquet"
        return self._once("sogs", lambda: pq.read_table(p, columns=_cols(p, SOG_COMPACT)))

    def loaded_rows(self):
        return {k: v.num_rows for k, v in self._c.items() if isinstance(v, pa.Table)}

    # ---- queries
    def og(self, og_id):
        t = self.ogs()
        m = pc.equal(t["og_id"], og_id)
        r = t.filter(m)
        return r.to_pylist()[0] if r.num_rows else None

    def og_sogs(self, og_id) -> list[dict]:
        t = self.sogs()
        rows = t.filter(pc.equal(t["og_id"], og_id)).to_pylist()
        return sorted(rows, key=lambda r: _sog_k(r["sog_id"]))

    def sog(self, sog_id):
        og_id = sog_id.rsplit(":", 1)[0]
        return next((r for r in self.og_sogs(og_id) if r["sog_id"] == sog_id), None)

    def _read_states(self, og_id, sog_ids=None, columns=None):
        f = self.dir / "sog_states" / f"bucket={bucket(og_id)}.parquet"
        if not f.exists() or (sog_ids is not None and not len(sog_ids)):
            return None
        filt = [("og_id", "=", og_id)]
        if sog_ids is not None:
            filt.append(("sog_id", "in", list(sog_ids)))
        return pq.read_table(f, columns=_cols(f, columns or STATE_COLS), filters=filt)

    def og_taxa(self, og_id) -> set:
        """Taxids that have at least one row in the OG (drives the pruned species tree); cached."""
        key = ("taxa", og_id)
        with self._lock:
            if key in self._states:
                return self._states[key]
        t = self._read_states(og_id, columns=["og_id", "taxid"])
        taxa = set(t["taxid"].unique().to_pylist()) if t is not None else set()
        self._remember(key, taxa)
        return taxa

    def _remember(self, key, val):
        with self._lock:
            self._states[key] = val
            while len(self._states) > 64:
                self._states.popitem(last=False)

    def og_states(self, og_id, sog_ids=None) -> list[dict]:
        """Member rows of the OG, optionally only those of `sog_ids` (a page of SOG columns). A big OG can have
        about 10^6 rows, so the web pages always pass the SOG ids of the page."""
        key = (og_id, tuple(sorted(sog_ids)) if sog_ids is not None else None)
        with self._lock:
            if key in self._states:
                self._states.move_to_end(key)
                return self._states[key]
        t = self._read_states(og_id, sog_ids)
        rows = t.to_pylist() if t is not None else []
        self._remember(key, rows)
        return rows

    def og_page(self, q="", flags="", half="", sort="og_id", desc=False, page=1, size=50):
        """Filtered, sorted, paginated orthogroup list -> (rows, total)."""
        t = self.ogs()
        m = None
        if q:
            ql = q.lower()
            ms = [pc.match_substring(pc.utf8_lower(t[c]), ql) for c in ("og_id", "label", "ref_accession")
                  if c in t.column_names]
            m = ms[0]
            for x in ms[1:]:
                m = pc.or_(m, x)
        if flags:
            fm = pc.match_substring(pc.fill_null(t["flags"], ""), flags)
            m = fm if m is None else pc.and_(m, fm)
        if half and "half" in t.column_names:
            hm = pc.equal(t["half"], half)
            m = hm if m is None else pc.and_(m, hm)
        if m is not None:
            t = t.filter(pc.fill_null(m, False))
        total = t.num_rows
        if sort in t.column_names and total:
            idx = pc.sort_indices(t, sort_keys=[(sort, "descending" if desc else "ascending")])
            t = t.take(idx.slice((page - 1) * size, size))
        else:
            t = t.slice((page - 1) * size, size)
        return t.to_pylist(), total

    def sog_page(self, q="", event="", case_only=False, min_score=None, sort="sog_id", desc=False, page=1, size=50):
        """Filtered, sorted, paginated SOG list over all SOGs of the version -> (rows, total). Rows have no states."""
        t = self.sogs()
        m = None

        def add(x):
            nonlocal m
            m = x if m is None else pc.and_(m, x)
        if q:
            add(pc.or_(pc.match_substring(pc.utf8_lower(t["sog_id"]), q.lower()),
                       pc.match_substring(pc.utf8_lower(pc.fill_null(t["label"], "")), q.lower())))
        if event == "!none":
            add(pc.not_equal(t["event"], "none"))
        elif event:
            add(pc.equal(t["event"], event))
        if case_only:
            add(pc.is_valid(t["case_id"]))
        if min_score is not None:
            add(pc.greater_equal(t["score"], min_score))
        if m is not None:
            t = t.filter(pc.fill_null(m, False))
        total = t.num_rows
        lo = (page - 1) * size
        if sort in t.column_names and sort != "sog_id" and total:
            idx = pc.sort_indices(t, sort_keys=[(sort, "descending" if desc else "ascending")])
            t = t.take(idx.slice(lo, size))
        else:
            t = t.slice(lo, size)
        return t.drop_columns([c for c in ("footprint_ref",) if c in t.column_names]).to_pylist(), total

    def event_counts(self) -> dict:
        t = self.sogs()
        return {r["values"]: r["counts"] for r in pc.value_counts(t["event"]).to_pylist()}

    def n_sogs(self):
        return pq.ParquetFile(self.dir / "sogs.parquet").metadata.num_rows

    def n_ogs(self):
        return pq.ParquetFile(self.dir / "orthogroups.parquet").metadata.num_rows


def _sog_k(sog_id):
    try:
        return int(sog_id.rsplit(":g", 1)[1])
    except (IndexError, ValueError):
        return 0


class Orthogroups:
    """Registry of orthology versions. `dirs` are `web/` directories (or fixture dirs)."""

    def __init__(self, dirs, index_dirs=()):
        self.dirs = [Path(d) for d in dirs]
        self.index_dirs = [Path(d) for d in index_dirs]  # roots holding og/<version>/bucket=XX.parquet
        self._versions, self._present = None, ()

    def versions(self) -> dict[str, Version]:
        present = tuple(d for d in self.dirs if (d / "orthogroups.parquet").exists())
        if self._versions is None or present != self._present or any(Version.signature(v.dir) != v.sig for v in self._versions.values()):
            out = {}
            for d in self.dirs:
                if not (d / "orthogroups.parquet").exists():
                    continue
                m = json.loads((d / "manifest.json").read_text()) if (d / "manifest.json").exists() else {}
                vid = "TOY" if m.get("toy") else (d.parent.name if d.name == "web" else d.name)
                while vid in out:
                    vid += "'"
                out[vid] = Version(vid, d, m)
            if self._versions:  # keep caches of unchanged versions
                for k, v in out.items():
                    if k in self._versions and self._versions[k].sig == v.sig:
                        out[k] = self._versions[k]
            self._versions, self._present = out, present
        return self._versions

    def default(self) -> Version | None:
        v = self.versions()
        return next(iter(v.values()), None)

    def pick(self, vid: str | None) -> Version | None:
        v = self.versions()
        return v.get(vid) if vid in v else self.default()

    def by_accession(self, acc: str) -> list[dict]:
        """Orthogroups/SOGs containing a member accession, from the sharded index written by
        scripts/webapp_build_og_index.py: [{version, og_id, label, sogs: [{sog_id, state, event}]}]."""
        acc = acc.strip().upper()
        out = []
        for vid, v in self.versions().items():
            for root in self.index_dirs:
                f = root / "og" / vid / f"bucket={bucket(acc)}.parquet"
                if not f.exists():
                    continue
                for r in pq.read_table(f, filters=[("accession", "=", acc)]).to_pylist():
                    og = v.og(r["og_id"]) or {}
                    ev = {s["sog_id"]: s.get("event") for s in v.og_sogs(r["og_id"])} if len(r["sog_k"]) < 200 else {}
                    out.append({"version": vid, "og_id": r["og_id"], "label": og.get("label"),
                                "n_sogs": og.get("n_sogs"),
                                "sogs": [{"sog_id": f'{r["og_id"]}:g{k}', "state": st,
                                          "event": ev.get(f'{r["og_id"]}:g{k}')}
                                         for k, st in sorted(zip(r["sog_k"], r["states"]))]})
                break
        return out

    def sources(self):
        return [(v.dir, v.manifest) for v in self.versions().values()]

    def find_version(self, og_id, preferred: Version | None):
        """Version holding og_id (preferred first)."""
        vs = list(self.versions().values())
        if preferred:
            vs = [preferred] + [v for v in vs if v is not preferred]
        for v in vs:
            if v.og(og_id):
                return v
        return None

    def get(self, og_id: str, ver: Version | None = None, sog_ids=None, sog_page=None) -> dict | None:
        """OG record with its SOGs and member rows. `sog_page=(page, size)` limits the SOG columns (and so the
        state rows read); `sog_ids` selects explicit SOGs. `taxa` = all taxids of the OG."""
        v = ver if ver is not None and ver.og(og_id) else self.find_version(og_id, ver)
        if v is None:
            return None
        sogs = v.og_sogs(og_id)
        n_total = len(sogs)
        if sog_ids is not None:
            sogs = [s for s in sogs if s["sog_id"] in set(sog_ids)]
        elif sog_page is not None:
            pg, size = sog_page
            sogs = sogs[(pg - 1) * size: pg * size]
        return {"manifest": v.manifest, "version": v.id, "orthogroup": v.og(og_id), "sogs": sogs, "n_sogs_total": n_total,
                "states": v.og_states(og_id, [s["sog_id"] for s in sogs]), "taxa": v.og_taxa(og_id),
                "species": v.species(), "tree": v.tree()}


def matrix(data: dict):
    """(pruned tree, leaf order [taxid str], {sog_id: {taxid: [rows]}})."""
    cells = {}
    taxa = {str(t) for t in data["taxa"]} if data.get("taxa") is not None else set()
    for r in data["states"]:
        cells.setdefault(r["sog_id"], {}).setdefault(str(r["taxid"]), []).append(r)
        taxa.add(str(r["taxid"]))
    t = prune(data["tree"], taxa) or data["tree"]
    return t, leaves(t), cells


def strongest(rows):
    return min((r["state"] for r in rows), key=STATE_ORDER.index)
