"""Per-protein binding *sites*: atlas pockets of one accession merged by residue overlap (DESIGN.md §3.1).

A site is a connected component (union-find) of the graph whose nodes are the atlas pockets of one UniProt
accession and whose edges join two pockets with lining-residue Jaccard >= ``jaccard`` (default 0.3), computed on
UniProt positions. Pockets from different AFDB fragments of the same accession are compared after converting the
fragment-local numbering to UniProt positions, so duplicated sites in overlapping fragment windows merge.

Fragment numbering (verified, see ``scripts/verify_fragment_offsets.py`` and tests/test_sites.py): AFDB long
proteins are split into 1400-residue windows with a 200-residue step; fragment k covers UniProt positions
200*(k-1)+1 .. 200*(k-1)+1400 and its author residue numbers restart at 1, i.e.
``uniprot_pos = author_resnum + 200*(k-1)`` (mmCIF ``_struct_ref_seq.db_align_beg`` = 200*(k-1)+1 with
``pdbx_auth_seq_align_beg`` = 1).

Everything here works on one pyarrow Table at a time (one atlas unit file, ~1e4 pockets), so the pipeline
streams over the atlas without ever holding more than one file's fingerprints.
"""
from __future__ import annotations

import gzip
import hashlib
import re
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components

SCHEMA_VERSION = "sites-v1"
DEFAULT_JACCARD = 0.3
FRAGMENT_STEP = 200      # AFDB: fragment k starts at UniProt position FRAGMENT_STEP*(k-1)+1
FRAGMENT_WINDOW = 1400   # AFDB: window length of each fragment (the last one may be shorter)
_POS_BITS = 24           # key = accession_code << 24 | uniprot_pos (positions < 16.7 M)

TIER_SOURCE_DEFAULT = {1: "afdb_v6_tar", 2: "afdb_swissprot_v6_tar"}  # = afdb_common.TIER_SOURCE_DEFAULT

#: pocket-table columns read by build_sites (fingerprint optional)
POCKET_COLUMNS = ["unit_id", "tier", "accession", "fragment", "model_version", "taxid", "organism", "pocket_id",
                  "burrow_score", "n_points", "volume_A3", "site_volume_A3", "residues", "plddt_mean",
                  "plddt_min", "plddt_frac_lt70"]


POCKET_SCHEMA = pa.schema([
    ("unit_id", pa.int32()), ("tier", pa.int8()), ("accession", pa.string()), ("fragment", pa.int16()),
    ("model_version", pa.int8()), ("taxid", pa.int32()), ("organism", pa.string()), ("pocket_id", pa.string()),
    ("burrow_score", pa.float64()), ("n_points", pa.int64()), ("volume_A3", pa.float64()),
    ("site_volume_A3", pa.float64()), ("residues", pa.string()), ("plddt_mean", pa.float64()),
    ("plddt_min", pa.float64()), ("plddt_frac_lt70", pa.float64()), ("source", pa.string()),
    ("fingerprint", pa.list_(pa.float16(), 384)),
])


def code_sha256() -> str:
    """sha256 of this module's source; part of the resume key so code changes invalidate outputs."""
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


# --------------------------------------------------------------------------------------------- numbering
def fragment_offset(fragment) -> np.ndarray | int:
    """UniProt offset to add to fragment-local author residue numbers (0 for fragment 1)."""
    return FRAGMENT_STEP * (np.asarray(fragment, dtype=np.int64) - 1) if not np.isscalar(fragment) \
        else FRAGMENT_STEP * (int(fragment) - 1)


_SRS_RE = re.compile(r"^_struct_ref_seq\.(\S+)\s+(\S+)", re.M)


def struct_ref_seq(cif: bytes | str | Path) -> dict:
    """Parse the (single-row, key-value form) ``_struct_ref_seq`` category of an AFDB mmCIF.

    Accepts a path (.cif or .cif.gz) or the raw / gzipped bytes. Returns a dict of the category items; the
    fragment offset is ``int(d['db_align_beg']) - int(d['pdbx_auth_seq_align_beg'])``."""
    if isinstance(cif, Path) or (isinstance(cif, str) and "\n" not in cif and len(cif) < 4096):
        cif = Path(cif).read_bytes()
    if isinstance(cif, bytes):
        if cif[:2] == b"\x1f\x8b":
            cif = gzip.decompress(cif)
        cif = cif.decode()
    d = dict(_SRS_RE.findall(cif))
    if not d:
        raise ValueError("no key-value _struct_ref_seq category found (loop form not supported)")
    return d


def parse_residues(residues: pa.Array | pa.ChunkedArray | list) -> tuple[np.ndarray, np.ndarray]:
    """Parse atlas ``residues`` strings ("A:53 A:79 ...") into CSR form.

    Returns (row_index int64[n_tokens], resnum int64[n_tokens]) with author residue numbers. Empty / null
    strings give no tokens. Raises on tokens that are not ``<chain>:<int>``."""
    arr = pa.chunked_array([residues]) if isinstance(residues, list) else residues
    if isinstance(arr, pa.ChunkedArray):
        arr = arr.combine_chunks()
    arr = pc.fill_null(arr, "")
    toks = pc.split_pattern(pc.utf8_trim_whitespace(arr), " ")
    lens = pc.list_value_length(toks).to_numpy(zero_copy_only=False).astype(np.int64)
    flat = pc.list_flatten(toks)
    nonempty = pc.not_equal(pc.utf8_length(flat), 0)
    row = np.repeat(np.arange(len(arr), dtype=np.int64), lens)
    row = row[nonempty.to_numpy(zero_copy_only=False)]
    flat = pc.filter(flat, nonempty)
    parts = pc.split_pattern(flat, ":")
    if len(flat) and pc.min(pc.list_value_length(parts)).as_py() != 2:
        bad = pc.filter(flat, pc.not_equal(pc.list_value_length(parts), 2))[0]
        raise ValueError(f"unparsable residue token {bad!r}")
    num = pc.list_element(parts, 1)
    try:
        resnum = pc.cast(num, pa.int64()).to_numpy(zero_copy_only=False)
    except pa.ArrowInvalid as e:  # insertion codes or garbage: fail loudly, never guess
        raise ValueError(f"non-integer residue number in atlas residues: {e}") from e
    return row, resnum


# --------------------------------------------------------------------------------------------- merging
def union_find(n: int, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Component labels (0..n_comp-1, ordered by smallest member index) of the graph with edges (i, j).

    Same result as a classic union-find over the edge list; implemented via scipy's connected_components and
    relabelled so that labels are a deterministic function of the edge set."""
    g = coo_matrix((np.ones(len(i), dtype=np.int8), (i, j)), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    # relabel by first occurrence -> independent of scipy's internal order
    _, first = np.unique(lab, return_index=True)
    order = np.argsort(first, kind="stable")
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    return remap[lab]


def jaccard_edges(row: np.ndarray, key: np.ndarray, n: int, jaccard: float):
    """Pairs (i<j) of rows whose key sets have Jaccard >= jaccard. ``(row, key)`` must be unique pairs."""
    if len(row) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0)
    ukey, kidx = np.unique(key, return_inverse=True)
    A = csr_matrix((np.ones(len(row), dtype=np.int32), (row, kidx)), shape=(n, len(ukey)))
    size = np.asarray(A.sum(axis=1)).ravel()
    inter = (A @ A.T).tocoo()
    m = inter.row < inter.col
    i, j, x = inter.row[m].astype(np.int64), inter.col[m].astype(np.int64), inter.data[m].astype(np.int64)
    union = size[i] + size[j] - x
    jac = x / union
    keep = x >= jaccard * union - 1e-9   # tolerant of float rounding at exactly-threshold fractions
    return i[keep], j[keep], jac[keep]


def _fill_source(t: pa.Table) -> pa.Table:
    tiers = t.column("tier").to_numpy(zero_copy_only=False)
    derived = pa.array([TIER_SOURCE_DEFAULT.get(int(x)) for x in tiers], pa.string())
    if "source" not in t.column_names:
        return t.append_column("source", derived)
    src = t.column("source").combine_chunks().cast(pa.string())
    return t.set_column(t.column_names.index("source"), "source", pc.coalesce(src, derived))


def build_sites(pockets: pa.Table, jaccard: float = DEFAULT_JACCARD, with_fingerprints: bool = True):
    """Merge the pockets of each accession in ``pockets`` into sites.

    ``pockets``: atlas pocket rows (columns POCKET_COLUMNS, optional ``source`` and ``fingerprint``). All
    pockets of an accession must be in this table (callers handle accessions split over files).
    Returns (sites, fingerprints, pocket_site) where ``sites`` has one row per site, ``fingerprints`` has
    (site_id, rep_pocket_id, fingerprint) for the representative pocket (or None), and ``pocket_site`` maps
    every input pocket_id to its site_id (for invariant checks).
    """
    t = _fill_source(pockets)
    n = t.num_rows
    if n == 0:
        return _empty_sites(), (_empty_fp() if with_fingerprints else None), \
            pa.table({"pocket_id": pa.array([], pa.string()), "site_id": pa.array([], pa.string())})
    # deterministic row order: accession, pocket_id
    order = pc.sort_indices(t, sort_keys=[("accession", "ascending"), ("pocket_id", "ascending")])
    t = t.take(order)
    pid = t.column("pocket_id").combine_chunks().cast(pa.string())
    if pc.count_distinct(pid).as_py() != n:
        raise ValueError("duplicate pocket_id in input")
    acc = t.column("accession").combine_chunks().cast(pa.string())
    acc_enc = pc.dictionary_encode(acc)
    acc_code = acc_enc.indices.to_numpy(zero_copy_only=False).astype(np.int64)
    frag = t.column("fragment").to_numpy(zero_copy_only=False).astype(np.int64)

    # residues -> UniProt positions -> per-accession keys
    row, resnum = parse_residues(t.column("residues"))
    pos = resnum + fragment_offset(frag)[row]
    if len(pos) and (pos.min() < 1 or pos.max() >= (1 << _POS_BITS)):
        raise ValueError(f"UniProt position out of range: {pos.min()}..{pos.max()}")
    key = (acc_code[row] << _POS_BITS) | pos
    rk = np.unique(np.stack([row, key], 1), axis=0) if len(row) else np.zeros((0, 2), np.int64)
    row, key = rk[:, 0], rk[:, 1]

    # union-find over Jaccard edges (edges only join rows that share a key -> same accession by construction)
    ei, ej, _ = jaccard_edges(row, key, n, jaccard)
    comp = union_find(n, ei, ej)
    ncomp = int(comp.max()) + 1

    score = t.column("burrow_score").to_numpy(zero_copy_only=False).astype(np.float64)
    # sort rows: component, score desc, pocket_id asc (rows are already pocket_id-sorted within an accession,
    # and components never span accessions, so row index breaks ties by pocket_id)
    ridx = np.arange(n)
    srt = np.lexsort((ridx, -score, comp))
    c_s = comp[srt]
    starts = np.flatnonzero(np.r_[True, c_s[1:] != c_s[:-1]])
    rep = srt[starts]                                   # representative row per component (highest score)
    nmem = np.diff(np.r_[starts, n])

    def red(col, ufunc):
        v = t.column(col).to_numpy(zero_copy_only=False).astype(np.float64)[srt]
        return ufunc.reduceat(v, starts)

    max_score = score[rep]
    mean_score = np.add.reduceat(score[srt], starts) / nmem

    # site numbering per accession: by max score desc, then rep pocket_id
    comp_acc = acc_code[rep]
    sorder = np.lexsort((rep, -max_score, comp_acc))   # rep row order == pocket_id order within accession
    k = np.empty(ncomp, np.int64)
    ca = comp_acc[sorder]
    gstart = np.r_[True, ca[1:] != ca[:-1]]
    grp_first = np.maximum.accumulate(np.where(gstart, np.arange(ncomp), 0))
    k[sorder] = np.arange(ncomp) - grp_first
    acc_rep = acc.take(pa.array(rep)).to_numpy(zero_copy_only=False)
    site_id = np.array([f"{a}:s{x}" for a, x in zip(acc_rep, k)], dtype=object)

    # residues per site: union of member keys
    ck = np.unique(np.stack([comp[row], key & ((1 << _POS_BITS) - 1)], 1), axis=0) if len(row) \
        else np.zeros((0, 2), np.int64)
    res_counts = np.bincount(ck[:, 0], minlength=ncomp)
    res_off = np.r_[0, np.cumsum(res_counts)].astype(np.int32)
    residues = pa.ListArray.from_arrays(pa.array(res_off), pa.array(ck[:, 1].astype(np.int32)))
    if len(ck):
        res_first = np.where(res_counts > 0, ck[np.minimum(res_off[:-1], len(ck) - 1), 1], 0)
        res_last = np.where(res_counts > 0, ck[np.maximum(res_off[1:] - 1, 0), 1], 0)
    else:
        res_first = res_last = np.zeros(ncomp, np.int64)
    # fragments per site
    cf = np.unique(np.stack([comp, frag], 1), axis=0)
    f_off = np.r_[0, np.cumsum(np.bincount(cf[:, 0], minlength=ncomp))].astype(np.int32)
    fragments = pa.ListArray.from_arrays(pa.array(f_off), pa.array(cf[:, 1].astype(np.int16)))
    # members (score desc, pocket_id asc)
    m_off = np.r_[0, np.cumsum(nmem)].astype(np.int32)
    members = pa.ListArray.from_arrays(pa.array(m_off), pid.take(pa.array(srt)))

    def take_rep(col):
        return t.column(col).combine_chunks().take(pa.array(rep))

    cols = {
        "site_id": pa.array(site_id, pa.string()),
        "accession": acc.take(pa.array(rep)),
        "site_index": pa.array(k.astype(np.int32)),
        "taxid": take_rep("taxid"),
        "organism": take_rep("organism").cast(pa.string()),
        "tier": take_rep("tier"),
        "unit_id": take_rep("unit_id"),
        "model_version": take_rep("model_version"),
        "source": take_rep("source"),
        "n_pockets": pa.array(nmem.astype(np.int32)),
        "member_pocket_ids": members,
        "fragments": fragments,
        "n_fragments": pa.array(np.diff(f_off).astype(np.int16)),
        "residues": residues,
        "n_res": pa.array(res_counts.astype(np.int32)),
        "res_min": pa.array(res_first.astype(np.int32)),
        "res_max": pa.array(res_last.astype(np.int32)),
        "max_burrow_score": pa.array(max_score),
        "mean_burrow_score": pa.array(mean_score),
        "plddt_mean_rep": pa.array(t.column("plddt_mean").to_numpy(zero_copy_only=False)[rep].astype(np.float64)),
        "plddt_mean_min": pa.array(red("plddt_mean", np.minimum)),
        "plddt_mean_max": pa.array(red("plddt_mean", np.maximum)),
        "plddt_min_min": pa.array(red("plddt_min", np.minimum)),
        "plddt_frac_lt70_max": pa.array(red("plddt_frac_lt70", np.maximum)),
        "volume_sum": pa.array(red("volume_A3", np.add)),
        "volume_max": pa.array(red("volume_A3", np.maximum)),
        "n_points_sum": pa.array(red("n_points", np.add).astype(np.int64)),
        "rep_pocket_id": pid.take(pa.array(rep)),
        "rep_fragment": pa.array(frag[rep].astype(np.int16)),
        "rep_site_volume_A3": take_rep("site_volume_A3"),
    }
    sites = pa.table(cols).cast(SITES_SCHEMA)
    # final row order: accession, site_index
    so = pc.sort_indices(sites, sort_keys=[("accession", "ascending"), ("site_index", "ascending")])
    sites = sites.take(so)
    so_np = so.to_numpy()
    fps = None
    if with_fingerprints:
        if "fingerprint" not in t.column_names:
            raise ValueError("with_fingerprints=True but no fingerprint column")
        fp = t.column("fingerprint").combine_chunks().take(pa.array(rep[so_np]))
        fps = pa.table({"site_id": sites.column("site_id"), "rep_pocket_id": sites.column("rep_pocket_id"),
                        "fingerprint": fp}).cast(FP_SCHEMA)
    pocket_site = pa.table({"pocket_id": pid, "site_id": pa.array(site_id[comp], pa.string())})
    return sites, fps, pocket_site


SITES_SCHEMA = pa.schema([
    ("site_id", pa.string()), ("accession", pa.string()), ("site_index", pa.int32()), ("taxid", pa.int32()),
    ("organism", pa.string()), ("tier", pa.int8()), ("unit_id", pa.int32()), ("model_version", pa.int8()),
    ("source", pa.string()), ("n_pockets", pa.int32()), ("member_pocket_ids", pa.list_(pa.string())),
    ("fragments", pa.list_(pa.int16())), ("n_fragments", pa.int16()), ("residues", pa.list_(pa.int32())),
    ("n_res", pa.int32()), ("res_min", pa.int32()), ("res_max", pa.int32()), ("max_burrow_score", pa.float64()),
    ("mean_burrow_score", pa.float64()), ("plddt_mean_rep", pa.float64()), ("plddt_mean_min", pa.float64()),
    ("plddt_mean_max", pa.float64()), ("plddt_min_min", pa.float64()), ("plddt_frac_lt70_max", pa.float64()),
    ("volume_sum", pa.float64()), ("volume_max", pa.float64()), ("n_points_sum", pa.int64()),
    ("rep_pocket_id", pa.string()), ("rep_fragment", pa.int16()), ("rep_site_volume_A3", pa.float64()),
])
FP_SCHEMA = pa.schema([("site_id", pa.string()), ("rep_pocket_id", pa.string()),
                       ("fingerprint", pa.list_(pa.float16(), 384))])


def _empty_sites() -> pa.Table:
    return SITES_SCHEMA.empty_table()


def _empty_fp() -> pa.Table:
    return FP_SCHEMA.empty_table()
