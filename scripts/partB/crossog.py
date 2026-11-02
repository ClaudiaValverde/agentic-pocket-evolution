#!/usr/bin/env python
"""Cross-orthogroup site matching (exploratory, not preregistered; DESIGN §12 D31). Method: docs/CROSSOG_METHOD.md
(frozen, committed 77db16c before any output). Results: CROSSOG.md, results/partB_crossog/<snapshot>/.

Modes
  validate   the 7 `known` items of scripts/cases/validation_controls.yaml on O1_r3 and O2_r1
  genome     all O1 same-family OG pairs (+ --null: seeded sample of same-family pairs, shuffled correspondence)

    export UV_CACHE_DIR=/tmp/pe_big/uvcache
    .venv/bin/python scripts/partB/crossog.py validate --out results/partB_crossog/20261030T085136Z
    .venv/bin/python scripts/partB/crossog.py genome --version O1 --out results/partB_crossog/20261030T085136Z
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
import time
import zlib
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")   # CPU only

import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(P / "src"))
sys.path.insert(0, str(P / "scripts"))

SNAP = "20261030T085136Z"
TAU = 0.2
ALNTM_MIN = 0.5
COVER_MIN = 0.8
MIN_SPECIES = 3
MAX_MEMBERS = 10
NULL_PAIRS = 200_000
NULL_BIN_MIN = 500
BINS = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0000001]
BH_Q = 0.05
FOLDSEEK = str(P / "envs/ancestra/bin/foldseek")
FSDB = "/tmp/pe_big/partA/foldseek/afdb_v6"
SITE_FP = f"/tmp/pe_big/sites/{SNAP}/site_fp"
RUNS = {"O1": "O1_r3", "O2": "O2_r1"}
GAP = ord("-")


# ============================================================================================ pure functions (tested)
def map_from_aln(qaln: str, taln: str, qstart: int, tstart: int) -> tuple[np.ndarray, np.ndarray]:
    """Aligned residue pairs (1-based q, t) from Foldseek qaln/taln strings."""
    qa = np.frombuffer(qaln.encode(), np.uint8)
    ta = np.frombuffer(taln.encode(), np.uint8)
    qn, tn = qa != GAP, ta != GAP
    qpos = qstart - 1 + np.cumsum(qn)
    tpos = tstart - 1 + np.cumsum(tn)
    m = qn & tn
    return qpos[m].astype(np.int64), tpos[m].astype(np.int64)


def map_from_cigar(cigar: str, qstart: int, tstart: int) -> tuple[np.ndarray, np.ndarray]:
    """Foldseek/MMseqs2 CIGAR (M; I = query residue not in target; D = target residue not in query)."""
    qs, ts, q, t, num = [], [], qstart, tstart, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num) if num else 1
        num = ""
        if ch == "M":
            qs.append(np.arange(q, q + n)); ts.append(np.arange(t, t + n)); q += n; t += n
        elif ch == "I":
            q += n
        elif ch == "D":
            t += n
    if not qs:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(qs).astype(np.int64), np.concatenate(ts).astype(np.int64)


def shift_targets(tpos: np.ndarray, L: int, offset: int) -> np.ndarray:
    """Null correspondence: t -> ((t - 1 + offset) mod L) + 1."""
    return ((tpos - 1 + offset) % L) + 1


def null_offset(rng: np.random.Generator, L: int) -> int:
    return int(rng.integers(L // 4, max(L // 4 + 1, (3 * L) // 4 + 1)))


def row_positions(aligned_row: np.ndarray) -> np.ndarray:
    """Per MSA column: 1-based residue position of the row, 0 at gaps."""
    ng = aligned_row != GAP
    return np.where(ng, np.cumsum(ng), 0).astype(np.int64)


def project_footprint(F_cols: np.ndarray, src_colpos: np.ndarray, amap: np.ndarray, dst_pos2col: np.ndarray):
    """Map footprint columns of the source OG into the target OG's columns.
    src_colpos: source rep residue per source column (0 = gap); amap: source residue -> target residue (0 = unaligned),
    index 0 unused; dst_pos2col: target rep residue (1-based) -> target column (0-based) as array index p-1.
    Returns (mapped target columns (unique, sorted), n_unmapped)."""
    res = src_colpos[F_cols]
    ok = res > 0
    res = res[ok]
    inb = res < len(amap)
    t = np.zeros(len(res), np.int64)
    t[inb] = amap[res[inb]]
    good = (t > 0) & (t <= len(dst_pos2col))
    cols = np.unique(dst_pos2col[t[good] - 1])
    n_unmapped = int((~ok).sum() + (~good).sum())
    return cols, n_unmapped


def directional_jaccard(mapped: np.ndarray, n_unmapped: int, F_dst: np.ndarray) -> float:
    inter = len(np.intersect1d(mapped, F_dst, assume_unique=True))
    union = len(np.union1d(mapped, F_dst)) + n_unmapped
    return inter / union if union else 0.0


def bin_index(identity: np.ndarray, null_identity: np.ndarray, bins=BINS, min_n=NULL_BIN_MIN) -> tuple[np.ndarray, np.ndarray]:
    """Bin assignment with sparse-bin merging (into the next lower bin; the lowest merges upward).
    Returns (bin of each query identity, bin of each null identity) as merged-bin labels."""
    nb = len(bins) - 1
    b_null = np.clip(np.digitize(null_identity, bins) - 1, 0, nb - 1)
    counts = np.bincount(b_null, minlength=nb)
    label = np.arange(nb)
    # merge downward from the top: a sparse bin takes the label of the next lower bin
    for k in range(nb - 1, 0, -1):
        if counts[k] < min_n:
            label[label == k] = label[k - 1]
            counts[k - 1] += counts[k]
            counts[k] = 0
    if counts[0] < min_n and nb > 1:   # lowest bin merges upward into the next existing label
        up = [l for l in np.unique(label) if l != 0]
        if up:
            label[label == 0] = min(up)
    b_q = np.clip(np.digitize(identity, bins) - 1, 0, nb - 1)
    return label[b_q], label[b_null]


def empirical_p(D: np.ndarray, qbin: np.ndarray, null_d: np.ndarray, null_bin: np.ndarray) -> np.ndarray:
    out = np.full(len(D), np.nan)
    for b in np.unique(qbin):
        nd = np.sort(null_d[null_bin == b])
        m = qbin == b
        ge = len(nd) - np.searchsorted(nd, D[m], side="left")
        out[m] = (1 + ge) / (1 + len(nd))
    return out


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def cos_dist(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    A = A / np.linalg.norm(A, axis=1, keepdims=True)
    B = B / np.linalg.norm(B, axis=1, keepdims=True)
    return 1.0 - A @ B.T


def call_b(n_P_species_x: int, states_y: pd.DataFrame) -> dict:
    """Presence/absence call. states_y: rows (taxid, state) of the absent-side OG at the projected footprint."""
    nPW = int(states_y.state.isin(["P", "W"]).sum())
    nA_sp = int(states_y.loc[states_y.state == "A", "taxid"].nunique())
    return {"n_P_species_x": n_P_species_x, "n_PW_y": nPW, "n_A_species_y": nA_sp,
            "call_b": bool(n_P_species_x >= MIN_SPECIES and nPW == 0 and nA_sp >= MIN_SPECIES)}


def family_labels(members: pd.DataFrame, ted: pd.DataFrame, pfam: pd.DataFrame, frac: float = 0.5) -> dict:
    """{og_id: (ted_labels frozenset, pfam_labels frozenset)}: labels held by >= frac of the OG's members.
    ted: accession, cath_label (H level); pfam: accession, pfam, clan."""
    n = members.groupby("og_id").accession.nunique()
    out = {}
    t = members.merge(ted.drop_duplicates(["accession", "cath_label"]), on="accession")
    tc = t.groupby(["og_id", "cath_label"]).accession.nunique().reset_index()
    tc = tc[tc.accession >= frac * tc.og_id.map(n)]
    tl = tc.groupby("og_id").cath_label.agg(frozenset).to_dict()
    pl_ = pd.concat([pfam[["accession", "pfam"]].rename(columns={"pfam": "label"}),
                     pfam[["accession", "clan"]].dropna().rename(columns={"clan": "label"})]).drop_duplicates()
    p = members.merge(pl_, on="accession")
    pc = p.groupby(["og_id", "label"]).accession.nunique().reset_index()
    pc = pc[pc.accession >= frac * pc.og_id.map(n)]
    plab = pc.groupby("og_id").label.agg(frozenset).to_dict()
    for og in n.index:
        out[og] = (tl.get(og, frozenset()), plab.get(og, frozenset()))
    return out


def has_any_ted(members: pd.DataFrame, ted: pd.DataFrame) -> set:
    return set(members[members.accession.isin(set(ted.accession))].og_id)


def same_family(la, lb, a_has_ted: bool, b_has_ted: bool) -> str | None:
    """Shared label or None. TED with TED when both OGs have TED labels; else Pfam/clan."""
    if a_has_ted and b_has_ted:
        s = la[0] & lb[0]
        return "TED:" + ",".join(sorted(s)) if s else None
    s = la[1] & lb[1]
    return "Pfam:" + ",".join(sorted(s)) if s else None


# ============================================================================================ data access
class Run:
    def __init__(self, version: str):
        self.version = version
        self.run_id = RUNS[version]
        self.bulk = Path(f"/tmp/pe_big/sogs/{SNAP}/{self.run_id}")
        self.web = P / f"results/sogs/{SNAP}/{self.run_id}/web"
        t0 = time.time()
        self.ogs = pd.read_parquet(self.bulk / "inputs/ogs.parquet", columns=["og_id", "shard", "level"])
        self.shard_of = dict(zip(self.ogs.og_id, self.ogs.shard))
        self.members = pd.read_parquet(self.bulk / "inputs/members.parquet")
        sg = pd.read_parquet(self.web / "sogs.parquet", columns=["og_id", "sog_id", "footprint_columns", "small_footprint",
                                                                   "empty_footprint", "n_P", "n_A"])
        self.sogs = sg
        import pyarrow.dataset as ds
        files = sorted(glob.glob(str(self.bulk / "states_t0.2" / "shard=*.parquet")))
        stP = ds.dataset(files, format="parquet").to_table(
            columns=["og_id", "sog_id", "accession", "taxid", "site_id", "footprint_plddt", "model_accession", "state"],
            filter=ds.field("state") == "P").to_pandas()
        self.P = stP
        nsp = stP.groupby("sog_id").taxid.nunique()
        bad = set(sg.sog_id[sg.small_footprint | sg.empty_footprint])
        el = nsp[(nsp >= MIN_SPECIES)]
        self.eligible = sorted(set(el.index) - bad)
        self.nPsp = nsp.to_dict()
        self.fp_cols = dict(zip(sg.sog_id, sg.footprint_columns))
        self._msa = {}
        self._models = {}
        self._sites = {}
        print(f"[{self.version}] loaded run {self.run_id}: {len(self.eligible)} eligible SOGs ({time.time() - t0:.0f}s)",
              flush=True)

    # -- per shard caches
    def load_shards(self, og_ids, what=("msa",)):
        shards = sorted({self.shard_of[o] for o in og_ids if o in self.shard_of})
        need = set(og_ids)
        for s in shards:
            if "msa" in what:
                m = pd.read_parquet(self.bulk / "msa" / f"shard={s:04d}.parquet")
                m = m[m.og_id.isin(need)]
                for og, g in m.groupby("og_id"):
                    ids = g.accession.tolist()
                    mat = np.frombuffer("".join(g.aligned).encode(), np.uint8).reshape(len(ids), -1)
                    self._msa[og] = (ids, mat)
            if "models" in what:
                d = self.bulk / "inputs/shard_models" / f"shard={s}"
                if d.exists():
                    m = pd.read_parquet(d)
                    m = m[m.og_id.isin(need)]
                    for og, g in m.groupby("og_id"):
                        self._models[og] = g
            if "sites" in what:
                d = self.bulk / "inputs/shard_sites" / f"shard={s}"
                if d.exists():
                    m = pd.read_parquet(d)
                    m = m[m.og_id.isin(need)]
                    for og, g in m.groupby("og_id"):
                        self._sites[og] = g

    def msa(self, og):
        return self._msa.get(og)


def load_fs_lookup():
    lk = pd.read_csv(FSDB + ".lookup", sep="\t", header=None, names=["key", "name", "file"], usecols=[0, 1])
    acc = lk.name.str.replace(r"^AF-", "", regex=True).str.replace(r"-F1-model_v\d+.*$", "", regex=True)
    keep = lk.name.str.contains("-F1-model_v")
    idx = pd.read_csv(FSDB + ".index", sep="\t", header=None, names=["key", "off", "len"])
    L = dict(zip(idx.key, idx.len - 2))   # entry = sequence + '\n' + '\0'
    return {a: (n, k, L.get(k)) for a, n, k in zip(acc[keep], lk.name[keep], lk.key[keep])}


def run_foldseek(pairs: list[tuple[str, str]], fsl: dict, workdir: Path, threads: int = 48) -> pd.DataFrame:
    """pairs of model accessions (q, t) -> alignment table (query/target = model accessions)."""
    from pocketevo import homology as hom
    workdir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(pairs, columns=["q", "t"]).drop_duplicates()
    df["qname"] = df.q.map(lambda a: fsl[a][0])
    df["tname"] = df.t.map(lambda a: fsl[a][0])
    name2key = {v[0]: v[1] for v in fsl.values()}
    fields = ["query", "target", "evalue", "bits", "fident", "alnlen", "qstart", "tstart", "qlen", "tlen", "cigar",
              "alntmscore", "qtmscore", "ttmscore"]
    out = []
    CH = 2_000_000
    for k in range(0, len(df), CH):
        sub = df.iloc[k:k + CH]
        r = hom.run_pair_search(sub[["qname", "tname"]], FSDB, name2key, workdir / f"c{k // CH}", "foldseek", FOLDSEEK,
                                threads=threads, evalue=10.0, fields=fields)
        r["q"] = sub.q.values
        r["t"] = sub.t.values
        out.append(r)
    return pd.concat(out, ignore_index=True)


# ============================================================================================ pipeline
class CrossOG:
    def __init__(self, run: Run, fsl: dict, work: Path, shuffle_seed: int | None = None):
        self.r, self.fsl, self.work = run, fsl, work
        self.shuffle_seed = shuffle_seed
        self.aln = {}

    # ---- representatives
    def sog_reps(self, sog_ids):
        Pr = self.r.P[self.r.P.sog_id.isin(set(sog_ids))]
        Pr = Pr[Pr.model_accession.isin(self.fsl.keys())]
        Pr = Pr.sort_values(["sog_id", "footprint_plddt", "accession"], ascending=[True, False, True])
        rep = Pr.drop_duplicates("sog_id").set_index("sog_id")[["accession", "model_accession"]]
        return {s: (a, m) for s, a, m in zip(rep.index, rep.accession, rep.model_accession)}

    def og_reps(self, og_ids):
        out = {}
        mem = self.r.members
        for og in og_ids:
            g = self.r._models.get(og)
            ms = self.r.msa(og)
            if g is None or ms is None:
                continue
            mm = mem[mem.og_id == og].set_index("accession")
            best = None
            ncol = ms[1].shape[1]
            for acc, macc, pq, nres in zip(g.accession, g.model_accession, g.plddt_q, g.n_res):
                if macc not in self.fsl or acc not in ms[0]:
                    continue
                st = mm.status.get(acc) if acc in mm.index else None
                if st != "ok":
                    continue
                score = float(np.sum(np.asarray(pq, np.float64)) / 2.0) / ncol
                key = (-score, acc)
                if best is None or key < best[0]:
                    best = (key, acc, macc)
            if best:
                out[og] = (best[1], best[2])
        return out

    # ---- alignment maps
    def amap(self, qm: str, tm: str):
        """(map q->t as array indexed by q residue, alntm, fident) for model accessions; reads either direction."""
        rec = self.aln.get((qm, tm))
        rev = False
        if rec is None:
            rec = self.aln.get((tm, qm))
            rev = True
        if rec is None or rec[0] is None:
            return None
        cigar, qs, ts, ql, tl, alntm, fid = rec
        qpos, tpos = map_from_cigar(cigar, int(qs), int(ts))
        if rev:
            qpos, tpos, ql, tl = tpos, qpos, tl, ql
        if self.shuffle_seed is not None:
            rng = np.random.default_rng([self.shuffle_seed, zlib.crc32(f"{qm}|{tm}".encode())])
            tpos = shift_targets(tpos, int(tl), null_offset(rng, int(tl)))
        m = np.zeros(int(ql) + 2, np.int64)
        m[qpos] = tpos
        return m, float(alntm), float(fid)

    def colpos(self, og, acc):
        ids, mat = self.r.msa(og)
        return row_positions(mat[ids.index(acc)])

    def pos2col(self, og, acc):
        ids, mat = self.r.msa(og)
        return np.flatnonzero(mat[ids.index(acc)] != GAP).astype(np.int64)

    def prepare(self, og_pairs: pd.DataFrame):
        """og_pairs: og_x, og_y. Loads data, picks reps, runs Foldseek on every needed rep pair."""
        ogs = sorted(set(og_pairs.og_x) | set(og_pairs.og_y))
        self.r.load_shards(ogs, what=("msa", "models", "sites"))
        el = [s for s in self.r.eligible if s.split(":")[0] in set(ogs)]
        self.el_by_og = defaultdict(list)
        for s in el:
            self.el_by_og[s.split(":")[0]].append(s)
        self.srep = self.sog_reps(el)
        self.orep = self.og_reps(ogs)
        need = set()
        for x, y in zip(og_pairs.og_x, og_pairs.og_y):
            rx = [self.srep[s][1] for s in self.el_by_og[x] if s in self.srep]
            ry = [self.srep[s][1] for s in self.el_by_og[y] if s in self.srep]
            for a in rx:
                for b in ry:
                    need.add(tuple(sorted((a, b))))
                if y in self.orep:
                    need.add(tuple(sorted((a, self.orep[y][1]))))
            for b in ry:
                if x in self.orep:
                    need.add(tuple(sorted((b, self.orep[x][1]))))
        need = sorted(p for p in need if p[0] != p[1])
        self.n_rep_pairs = len(need)
        print(f"[{self.r.version}] {len(og_pairs)} OG pairs, {len(el)} eligible SOGs, {len(self.srep)} with rep, "
              f"{len(self.orep)} OG reps, {len(need)} rep pairs to align", flush=True)
        t0 = time.time()
        res = run_foldseek(need, self.fsl, self.work / "foldseek") if need else pd.DataFrame()
        self.t_foldseek = time.time() - t0
        for r in res.itertuples(index=False):
            if r.status == "hit":
                self.aln[(r.q, r.t)] = (r.cigar, r.qstart, r.tstart, r.qlen, r.tlen, r.alntmscore, r.fident)
            else:
                self.aln[(r.q, r.t)] = (None,) * 7
        print(f"[{self.r.version}] foldseek {len(res)} pairs, {sum(v[0] is not None for v in self.aln.values())} aligned "
              f"({self.t_foldseek:.0f}s)", flush=True)

    # ---- links
    def links(self, x, y):
        rows = []
        for sx in self.el_by_og[x]:
            if sx not in self.srep:
                continue
            ax, mx = self.srep[sx]
            Fx = np.asarray(self.r.fp_cols[sx], np.int64)
            for sy in self.el_by_og[y]:
                if sy not in self.srep:
                    continue
                ay, my = self.srep[sy]
                Fy = np.asarray(self.r.fp_cols[sy], np.int64)
                fw = self.amap(mx, my)
                bw = self.amap(my, mx)
                if fw is None or bw is None:
                    rows.append({"og_x": x, "og_y": y, "sog_x": sx, "sog_y": sy, "rep_x": ax, "rep_y": ay,
                                 "aligned": False})
                    continue
                c, nu = project_footprint(Fx, self.colpos(x, ax), fw[0], self.pos2col(y, ay))
                jxy = directional_jaccard(c, nu, Fy)
                c2, nu2 = project_footprint(Fy, self.colpos(y, ay), bw[0], self.pos2col(x, ax))
                jyx = directional_jaccard(c2, nu2, Fx)
                rows.append({"og_x": x, "og_y": y, "sog_x": sx, "sog_y": sy, "rep_x": ax, "rep_y": ay, "aligned": True,
                             "alntm": fw[1], "fident": fw[2], "J_xy": jxy, "J_yx": jyx, "J": (jxy + jyx) / 2})
        return rows

    # ---- projections for test (b)
    def projections(self, x, y):
        """SOGs of x projected into y via rep(s_x) -> rep(OG_y)."""
        out = []
        if y not in self.orep:
            return out
        ay, my = self.orep[y]
        for sx in self.el_by_og[x]:
            if sx not in self.srep:
                continue
            ax, mx = self.srep[sx]
            fw = self.amap(mx, my)
            if fw is None:
                out.append({"og_x": x, "og_y": y, "sog_x": sx, "rep_x": ax, "rep_y": ay, "aligned": False})
                continue
            Fx = np.asarray(self.r.fp_cols[sx], np.int64)
            c, nu = project_footprint(Fx, self.colpos(x, ax), fw[0], self.pos2col(y, ay))
            out.append({"og_x": x, "og_y": y, "sog_x": sx, "rep_x": ax, "rep_y": ay, "aligned": True, "alntm": fw[1],
                        "coverage": 1 - nu / max(len(Fx), 1), "proj_cols": c.astype(np.int32)})
        return out


# ---------------------------------------------------------------------------------------------- states at projections
_SPTM = {}


def sptm_masks(accs, n_res):
    from pocketevo import states
    if not _SPTM:
        u = pd.read_parquet(P / "data/annotations/uniprot_proteins.parquet", columns=["accession", "signal", "transmem"])
        _SPTM.update({r.accession: (r.signal, r.transmem) for r in u.itertuples(index=False)})
    out = []
    for a_, n in zip(accs, n_res):
        if a_ not in _SPTM:
            out.append(None)
            continue
        sig, tm = _SPTM[a_]
        segs = ";".join(x for x in [sig, tm] if isinstance(x, str) and x)
        out.append(states.segments_mask(segs, int(n)))
    return out


def states_at(run: Run, og: str, F_list: list[np.ndarray], params: dict) -> pd.DataFrame:
    """compute_states for projected footprints (rows of F) on OG `og`; returns sog_index (= index in F_list) x member."""
    import scipy.sparse as sps
    from build_sogs import _member_table
    from pocketevo import sog, states
    ids, mat = run.msa(og)
    aln = {a: bytes(r).decode() for a, r in zip(ids, mat)}
    msa = sog.msa_from_aligned(aln, order=sorted(aln))
    mem = _member_table(None, og, run.members[run.members.og_id == og], run._models.get(og, pd.DataFrame(
        columns=["og_id", "accession", "n_res", "sequence", "seq_md5", "crc64", "crc_ok", "fragment_conflict", "plddt_q",
                 "parse_error"])))
    mw = mem[mem.has_model].set_index("accession").loc[msa.ids].reset_index()
    mw["og_median_len"] = float(mw.n_res.median())
    mw["sptm"] = sptm_masks(mw.accession, mw.n_res)
    s = run._sites.get(og, pd.DataFrame(columns=["og_id", "accession", "site_id", "residues", "site_index"]))
    s = s[s.accession.isin(set(msa.ids))].sort_values(["accession", "site_index"])
    row_of = {k: i for i, k in enumerate(msa.ids)}
    sm = s.accession.map(row_of).to_numpy(np.int64)
    S, _ = sog.site_columns(msa, sm, list(s.residues))
    rows = np.repeat(np.arange(len(F_list)), [len(f) for f in F_list])
    cols = np.concatenate(F_list) if F_list else np.zeros(0, np.int64)
    F = sps.csr_matrix((np.ones(len(cols), np.int32), (rows, cols)), shape=(len(F_list), msa.n_cols))
    st = states.compute_states(msa, mw, F, S, sm, s.site_id.to_numpy(), np.zeros(0, np.int64), TAU, params)
    tax = dict(zip(mem.accession, mem.taxid))
    st["taxid"] = st.accession.map(tax)
    st["model_accession"] = st.accession.map(dict(zip(mem.accession, mem.model_accession)))
    npk = dict(zip(mem.accession, (mem.n_pockets_all.fillna(0) - mem.n_pockets_kept.fillna(0)).astype(int)))
    st["n_filtered_pockets"] = st.accession.map(npk)
    # member footprint residues at the projected columns (for rescue)
    pos = {a: row_positions(mat[ids.index(a)]) for a in ids}
    st["fp_residues"] = [pos[a][F_list[i]][pos[a][F_list[i]] > 0] for i, a in zip(st.sog_index, st.accession)]
    return st


# ---------------------------------------------------------------------------------------------- rescue (k = 1)
def rescue_rows(rows: pd.DataFrame, work: Path, workers: int = 40) -> pd.DataFrame:
    """rows: key, model_accession, fp_residues. Returns rows + n_settings_hit, n_settings_run, rescued_state."""
    from pocketevo import rescue as R
    from rescue_sogs import extract
    from build_sogs import model_index
    if len(rows) == 0:
        return rows.assign(n_settings_hit=[], n_settings_run=[], rescued_state=[])

    class _C:
        pass
    ctx = _C()
    ctx.snapshot = SNAP
    ctx.a = argparse.Namespace(tierl_entries="/tmp/pe_big/tierL/fc_entries_tierL.parquet")
    ix = model_index(ctx, [1, 2, 3, 4])
    ix = ix[ix.accession.isin(set(rows.model_accession))]
    nfrag = ix.groupby("accession").fragment.max()
    multi = set(nfrag.index[nfrag > 1])
    one = ix[(ix.fragment == 1) & ~ix.accession.isin(multi)].rename(columns={"accession": "model_accession"})
    paths = extract(one, work / "models")
    todo = rows[rows.model_accession.isin(paths.keys())]
    hits = {k: 0 for k in todo.key}
    nrun = {k: 0 for k in todo.key}
    open_keys = set(todo.key)
    fps = dict(zip(rows.key, rows.fp_residues))
    for name in R.SETTING_NAMES:
        g = todo[todo.key.isin(open_keys)]
        models = sorted(set(g.model_accession))
        if not models:
            break
        t2 = time.time()
        outs = R.run_settings([str(paths[m]) for m in models], work / "burrow", workers=workers,
                              settings={name: R.SETTINGS[name]})
        n_hit = 0
        for key, macc in zip(g.key, g.model_accession):
            h = R.setting_hit(R.read_burrow_csv(Path(outs[name]) / f"{R.model_stem(paths[macc])}.csv"),
                              set(int(x) for x in fps[key]), 0.3, "jaccard")
            nrun[key] += 1
            if h["hit"]:
                hits[key] += 1
                open_keys.discard(key)
                n_hit += 1
        print(f"  rescue {name}: {len(models)} models, {n_hit}/{len(g)} hit, {time.time() - t2:.0f}s", flush=True)
    out = rows.copy()
    out["n_settings_hit"] = out.key.map(hits)
    out["n_settings_run"] = out.key.map(nrun)
    out["rescued_state"] = np.where(out.n_settings_hit.isna(), "U", np.where(out.n_settings_hit.fillna(0) >= 1, "W", "A"))
    return out


# ---------------------------------------------------------------------------------------------- fingerprints
def load_fingerprints(site_ids: set) -> dict:
    import pyarrow as pa
    import pyarrow.dataset as ds
    files = sorted(glob.glob(f"{SITE_FP}/tier*/*.parquet"))
    d = ds.dataset(files, format="parquet")
    t = d.to_table(columns=["site_id", "fingerprint"], filter=ds.field("site_id").isin(pa.array(sorted(site_ids))))
    sid = t.column("site_id").to_pylist()
    fp = t.column("fingerprint").combine_chunks()
    vals = fp.values.to_numpy(zero_copy_only=False).astype(np.float32).reshape(len(sid), -1)
    return dict(zip(sid, vals))


def select_members(Prow: pd.DataFrame) -> pd.DataFrame:
    """<= MAX_MEMBERS P members per SOG, one per species, highest footprint pLDDT."""
    g = Prow.sort_values(["sog_id", "footprint_plddt", "accession"], ascending=[True, False, True])
    g = g.drop_duplicates(["sog_id", "taxid"])
    return g.groupby("sog_id").head(MAX_MEMBERS)


def msa_identity(mat: np.ndarray, i: int, j: int) -> float:
    a, b = mat[i], mat[j]
    both = (a != GAP) & (b != GAP)
    n = both.sum()
    return float((a[both] == b[both]).sum() / n) if n else np.nan


def build_null_a(run: Run, cache: Path) -> pd.DataFrame:
    if cache.exists():
        return pd.read_parquet(cache)
    rng = np.random.default_rng(0)
    el = np.array(run.eligible)
    pick = rng.choice(len(el), size=min(30000, len(el)), replace=False)
    sogs = set(el[pick])
    sel = select_members(run.P[run.P.sog_id.isin(sogs)])
    sel = sel[sel.groupby("sog_id").accession.transform("size") >= 2]
    pairs = []
    for s, g in sel.groupby("sog_id"):
        a = g.accession.tolist(); sid = g.site_id.tolist()
        for i in range(len(a)):
            for j in range(i + 1, len(a)):
                pairs.append((s, a[i], a[j], sid[i], sid[j]))
    pairs = pd.DataFrame(pairs, columns=["sog_id", "acc_i", "acc_j", "site_i", "site_j"])
    if len(pairs) > NULL_PAIRS:
        pairs = pairs.sample(n=NULL_PAIRS, random_state=0)
    pairs["og_id"] = pairs.sog_id.str.split(":").str[0]
    run.load_shards(set(pairs.og_id), what=("msa",))
    ident = []
    for og, ai, aj in zip(pairs.og_id, pairs.acc_i, pairs.acc_j):
        ids, mat = run.msa(og)
        ident.append(msa_identity(mat, ids.index(ai), ids.index(aj)))
    pairs["identity"] = ident
    fps = load_fingerprints(set(pairs.site_i) | set(pairs.site_j))
    ok = pairs.site_i.isin(fps.keys()) & pairs.site_j.isin(fps.keys())
    pairs = pairs[ok].copy()
    A = np.stack([fps[s] for s in pairs.site_i]); B = np.stack([fps[s] for s in pairs.site_j])
    A /= np.linalg.norm(A, axis=1, keepdims=True); B /= np.linalg.norm(B, axis=1, keepdims=True)
    pairs["d"] = 1 - (A * B).sum(1)
    pairs = pairs.dropna(subset=["identity"])
    cache.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_parquet(cache, index=False)
    return pairs


def test_a(run: Run, links: pd.DataFrame, null: pd.DataFrame, fps: dict, msa_cons=None) -> pd.DataFrame:
    """Test (a) on tested links (best-partner union). Adds D, W_x, W_y, n_x, n_y, p_a."""
    sel = select_members(run.P[run.P.sog_id.isin(set(links.sog_x) | set(links.sog_y))])
    sel = sel[sel.site_id.isin(fps.keys())]
    mem = {s: np.stack([fps[x] for x in g.site_id]) for s, g in sel.groupby("sog_id")}
    rows = []
    for r in links.itertuples(index=False):
        A, B = mem.get(r.sog_x), mem.get(r.sog_y)
        if A is None or B is None or len(A) < 2 or len(B) < 2:
            rows.append((np.nan, np.nan, np.nan, 0 if A is None else len(A), 0 if B is None else len(B)))
            continue
        D = float(np.median(cos_dist(A, B)))
        wa = cos_dist(A, A)[np.triu_indices(len(A), 1)]
        wb = cos_dist(B, B)[np.triu_indices(len(B), 1)]
        rows.append((D, float(np.median(wa)), float(np.median(wb)), len(A), len(B)))
    out = links.copy()
    out[["D", "W_x", "W_y", "n_x", "n_y"]] = pd.DataFrame(rows, index=out.index)
    m = out.D.notna()
    qb, nb = bin_index(out.loc[m, "fident"].to_numpy(), null.identity.to_numpy())
    out["p_a"] = np.nan
    out.loc[m, "p_a"] = empirical_p(out.loc[m, "D"].to_numpy(), qb, null.d.to_numpy(), nb)
    out["id_bin"] = np.nan
    out.loc[m, "id_bin"] = qb
    return out


def finalize_a(out: pd.DataFrame, p_cut: float | None = None) -> tuple[pd.DataFrame, float]:
    m = out.p_a.notna()
    out["q_a"] = np.nan
    out.loc[m, "q_a"] = bh(out.loc[m, "p_a"].to_numpy())
    if p_cut is None:
        passed = out.loc[m & (out.q_a <= BH_Q), "p_a"]
        p_cut = float(passed.max()) if len(passed) else 0.0
    out["call_a"] = m & (out.p_a <= p_cut) & (out.D > out[["W_x", "W_y"]].max(axis=1))
    return out, p_cut


def chemistry(run: Run, co: CrossOG, r) -> int:
    """# mapped footprint positions where each side's members are >= 80% identical and consensus residues differ."""
    try:
        x, y = r.og_x, r.og_y
        ids_x, mat_x = run.msa(x); ids_y, mat_y = run.msa(y)
        fw = co.amap(co.srep[r.sog_x][1], co.srep[r.sog_y][1])
        cpx = co.colpos(x, r.rep_x)
        p2c_y = co.pos2col(y, r.rep_y)
        n = 0
        for c in np.asarray(run.fp_cols[r.sog_x], np.int64):
            rx = cpx[c]
            if rx == 0 or rx >= len(fw[0]) or fw[0][rx] == 0:
                continue
            cy = p2c_y[fw[0][rx] - 1]
            colx, coly = mat_x[:, c], mat_y[:, cy]
            colx, coly = colx[colx != GAP], coly[coly != GAP]
            if len(colx) == 0 or len(coly) == 0:
                continue
            vx, cx_ = np.unique(colx, return_counts=True); vy, cy_ = np.unique(coly, return_counts=True)
            if cx_.max() / len(colx) >= 0.8 and cy_.max() / len(coly) >= 0.8 and vx[cx_.argmax()] != vy[cy_.argmax()]:
                n += 1
        return n
    except Exception:
        return -1


# ---------------------------------------------------------------------------------------------- one full run
def run_pairs(run: Run, og_pairs: pd.DataFrame, fsl: dict, work: Path, null: pd.DataFrame, params: dict,
              shuffle_seed=None, do_rescue=True, p_cut=None, workers=40, tag="real"):
    t0 = time.time()
    co = CrossOG(run, fsl, work / tag, shuffle_seed=shuffle_seed)
    co.prepare(og_pairs)
    link_rows, proj_rows = [], []
    for x, y in zip(og_pairs.og_x, og_pairs.og_y):
        link_rows += co.links(x, y)
        proj_rows += co.projections(x, y)
        proj_rows += co.projections(y, x)
    L = pd.DataFrame(link_rows)
    if len(L) == 0:
        L = pd.DataFrame(columns=["og_x", "og_y", "sog_x", "sog_y", "rep_x", "rep_y", "aligned", "alntm", "fident", "J_xy",
                                  "J_yx", "J"])
    for c in ["alntm", "fident", "J_xy", "J_yx", "J"]:
        if c not in L:
            L[c] = np.nan
    L["supported"] = L.aligned.astype(bool) & (L.alntm >= ALNTM_MIN)
    L["linked"] = L.supported & (L.J >= TAU - 1e-12)
    # best partner per SOG (both sides), among supported pairs with J >= tau
    lk = L[L.linked].copy()
    lk["pair"] = lk.og_x + "|" + lk.og_y
    bx = lk.sort_values(["sog_x", "pair", "J", "sog_y"], ascending=[True, True, False, True]).drop_duplicates(["sog_x", "pair"])
    by = lk.sort_values(["sog_y", "pair", "J", "sog_x"], ascending=[True, True, False, True]).drop_duplicates(["sog_y", "pair"])
    tested = pd.concat([bx, by]).drop_duplicates(["sog_x", "sog_y"]).drop(columns=["pair"])
    t_link = time.time()
    # test (a)
    fps = load_fingerprints(set(select_members(run.P[run.P.sog_id.isin(set(tested.sog_x) | set(tested.sog_y))]).site_id)) \
        if len(tested) else {}
    A = test_a(run, tested, null, fps) if len(tested) else tested.assign(D=[], W_x=[], W_y=[], n_x=[], n_y=[], p_a=[], id_bin=[])
    A, p_cut = finalize_a(A, p_cut)
    if len(A):
        A["chem_n_sdp"] = [chemistry(run, co, r) for r in A.itertuples(index=False)]
    # test (b): projections
    PJ = pd.DataFrame(proj_rows)
    if len(PJ) == 0:
        PJ = pd.DataFrame(columns=["og_x", "og_y", "sog_x", "rep_x", "rep_y", "aligned", "alntm", "coverage", "proj_cols"])
    for c in ["alntm", "coverage"]:
        if c not in PJ:
            PJ[c] = np.nan
    PJ["testable"] = PJ.aligned.astype(bool) & (PJ.alntm >= ALNTM_MIN) & (PJ.coverage >= COVER_MIN - 1e-12)
    PJ["n_P_species_x"] = PJ.sog_x.map(run.nPsp).fillna(0).astype(int)
    st_all = []
    T = PJ[PJ.testable].reset_index(drop=True)
    for og, g in T.groupby("og_y"):
        st = states_at(run, og, [np.asarray(c, np.int64) for c in g.proj_cols], params)
        st["proj_idx"] = g.index.to_numpy()[st.sog_index.to_numpy()]
        st_all.append(st)
    ST = pd.concat(st_all, ignore_index=True) if st_all else pd.DataFrame(columns=["proj_idx", "state", "taxid"])
    pre = []
    for i, g in ST.groupby("proj_idx"):
        c = call_b(int(T.at[i, "n_P_species_x"]), g)
        pre.append({"proj_idx": i, **{k + "_pre": v for k, v in c.items() if k != "n_P_species_x"}})
    pre = pd.DataFrame(pre) if pre else pd.DataFrame(columns=["proj_idx", "n_PW_y_pre", "n_A_species_y_pre", "call_b_pre"])
    T = T.merge(pre, left_index=True, right_on="proj_idx", how="left").set_index("proj_idx", drop=False)
    t_b = time.time()
    # rescue A rows of pre-rescue passes
    resc = None
    if do_rescue and len(T) and T.call_b_pre.fillna(False).any():
        passes = set(T.index[T.call_b_pre.fillna(False).astype(bool)])
        Ar = ST[ST.proj_idx.isin(passes) & (ST.state == "A")].copy()
        Ar["key"] = list(zip(Ar.proj_idx, Ar.accession))
        print(f"[{run.version}/{tag}] rescue: {len(passes)} pre-rescue (b) passes, {len(Ar)} A rows, "
              f"{Ar.model_accession.nunique()} models", flush=True)
        resc = rescue_rows(Ar[["key", "model_accession", "fp_residues"]], work / "rescue", workers=workers)
        rs = dict(zip(resc.key, resc.rescued_state))
        m = ST.proj_idx.isin(passes) & (ST.state == "A")
        ST.loc[m, "state"] = [rs.get(k, "U") for k in zip(ST.loc[m, "proj_idx"], ST.loc[m, "accession"])]
    post = []
    for i, g in ST.groupby("proj_idx"):
        c = call_b(int(T.at[i, "n_P_species_x"]), g)
        post.append({"proj_idx": i, "n_PW_y": c["n_PW_y"], "n_A_species_y": c["n_A_species_y"], "call_b": c["call_b"]})
    post = pd.DataFrame(post) if post else pd.DataFrame(columns=["proj_idx", "n_PW_y", "n_A_species_y", "call_b"])
    T = T.drop(columns=["proj_idx"]).reset_index().merge(post, on="proj_idx", how="left")
    if not do_rescue:
        T["call_b"] = T.call_b_pre
    T["rescued"] = bool(do_rescue)
    timing = {"foldseek_s": co.t_foldseek, "links_s": t_link - t0 - co.t_foldseek, "b_states_s": t_b - t_link,
              "rescue_s": time.time() - t_b, "total_s": time.time() - t0, "n_rep_pairs": co.n_rep_pairs}
    return {"links": L, "tested": A, "proj": PJ.drop(columns=["proj_cols"]), "b": T.drop(columns=["proj_cols"]),
            "states_b": ST.drop(columns=["fp_residues"]), "rescue": resc, "p_cut": p_cut, "timing": timing, "co": co}


def funnel(res, og_pairs) -> dict:
    L, A, PJ, B = res["links"], res["tested"], res["proj"], res["b"]
    return {"og_pairs": int(len(og_pairs)),
            "sog_pairs_considered": int(len(L)), "sog_pairs_aligned": int(L.aligned.astype(bool).sum()) if len(L) else 0,
            "sog_pairs_supported_alntm": int(L.supported.sum()) if len(L) else 0,
            "sog_pairs_linked_J": int(L.linked.sum()) if len(L) else 0,
            "links_tested_best_partner": int(len(A)), "links_with_fingerprints": int(A.p_a.notna().sum()) if len(A) else 0,
            "call_a": int(A.call_a.sum()) if len(A) else 0,
            "projections": int(len(PJ)), "projections_testable": int(PJ.testable.sum()) if len(PJ) else 0,
            "b_pre_rescue": int(B.call_b_pre.fillna(False).astype(bool).sum()) if len(B) else 0,
            "call_b": int(B.call_b.fillna(False).astype(bool).sum()) if len(B) else 0,
            "og_pairs_with_call_a": int(A[A.call_a].groupby(["og_x", "og_y"]).ngroups) if len(A) else 0,
            "og_pairs_with_call_b": int(B[B.call_b.fillna(False).astype(bool)].assign(
                k=lambda d: [tuple(sorted(t)) for t in zip(d.og_x, d.og_y)]).k.nunique()) if len(B) else 0,
            "p_cut_a": res["p_cut"], **{k: v for k, v in res["timing"].items()}}


def git_commit():
    try:
        c = subprocess.run(["git", "-C", str(P), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        d = subprocess.run(["git", "-C", str(P), "status", "--porcelain", "--", "scripts/partB/crossog.py"],
                           capture_output=True, text=True).stdout.strip()
        return c + ("-dirty" if d else "")
    except Exception:
        return "unknown"


def load_params():
    import yaml
    from pocketevo import sog
    pre = yaml.safe_load(open(P / "config/prereg.yaml"))
    p = sog.prereg_params(pre)
    p["A_min_footprint_residues"] = 10   # D11
    return p


def family_tables(run: Run):
    ted = pd.read_parquet(P / "data/annotations/ted_domains.parquet", columns=["accession", "cath_label", "cath_level"])
    ted = ted[ted.cath_level == "H"][["accession", "cath_label"]]
    pf = pd.read_parquet(P / "data/annotations/pfam_regions.parquet", columns=["accession", "pfam", "clan"])
    mem = run.members[["og_id", "accession"]]
    labs = family_labels(mem, ted, pf)
    with_ted = {og for og, (t, _) in labs.items() if t}
    return labs, with_ted


def write(res: dict, out: Path, prefix: str):
    out.mkdir(parents=True, exist_ok=True)
    for k in ["links", "tested", "b"]:
        df = res[k]
        if k == "links":
            df = df[df.linked] if len(df) else df
        df.to_parquet(out / f"{prefix}_{k}.parquet", index=False)


# ============================================================================================ modes
def mode_validate(a):
    import yaml
    out = Path(a.out)
    work = Path(a.work) / ("validate" + a.subdir)
    fsl = load_fs_lookup()
    params = load_params()
    items = [i for i in yaml.safe_load(open(P / "scripts/cases/validation_controls.yaml"))["items"] if i["role"] == "known"]
    vt = pd.read_parquet(P / f"results/cases/{SNAP}/validation_table.parquet")
    summary, all_rows = [], []
    for ver in ["O1", "O2"]:
        run = Run(ver)
        labs, with_ted = family_tables(run)
        null = build_null_a(run, Path(a.work) / f"null_a_{ver}.parquet")
        print(f"[{ver}] null (a): {len(null)} within-OG pairs", flush=True)
        pairs, meta = [], {}
        for it in items:
            row = vt[(vt.item == it["id"]) & (vt.orthology == ver)].iloc[0]
            ogs = json.loads(row.ogs)
            dog = sorted({ogs[p] for p in it["derived"] if p in ogs})
            aog = sorted({ogs[p] for p in it["ancestral"] if p in ogs})
            prs = sorted({(d, x) for d in dog for x in aog if d != x})
            meta[it["id"]] = {"derived_ogs": dog, "ancestral_ogs": aog, "pairs": prs, "anchor": row.sog_id}
            pairs += prs
        og_pairs = pd.DataFrame(sorted(set(pairs)), columns=["og_x", "og_y"])
        og_pairs["family"] = [same_family(labs.get(x, (frozenset(), frozenset())), labs.get(y, (frozenset(), frozenset())),
                                          x in with_ted, y in with_ted) for x, y in zip(og_pairs.og_x, og_pairs.og_y)]
        res = run_pairs(run, og_pairs, fsl, work / ver, null, params, workers=a.workers, tag="real")
        write(res, out / ("validation" + a.subdir), f"{ver}")
        og_pairs.to_csv(out / ("validation" + a.subdir) / f"{ver}_og_pairs.csv", index=False)
        A, B = res["tested"], res["b"]
        for it in items:
            m = meta[it["id"]]
            prset = {tuple(p) for p in m["pairs"]}
            und = prset | {(y, x) for x, y in prset}
            in_pair = A.apply(lambda r: (r.og_x, r.og_y) in und, axis=1) if len(A) else pd.Series(dtype=bool)
            Ai = A[in_pair] if len(A) else A
            anchor = m["anchor"]
            anc = Ai[(Ai.sog_x == anchor) | (Ai.sog_y == anchor)] if len(Ai) else Ai
            rec_spec = bool(anc.call_a.any()) if len(anc) else False
            rec_spec_len = bool(Ai.call_a.any()) if len(Ai) else False
            Bc = B[B.call_b.fillna(False).astype(bool)] if len(B) else B
            dset, aset = set(m["derived_ogs"]), set(m["ancestral_ogs"])
            rec_loss = bool(((Bc.og_y.isin(dset)) & (Bc.og_x.isin(aset))).any()) if len(Bc) else False
            rec_gain = bool(((Bc.og_x.isin(dset)) & (Bc.og_y.isin(aset))).any()) if len(Bc) else False
            exp = it["expected_type"]
            rec = (("specificity" in exp and rec_spec) or ("loss" in exp and rec_loss) or ("gain" in exp and rec_gain))
            Bi = B[B.apply(lambda r: (r.og_x, r.og_y) in und, axis=1)] if len(B) else B
            fam = og_pairs[og_pairs.apply(lambda r: (r.og_x, r.og_y) in prset, axis=1)].family.tolist()
            summary.append({
                "item": it["id"], "orthology": ver, "expected": ";".join(exp), "in_scope": bool(prset),
                "derived_ogs": ";".join(m["derived_ogs"]), "ancestral_ogs": ";".join(m["ancestral_ogs"]),
                "family_rule": ";".join(str(f) for f in fam), "anchor_sog": anchor,
                "anchor_eligible": anchor in set(run.eligible),
                "anchor_has_rep": anchor in res["co"].srep,
                "n_links_tested": int(len(Ai)), "n_anchor_links": int(len(anc)),
                "anchor_best_J": float(anc.J.max()) if len(anc) else np.nan,
                "anchor_min_q_a": float(anc.q_a.min()) if len(anc) and anc.q_a.notna().any() else np.nan,
                "n_call_a_any_link": int(Ai.call_a.sum()) if len(Ai) else 0,
                "n_projections_testable": int(Bi.shape[0]),
                "n_b_pre_rescue": int(Bi.call_b_pre.fillna(False).astype(bool).sum()) if len(Bi) else 0,
                "n_call_b": int(Bi.call_b.fillna(False).astype(bool).sum()) if len(Bi) else 0,
                "rec_specificity_anchor": rec_spec, "rec_specificity_lenient": rec_spec_len,
                "rec_loss": rec_loss, "rec_gain": rec_gain, "recovered": bool(rec)})
        fn = funnel(res, og_pairs)
        json.dump(fn, open(out / ("validation" + a.subdir) / f"{ver}_funnel.json", "w"), indent=1, default=float)
        print(json.dumps(fn, default=float), flush=True)
    S = pd.DataFrame(summary)
    S.to_csv(out / ("validation" + a.subdir) / "validation_recall.csv", index=False)
    per_item = S.groupby("item").recovered.any()
    rec = {"recall_k": int(per_item.sum()), "n_items": int(len(per_item)),
           "recall_O1": int(S[S.orthology == "O1"].recovered.sum()), "recall_O2": int(S[S.orthology == "O2"].recovered.sum()),
           "in_scope_items": int(S.groupby("item").in_scope.any().sum()),
           "lenient_specificity_items": int(S.groupby("item").rec_specificity_lenient.any().sum()),
           "git_commit": git_commit(), "method": "docs/CROSSOG_METHOD.md @ 77db16c", "snapshot": SNAP,
           "label": "exploratory, not preregistered (D31)"}
    json.dump(rec, open(out / ("validation" + a.subdir) / "validation_recall.json", "w"), indent=1)
    print(S.to_string())
    print(json.dumps(rec, indent=1))


def same_family_pairs(run: Run, labs, with_ted, ogs_set) -> pd.DataFrame:
    rows = defaultdict(set)
    ogs = sorted(ogs_set)
    by_label = defaultdict(list)
    for og in ogs:
        t, pf = labs.get(og, (frozenset(), frozenset()))
        if og in with_ted:
            for l in t:
                by_label["TED:" + l].append(og)
        for l in pf:
            by_label["PF:" + l].append(og)
    for lab, lst in by_label.items():
        lst = sorted(set(lst))
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                x, y = lst[i], lst[j]
                if lab.startswith("TED:") or not (x in with_ted and y in with_ted):
                    rows[(x, y)].add(lab)
    return pd.DataFrame([(x, y, ";".join(sorted(l))) for (x, y), l in sorted(rows.items())],
                        columns=["og_x", "og_y", "family"])


def mode_genome(a):
    name = ("genome" if not a.max_pairs else f"characterization_sample{a.max_pairs}") + a.subdir
    out = Path(a.out) / name
    work = Path(a.work) / name
    fsl = load_fs_lookup()
    params = load_params()
    run = Run(a.version)
    labs, with_ted = family_tables(run)
    null = build_null_a(run, Path(a.work) / f"null_a_{a.version}.parquet")
    el_ogs = {s.split(":")[0] for s in run.eligible}
    pairs = same_family_pairs(run, labs, with_ted, el_ogs)
    if a.max_pairs:
        pairs = pairs.sample(n=min(a.max_pairs, len(pairs)), random_state=0).sort_values(["og_x", "og_y"])
    print(f"[genome] {len(pairs)} same-family OG pairs (both with >= 1 eligible SOG)", flush=True)
    out.mkdir(parents=True, exist_ok=True)
    pairs.to_parquet(out / f"{a.version}_og_pairs.parquet", index=False)
    rng_pairs = pairs.sample(n=min(a.null_pairs, len(pairs)), random_state=0).sort_values(["og_x", "og_y"])
    res = run_pairs(run, pairs, fsl, work, null, params, workers=a.workers, tag="real", do_rescue=not a.no_rescue)
    write(res, out, f"{a.version}_real")
    fn_real = funnel(res, pairs)
    print(json.dumps(fn_real, default=float), flush=True)
    # null: shuffled correspondence on the sampled pairs; real restricted to the same pairs for the per-pair rates
    resn = run_pairs(run, rng_pairs, fsl, work, null, params, shuffle_seed=0, workers=a.workers, tag="null",
                     do_rescue=not a.no_rescue, p_cut=res["p_cut"])
    write(resn, out, f"{a.version}_null")
    fn_null = funnel(resn, rng_pairs)
    key = set(zip(rng_pairs.og_x, rng_pairs.og_y))
    sub = {"links": res["links"][[k in key for k in zip(res["links"].og_x, res["links"].og_y)]] if len(res["links"]) else res["links"],
           "tested": res["tested"][[k in key for k in zip(res["tested"].og_x, res["tested"].og_y)]] if len(res["tested"]) else res["tested"],
           "proj": res["proj"][[tuple(sorted(k)) in key for k in zip(res["proj"].og_x, res["proj"].og_y)]] if len(res["proj"]) else res["proj"],
           "b": res["b"][[tuple(sorted(k)) in key for k in zip(res["b"].og_x, res["b"].og_y)]] if len(res["b"]) else res["b"],
           "p_cut": res["p_cut"], "timing": {}}
    fn_real_sub = funnel(sub, rng_pairs)
    json.dump({"real_all_pairs": fn_real, "real_on_null_sample": fn_real_sub, "null_shuffled_on_sample": fn_null,
               "git_commit": git_commit(), "method": "docs/CROSSOG_METHOD.md @ 77db16c", "snapshot": SNAP,
               "orthology": a.version, "label": "exploratory, not preregistered (D31)"},
              open(out / f"{a.version}_funnel.json", "w"), indent=1, default=float)
    print(json.dumps({"real_sub": fn_real_sub, "null": fn_null}, default=float, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["validate", "genome"])
    ap.add_argument("--out", default=str(P / f"results/partB_crossog/{SNAP}"))
    ap.add_argument("--work", default=f"/tmp/pe_big/partB_crossog/{SNAP}")
    ap.add_argument("--version", default="O1")
    ap.add_argument("--workers", type=int, default=40)
    ap.add_argument("--max-pairs", type=int, default=0)
    ap.add_argument("--null-pairs", type=int, default=2000)
    ap.add_argument("--no-rescue", action="store_true")
    ap.add_argument("--alntm-min", type=float, default=ALNTM_MIN,
                    help="POST HOC diagnostic only (frozen method: 0.5); output goes to a *_posthoc_alntm<x> dir")
    a = ap.parse_args()
    if a.alntm_min != ALNTM_MIN:
        globals()["ALNTM_MIN"] = a.alntm_min
        a.subdir = f"_posthoc_alntm{a.alntm_min:g}"
    else:
        a.subdir = ""
    if a.mode == "validate":
        mode_validate(a)
    else:
        mode_genome(a)


if __name__ == "__main__":
    main()
