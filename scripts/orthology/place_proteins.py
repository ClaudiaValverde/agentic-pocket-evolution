#!/usr/bin/env python
"""Place atlas proteins that are not O2 HOG members into the O2 HOG hierarchy by homology (OMAmer-like).

For each query: DIAMOND vs all O2 input proteins (HOG members and singletons) -> best hit (max bitscore) among
hits with e-value <= evalue and query coverage >= min_query_cov -> root HOG of the hit; deepest compatible HOG =
the deepest group on the hit's HOG path whose level (NCBI taxid) is an ancestor of (or equal to) the query's
taxon. If no group on the path is compatible (e.g. a eukaryotic query hitting a bacterial-only root HOG) the root
is kept with compatible=False. second_root_ratio = best bitscore in a *different* root HOG / best bitscore
(near 1 = ambiguous family assignment). Placement is presence-only evidence (never an absence).

Query sets (--set):
  tier2     snapshot tier-2 (Swiss-Prot) accessions; sequences from the release's uniprot_sprot.fasta.gz
  context   atlas accessions (any tier) of the O2 context taxa that are not O2 members; sequences from UniProt
            REST (uniprotkb/stream, organism_id)
  tier3     snapshot tier-3 accessions, streamed from a FASTA (--fasta, e.g. uniprot_trembl.fasta.gz) in chunks
            of placement.chunk_size; each chunk -> one parquet part (resumable; --max-chunks for testing)
Accessions that are O2 members are reported with method 'member' (no search).

Output: <out_dir>/placement/<set>/part-*.parquet  columns: accession, tier, taxid, method (member|diamond|no_hit),
  hit, pident, qcov, evalue, bitscore, root, root_name, deepest, deepest_level, deepest_level_taxid, compatible,
  n_roots_hit, second_root_ratio
Usage: .venv/bin/python scripts/orthology/place_proteins.py --set tier2 [--config config/orthology.yaml]
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import subprocess
import time
import urllib.parse
import urllib.request

import duckdb
import numpy as np
import pandas as pd

from _common import ANCESTRA_BIN, acc_of_header, atlas_accessions, fetch, load_config, log, read_fasta, snapshot_id
from pocketevo.taxonomy import Taxonomy

COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "evalue", "bitscore"]


def build_db(cfg, tdir) -> str:
    """DIAMOND db of all O2 input proteins (cached)."""
    W = os.path.join(cfg["work_dir"], "placement")
    os.makedirs(W, exist_ok=True)
    db = os.path.join(W, "o2_members")
    if os.path.exists(db + ".dmnd"):
        return db
    o2 = pd.read_csv(os.path.join(cfg["out_dir"], "o2_proteomes.tsv"), sep="\t")
    fa = os.path.join(W, "o2_members.fa")
    with open(fa, "w") as fh:
        for f in o2.file:
            for h, s in read_fasta(f):
                fh.write(f">{acc_of_header(h)}\n{s}\n")
    subprocess.run([os.path.join(ANCESTRA_BIN, "diamond"), "makedb", "--in", fa, "-d", db, "--quiet",
                    "--threads", str(cfg["placement"]["threads"])], check=True)
    os.remove(fa)
    return db


class Placer:
    def __init__(self, cfg, tdir):
        self.cfg = cfg
        h = pd.read_parquet(os.path.join(tdir, "hogs.parquet"))
        p = pd.read_parquet(os.path.join(tdir, "protein_hogs.parquet"))
        self.hog_name = h.hog_name.values
        self.hog_level = h.level.values
        self.hog_ltax = h.level_taxid.values
        self.path = dict(zip(p.protein, p.path))
        self.root = dict(zip(p.protein, p.root))
        self.tax = Taxonomy(cfg["taxdump_dir"], load_names=False)

    def assign_member(self, acc, taxid):
        path = self.path[acc]
        r = int(path[0]) if len(path) else -1
        d = int(path[-1]) if len(path) else -1
        return self._row(acc, taxid, "member", acc, 100.0, 1.0, 0.0, np.nan, r, d, True, 1, np.nan)

    def _row(self, acc, taxid, method, hit, pident, qcov, evalue, bits, r, d, comp, nroots, ratio):
        return dict(accession=acc, taxid=int(taxid), method=method, hit=hit, pident=pident, qcov=qcov,
                    evalue=evalue, bitscore=bits, root=r, root_name=self.hog_name[r] if r >= 0 else None,
                    deepest=d, deepest_level=self.hog_level[d] if d >= 0 else None,
                    deepest_level_taxid=int(self.hog_ltax[d]) if d >= 0 else -1, compatible=comp,
                    n_roots_hit=nroots, second_root_ratio=ratio)

    def assign_hits(self, hits: pd.DataFrame, qtax: dict) -> list[dict]:
        """Vectorised: best hit per query, its root, 2nd-root ratio, #roots; per-query loop only for the
        (short) taxonomy walk along the best hit's HOG path."""
        pc = self.cfg["placement"]
        h = hits[(hits.evalue <= pc["evalue"]) & (100 * hits.qcov >= pc["min_query_cov"])].copy()
        if h.empty:
            return []
        h = h.sort_values(["qseqid", "bitscore", "sseqid"], ascending=[True, False, True], kind="mergesort")
        h["root"] = h.sseqid.map(self.root).fillna(-1).astype("int64")
        # singletons count as their own "root" (keyed by accession)
        h["rkey"] = np.where(h.root >= 0, h.root.astype(str), "s:" + h.sseqid)
        best = h.drop_duplicates("qseqid", keep="first").set_index("qseqid")
        h["bkey"] = h.qseqid.map(best.rkey)
        oth = h[h.rkey != h.bkey].drop_duplicates("qseqid", keep="first").set_index("qseqid").bitscore
        nroots = h.drop_duplicates(["qseqid", "rkey"]).groupby("qseqid").size()
        out = []
        for q, b in best.iterrows():
            r0 = int(b.root)
            ratio = float(oth.get(q, 0.0) / b.bitscore)
            t = qtax.get(q, -1)
            d, comp = -1, False
            if r0 >= 0:
                lin = set(self.tax.lineage(t)) if t > 0 else set()
                path = self.path[b.sseqid]
                d = int(path[0])
                for x in reversed(list(path)):
                    if int(self.hog_ltax[x]) in lin:
                        d, comp = int(x), True
                        break
            out.append(self._row(q, t, "diamond", b.sseqid, float(b.pident), float(b.qcov), float(b.evalue),
                                 float(b.bitscore), r0, d, comp, int(nroots[q]), ratio))
        return out


def run_diamond(cfg, db, qfa, out_tsv):
    """DIAMOND blastp; the hit table is cached (out_tsv) so changing the HOG tables only re-assigns."""
    pc = cfg["placement"]
    t0 = time.time()
    if not os.path.exists(out_tsv):
        cmd = [os.path.join(ANCESTRA_BIN, "diamond"), "blastp", "-q", qfa, "-d", db, "-o", out_tsv + ".part",
               "--quiet", "--threads", str(pc["threads"]), "-e", str(pc["evalue"]), "-k", "25", "-b", "4", "-c", "1",
               "--outfmt", "6", *COLS]
        if pc.get("diamond_mode"):
            cmd.append(pc["diamond_mode"])
        subprocess.run(cmd, check=True)
        os.replace(out_tsv + ".part", out_tsv)
    h = pd.read_csv(out_tsv, sep="\t", header=None, names=COLS)
    h["qcov"] = (h.qend - h.qstart + 1) / h.qlen
    return h, time.time() - t0


def place_chunk(cfg, placer, db, recs: list[tuple], qtax: dict, qtier: dict, part_path: str, hits_path: str):
    """recs: (accession, sequence) of non-members."""
    os.makedirs(os.path.dirname(hits_path), exist_ok=True)
    qfa = hits_path + ".q.fa"
    if not os.path.exists(hits_path):
        with open(qfa, "w") as fh:
            for a, s in recs:
                fh.write(f">{a}\n{s}\n")
    h, sec = run_diamond(cfg, db, qfa, hits_path)
    if os.path.exists(qfa):
        os.remove(qfa)
    rows = placer.assign_hits(h, qtax)
    placed = {r["accession"] for r in rows}
    for a, _ in recs:
        if a not in placed:
            rows.append(placer._row(a, qtax.get(a, -1), "no_hit", None, np.nan, np.nan, np.nan, np.nan, -1, -1,
                                    False, 0, np.nan))
    df = pd.DataFrame(rows)
    df.insert(1, "tier", df.accession.map(qtier).astype("int8"))
    df.to_parquet(part_path, index=False)
    log(f"  chunk {os.path.basename(part_path)}: {len(recs)} queries, diamond {sec:.0f}s, "
        f"placed {len(placed)} ({len(placed) / max(len(recs), 1):.3f})")
    return sec


def stream_uniprot_organism(taxid: int, dest: str) -> str:
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {"query": f"organism_id:{taxid}", "format": "fasta", "compressed": "true"})
    return fetch(url, dest)


def fetch_uniparc_by_accession(accs: list[str], dest: str, threads: int = 4) -> str:
    """Sequences of (possibly deleted/merged) UniProtKB accessions from UniParc (dbid: query), cached FASTA."""
    if os.path.exists(dest):
        return dest
    from concurrent.futures import ThreadPoolExecutor

    def one(batch):
        url = "https://rest.uniprot.org/uniparc/stream?" + urllib.parse.urlencode(
            {"query": " OR ".join(f"dbid:{x}" for x in batch), "fields": "upi,accession,sequence", "format": "tsv"})
        for i in range(5):
            try:
                with urllib.request.urlopen(url, timeout=300) as r:
                    d = pd.read_csv(io.StringIO(r.read().decode()), sep="\t", dtype=str)
                want = set(batch)
                out = []
                for u, al, sq in zip(d.Entry, d.UniProtKB.fillna(""), d.Sequence):
                    for x in al.replace(",", ";").split(";"):
                        x = x.strip().split(".")[0]
                        if x in want:
                            out.append((x, sq))
                return out
            except Exception as e:  # noqa: BLE001
                log(f"uniparc retry {i}: {e}")
                time.sleep(3 * (i + 1))
        return []
    batches = [accs[i:i + 100] for i in range(0, len(accs), 100)]
    with ThreadPoolExecutor(threads) as ex:
        res = [r for b in ex.map(one, batches) for r in b]
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with gzip.open(dest + ".part", "wt") as fh:
        seen = set()
        for a, sq in res:
            if a not in seen:
                seen.add(a)
                fh.write(f">{a}\n{sq}\n")
    os.replace(dest + ".part", dest)
    return dest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config")
    ap.add_argument("--set", required=True, choices=["tier2", "context", "tier3"])
    ap.add_argument("--fasta", help="tier3: FASTA (gz ok) to stream sequences from")
    ap.add_argument("--tables", help="HOG table dir (default <out_dir>/tables)")
    ap.add_argument("--max-chunks", type=int, default=0)
    ap.add_argument("--force", action="store_true", help="re-assign chunks whose parquet exists (hits are reused)")
    ap.add_argument("--odir", help="output dir (default <out_dir>/placement/<set>)")
    ap.add_argument("--chunk-size", type=int, help="override placement.chunk_size")
    a = ap.parse_args()
    cfg = load_config(a.config)
    tdir = a.tables or os.path.join(cfg["out_dir"], "tables")
    odir = a.odir or os.path.join(cfg["out_dir"], "placement", a.set)
    if a.chunk_size:
        cfg["placement"]["chunk_size"] = a.chunk_size
    os.makedirs(odir, exist_ok=True)
    t0 = time.time()
    db = build_db(cfg, tdir)
    placer = Placer(cfg, tdir)
    members = set(placer.path)
    acc_pq = atlas_accessions(cfg)
    con = duckdb.connect()
    if a.set == "context":
        ctx = pd.read_csv(os.path.join(cfg["out_dir"], "context_species.tsv"), sep="\t")
        q = con.execute(f"""SELECT accession, min(tier) tier, any_value(taxid) taxid FROM '{acc_pq}'
                            WHERE taxid IN ({",".join(map(str, ctx.taxid))}) GROUP BY 1""").df()
    else:
        tier = 2 if a.set == "tier2" else 3
        q = con.execute(f"SELECT accession, min(tier) tier, any_value(taxid) taxid FROM '{acc_pq}' "
                        f"WHERE tier={tier} GROUP BY 1").df()
    qtax = dict(zip(q.accession, q.taxid))
    qtier = dict(zip(q.accession, q.tier))
    is_mem = q.accession.isin(members)
    log(f"{a.set}: {len(q)} atlas accessions, {int(is_mem.sum())} are O2 members")
    mem_rows = [placer.assign_member(r.accession, r.taxid) for r in q[is_mem].itertuples()]
    if mem_rows:
        df = pd.DataFrame(mem_rows)
        df.insert(1, "tier", df.accession.map(qtier).astype("int8"))
        df.to_parquet(os.path.join(odir, "part-members.parquet"), index=False)
    want = set(q.accession[~is_mem])

    # sequence source
    W = os.path.join(cfg["work_dir"], "placement")
    prov = {"set": a.set, "snapshot": snapshot_id(cfg["snapshot_manifest"]), "n_atlas": int(len(q)),
            "n_members": int(is_mem.sum()), "n_queries": len(want)}
    if a.set == "tier2":
        base = cfg["uniprot_refprot_base"].rsplit("/reference_proteomes", 1)[0] + "/complete"
        src = fetch(base + "/uniprot_sprot.fasta.gz", os.path.join(cfg["work_dir"], "raw", "sprot",
                                                                    "uniprot_sprot.fasta.gz"))
        rel = fetch(base + "/reldate.txt", os.path.join(cfg["work_dir"], "raw", "sprot", "reldate.txt"))
        prov["source"] = src
        prov["release"] = open(rel).read().strip().split("\n")[1]
        sources = [src]
    elif a.set == "context":
        sources = []
        for t in sorted(set(q.taxid[~is_mem])):
            sources.append(stream_uniprot_organism(int(t), os.path.join(W, "context_fa", f"{t}.fasta.gz")))
        prov["source"] = "rest.uniprot.org/uniprotkb/stream?query=organism_id:<taxid>&format=fasta"
        prov["retrieved"] = time.strftime("%Y-%m-%d")
    else:
        if not a.fasta:
            raise SystemExit("--set tier3 needs --fasta (e.g. uniprot_trembl.fasta.gz)")
        sources = [a.fasta]
        prov["source"] = a.fasta
    chunk = cfg["placement"]["chunk_size"]
    hkey = os.path.basename(a.fasta).split(".")[0] if a.fasta else a.set

    def hits_of(k, recs):   # cached DIAMOND hits (bulk, /tmp); keyed by set, source, chunk and query content
        import hashlib
        hq = hashlib.md5("\n".join(sorted(x for x, _ in recs)).encode()).hexdigest()[:12]
        return os.path.join(W, "hits", f"{a.set}_{hkey}_c{chunk}", f"part-{k:05d}_{hq}.tsv")
    buf, k, seen, dsec = [], 0, set(), 0.0
    for src in sources:
        for h, s in read_fasta(src):
            acc = acc_of_header(h)
            if acc in want and acc not in seen:
                seen.add(acc)
                buf.append((acc, s))
                if len(buf) >= chunk:
                    part = os.path.join(odir, f"part-{k:05d}.parquet")
                    if a.force or not os.path.exists(part):
                        dsec += place_chunk(cfg, placer, db, buf, qtax, qtier, part, hits_of(k, buf))
                    buf, k = [], k + 1
                    if a.max_chunks and k >= a.max_chunks:
                        break
        if a.max_chunks and k >= a.max_chunks:
            break
    if a.set != "tier3" and not (a.max_chunks and k >= a.max_chunks):
        # accessions without a sequence in the primary source (deleted / merged / not in this release):
        # sequence from UniParc (the sequence that was modelled in AFDB)
        rest = sorted(want - seen)
        if rest:
            up = fetch_uniparc_by_accession(rest, os.path.join(W, f"uniparc_{a.set}.fasta.gz"))
            prov["uniparc_fallback"] = "rest.uniprot.org/uniparc/stream?query=dbid:<acc>"
            for h, s in read_fasta(up):
                acc = acc_of_header(h)
                if acc in want and acc not in seen:
                    seen.add(acc)
                    buf.append((acc, s))
    if buf:
        part = os.path.join(odir, f"part-{k:05d}.parquet")
        if a.force or not os.path.exists(part):
            dsec += place_chunk(cfg, placer, db, buf, qtax, qtier, part, hits_of(k, buf))
    missing = sorted(want - seen)
    if a.set != "tier3":   # tier 3: only counted (a streamed FASTA covers a subset by design)
        pd.DataFrame({"accession": missing}).to_csv(os.path.join(odir, "no_sequence.tsv"), sep="\t", index=False)
    prov.update(n_with_sequence=len(seen), n_no_sequence=len(missing), diamond_s=round(dsec, 1),
                wall_s=round(time.time() - t0, 1), config=cfg["placement"])
    import glob
    df = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(os.path.join(odir, "part-*.parquet")))])
    prov["methods"] = df.method.value_counts().to_dict()
    d = df[df.method == "diamond"]
    prov["frac_compatible_of_diamond"] = float(d.compatible.mean()) if len(d) else None
    prov["frac_placed_in_hog"] = float((df.root >= 0).mean())
    json.dump(prov, open(os.path.join(odir, "placement.json"), "w"), indent=1, default=str)
    log(json.dumps(prov, default=str))


if __name__ == "__main__":
    main()
