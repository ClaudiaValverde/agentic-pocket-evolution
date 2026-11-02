#!/usr/bin/env python
"""Sequence-based orthology check for Part B candidates (PI D28 follow-up): reciprocal best hits (DIAMOND).

For every candidate (primary rows of results/partB/<snap>/<run>/candidates.parquet; gains and losses) and every
outgroup supporting member q (gain: A outgroup members; loss: P outgroup members):
  forward  best hit of q among the proteomes of the derived-clade species (whole proteomes of the analysis:
           O1 = 46 tier-1 proteomes, O2 = 157) must be a member of the derived clade of the same HOG tree;
  reverse  that member's best hit in q's own proteome must be q.
RBH pass = forward and reverse. q's own species is never a target (paralogous outgroups in a derived species are
not testable: rbh = None, excluded from the fraction). Domain architecture: ordered Pfam domains (data/annotations/pfam_regions.parquet)
of q vs its forward best hit; `arch_match` = identical ordered list (None when either is unannotated).
DIAMOND 2.2.8 (envs/ancestra), blastp --sensitive, --max-target-seqs 1, one DB per proteome
(/tmp/pe_big/partB/diamond_db). Output: results/partB/<snap>/<run>/rbh_outgroups.parquet (one row per candidate x
q) used by run_artifacts.py (orthology_split fails if < 50 % of the outgroup supporters pass RBH).

    python scripts/partB/rbh_check.py --runs O1_r3:O1 O2_r1:O2 --threads 16
"""
from __future__ import annotations

import argparse
import glob
import gzip
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()
DIAMOND = str(ROOT / "envs/ancestra/bin/diamond")
DBDIR = Path("/tmp/pe_big/partB/diamond_db")
WORK = Path("/tmp/pe_big/partB/rbh")


def acc_of(h: str) -> str:
    h = h.split()[0]
    return h.split("|")[1] if h.count("|") >= 2 else h


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur, buf = acc_of(line[1:]), []
            else:
                buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def run_diamond(pairs: dict, seq: dict, prot_file: dict, threads: int, tag: str) -> dict:
    """pairs: {species: set(query accessions)} -> {(query, species): (best subject, bitscore)}."""
    out = {}
    DBDIR.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    for sp, qs in sorted(pairs.items()):
        qs = sorted(q for q in qs if q in seq)
        if not qs:
            continue
        db = DBDIR / f"{sp}.dmnd"
        if not db.exists():
            subprocess.run([DIAMOND, "makedb", "--in", prot_file[sp], "-d", str(db), "--threads", str(threads),
                            "--quiet"], check=True)
        qf = WORK / f"{tag}_{sp}.faa"
        qf.write_text("".join(f">{q}\n{seq[q]}\n" for q in qs))
        of = WORK / f"{tag}_{sp}.tsv"
        subprocess.run([DIAMOND, "blastp", "-q", str(qf), "-d", str(db), "-o", str(of), "--sensitive",
                        "--max-target-seqs", "1", "--threads", str(threads), "--quiet",
                        "--outfmt", "6", "qseqid", "sseqid", "bitscore", "evalue"], check=True)
        if of.stat().st_size:
            t = pd.read_csv(of, sep="\t", header=None, names=["q", "s", "bits", "evalue"])
            t = t.sort_values("bits", ascending=False).drop_duplicates("q")
            for q, s, b in zip(t["q"], t["s"], t["bits"]):
                out[(acc_of(q), sp)] = (acc_of(s), float(b))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", default=["O1_r3:O1", "O2_r1:O2"])
    ap.add_argument("--threads", type=int, default=16)
    recon.add_partb_args(ap, skip=("sog_run", "orthology"))
    a = ap.parse_args()
    t0 = time.time()
    prot = pd.read_csv(ROOT / "data/orthology/ancestra_O2/o2_proteomes.tsv", sep="\t")
    prot["stem"] = [Path(f).name.split(".")[0] for f in prot["file"]]
    prot_file = dict(zip(prot["stem"], prot["file"]))
    tier1 = set(prot.loc[prot["set"] == "tier1", "stem"])
    tasks = []      # (run, tau, sog, node, type, q, q_species, derived set, derived species set)
    for spec in a.runs:
        run, orth = spec.split(":")
        cand = pd.read_parquet(ROOT / f"results/partB/{SNAP}/{run}/candidates.parquet")
        cand = cand[cand["primary"].astype(bool)].drop_duplicates(["tau", "sog_id", "node", "type"])
        ox = ROOT / ("data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz" if orth == "O2"
                     else "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz")
        trees = {t.tree_id: t for t in recon.iter_orthoxml_trees(ox, keep=set(cand["og_id"]))}
        stt = {}
        for tau in sorted(set(cand["tau"])):
            if orth == "O2":
                f = sorted(glob.glob(f"/tmp/pe_big/partB/{SNAP}/{run}/states_resc_t{tau}/shard=*.parquet"))
                st = ds.dataset(f, format="parquet").to_table(columns=["sog_id", "accession", "state"]).to_pandas()
            else:
                st = pd.read_parquet(ROOT / f"results/partB/events/{SNAP}/{run}/tau{tau}/primary/rescue/sog_states/bucket=all.parquet",
                                     columns=["sog_id", "accession", "state"])
            st = st[st["sog_id"].isin(set(cand["sog_id"]))]
            stt[tau] = dict(zip(zip(st["sog_id"], st["accession"]), st["state"]))
        for r in cand.itertuples():
            t = trees[r.og_id]
            names = t.leaf_names()
            lo, hi = t.leaf_lo[r.node], t.leaf_hi[r.node]
            der = set(names[lo:hi])
            dsp = {t.species[t.leaves[k]] for k in range(lo, hi)}
            sup_state = "A" if r.type == "gain" else "P"
            for k in list(range(0, lo)) + list(range(hi, t.n_leaves)):
                q = names[k]
                if stt[r.tau].get((r.sog_id, q)) == sup_state:
                    tasks.append((run, orth, r.tau, r.sog_id, int(r.node), r.type, q, t.species[t.leaves[k]], der, dsp))
    print(f"[{time.time() - t0:.0f}s] {len(tasks)} candidate x outgroup-supporter rows", flush=True)
    need = {t[6] for t in tasks} | {x for t in tasks for x in t[8]}
    seq = {}
    for stem, f in prot_file.items():
        s = read_fasta(f)
        seq.update({k: v for k, v in s.items() if k in need})
    print(f"[{time.time() - t0:.0f}s] sequences {len(seq)} / {len(need)}", flush=True)
    fwd_pairs = {}
    for t in tasks:
        for sp in t[9]:
            if (t[1] == "O1" and sp not in tier1) or sp == t[7]:   # q's own species is not a test target
                continue
            fwd_pairs.setdefault(sp, set()).add(t[6])
    fwd = run_diamond(fwd_pairs, seq, prot_file, a.threads, "fwd")
    print(f"[{time.time() - t0:.0f}s] forward done ({len(fwd)} hits)", flush=True)
    best = {}
    for t in tasks:
        cands = [(fwd[(t[6], sp)], sp) for sp in t[9] if (t[6], sp) in fwd and sp != t[7]]
        if cands:
            (s, b), sp = max(cands, key=lambda x: x[0][1])
            best[(t[6], frozenset(t[9]))] = (s, b, sp)
    rev_pairs = {}
    for t in tasks:
        bh = best.get((t[6], frozenset(t[9])))
        if bh:
            rev_pairs.setdefault(t[7], set()).add(bh[0])
    rev = run_diamond(rev_pairs, seq, prot_file, a.threads, "rev")
    print(f"[{time.time() - t0:.0f}s] reverse done ({len(rev)} hits)", flush=True)
    pf = pd.read_parquet(ROOT / "data/annotations/pfam_regions.parquet", columns=["accession", "pfam", "seq_start"])
    pf = pf[pf["accession"].isin(need | {b[0] for b in best.values()})]
    arch = pf.sort_values("seq_start").groupby("accession")["pfam"].apply(tuple).to_dict()
    rows = []
    for t in tasks:
        run, orth, tau, sog, node, typ, q, qsp, der, dsp = t
        bh = best.get((q, frozenset(dsp)))
        fwd_ok = bool(bh and bh[0] in der)
        rv = rev.get((bh[0], qsp)) if bh else None
        rev_ok = bool(rv and rv[0] == q)
        aq, ab = arch.get(q), arch.get(bh[0]) if bh else None
        testable = any(sp != qsp for sp in dsp)
        rows.append(dict(run=run, orthology=orth, tau=tau, sog_id=sog, node=node, type=typ, outgroup_member=q,
                         testable=testable,
                         outgroup_species=qsp, best_hit=bh[0] if bh else None, best_hit_bits=bh[1] if bh else None,
                         best_hit_species=bh[2] if bh else None, fwd_in_derived=fwd_ok,
                         reverse_best=rv[0] if rv else None, rbh=(fwd_ok and rev_ok) if testable else None,
                         arch_match=(aq == ab) if (aq and ab) else None, arch_q=";".join(aq or ()), arch_hit=";".join(ab or ())))
    df = pd.DataFrame(rows)
    for run in df["run"].unique():
        f = ROOT / f"results/partB/{SNAP}/{run}/rbh_outgroups.parquet"
        df[df["run"] == run].to_parquet(f, index=False)
        recon.write_sidecar(f, SNAP, 0, {"script": "scripts/partB/rbh_check.py", "argv": sys.argv[1:],
                                         "diamond": subprocess.run([DIAMOND, "version"], capture_output=True,
                                                                   text=True).stdout.strip(),
                                         "runtime_s": time.time() - t0})
    print(df.dropna(subset=["rbh"]).groupby(["run", "type"])["rbh"].apply(lambda x: x.astype(bool).mean()).round(3).to_dict(),
          f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
