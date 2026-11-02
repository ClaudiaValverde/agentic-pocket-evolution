"""Shared CLI plumbing for scripts/partA/*: snapshot, bulk work dir, compact results dir, logging."""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from pocketevo import pairs  # noqa: E402


def parser(desc: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=desc)
    ap.add_argument("--snapshot", default="20261030T085136Z")
    ap.add_argument("--sites", default=None, help="site dir (default /tmp/pe_big/sites/<snapshot>)")
    ap.add_argument("--work", default=None, help="bulk dir (default /tmp/pe_big/partA/<snapshot>)")
    ap.add_argument("--results", default=None, help="compact results (default results/partA/<snapshot>)")
    ap.add_argument("--orthology", default=str(REPO / "data/orthology/ancestra_tier1_v0/accurate"),
                    help="Ancestra run dir with hogs.orthoxml.gz (O1 default; O2 when available)")
    ap.add_argument("--annotations", default=str(REPO / "data/annotations"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="", help="suffix for classification-dependent outputs (e.g. _t0.3 sensitivity)")
    return ap


def resolve(a):
    a.sites = Path(a.sites or f"/tmp/pe_big/sites/{a.snapshot}")
    a.work = Path(a.work or f"/tmp/pe_big/partA/{a.snapshot}")
    a.results = Path(a.results or REPO / "results" / "partA" / a.snapshot)
    a.work.mkdir(parents=True, exist_ok=True)
    a.results.mkdir(parents=True, exist_ok=True)
    a.prereg_sha = pairs.check_prereg()
    a.prereg = pairs.load_prereg()
    return a


def logger(path: Path):
    fh = open(path, "a")

    def log(*msg):
        s = time.strftime("%H:%M:%S ") + " ".join(str(m) for m in msg)
        print(s, flush=True)
        fh.write(s + "\n")
        fh.flush()
    return log


def build_seqdb(a, acc, log=None) -> dict:
    """The one MMseqs sequence DB of a Part A work dir (<work>/mmseqs/seqdb), shared by 04b, 06, 09, 12.

    Model sequence (F1 of the AFDB v6 file the pockets were computed on, from <foldseek_db>.fasta) for every protein,
    UniProt sequence only where no single-file model exists (fragmented > 2,700 aa) or the model is missing.
    Built only if absent. Every stage must build it the same way: the first stage to run creates it and the others
    reuse it (G5: 04b used to build a UniProt-only DB first, which dropped 7,028 proteins without a UniProt
    sequence and changed the 06 classification)."""
    import subprocess
    import pandas as pd
    seqdb = Path(a.work) / "mmseqs" / "seqdb"
    if Path(str(seqdb) + ".dbtype").exists():
        return {}
    seqdb.parent.mkdir(exist_ok=True)
    model_seq = {}
    with open(str(a.foldseek_db) + ".fasta") as fh:
        name = None
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().replace("AF-", "").split("-F1")[0]
            else:
                model_seq[name] = line.strip()
    sq = pd.read_parquet(a.sequences)
    upm = pd.read_parquet(Path(a.annotations) / "uniprot_proteins.parquet", columns=["accession", "uniprot_acc"])
    sq = upm.merge(sq, on="uniprot_acc").drop_duplicates("accession")
    up_seq = dict(zip(sq.accession.values, sq.sequence.values))
    n_model = n_up = 0
    tmp = seqdb.parent / "all.fasta"
    with open(tmp, "w") as fh:
        for acc_ in acc:
            s_ = model_seq.get(acc_)
            u_ = up_seq.get(acc_)
            if isinstance(u_, str) and len(u_) > 2700 and (s_ is None or len(s_) == 1400):
                s_ = u_
            if s_ is None:
                s_ = u_
                n_up += isinstance(u_, str)
            else:
                n_model += 1
            if isinstance(s_, str) and s_:
                fh.write(f">{acc_}\n{s_}\n")
    subprocess.run([a.mmseqs, "createdb", tmp, seqdb], check=True, capture_output=True)
    if log:
        log(f"seqdb built: {n_model} model sequences, {n_up} UniProt fallback")
    return {"seqdb_model_seq": n_model, "seqdb_uniprot_fallback": n_up}
