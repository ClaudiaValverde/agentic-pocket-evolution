#!/usr/bin/env python3
"""Fetch unpaired MSAs from the ColabFold MMseqs2 server and convert them to AF3 JSON.

Route for proteins that have no cached AF3 (jackhmmer) MSA: non-human/non-yeast
orthologs, ancestral reconstructions, engineered variants.

Protocol (same as colabfold.mmseqs2.run_mmseqs2):
  POST {host}/ticket/msa   q=<fasta> mode=env|nomsa|...   -> {"id":..,"status":..}
  GET  {host}/ticket/{id}                                  -> poll until COMPLETE
  GET  {host}/result/download/{id}                         -> tar.gz with
       uniref.a3m, bfd.mgnify30.metaeuk30.smag30.a3m (env mode), pdb100.m8
The unpaired MSA handed to AF3 is uniref.a3m + envdb a3m concatenated, query first.
pairedMsa is "" (monomer + ligand) and templates [] (no-templates mode), matching
the primary co-folding protocol in DESIGN.md Sec.6.

Politeness: one submission per --min-interval seconds (default 20) enforced with a
lock file shared by all processes on the node; polling backs off 5 -> 30 s.
Everything is cached by sequence SHA-256 under --cache-dir, so a re-run is free.

Usage
  colabfold_msa.py --acc P24941                          # fetch (UniProt seq), cache
  colabfold_msa.py --acc ANC1 --seq-fasta anc.fa
  colabfold_msa.py --acc P24941 --ligand CCD:ATP --ligand CCD:MG \
      --name cdk2_atp_cf --out /tmp/af3_inputs/cdk2_atp_cf.json
  colabfold_msa.py --fasta many.fa --msa-only            # batch, cache only
  colabfold_msa.py --export-compact                      # copy gz a3m into data/msas/
  colabfold_msa.py --record-data-md                      # append provenance rows to DATA.md

Exit codes: 0 ok, 1 error, 4 server error/timeout.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import random
import sys
import tarfile
import threading
import time
import urllib.parse
import urllib.request
from pathlib import Path

P = Path(__file__).resolve().parents[2]  # repo root
HOST = "https://api.colabfold.com"
CACHE = Path("/tmp/pe_big/msas")
COMPACT = P / "data/msas"
COMPACT_BUDGET_MB = 50.0
UA = "pocket-evolution/0.1 (cofold MSA; contact: " + os.environ.get("COLABFOLD_CONTACT_EMAIL", "unset") + ")"
LOCKFILE = Path("/tmp/pe_big/msas/.colabfold_rate.lock")


def log(m: str) -> None:
    print(f"[colabfold_msa] {m}", file=sys.stderr, flush=True)


def sha(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def read_fasta_multi(path: str) -> list[tuple[str, str]]:
    out, name, cur = [], None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(cur).upper()))
            name, cur = line[1:].split()[0], []
        elif line:
            cur.append(line)
    if name is not None:
        out.append((name, "".join(cur).upper()))
    return out


def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc.upper()}.fasta"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}),
                                        timeout=30) as r:
                txt = r.read().decode()
            return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">")).upper()
        except Exception:  # noqa: BLE001
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


# ---------------------------------------------------------------- rate limiting
_SUBMIT_LOCK = threading.Lock()


def rate_limit(min_interval: float) -> None:
    """Block until at least min_interval s since the last submission on this node.

    The lock file spaces submissions across processes; the in-process lock keeps
    concurrent threads from racing past it.
    """
    LOCKFILE.parent.mkdir(parents=True, exist_ok=True)
    with _SUBMIT_LOCK:
        _rate_limit_locked(min_interval)


def _rate_limit_locked(min_interval: float) -> None:
    while True:
        last = float(LOCKFILE.read_text().strip()) if LOCKFILE.exists() else 0.0
        wait = last + min_interval - time.time()
        if wait <= 0:
            break
        log(f"rate limit: sleeping {wait:.0f}s")
        time.sleep(min(wait, min_interval))
    LOCKFILE.write_text(str(time.time()))


def post(url: str, fields: dict) -> dict:
    data = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(url, data=data, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def submit(seq: str, mode: str, host: str, min_interval: float) -> str:
    """Submit one query, retrying on RATELIMIT/MAINTENANCE. Returns the ticket id."""
    q = f">query\n{seq}\n"
    for attempt in range(12):
        rate_limit(min_interval)
        try:
            r = post(f"{host}/ticket/msa", {"q": q, "mode": mode})
        except Exception as e:  # noqa: BLE001
            log(f"submit failed ({e}); retrying")
            time.sleep(30)
            continue
        st = r.get("status")
        if st in ("PENDING", "RUNNING", "COMPLETE"):
            return r["id"]
        if st in ("RATELIMIT", "MAINTENANCE"):
            back = 60 + 30 * attempt + random.uniform(0, 10)
            log(f"server says {st}; backing off {back:.0f}s")
            time.sleep(back)
            continue
        raise RuntimeError(f"submit returned {r}")
    raise RuntimeError("submit: gave up after 12 attempts")


def poll(tid: str, host: str, timeout: float = 1800) -> None:
    t0, wait = time.time(), 5.0
    while True:
        try:
            st = get(f"{host}/ticket/{tid}").get("status")
        except Exception as e:  # noqa: BLE001
            log(f"poll error {e}")
            st = None
        if st == "COMPLETE":
            return
        if st in ("ERROR", "UNKNOWN"):
            raise RuntimeError(f"ticket {tid} status {st}")
        if time.time() - t0 > timeout:
            raise TimeoutError(f"ticket {tid} not complete after {timeout}s")
        time.sleep(wait)
        wait = min(wait * 1.4, 30.0)


def download(tid: str, host: str) -> bytes:
    req = urllib.request.Request(f"{host}/result/download/{tid}", headers={"User-Agent": UA})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                raise
            log(f"download failed ({e}); retry")
            time.sleep(15)
    raise RuntimeError("unreachable")


# ---------------------------------------------------------------- a3m handling
def merge_a3m(members: dict[str, str], query: str) -> tuple[str, dict]:
    """Concatenate ColabFold a3m files into one unpaired A3M with the query first.

    Keeps lowercase insertion columns (AF3 parses A3M). Drops the per-file query
    rows after the first, drops null-byte separators, and de-duplicates headers.
    """
    out_h, out_s, seen = [], [], set()
    out_h.append(">query")
    out_s.append(query)
    counts = {}
    for name in sorted(members):
        n = 0
        h, s = None, []
        lines = members[name].replace("\x00", "").split("\n")
        for line in lines:
            if line.startswith(">"):
                if h is not None:
                    n += _emit(h, "".join(s), out_h, out_s, seen, query)
                h, s = line.rstrip(), []
            elif line.strip():
                s.append(line.strip())
        if h is not None:
            n += _emit(h, "".join(s), out_h, out_s, seen, query)
        counts[name] = n
    a3m = "\n".join(f"{h}\n{s}" for h, s in zip(out_h, out_s)) + "\n"
    return a3m, counts


def _emit(h: str, s: str, out_h: list, out_s: list, seen: set, query: str) -> int:
    if not s:
        return 0
    ungapped = s.replace("-", "").upper()
    if ungapped == query.upper():  # per-file query row
        return 0
    key = (h, s)
    if key in seen:
        return 0
    seen.add(key)
    out_h.append(h)
    out_s.append(s)
    return 1


def neff_proxy(a3m: str, qlen: int) -> float:
    """Cheap depth proxy: number of sequences / sqrt(length). Not a true Neff."""
    n = a3m.count("\n>") + (1 if a3m.startswith(">") else 0)
    return round(n / max(1.0, qlen ** 0.5), 2)


# ---------------------------------------------------------------- main fetch
def fetch_msa(acc: str, seq: str, *, host: str, mode: str, cache: Path,
              min_interval: float, force: bool = False) -> dict:
    """Return the provenance dict for acc/seq, fetching from the server if needed."""
    h = sha(seq)
    d = cache / h
    prov_p = d / "provenance.json"
    if prov_p.exists() and not force:
        prov = json.loads(prov_p.read_text())
        log(f"{acc}: cache hit {d}")
        return prov
    d.mkdir(parents=True, exist_ok=True)
    log(f"{acc}: submitting {len(seq)} aa to {host} (mode={mode})")
    t0 = time.time()
    tid = submit(seq, mode, host, min_interval)
    log(f"{acc}: ticket {tid}; polling")
    poll(tid, host)
    blob = download(tid, host)
    (d / "raw.tar.gz").write_bytes(blob)
    members: dict[str, str] = {}
    extras: dict[str, int] = {}
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            data = tf.extractfile(m).read().decode(errors="replace")
            if m.name.endswith(".a3m"):
                members[m.name] = data
            else:
                extras[m.name] = m.size
    if not members:
        raise RuntimeError(f"{acc}: no .a3m in server tar (members: {list(extras)})")
    a3m, counts = merge_a3m(members, seq)
    with gzip.open(d / "unpaired.a3m.gz", "wt") as f:
        f.write(a3m)
    prov = {
        "acc": acc, "seq_sha256": h, "length": len(seq),
        "server": host, "endpoint": "/ticket/msa", "mode": mode, "ticket": tid,
        "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_s": round(time.time() - t0, 1),
        # The tar member names ARE the ColabFold DB version strings the server uses.
        "db_files": {k: counts[k] for k in sorted(counts)},
        "other_members": extras,
        "n_seqs": sum(counts.values()) + 1,
        "neff_proxy": neff_proxy(a3m, len(seq)),
        "a3m_gz_bytes": (d / "unpaired.a3m.gz").stat().st_size,
    }
    prov_p.write_text(json.dumps(prov, indent=1))
    log(f"{acc}: {prov['n_seqs']} seqs from {list(counts)} in {prov['elapsed_s']}s")
    return prov


def load_a3m(seq: str, cache: Path) -> str:
    with gzip.open(cache / sha(seq) / "unpaired.a3m.gz", "rt") as f:
        return f.read()


# ---------------------------------------------------------------- AF3 JSON
def parse_ligand(spec: str) -> dict:
    kind, _, val = spec.partition(":")
    if kind.upper() == "CCD":
        return {"ccdCodes": [c.strip().upper() for c in val.split("+")]}
    if kind.upper() == "SMILES":
        return {"smiles": val}
    raise SystemExit(f"ligand must be CCD:<code> or SMILES:<smiles>, got {spec}")


def write_af3_json(out: str, name: str, seq: str, a3m: str, ligands: list[str],
                   seeds: list[int], prov: dict, copies: int = 1) -> None:
    import string
    ids = iter(list(string.ascii_uppercase))
    prot_ids = [next(ids) for _ in range(copies)]
    seqs = [{"protein": {"id": prot_ids if copies > 1 else prot_ids[0], "sequence": seq,
                         "unpairedMsa": a3m, "pairedMsa": "", "templates": []}}]
    for spec in ligands:
        seqs.append({"ligand": {"id": next(ids), **parse_ligand(spec)}})
    job = {"dialect": "alphafold3", "version": 2, "name": name, "modelSeeds": seeds,
           "sequences": seqs, "bondedAtomPairs": None, "userCCD": None}
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(job))
    meta = {"acc": prov["acc"], "name": name, "length": len(seq), "copies": copies,
            "ligands": ligands, "seeds": seeds, "templates": 0,
            "msa": "colabfold_mmseqs2", "msa_provenance": prov,
            "sequence_source": "see caller"}
    Path(out).with_suffix(".meta.json").write_text(json.dumps(meta, indent=1))
    log(f"wrote {out} (name={name}, {len(a3m)//1024} kB MSA)")


# ---------------------------------------------------------------- side commands
def export_compact(cache: Path, compact: Path, budget_mb: float) -> None:
    """Write provenance JSON + sha256 for each cached MSA into the repo; payloads stay in /tmp.

    Home quota alert (2026-10-30): A3M payloads are no longer copied into data/msas/.
    They live at /tmp/pe_big/msas/<seq_sha256>/unpaired.a3m.gz and can be refetched from
    the server ticket recorded in DATA.md. `budget_mb` is kept for CLI compatibility.
    """
    import hashlib as _hl
    compact.mkdir(parents=True, exist_ok=True)
    man = compact / "MANIFEST.sha256"
    lines = set(man.read_text().splitlines()) if man.exists() else set()
    n = 0
    for pp in sorted(cache.glob("*/provenance.json")):
        prov = json.loads(pp.read_text())
        stem = f"{prov['acc']}_{prov['seq_sha256'][:12]}"
        dst = compact / f"{stem}.provenance.json"
        if not dst.exists():
            dst.write_text(json.dumps(prov, indent=1))
            n += 1
        src = pp.parent / "unpaired.a3m.gz"
        if src.exists():
            lines.add(f"{_hl.sha256(src.read_bytes()).hexdigest()}  {stem}.a3m.gz")
    man.write_text("\n".join(sorted(lines)) + "\n")
    log(f"provenance: {n} new records in {compact}; payloads stay in {cache}")


DATA_MD_MARK = "## ColabFold MMseqs2 MSAs (server provenance)"


def record_data_md(cache: Path, data_md: Path) -> None:
    """Append one row per fetched MSA to DATA.md (idempotent on ticket id)."""
    txt = data_md.read_text()
    if DATA_MD_MARK not in txt:
        txt += (f"\n{DATA_MD_MARK}\n\nUnpaired MSAs for proteins without a cached AF3 "
                "(jackhmmer) MSA. Fetched with `scripts/cofold/colabfold_msa.py`. The "
                "ColabFold server does not expose a version endpoint; the returned tar "
                "member names are the server's own database version strings and are "
                "recorded verbatim. License: ColabFold/MMseqs2 server is free for "
                "non-commercial and academic use (Mirdita et al. 2022, Nat Methods); "
                "underlying DBs UniRef30 (CC BY 4.0) and ColabFoldDB envdb.\n\n"
                "| acc | len | seq sha256[:12] | server | mode | ticket | fetched (UTC) | "
                "db files returned (n seqs) | n seqs | depth/sqrt(L) |\n"
                "|---|---|---|---|---|---|---|---|---|---|\n")
    rows = []
    for pp in sorted(cache.glob("*/provenance.json")):
        prov = json.loads(pp.read_text())
        if f"| {prov['ticket']} |" in txt:
            continue
        dbs = "; ".join(f"`{k}` ({v})" for k, v in prov["db_files"].items())
        rows.append(f"| {prov['acc']} | {prov['length']} | {prov['seq_sha256'][:12]} | "
                    f"{prov['server']} | {prov['mode']} | {prov['ticket']} | "
                    f"{prov['fetched_utc']} | {dbs} | {prov['n_seqs']} | {prov['neff_proxy']} |")
    if rows:
        txt = txt.rstrip("\n") + "\n" + "\n".join(rows) + "\n"
    data_md.write_text(txt)
    log(f"DATA.md: {len(rows)} new rows")


# ---------------------------------------------------------------- CLI
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--acc", help="UniProt accession or free label for the query")
    ap.add_argument("--seq")
    ap.add_argument("--seq-fasta")
    ap.add_argument("--fasta", help="multi-FASTA; batch mode, MSAs only")
    ap.add_argument("--ligand", action="append", default=[])
    ap.add_argument("--copies", type=int, default=1)
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("--name")
    ap.add_argument("--out", help="AF3 input JSON to write")
    ap.add_argument("--msa-only", action="store_true")
    ap.add_argument("--mode", default="env", choices=["env", "all", "nofilter", "noenv"],
                    help="env = uniref30 + ColabFoldDB envdb (default)")
    ap.add_argument("--host", default=HOST)
    ap.add_argument("--cache-dir", default=str(CACHE))
    ap.add_argument("--min-interval", type=float, default=20.0,
                    help="min seconds between submissions to the server (politeness)")
    ap.add_argument("--force", action="store_true", help="ignore the cache")
    ap.add_argument("--export-compact", action="store_true")
    ap.add_argument("--record-data-md", action="store_true")
    a = ap.parse_args()

    cache = Path(a.cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    mode = {"noenv": ""}.get(a.mode, a.mode)

    if a.fasta:
        for name, seq in read_fasta_multi(a.fasta):
            try:
                fetch_msa(name, seq, host=a.host, mode=mode, cache=cache,
                          min_interval=a.min_interval, force=a.force)
            except Exception as e:  # noqa: BLE001
                log(f"{name}: FAILED {e}")
    elif a.acc:
        if a.seq:
            seq = a.seq.strip().upper()
        elif a.seq_fasta:
            seq = read_fasta_multi(a.seq_fasta)[0][1]
        else:
            seq = fetch_uniprot(a.acc)
        prov = fetch_msa(a.acc, seq, host=a.host, mode=mode, cache=cache,
                         min_interval=a.min_interval, force=a.force)
        if a.out and not a.msa_only:
            name = a.name or f"{a.acc.lower()}_cf"
            write_af3_json(a.out, name, seq, load_a3m(seq, cache), a.ligand,
                           [int(s) for s in a.seeds.split(",") if s.strip()], prov, a.copies)

    if a.export_compact:
        export_compact(cache, COMPACT, COMPACT_BUDGET_MB)
    if a.record_data_md:
        record_data_md(cache, P / "DATA.md")
    if not (a.fasta or a.acc or a.export_compact or a.record_data_md):
        ap.error("give --acc, --fasta, --export-compact or --record-data-md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
