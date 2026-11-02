#!/usr/bin/env python3
"""Build an AlphaFold3 inference-only input JSON, reusing cached AF3 data JSONs.

The cached files (output of the AF3 data pipeline, one per UniProt accession)
live in read-only shared directories:
  /cluster/project/beltrao/shared/25.05_alphafold3_msas_human/<acc>_data.json.gz
  /cluster/project/beltrao/shared/25.06_alphafold3_msas_yeast/<acc>_data.json.gz
They contain the query sequence, unpairedMsa, pairedMsa and templates (mmCIF
text + query/template index mappings), so --norun_data_pipeline works directly.

Sequence safety: the sequence to fold is (in priority order) --seq / --seq-fasta,
else the current UniProt sequence (REST), else (only with --trust-cache-seq) the
cached sequence. The cached MSA/templates are reused ONLY when the cached query
sequence is identical to the sequence to fold. On mismatch the script exits
with code 2 and prints the differences; it never silently folds another sequence.

Exit codes: 0 ok, 2 sequence mismatch vs cache, 3 no cached MSA (MSA needed),
1 other error.

Examples:
  make_af3_input.py --acc P24941 --ligand CCD:ATP --ligand CCD:MG --name cdk2_atp
  make_af3_input.py --acc P00918 --ligand CCD:ZN --ligand 'SMILES:CC(=O)Nc1nnc(s1)S(N)(=O)=O'
  make_af3_input.py --acc P24941 --seq-fasta anc.fa --reuse-msa-for-variant ...
"""
from __future__ import annotations

import argparse
import gzip
import json
import string
import sys
import urllib.request
from pathlib import Path

DEFAULT_CACHE_DIRS = [
    "/cluster/project/beltrao/shared/25.05_alphafold3_msas_human",
    "/cluster/project/beltrao/shared/25.06_alphafold3_msas_yeast",
]
DEFAULT_SEEDS = [1, 2, 3]


def log(msg: str) -> None:
    print(f"[make_af3_input] {msg}", file=sys.stderr)


def read_fasta(path: str) -> str:
    lines = [l.strip() for l in Path(path).read_text().splitlines()]
    seqs, cur = [], []
    for l in lines:
        if l.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        elif l:
            cur.append(l)
    if cur:
        seqs.append("".join(cur))
    if len(seqs) != 1:
        raise SystemExit(f"{path}: expected exactly 1 sequence, got {len(seqs)}")
    return seqs[0].upper()


def fetch_uniprot(acc: str) -> str:
    import time
    url = f"https://rest.uniprot.org/uniprotkb/{acc.upper()}.fasta"
    for attempt in range(3):  # UniProt REST has occasional transient failures
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                txt = r.read().decode()
            break
        except Exception:  # noqa: BLE001
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
    if not txt.strip():
        raise RuntimeError(f"empty UniProt response for {acc}")
    return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">")).upper()


def find_cache(acc: str, cache_dirs: list[str]) -> Path | None:
    for d in cache_dirs:
        p = Path(d) / f"{acc.lower()}_data.json.gz"
        if p.exists():
            return p
    return None


def load_cached_protein(path: Path) -> dict:
    with gzip.open(path, "rt") as f:
        d = json.load(f)
    prots = [s["protein"] for s in d["sequences"] if "protein" in s]
    if len(prots) != 1:
        raise SystemExit(f"{path}: expected 1 protein chain, found {len(prots)}")
    return prots[0]


def describe_mismatch(a: str, b: str) -> str:
    out = [f"cached len={len(a)}, wanted len={len(b)}"]
    if len(a) == len(b):
        diffs = [f"{x}{i+1}{y}" for i, (x, y) in enumerate(zip(a, b)) if x != y]
        out.append(f"{len(diffs)} substitutions: {' '.join(diffs[:30])}"
                   + (" ..." if len(diffs) > 30 else ""))
    else:
        import difflib
        sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
        ops = [op for op in sm.get_opcodes() if op[0] != "equal"]
        out.append(f"identity-ish ratio={sm.ratio():.3f}; first edits: "
                   + "; ".join(f"{t} cached[{i1}:{i2}] -> wanted[{j1}:{j2}]"
                               for t, i1, i2, j1, j2 in ops[:10]))
    return "\n  ".join(out)


def replace_query_row(a3m: str, new_query: str) -> str:
    """Replace the first (query) sequence of an A3M string. Same length only."""
    lines = a3m.split("\n")
    assert lines[0].startswith(">"), "A3M must start with a header"
    j = 1
    while j < len(lines) and not lines[j].startswith(">"):
        j += 1
    old = "".join(lines[1:j])
    if len(old) != len(new_query):
        raise ValueError("query row length differs")
    return "\n".join([lines[0], new_query] + lines[j:])


def chain_ids():
    letters = string.ascii_uppercase
    for c in letters:
        yield c
    for c1 in letters:
        for c2 in letters:
            yield c1 + c2


def parse_ligand(spec: str) -> dict:
    kind, _, val = spec.partition(":")
    kind = kind.upper()
    if kind == "CCD":
        return {"ccdCodes": [c.strip().upper() for c in val.split("+")]}
    if kind == "SMILES":
        try:
            from rdkit import Chem  # available in the AF3 env
            if Chem.MolFromSmiles(val) is None:
                raise SystemExit(f"invalid SMILES: {val}")
        except ImportError:
            pass
        return {"smiles": val}
    raise SystemExit(f"ligand spec must be CCD:<code> or SMILES:<smiles>, got {spec}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--acc", required=True, help="UniProt accession")
    ap.add_argument("--seq", help="sequence override (to fold)")
    ap.add_argument("--seq-fasta", help="FASTA with the sequence to fold")
    ap.add_argument("--trust-cache-seq", action="store_true",
                    help="if no --seq and UniProt unreachable, fold the cached sequence")
    ap.add_argument("--no-uniprot", action="store_true",
                    help="do not query UniProt; requires --seq/--seq-fasta or --trust-cache-seq")
    ap.add_argument("--ligand", action="append", default=[],
                    help="CCD:ATP, CCD:ZN, CCD:NAP+..., or SMILES:<smiles>; repeatable, one chain each")
    ap.add_argument("--copies", type=int, default=1, help="homo-oligomer copies of the protein")
    ap.add_argument("--seeds", default=",".join(map(str, DEFAULT_SEEDS)),
                    help="comma-separated model seeds (default 1,2,3)")
    ap.add_argument("--no-templates", action="store_true", help="drop cached templates")
    ap.add_argument("--msa-free", action="store_true",
                    help="if no cache: write single-sequence input (empty MSA, no templates)")
    ap.add_argument("--reuse-msa-for-variant", action="store_true",
                    help="same-length variant (e.g. ancestral): reuse cached MSA/templates with "
                         "the query row replaced. Approximation, recorded in the output name/meta.")
    ap.add_argument("--cache-dir", action="append", help="override cache dirs (repeatable)")
    ap.add_argument("--name", help="job name (default: <acc>_<ligands>)")
    ap.add_argument("--out", required=True, help="output JSON path")
    args = ap.parse_args()

    acc = args.acc.upper()
    cache_dirs = args.cache_dir or DEFAULT_CACHE_DIRS
    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]

    # 1) sequence to fold
    wanted, src = None, None
    if args.seq:
        wanted, src = args.seq.strip().upper(), "override(--seq)"
    elif args.seq_fasta:
        wanted, src = read_fasta(args.seq_fasta), f"override({args.seq_fasta})"
    elif not args.no_uniprot:
        try:
            wanted, src = fetch_uniprot(acc), "uniprot_rest"
        except Exception as e:  # noqa: BLE001
            log(f"UniProt fetch failed: {e}")
    cache = find_cache(acc, cache_dirs)
    cached = load_cached_protein(cache) if cache else None
    if wanted is None:
        if args.trust_cache_seq and cached:
            wanted, src = cached["sequence"], "cache(--trust-cache-seq)"
        else:
            log("no sequence to fold: give --seq/--seq-fasta or --trust-cache-seq")
            return 1
    log(f"{acc}: sequence source={src}, length={len(wanted)}")

    meta = {"acc": acc, "sequence_source": src, "cache_file": str(cache) if cache else None}

    # 2) MSA/templates
    if cached is not None and cached["sequence"] == wanted:
        protein = {"sequence": wanted, "unpairedMsa": cached["unpairedMsa"],
                   "pairedMsa": cached["pairedMsa"],
                   "templates": [] if args.no_templates else cached["templates"]}
        meta["msa"] = "cached_exact_match"
        log(f"cache hit, exact sequence match ({cache}); templates={len(protein['templates'])}")
    elif cached is not None:
        msg = describe_mismatch(cached["sequence"], wanted)
        if args.reuse_msa_for_variant and len(cached["sequence"]) == len(wanted):
            protein = {"sequence": wanted,
                       "unpairedMsa": replace_query_row(cached["unpairedMsa"], wanted),
                       "pairedMsa": replace_query_row(cached["pairedMsa"], wanted)
                       if cached["pairedMsa"] else cached["pairedMsa"],
                       "templates": [] if args.no_templates else cached["templates"]}
            meta["msa"] = "cached_variant_query_replaced"
            log(f"WARNING: reusing MSA of cached sequence for a variant:\n  {msg}")
        else:
            log(f"SEQUENCE MISMATCH vs cache {cache}:\n  {msg}\n"
                "Refusing to reuse MSA. Compute a new MSA, or use --reuse-msa-for-variant "
                "(same length only).")
            return 2
    else:
        if args.msa_free:
            protein = {"sequence": wanted, "unpairedMsa": "", "pairedMsa": "", "templates": []}
            meta["msa"] = "none_single_sequence"
            log("WARNING: no cached MSA; writing single-sequence (MSA-free) input")
        else:
            log(f"NO CACHED MSA for {acc} in {cache_dirs}. An MSA must be computed "
                "(see docs/COFOLD.md: ColabFold/MMseqs2 server or local jackhmmer).")
            Path(args.out).with_suffix(".NEEDS_MSA").write_text(f"{acc}\n{wanted}\n")
            return 3

    # 3) assemble
    ids = chain_ids()
    seqs = []
    prot_ids = [next(ids) for _ in range(args.copies)]
    seqs.append({"protein": {"id": prot_ids if len(prot_ids) > 1 else prot_ids[0], **protein}})
    lig_tags = []
    for spec in args.ligand:
        lig = parse_ligand(spec)
        seqs.append({"ligand": {"id": next(ids), **lig}})
        lig_tags.append("-".join(lig["ccdCodes"]) if "ccdCodes" in lig else "smi")
    name = args.name or "_".join([acc.lower()] + [t.lower() for t in lig_tags])
    if meta["msa"] == "cached_variant_query_replaced" and "variantmsa" not in name:
        name += "_variantmsa"
    job = {"dialect": "alphafold3", "version": 2, "name": name, "modelSeeds": seeds,
           "sequences": seqs, "bondedAtomPairs": None, "userCCD": None}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(job, f)
    meta.update({"name": name, "length": len(wanted), "copies": args.copies,
                 "ligands": args.ligand, "seeds": seeds, "templates": len(protein["templates"])})
    Path(args.out).with_suffix(".meta.json").write_text(json.dumps(meta, indent=1))
    log(f"wrote {args.out} (name={name}, seeds={seeds})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
