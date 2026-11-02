#!/usr/bin/env python
"""Verify prereg freeze, sidecar prereg hashes, and deviation-ledger references.

Usage: check_prereg.py [--root P] [--docs-glob 'results/**/*.md ...']
Exit 0 only if: prereg hash == frozen, no sidecar mismatches/missing, every D-id in results docs is in DESIGN ledger.
"""
import argparse, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import provlib as L

SIDE_EXT = re.compile(r"\.(parquet|csv|tsv|png|pdf|svg|npz|npy|pkl|txt|md|html)\.json$")


def check(root, frozen=L.FROZEN_PREREG_SHA256, doc_roots=("results",)):
    root = Path(root)
    res = {"fail": []}
    actual = L.sha256_file(root / "config/prereg.yaml")
    res["prereg_ok"] = actual == frozen
    res["prereg_actual"] = actual
    if not res["prereg_ok"]:
        res["fail"].append(f"config/prereg.yaml sha256 {actual} != frozen {frozen}")
    stats = {"n": 0, "ok": 0, "differ": [], "missing": []}
    for p in L.iter_json(root / "results", root / "cases"):
        d = L.load_json(p)
        if d is None:
            continue
        stats["n"] += 1
        h = L.prereg_of(d)
        rel = str(p.relative_to(root))
        if h is None:
            stats["missing"].append(rel)
        elif h != frozen:
            stats["differ"].append((rel, h))
        else:
            stats["ok"] += 1
    res["sidecars"] = stats
    if stats["differ"] or stats["missing"]:
        res["fail"].append(f"{len(stats['differ'])} json outputs with differing and {len(stats['missing'])} with missing prereg hash")
    ledger = L.ledger_ids(root / "record/DESIGN.md")
    refs = {}
    for dr in doc_roots:
        for md in sorted((root / dr).rglob("*.md")):
            for i, ln in enumerate(L.strip_code(md.read_text(errors="replace")), 1):
                for m in re.finditer(r"(?<![A-Za-z0-9_])D(\d+)(?![A-Za-z0-9_])", ln):
                    refs.setdefault("D" + m.group(1), []).append(f"{md.relative_to(root)}:{i}")
    res["ledger"] = sorted(ledger)
    res["refs"] = refs
    res["unknown_ids"] = {k: v for k, v in refs.items() if k not in ledger}
    if res["unknown_ids"]:
        res["fail"].append(f"deviation ids referenced but not in ledger: {sorted(res['unknown_ids'])}")
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--max-list", type=int, default=15)
    a = ap.parse_args(argv)
    r = check(a.root)
    s = r["sidecars"]
    print("| check | result | detail |\n|---|---|---|")
    print(f"| prereg.yaml hash | {'PASS' if r['prereg_ok'] else 'FAIL'} | {r['prereg_actual'][:12]}... |")
    print(f"| json outputs scanned | {s['n']} | ok={s['ok']} differ={len(s['differ'])} missing={len(s['missing'])} |")
    print(f"| ledger ids | {len(r['ledger'])} | {', '.join(r['ledger'])} |")
    print(f"| referenced ids | {len(r['refs'])} | unknown={sorted(r['unknown_ids'])} |")
    for title, items in (("DIFFERING prereg hash", [f"{a_} ({b[:12]})" for a_, b in s["differ"]]),
                         ("MISSING prereg hash", s["missing"])):
        if items:
            print(f"\n{title} ({len(items)}), first {a.max_list}:")
            for x in items[:a.max_list]:
                print("  ", x)
    for k, v in r["unknown_ids"].items():
        print(f"\nunknown ledger id {k}: {', '.join(v[:5])}")
    print("\nRESULT:", "FAIL" if r["fail"] else "PASS")
    for f in r["fail"]:
        print(" -", f)
    return 1 if r["fail"] else 0


if __name__ == "__main__":
    sys.exit(main())
