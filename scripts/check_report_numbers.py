#!/usr/bin/env python
"""Check that every notable number in a markdown report has a producer.

A number is 'notable' if it has >=2 significant digits, a decimal point, or a percent sign (code blocks, dates,
times, section numbers, DOIs, accession-like tokens, years, URLs/link targets are skipped). It matches if it equals
(within rounding to the shown precision) a numeric leaf of results/numbers.json, or a number in a JSON/CSV file linked
in the same paragraph/table (+-3 lines), or it is listed in config/report_numbers_allowlist.txt.
Percentages match v*100 (or v itself). Comparison is on absolute value.
Usage: check_report_numbers.py [report.md] [--root P]
"""
import argparse, csv, json, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import provlib as L

NUM = re.compile(r"(?<![\w.\-/])[-−+]?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?\s?%?|(?<![\w.\-/])[-−+]?\.\d+%?")
SKIP_CTX = [
    re.compile(r"\d{4}-\d{2}-\d{2}(?:[T ]\d{2}:\d{2}(?::\d{2})?)?"),   # dates
    re.compile(r"\b\d{1,2}:\d{2}(?::\d{2})?\b"),                          # times
    re.compile(r"\b\d{8}T\d{6}Z?\b"),                                     # snapshot ids
    re.compile(r"\b10\.\d{4,9}/\S+"),                                     # DOIs
    re.compile(r"[§§]\s?\d+(?:\.\d+)*"),                                  # section refs
    re.compile(r"\b(?:PMID|PMC|PDB|CHEMBL|ChEMBL|UniProt|EC|GO|ECO|HOG|OMA|arXiv)[:\s_-]*[A-Za-z]*\d[\w.]*", re.I),
    re.compile(r"\b[A-Za-z]+[_-]?\d+[A-Za-z0-9_.-]*\b"),                  # identifiers with letters (D7, G2, AF3, p53, v1)
    re.compile(r"\b\d+[A-Za-z_][\w.-]*\b"),                               # 20261030T..., 1abc, 3D, 5x
    re.compile(r"\b[0-9a-f]{7,}\b"),                                      # hashes
    re.compile(r"https?://\S+"),
    re.compile(r"\]\([^)]*\)"),                                           # link targets
    re.compile(r"\(\d{4}\)|\b(?:19|20)\d{2}\b(?!\.\d)"),                  # years
]
LINK = re.compile(r"\[[^\]]*\]\(([^)#\s]+)")


def flatten(o, out):
    if isinstance(o, bool):
        return
    if isinstance(o, (int, float)):
        out.append(float(o))
    elif isinstance(o, dict):
        for v in o.values():
            flatten(v, out)
    elif isinstance(o, list):
        for v in o:
            flatten(v, out)


def file_numbers(p):
    out = []
    try:
        if p.suffix == ".json":
            flatten(json.load(open(p)), out)
        elif p.suffix in (".csv", ".tsv"):
            with open(p, newline="") as f:
                for row in csv.reader(f, delimiter="\t" if p.suffix == ".tsv" else ","):
                    for c in row:
                        try:
                            out.append(float(c))
                        except ValueError:
                            pass
    except (OSError, ValueError, csv.Error):
        pass
    return out


def parse(tok):
    t = tok.strip().replace("−", "-").replace(",", "")
    pct = t.endswith("%")
    t = t.rstrip("%").strip().lstrip("+-")
    mant = re.split(r"[eE]", t)[0]
    dec = len(mant.split(".")[1]) if "." in mant else 0
    return float(t), dec, pct, mant


def notable(tok):
    v, dec, pct, mant = parse(tok)
    digits = mant.replace(".", "").lstrip("0")
    digits = digits.rstrip("0") if "." not in mant else digits
    return pct or dec > 0 or len(digits) >= 2


def matches(tok, pool):
    x, dec, pct, _ = parse(tok)
    tol = 0.5 * 10 ** (-dec) * (1 + 1e-9)
    exps = [(x, tol)]
    if pct:
        exps.append((x / 100, tol / 100))
    for target, t in exps:
        for v in pool:
            if abs(abs(v) - target) <= t:
                return True
    return False


def extract(lines):
    """Yield (lineno, token) for notable numbers outside code, dates, ids, etc."""
    for i, ln in enumerate(lines, 1):
        if re.match(r"^\s*#+\s", ln) or re.match(r"^\s*(\d+\.)+\s", ln):
            ln = re.sub(r"^\s*(#+\s*)?(\d+\.)+\d*\s", " ", ln)   # heading / list numbering
        ln = re.sub(r"^\s*\|?\s*:?-{3,}.*$", "", ln)
        for rx in SKIP_CTX:
            ln = rx.sub(lambda m: " " * len(m.group(0)), ln)
        for m in NUM.finditer(ln):
            tok = m.group(0).strip()
            try:
                if notable(tok):
                    yield i, tok
            except ValueError:
                continue


def load_allow(path):
    toks = []
    if Path(path).is_file():
        for ln in Path(path).read_text().splitlines():
            ln = ln.split("#")[0].strip()
            if ln:
                toks.append(ln)
    return toks


def check(report, root, numbers=None, allow=None):
    root, report = Path(root), Path(report)
    raw = report.read_text(errors="replace")
    lines = L.strip_code(raw)
    full = raw.splitlines()
    numbers = Path(numbers) if numbers else root / "results/numbers.json"
    pool = []
    if numbers.is_file():
        flatten(json.load(open(numbers)), pool)
    allow = load_allow(allow if allow else root / "config/report_numbers_allowlist.txt")
    allow_vals = []
    for a in allow:
        try:
            allow_vals.append(parse(a)[0])
        except ValueError:
            pass
    cache, unmatched, n = {}, [], 0
    for i, tok in extract(lines):
        n += 1
        if matches(tok, pool) or any(abs(parse(tok)[0] - v) < 1e-12 for v in allow_vals):
            continue
        lo, hi = max(0, i - 4), min(len(full), i + 3)
        # extend to paragraph edges not needed; use +-3 lines window
        ok = False
        for tgt in LINK.findall(" ".join(full[lo:hi])):
            for base in (report.parent, root):
                p = (base / tgt).resolve()
                if p.is_file() and p.suffix in (".json", ".csv", ".tsv"):
                    if p not in cache:
                        cache[p] = file_numbers(p)
                    if matches(tok, cache[p]):
                        ok = True
                    break
            if ok:
                break
        if not ok:
            unmatched.append((i, tok, full[i - 1].strip()[:100]))
    return n, unmatched, len(pool)


def main(argv=None):
    root0 = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser()
    ap.add_argument("report", nargs="?", default=str(root0 / "record/REPORT.md"))
    ap.add_argument("--root", default=str(root0))
    ap.add_argument("--numbers")
    ap.add_argument("--allowlist")
    a = ap.parse_args(argv)
    if not Path(a.report).is_file():
        print(f"FAIL: report {a.report} not found")
        return 2
    n, un, npool = check(a.report, a.root, a.numbers, a.allowlist)
    print(f"{a.report}: {n} notable numbers, {npool} values in numbers.json, {len(un)} unmatched")
    for i, tok, ctx in un:
        print(f"  line {i}: {tok}    | {ctx}")
    print("RESULT:", "FAIL" if un else "PASS")
    return 1 if un else 0


if __name__ == "__main__":
    sys.exit(main())
