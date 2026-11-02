#!/usr/bin/env python
"""Provenance checks for results/: sidecar commits exist in git, dirty trees, uncommitted modifications, checksum manifests.

Exit non-zero on any failure. Missing bulk files are reported as 'bulk-missing (node-local)', not failures.
"""
import argparse, re, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import provlib as L

MANIFEST = re.compile(r"(BULK_MANIFEST|checksums?)[^/]*$", re.I)


def git(root, *args):
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True)


def commit_exists(root, sha):
    return git(root, "cat-file", "-e", sha + "^{commit}").returncode == 0


def verify_manifest(mpath, results_root):
    """Return dict(ok, bad[list], missing[list], n). Paths resolved vs manifest dir, then its ancestors up to results."""
    ok, bad, missing = 0, [], []
    for ln in Path(mpath).read_text(errors="replace").splitlines():
        m = re.match(r"^([0-9a-fA-F]{64})\s+\*?(.+)$", ln.strip())
        if not m:
            continue
        h, rel = m.group(1).lower(), m.group(2).strip()
        cand, base = None, Path(mpath).parent
        while True:
            if (base / rel).is_file():
                cand = base / rel
                break
            if base == results_root or base == base.parent:
                break
            base = base.parent
        if cand is None:
            missing.append(rel)
        elif L.sha256_file(cand) != h:
            bad.append(rel)
        else:
            ok += 1
    return {"ok": ok, "bad": bad, "missing": missing}


def check(root):
    root = Path(root)
    res = {"fail": []}
    commits, dirty, nocommit, notfound = {}, [], [], []
    for p in L.iter_json(root / "results", root / "cases"):
        d = L.load_json(p)
        if not isinstance(d, dict) or (L.prereg_of(d) is None and not L.first_key(d, L.COMMIT_KEYS)[0]):
            continue
        rel = str(p.relative_to(root))
        c = L.commit_of(d)
        if not c:
            nocommit.append(rel)
            continue
        is_dirty = bool(re.search(r"[-+.]dirty$", c))
        sha = re.sub(r"[-+.]dirty$", "", c)
        if is_dirty:
            dirty.append((rel, sha[:10]))
        if sha not in commits:
            commits[sha] = commit_exists(root, sha) if re.fullmatch(r"[0-9a-f]{7,40}", sha) else False
        if not commits[sha]:
            notfound.append((rel, sha))
    res.update(commits=commits, dirty=dirty, nocommit=nocommit, notfound=notfound)
    if notfound:
        res["fail"].append(f"{len(notfound)} sidecars cite a commit not in git history")
    if dirty:
        res["fail"].append(f"{len(dirty)} sidecars were produced from a dirty tree")
    # tracked results files modified vs HEAD
    st = git(root, "status", "--porcelain", "--untracked-files=no", "--", "results", "cases").stdout.splitlines()
    res["modified"] = [l[3:] for l in st if l.strip()]
    if res["modified"]:
        res["fail"].append(f"{len(res['modified'])} tracked results files differ from HEAD")
    un = git(root, "status", "--porcelain", "--", "results", "cases").stdout.splitlines()
    res["untracked"] = sorted({l[3:] for l in un if l.startswith("??")})
    # manifests
    mans = []
    for p in sorted((root / "results").rglob("*")):
        if p.is_file() and MANIFEST.search(p.name) and not p.name.endswith(".json"):
            v = verify_manifest(p, root / "results")
            v["path"] = str(p.relative_to(root))
            mans.append(v)
            if v["bad"]:
                res["fail"].append(f"{v['path']}: {len(v['bad'])} files fail checksum")
    res["manifests"] = mans
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--max-list", type=int, default=12)
    a = ap.parse_args(argv)
    r = check(a.root)
    n = a.max_list
    print("| check | count | detail |\n|---|---|---|")
    print(f"| distinct sidecar commits | {len(r['commits'])} | missing from git: {sum(not v for v in r['commits'].values())} |")
    print(f"| sidecars w/o commit field | {len(r['nocommit'])} | |")
    print(f"| sidecars citing unknown commit | {len(r['notfound'])} | |")
    print(f"| dirty-tree sidecars | {len(r['dirty'])} | commits: {sorted({c for _, c in r['dirty']})} |")
    print(f"| tracked results modified vs HEAD | {len(r['modified'])} | |")
    print(f"| untracked results paths | {len(r['untracked'])} | (info) |")
    for m in r["manifests"]:
        status = "FAIL" if m["bad"] else "ok"
        extra = f"; bulk-missing (node-local): {len(m['missing'])}" if m["missing"] else ""
        print(f"| manifest {m['path']} | {status} | verified={m['ok']} bad={len(m['bad'])}{extra} |")
    for title, items in (("UNKNOWN COMMIT", [f"{a_} {b}" for a_, b in r["notfound"]]),
                         ("DIRTY", [f"{a_} {b}" for a_, b in r["dirty"]]),
                         ("MODIFIED vs HEAD", r["modified"]),
                         ("NO COMMIT FIELD", r["nocommit"]),
                         ("UNTRACKED (info)", r["untracked"])):
        if items:
            print(f"\n{title} ({len(items)}), first {n}:")
            for x in items[:n]:
                print("  ", x)
    for m in r["manifests"]:
        if m["bad"]:
            print(f"\nCHECKSUM MISMATCH in {m['path']} ({len(m['bad'])}), first {n}:")
            for x in m["bad"][:n]:
                print("  ", x)
    print("\nRESULT:", "FAIL" if r["fail"] else "PASS")
    for f in r["fail"]:
        print(" -", f)
    return 1 if r["fail"] else 0


if __name__ == "__main__":
    sys.exit(main())
