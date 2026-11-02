"""Shared helpers for check_prereg / check_provenance / check_report_numbers (stdlib only)."""
import hashlib, json, re
from pathlib import Path

FROZEN_PREREG_SHA256 = "52f9461890984c210b730af86c40de221307a668e7d4b38098b3316059690d83"
PREREG_KEYS = ("prereg_sha256", "prereg_hash", "prereg_sha", "prereg")
COMMIT_KEYS = ("git_commit", "code_commit", "git_sha", "commit", "git_hash")
SEED_KEYS = ("seed", "seeds")
SNAP_KEYS = ("snapshot_id", "snapshot")
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def sha256_file(p, chunk=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def load_json(p):
    try:
        with open(p) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def iter_json(*roots):
    for r in roots:
        r = Path(r)
        if r.is_dir():
            yield from sorted(r.rglob("*.json"))


def first_key(d, keys, nested=True):
    """Value of first key found at top level (then one level down). Returns (found, value)."""
    if not isinstance(d, dict):
        return False, None
    for k in keys:
        if k in d:
            return True, d[k]
    if nested:
        for v in d.values():
            if isinstance(v, dict):
                for k in keys:
                    if k in v:
                        return True, v[k]
    return False, None


def prereg_of(d):
    """Return prereg hash string or None. Handles str, {'sha256':..} dicts, and 'sha256:'-prefixed."""
    ok, v = first_key(d, PREREG_KEYS)
    if not ok:
        return None
    if isinstance(v, dict):
        v = v.get("sha256") or v.get("hash")
    if not isinstance(v, str):
        return None
    return v.lower().replace("sha256:", "").strip()


def commit_of(d):
    ok, v = first_key(d, COMMIT_KEYS)
    return v if ok and isinstance(v, str) else None


def ledger_ids(design_md):
    """Ids (D7, C1, ...) defined as first cell of a table row in the deviation ledger section."""
    txt = Path(design_md).read_text()
    m = re.search(r"^#+\s*Deviation ledger.*$", txt, re.M | re.I)
    if not m:
        return set()
    return set(re.findall(r"^\|\s*([A-Z]\d+)\s*\|", txt[m.end():], re.M))


def strip_code(text):
    """Blank out fenced code blocks (keeping line count) and inline code."""
    out, fence = [], False
    for ln in text.splitlines():
        if re.match(r"^\s*(```|~~~)", ln):
            fence = not fence
            out.append("")
            continue
        out.append("" if fence else re.sub(r"`[^`]*`", " ", ln))
    return out
