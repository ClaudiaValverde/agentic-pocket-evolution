"""Run the research browser:  uv run --project . python -m webapp --port 8765   (from the repo root)."""
from __future__ import annotations

import argparse


def build_index(argv):
    """`python -m webapp build-index [--tiers 1,2,3] [--snapshot S] ...`: (re)generate the protein index by running
    scripts/webapp_build_index.py (tiers 1-2 -> results/webapp_index, tier 3 -> /tmp/pe_big/webapp_index, ~2 min)."""
    import subprocess
    import sys

    from .config import DEFAULT_SNAPSHOT, REPO
    argv = list(argv)
    if "--tiers" not in argv:
        argv += ["--tiers", "1,2,3"]
    if "--snapshot" not in argv:
        argv += ["--snapshot", DEFAULT_SNAPSHOT]
    raise SystemExit(subprocess.call([sys.executable, "-I", str(REPO / "scripts" / "webapp_build_index.py")] + argv))


def main(argv=None):
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "build-index":
        return build_index(argv[1:])
    if argv and argv[0] == "build-og-index":
        import subprocess
        from .config import REPO
        raise SystemExit(subprocess.call([sys.executable, "-I", str(REPO / "scripts" / "webapp_build_og_index.py")] + argv[1:]))
    ap = argparse.ArgumentParser(prog="python -m webapp", description="pocket-evolution research browser")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--snapshot", default=None, help="atlas snapshot id (default: config.DEFAULT_SNAPSHOT)")
    ap.add_argument("--index-dirs", default=None, help="colon-separated index roots (default: results/webapp_index/<snap>:/tmp/pe_big/webapp_index/<snap>)")
    a = ap.parse_args(argv)

    import uvicorn

    from .app import create_app
    from .config import Config

    kw = {}
    if a.snapshot:
        kw["snapshot"] = a.snapshot
    if a.index_dirs:
        from pathlib import Path
        kw["index_dirs"] = [Path(p) for p in a.index_dirs.split(":")]
    app = create_app(Config(**kw))
    tiers = sorted(app.state.atlas.tier_dirs())
    print(f"pocket-evolution browser: http://{a.host}:{a.port}/  (index tiers: {tiers}; "
          f"orthology versions: {list(app.state.ogs.versions())})", flush=True)
    if tiers != [1, 2, 3]:
        print("  note: protein index incomplete (tier 3 lives in node-local /tmp). Regenerate: "
              "uv run --project . python -m webapp build-index", flush=True)
    uvicorn.run(app, host=a.host, port=a.port, log_level="info")


if __name__ == "__main__":
    main()
