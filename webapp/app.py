"""FastAPI app: protein / orthogroup / case views (server-rendered HTML + small JSON API)."""
from __future__ import annotations

import html
import json
import mimetypes
import re
import threading
from urllib.parse import quote, urlencode

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles

from .atlas import Atlas, fragment_offset
from .cases import Cases
from .config import REPO, WEBAPP, Config
from .orthogroups import STATE_COLORS, STATE_NAMES, STATE_ORDER, Orthogroups, layout, leaves, matrix, strongest

e = lambda x: html.escape("" if x is None else str(x))  # noqa: E731


def fmt(x, nd=2):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def page(title: str, body: str, scripts: str = "", vsel: str = "") -> HTMLResponse:
    return HTMLResponse(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>{e(title)} · pocket-evolution</title>
<link rel="stylesheet" href="/static/style.css"></head><body>
<header><a class="brand" href="/">pocket-evolution browser</a>
<form action="/search" method="get"><input name="q" placeholder="UniProt accession or gene" aria-label="search">
<button>Go</button></form><nav><a href="/og">orthogroups</a> <a href="/sog">SOGs</a> <a href="/case">cases</a> <a href="/validation">validation</a> {vsel}</nav></header>
<main>{body}</main>{scripts}</body></html>""")


def table(rows, cols, cls="", row_attrs=None):
    head = "".join(f"<th>{e(c[1])}</th>" for c in cols)
    body = []
    for i, r in enumerate(rows):
        attrs = row_attrs(r, i) if row_attrs else ""
        body.append(f"<tr {attrs}>" + "".join(f"<td>{c[2](r) if len(c) > 2 else e(fmt(r.get(c[0])))}</td>"
                                             for c in cols) + "</tr>")
    return f'<table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>'


def jclean(o):
    """NaN/Inf -> None so that FastAPI's strict JSON encoder accepts parquet rows."""
    if isinstance(o, float):
        return o if o == o and o not in (float("inf"), float("-inf")) else None
    if isinstance(o, dict):
        return {k: jclean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jclean(v) for v in o]
    return o


def badge(status):
    s = (status or "").lower()
    return f'<span class="badge b-{e(s)}">{e(status)}</span>'


def create_app(cfg: Config | None = None) -> FastAPI:
    cfg = cfg or Config()
    atlas, ogs, cases = Atlas(cfg), Orthogroups(cfg.sog_dirs, cfg.index_dirs), Cases(cfg.case_dirs)
    app = FastAPI(title="pocket-evolution browser", docs_url="/api/docs", redoc_url=None)
    threading.Thread(target=atlas.genes, daemon=True).start()  # warm the gene table (about 3 s)
    app.state.atlas, app.state.cfg, app.state.ogs, app.state.cases = atlas, cfg, ogs, cases
    app.mount("/static", StaticFiles(directory=WEBAPP / "static"), name="static")

    # ------------------------------------------------------------------ home / search
    @app.get("/", response_class=HTMLResponse)
    def home(request: Request):
        ver = getv(request)
        st = atlas.status()
        tiers = ", ".join(f"tier {t}" for t in sorted(st["index_tiers"])) or "none (run: python -m webapp build-index)"
        gene = "gene names from " + e(st["gene_table"]) if st["gene_table"] else \
            "accession only (no gene table in data/annotations/ yet)"
        return page("home", f"""<h1>Pocket-evolution research browser</h1>
<p>Atlas snapshot <code>{e(cfg.snapshot)}</code>; indexed: {tiers}. Search: {gene}.</p>
<p>Examples: <a href="/protein/P24941">P24941</a> (CDK2, tier 1), <a href="/protein/P04275">P04275</a>
(VWF, 9 fragments), <a href="/protein/A0A633RIM1">A0A633RIM1</a> (tier 3, Foldcomp).</p>
<h2>Orthogroups{(" · " + e(ver.id)) if ver else ""}</h2>{og_list_html(request, ver, limit=15)}<h2>Cases</h2>{case_list_html()}""", vsel=vsel_html(ver))

    @app.get("/search")
    def search(q: str = ""):
        r = atlas.search(q)
        if r["by"] == "accession" or len(r["accessions"]) == 1:
            return RedirectResponse(f"/protein/{quote(r['accessions'][0])}", status_code=303)
        if r["hits"]:
            rows = table(r["hits"], [("accession", "accession", lambda h: f'<a href="/protein/{e(h["accession"])}">{e(h["accession"])}</a>'),
                                     ("gene", "gene"), ("protein_name", "protein"), ("match", "match")])
            return page("search", f"<h1>{len(r['hits'])} hits for “{e(q)}”</h1>{rows}")
        hint = "" if r.get("gene_search") else " Gene-name search needs data/annotations/uniprot_proteins.parquet."
        return page("search", f"<h1>No atlas entry for “{e(q)}”</h1><p>Not an accession in the indexed tiers "
                              f"{e(sorted(atlas.tier_dirs()))}, and no gene or protein name matches.{hint}</p>")

    @app.get("/api/status")
    def api_status():
        return {**atlas.status(), "og_sources": [str(d) for d, _ in ogs.sources()],
                "case_dirs": [str(d) for d in cfg.case_dirs]}

    @app.get("/api/search")
    def api_search(q: str):
        return atlas.search(q)

    # ------------------------------------------------------------------ protein
    @app.get("/api/protein/{acc}")
    def api_protein(acc: str):
        structs = atlas.lookup(acc)
        if not structs:
            raise HTTPException(404, f"{acc} not in the atlas index")
        _, a2g = atlas.genes()
        g = a2g.get(acc.upper(), (None, None))
        return jclean({"accession": acc.upper(), "gene": g[0], "protein_name": g[1], "structures": structs,
                       "pockets": atlas.pockets(acc), "orthogroups": ogs.by_accession(acc)})

    @app.get("/api/structure/{acc}")
    def api_structure(acc: str, tier: int | None = None, fragment: int | None = None):
        try:
            txt, s = atlas.structure(acc, tier, fragment)
        except KeyError:
            raise HTTPException(404, f"{acc} (tier={tier}, fragment={fragment}) not in the atlas index")
        except FileNotFoundError as ex:
            raise HTTPException(503, f"structure source not readable: {ex}")
        return Response(txt, media_type="chemical/x-mmcif", headers={
            "X-Atlas-Tier": str(s["tier"]), "X-Atlas-Fragment": str(s["fragment"]),
            "X-Atlas-Source": f"{s['blob_kind']}:{s['blob_file']}",
            "Content-Disposition": f'inline; filename="AF-{s["accession"]}-F{s["fragment"]}-model_v{s["model_version"]}.cif"'})

    @app.get("/protein/{acc}", response_class=HTMLResponse)
    def protein(acc: str, tier: int | None = None, fragment: int | None = None, pocket: str | None = None):
        acc = acc.strip().upper()
        mm = re.search(r"-F(\d+)-", pocket or "")
        if mm and fragment is None:
            fragment = int(mm.group(1))
        structs = atlas.lookup(acc)
        if not structs:
            raise HTTPException(404, f"{acc} not in the atlas index (tiers {sorted(atlas.tier_dirs())})")
        sel = [s for s in structs if (tier is None or s["tier"] == tier) and (fragment is None or s["fragment"] == fragment)]
        sel = sel[0] if sel else structs[0]
        pockets = atlas.pockets(acc)
        _, a2g = atlas.genes()
        gene, pname = a2g.get(acc, (None, None))
        st_rows = table(structs, [
            ("tier", "tier"), ("fragment", "fragment"),
            ("model", "model", lambda s: f'<a href="?tier={s["tier"]}&fragment={s["fragment"]}">AF-{e(s["accession"])}-F{s["fragment"]}-model_v{s["model_version"]}</a>'
             + (" (shown)" if s is sel else "")),
            ("organism", "organism"), ("taxid", "taxid"), ("n_res", "residues"),
            ("plddt_mean", "pLDDT mean", lambda s: fmt(s["plddt_mean"], 1)),
            ("n_pockets_kept", "pockets kept"), ("status", "status"),
            ("src", "source", lambda s: e(s["blob_file"]))])
        sel_p = [p for p in pockets if p["tier"] == sel["tier"] and p["fragment"] == sel["fragment"]]
        other = len(pockets) - len(sel_p)

        def attrs(p, i):
            return (f'class="pocket" data-pid="{e(p["pocket_id"])}" data-residues="{e(p["residues"])}" '
                    f'data-center="{p["center_x"]:.3f},{p["center_y"]:.3f},{p["center_z"]:.3f}"')
        pk = table(sel_p, [
            ("pocket_id", "pocket", lambda p: f'<code>{e(p["pocket_id"].split(":")[-1])}</code>'),
            ("burrow_rank", "BURROW rank"), ("burrow_score", "score", lambda p: fmt(p["burrow_score"], 3)),
            ("volume_A3", "volume Å³", lambda p: fmt(p["volume_A3"], 0)),
            ("site_volume_A3", "site vol Å³", lambda p: fmt(p["site_volume_A3"], 0)),
            ("plddt_mean", "pLDDT mean", lambda p: fmt(p["plddt_mean"], 1)),
            ("plddt_min", "pLDDT min", lambda p: fmt(p["plddt_min"], 1)),
            ("n", "n res", lambda p: str(len(p["residues"].split()))),
            ("residues", "residues (UniProt numbering)", lambda p: f'<span class="res">{e(p["uniprot_positions"])}</span>')],
            cls="pockets", row_attrs=attrs)
        off = fragment_offset(sel["fragment"])
        note = f" Fragment {sel['fragment']}: model residue r = UniProt position r + {off}." if off else ""
        other_note = f"<p class='muted'>{other} more pockets in other fragments/tiers (select the model above).</p>" if other else ""
        ogl = ogs.by_accession(acc)
        if ogl:
            def ogrow(r):
                cells = " ".join(
                    f'<a class="st s{e(x["state"])}" href="/og/{e(r["og_id"])}?v={quote(r["version"])}&sog={quote(x["sog_id"])}" '
                    f'title="{e(x["sog_id"])} {e(x.get("event") or "")}">{e(x["sog_id"].split(":")[-1])}:{e(x["state"])}</a>'
                    for x in r["sogs"][:60]) + (f" … ({len(r['sogs'])} SOGs)" if len(r["sogs"]) > 60 else "")
                return cells
            og_html = "<h2>Orthogroups and site orthogroups</h2>" + table(ogl, [
                ("version", "version"),
                ("og_id", "OG", lambda r: f'<a href="/og/{e(r["og_id"])}?v={quote(r["version"])}">{e(r["og_id"])}</a>'),
                ("label", "label"), ("n_sogs", "SOGs in OG"), ("sogs", "this protein's SOGs (state)", ogrow)])
        else:
            og_html = ""
        struct_url = f"/api/structure/{quote(acc)}?tier={sel['tier']}&fragment={sel['fragment']}"
        title = f"{acc}" + (f" · {gene}" if gene else "")
        body = f"""<h1>{e(title)}</h1><p>{e(pname or '')} <a href="https://www.uniprot.org/uniprotkb/{e(acc)}">UniProt</a>
· <a href="https://alphafold.ebi.ac.uk/entry/{e(acc)}">AFDB</a> · <a href="/api/protein/{e(acc)}">JSON</a></p>
<h2>Atlas structures</h2>{st_rows}
<div class="split"><div class="left"><h2>Pockets ({len(sel_p)}) in AF-{e(acc)}-F{sel['fragment']}</h2>
<p class="muted">Atlas filters: score ≥ 0.03, lining pLDDT mean ≥ 70. Click a row to highlight its lining residues.{e(note)}</p>
{pk}{other_note}</div>
<div class="right"><div class="vbar"><label><input type="checkbox" id="plddt" checked> colour by pLDDT</label>
<label><input type="checkbox" id="surf"> pocket surface</label> <button id="reset">reset view</button>
<a href="{struct_url}" download>mmCIF</a></div>
<div id="viewer" class="viewer" data-url="{struct_url}"></div>
<div class="legend"><span style="background:#0053D6"></span>&gt;90 <span style="background:#65CBF3"></span>70–90
<span style="background:#FFDB13"></span>50–70 <span style="background:#FF7D45"></span>&lt;50 · pocket residues: orange sticks</div>
<div id="vmsg" class="muted"></div></div></div>{og_html}"""
        return page(acc, body, '<script src="/static/vendor/3Dmol-min.js"></script><script src="/static/viewer.js"></script>')

    # ------------------------------------------------------------------ orthogroups
    def getv(request: Request):
        """Selected orthology version: ?v= wins, then the pe_v cookie (set by the selector), then the default."""
        return ogs.pick(request.query_params.get("v") or request.cookies.get("pe_v"))

    def vsel_html(ver):
        vs = ogs.versions()
        if not vs:
            return ""
        opts = "".join(f'<option value="{e(k)}"{" selected" if ver is not None and k == ver.id else ""}>'
                       f'{e(k)}{" (toy)" if x.toy else ""}</option>' for k, x in vs.items())
        return ('<label class="vsel">orthology <select id="vsel" aria-label="orthology version">' + opts + '</select></label>'
                '<script>document.getElementById("vsel").addEventListener("change",function(){'
                'document.cookie="pe_v="+encodeURIComponent(this.value)+"; path=/; max-age=31536000; SameSite=Lax";'
                'var u=new URL(location.href);u.searchParams.set("v",this.value);u.searchParams.delete("page");'
                'u.searchParams.delete("sp");location.href=u.toString();});</script>')

    def P(title, body, request, scripts=""):
        return page(title, body, scripts, vsel_html(getv(request)))

    def vq(ver):
        return f"v={quote(ver.id)}" if ver is not None else ""

    def pager(request: Request, key, total, size, cur):
        n = max(1, -(-total // size))
        if n <= 1:
            return f"<p class='muted'>{total} rows</p>"

        def link(p, txt):
            qp = dict(request.query_params)
            qp[key] = str(p)
            return f'<a href="?{urlencode(qp)}">{txt}</a>'
        parts = [f"{total} rows · page {cur} / {n}"]
        if cur > 1:
            parts = [link(1, "«"), link(cur - 1, "‹ prev")] + parts
        if cur < n:
            parts += [link(cur + 1, "next ›"), link(n, "»")]
        return f'<p class="pager">{" ".join(parts)}</p>'

    def qint(request, key, default, lo=1, hi=10**9):
        try:
            return max(lo, min(hi, int(request.query_params.get(key, default))))
        except ValueError:
            return default

    def og_filter_form(request, extra=""):
        qp = request.query_params
        return (f'<form class="filters" method="get"><input name="q" value="{e(qp.get("q", ""))}" placeholder="OG id, label or reference accession">'
                f'<input name="flags" value="{e(qp.get("flags", ""))}" placeholder="flag contains" size="12">'
                f'<select name="half"><option value="">tune+test</option>'
                + "".join(f'<option{" selected" if qp.get("half") == h else ""}>{h}</option>' for h in ("tune", "test"))
                + f'</select>{extra}<input type="hidden" name="v" value="{e(qp.get("v") or "")}"><button>filter</button></form>')

    def og_list_html(request: Request, ver, limit=None):
        if ver is None:
            return "<p class='muted'>No orthogroup data yet (results/sogs/&lt;snapshot&gt;/&lt;version&gt;/web/).</p>"
        size = limit or qint(request, "size", 50, 5, 500)
        pg = qint(request, "page", 1)
        qp = request.query_params
        rows, total = ver.og_page(q=qp.get("q", ""), flags=qp.get("flags", ""), half=qp.get("half", ""),
                                  sort=qp.get("sort", "og_id"), desc=qp.get("desc") == "1", page=pg, size=size)
        t = table(rows, [("og_id", "OG", lambda r: f'<a href="/og/{e(r["og_id"])}?{vq(ver)}">{e(r["og_id"])}</a>'),
                         ("label", "label"), ("n_species", "species"), ("n_members", "members"), ("n_sogs", "SOGs"),
                         ("ref_accession", "reference", lambda r: acc_link(r.get("ref_accession"))),
                         ("level", "root level"), ("flags", "flags")])
        return t + (pager(request, "page", total, size, pg) if not limit else
                    f"<p class='muted'>first {limit} of {total}; <a href='/og?{vq(ver)}'>all orthogroups</a></p>")

    def acc_link(a, tier=None, pocket=None):
        if not a:
            return ""
        a = str(a)
        if a.startswith("TOY") or (tier is not None and tier not in (1, 2, 3)):
            return e(a)
        return f'<a href="/protein/{e(a)}{("?pocket=" + quote(pocket)) if pocket else ""}">{e(a)}</a>'

    @app.get("/og", response_class=HTMLResponse)
    def og_index(request: Request):
        ver = getv(request)
        head = ""
        if ver is not None:
            head = (f"<p>Version <b>{e(ver.id)}</b>: {ver.n_ogs():,} orthogroups, {ver.n_sogs():,} site orthogroups "
                    f"(SOGs) · <a href='/sog?{vq(ver)}'>browse SOGs</a>"
                    + (" · <span class='badge b-toy'>toy</span>" if ver.toy else "") + "</p>")
        return P("orthogroups", f"<h1>Orthogroups</h1>{head}{og_filter_form(request)}{og_list_html(request, ver)}", request)

    @app.get("/sog", response_class=HTMLResponse)
    def sog_index(request: Request):
        ver = getv(request)
        if ver is None:
            return P("SOGs", "<h1>SOGs</h1>" + og_list_html(request, None), request)
        qp = request.query_params
        size, pg = qint(request, "size", 50, 5, 500), qint(request, "page", 1)
        try:
            ms = float(qp["min_score"]) if qp.get("min_score") else None
        except ValueError:
            ms = None
        ev = qp.get("event", "")
        sort = qp.get("sort", "sog_id")
        rows, total = ver.sog_page(q=qp.get("q", ""), event=ev, case_only=qp.get("case") == "1", min_score=ms,
                                   sort=sort, desc=qp.get("desc", "1" if sort == "score" else "0") == "1",
                                   page=pg, size=size)
        evc = ver.event_counts()
        evopts = '<option value="">any event</option><option value="!none"%s>any non-none event</option>' % (
            " selected" if ev == "!none" else "") + "".join(
            f'<option value="{k}"{" selected" if ev == k else ""}>{k} ({n:,})</option>' for k, n in sorted(evc.items()))
        form = (f'<form class="filters" method="get"><input name="q" value="{e(qp.get("q", ""))}" placeholder="SOG id or label">'
                f'<select name="event">{evopts}</select><input name="min_score" value="{e(qp.get("min_score", ""))}" placeholder="min S" size="6">'
                f'<label><input type="checkbox" name="case" value="1"{" checked" if qp.get("case") == "1" else ""}> has case</label>'
                f'<select name="sort">' + "".join(f'<option value="{c}"{" selected" if sort == c else ""}>sort: {c}</option>'
                                                  for c in ("sog_id", "score", "n_P", "n_A", "n_species_sites", "footprint_size"))
                + f'</select><input type="hidden" name="v" value="{e(ver.id)}"><button>filter</button></form>')
        t = table(rows, [("sog_id", "SOG", lambda s: f'<a href="/og/{e(s["og_id"])}?{vq(ver)}&sog={quote(s["sog_id"])}">{e(s["sog_id"])}</a>'),
                         ("label", "label"), ("footprint_size", "footprint"),
                         ("counts", "P/W/A/U/G", lambda s: "/".join(str(s.get(f"n_{k}", "")) for k in STATE_ORDER)),
                         ("event", "event"), ("event_node", "event node"), ("score", "S", lambda s: fmt(s.get("score"), 3)),
                         ("case_id", "case", lambda s: f'<a href="/case/{e(s["case_id"])}">{e(s["case_id"])}</a>' if s.get("case_id") else "")])
        summ = ", ".join(f"{k}: {n:,}" for k, n in sorted(evc.items()))
        return P("SOGs", f"<h1>Site orthogroups · {e(ver.id)}</h1><p class='muted'>All SOGs of the version "
                         f"({summ}). Reads a column subset of sogs.parquet; states are read per orthogroup only.</p>"
                         f"{form}{t}{pager(request, 'page', total, size, pg)}", request)

    @app.get("/api/versions")
    def api_versions():
        return [{"id": v.id, "dir": str(v.dir), "toy": v.toy, "manifest": v.manifest} for v in ogs.versions().values()]

    @app.get("/api/sogs")
    def api_sogs(request: Request, q: str = "", event: str = "", case: bool = False, min_score: float | None = None,
                 sort: str = "sog_id", desc: bool = False, page: int = 1, size: int = 50):
        ver = getv(request)
        if ver is None:
            raise HTTPException(404, "no orthogroup data")
        size = max(1, min(size, 500))
        rows, total = ver.sog_page(q=q, event=event, case_only=case, min_score=min_score, sort=sort, desc=desc,
                                   page=max(1, page), size=size)
        return jclean({"version": ver.id, "total": total, "page": page, "size": size, "rows": rows})

    @app.get("/api/ogs")
    def api_ogs(request: Request, q: str = "", flags: str = "", half: str = "", sort: str = "og_id",
                desc: bool = False, page: int = 1, size: int = 50):
        ver = getv(request)
        if ver is None:
            raise HTTPException(404, "no orthogroup data")
        rows, total = ver.og_page(q=q, flags=flags, half=half, sort=sort, desc=desc, page=max(1, page),
                                  size=max(1, min(size, 500)))
        return jclean({"version": ver.id, "total": total, "rows": rows})

    @app.get("/api/og/{og_id}")
    def api_og(og_id: str, request: Request, sp: int = 1, n: int = 50, sog: str | None = None):
        """OG record + a page of SOG columns (`sp`, `n`) with their member rows, or one SOG (`sog`)."""
        ver = getv(request)
        d = ogs.get(og_id, ver, sog_ids=[sog] if sog else None, sog_page=None if sog else (max(1, sp), max(1, min(n, 200))))
        if d is None:
            raise HTTPException(404, f"orthogroup {og_id} not found")
        d["taxa"] = sorted(d["taxa"])
        return JSONResponse(jclean(json.loads(json.dumps({k: v for k, v in d.items() if k != "species"}
                                                  | {"species": list(d["species"].values())}, default=str))))

    SOG_PER_PAGE = 30

    @app.get("/og/{og_id}", response_class=HTMLResponse)
    def og_view(og_id: str, request: Request):
        ver = ogs.find_version(og_id, getv(request))
        if ver is None:
            raise HTTPException(404, f"orthogroup {og_id} not found")
        og, all_sogs = ver.og(og_id), ver.og_sogs(og_id)
        qp = request.query_params
        focus = qp.get("sog")
        evf = qp.get("ev", "")
        sogs_f = [s for s in all_sogs if (evf != "!none" or s.get("event") not in (None, "none"))]
        if focus:
            sogs = [s for s in all_sogs if s["sog_id"] == focus]
            sp_total, sp = 1, 1
        else:
            per = qint(request, "n", SOG_PER_PAGE, 5, 100)
            sp = qint(request, "sp", 1)
            sp_total = max(1, -(-len(sogs_f) // per))
            sp = min(sp, sp_total)
            sogs = sogs_f[(sp - 1) * per: sp * per]
        d = ogs.get(og_id, ver, sog_ids=[s["sog_id"] for s in sogs])
        full = qp.get("full") == "1"
        tree, order, cells = matrix(d)
        if full:
            tree = d["tree"]
            order = leaves(tree)
        dy, cw = 20, 26
        segs, leaf_y, tw, labels = layout(tree, dx=16, dy=dy)
        name_w, x0 = 190, tw + 8
        mx = x0 + name_w
        W = mx + cw * len(sogs) + 20
        H = dy * len(order) + 70
        top = 60
        evmark = {"gain": "▲", "loss": "▼", "specificity": "◆", "convergence": "✚"}
        svg = [f'<svg class="ogm" width="{W}" height="{H}" viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg">']
        for i, s in enumerate(sogs):
            x = mx + i * cw + cw / 2
            mk = evmark.get(s.get("event"), "")
            svg.append(f'<a href="/og/{e(og_id)}?{vq(ver)}&sog={quote(s["sog_id"])}"><text x="{x}" y="{top - 6}" '
                       f'transform="rotate(-45 {x} {top - 6})" class="colh"><title>{e(s["sog_id"])}: {e(s.get("label") or "")} '
                       f'[{e(s.get("event"))}]</title>{e(s["sog_id"].split(":")[-1])}{mk}</text></a>')
        g = [f'<g transform="translate(0,{top})">']
        for x1, y1, x2, y2 in segs:
            g.append(f'<line x1="{x1 + 4}" y1="{y1}" x2="{x2 + 4}" y2="{y2}" class="br"/>')
        for x, y, lab in labels:
            g.append(f'<text x="{x + 6}" y="{y - 3}" class="clade">{e(lab)}</text>')
        for tx in order:
            y = leaf_y[tx]
            sp_ = d["species"].get(int(tx)) if tx.lstrip("-").isdigit() else None
            nm = sp_["name"] if sp_ else tx
            g.append(f'<text x="{x0}" y="{y + 4}" class="sp"><title>taxid {e(tx)}</title>{e(nm)}</text>')
            for i, s in enumerate(sogs):
                rows = cells.get(s["sog_id"], {}).get(tx)
                cx = mx + i * cw + 3
                if not rows:
                    g.append(f'<rect x="{cx}" y="{y - dy / 2 + 2}" width="{cw - 6}" height="{dy - 4}" class="norow"/>')
                    continue
                stt = strongest(rows)
                tip = "\n".join(f'{r["accession"] or "-"}: {r["state"]} ({STATE_NAMES[r["state"]]})'
                                + (f' J={r["jaccard"]:.2f}' if r.get("jaccard") is not None else "")
                                + (f' pLDDT={r["footprint_plddt"]:.0f}' if r.get("footprint_plddt") is not None else "")
                                + (f' [{r["reason"]}]' if r.get("reason") else "") for r in rows[:12]) + (f"\n+{len(rows) - 12} more" if len(rows) > 12 else "")
                rect = (f'<rect x="{cx}" y="{y - dy / 2 + 2}" width="{cw - 6}" height="{dy - 4}" '
                        f'fill="{STATE_COLORS[stt]}" class="cell s{stt}"><title>{e(s["sog_id"])} · {e(nm)}\n{e(tip)}</title></rect>'
                        + (f'<text x="{cx + (cw - 6) / 2}" y="{y + 4}" class="ct s{stt}">{stt}{len(rows)}</text>' if len(rows) > 1 else ""))
                g.append(rect)
        g.append("</g>")
        svg += g + ["</svg>"]
        legend = " ".join(f'<span class="lg s{s}" style="background:{STATE_COLORS[s]}"></span>{s} {STATE_NAMES[s]}'
                          for s in STATE_ORDER)
        toy = "<p class='toy'>TOY fixture data — synthetic, not a result.</p>" if d["manifest"].get("toy") else ""
        sog_t = table(sogs, [("sog_id", "SOG", lambda s: f'<a href="?{vq(ver)}&sog={quote(s["sog_id"])}">{e(s["sog_id"])}</a>'),
                             ("label", "label"), ("footprint_ref", "footprint (ref numbering)"),
                             ("counts", "P/W/A/U/G", lambda s: "/".join(str(s.get(f"n_{k}", "")) for k in STATE_ORDER)),
                             ("event", "event"), ("event_node", "event node"), ("score", "S"),
                             ("case_id", "case", lambda s: f'<a href="/case/{e(s["case_id"])}">{e(s["case_id"])}</a>' if s.get("case_id") else "")])
        ref = og.get("ref_accession")
        refl = f" · reference {acc_link(ref)}" if ref else ""
        # navigation: SOG column pages, event filter, full-tree toggle
        nav = []
        if focus:
            nav.append(f'<a href="?{vq(ver)}">« all SOGs of {e(og_id)}</a>')
        else:
            def lnk(txt, **kw):
                q2 = {"v": ver.id, **{k: v for k, v in qp.items() if k in ("n", "ev", "full")}, **kw}
                return f'<a href="?{urlencode({k: v for k, v in q2.items() if v not in (None, "")})}">{txt}</a>'
            nav.append(f"SOGs {len(sogs_f)}" + (f" of {len(all_sogs)}" if evf else "") + f" · page {sp}/{sp_total}")
            if sp > 1:
                nav.append(lnk("‹ prev", sp=sp - 1))
            if sp < sp_total:
                nav.append(lnk("next ›", sp=sp + 1))
            nav.append(lnk("only SOGs with an event" if evf != "!none" else "show all SOGs", ev="!none" if evf != "!none" else "", sp=1))
        nav.append(f'<a href="?{urlencode({k: v for k, v in {**dict(qp), "full": "0" if full else "1"}.items()})}">'
                   f'{"pruned tree (species with rows)" if full else "full species tree"}</a>')
        members = ""
        if focus and sogs:
            s = sogs[0]
            rows = sorted([r for r in d["states"] if r["sog_id"] == focus],
                          key=lambda r: (STATE_ORDER.index(r["state"]), r["taxid"]))

            def spn(r):
                sp_ = d["species"].get(r["taxid"])
                return e(sp_["name"] if sp_ else r["taxid"])
            nrows = len(rows)
            mt = table(rows[:500], [("taxid", "species", spn),
                              ("state", "state", lambda r: f'<span class="st s{e(r["state"])}">{e(r["state"])}</span>'),
                              ("accession", "member", lambda r: acc_link(r.get("accession"), r.get("model_tier"), r.get("pocket_id"))
                               + ("" if r.get("model_tier") in (None, 1, 2, 3) else " <span class='muted'>(tier 4: no atlas model)</span>")),
                              ("site_id", "site"), ("pocket_id", "pocket", lambda r: f'<code>{e((r.get("pocket_id") or "").split(":")[-1])}</code>'),
                              ("jaccard", "Jaccard", lambda r: fmt(r.get("jaccard"), 2)),
                              ("footprint_plddt", "footprint pLDDT", lambda r: fmt(r.get("footprint_plddt"), 0)),
                              ("reason", "reason")])
            members = (f"<h2>Members of {e(focus)}</h2><p>footprint (reference numbering): <code>{e(s.get('footprint_ref'))}</code></p>{mt}" + (f"<p class='muted'>first 500 of {nrows} member rows (strongest states first); full list: JSON API.</p>" if nrows > 500 else ""))
        body = f"""<h1>{e(og["og_id"])} <small>{e(og.get("label"))}</small> <small class="muted">· {e(ver.id)}</small></h1>{toy}
<p>{og.get("n_species")} species · {og.get("n_members")} members · {og.get("n_sogs")} SOGs{refl}
{(" · flags: " + e(og["flags"])) if og.get("flags") else ""} · <a href="/api/og/{e(og_id)}?{vq(ver)}">JSON</a></p>
<h2>SOG × species</h2><p class="legend">{legend} · blank = no row · ▲ gain ▼ loss ◆ specificity ✚ convergence</p>
<p class="pager">{" · ".join(nav)}</p><div class="scroll">{"".join(svg)}</div>
{members}<h2>Site orthogroups{" (this page)" if not focus else ""}</h2>{sog_t}"""
        return P(og_id, body, request)

    # ------------------------------------------------------------------ cases
    def calib_badge(f):
        f = f or {}
        k = f.get("kind")
        if not f.get("text"):
            return ""
        return f'<span class="badge b-{"pass" if k == "fdr" else "pending"}">{e(f["text"])}</span>'

    def case_list_html():
        rows = cases.list()
        if not rows:
            return "<p class='muted'>No cases yet (cases/&lt;case_id&gt;/case.yaml).</p>"
        return table(rows, [("case_id", "case", lambda r: f'<a href="/case/{e(r["case_id"])}">{e(r["case_id"])}</a>' + (" <span class='badge b-toy'>toy</span>" if r["toy"] else "")),
                            ("title", "title"), ("category", "category"), ("label", "label", lambda r: badge(r["label"])),
                            ("fdr", "calibration", lambda r: calib_badge(r.get("fdr")))])

    @app.get("/case", response_class=HTMLResponse)
    def case_index():
        return page("cases", f"<h1>Cases</h1>{case_list_html()}")

    def _md_table(text):
        """Parse the first pipe table of a markdown text -> (preamble lines, header, rows, trailing lines)."""
        lines = text.splitlines()
        idx = [i for i, ln in enumerate(lines) if ln.startswith("|")]
        if not idx:
            return lines, [], [], []
        a, b = idx[0], idx[-1]
        cells = lambda ln: [x.strip() for x in ln.strip().strip("|").split(" | ")]  # noqa: E731
        return lines[:a], cells(lines[a]), [cells(ln) for ln in lines[a + 2:b + 1]], lines[b + 1:]

    @app.get("/validation", response_class=HTMLResponse)
    def validation_index():
        rows = [r for r in cases.list() if r["case_id"].startswith("val_")]
        crow = table(rows, [("case_id", "case", lambda r: f'<a href="/case/{e(r["case_id"])}">{e(r["case_id"])}</a>'),
                            ("title", "title"), ("label", "pipeline outcome", lambda r: badge(r["label"]))]) if rows else \
            "<p class='muted'>No val_* cases found.</p>"
        md = REPO / "results" / "cases" / cfg.snapshot / "validation_table.md"
        if md.is_file():
            pre, head, trs, post = _md_table(md.read_text())
            doi = lambda t: re.sub(r"(10\.\d{4,9}/[^\s;]+)", lambda m: f'<a href="https://doi.org/{m.group(1)}">{m.group(1)}</a>', e(t))  # noqa: E731
            tab = ('<div class="scroll"><table class="valtab"><thead><tr>' + "".join(f"<th>{e(h)}</th>" for h in head) + "</tr></thead><tbody>"
                   + "".join("<tr>" + "".join(f"<td>{doi(x) if head[i:i+1] == ['DOIs'] else e(x)}</td>" for i, x in enumerate(r)) + "</tr>" for r in trs)
                   + "</tbody></table></div>")
            intro = "".join(f"<p>{e(x)}</p>" for x in pre if x.strip() and not x.startswith("#"))
            outro = "".join(f"<p class='muted'>{e(x)}</p>" for x in post if x.strip())
            tab = f"{intro}{tab}{outro}<p class='muted'>Source: <code>{e(md.relative_to(REPO) if REPO in md.parents else md)}</code></p>"
        else:
            tab = f"<p class='muted'>Recall table not found at <code>{e(md)}</code>.</p>"
        body = ("<h1>Validation (known events)</h1><p class='toy'>Annotation-first validation of the unchanged pipeline against published "
                "events. Not discoveries, not FDR-calibrated.</p>"
                f"<h2>Cases</h2>{crow}<h2>Recall table</h2>{tab}")
        return page("validation", body)

    @app.get("/api/validation")
    def api_validation():
        md = REPO / "results" / "cases" / cfg.snapshot / "validation_table.md"
        pre, head, trs, post = _md_table(md.read_text()) if md.is_file() else ([], [], [], [])
        return JSONResponse(jclean({"cases": [r for r in cases.list() if r["case_id"].startswith("val_")], "header": head,
                                    "rows": [dict(zip(head, r)) for r in trs]}))

    @app.get("/api/case/{case_id}")
    def api_case(case_id: str):
        c = cases.get(case_id)
        if c is None:
            raise HTTPException(404, f"case {case_id} not found")
        return JSONResponse(json.loads(json.dumps({k: v for k, v in c.items() if k != "_dir"}, default=str)))

    @app.get("/case/{case_id}/file/{rel:path}")
    def case_file(case_id: str, rel: str):
        p = cases.file(case_id, rel)
        if p is None:
            raise HTTPException(404, "file not found in case directory")
        mt = mimetypes.guess_type(p.name)[0]
        if p.suffix == ".gz" and p.name[:-3].endswith((".cif", ".mmcif")):
            import gzip
            return Response(gzip.decompress(p.read_bytes()), media_type="chemical/x-mmcif")
        if p.suffix in (".cif", ".mmcif"):
            mt = "chemical/x-mmcif"
        return FileResponse(p, media_type=mt or "application/octet-stream")

    @app.get("/case/{case_id}", response_class=HTMLResponse)
    def case_view(case_id: str):
        c = cases.get(case_id)
        if c is None:
            raise HTTPException(404, f"case {case_id} not found")
        cid = e(case_id)
        none = "<p class='muted'>none recorded</p>"

        def lst(key):
            v = c.get(key) or []
            return v if isinstance(v, list) else [v]

        def acc_link(a):
            if not a:
                return ""
            return e(a) if str(a).startswith("TOY") else f'<a href="/protein/{e(a)}">{e(a)}</a>'
        prot = table(lst("proteins"), [("accession", "accession", lambda r: acc_link(r.get("accession"))),
                                       ("species", "species"), ("taxid", "taxid"), ("role", "role"),
                                       ("state", "state", lambda r: f'<span class="st s{e(r.get("state"))}">{e(r.get("state"))}</span>'),
                                       ("site_id", "site"), ("residues", "residues"), ("note", "note")]) if lst("proteins") else none
        art = table(lst("artifact_checks"), [("name", "check"), ("status", "status", lambda r: badge(r.get("status"))),
                                             ("detail", "detail")]) if lst("artifact_checks") else none

        def ev_link(r):
            t, i = (r.get("type") or "").lower(), r.get("id") or ""
            url = {"pdb": f"https://www.rcsb.org/structure/{i}", "uniprot": f"https://www.uniprot.org/uniprotkb/{i}",
                   "chembl": f"https://www.ebi.ac.uk/chembl/explore/compound/{i}", "doi": f"https://doi.org/{i}",
                   "pmid": f"https://pubmed.ncbi.nlm.nih.gov/{i}", "mcsa": f"https://www.ebi.ac.uk/thornton-srv/m-csa/entry/{i}"}.get(t)
            return f'<a href="{e(url)}">{e(i)}</a>' if url else e(i)

        def cfile(r):
            f = r.get("file")
            return f'<a href="/case/{cid}/file/{e(f)}">{e(f)}</a>' if f else ""
        def filter_evidence(r):
            ids = {str(x) for x in (r.get("evidence") or [])}
            hits = [x for x in lst("evidence") if x.get("filter") == r.get("name") or str(x.get("id")) in ids]
            return "<br>".join(f'{e(x.get("type"))} {ev_link(x)}' + (f' <span class="muted">{e(x.get("detail"))}</span>' if x.get("detail") else "")
                               + (f' · {cfile(x)}' if x.get("file") else "") for x in hits) or "<span class='muted'>none linked</span>"
        filt = table(lst("filters"), [("name", "filter"), ("status", "status", lambda r: badge(r.get("status"))),
                                      ("detail", "detail"), ("ev", "evidence", filter_evidence)]) if lst("filters") else none
        ev = table(lst("evidence"), [("type", "type"), ("id", "id", ev_link), ("about_change", "about the change"),
                                     ("filter", "supports filter"), ("verified", "verified"), ("source", "source"),
                                     ("detail", "detail"), ("file", "cached record", cfile)]) if lst("evidence") else none
        cof_rows = lst("cofold")

        def cof_btn(r):
            f = r.get("file")
            return cfile(r) if f else ""
        def kv(d):
            return "".join(f"<li><b>{e(k)}</b>: {e(v if not isinstance(v, (dict, list)) else json.dumps(v, default=str))}</li>"
                           for k, v in d.items())

        def cof_card(r):
            if not isinstance(r, dict):
                return f"<p>{e(r)}</p>"
            st = r.get("status")
            bl = r.get("blinding")
            bls = (bl.get("status") if isinstance(bl, dict) else bl) or ""
            head = f'<h3>{e(r.get("model_id"))} <span class="muted">({e(r.get("role"))})</span> verdict {badge(st)}'
            if bls:
                head += f' blinding <span class="badge b-{"fail" if str(bls).upper() == "UNBLINDED" else "pass"}">{e(bls)}</span>'
            head += "</h3>"
            out = [head]
            if st and "uninformative" in str(st):
                out.append("<p class='toy'>Uninformative: the co-folding ligand does not sit at the case footprint, so this result "
                           "says nothing for or against the gain/loss.</p>")
            if str(st) == "pending":
                out.append("<p class='muted'>Co-folding not run yet (placeholder).</p>")
            if isinstance(bl, dict):
                out.append(f"<p><b>Blinding {e(bls)}.</b> {e(bl.get('honest_claim', ''))}</p>"
                           f"<details><summary>what happened</summary><p>{e(bl.get('what_actually_happened', ''))}</p>"
                           f"<p>{e(bl.get('remedy', ''))}</p></details>")
            if r.get("note"):
                out.append(f"<p>{e(r['note'])}</p>")
            if r.get("decision_rule"):
                out.append(f"<p>Decision rule: {e(r['decision_rule'])}"
                           + (f" · seeds gain {e(r.get('seeds_gain'))} · seeds loss {e(r.get('seeds_loss'))}" if r.get("seeds_gain") is not None else "") + "</p>")
            if r.get("interpretation"):
                out.append("<ul>" + "".join(f"<li>{e(x)}</li>" for x in (r["interpretation"] if isinstance(r["interpretation"], list) else [r["interpretation"]])) + "</ul>")
            fvc = r.get("footprint_vs_cognate_site")
            if isinstance(fvc, dict) and fvc.get("conclusion"):
                out.append(f"<p><b>Footprint vs cognate ligand site:</b> {e(fvc['conclusion'])}</p>")
            if r.get("metrics"):
                out.append(f"<p>metrics: {e(', '.join(f'{k}={v}' for k, v in r['metrics'].items()))}</p>")
            if r.get("file"):
                out.append(f"<p>model file: {cfile(r)}</p>")
            shown = {"model_id", "role", "status", "blinding", "note", "decision_rule", "seeds_gain", "seeds_loss", "interpretation",
                     "metrics", "file"}
            rest = {k: v for k, v in r.items() if k not in shown}
            if rest:
                out.append(f"<details><summary>all other co-fold fields ({len(rest)})</summary><ul>{kv(rest)}</ul></details>")
            return f'<div class="cofcard">{"".join(out)}</div>'
        cof = "".join(cof_card(r) for r in cof_rows) if cof_rows else none
        # pocket residues per accession from the case's protein list ("A:10 A:18" or "10 18")
        res_by_acc = {}
        for pr in lst("proteins"):
            nums = [int(x.split(":")[-1]) for x in str(pr.get("residues") or "").split() if x.split(":")[-1].isdigit()]
            if nums:
                res_by_acc.setdefault((pr.get("accession"), pr.get("role")), nums)
                res_by_acc.setdefault((pr.get("accession"), None), nums)
        cof_json = [{"model_id": r.get("model_id"), "accession": r.get("accession"), "ligand": r.get("ligand"),
                     "seed": r.get("seed"), "role": r.get("role"), "metrics": r.get("metrics") or {},
                     "url": f"/case/{cid}/file/{r['file']}", "residues": res_by_acc.get((r.get("accession"), r.get("role")))
                     or res_by_acc.get((r.get("accession"), None)) or []}
                    for r in cof_rows if r.get("file")]
        figs = []
        for f in lst("figures"):
            fp = f.get("file", "")
            url = f"/case/{cid}/file/{e(fp)}"
            inner = f'<a href="{url}">{e(fp)}</a>' if fp.endswith(".pdf") else f'<img src="{url}" alt="{e(f.get("caption"))}">'
            figs.append(f"<figure>{inner}<figcaption>{e(f.get('caption'))}</figcaption></figure>")

        def lit(r):
            if r.get("doi"):
                return f'<a href="https://doi.org/{e(r["doi"])}">{e(r["doi"])}</a>'
            if r.get("pmid"):
                return f'<a href="https://pubmed.ncbi.nlm.nih.gov/{e(r["pmid"])}">PMID {e(r["pmid"])}</a>'
            return e(r.get("url", ""))
        litt = table(lst("literature"), [("ref", "reference", lit), ("title", "title"), ("note", "note")]) if lst("literature") else none
        og_link = ""
        if c.get("og_id"):
            og_link = f'Orthogroup <a href="/og/{e(c["og_id"])}">{e(c["og_id"])}</a>' + (f' · SOG <code>{e(c.get("sog_id"))}</code>' if c.get("sog_id") else "")
        warn = f"<p class='toy'>case.yaml is missing required keys: {e(c['_missing_required'])}</p>" if c["_missing_required"] else ""
        toy = "<p class='toy'>TOY fixture case — synthetic, not a result.</p>" if c.get("toy") else ""
        extra = f"<h2>Other fields</h2><pre>{e(json.dumps(c['_extra'], indent=1, default=str))}</pre>" if c["_extra"] else ""
        def panel(role, title):
            return (f'<div class="cofpanel" data-role="{role}"><h3>{title}</h3><div class="vbar">'
                    f'<select aria-label="{role} model"></select> <label><input type="checkbox" class="plddt" checked> pLDDT</label> '
                    f'<label><input type="checkbox" class="lig" checked> zoom to ligand</label></div>'
                    f'<div class="viewer small"></div><div class="vmsg muted"></div><div class="cofinfo muted"></div></div>')
        viewer = ("<div class='sbs'>" + panel("ancestral", "Ancestral state") + panel("derived", "Derived state") + "</div>"
                  "<p class='legend'>Ligand: green sticks · residues within 4.5 Å of the ligand: cyan · case pocket residues: orange · "
                  "cartoon coloured by pLDDT (B-factor).</p>"
                  f'<script type="application/json" id="cofold-data">{json.dumps(cof_json).replace("</", "<\\/")}</script>'
                  if cof_json else "")
        fd = c.get("_fdr") or {}
        banner = ""
        if fd.get("text"):
            banner = (f"<p class='{'toy' if fd.get('kind') != 'fdr' else 'calib'}'><b>{e(fd['text'])}</b>"
                      + (f" · selection: {e(fd['selection'])}" if fd.get("selection") else "") + "</p>")
        if c.get("headline") or c.get("filter_label_with_literature_post_hoc"):
            banner += ("<p>Label from the preregistered filters only: <b>" + e(c.get("headline") or c.get("label"))
                       + "</b> (headline) · label after adding literature evidence: <b>"
                       + e(c.get("filter_label_with_literature_post_hoc")) + "</b> <span class='badge b-pending'>post hoc</span></p>")
        prov = (f"<h2>Provenance</h2><pre>{e(json.dumps(c['provenance'], indent=1, default=str))}</pre>"
                if c.get("provenance") else "")
        body = f"""<h1>{e(c.get("title"))}</h1>{toy}{warn}{banner}
<p><code>{cid}</code> · {e(c.get("category"))} · {e(c.get("track", ""))} · label {badge(c.get("label"))}
{(" · S = " + e(c.get("score"))) if c.get("score") is not None else ""} · {og_link} · <a href="/api/case/{cid}">JSON</a></p>
<p>{e(c.get("summary", ""))}</p>
<h2>Proteins</h2>{prot}<h2>Evidence filters</h2>{filt}<h2>Artifact checks</h2>{art}
<h2>Evidence</h2>{ev}<h2>Co-folding models</h2>{cof}{viewer}
<h2>Figures</h2>{"".join(figs) or none}<h2>Literature</h2>{litt}
<h2>Notes</h2><p>{e(c.get("notes", ""))}</p>{prov}{extra}"""
        sc = '<script src="/static/vendor/3Dmol-min.js"></script><script src="/static/viewer.js"></script>' if cof_rows else ""
        return page(case_id, body, sc)

    return app
