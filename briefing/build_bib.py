"""Build briefing/references.bib from Crossref metadata, one DOI per citation key.

    python briefing/build_bib.py          # needs network; prints each key with its Crossref title

Most DOIs come from record/LITERATURE.md (checked against Crossref during the project); the standard-method
references (Sankoff, Leiden, curveball, FAMSA, DIAMOND, Benjamini-Hochberg, ColabFold) were added for the briefing.
"""
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

KEYS = {
    "thornton2003": "10.1126/science.1086185", "bridgham2006": "10.1126/science.1123348",
    "babbitt1997": "10.1074/jbc.272.49.30591", "varadi2024": "10.1093/nar/gkad1011",
    "abramson2024": "10.1038/s41586-024-07487-w", "akdel2022": "10.1038/s41594-022-00849-w",
    "deorowicz2016": "10.1038/srep33964", "traag2019": "10.1038/s41598-019-41695-z",
    "zhang2024biolip": "10.1093/nar/gkad630", "utges2024": "10.1186/s13321-024-00923-z",
    "gao2013": "10.1371/journal.pcbi.1003302", "skolnick2013": "10.1073/pnas.1300011110",
    "vankempen2024": "10.1038/s41587-023-01773-0", "altenhoff2012": "10.1371/journal.pcbi.1002514",
    "nehrt2011": "10.1371/journal.pcbi.1002073", "sankoff1975": "10.1137/0128004",
    "strona2014": "10.1038/ncomms5114", "tunyasuvunakool2021": "10.1038/s41586-021-03828-1",
    "meyer2020": "10.1186/s12859-020-03855-1", "barriohernandez2023": "10.1038/s41586-023-06510-w",
    "thomas2022": "10.1002/pro.4218", "buchfink2021": "10.1038/s41592-021-01101-x",
    "masters2025": "10.1038/s41467-025-63947-5", "skrinjar2025": "10.1101/2025.02.03.636309",
    "lewis2001": "10.1080/106351501753462876", "doolittle1994": "10.1016/0968-0004(94)90167-8",
    "riziotis2025": "10.1111/febs.17332", "saldano2022": "10.1093/bioinformatics/btac202",
    "cohen2010": "10.1093/bioinformatics/btq549", "kim2026": "10.7554/eLife.110475.2",
    "harms2013": "10.1038/nrg3540", "uniprot2025": "10.1093/nar/gkae1010", "ribeiro2018": "10.1093/nar/gkx1012",
    "lau2024": "10.1126/science.adq4946", "benjamini1995": "10.1111/j.2517-6161.1995.tb02031.x",
    "mirdita2022": "10.1038/s41592-022-01488-1",
}


def esc(s):
    s = re.sub(r"<[^>]+>", "", s or "")
    s = s.replace("{", "").replace("}", "").replace("&amp;", "and").replace("&", "and")
    return " ".join(s.split())


def entry(key, doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers={"User-Agent": "agentic-pocket-evolution briefing (bib builder)"})
    m = json.load(urllib.request.urlopen(req, timeout=30))["message"]
    authors = []
    for a in m.get("author", []):
        if "family" in a:
            authors.append(f"{esc(a['family'])}, {esc(a.get('given', ''))}".rstrip(", "))
        elif "name" in a:
            authors.append("{" + esc(a["name"]) + "}")
    year = (m.get("published-print") or m.get("published-online") or m.get("issued"))["date-parts"][0][0]
    venue = esc((m.get("container-title") or [""])[0]) or esc(m.get("institution", [{}])[0].get("name", "")) \
        if isinstance(m.get("institution"), list) else esc((m.get("container-title") or [""])[0])
    if not venue and m.get("type") == "posted-content":
        venue = "bioRxiv (preprint)" if doi.startswith("10.1101/") else "Preprint"
    f = {"title": esc(m["title"][0]), "author": " and ".join(authors), "journal": venue, "year": str(year),
         "volume": m.get("volume", ""), "pages": m.get("page", ""), "doi": doi}
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in f.items() if v)
    return f"@article{{{key},\n{body}\n}}\n", f["title"]


if __name__ == "__main__":
    out = []
    for k, d in KEYS.items():
        e, t = entry(k, d)
        out.append(e)
        print(f"{k:22s} {t[:90]}")
    (Path(__file__).parent / "references.bib").write_text("\n".join(out))
