#!/usr/bin/env python
"""Select context UniProt reference proteomes for Ancestra run O2 (tier-1 + context).

Rule (deterministic given the snapshot manifest, the UniProt reference-proteome README of the release, the
Foldcomp v4 entry table and the NCBI taxdump; all parameters in config/orthology.yaml `context:`):

Common candidate filter: UniProt reference proteome (Bacteria/Archaea/Eukaryota), taxid not in tier 1 and not
the same NCBI species as a tier-1 proteome, min_proteins <= #canonical proteins <= max_proteins, no MAG-like
placeholder name.

(a) Eukaryotic densification (pockets_source = local unless already atlas-complete). Slots (clade taxid, n)
    are filled in config order. Inside a slot: priority = has a UniProt mnemonic (OSCODE; established
    reference organism) > #additional isoforms in the release (annotation depth) > smaller taxid; at most one
    proteome per NCBI order and per genus within the slot (diversity); accepted only if >= min_fc_coverage of
    its accessions are in fc_entries_afdb_uniprot_v4 (structures exist -> pockets can be computed locally
    with the atlas worker) or already in the atlas. Chosen regardless of atlas coverage (PI decision: atlas-
    complete eukaryotes are only protists/early fungi and cannot densify clades near human/fly/plants).
(b) Atlas-complete taxa (pockets_source = atlas): prefilter (#atlas accessions of the taxid)/(#proteins) >=
    min_size_ratio; groups protists/bacteria/archaea with quotas; priority d = #standard-rank taxa
    (domain..genus) shared with the closest tier-1 species (densify), then OSCODE, then taxid; round-robin over
    phyla (ordered by their best candidate); caps max_per_class per class, max_per_genus per genus;
    accepted only if >= min_atlas_overlap of the proteome accessions are present in the snapshot (any tier).
    Stops at atlas_max_total or when tier-1 + context proteins would exceed max_total_proteins.

Outputs (out_dir): context_species.tsv (accepted; column pockets_source = atlas|local), context_considered.tsv
(every candidate evaluated, with status), o2_proteomes.tsv (tier 1 + context with md5), context_selection.json.

Usage: .venv/bin/python scripts/orthology/select_context_species.py [--config config/orthology.yaml]
"""
from __future__ import annotations

import argparse
import json
import os
import re

import duckdb
import pandas as pd

from _common import P, atlas_accessions, fasta_accessions, fetch, load_config, log, md5, snapshot_id  # noqa: F401
from pocketevo.taxonomy import Taxonomy, taxdump_info

STD_RANKS = ("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
PLACEHOLDER = re.compile(r"(^uncultured )|( bacterium\b)|( archaeon\b)|(^unclassified )|( sp\.( |$))|"
                         r"(environmental)|(unidentified)|(^Candidatus .* (bacterium|archaeon))", re.I)
KINGDOM_DIR = {"bacteria": "Bacteria", "archaea": "Archaea", "eukaryota": "Eukaryota"}


def parse_readme(path: str) -> tuple[pd.DataFrame, str]:
    lines = open(path).read().split("\n")
    rel = next((l.strip() for l in lines if l.startswith("Release ")), "")
    i = next(k for k, l in enumerate(lines) if l.startswith("Proteome_ID\tTax_ID"))
    rows = []
    for l in lines[i + 1:]:
        f = l.split("\t")
        if len(f) >= 8 and f[0].startswith("UP"):
            rows.append(f[:8])
        elif rows and l.strip() and not l.startswith("UP"):
            break
    r = pd.DataFrame(rows, columns=["proteome", "taxid", "oscode", "superregnum", "n1", "n2", "n3", "species"])
    for c in ("taxid", "n1", "n2", "n3"):
        r[c] = r[c].astype(int)
    r["oscode"] = r["oscode"].replace({"None": ""})
    return r, rel


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config")
    args = ap.parse_args()
    cfg = load_config(args.config)
    cc = cfg["context"]
    W = cfg["work_dir"]
    out_dir = cfg["out_dir"]
    os.makedirs(out_dir, exist_ok=True)
    pdir = os.path.join(W, "proteomes")
    fce = os.path.join(cfg["atlas_root"], "fc_entries_afdb_uniprot_v4.parquet")

    readme = fetch(cfg["uniprot_refprot_base"] + "/README", os.path.join(W, "raw", "README"))
    ref, release = parse_readme(readme)
    log(f"{len(ref)} reference proteomes in README ({release})")

    tax = Taxonomy(cfg["taxdump_dir"])
    t1 = pd.read_csv(cfg["tier1_proteomes"], sep="\t", header=None,
                     names=["taxid", "proteome", "ptype", "species", "superregnum", "genes", "proteins"])
    t1_taxids = set(t1.taxid.astype(int))
    t1_species = {tax.at_rank(t, "species") for t in t1_taxids} - {None}
    t1_files = {r.taxid: os.path.join(cfg["tier1_fasta_dir"], f"{r.proteome}_{r.taxid}.fasta.gz")
                for r in t1.itertuples(index=False)}
    t1_n = {t: len(fasta_accessions(f)) for t, f in t1_files.items()}
    budget_used = sum(t1_n.values())

    def stdlin(t):
        return [x for x in tax.lineage(t) if tax.rank.get(x) in STD_RANKS]
    t1_lin = {t: stdlin(t) for t in sorted(t1_taxids)}

    acc_pq = atlas_accessions(cfg)
    con = duckdb.connect()
    counts = con.execute(f"SELECT taxid, count(DISTINCT accession) n_atlas FROM '{acc_pq}' GROUP BY 1").df()
    ref = ref[ref.superregnum.isin(KINGDOM_DIR)].merge(counts, on="taxid", how="left").fillna({"n_atlas": 0})
    ref["size_ratio"] = ref.n_atlas / ref.n1
    con.execute(f"CREATE TEMP TABLE atl AS SELECT accession, min(tier) tier, max(n_pockets_kept) npk "
                f"FROM '{acc_pq}' GROUP BY 1")

    recs = []
    for r in ref.itertuples(index=False):
        if r.taxid in t1_taxids:
            continue
        lin = tax.lineage(r.taxid)
        if not lin or tax.at_rank(r.taxid, "species") in t1_species:
            continue
        if not (cc["min_proteins"] <= r.n1 <= cc["max_proteins"]):
            continue
        if cc.get("exclude_placeholder_names") and PLACEHOLDER.search(r.species):
            continue
        sl = stdlin(r.taxid)
        best_d, anchor = -1, None
        for t, l in t1_lin.items():
            d = 0
            for a, b in zip(sl, l):
                if a != b:
                    break
                d += 1
            if d > best_d:
                best_d, anchor = d, t
        recs.append(dict(taxid=r.taxid, proteome=r.proteome, species=r.species, superregnum=r.superregnum,
                         oscode=r.oscode, n_isoforms=r.n2, n_proteins_readme=r.n1, n_atlas_taxid=int(r.n_atlas),
                         size_ratio=round(r.size_ratio, 3), lin=set(lin),
                         phylum=tax.name_at_rank(r.taxid, "phylum") or "(none)",
                         class_=tax.name_at_rank(r.taxid, "class") or tax.name_at_rank(r.taxid, "order") or "(none)",
                         order_=tax.name_at_rank(r.taxid, "order") or "(none)",
                         genus=tax.name_at_rank(r.taxid, "genus") or r.species.split()[0],
                         shared_depth=best_d, nearest_tier1=anchor, nearest_tier1_name=tax.name.get(anchor, ""),
                         lineage=tax.lineage_names(r.taxid)))
    cand = pd.DataFrame(recs)
    cand["has_oscode"] = (cand.oscode != "").astype(int)
    log(f"{len(cand)} candidates after common filters")

    evaluated: dict[int, dict] = {}

    def check(row) -> dict:
        if row.taxid in evaluated:
            return evaluated[row.taxid]
        url = (f"{cfg['uniprot_refprot_base']}/{KINGDOM_DIR[row.superregnum]}/{row.proteome}/"
               f"{row.proteome}_{row.taxid}.fasta.gz")
        dest = os.path.join(pdir, f"{row.proteome}_{row.taxid}.fasta.gz")
        try:
            fetch(url, dest)
        except RuntimeError:
            evaluated[row.taxid] = dict(n_proteins=0, atlas_overlap=0.0, fc_coverage=0.0, frac_tier3=0.0,
                                        frac_with_pockets=0.0, status="download failed")
            return evaluated[row.taxid]
        accs = pd.DataFrame({"accession": sorted(set(fasta_accessions(dest)))})
        con.register("q", accs)
        n, ov, t3, wp = con.execute("""SELECT count(*), count(a.accession), count(*) FILTER (WHERE a.tier=3),
                                       count(*) FILTER (WHERE a.npk>0)
                                       FROM q LEFT JOIN atl a USING (accession)""").fetchall()[0]
        fc = con.execute(f"""WITH f AS (SELECT DISTINCT e.accession FROM read_parquet('{fce}') e
                                        WHERE e.accession IN (SELECT accession FROM q))
                             SELECT count(*) FROM q WHERE accession IN (SELECT accession FROM atl)
                                OR accession IN (SELECT accession FROM f)""").fetchall()[0][0]
        con.unregister("q")
        n = max(n, 1)
        res = dict(n_proteins=n, atlas_overlap=round(ov / n, 4), fc_coverage=round(fc / n, 4),
                   frac_tier3=round(t3 / n, 4), frac_with_pockets=round(wp / n, 4), fasta=dest)
        evaluated[row.taxid] = res
        return res

    chosen: list[tuple] = []  # (row, part)

    def accept(row, part, res):
        nonlocal budget_used
        chosen.append((row, part))
        budget_used += res["n_proteins"]

    # ---------------- (a) eukaryotic densification slots ----------------
    euk = cand[cand.superregnum == "eukaryota"].sort_values(["has_oscode", "n_isoforms", "taxid"],
                                                            ascending=[False, False, True])
    for slot in cc["eukaryote_slots"]:
        sub = euk[euk.lin.map(lambda s: slot["taxid"] in s)]
        got, orders, genera = 0, set(), set()
        taken = {r.taxid for r, _ in chosen}
        for row in sub.itertuples(index=False):
            if got >= slot["n"]:
                break
            if row.taxid in taken or row.order_ in orders or row.genus in genera:
                continue
            if row.n_proteins_readme > cc.get("max_proteins_slot", cc["max_proteins"]):
                continue
            res = check(row)
            ok = res.get("fc_coverage", 0) >= cc["min_fc_coverage"]
            res["status"] = "accepted" if ok else "low fc_entries coverage"
            res["part"] = f"a:{slot['name']}"
            log(f"  a {slot['name']:16s} {row.proteome}_{row.taxid} {row.species[:45]:45s} n={res['n_proteins']} "
                f"fc={res['fc_coverage']:.3f} atlas={res['atlas_overlap']:.3f} {res['status']}")
            if ok:
                accept(row, f"a:{slot['name']}", res)
                got += 1
                orders.add(row.order_)
                genera.add(row.genus)
        log(f"slot {slot['name']}: {got}/{slot['n']} (candidates {len(sub)})")

    # ---------------- (b) atlas-complete taxa ----------------
    def group_of(lin):
        for g in cc["atlas_groups"]:
            if g["taxid"] in lin and not any(x in lin for x in g.get("exclude", [])):
                return g["name"]
        return None
    pb = cand[cand.size_ratio >= cc["min_size_ratio"]].copy()
    pb["group"] = pb.lin.map(group_of)
    pb = pb[pb.group.notna()].sort_values(["shared_depth", "has_oscode", "taxid"], ascending=[False, False, True])
    n_class: dict = {}
    n_genus: dict = {}

    def picker(sub):
        """Generator over accepted rows of one group: round-robin over phyla (ordered by best candidate)."""
        phyla = list(dict.fromkeys(sub.phylum))
        queues = {ph: list(sub[sub.phylum == ph].itertuples(index=False)) for ph in phyla}
        while any(queues.values()):
            for ph in phyla:
                q = queues[ph]
                while q:
                    row = q.pop(0)
                    if row.taxid in {r.taxid for r, _ in chosen}:
                        continue
                    if n_class.get(row.class_, 0) >= cc["max_per_class"] or \
                            n_genus.get(row.genus, 0) >= cc["max_per_genus"]:
                        continue
                    if budget_used + row.n_proteins_readme > cc["max_total_proteins"]:
                        continue
                    res = check(row)
                    ok = res["atlas_overlap"] >= cc["min_atlas_overlap"]
                    res["status"] = "accepted" if ok else "low atlas overlap"
                    res["part"] = f"b:{row.group}"
                    if ok:
                        yield row, res
                        break

    # groups are interleaved (one pick per group per round) so the protein budget is shared fairly;
    # a group that reaches its quota or runs dry passes its turns to the others.
    gens = {g["name"]: picker(pb[pb.group == g["name"]]) for g in cc["atlas_groups"]}
    quota = {g["name"]: g["quota"] for g in cc["atlas_groups"]}
    got = {g: 0 for g in gens}
    n_b = 0
    active = list(gens)
    while active and n_b < cc["atlas_max_total"]:
        under = [g for g in active if got[g] < quota[g]] or active   # after quotas: redistribute
        for g in under:
            if n_b >= cc["atlas_max_total"]:
                break
            nxt = next(gens[g], None)
            if nxt is None:
                active.remove(g)
                continue
            row, res = nxt
            accept(row, f"b:{g}", res)
            n_class[row.class_] = n_class.get(row.class_, 0) + 1
            n_genus[row.genus] = n_genus.get(row.genus, 0) + 1
            got[g] += 1
            n_b += 1
    for g in gens:
        log(f"group {g}: {got[g]}/{quota[g]} (candidates {(pb.group == g).sum()})")

    rows = []
    for i, (r, part) in enumerate(chosen):
        d = r._asdict()
        d.update(evaluated[r.taxid])
        d["order_selected"] = i
        d["part"] = part
        d["pockets_source"] = "atlas" if evaluated[r.taxid]["atlas_overlap"] >= cc["min_atlas_overlap"] else "local"
        rows.append(d)
    sel = pd.DataFrame(rows).rename(columns={"class_": "class", "order_": "order"})
    cols = ["taxid", "proteome", "species", "part", "pockets_source", "superregnum", "phylum", "class", "order",
            "genus", "n_proteins", "atlas_overlap", "fc_coverage", "frac_tier3", "frac_with_pockets", "shared_depth",
            "nearest_tier1", "nearest_tier1_name", "oscode", "lineage"]
    sel[cols].to_csv(os.path.join(out_dir, "context_species.tsv"), sep="\t", index=False)
    ev = cand[cand.taxid.isin(evaluated)].drop(columns=["lin", "has_oscode"]).copy()
    for k in ("part", "n_proteins", "atlas_overlap", "fc_coverage", "frac_tier3", "frac_with_pockets", "status"):
        ev[k] = ev.taxid.map(lambda t: evaluated[t].get(k))
    ev.rename(columns={"class_": "class", "order_": "order"}).to_csv(
        os.path.join(out_dir, "context_considered.tsv"), sep="\t", index=False)

    o2 = []
    for r in t1.itertuples(index=False):
        f = t1_files[r.taxid]
        o2.append(dict(set="tier1", proteome=r.proteome, taxid=r.taxid, species=r.species, file=f, md5=md5(f),
                       n_proteins=t1_n[r.taxid]))
    for r in sel.itertuples(index=False):
        f = evaluated[r.taxid]["fasta"]
        o2.append(dict(set="context", proteome=r.proteome, taxid=r.taxid, species=r.species, file=f, md5=md5(f),
                       n_proteins=r.n_proteins))
    o2 = pd.DataFrame(o2)
    o2.to_csv(os.path.join(out_dir, "o2_proteomes.tsv"), sep="\t", index=False)
    prov = dict(snapshot=snapshot_id(cfg["snapshot_manifest"]), uniprot_refprot_release=release,
                readme_md5=md5(readme), taxdump=taxdump_info(cfg["taxdump_dir"]), config=cc,
                n_candidates=int(len(cand)), n_evaluated=len(evaluated), n_selected=int(len(sel)),
                n_selected_by_source=sel.pockets_source.value_counts().to_dict(),
                n_proteins_context=int(sel.n_proteins.sum()), n_proteins_o2=int(o2.n_proteins.sum()),
                per_part=sel.part.str.split(":").str[0].value_counts().to_dict())
    json.dump(prov, open(os.path.join(out_dir, "context_selection.json"), "w"), indent=1, default=str)
    log(json.dumps({k: prov[k] for k in ("n_candidates", "n_evaluated", "n_selected", "n_selected_by_source",
                                          "n_proteins_context", "n_proteins_o2", "per_part")}, default=str))


if __name__ == "__main__":
    main()
