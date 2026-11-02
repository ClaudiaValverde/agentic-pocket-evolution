#!/usr/bin/env python
"""Coverage / bias report for one atlas snapshot.

Reads ONLY the units listed in a manifest written by scripts/atlas_snapshot.py. All aggregation is done
in DuckDB over column projections (no fingerprints are read), so it scales to the complete atlas.

Outputs in <out>/<snapshot_id>/:
  tier_summary.csv               structures / pockets / quality per tier
  pocket_quality_by_tier.csv     pocket-level score / pLDDT / volume quantiles per tier
  tier1_species.csv              per-species counts for tier 1
  tier3_<rank>.csv               tier-3 done counts by superkingdom / phylum / class, with done-vs-AFDB shares
  tier3_taxa_thresholds.csv      number of tier-3 taxa (taxid and species level) with >= N done proteins
  tier3_completeness.csv         taxa x (AFDB size bin) x (fraction of the taxon's AFDB entries present)
  tier3_unclassified.csv         unmapped / unclassified / metagenomic / placeholder-name shares
  tier3_done_vs_pending_units.csv unit-level size / taxa-per-unit comparison
  tier3_entry_size.csv           Foldcomp entry length (proxy for protein length) done vs pending
  quality_by_superkingdom.csv    tier-3 structure quality by superkingdom
  COVERAGE.md                    human-readable report incl. "Biases & caveats"

Usage:
  uv run --with duckdb --with pyarrow --with pandas python scripts/coverage_report.py \
      --manifest data/snapshot/<id>/manifest.parquet [--out results/coverage]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import duckdb
import pandas as pd
import pyarrow.parquet as pq

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIER_SOURCE = {1: "afdb_v6_tar", 2: "afdb_swissprot_v6_tar"}
# placeholder organism names typical of MAG / environmental entries ("Acidobacteriota bacterium", "uncultured X")
PLACEHOLDER_SQL = ("(regexp_matches(lower(tax_name), '(^uncultured )|( bacterium$)|( archaeon$)|"
                   "(^unclassified )|( sp\\.$)|(environmental)|(unidentified)'))")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, file=sys.stderr, flush=True)


def md_table(df: pd.DataFrame, floatfmt: str = "{:,.3g}", max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                cells.append("" if pd.isna(v) else (f"{v:,.0f}" if abs(v) >= 1000 else floatfmt.format(v)))
            elif isinstance(v, (int,)) and not isinstance(v, bool):
                cells.append(f"{v:,}")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--atlas-root", default=None, help="default: the root recorded in the manifest")
    ap.add_argument("--out", default=os.path.join(P, "results", "coverage"))
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--memory-limit", default="96GB")
    ap.add_argument("--tmp", default="/tmp/duckdb_coverage")
    args = ap.parse_args(argv)

    man_path = os.path.abspath(args.manifest)
    snap_dir = os.path.dirname(man_path)
    summary = json.load(open(os.path.join(snap_dir, "manifest.json")))
    atlas = args.atlas_root or summary["atlas_root"]
    snap_id = summary["snapshot_id"]
    out = os.path.join(os.path.abspath(args.out), snap_id)
    os.makedirs(out, exist_ok=True)
    os.makedirs(args.tmp, exist_ok=True)

    con = duckdb.connect()
    con.execute(f"SET threads={args.threads}; SET memory_limit='{args.memory_limit}'; "
                f"SET temp_directory='{args.tmp}'; SET preserve_insertion_order=false;")
    con.execute(f"CREATE TABLE man AS SELECT * FROM read_parquet('{man_path}')")
    man = con.execute("SELECT unit_id, tier, pockets_path, structures_path FROM man ORDER BY unit_id").df()
    sfiles = [os.path.join(atlas, p) for p in man.structures_path]
    pfiles = [os.path.join(atlas, p) for p in man.pockets_path]
    log(f"snapshot {snap_id}: {len(man)} units")

    # ---------------- lineage table (NCBI rankedlineage + merged) ----------------
    ncbi = os.path.join(atlas, "ncbi")
    con.execute(f"""
      CREATE TABLE lin AS SELECT
        CAST(column0 AS BIGINT) taxid, column1 tax_name, NULLIF(column2,'') species, NULLIF(column3,'') genus,
        NULLIF(column4,'') family, NULLIF(column5,'') "order", NULLIF(column6,'') class,
        NULLIF(column7,'') phylum, NULLIF(column8,'') kingdom,
        NULLIF(regexp_replace(column9, '[\t|]', '', 'g'),'') superkingdom
      FROM read_csv('{ncbi}/rankedlineage.dmp', sep='\t|\t', header=false, quote='', all_varchar=true)""")
    con.execute(f"""
      CREATE TABLE merged AS SELECT CAST(column0 AS BIGINT) old, CAST(regexp_replace(column1, '[\t|]', '', 'g') AS BIGINT) new
      FROM read_csv('{ncbi}/merged.dmp', sep='\t|\t', header=false, quote='', all_varchar=true)""")

    # ---------------- structures (metadata only) ----------------
    t0 = time.time()
    con.execute("""
      CREATE TABLE s AS SELECT CAST(tier AS TINYINT) tier, CAST(unit_id AS INTEGER) unit_id,
        CAST(taxid AS BIGINT) taxid, accession, CAST(fragment AS SMALLINT) fragment, n_res, plddt_mean,
        plddt_frac_lt70, plddt_frac_lt50, status, n_pockets_all, n_pockets_kept
      FROM read_parquet(?, union_by_name=true)""", [sfiles])
    log(f"structures loaded: {con.execute('select count(*) from s').fetchone()[0]:,} rows, {time.time()-t0:.0f}s")
    # taxid -> lineage (resolve merged ids)
    # tier-3 taxa: the pipeline's own lineage table (fc_taxa) has priority; NCBI rankedlineage (after merged.dmp)
    # fills taxa not in it (tiers 1-2) and supplies the species name used for species-level collapsing.
    fct = os.path.join(atlas, "fc_taxa_afdb_uniprot_v4.parquet")
    con.execute(f"""CREATE TABLE fct AS SELECT taxid, NULLIF(tax_name,'') tax_name, NULLIF(genus,'') genus,
        NULLIF(family,'') family, NULLIF("order",'') "order", NULLIF(class,'') class, NULLIF(phylum,'') phylum,
        NULLIF(superkingdom,'') superkingdom FROM read_parquet('{fct}')""")
    con.execute("""
      CREATE MACRO resolve_lineage(tbl) AS TABLE
      WITH r AS (SELECT x.*, coalesce(m.new, x.taxid) eff FROM query_table(tbl) x LEFT JOIN merged m ON m.old = x.taxid)
      SELECT r.*, coalesce(f.tax_name, l.tax_name, 'unmapped accession') tax_name,
             coalesce(l.species, f.tax_name, l.tax_name, 'unmapped accession') species_key,
             coalesce(f.genus, l.genus) genus, coalesce(f.family, l.family) family,
             coalesce(f."order", l."order") "order", coalesce(f.class, l.class) class,
             coalesce(f.phylum, l.phylum) phylum,
             coalesce(f.superkingdom, l.superkingdom, 'unclassified') superkingdom,
             (l.taxid IS NULL) not_in_ncbi
      FROM r LEFT JOIN fct f ON f.taxid = r.taxid LEFT JOIN lin l ON l.taxid = r.eff""")
    con.execute("CREATE TABLE s_taxa AS SELECT DISTINCT taxid FROM s")
    con.execute("CREATE TABLE taxmap AS SELECT * FROM resolve_lineage(s_taxa)")

    # ---------------- 1. tier summary ----------------
    tier_sum = con.execute("""
      SELECT s.tier, count(DISTINCT unit_id) units, count(*) structures, count(DISTINCT accession) accessions,
        sum((status='ok')::INT) status_ok, sum((status='error')::INT) status_error,
        sum((status='missing')::INT) status_missing,
        sum(n_pockets_all) pockets_all, sum(n_pockets_kept) pockets_kept,
        sum(n_pockets_kept)::DOUBLE / count(*) pockets_per_structure,
        avg((n_pockets_kept>0)::INT) frac_with_kept_pocket,
        sum(n_pockets_kept)::DOUBLE / nullif(sum(n_pockets_all),0) frac_pockets_kept,
        quantile_cont(n_pockets_kept, 0.5) median_kept, quantile_cont(n_pockets_all, 0.5) median_all,
        quantile_cont(plddt_mean, [0.1,0.25,0.5,0.75,0.9]) plddt_mean_q,
        avg(plddt_mean) plddt_mean_avg,
        quantile_cont(plddt_frac_lt70, [0.1,0.5,0.9]) frac_lt70_q,
        avg((plddt_mean<70)::INT) frac_struct_plddt_lt70,
        quantile_cont(n_res, [0.1,0.5,0.9]) n_res_q,
        count(DISTINCT s.taxid) n_taxids, count(DISTINCT t.species_key) n_species
      FROM s JOIN taxmap t USING (taxid) GROUP BY s.tier ORDER BY s.tier""").df()
    for c in ["plddt_mean_q", "frac_lt70_q", "n_res_q"]:
        tier_sum[c] = tier_sum[c].apply(lambda v: "/".join(f"{x:.3g}" for x in v))
    tier_sum.to_csv(os.path.join(out, "tier_summary.csv"), index=False)
    errs = con.execute("""SELECT tier, status, count(*) n FROM s WHERE status<>'ok' GROUP BY ALL ORDER BY ALL""").df()
    errs.to_csv(os.path.join(out, "status_not_ok.csv"), index=False)
    # consistency with manifest
    mcheck = con.execute("""SELECT tier, sum(structures_rows) m_structures, sum(pockets_rows) m_pockets
                            FROM man GROUP BY tier ORDER BY tier""").df()
    log("tier summary done")

    # pocket-level quality per tier (column projection only)
    t0 = time.time()
    pq_tier = con.execute("""
      SELECT tier, count(*) pockets, count(DISTINCT pocket_id) FILTER (WHERE tier<=2) distinct_ids_t12,
        quantile_cont(burrow_score,[0.1,0.5,0.9]) burrow_score_q, quantile_cont(plddt_mean,[0.1,0.5,0.9]) plddt_q,
        quantile_cont(site_volume_A3,[0.1,0.5,0.9]) site_volume_q, avg(sph8_fallback::INT) frac_sph8_fallback,
        avg((fragment>1)::INT) frac_fragment_gt1
      FROM read_parquet(?, union_by_name=true) GROUP BY tier ORDER BY tier""", [pfiles]).df()
    for c in ["burrow_score_q", "plddt_q", "site_volume_q"]:
        pq_tier[c] = pq_tier[c].apply(lambda v: "/".join(f"{x:.3g}" for x in v))
    pq_tier.to_csv(os.path.join(out, "pocket_quality_by_tier.csv"), index=False)
    log(f"pocket quality done {time.time()-t0:.0f}s")

    # ---------------- 2. tier 1 per species ----------------
    wl = os.path.join(atlas, "worklist.parquet")
    con.execute(f"CREATE TABLE w AS SELECT * FROM read_parquet('{wl}')")
    t1 = con.execute("""
      WITH wsp AS (SELECT taxid, any_value(code) code, sum(est_structures) worklist_structures FROM w WHERE tier=1
                   GROUP BY taxid)
      SELECT s.taxid, any_value(t.tax_name) organism, any_value(wsp.code) code, any_value(t.superkingdom) superkingdom,
        count(*) structures, count(DISTINCT accession) accessions, any_value(wsp.worklist_structures) worklist_structures,
        sum(n_pockets_kept) pockets, sum(n_pockets_kept)::DOUBLE/count(*) pockets_per_structure,
        avg(plddt_mean) plddt_mean_avg, avg((plddt_mean<70)::INT) frac_struct_plddt_lt70,
        sum((status<>'ok')::INT) n_not_ok
      FROM s JOIN taxmap t USING (taxid) LEFT JOIN wsp USING (taxid)
      WHERE s.tier=1 GROUP BY s.taxid ORDER BY structures DESC""").df()
    t1.to_csv(os.path.join(out, "tier1_species.csv"), index=False)

    # ---------------- 3. tier 3 lineage + AFDB totals ----------------
    # entry-level done flag over ALL tier-3 worklist entries (fc_entries: fc_unit = worklist member_start)
    t0 = time.time()
    fce = os.path.join(atlas, "fc_entries_afdb_uniprot_v4.parquet")
    con.execute("""CREATE TABLE w3 AS SELECT unit_id, member_start fc_unit, est_structures, len(source_taxids) n_taxa,
                     (unit_id IN (SELECT unit_id FROM man)) done FROM w WHERE tier=3""")
    con.execute(f"""
      CREATE TABLE tax3 AS
      SELECT e.taxid, sum(w3.done::INT) done_entries, count(*) afdb_entries,
             sum(CASE WHEN w3.done THEN e.len ELSE 0 END) done_bytes, sum(e.len) afdb_bytes
      FROM read_parquet('{fce}') e JOIN w3 ON w3.fc_unit = e.fc_unit GROUP BY e.taxid""")
    # entry size proxy: done vs pending
    size = con.execute(f"""
      SELECT w3.done, count(*) entries, quantile_cont(e.len,[0.1,0.25,0.5,0.75,0.9]) len_bytes_q, avg(e.len) len_mean
      FROM read_parquet('{fce}') e JOIN w3 ON w3.fc_unit = e.fc_unit GROUP BY w3.done ORDER BY w3.done""").df()
    size["len_bytes_q"] = size["len_bytes_q"].apply(lambda v: "/".join(f"{x:.0f}" for x in v))
    # calibrate bytes -> residues on done entries (median n_res / len)
    calib = con.execute(f"""
      WITH d AS (SELECT accession, fragment, len FROM read_parquet('{fce}') e JOIN w3 ON w3.fc_unit=e.fc_unit
                 WHERE w3.done USING SAMPLE 200000 ROWS)
      SELECT median(s.n_res::DOUBLE / d.len) res_per_byte, corr(s.n_res, d.len) r
      FROM d JOIN s ON s.accession=d.accession AND s.fragment=d.fragment WHERE s.tier=3""").df()
    size["median_res_per_byte(done sample)"] = calib.res_per_byte[0]
    size["corr_nres_len"] = calib.r[0]
    size.to_csv(os.path.join(out, "tier3_entry_size.csv"), index=False)
    log(f"fc_entries aggregation {time.time()-t0:.0f}s")

    con.execute("CREATE TABLE tax3l AS SELECT * FROM resolve_lineage(tax3)")
    tot_done = con.execute("SELECT sum(done_entries), sum(afdb_entries) FROM tax3l").fetchone()
    # cross-check: done entries from fc_entries vs structures rows in tier 3
    s3 = con.execute("SELECT count(*) FROM s WHERE tier=3").fetchone()[0]

    rank_tables = {}
    for rank in ["superkingdom", "phylum", "class"]:
        grp = "superkingdom, phylum, class" if rank == "class" else ("superkingdom, phylum" if rank == "phylum" else "superkingdom")
        df = con.execute(f"""
          WITH st AS (SELECT {', '.join('t.' + g.strip() for g in grp.split(','))}, count(*) structures,
                             sum(n_pockets_kept) pockets, avg(plddt_mean) plddt_mean_avg,
                             count(DISTINCT s.taxid) taxa_done
                      FROM s JOIN taxmap t USING (taxid) WHERE s.tier=3 GROUP BY ALL),
          af AS (SELECT {grp}, sum(done_entries) done_entries, sum(afdb_entries) afdb_entries,
                        count(*) FILTER (WHERE done_entries>0) taxa_with_done, count(*) taxa_afdb
                 FROM tax3l GROUP BY ALL)
          SELECT coalesce(af.{rank}, '(none)') {rank}, {', '.join('af.' + g.strip() for g in grp.split(',') if g.strip() != rank)}{',' if rank != 'superkingdom' else ''}
            coalesce(st.structures,0) structures, coalesce(st.pockets,0) pockets,
            coalesce(st.pockets,0)::DOUBLE/nullif(st.structures,0) pockets_per_structure, st.plddt_mean_avg,
            af.taxa_with_done, af.taxa_afdb, af.done_entries, af.afdb_entries,
            af.done_entries::DOUBLE/af.afdb_entries frac_done,
            af.done_entries::DOUBLE/{tot_done[0]} share_done, af.afdb_entries::DOUBLE/{tot_done[1]} share_afdb,
            (af.done_entries::DOUBLE/{tot_done[0]}) / (af.afdb_entries::DOUBLE/{tot_done[1]}) enrichment
          FROM af LEFT JOIN st ON {' AND '.join(f'st.{g.strip()} IS NOT DISTINCT FROM af.{g.strip()}' for g in grp.split(','))}
          ORDER BY af.afdb_entries DESC""").df()
        df.to_csv(os.path.join(out, f"tier3_{rank}.csv"), index=False)
        rank_tables[rank] = df
    log("lineage tables done")

    # taxa thresholds (done proteins per taxon) at taxid and species level
    thr = []
    for N in [1, 10, 100, 1000]:
        r = con.execute(f"""SELECT
            (SELECT count(*) FROM tax3l WHERE done_entries>={N}) taxids_done,
            (SELECT count(*) FROM (SELECT species_key FROM tax3l GROUP BY 1 HAVING sum(done_entries)>={N})) species_done,
            (SELECT count(*) FROM tax3l WHERE afdb_entries>={N}) taxids_afdb,
            (SELECT count(*) FROM (SELECT species_key FROM tax3l GROUP BY 1 HAVING sum(afdb_entries)>={N})) species_afdb
        """).df()
        r.insert(0, "N", N)
        thr.append(r)
    thr = pd.concat(thr)
    thr.to_csv(os.path.join(out, "tier3_taxa_thresholds.csv"), index=False)

    # completeness: per taxid, fraction of its tier-3 AFDB entries present (tiers 1-2 entries counted as present)
    con.execute("""CREATE TABLE t12tax AS SELECT taxid, count(*) n12 FROM s WHERE tier<=2 GROUP BY taxid""")
    comp = con.execute("""
      WITH x AS (SELECT t.taxid, t.done_entries + coalesce(a.n12,0) present, t.afdb_entries + coalesce(a.n12,0) total
                 FROM tax3l t LEFT JOIN t12tax a USING (taxid)),
      b AS (SELECT *, CASE WHEN total<10 THEN 'a: 1-9' WHEN total<100 THEN 'b: 10-99' WHEN total<1000 THEN 'c: 100-999'
                           WHEN total<5000 THEN 'd: 1k-5k' ELSE 'e: >=5k' END size_bin,
                      present::DOUBLE/total f FROM x)
      SELECT size_bin, count(*) taxa_afdb, sum(total) entries_afdb,
        count(*) FILTER (WHERE f=0) f_0, count(*) FILTER (WHERE f>0 AND f<0.5) f_lt50,
        count(*) FILTER (WHERE f>=0.5 AND f<0.9) f_50_90, count(*) FILTER (WHERE f>=0.9 AND f<1) f_90_100,
        count(*) FILTER (WHERE f=1) f_complete, sum(present) entries_present
      FROM b GROUP BY size_bin ORDER BY size_bin""").df()
    comp.to_csv(os.path.join(out, "tier3_completeness.csv"), index=False)
    # list of (near-)complete large taxa
    big_complete = con.execute("""
      SELECT t.taxid, t.tax_name, t.superkingdom, t.phylum, t.class, t.done_entries, t.afdb_entries,
             coalesce(a.n12,0) tier12_entries,
             (t.done_entries+coalesce(a.n12,0))::DOUBLE/(t.afdb_entries+coalesce(a.n12,0)) frac_present
      FROM tax3l t LEFT JOIN t12tax a USING (taxid)
      WHERE t.afdb_entries+coalesce(a.n12,0) >= 1000
        AND (t.done_entries+coalesce(a.n12,0))::DOUBLE/(t.afdb_entries+coalesce(a.n12,0)) >= 0.9
      ORDER BY t.afdb_entries DESC""").df()
    big_complete.to_csv(os.path.join(out, "tier3_nearcomplete_taxa_ge1000.csv"), index=False)
    # split taxa: taxa whose entries span >1 unit, and done fraction for them
    split = con.execute(f"""
      WITH tu AS (SELECT e.taxid, count(DISTINCT e.fc_unit) n_units FROM read_parquet('{fce}') e GROUP BY e.taxid)
      SELECT CASE WHEN n_units=1 THEN '1 unit' WHEN n_units<=10 THEN '2-10 units' ELSE '>10 units' END units_per_taxon,
             count(*) taxa, sum(t.afdb_entries) entries, sum(t.done_entries) done_entries,
             count(*) FILTER (WHERE t.done_entries>0 AND t.done_entries<t.afdb_entries) taxa_partially_done,
             count(*) FILTER (WHERE t.done_entries=t.afdb_entries) taxa_fully_done
      FROM tu JOIN tax3l t USING (taxid) GROUP BY 1 ORDER BY 1""").df()
    split.to_csv(os.path.join(out, "tier3_taxon_splitting.csv"), index=False)
    log("completeness done")

    # unclassified / metagenomic
    unc = con.execute(f"""
      WITH c AS (SELECT *, CASE
          WHEN taxid=0 THEN 'unmapped (taxid 0)'
          WHEN not_in_ncbi THEN 'taxid not in NCBI dump'
          WHEN lower(tax_name) LIKE '%metagenome%' THEN 'metagenome'
          WHEN superkingdom = 'unclassified' THEN 'superkingdom unclassified (other)'
          WHEN {PLACEHOLDER_SQL} THEN 'placeholder name (MAG-like: "X bacterium", "uncultured", "sp.")'
          ELSE 'classified, named' END cat FROM tax3l)
      SELECT cat, count(*) FILTER (WHERE done_entries>0) taxa_done, sum(done_entries) done_entries,
             sum(done_entries)::DOUBLE/{tot_done[0]} frac_done_entries,
             count(*) taxa_afdb, sum(afdb_entries) afdb_entries, sum(afdb_entries)::DOUBLE/{tot_done[1]} frac_afdb_entries
      FROM c GROUP BY cat ORDER BY afdb_entries DESC""").df()
    unc.to_csv(os.path.join(out, "tier3_unclassified.csv"), index=False)

    # unit-level done vs pending
    units = con.execute("""
      SELECT done, count(*) units, sum(est_structures) structures, avg(est_structures) mean_unit_size,
        avg((est_structures=4000)::INT) frac_full_units, avg(n_taxa) mean_taxa_per_unit,
        quantile_cont(n_taxa,0.5) median_taxa_per_unit, min(unit_id) min_unit, max(unit_id) max_unit
      FROM w3 GROUP BY done ORDER BY done""").df()
    units.to_csv(os.path.join(out, "tier3_done_vs_pending_units.csv"), index=False)
    first_k = con.execute("""
      WITH r AS (SELECT unit_id, done, row_number() OVER (ORDER BY unit_id) rk FROM w3)
      SELECT max(rk) FILTER (WHERE done) last_done_rank,
             min(rk) FILTER (WHERE NOT done) first_pending_rank,
             count(*) FILTER (WHERE done) n_done FROM r""").df()

    # quality by superkingdom (tier 3) and all tiers
    qsk = con.execute("""
      SELECT s.tier, coalesce(t.superkingdom,'(none)') superkingdom, count(*) structures,
        avg(plddt_mean) plddt_mean_avg, quantile_cont(plddt_mean,0.5) plddt_mean_median,
        avg(plddt_frac_lt70) frac_lt70_avg, sum(n_pockets_kept)::DOUBLE/count(*) pockets_per_structure,
        sum(n_pockets_kept)::DOUBLE/nullif(sum(n_pockets_all),0) frac_pockets_kept, avg(n_res) n_res_avg
      FROM s JOIN taxmap t USING (taxid) GROUP BY ALL ORDER BY s.tier, structures DESC""").df()
    qsk.to_csv(os.path.join(out, "quality_by_superkingdom.csv"), index=False)
    log("all tables done; writing markdown")

    # ---------------- markdown ----------------
    sk = rank_tables["superkingdom"]
    ph = rank_tables["phylum"].copy()
    ph_big = ph[ph.afdb_entries >= 200_000].copy()
    over = ph_big.sort_values("enrichment", ascending=False).head(8)
    under = ph_big.sort_values("enrichment").head(8)
    cls = rank_tables["class"]
    n_phyla_done = int((ph.done_entries > 0).sum())
    n_class_done = int((cls.done_entries > 0).sum())
    meta_row = unc.set_index("cat")
    def ufrac(cat, col):
        return float(meta_row.loc[cat, col]) if cat in meta_row.index else 0.0

    L = []
    L.append(f"# Atlas coverage — snapshot `{snap_id}`\n")
    L.append(f"Generated by `scripts/coverage_report.py` from `{os.path.relpath(man_path, P)}` "
             f"(manifest sha256 `{summary['manifest_sha256'][:16]}…`, atlas `{atlas}`, worklist sha256 "
             f"`{(summary.get('worklist_sha256') or '')[:16]}…`). Snapshot created {summary['created_utc']}. "
             "All numbers refer only to the units in this manifest; the atlas was still being written. "
             "Treat every tier-3 statement as preliminary.\n")
    L.append("## 1. Units, structures, pockets per tier\n")
    pt = summary["per_tier"]
    ut = pd.DataFrame([{"tier": int(t), "units_done": v.get("units", 0), "worklist_units": v.get("worklist_units"),
                        "frac_units_done": v.get("units", 0) / v["worklist_units"] if v.get("worklist_units") else None}
                       for t, v in sorted(pt.items())])
    L.append(md_table(ut))
    L.append("")
    cols = ["tier", "structures", "accessions", "pockets_kept", "pockets_all", "pockets_per_structure",
            "frac_with_kept_pocket", "frac_pockets_kept", "n_taxids", "n_species"]
    L.append(md_table(tier_sum[cols]))
    L.append("\n`pockets_all` = merged BURROW pockets before the score (≥ 0.03) / lining-pLDDT (≥ 70) filters; "
             "`frac_pockets_kept` = kept/all. Manifest row counts agree with the structure tables: "
             + ", ".join(f"tier {int(r.tier)} {int(r.m_structures):,} structures / {int(r.m_pockets):,} pockets"
                         for r in mcheck.itertuples()) + ".\n")
    L.append("### Structure quality per tier\n")
    L.append(md_table(tier_sum[["tier", "status_ok", "status_error", "status_missing", "plddt_mean_avg",
                                "plddt_mean_q", "frac_struct_plddt_lt70", "frac_lt70_q", "n_res_q",
                                "median_all", "median_kept"]]))
    L.append("\nQuantiles: pLDDT mean 10/25/50/75/90 %; frac_lt70 and n_res 10/50/90 %. "
             "Tier-3 pLDDT comes from CA B-factors of Foldcomp-decompressed v4 models (±0.03).\n")
    if len(errs):
        L.append("Non-ok structures:\n\n" + md_table(errs) + "\n")
    else:
        L.append("No structures with status ≠ ok in this snapshot.\n")
    L.append("### Pocket-level properties per tier\n")
    L.append(md_table(pq_tier))
    L.append("\n(quantiles 10/50/90 %; `distinct_ids_t12` is computed for tiers 1–2 only.)\n")
    L.append("### Tier-3 quality by superkingdom\n")
    L.append(md_table(qsk))
    L.append("")

    L.append("## 2. Tier 1: per-species counts\n")
    L.append(md_table(t1[["taxid", "organism", "code", "superkingdom", "structures", "accessions",
                          "worklist_structures", "pockets", "pockets_per_structure", "plddt_mean_avg",
                          "frac_struct_plddt_lt70"]]))
    L.append("")

    L.append("## 3. Tier 3: taxonomic breakdown\n")
    L.append(f"Tier-3 entries in this snapshot: {tot_done[0]:,} of {tot_done[1]:,} AFDB-v4 tier-3 entries "
             f"({tot_done[0]/tot_done[1]:.1%}); tier-3 structure rows read: {s3:,} "
             f"({'consistent' if s3 == tot_done[0] else 'MISMATCH'} with fc_entries). "
             f"Phyla with ≥1 done entry: {n_phyla_done} of {int((ph.afdb_entries>0).sum())}; "
             f"classes: {n_class_done} of {int((cls.afdb_entries>0).sum())}. "
             "`enrichment` = share of done entries / share of all tier-3 AFDB entries "
             "(1 = proportional; lineages from NCBI rankedlineage.dmp after merged.dmp).\n")
    L.append("### Superkingdom\n")
    L.append(md_table(sk[["superkingdom", "structures", "pockets", "pockets_per_structure", "plddt_mean_avg",
                          "taxa_with_done", "taxa_afdb", "frac_done", "share_done", "share_afdb", "enrichment"]]))
    L.append("\n### Phylum (top 25 by AFDB size; full table `tier3_phylum.csv`)\n")
    L.append(md_table(ph[["phylum", "superkingdom", "structures", "pockets", "taxa_with_done", "taxa_afdb",
                          "frac_done", "share_done", "share_afdb", "enrichment"]], max_rows=25))
    L.append("\n### Class (top 25 by AFDB size; full table `tier3_class.csv`)\n")
    L.append(md_table(cls[["class", "phylum", "structures", "taxa_with_done", "taxa_afdb", "frac_done",
                           "share_done", "share_afdb", "enrichment"]], max_rows=25))
    L.append("\n### Tier-3 taxa with ≥ N done proteins\n")
    L.append(md_table(thr))
    L.append("\n`taxids_*` counts NCBI taxids as stored (strains count separately); `species_*` collapses to the "
             "NCBI species name (or the taxon name for taxa at/above species). `*_afdb` = all of tier-3 AFDB v4.\n")
    L.append("### Unmapped / unclassified / metagenomic\n")
    L.append(md_table(unc))
    L.append("")

    L.append("## 4. Done vs pending (is the done subset biased?)\n")
    L.append("Unit level:\n\n" + md_table(units))
    L.append(f"\nDone units are {'a contiguous prefix' if int(first_k.last_done_rank[0]) < int(first_k.first_pending_rank[0]) else 'NOT a contiguous prefix'} "
             f"of the round-robin worklist order (last done rank {int(first_k.last_done_rank[0]):,}, first pending rank "
             f"{int(first_k.first_pending_rank[0]):,}, {int(first_k.n_done[0]):,} done).\n")
    L.append("Entry size (Foldcomp entry bytes, proxy for protein length; median residues/byte and Pearson r with "
             "n_res measured on a 200k sample of done entries):\n\n" + md_table(size))
    L.append("\nMost over-represented phyla among done entries (phyla with ≥ 200k AFDB entries):\n\n"
             + md_table(over[["phylum", "superkingdom", "done_entries", "afdb_entries", "frac_done", "enrichment"]]))
    L.append("\nMost under-represented:\n\n"
             + md_table(under[["phylum", "superkingdom", "done_entries", "afdb_entries", "frac_done", "enrichment"]]))
    L.append("")

    L.append("## 5. Per-species proteome completeness (tier 3; matters for orthology)\n")
    L.append("A tier-3 work unit holds ≤ 4,000 Foldcomp entries: large taxa are split over many units, and many "
             "small taxa are packed into one unit (`source_taxids`). In round-robin order a large taxon therefore "
             "receives one slice per round, so early snapshots contain *slices* of big proteomes. "
             "Below, for every tier-3 taxid, `present` = its done tier-3 entries + its tier-1/2 entries, "
             "`total` = all its tier-3 AFDB-v4 entries + its tier-1/2 entries (AFDB entries ≈ UniProt entries for "
             "the taxon, which is not the same as a reference proteome).\n")
    L.append(md_table(comp))
    L.append(f"\nTaxa with ≥ 1,000 AFDB entries that are ≥ 90 % present: **{len(big_complete):,}** "
             "(`tier3_nearcomplete_taxa_ge1000.csv`).\n")
    if len(big_complete):
        L.append(md_table(big_complete.head(20)))
        L.append("")
    L.append("Taxon splitting across units:\n\n" + md_table(split))
    L.append("")

    # ---------------- biases & caveats ----------------
    def sk_enrich(name):
        r = sk[sk.superkingdom == name]
        return (float(r.share_done.iloc[0]), float(r.share_afdb.iloc[0]), float(r.enrichment.iloc[0])) if len(r) else (0, 0, float("nan"))
    L.append("## 6. Biases & caveats\n")
    b = []
    for name in ["Bacteria", "Eukaryota", "Archaea", "Viruses", "unclassified"]:
        sd, sa, e = sk_enrich(name)
        b.append(f"{name}: {sd:.1%} of done vs {sa:.1%} of all tier-3 entries (×{e:.2f})")
    L.append("- **Superkingdom composition of done tier-3 entries vs all tier 3:** " + "; ".join(b) + ".")
    exhausted = sk[sk.frac_done >= 0.99]
    if len(exhausted):
        L.append("- **Whole superkingdoms already exhausted:** " + ", ".join(
            f"{r.superkingdom} ({int(r.done_entries):,}/{int(r.afdb_entries):,} entries done)" for r in exhausted.itertuples())
            + ". The round-robin interleaves superkingdoms, so small superkingdoms finish first and later units are "
              "almost all Bacteria/Eukaryota: the composition of the done subset will keep shifting until the run ends; "
              "any done/pending comparison is a snapshot property, not a stable one. "
            + "Remaining fraction done: " + ", ".join(
            f"{r.superkingdom} {r.frac_done:.1%}" for r in sk[sk.frac_done < 0.99].itertuples()) + ".")
    L.append(f"- **Round-robin over the tree flattens clade sizes:** the done subset is deliberately broad, not "
             f"proportional. Over-represented phyla (≥200k entries): "
             + ", ".join(f"{r.phylum} ×{r.enrichment:.1f}" for r in over.head(5).itertuples())
             + "; under-represented: "
             + ", ".join(f"{r.phylum} ×{r.enrichment:.2f}" for r in under.head(5).itertuples())
             + ". Small/rare clades are over-sampled relative to the giant bacterial phyla. Any statistic pooled "
               "over tier 3 must be weighted (per taxon / per clade) or compared against the matching pending set.")
    unit_done = units[units.done == True]  # noqa: E712
    unit_pend = units[units.done == False]  # noqa: E712
    if len(unit_done) and len(unit_pend):
        L.append(f"- **Unit composition:** done units average {float(unit_done.mean_unit_size.iloc[0]):,.0f} structures and "
                 f"{float(unit_done.mean_taxa_per_unit.iloc[0]):.1f} taxa/unit vs {float(unit_pend.mean_unit_size.iloc[0]):,.0f} "
                 f"and {float(unit_pend.mean_taxa_per_unit.iloc[0]):.1f} for pending units.")
    if len(size) == 2:
        L.append(f"- **Protein size:** median Foldcomp entry size done {size.len_bytes_q.iloc[1].split('/')[2]} B vs pending "
                 f"{size.len_bytes_q.iloc[0].split('/')[2]} B (≈ {calib.res_per_byte[0]:.3f} residues/byte, r = {calib.r[0]:.2f}).")
    L.append(f"- **Partial proteomes:** of tier-3 taxa with ≥1 done entry, most are represented by a slice; only "
             f"{len(big_complete):,} taxa with ≥1,000 entries are ≥90 % present. Presence/absence (gain/loss) calls "
             "across tier-3 species are not interpretable until the species' full entry set is done — restrict such "
             "analyses to tier 1 (complete reference proteomes) + tier-3 taxa flagged complete, or condition on "
             "`frac_present`.")
    L.append("- **Slices are contiguous, not random:** within a taxon, units are consecutive runs of Foldcomp-DB "
             "entries (DB offset order, which follows accession/key order and hence often one UniProt proteome "
             "submission). A done slice of a species is therefore not a random sample of its proteome.")
    L.append(f"- **Unclassified / metagenomic:** in done tier 3, {ufrac('metagenome','frac_done_entries'):.1%} of entries are "
             f"from metagenome taxa (vs {ufrac('metagenome','frac_afdb_entries'):.1%} of all tier 3), "
             f"{ufrac('superkingdom unclassified (other)','frac_done_entries'):.1%} are other non-metagenome unclassified taxa, "
             f"{ufrac('unmapped (taxid 0)','frac_done_entries'):.2%} are unmapped (taxid 0), and "
             f"{ufrac('placeholder name (MAG-like: \"X bacterium\", \"uncultured\", \"sp.\")','frac_done_entries'):.1%} carry "
             "placeholder (MAG-like) names whose lineage stops at phylum/class.")
    L.append("- **Tiers are different populations:** tier 1 = 46 complete reference proteomes (model organisms and "
             "pathogens, eukaryote-heavy), tier 2 = Swiss-Prot (curated, well-studied, higher pLDDT), tier 3 = TrEMBL "
             "remainder. Tier 1/2 are v6 models from AFDB files; tier 3 are v4 Foldcomp-decompressed models — model "
             "version/source is a covariate (POCKETOME §6.3, §6.10); never compare the same protein across sources.")
    L.append("- **Accession dedup across tiers:** each accession appears in its first tier only, so a tier-3 taxon "
             "that also has Swiss-Prot entries has those proteins in tier 2 (counted as present in §5).")
    L.append("- **Filters shape content:** pockets with BURROW score < 0.03 or lining pLDDT < 70 are dropped and "
             "pLDDT < 50 residues pruned; low-confidence/disordered proteins and clades with lower pLDDT have fewer "
             "kept pockets per structure (see quality tables). `frac_pockets_kept` differs by tier and superkingdom.")
    L.append("- **`source` column vs the guide:** tier-1 structure tables store the proteome tar name "
             "(`UP…_v6.tar`) in `source`, not `afdb_v6_tar` as POCKETOME §3 says; tier-1 pocket tables lack the "
             "column. Normalise by tier, never by string equality on `source`.")
    L.append("- **Empty units:** units with 0 kept pockets write a pocket file without a `source` column even in "
             "tier 3; read with `union_by_name=true` and fill `source` from the tier.")
    L.append("- **Fragments:** long human proteins are split into overlapping fragments in tier 1 (structures > "
             "accessions); deduplicate by residue overlap before per-protein counts.")
    L.append("- **Strain-level taxids:** many tier-3 taxids are strains/isolates of the same species; species-level "
             "counts collapse them by NCBI species name.")
    with open(os.path.join(out, "COVERAGE.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    json.dump({"snapshot_id": snap_id, "manifest": man_path, "manifest_sha256": summary["manifest_sha256"],
               "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z")},
              open(os.path.join(out, "provenance.json"), "w"), indent=1)
    log(f"wrote {out}")
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
