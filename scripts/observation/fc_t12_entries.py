"""Foldcomp afdb_uniprot_v4 offsets (.lookup + .index) of all tier-1/2 accessions -> /tmp/pe_big/obs/fc_t12_entries.parquet (17 s)."""
import duckdb, time
t=time.time()
F='/cluster/scratch/gmagnusson/pocketcompare_afdb'
con=duckdb.connect(); con.execute("SET threads=6"); con.execute("SET memory_limit='40GB'")
con.execute(f"""CREATE TABLE skip AS SELECT * FROM read_csv('{F}/skip_accessions.txt', delim='\t', header=false, columns={{'acc':'VARCHAR','tier':'INTEGER'}})""")
con.execute(f"""CREATE TABLE lk AS SELECT l.id, l.name, l.tarid, s.acc, s.tier FROM read_csv('{F}/foldcomp/afdb_uniprot_v4.lookup', delim='\t', header=false, columns={{'id':'UBIGINT','name':'VARCHAR','tarid':'UBIGINT'}}) l
  JOIN skip s ON split_part(l.name,'-',2)=s.acc""")
print(con.execute("select count(*), count(distinct acc) from lk").fetchall(), time.time()-t)
con.execute(f"""COPY (SELECT lk.*, i.off, i.len FROM lk JOIN read_csv('{F}/foldcomp/afdb_uniprot_v4.index', delim='\t', header=false, columns={{'id':'UBIGINT','off':'UBIGINT','len':'UBIGINT'}}) i USING (id)) TO '/tmp/pe_big/obs/fc_t12_entries.parquet' (FORMAT parquet)""")
print(time.time()-t)
