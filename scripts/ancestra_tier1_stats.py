import sys, collections, pandas as pd
from ancestra.orthoxml import read_orthoxml, hogs_to_pairs, Group
ox=read_orthoxml(sys.argv[1])
sp={k:v.species for k,v in ox.genes.items()}
ingrp=set(); 
for g in ox.groups:
    if isinstance(g,Group): ingrp.update(g.iter_genes())
c=collections.Counter(sp[x] for x in ingrp); tot=collections.Counter(sp.values())
cov=pd.Series({s:c[s]/tot[s] for s in tot}).sort_values()
print('genes in multi-gene HOGs: %d/%d'%(len(ingrp),len(sp)),'per-species coverage min/median/max',cov.min().round(3),cov.median().round(3),cov.max().round(3)); print(cov.head(5).round(3).to_string())
orth,para=hogs_to_pairs(ox,with_paralogs=True)
H,M='UP000005640_9606','UP000000589_10090'
hm=[p for p in orth if {sp[p[0]],sp[p[1]]}=={H,M}]
deg=collections.Counter(); 
for a,b in hm: deg[a]+=1; deg[b]+=1
one=sum(1 for a,b in hm if deg[a]==1 and deg[b]==1)
print('human-mouse ortholog pairs',len(hm),'1:1',one)
for a,b in [('P07900','P02829'),('P07900','P08238'),('P00533','P04626'),('P00533','Q01279'),('P68871','P02088')]:
    p=(min(a,b),max(a,b)); print(a,b,'ortholog' if p in orth else 'paralog' if p in para else 'unrelated/none')
t=pd.read_csv('tier1_acc_taxid.tsv',sep='\t'); t['in_run']=t.accession.isin(sp.keys()); t['in_hog']=t.accession.isin(ingrp)
print('atlas tier1 accessions in run: %.3f, in multi-gene HOG: %.3f'%(t.in_run.mean(),t.in_hog.mean()))
q=pd.read_csv('qfo_genes.tsv',sep='\t'); 
print('ok')
