import sys, time
from ancestra.orthoxml import read_orthoxml, hogs_to_pairs, hogs_at_level, Group
from ancestra.newick import *  # noqa
ox_path, tree_path = sys.argv[1], sys.argv[2]
t=time.time(); ox = read_orthoxml(ox_path); print('read', round(time.time()-t,2),'s', len(ox.genes),'genes', len(ox.groups),'top groups')
# 1. orthogroup (root HOG) of each protein
root_of = {}
for i, g in enumerate(ox.groups):
    if isinstance(g, Group):
        for p in g.iter_genes(): root_of[p] = i
# 2. pair relations
orth, para = hogs_to_pairs(ox, with_paralogs=True)
a, b = next(iter(orth)); print('ortholog example', a, b, ox.species_of(a), ox.species_of(b))
a, b = next(iter(para)); print('paralog example', a, b)
# 3. labelled gene tree of one HOG: nested Newick with D/S labels and level
def to_nwk(x):
    if not isinstance(x, Group): return x
    lab = ('S_' if x.kind=='ortholog' else 'D_') + str(x.level)
    return '(' + ','.join(to_nwk(c) for c in x.children) + ')' + lab.replace(' ','_')
big = max((g for g in ox.groups if isinstance(g, Group) and any(h.kind=='paralog' for h in g.iter_groups())), key=lambda g: len(g.genes()))
print('HOG as labelled tree:', to_nwk(big)[:300])
# 4. HOGs at an internal level
lv = [n.name for n in ox.species_tree.preorder() if n.children][1] if ox.species_tree else None
print('level', lv, 'n HOGs', len(hogs_at_level(ox, lv)))
# 5. pyham: per-branch gains / losses / duplications
import pyham
ham = pyham.Ham(tree_path, ox_path, use_internal_name=True)
anc = ham.get_ancestral_genome_by_name(lv)
print('pyham ancestral genome', lv, 'n HOGs', len(anc.genes))
ext = ham.get_extant_genome_by_name(ham.get_list_extant_genomes()[0].name)
vh = ham.compare_genomes_vertically(ext, anc)
print('vertical', ext.name, 'vs', lv, 'lost', len(vh.get_lost()), 'gained', len(vh.get_gained()), 'duplicated', len(vh.get_duplicated()))
