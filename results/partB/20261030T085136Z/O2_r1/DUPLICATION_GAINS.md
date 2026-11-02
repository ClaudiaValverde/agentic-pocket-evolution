# Duplication-associated gains — O2_r1, snapshot 20261030T085136Z (EXPLORATORY, NOT PREREGISTERED)

Generated 2026-10-30 16:06 UTC by `scripts/partB/dup_gain_test.py`; git `5cd4fd2f967a3673db5ac7fead77b9237bf64853+dirty`. Scope: primary gains on a paralog branch that fail `orthology_split` (55 of 104 primary gains).

Test: species with the pocket (P) in the derived paralog copy and an observed absence (A only; pLDDT >= 90, D11, after rescue) in the sister copy; category if >= 3 such species and more than under the label-swap null (sister P, derived A only). No FDR control; labelled "duplication-associated gain (exploratory, not preregistered)".

**18 of 55 pass.** Label-swap null: total pattern species 129 vs swapped 1; SOGs with >= 3 swapped species: 0.

Background null (no selection): all other duplication-node copies of the same SOG trees, both orientations: 5276 copies, 4 (0.001) reach >= 3 pattern species. The label-swap null is conditioned on the gain selection (derived copy P, outgroups A) and is therefore weak; read the background rate instead.

| sog_id | level | duplication_level | n_species_both | n_pattern | n_swapped | dup-associated gain | pass_artifacts |
|---|---|---|---|---|---|---|---|
| HOG0054963:g20 | Eumetazoa | Eukaryota | 17 | 12 | 0 | True | False |
| HOG0050325:g11 | Bilateria | Bilateria | 32 | 10 | 0 | True | False |
| HOG0054257:g2 | Opisthokonta | Eukaryota | 19 | 9 | 1 | True | False |
| HOG0008302:g0 | Pseudomonadati | Pseudomonadati | 9 | 8 | 0 | True | False |
| HOG0056321:g9 | Embryophyta | Embryophyta | 8 | 6 | 0 | True | False |
| HOG0046679:g11 | Gnathostomata | Gnathostomata | 11 | 5 | 0 | True | False |
| HOG0057094:g6 | Thermoproteati | Thermoproteati | 10 | 5 | 0 | True | False |
| HOG0054166:g12 | Embryophyta | Embryophyta | 9 | 5 | 0 | True | False |
| HOG0054071:g9 | Rhabditida | Nematoda | 7 | 4 | 0 | True | False |
| HOG0050569:g8 | leotiomyceta | leotiomyceta | 9 | 4 | 0 | True | False |
| HOG0054282:g25 | Bilateria | Bilateria | 32 | 4 | 0 | True | False |
| HOG0046258:g14 | Gnathostomata | Gnathostomata | 12 | 3 | 0 | True | False |
| HOG0050209:g8 | Rhabditida | Eumetazoa | 7 | 3 | 0 | True | False |
| HOG0051106:g12 | Viridiplantae | Viridiplantae | 10 | 3 | 0 | True | False |
| HOG0036074:g7 | Euarchontoglires | Tetrapoda | 4 | 3 | 0 | True | False |
| HOG0048923:g22 | Spirurina | Ecdysozoa | 3 | 3 | 0 | True | False |
| HOG0002084:g4 | Sulfolobaceae | Thermoprotei | 3 | 3 | 0 | True | False |
| HOG0043090:g10 | Eumetazoa | Eumetazoa | 13 | 3 | 0 | True | False |
| HOG0042176:g10 | leotiomyceta | saccharomyceta | 3 | 2 | 0 | False | False |
| HOG0044760:g1 | Rhabditida | Rhabditida | 5 | 2 | 0 | False | False |
| HOG0045727:g8 | Mesangiospermae | Embryophyta | 5 | 2 | 0 | False | False |
| HOG0005220:g5 | Enterobacteriaceae | Enterobacteriaceae | 4 | 2 | 0 | False | False |
| HOG0008276:g3 | Gammaproteobacteria | Pseudomonadati | 6 | 2 | 0 | False | False |
| HOG0053605:g5 | Embryophyta | Embryophyta | 9 | 2 | 0 | False | False |
| HOG0050360:g23 | Rhabditida | Bilateria | 7 | 2 | 0 | False | False |
| HOG0050872:g14 | Nematoda | Nematoda | 6 | 2 | 0 | False | False |
| HOG0052311:g7 | Mesangiospermae | Embryophyta | 6 | 2 | 0 | False | False |
| HOG0056718:g5 | Pseudomonadati | Pseudomonadati | 1 | 1 | 0 | False | False |
| HOG0056709:g10 | Protostomia | Bilateria | 3 | 1 | 0 | False | False |
| HOG0045534:g15 | Eurotiomycetes | Eurotiomycetes | 3 | 1 | 0 | False | False |
| HOG0048636:g12 | Rhabditida | Rhabditida | 6 | 1 | 0 | False | False |
| HOG0045692:g8 | Boreoeutheria | Boreoeutheria | 2 | 1 | 0 | False | False |
| HOG0006772:g2 | Gammaproteobacteria | Gammaproteobacteria | 2 | 1 | 0 | False | False |
| HOG0046728:g6 | saccharomyceta | saccharomyceta | 8 | 1 | 0 | False | False |
| HOG0047083:g6 | Eumetazoa | Eumetazoa | 1 | 1 | 0 | False | False |
| HOG0046638:g5 | leotiomyceta | saccharomyceta | 5 | 1 | 0 | False | False |
| HOG0042807:g13 | Nematoda | Nematoda | 8 | 1 | 0 | False | False |
| HOG0050872:g9 | Nematoda | Nematoda | 6 | 1 | 0 | False | False |
| HOG0053601:g6 | Thermoplasmatota | Thermoplasmatota | 1 | 1 | 0 | False | False |
| HOG0051112:g18 | Mesangiospermae | Embryophyta | 7 | 1 | 0 | False | False |
| HOG0056107:g8 | Pseudomonadati | Pseudomonadati | 4 | 1 | 0 | False | False |
| HOG0051923:g8 | Gammaproteobacteria | Gammaproteobacteria | 1 | 1 | 0 | False | False |
| HOG0053806:g5 | leotiomyceta | Ascomycota | 8 | 1 | 0 | False | False |
| HOG0056107:g21 | Pseudomonadati | Pseudomonadati | 4 | 1 | 0 | False | False |
| HOG0056094:g11 | Pseudomonadati | Bacteria | 9 | 1 | 0 | False | False |
| HOG0046486:g12 | Mesangiospermae | Mesangiospermae | 2 | 0 | 0 | False | False |
| HOG0042720:g1 | Amniota | Eumetazoa | 0 | 0 | 0 | False | False |
| HOG0053895:g12 | Ecdysozoa | Protostomia | 3 | 0 | 0 | False | False |
| HOG0051221:g20 | leotiomyceta | Ascomycota | 5 | 0 | 0 | False | False |
| HOG0055139:g9 | Fungi | Fungi | 6 | 0 | 0 | False | False |
| HOG0054501:g10 | Rhabditida | Rhabditida | 4 | 0 | 0 | False | False |
| HOG0055724:g2 | Bilateria | Eukaryota | 0 | 0 | 0 | False | False |
| HOG0055628:g19 | leotiomyceta | leotiomyceta | 5 | 0 | 0 | False | False |
| HOG0056203:g8 | Eukaryota | Eukaryota | 1 | 0 | 0 | False | False |
| HOG0056733:g9 | Eumetazoa | Eukaryota | 0 | 0 | 0 | False | False |
