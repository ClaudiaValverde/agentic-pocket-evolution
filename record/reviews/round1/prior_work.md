# Round 1 review: prior work, novelty, convergence definitions, positive controls
Reviewer lens: binding-site / enzyme-function evolution and active-site convergence. Reviewed DESIGN.md v0, LITERATURE.md, PROMPT.md, 2026-10-30.
Method: references checked on the web or Crossref. Every control accession below was looked up in the snapshot atlas
(`pockets/tier{1,2}`, script in scratchpad `chk.py`/`cov.py`). For each one I checked whether the known UniProt
active-site or binding residues fall inside a BURROW pocket. These are single-protein spot checks, not a benchmark.

## Verdict
The engagement with prior work is broadly good: Babbitt/Gerlt, Furnham, Riziotis, Gao & Skolnick and Masters are all used
correctly. The novelty claim mostly holds. A search found no study that maps individual pockets onto gene and species
trees with a gain/loss model. The closest work is Zillmer & Walther 2025, which works at species level. Two things in the
design need to change. (a) Its convergence vocabulary does not match how the field uses the words. (b) Its main
convergence positive set (Riziotis catalytic sites) is largely invisible to BURROW on monomer models. I verified (b) in the atlas.

## Numbered changes
1. **Use pattern-based terms the way the field uses them** (Stayton 2015 Evolution doi:10.1111/evo.12729; Storz 2016 Nat Rev
   Genet doi:10.1038/nrg.2016.11; Doolittle 1994). Define three labels in §4.2/§5.2:
   - *analogous similarity*: similar sites in N pairs, with no claim about direction;
   - *convergence*: independently evolved similarity, meaning the derived states are closer than their ancestors (the Δ test);
   - *parallelism/homoplasy within a family*: the same SOG gained, or the same lining substitution made, more than once
     from the same ancestral state inside one HOG tree.
   Reversal (loss then regain) is homoplasy too. Rationale: F_conv as currently defined measures analogy, not convergence.
   Calling it convergence would contradict Gao & Skolnick and Stayton's warning that stochastic processes alone produce
   apparent convergence.
2. **Rename and split the Part A metric.** Rename F_conv to **F_analog** (excess similarity mass in N pairs). Add
   **F_conv,dir**, the fraction of N-pair excess mass on pairs with Δ>0 that is significant after the symmetric-swap null.
   Report only F_conv,dir as "convergence". Also run the Δ test on H2/H3 pairs. Riziotis found 9–15% of
   function-matched pairs are homologues, and apicomplexan LDH (below) shows convergence can happen inside a superfamily.
   Rationale: the project goal asks "what fraction points to convergence". Excess similarity alone cannot answer that.
3. **Strengthen the Gao–Skolnick null.** A null matched only on volume decile is too weak. The prediction from
   Gao & Skolnick 2013 is that non-homologous pockets of matching size, buriedness and polarity are similar *without
   selection*. The null should match on joint strata of volume × burrow_score × lining polarity/aromatic fraction ×
   mean pLDDT. Report F_analog as excess over that stratified null. State in advance that a non-zero F_analog is
   expected under the "small pocket space" hypothesis and is not, by itself, evidence of convergence. Keep the
   shared-ligand PPV step. It is the right functional test (Barelier/Kahraman precedent).
4. **Keep H4 (remote-structural, no sequence signal) out of F_anc.** Report it as its own "ancestry-uncertain" class.
   Rationale: at the remote-structure level, homology and analogy are exactly what the field disputes, so folding H4
   into F_anc biases F_anc upward.
5. **Add a site-detectability preflight to §8.2 and to every control.** For each benchmark or control site, record
   whether the known residues lie in a BURROW pocket (≥50% residue overlap). Report convergence recall *conditional on
   detection*. Atlas evidence:
   - Chymotrypsin P00766: H57/D102/S195 are in no pocket; the best pocket contains only S189.
   - Trypsin P07477: 2 of 7 annotated sites in a pocket.
   - CA2 P00918: only T199/T200 of the Zn cone.
   - β-CAs P61517 and P53615: 0 of 4 Zn ligands.
   - γ-CA P40881: 1 of 8. Its site lies at the trimer interface.
   - SOD1/SOD2/SodA: 0 metal ligands.
   - Subtilisin P00782: H171 and S328 are in one pocket (score 0.34); D139 is not.

   Shallow catalytic triads and metal sites, which make up most of Riziotis's 34 cases, therefore mostly fail at
   *detection*, before similarity is ever tested. Add ligand-pocket convergence positives: cofactor sites (adenine/NAD/FAD
   recognition in different folds; Moodie 1996, Denessiouk 2001, Dym & Eisenberg 2001, all already in LITERATURE.md) and
   ProSPECCTs Barelier pairs.
6. **Treat oligomeric-interface sites as U, not A.** Monomer AF models lack sites that form at an interface:
   - γ-CA trimer;
   - uricase dimer (P25688: 4 of 9 site residues);
   - parts of LDH/MDH tetramer contacts;
   - ClpP rings.

   If the footprint overlaps a PDB biological-assembly interface (PDBe/SIFTS) or a UniProt "homodimer" note, set A→U
   and flag it. Rationale: otherwise a site absent in a monomer model gets called a "loss", which is an artifact class
   missing from §5.3.
7. **Specificity change should be residue-first.** All the classic cases were caused by 1–3 lining or second-shell
   residues, or by distant permissive substitutions:
   - Ortlund 2007;
   - Boucher 2014 eLife doi:10.7554/eLife.02304 (R102→W107f, moved by an *insertion*);
   - Zhu 2005 Science doi:10.1126/science.1106974 (IDH NAD→NADP).

   A whole-pocket fingerprint cosine will often not separate these. Make the clade-consistent lining-substitution test
   the primary criterion and the fingerprint split supporting evidence. Run a kinase gatekeeper-type single substitution
   as a deliberately small effect, to measure the sensitivity floor and report it.
8. **Let insertions contribute footprint columns.** Both GID1 (N-terminal lid) and apicomplexan LDH (loop insertion)
   build their new pocket from residues that are gap columns in most orthologs. The rule "footprint = columns present in
   ≥50% of member sites" will drop them. Use structure-based correspondence (US-align/Foldseek to a reference member) for
   case candidates, and allow lineage-specific columns in gain calls. Rationale: in documented "gains", new specificity
   comes from new structure on an old scaffold far more often than a pocket appears de novo (Bridgham 2010 PLoS Biol
   doi:10.1371/journal.pbio.1000497, "tinkering").
9. **Tone down M1/M4 as a discovery.** Homology transfer of sites is standard practice (AlphaFill, FINDSITE-style
   transfer, SFLD), and Panjkovich & Daura reported conserved pockets per family. The new contributions are: *rates* of
   turnover per unit of divergence, with an explicit observation-error model; the ortholog-vs-paralog contrast; and the
   exceptions. Frame M1 that way. M1 near 1 is expected, and the result is how far it falls below 1 and where.
10. **Gene loss vs pocket loss, and decaying genomes, as built-in negative controls.** Add two controls to §5.3:
    - Uricase (mouse P25688, rat) and GULO (mouse P58710, rat P10867) are lost as genes in hominoids/primates (Oda 2002
      MBE doi:10.1093/oxfordjournals.molbev.a004123). The pipeline must report **G** for human, never A.
    - M. leprae vs M. tuberculosis, both in tier 1: M. leprae has 1,116 pseudogenes against 6 (Cole 2001 Nature
      doi:10.1038/35059006). It is the natural stress test for fragment and truncation artifacts. Report the A-call rate
      for M. leprae orthologs before and after the §5.3 filters.
11. **Do not use olfactory or bitter-taste receptors as controls.** As monomers their AF models are adequate, but the
    families turn over fast by birth–death. That gives many-to-many orthology, frequent pseudogenes and few holo
    structures, so persistence and orthology filters cannot be checked cleanly. Use them only as an orthology stress
    group, if at all.

## Positive-control seed list (all present in snapshot; tier 1 unless noted; "atlas" = my spot check)
| # | Proteins (UniProt, species) | Expected event | Reference | Atlas |
|---|---|---|---|---|
| 1 | ESR1 P03372 (human), P57717 (zebrafish); AR P10275, GR P04150, MR P08235, PR P06401 (human) | specificity change, conserved LBD pocket (estrogen → 3-ketosteroid receptors) | Thornton 2003 doi:10.1126/science.1086185; Bridgham 2006; Ortlund 2007 | LBD pocket score 0.60–0.62; ESR1 E353/R394/H524 in pocket |
| 2 | Nurr1 P43354 (human), Q06219 (mouse); DHR38 P49869 (fly); NHR-6 P41829 (worm) | pocket loss after duplication (NR4A LBD has no cavity) from a ligand-binding NR ancestor | Wang 2003 Nature doi:10.1038/nature01645 (1OVL); Bridgham 2010 | human/mouse/fly: only weak LBD pockets ≤0.058. **NHR-6 has a 0.475 LBD pocket. Inspect it (regain or artifact?)** |
| 3 | PfLDH Q76NM3, PfMDH C6KT25 (PF3D7_0618500; UniProt mislabels it as LDH); human LDHA P00338, MDH1 P40925 | convergent/parallel specificity change (LDH from MDH in Apicomplexa) | Boucher 2014 doi:10.7554/eLife.02304 | substrate pockets 0.53/0.59 |
| 4 | E. coli IDH P08200 (NADP) vs human IDH3A P50213 (NAD) | coenzyme specificity change | Zhu 2005 doi:10.1126/science.1106974 | 9/12 site residues in one pocket |
| 5 | GID1 Q6L545 (rice), GID1A Q9MAA7 (Arabidopsis) vs HSL Q05469 (human) | function change / receptor-pocket gain on the HSL scaffold | Ueguchi-Tanaka 2005 doi:10.1038/nature04028 | 6/6 site residues, score 0.70 |
| 6 | D14 Q9SQR3, KAI2 Q9SZU7 (Arabidopsis), D14 Q10QA5 (rice) | specificity change after duplication (KAI2 ancestral) | Waters 2012 doi:10.1242/dev.074567 | 2/3 triad residues in pocket, 0.60–0.65 |
| 7 | Lysozyme C P61626 vs α-lactalbumin P00709 (human), P29752 (mouse) | catalytic-cleft function loss + Ca-site gain after a mammalian duplication | Qasba & Kumar 1997 doi:10.3109/10409239709082574 | LYZ cleft 2/2 (0.42); LALBA Ca site not in a pocket |
| 8 | Chymotrypsin P00766 (cow, tier 2), trypsin P07477, subtilisin BPN′ P00782 (tier 2), subtilisin E P04189 (tier 2), ClpP P0A6G7 | triad convergence across 3 folds | Wallace 1996; Riziotis 2025 | **detection-limited**: 0/3, 2/7, 2/3, –, 2/3 |
| 9 | α-CA P00918; β-CA P61517, P0ABE9 (E. coli), P53615 (yeast), P27140 (Arabidopsis); γ-CA P40881 (tier 2) | convergence (CA classes) | Hewett-Emmett & Tashian 1996 doi:10.1006/mpev.1996.0006 | **Zn sites largely not in pockets; γ at an interface. Detection negative.** |
| 10 | SOD1 P00441, SodC P0AGD1 vs SOD2 P04179, SodA P00448 | analogous enzymes (convergence) | Galperin 1998 doi:10.1101/gr.8.8.779 | **0 metal ligands in pockets. Detection negative.** |
| 11 | Uricase P25688, GULO P58710 (mouse), P10867 (rat); human has no protein | gene loss → must yield G, not A | Oda 2002 | uricase site partly at the dimer interface (see change 6) |
| 12 | M. leprae vs M. tuberculosis orthologs (tier 1) | artifact negative control for truncation-driven "loss" | Cole 2001 | – |

Recommended use:
- Cases 1–7 are the real Part B positives. They span gain, loss, specificity change and convergence, and the atlas sees them.
- Cases 8–10 are convergence controls. Report them honestly as detection-limited and do not drop them silently.
- Cases 11–12 are negative controls.
