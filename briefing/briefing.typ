// Research briefing: pocket evolution across the AlphaFold DB pocket atlas, run by a team of AI agents.
// Build: python -m typst compile briefing/briefing.typ briefing/briefing.pdf   (pip install typst)
// Figures: .venv/bin/python briefing/make_figures.py

#let accent = rgb("#2a78d6")
#let warn = rgb("#eb6834")
#let soft = rgb("#f3f6fb")
#let ink2 = rgb("#52514e")

#set document(title: "How do binding pockets evolve? A test of agentic science on the AlphaFold pocket atlas",
              author: "Project lead agent")
#set page(paper: "a4", margin: (x: 2.1cm, top: 2.2cm, bottom: 2.2cm),
  header: context { if counter(page).get().first() > 1 [
    #set text(8pt, fill: ink2)
    Pocket evolution across the AlphaFold pocket atlas: a research briefing #h(1fr) agentic-pocket-evolution
  ] },
  numbering: "1")
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 0.95em)
#set heading(numbering: "1.")
#show heading.where(level: 1): it => { v(0.6em); set text(13pt, fill: accent); it; v(0.15em) }
#show heading.where(level: 2): it => { v(0.3em); set text(11pt); it }
#show figure.caption: it => { set text(9pt); set par(justify: true); align(left, it) }
#show figure: set block(spacing: 1.4em)
#set table(stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.25pt + luma(220)) }, inset: (x: 5pt, y: 3.5pt))
#show table: set text(9pt)
#show link: set text(fill: accent)

#let box-note(title, body, fill: soft, stroke: accent) = block(
  fill: fill, stroke: (left: 2.5pt + stroke), inset: (x: 10pt, y: 8pt), radius: 2pt, width: 100%,
  [#text(weight: "bold", fill: stroke)[#title] \ #set text(9.6pt); #body])

#let term(t) = text(weight: "bold")[#t]

// ------------------------------------------------------------------------------------------------------------
#align(center)[
  #v(0.4cm)
  #text(18pt, weight: "bold")[How do binding pockets evolve?]
  #v(-0.25cm)
  #text(13pt)[A test of agentic science on the AlphaFold DB pocket atlas]
  #v(0.15cm)
  #text(9.5pt, fill: ink2)[Research briefing · development subset, preliminary results · November 2026 \
  Written by the project's lead AI agent ("PI") from the committed research record; not yet peer reviewed]
  #v(0.2cm)
]

#box-note("In brief")[
A team of AI agents was given 24 hours, a cluster with four GPUs, three newly built tools and a pocket atlas
computed over the AlphaFold Database, and asked an open question: how do the binding pockets of proteins evolve
across the tree of life? The agents designed the study, had the design criticised by other agents, froze it,
built and ran an end-to-end pipeline, checked their own results with independent "verifier" agents, and wrote a
report. This briefing explains what they did and found, and what the exercise says about using AI agents for
bioinformatics.

- *Pocket similarity is strongly tied to ancestry pair by pair, yet much of it is not explained by detectable
  homology.* Pockets at the same position in orthologous proteins are about 14,000 times more likely to be
  significantly similar than pockets in unrelated proteins. Under the preregistered weighting (each pair of protein
  families counts once), 82 % (O1) to 94 % (O2) of the similarity has no detectable common ancestry; counting each
  pocket pair once (a framing chosen after the results) gives 37 %. These shares are upper bounds on convergent
  evolution, not estimates of it.
- *Convergent evolution of pockets is weak at this scale.* We found a small functional excess (about 8 %, with a 95 %
  interval up to 18 %; about 2 % for specific ligands), and the directional test found no convincing case of
  unrelated pockets becoming more alike, although it could miss weak, lineage-level convergence.
- *No individual gain, loss or change of a pocket passed all the evidence filters.* Gains are a real signal in
  aggregate, but one in five to one in four calls is expected to be false. Losses could not be told apart from noise, because even
  confident structure models miss 40 % of experimentally known binding sites. The pipeline also recovered none of
  seven textbook events from the literature.
- *As an agentic exercise, the team was fast and unusually self-critical, but unreliable as a project manager.* The
  design, the pipeline, and checks that caught real errors were in place within hours. The lead agent then sat idle
  for about 21 hours, missed the deadline by more than 11 hours, and left several decisions to be made after the
  results were known. All of these are documented.
]

= Why this question, and why agents?

Most proteins act by binding small molecules: substrates, cofactors, signalling molecules, drugs. Binding happens in
#term[pockets], concave patches on the protein surface whose shape and chemistry decide what fits. When proteins
evolve, their pockets can be kept, lost, gained, or reshaped to bind something new. Biologists have studied this
closely in a few families, for example how steroid hormone receptors acquired their specificities
@thornton2003 @bridgham2006, or how enzyme superfamilies keep their chemistry while changing their substrates
@babbitt1997. It has never been mapped systematically across whole proteomes, because doing so needs three things that
did not exist together: protein structures for everything, pockets found consistently on all of them, and reliable
family trees of the genes.

The AlphaFold Database (AFDB) now provides predicted structures for hundreds of millions of proteins @varadi2024. Three
tools built by parallel teams supplied the rest: #term[BURROW], a pocket finder; #term[pocketcompare], which turns
every pocket into a fingerprint vector so that pockets can be compared by a dot product; and #term[Ancestra], which
groups genes into families and reconstructs their history. A pocket atlas had been precomputed with BURROW and
pocketcompare over AFDB, and was still being generated when the project began.

The question was split into two parts. *Part A* asks about global patterns: do pockets descend along gene family trees,
and how much of the similarity between pockets is due to shared ancestry versus #term[convergence] (unrelated proteins
evolving similar pockets independently)? *Part B* asks for convincing individual events: pockets #term[gained] in a
lineage, #term[lost] in a lineage, or #term[changed in specificity] while staying in the same place. A case counts as
convincing only if (1) the change persists across several species of the lineage, (2) experimental evidence supports
it, and (3) an independent computational line supports it, after #term[artifacts] such as poorly predicted regions or
wrong gene models have been ruled out. Convincing cases were then to be checked by #term[co-folding]: predicting the
protein together with its ligand using AlphaFold 3 @abramson2024.

The second reason for the project was to test #term[agentic science]. One AI agent acted as principal investigator
(PI). It spawned other agents to do the work, to review the design and to verify results; it kept a research tree, a
log of decisions and a list of process lessons; and it worked largely without human input (the user intervened once,
to stop a duplicate lead-agent session). The instructions demanded fixed random seeds, positive and negative controls,
known-answer tests, review of every result by an agent that did not produce it, and an honest final report;
preregistration, null-input pipelines and the deviation ledger were proposed by the agents' own design reviewers. This briefing therefore reports two things at once: what was found about pockets, and how well an agent
team does science of this kind.

#box-note("Glossary", fill: rgb("#fafaf7"), stroke: ink2)[
#set par(justify: false)
#grid(columns: (auto, 1fr), column-gutter: 8pt, row-gutter: 4.5pt,
[*pLDDT*], [AlphaFold's per-residue confidence (0–100). Above 90 is "very high"; pockets in such regions match
experimental pockets about as well as experimental structures do @akdel2022.],
[*Ortholog / paralog*], [Genes in two species that descend from one ancestral gene through speciation (orthologs),
or through a gene duplication (paralogs).],
[*HOG, orthogroup*], [A hierarchical orthologous group: all genes descending from one ancestral gene at a given level
of the species tree. The project used two orthology runs, *O1* (46 reference proteomes) and *O2* (157 species).],
[*Site, site orthogroup (SOG)*], [A site is one pocket on one protein. A site orthogroup collects the sites that sit
at the same aligned position in the members of one gene family, i.e. "the same pocket" across species.],
[*Observation state*], [For each SOG and each protein: present, weakly present, absent, uninformative (model too
poor to tell), or no gene. Only confident absences count as evidence of loss.],
[*Null model*], [A version of the data in which the effect of interest has been destroyed on purpose, e.g. by
shuffling. Whatever the pipeline still "finds" there is the false-positive rate.],
[*FDR*], [False discovery rate: the expected fraction of reported hits that are false.],
[*Preregistration*], [Fixing the analysis plan, thresholds and success criteria before seeing the results. Later
changes are logged as deviations, with the direction in which they bias the result.],
[*Co-folding*], [Predicting a protein structure together with a ligand (here with AlphaFold 3) to see whether and
where the ligand is placed.],
)]

= How the work was organised

#figure(image("figures/fig1_pipeline.svg", width: 100%),
  caption: [*The pipeline.* Pockets from the atlas are mapped onto gene families to form site orthogroups. Each
  protein then gets an observation state per site. Part A uses all pocket pairs to ask how much similarity is
  inherited; Part B reconstructs gains, losses and specificity changes on gene trees and filters them hard.
  Co-folding and a web browser sit downstream. Known answers (experimental binding sites, literature events) are
  used throughout to measure error rates.]) <fig-pipeline>

The PI agent first read the documentation of the three tools and the atlas, surveyed the literature (127 DOIs of papers and
datasets, each checked against Crossref or Europe PMC) and wrote a design. Six reviewer
agents, each with a different brief (an external PI, statistics, data and compute, prior work, structural biology and
co-folding, and failure modes of AI agents), proposed concrete changes; a second round with two combined reviewers
followed. The design was then frozen together with a machine-readable file of 42 parameter values, and its checksum
was stamped into every output. Later changes were entered in a deviation ledger (31 entries) with the direction of
bias each introduces; one missing entry was caught by a verifier.

Work was split into stages, each run by one command on a frozen snapshot of the still-growing atlas, so that the
whole pipeline can be rerun unchanged on a later snapshot. Most headline numbers were reproduced by a fresh agent that had
not produced them before they were used (eight verification reports are part of the record); the O2 Part B results
were used for about five hours before they were verified. Results, figures and the
reasoning behind every decision were committed to git as they happened (223 commits in 36 hours).

== The data and its biases

The atlas came in tiers (@tab-data). Tiers 1 and 2 were complete. Tier 3, the whole AFDB, was about 10 %
done at the start, and the finished part was not a random sample: archaea and metagenomic sequences were over-represented about ninefold,
bacteria and eukaryotes under-represented, and within each taxon the finished part was a contiguous block. Because
persistence across species needed dense sampling of eukaryotes, the team computed its own pockets for 39 further
eukaryotic proteomes (tier L) in 87 minutes. Many eukaryotic models are of modest quality (mean pLDDT 73.7 for tier-1 eukaryotes),
so many observations had to be classed as uninformative. All findings are therefore *preliminary*. A second snapshot,
taken when tier 3 was complete (80,304 work units, 690 million pockets), was used to test reproducibility (§8); only
the analyses that use tier 3 saw new data in it.

#figure(
  table(columns: (auto, 1fr, auto, auto), align: (left, left, right, right),
    [*Tier*], [*Content*], [*Structures*], [*Pockets kept*],
    [1], [46 AFDB reference proteomes (model organisms)], [555,846], [1,788,750],
    [2], [Swiss-Prot (manually reviewed proteins)], [436,886], [1,442,678],
    [3], [whole AFDB (v4); 7,841 of 80,304 units done at snapshot 1], [22,962,643], [68,104,408],
    [L], [39 extra eukaryotic proteomes, computed by the team], [637,503], [1,837,831],
  ),
  caption: [*Data used.* Tiers 1 and 2 give 3,221,702 sites after merging fragments of long proteins. Sources,
  versions, checksums and licences are listed in `DATA.md`.]) <tab-data>

= Step 1: Matching the same pocket across species

To ask whether a pocket is kept or lost we must first know which pocket in a mouse protein corresponds to which pocket
in its human ortholog. The pipeline aligns the sequences of each gene family (FAMSA @deorowicz2016), maps every
pocket's lining residues onto alignment columns, and calls two pockets #term[positional homologs] if their columns
overlap enough (Jaccard index ≥ $tau$). Pockets are then grouped into site orthogroups by community detection (Leiden
@traag2019).

No benchmark existed for this task, so the team built one: pairs of ortholog sites that both have experimental ligand
evidence (from BioLiP2 @zhang2024biolip and UniProt). The overlap threshold was tuned on one half of the families and
tested on the other. At $tau = 0.2$ the matcher reached F1 = 0.68 [0.60, 0.71] on O1 and 0.70 [0.62, 0.78] on O2
(precision ≈ 0.8, recall ≈ 0.6). The tuned value sat at the edge of the search grid, so every headline was also
recomputed at $tau = 0.3$ (no conclusion changed). A leak of families between the two halves of the O2 split was found and
fixed before the O2 test numbers above were computed.

The result was 260,991 site orthogroups on O1 and 463,014 on O2. Matching is imperfect: about two in five known
correspondences are missed even when both pockets were detected, and counting pockets the finder did not detect, only
about a third (0.31–0.36) of known correspondences are recovered. This matters for everything that follows.

= Step 2: How often is a real pocket missed?

A claim that a pocket was "lost" rests on not seeing it. The design therefore required measuring how often the
pipeline fails to see a pocket that is known to exist, before any event was called. The team collected 73,107
experimentally known binding sites (BioLiP2, UniProt with experimental evidence codes, LIGYSIS @utges2024), kept only
those in high-confidence model regions (pLDDT ≥ 90, a threshold taken from @akdel2022) and asked how often the site
would be called "absent".

#figure(image("figures/fig2_absence.svg", width: 88%),
  caption: [*Absence is weak evidence.* Fraction of experimentally known binding sites, in confident AlphaFold
  regions, that the pipeline would call absent. Metal sites are missed most often. "Rescue" re-runs the pocket finder
  with relaxed settings and other random seeds and turns an absence into "weakly present" if any run finds the pocket;
  it costs specificity (it also rescues 35 % of regions where no pocket exists). Sites larger than 10 residues come
  closest (0.054), but no setting reaches the preregistered 5 % target. Tune half of the known sites; the test half gives 0.14 after rescue.])
  <fig-absence>

The answer was sobering: 40 % of known sites would be called absent (@fig-absence). The same structure with the same seed
gives the same pockets, but a different seed re-finds only 81 % of sites, and whether a site forms a pocket depends on
details of the predicted conformation. A likely further reason, known from the literature but not measured here, is
that many real binding sites are shallow or form only when the ligand is bound. Even comparing the same protein in two AFDB versions (v4 and v6)
re-detects only 65 % of sites (72 % in very-confident models).

The team responded with three conservative rules: absences count only where the model is very confident; absences
count only for footprints larger than 10 residues; and every absence must survive a rescue attempt. No setting met the
preregistered 5 % false-absence target overall, and this was logged as a deviation rather than hidden. The measured
rates were then fed into simulations used to set the thresholds for calling events. This step shaped every later
result: in this data a single missing pocket says little.

= Step 3: Part A — how much pocket similarity is inherited?

== The approach

Pocket shapes are not infinitely varied. Earlier work showed that a few hundred to a thousand representative shapes
cover most natural pockets, and that pockets resembling natural ones appear even in computer-generated proteins that
were never selected for anything @gao2013 @skolnick2013. Pocket similarity between unrelated proteins is therefore
expected even without any convergent evolution. The design took this as the null hypothesis.

For all 3.2 million tier-1/2 sites the team found each site's 1,000 nearest neighbours by fingerprint (on GPU) and
called a pair significant only if it was more similar than expected for pockets of the same size, buriedness,
chemistry and model quality from unrelated proteins (expected false hits per query < 0.01, with the tail of the null
distribution fitted). This gave 9.8 million significant pairs. Each pair was then classified by the relationship of
the two proteins: orthologs or paralogs at the same site, the same proteins at different sites, members of the same
domain family, only remotely similar in structure (Foldseek @vankempen2024), or #term[no detectable homology]. From this
the team estimated $F_"analog"$, the share of significant similarity that is not explained by any detectable homology.

== Results

#figure(image("figures/fig3_partA.svg", width: 100%),
  caption: [*Similarity is strongly tied to ancestry per pair, but similarity without detectable homology is common
  in aggregate.* (a) Probability that a pocket pair is significantly similar, by the relationship of the two proteins
  (O1; 95 % bootstrap intervals over families). Same-site orthologs: 0.010; unrelated proteins: $7 times 10^(-7)$.
  (b) $F_"analog"$ under two weightings. Counting each pocket pair once (post hoc framing) gives 0.37 on both
  orthologies; counting each pair of families once (the preregistered headline) gives 0.82 on O1 and 0.94 on O2, because that weighting depends on
  how finely the families are cut. The per-pair value also differs between the held-out halves (0.29 on O1, 0.14 on
  O2). Restricting to very confident models lowers the per-pair value to 0.16.])
  <fig-partA>

Three findings stand out (@fig-partA).

*Ancestry dominates at the level of individual pairs.* Two pockets at the same position in orthologous proteins are
about 14,000 times more likely to be significantly similar than two pockets from unrelated proteins.

*But unrelated pairs are so numerous that they still contribute a lot.* Of all significant pairs, 37 % involve proteins
with no detectable homology. This drops to 16 % when both models are very confident, which suggests that part of the
"analogous" similarity comes from poorly modelled regions. The preregistered headline, which weights each pair of
protein families equally, is 0.82 on O1 but 0.94 on O2. The difference is not biological: O2 groups are coarser on
average, which changes how many family pairs there are. The per-pair estimate agrees between the two orthologies on all
pairs (0.37 on both) but not on the held-out halves (0.29 vs 0.14), and it was given prominence only after the results
were seen; neither weighting is a stable single number.

*$F_"analog"$ is an upper bound on convergence, not a measurement of it.* Two further tests asked whether the analogous
similarity looks like convergent evolution. A functional test compared how often significant unrelated pairs bind the
same class of ligand with matched random pairs: there is a small excess (2.6-fold, giving an estimated convergent
fraction of 0.08 [0.01, 0.18]), but it comes almost entirely from nucleotides and metals; with specific, non-metal
ligands it is 0.02 [0, 0.06]. A #term[directional] test asked whether a pocket in one family became _more_ similar to a
pocket in an unrelated family along a lineage than its sister lineages are. Of 2.6 million testable pairs, 3 passed,
all with artifact flags; this is below the false-pass rate measured on a simulated world without convergence. The
test has power for strong effects (0.96 in simulation) but could miss weak convergence spread over a lineage. Our
reading: most similarity between unrelated pockets is consistent with the limited number of possible pocket shapes,
with a small excess that could be convergence. One caveat limits power: in controls with known convergent active sites
(the serine-protease triad, carbonic anhydrases), the catalytic residues mostly lie outside BURROW pockets, so this
pipeline cannot see some of the classic convergence cases.

== Do pockets descend along the gene tree?

#figure(image("figures/fig4_descent.svg", width: 72%),
  caption: [*Pocket retention falls with sequence divergence, and orthologs keep pockets no more often than paralogs.*
  Probability that a site present in one protein is also found at the aligned position of the other, given that the
  other model is good enough to tell, by sequence identity (O1, preregistered quality rule, all pairs, 95 %
  intervals).]) <fig-descent>

Yes, broadly. Where the other protein's model is good enough to judge, a site is found again in 77 % of ortholog
pairs under the preregistered rule (all pairs, O1), and in 89 % when the same quality rules are applied to both proteins
(held-out half, O1); for comparison, re-detecting a large site in two models of the same protein succeeds 87 % of the
time. Retention falls steeply only below
about 20 % sequence identity (@fig-descent). Presence patterns across species also fit the gene trees better than a
shuffled control that keeps each species' and each pocket's totals (0.933 vs 0.919, p = 0.001), a modest but real
effect.

The #term[ortholog conjecture] holds that orthologs keep their function better than paralogs; for functions in general
it was contested @nehrt2011 and later found to be weak but significant @altenhoff2012. For pocket presence we found no such advantage: at matched
sequence identity, orthologs keep a pocket slightly _less_ often than paralogs (−1.1 percentage points [−1.8, −0.5] on
O1, −0.2 [−0.8, +0.5] on O2; held-out half, same quality rules on both sides). We report this as "not supported", not as "orthologs are less conserved": the age of a
gene duplication is unknown, so the comparison is confounded by time.

= Step 4: Part B — individual gains, losses and changes

== The approach

For each site orthogroup, the presence or absence of the pocket was mapped onto the gene tree of its family, and
Sankoff parsimony @sankoff1975 found the cheapest history of gains and losses (a gain costs twice a loss, so gains must
be well supported). The decision thresholds were set by simulation, using the measured error rates of Step 2, to give
a 10 % false discovery rate. Two independent safeguards were then added.

*Null-input pipelines.* The whole pipeline was rerun 20 times each on two scrambled versions of the data: with
presence and absence shuffled among species while keeping row and column totals ("curveball" randomisation
@strona2014), and with random groups of species posing as lineages ("pseudo-clades"). Whatever these runs call is false
by construction.

*An artifact cascade.* Every candidate was checked for the known ways a pocket can appear or vanish without
evolution @tunyasuvunakool2021 @meyer2020 @barriohernandez2023: low model confidence at the pocket, disordered
segments, fragmentary or truncated gene models on either side, wrong orthology (checked against a second orthology and
PANTHER @thomas2022, and by reciprocal best hits with DIAMOND @buchfink2021), and whether the lineage makes sense on
the species tree. Cases that survived were then tested against the three evidence filters. Specificity changes used a
separate test: at least two lineage-consistent, chemically non-conservative substitutions in the pocket plus a
fingerprint shift beyond the noise floor. In parallel, an *annotation-first validation track* followed seven
well-documented events from the literature (for example the mineralocorticoid vs glucocorticoid receptor, lactate vs
malate dehydrogenase, NADP- vs NAD-dependent isocitrate dehydrogenase, lysozyme vs α-lactalbumin) through the same
outputs, to measure recall.

== Results

#figure(image("figures/fig5_partB.svg", width: 100%),
  caption: [*Gains carry a statistical signal; losses do not; no case passes everything.* (a) Event calls on real data
  and the mean over 20 runs of each null (O1 calls before rescue, O2 after). Gains exceed both nulls; losses do not
  exceed the shuffled-state null. (b) Candidates after rescue and the checks they pass. "Persistence + 1 other filter"
  and "pass artifact checks" are assessed separately (not nested); "convincing" requires all three filters and the
  artifact cascade.]) <fig-partB>

*Gains are real in aggregate but individually unreliable.* On O1 the pipeline called 175 gains, against 35 and 14 on
average in the two nulls (@fig-partB a). The simulated FDR was 5 %, but the empirical FDR from the shuffled null is
20 % (27 % on O2). The simulation was too optimistic, and no threshold was changed after this was seen.

*Losses carry no signal.* The shuffled null calls as many losses as the real data (66 vs 61 on O1). Given the 40 %
false-absence rate, this is what Step 2 predicts. Losses were kept only at "audit" level.

*Most gain candidates are artifacts of orthology.* After rescue, 50 O1 gains remained; only 7 passed the artifact
cascade. The commonest failure was orthology: for 79.5 % of the outgroup species whose absence made an O1 gain, the
O2 orthology places the corresponding proteins in a different family, so the "gain" compared different genes (O1
had merged families that O2 separates, although O2 groups are coarser on average).
Even within one O2 family, the absent outgroup was usually a different paralog copy rather than a true ortholog (only
20 % were reciprocal best hits). Fragmentary gene models on the present side were the second most common cause. (The orthology check used on O2,
reciprocal best hits, was chosen after the planned check turned out to be uninformative on O2; the planned rule would
have passed 34 O2 gains instead of 17.) About 23 % of the absences behind O2 gains were never put through the rescue
step.

*Zero convincing cases.* No gain, loss, specificity change or convergence candidate passed all three filters
(@tab-cases). The experimental filter failed almost everywhere: no experimental record lies at the gained or lost
pocket in the relevant state; a verifier searched UniProt, PDBe and RCSB independently and found nothing that had been
missed. The specificity test produced no individual call (the null resolution was too coarse for the multiple-testing
correction to ever fire), although at an earlier stage of the test real lineage pairs outnumber null pairs about
6-fold (513 vs 71–88 on O1, 3 null seeds). Most of that enrichment comes from the residue filter; the fingerprint shift
adds only 1.4–2-fold. It is a population-level signal without nameable cases, found after the fact.

#figure(
  table(columns: (auto, auto, auto, auto, auto, auto), align: (left, right, right, right, right, right),
    [*Category*], [*Candidates*], [*Persistence*], [*Experimental*], [*Pass artifacts*], [*Convincing*],
    [O1 gains], [50], [42], [1], [7], [*0*],
    [O1 losses], [17], [8], [2], [1], [*0*],
    [O2 gains], [104], [76], [0], [17], [*0*],
    [O2 losses], [66], [19], [3], [6], [*0*],
    [Specificity (top 5, O1)], [5], [2], [0], [0], [*0*],
    [Convergence (directional passes)], [3], [3], [0], [0], [*0*],
  ),
  caption: [*Evidence filters per category* (frozen case list). The computational filter passed often (47 of 50 O1
  gains), the experimental filter almost never.]) <tab-cases>

*Recall of known events: 0 of 7.* This is the most informative negative result. Several of the textbook events happen
between paralogous families that the orthology places in different groups (MR/GR on O1, IDH, GID1/HSL), which a
pipeline that compares pockets within one group cannot see. The others are within one group but are missed because the
pocket is detected in both states, because absences are not confident enough, or because the specificity test is not
significant. An exploratory extension that matched sites across related groups also recovered 0 of 7 and was stopped
by its own pre-committed gate. The negative controls behaved: uricase (and GULO on O1), lost in humans, appear as missing
genes rather than pocket losses, and truncated _M. leprae_ proteins gave no spurious losses.

= Step 5: Co-folding as an independent check

Recent work shows that co-folding models often place a ligand in a familiar-looking pocket even after every residue
of the site has been replaced by glycine @masters2025, and that they largely reproduce poses seen in training
@skrinjar2025. A confident co-folded complex is therefore weak evidence that a pocket exists. The team calibrated AF3
before folding any case: 24 experimentally confirmed complexes (cofactors, drug-like molecules, metals, nucleotides)
against 72 negatives (48 property-matched decoy ligands and 24 site-ablated proteins). Thresholds and a decision rule
were frozen before any case was folded.

#figure(image("figures/fig6_cofold.svg", width: 100%),
  caption: [*Co-folding is informative only one way.* (a) How well single AF3 readouts separate true complexes from
  decoys and ablations (AUC with 95 % intervals). For metal ions, ligand confidence is at chance (geometry
  still works). The frozen combined rule has specificity
  0.96 but sensitivity 0.42: a "binds" call means something, a "does not bind" call much less. (b) Ligand-free AF3
  refolds of both sides of each candidate: fraction of the pocket's residues re-found as a pocket (median of three
  models). For gains, filled = the lineage with the gained pocket. For false-loss controls (losses whose absent side
  had poor AFDB models), the prediction was that a fresh model would restore the pocket on the open-circle side; in
  HOG0030328 both sides give 0.6.]) <fig-cofold>

The interface-confidence score (ipTM) was dropped because decoys and site-destroyed proteins kept high ipTM. Ligand
confidence and geometry worked moderately well (AUC 0.73–0.79); for metal ions ligand confidence is at chance, although
the geometry readouts still separate them (@fig-cofold a). The results on
cases were mostly negative or ambiguous:

- In the three gain cases whose positive control passed, AF3 placed the family's known ligand 12–34 Å away from the
  gained pocket, in another pocket. The gained pockets are *secondary sites* of unknown function; co-folding cannot
  confirm them at the ligand level.
- Ligand-free refolds reproduced 3 of 4 testable gains: the pocket appears in fresh AF3 models of the gaining lineage
  and not of its relatives, so the gains are not artifacts of one AFDB model.
- The prediction for the false-loss controls was mostly refuted: in 3 of 4 testable cases a fresh model also lacked the
  pocket. Either some of these are real losses discarded by the quality filter, or AlphaFold 2 and 3 share the same
  error; the data cannot tell which.
- *Deep dives with known answers.* For the steroid receptors (aldosterone vs cortisol), AF3 gave the expected
  receptor-by-ligand interaction in all 3 seeds, but the effect is small (about 3 pLDDT units), all four combinations
  passed the binding rule, and these complexes are likely in AF3's training data. For lactate vs malate dehydrogenase
  the positive control failed, so the result was reported as uninterpretable rather than rescued with a different
  control after the fact.

Case scoring was meant to be blinded. A post-processing bug regenerated the blinding key, and the mapping was recovered
by sequence before the final scoring pass, so case scoring is reported as unblinded (the folds themselves ran blind).
Runtimes were modest: one AF3 job with three seeds took 80–310 s on an RTX 4090, and the whole calibration panel took
1.7 GPU-hours.

= What the agents did well, and what went wrong

#figure(image("figures/fig7_process.svg", width: 100%),
  caption: [*Activity over the run* (commits per hour, CET). Most of the science was done in the first nine hours. The
  shaded periods are documented failures: a co-folding queue that stalled without anyone noticing, and about 21 hours in
  which the lead agent did nothing because it had not scheduled a wake-up. An unattended rerun on the second snapshot
  ran in that window and failed in three stages. Source: `record/commit_timeline.tsv`.]) <fig-process>

== Strengths

*Speed with structure.* Within about one hour the agents had read three new tools, audited the atlas for bias,
reviewed the literature, drafted a design and frozen it after two review rounds (several agents worked in
parallel). Within nine
hours the pipeline had produced Part A, Part B on both orthologies, measured error rates and a calibrated co-folding
protocol (@fig-process).

*Self-correction that changed conclusions.* The most valuable work was negative. Measuring the false-absence rate
before calling events (40 % against a placeholder of 5 %) turned "loss" from a headline into an audit category.
Rechecking O1 gains against the finer O2 orthology removed most of them. Verifier agents reproduced most headline
numbers before use and found mislabelled quantities, a quality rule applied to only one side that inflated pocket
retention, a simulated FDR that was too optimistic, and a missing ledger entry. Other checks during the run found a
leaky train/test split and a bug in the co-folding decision function that had briefly recorded a wrong verdict.

*Honest reporting.* The final report leads with "0 convincing cases" and "recall 0/7", reports the preregistered
headline next to the post hoc one, gives every post-freeze change in a ledger with its bias direction, and lists its
own process failures. Every number in the report is checked automatically against a machine-readable results file.

*Reproducibility.* A second, larger atlas snapshot was processed with the same commands. After at least six
re-runnability defects were fixed, Part A and the O1 analysis, whose inputs were identical, were reproduced exactly
(byte-identical site orthogroups and identical headline numbers): a test of reproducibility, not of robustness to new
data. The O2 results, which did see new data, were stable in their counts and conclusions, although individual O2
candidates changed.

== Weaknesses

*Project management.* The lead agent went quiet for about 21 hours without scheduling a wake-up; the deadline was
missed by more than 11 hours. Earlier, the co-folding queue stalled for about four hours without detection. After a client
restart two lead-agent sessions ran in parallel, and an agent wrote a "PI review" into the record under false
attribution; it was reverted.

*Decisions after the results.* Several choices were made after seeing data: the per-pair framing of Part A, the
reciprocal-best-hit orthology rule in the cascade, the population-level specificity analysis and the cross-family
extension. They are labelled post hoc and none changes the count of convincing cases, but some shape how results are
presented, and they show how easily an agent team drifts from
its plan when results disappoint.

*Fragile engineering.* The first unattended rerun failed on inputs that were tied to the first snapshot and overwrote
a section of the results file (restored from git). Several outputs carry code hashes from uncommitted working trees.
The web app's 3D viewer was never opened in a browser by an agent; the container had no headless browser.

*Limits of the evidence filters.* The experimental filter is very hard to satisfy for new events: evidence about a
pocket in the lineage that has it rarely exists. The bar was set by the task; it produced a rigorous null result
rather than a list of discoveries.

= Take-away messages

+ *For pocket biology.* Pockets are largely inherited and tend to be retained along gene trees until sequences diverge
  below about 20 % identity. Similarity between unrelated pockets is common and is consistent with the limited
  number of pocket shapes; convergence is a small effect at most at this scale, as far as these tests can see. These are development-subset results
  and need the full atlas.
+ *For proteome-scale inference of gains and losses.* Absence is the weak link. When 40 % of known sites are missed
  even in confident models, losses cannot be called from single structures, and gains depend on correct orthology more
  than on the pocket finder. Better orthology (paralog-aware comparisons), a measured false-presence rate, likelihood
  models of noisy traits @lewis2001 and cross-family event calling are the next steps; the 0/7 recall shows where the
  current design is blind.
+ *For co-folding as evidence.* Calibrate before use, include site-ablation negatives, and treat a "no binding" result
  as weak. AF3 reproduced gained pockets in fresh models, but could not tell what binds them.
+ *For agentic science.* An agent team can design, build and stress-test a nontrivial bioinformatics study in hours.
  What made it trustworthy were established scientific practices: the instructions asked for controls, known-answer
  tests and independent verification, and the agents' own design reviewers added preregistration, null-input pipelines
  and a deviation ledger. Its failures were both organisational (attention, scheduling, parallel sessions, post hoc
  drift) and engineering (bugs in reruns, blinding and bookkeeping), and every one of them was visible because the
  record was kept in git. Human oversight remains essential at the points where the
  record shows the agents failed: keeping the work moving, and resisting decisions made after the results were seen.

= What worked, and whose ideas it rests on

The parts of the method that held up under verification are listed below with the published work they draw on. These
ideas were extracted by the agents during the literature review (`record/LITERATURE.md`), before the design was written.
"Held up" means the component did its job, which here often meant exposing a limitation (the validation track, for
example, worked as a recall test and returned 0 of 7).

#figure(
  table(columns: (1.2fr, 1.5fr, 1.6fr), align: left,
    [*Component*], [*Idea taken from*], [*How it was used*],
    [Absence only in very confident models], [Pocket accuracy depends on pLDDT > 90 @akdel2022], [Hard gate for "absent";
      pLDDT ≥ 90 stratum reported for Part A],
    [Known-site error rates], [Pocket-finder recall varies widely between tools @utges2024], [False-absence rate
      measured on 73,107 sites before any event was called (the rescue step came from the design reviewers)],
    [Null for pocket similarity], [Pocket-shape space is small and reachable without selection @gao2013 @skolnick2013],
      [Covariate-matched null; $F_"analog"$ framed as an upper bound],
    [What counts as convergence], [Define convergence explicitly @doolittle1994; check site, mechanism and ancestry
      @riziotis2025], [Separate functional and directional tests; homology classes],
    [Ortholog vs paralog contrast], [Ortholog conjecture and its small effect size @nehrt2011 @altenhoff2012],
      [Matched-identity comparison, preregistered],
    [Artifact cascade], [Disorder and length limits @tunyasuvunakool2021; gene-model errors @meyer2020; species-specific
      clusters @barriohernandez2023; apo/holo bias @saldano2022], [Required checks before any gain or loss],
    [Co-folding controls], [Site-ablation challenge @masters2025; memorisation @skrinjar2025; confidence ≠ binding
      @kim2026], [Calibration panel with decoys and ablations; relative evidence only (dropping ipTM followed from
      the team's own control runs)],
    [Known-answer events], [Resurrected steroid receptors @thornton2003 @bridgham2006; ancestral vs derived design
      @harms2013; conserved chemistry, shifting specificity @babbitt1997], [Validation track and deep dives],
  ),
  caption: [*Credits.* Experimental annotations came from BioLiP2 @zhang2024biolip, UniProt @uniprot2025, M-CSA
  @ribeiro2018, LIGYSIS @utges2024, PANTHER @thomas2022 and TED @lau2024. Standard algorithms (Sankoff parsimony
  @sankoff1975, Leiden clustering @traag2019, curveball randomisation @strona2014, FAMSA @deorowicz2016, DIAMOND
  @buchfink2021, Foldseek @vankempen2024, Benjamini–Hochberg @benjamini1995, ColabFold MSAs @mirdita2022, AlphaFold 3
  @abramson2024) are cited here for the reader; they were used as tools rather than taken from the literature
  review.]) <tab-credits>

One idea from the literature did not hold up as used. The literature review stressed that a structurally inferred
trait is observed with error, and pointed to likelihood models for binary traits on trees @lewis2001 @cohen2010. The
team deferred these and calibrated parsimony thresholds by simulation with the measured error rates instead; the
simulated FDR for gains (0.05) proved too optimistic against the null-input pipelines (0.20–0.27). A likelihood
observation model remains the main methodological follow-up.

The three project tools were built by parallel teams with AI agents and are described in their repositories:
BURROW (#link("https://github.com/zoidpuff/burrow")), pocketcompare (#link("https://github.com/zoidpuff/pocketcompare"))
and Ancestra (#link("https://github.com/evocellnet/ancestra-mvp-ai")). Their own documentation states that their
claims are not yet human-verified; every result here inherits that caveat.

= Limitations

All results come from a development subset of the atlas with known taxonomic biases. Tier 3 was about 10 % complete at
the main snapshot; the second snapshot changed no conclusion, but only the O2-based analyses saw new data and
orthology was not recomputed for the new species. Pocket models come from two AFDB versions (v6 files and compressed
v4 models), and model source is a covariate.
Orthology at deep (eukaryote-wide) splits is uncertain, and 16 % of human–mouse pairs that PANTHER places in one
subfamily fall into different O2 groups. The ortholog–paralog contrast is confounded by the unknown age of
duplications. Co-folding positive controls are probably in AF3's training data. The preregistered Part A headline
depends on how finely families are defined. Many choices after the design freeze are logged as deviations; each is
labelled with the direction of its bias, but together they leave more room for the analyst's choices than a strict
preregistration would. This briefing was written by the same agent system that ran the project and has not been
reviewed by a human expert.

= Reproducing this work and browsing the results

The repository contains the code, configuration, small result files, case records and the full research record
(prompt, design, research tree, lessons, reviews, verification reports and report). Bulk intermediates (site tables,
SOG tables, orthology XML, structures, AF3 models) are excluded and are regenerated by
`scripts/run_pipeline.py` on a cluster with access to the atlas and AFDB. The research browser (a small web app) shows
proteins with their pockets in 3D, the SOG × species matrix of each gene family, and every case with its evidence; the
repository README explains how to open it on the cluster and how to run its demo mode anywhere.

#v(0.6em)
#set text(9pt)
#bibliography("references.bib", title: "References", style: "nature")
