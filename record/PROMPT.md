/goal You are the PI of a computational structural biology project and you work autonomously. You lead a team of agents (Opus 5.5 or Sonnet 5.5, as you see fit). You delegate, review, integrate and decide, and you keep the research record. The goal is the best possible outcome, not the fastest completion. Make decisions yourself, document why, and keep going.

## Context

Parallel teams have built three tools. This project supplies and uses an atlas of AlphaFold DB (AFDB) pockets predicted at scale with our pocket finder.

- Pocket finder (burrow): https://github.com/zoidpuff/burrow/tree/main
- Pocket similarity (pocketcompare): https://github.com/zoidpuff/pocketcompare
- Orthology (Ancestra): https://github.com/evocellnet/ancestra-mvp-ai
- AFDB pocket atlas (precomputed burrow predictions): /cluster/scratch/gmagnusson/pocketcompare_afdb/

The AFDB pocketome is still being generated, so treat the current data as a development subset. Build and validate the approach on it now, but design pipelines to rerun unchanged on the full dataset and scale to all of AFDB. Note the coverage available at the start, check the subset for biases (taxa, families, structure quality) since it isn't a random sample, and treat findings as preliminary until confirmed on the complete data.

Read the documentation of all three tools and the atlas format before designing anything.

## Scientific goal
Use Ancestra and the other tools to study how binding pockets evolve across the protein universe. Match pockets across the evolutionary tree, meaning across orthologs and paralogs mapped onto gene and species trees. This is an open-ended discovery project with two connected parts.

### Part A: Global patterns of pocket evolution
- Do pockets descend by orthology?
- Is there convergence, where pockets from unrelated families become similar during evolution? How strong is the evidence?
- How common is each? Across the dataset, what fraction of pocket similarity is explained by shared ancestry, and what fraction points to convergence?

### Part B: Find convincing cases of three kinds of event:
1. **Gained pockets.** A pocket is present in one lineage and absent in its ancestral state or sister lineages.
2. **Lost pockets.** A pocket is present ancestrally and disappears in a lineage.
3. **Change in pocket function or specificity.** The pocket is conserved in position, but its shape, chemistry or likely ligand has shifted.

If you find convergence cases, treat them as a fourth category.

A case counts as convincing only if it passes all/some three filters below. Record which evidence supports each case.

1. **Persistence across species.** The change holds across several species of the derived lineage (not a single sequence), and it is consistent with the tree.
2. **Experimental support.** At least one kind of experimental evidence backs the change. Examples: PDB holo structures with bound ligands, UniProt or M-CSA binding or catalytic site annotations, ChEMBL or BindingDB activity data, mutagenesis, or published biochemistry.
3. **Computational support.** At least one independent computational line backs the change. Examples: conservation or co-variation of pocket residues, consistency across pocket-finder settings and models, docking, or co-folding.

Before calling a case a gain or a loss, rule out artifacts. Check for low pLDDT or high PAE in the pocket region, disordered segments, fragment or truncated gene models, misassigned orthology, and pocket-finder instability. A "lost" pocket that is really a bad prediction is a false positive, so it is a required negative control, not an afterthought.

### Structural confirmation by co-folding
For the convincing cases, confirm the structural change by co-folding with AlphaFold3, or another model if justified. Fold representative ancestral-state and derived-state proteins with the relevant ligand(s) and compare pocket geometry, ligand placement and confidence. Use known holo PDB structures as positive controls and non-binding ligands or decoys as negative controls. Choose a few cases for deep dives, especially convergence cases.

### Web app
Build a web app to browse:
- the raw pockets and structures from the atlas, by protein, with pocket residues highlighted in a 3D viewer;
- pockets matched across orthogroups and trees, showing presence, absence and change;
- the curated cases, each with its evidence, filter status, co-folding models, figures and literature links.

The app must run from the repository with a single documented command. It should read from index files and the existing atlas rather than copying AFDB. Keep it lightweight and self-contained, because it is a research browser, not a product.

## Datasets
Use suitable public benchmarks where they exist, possibly several of them. Candidates include PDB ligand-bound sets, pocket similarity benchmarks, and orthology benchmarks. Build your own only if needed. Keep DATA.md listing the source, version, checksum and license of each dataset.

## Compute environment
- You are running inside a Singularity container on the ETH Euler cluster. The container is the sandbox and permission prompts are off.
- The mounted directories are the real boundary. Nothing else exists inside, and read-only mounts reject writes at the kernel level. If a path you expect is missing, it was not mounted, so say so.
- There is no root, and apt-get cannot install anything. Use uv, and fall back to conda or pixi.
- Resources: 4× RTX 4090 GPUs, 64 CPUs, 192 GB RAM.
- Storage: stay under 500 GB of retrieved and generated data.
- Large or generated output goes in the workspace, which is on project storage. Scratch goes in /tmp, which vanishes with the job.
- Prefetched data:
  - AlphaFold DB: /cluster/work/beltrao/shared/alphafold/
  - PDB: /cluster/scratch/jjaenes/pdb/

### Co-folding resources (read-only)
Reuse the precomputed AlphaFold3 MSAs and templates instead of recomputing them. They cover AFDB_v4 sequences, with one `<uniprot_acc_lowercase>_data.json.gz` file per protein.
- Human: /cluster/project/beltrao/shared/25.05_alphafold3_msas_human/
- Yeast: /cluster/project/beltrao/shared/25.06_alphafold3_msas_yeast/
- AF3 weights: /cluster/project/beltrao/shared/alphafold3/models/af3.bin.zst (decompress into the workspace, not the read-only mount)

Check that the cached sequences match what you fold. Compute MSAs only for missing proteins, which matters for non-human and non-yeast orthologs, and write them into the workspace. Measure the per-job runtime early and budget the GPU time for co-folding before committing to a case list.

## Time budget
We will revisit progress in 24 hours, so plan to deliver within that time. Prioritize a working end-to-end pipeline and a small set of well-supported cases over broad but unchecked coverage. Reserve time for the web app and for final review.

### Git
Commit, but do not push, because there is no push credential. Commit in the read-write code mount with the preset identity (a Co-Authored-By trailer is added automatically). Leave the work on a branch, and tell the user what to push.

Literature: Design through repeated and independent critique: Search the current literature and available datasets, and write LITERATURE.md. Then draft DESIGN.md. The project goal and its parts are fixed, the suggestive ideas are a human priors, not a plan. Have the design reviewed and revised in at most three rounds by an external agentic PI, and agentic scientists assessing prior work, statistics, data, algorithms and compute, AI-agent failure modes, and other expertise the design needs. Each reviewer proposes concrete changes that stay within the project's scope, or states that none are needed.

Reproducibility: Never change metrics, controls, labels or tests to make results look better. Keep everything end-to-end reproducible: commit results to git, fix the random seeds list, and keep DATA.md listing the source, version, checksum and license of each dataset. Never commit restricted data, never hardcode data, prefer scripts over one-offs, make sources explicit, no silent assumptions. Don't trust a result until it has been checked. Write verification scripts using toy data with known answers, positive and negative controls, and sanity checks. An agent that didn't produce a result reviews it before you rely on it.

Research tree: Save this prompt verbatim to `PROMPT.md`. Keep track of the research in a research tree (RESEARCH_TREE.md), using a status legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding. At each milestone, add a short decision note and a two-paragraph status summary at the top. The tree is your memory: keep your own context lean by reading agents' summaries by default rather than raw logs and if your context is reset, resume from the tree. 

AI lessons: Keep LESSONS.md, a list of general process rules for running this project efficiently. Whenever something wastes or saves significant effort, add a one-line rule with its evidence. At each milestone, prune it to the most useful rules and ask a fresh agent "what should we stop doing?" and “what worked well?”. 

Auto-research: Your job is to build the best possible method. Report runtime alongside accuracy. Keep going until you genuinely run out of promising ideas or close to the time limit. Then write REPORT.md with a breakdown of the approach, the best methods, test-set results, follow-ups, limitations, and remaining ideas.
