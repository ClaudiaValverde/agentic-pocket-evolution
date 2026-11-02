# How do binding pockets evolve? An exercise in agentic AI for bioinformatics

This repository is an **experiment in agentic science**. A team of AI agents (one lead agent acting as principal
investigator, plus agents it spawned to build, review and verify) was given 24 hours, a GPU node, three newly built
tools and a pocket atlas computed over the AlphaFold Database, and asked an open research question:
*how do protein binding pockets evolve across the tree of life?* The agents designed the study, had it reviewed,
froze it, built and ran the pipeline, verified their own results and wrote the report, without a human in the loop.

The repository holds what is needed to understand, check and rerun that work: the code, the configuration, the small
result files behind every reported number, the curated case records, and the full research record (the prompt, the
design and its reviews, the decision log, the verification reports and the final report). It is meant to be read as
much for **how the agents worked** as for what they found.

**Start here:** [`briefing/briefing.pdf`](briefing/briefing.pdf), a research briefing for a broad scientific
audience: what was done, what was found, what went wrong, and whose ideas the methods rest on.

> **Status.** All results come from a development subset of the atlas and are preliminary. The briefing and report
> were written by the lead agent and have not yet been reviewed by a human expert. The three tools the project builds
> on were themselves built with AI agents and describe their own claims as not yet human-verified.

## What was found, in short

| Question | Answer (development subset) |
|---|---|
| Do pockets descend along gene trees? | Broadly yes: a site is re-found in 77 % of informative ortholog pairs (preregistered rule; 89 % with symmetric quality rules on a held-out half), retention falls with sequence divergence, and presence patterns fit gene trees modestly better than shuffled data. Orthologs keep pockets no better than paralogs at matched identity. |
| How much pocket similarity is inherited? | Same-site orthologs are ~14,000× more likely to be significantly similar than unrelated pockets. Yet under the preregistered weighting (each family pair once) 82 % (O1) to 94 % (O2) of similar pairs have no detectable homology; counting each pocket pair once (post hoc) gives 37 %. Both are upper bounds on convergence, not estimates of it. |
| Is there convergent pocket evolution? | Weak: a small functional excess (0.08, 95 % CI 0.01–0.18; ≈ 0.02 for specific ligands); the directional test found no convincing case, though it could miss weak lineage-level convergence. |
| Convincing gains, losses, specificity changes? | **None** passed all evidence filters. Gains are a population signal with an empirical FDR of 0.20–0.27; losses cannot be separated from noise, because 40 % of known binding sites are missed even in confident models. Recall of 7 literature events: 0/7. |
| Does co-folding (AlphaFold 3) help? | After calibration it is informative one way (specificity 0.96, sensitivity 0.42). Gained pockets reappear in fresh ligand-free models but are secondary sites, not where known ligands bind. |
| How did the agents do? | Fast, structured and self-critical (preregistration, null models, independent verifiers, a deviation ledger), but the lead agent idled ~21 h, missed the deadline by more than 11 h, and made several post hoc choices, all documented. |

Details and numbers: [`briefing/briefing.pdf`](briefing/briefing.pdf) and the agents' own final report,
[`record/REPORT.md`](record/REPORT.md).

## Repository layout

```
briefing/        research briefing (PDF + Typst source), figure script, bibliography builder
record/          the agents' research record
  PROMPT.md        the task exactly as given to the lead agent
  LITERATURE.md    literature and dataset review (127 DOIs, each checked against Crossref/Europe PMC)
  DESIGN.md        preregistered design, review decisions and the deviation ledger (§12)
  RESEARCH_TREE.md research tree, status summaries and decision log
  LESSONS.md       process rules the agents wrote for themselves
  REPORT.md        final report (every number machine-checked against results/numbers.json)
  CROSSOG.md       exploratory cross-orthogroup extension (negative)
  reviews/         design reviews (2 rounds), retrospectives and independent verification reports
  commit_timeline.tsv  time and subject of the 223 commits of the original run (no author data)
src/pocketevo/   pipeline library (sites, site orthogroups, states, Part A, Part B engine)
scripts/         one script per stage; scripts/run_pipeline.py runs everything on an atlas snapshot
config/          preregistered parameters (prereg.yaml; its sha256 is stamped into every output) and orthology config
results/         small result files: results/numbers.json (all headline numbers) and the JSON/CSV/MD summaries
cases/           curated case records (case.yaml with evidence, filter status and co-folding readouts)
webapp/          research browser (FastAPI + 3Dmol.js), with a tiny synthetic fixture for demo mode
tests/           unit tests, toy datasets with known answers, and data-dependent tests that skip without data
docs/            tool briefs, co-folding protocol, evidence rules, environment notes
DATA.md          every dataset with source, version, checksum and licence; pinned tool commits
```

Bulk data are **not** in the repository: site and SOG tables, orthology XML, annotation tables, AFDB structures,
pocket fingerprints, MSAs and AlphaFold 3 models (see `.gitignore`). They are regenerated by the pipeline or stay at
their source (see `DATA.md`). AlphaFold 3 weights and outputs fall under DeepMind's terms of use and are never
included.

## Setup

Requires Linux, Python 3.12 and [uv](https://docs.astral.sh/uv/).

```bash
git clone <this repository> agentic-pocket-evolution && cd agentic-pocket-evolution
bash scripts/setup_external.sh            # clones BURROW, pocketcompare and Ancestra at the pinned commits
uv sync                                   # core environment in .venv
uv pip install --python .venv/bin/python -r requirements.lock -r webapp/requirements.txt
.venv/bin/python -m pytest -q tests/      # 170 pass; tests that need bulk data or the cluster are skipped
```

## The web app (research browser)

The browser shows (1) atlas proteins with their pockets highlighted in a 3D viewer, (2) each gene family as a
site-orthogroup × species matrix of present / weak / absent / uninformative states on the species tree, and (3) every
curated case with its evidence, filter status, co-folding readouts and side-by-side 3D models. It runs locally and is
opened in your own browser. How much it can show depends on where you run it.

### Option A: demo mode, anywhere (no cluster needed)

After [Setup](#setup):

```bash
uv run --project . python -m webapp --port 8765
```

Then open <http://127.0.0.1:8765/>. Available: all case pages (`/case`) with their evidence and co-folding
tables, the validation table of known events (`/validation`), and a small synthetic orthogroup (labelled TOY) to try
the SOG × species matrix. Protein pages and real orthogroups need the bulk data, so they show as unavailable.

### Option B: the full browser on the ETH Euler cluster

This needs access to the AFDB and atlas mounts listed in `DATA.md` and to a full copy of the project workspace (a
checkout in which the pipeline was run, so that it holds the bulk results). Start the app on the cluster and point it
at that workspace:

```bash
bash scripts/webapp_with_workspace.sh /path/to/full/workspace 8765
```

The server listens on `127.0.0.1:8765` of the machine it runs on. To open it from your laptop, forward the port
over SSH, in a second terminal on your laptop:

```bash
# app running on a login node (note which one, e.g. eu-login-12):
ssh -L 8765:localhost:8765 <user>@eu-login-12.euler.ethz.ch
# app running on a compute node inside a job (find the node with `squeue -u $USER`, e.g. eu-a2p-123):
ssh -L 8765:eu-a2p-123:8765 <user>@euler.ethz.ch
```

then open <http://localhost:8765/>. If port 8765 is taken, choose another on both sides. The pages:

| Page | Shows |
|---|---|
| `/` | search box, index status, orthogroups and cases |
| `/search?q=CDK2` | search by UniProt accession, gene name or protein name |
| `/protein/<accession>` | structure and pocket table; click a pocket to highlight its lining residues in 3D |
| `/og` and `/og/<id>` | gene families, and one family as a SOG × species state matrix on the species tree |
| `/sog` | all site orthogroups, filterable by event type and score |
| `/case` and `/case/<id>` | curated cases: evidence by filter, artifact checks, co-folding table, ancestral vs derived 3D panels |
| `/validation` | the seven known literature events and negative controls, traced through the pipeline |
| `/api/docs` | the JSON API |

The header has a selector for the two orthology versions (O1_r2, O2_r1). Tier-3 proteins need a node-local index:
`uv run --project . python -m webapp build-index --tiers 3` (about 2 min); without it the app serves tiers 1–2.
More detail (endpoints, environment variables, caching, latency): [`webapp/README.md`](webapp/README.md). Note: the
3D viewer was tested through its API but never opened in a browser by the agents (see the briefing).

### Rebuilding the data from scratch

On a cluster with the atlas and AFDB mounted, `scripts/run_pipeline.py` regenerates every stage on a snapshot of
the atlas (print the plan with `--dry-run` first); `python -m webapp build-index` and `build-og-index` then build the
browser indexes. The full run took several hours on 64 CPUs and 4 GPUs; runtimes per stage are in
[`record/REPORT.md`](record/REPORT.md) §1.

```bash
.venv/bin/python scripts/run_pipeline.py --snapshot new --run-id myrun --dry-run
.venv/bin/python scripts/run_pipeline.py --snapshot new --run-id myrun
```

## Checking the reported numbers and rebuilding the briefing

```bash
.venv/bin/python scripts/check_report_numbers.py      # every number in record/REPORT.md vs results/numbers.json
.venv/bin/python briefing/make_figures.py             # briefing figures, from committed results only
uv pip install --python .venv/bin/python typst && .venv/bin/python -c "import typst; typst.compile('briefing/briefing.typ', output='briefing/briefing.pdf')"
```

`briefing/build_bib.py` regenerates the bibliography from Crossref (needs network).

## Tools and data this work builds on

- **BURROW** pocket finder: <https://github.com/zoidpuff/burrow>
- **pocketcompare** pocket fingerprints and the AFDB pocket atlas: <https://github.com/zoidpuff/pocketcompare>
- **Ancestra** hierarchical orthology: <https://github.com/evocellnet/ancestra-mvp-ai>
- AlphaFold DB, AlphaFold 3, BioLiP2, UniProt, M-CSA, LIGYSIS, PANTHER, TED, ChEMBL and others, with versions,
  checksums and licences in [`DATA.md`](DATA.md). The literature that shaped the methods is credited in the
  briefing (section 10) and reviewed in [`record/LITERATURE.md`](record/LITERATURE.md).
