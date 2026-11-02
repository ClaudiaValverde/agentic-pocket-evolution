"""Make the briefing figures from committed results only (no cluster data needed).

    .venv/bin/python briefing/make_figures.py        # writes briefing/figures/*.svg and *.png

Every number drawn here is read from results/numbers.json or from a committed JSON/CSV file named next to the
figure function; nothing is typed in by hand except the pipeline schematic (Fig. 1) and the gate times of Fig. 7,
which come from record/RESEARCH_TREE.md.
"""
import csv
import json
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

R = Path(__file__).resolve().parents[1]  # repo root
OUT = R / "briefing/figures"
SNAP = "20261030T085136Z"
NUM = json.loads((R / "results/numbers.json").read_text())

# Reference palette (dataviz skill, light mode): categorical slots 1-3, plus neutral inks.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
BLUE_LIGHT = "#9ec5f4"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "legend.frameon": False, "svg.fonttype": "path",
})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / f"{name}.svg")


def panel(ax, letter):
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK,
            ha="right", va="bottom")


def hgrid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------------------------------------------------------
def fig1_pipeline():
    """Schematic of the pipeline (boxes and arrows only; counts from REPORT.md §1)."""
    fig, ax = plt.subplots(figsize=(8.6, 4.3))
    ax.set_xlim(0, 100), ax.set_ylim(0, 50), ax.axis("off")

    def box(x, y, w, h, title, sub, fc="#f3f6fb", ec=BLUE):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=0.9))
        ax.text(x + w / 2, y + h - 2.2, title, ha="center", va="top", fontsize=7.6, fontweight="bold", color=INK)
        ax.text(x + w / 2, y + h - 5.6, sub, ha="center", va="top", fontsize=6.4, color=INK2, linespacing=1.25)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9,
                                                              shrinkA=0, shrinkB=0, mutation_scale=8))

    # inputs
    box(1, 35, 20, 13, "Pocket atlas", "BURROW pockets on\nAlphaFold DB models\n+ pocketcompare vectors")
    box(1, 18, 20, 13, "Orthology", "Ancestra HOGs\nO1: 46 proteomes\nO2: 157 proteomes")
    box(1, 1, 20, 13, "Known answers", "BioLiP2, UniProt, M-CSA,\nLIGYSIS sites; 7 literature\nevents; holo PDB")
    # core
    box(27, 26, 21, 22, "Site orthogroups", "pockets mapped onto a\nper-family alignment;\nmatched if columns overlap\n(Jaccard >= 0.2), grouped\nby community detection")
    box(27, 1, 21, 20, "Observation states", "present / weak / absent /\nuninformative / no gene;\nabsence only in confident\nmodels; error rates measured\non known sites")
    # analyses
    box(54, 30, 21, 18, "Part A: global", "how much pocket similarity\nis explained by homology?\ndo pockets follow the\ngene tree?")
    box(54, 1, 21, 25, "Part B: events", "gain / loss / specificity\non gene trees (parsimony)\n-> null-input pipelines\n-> artifact cascade\n-> 3 evidence filters")
    box(80, 18, 19, 30, "Co-folding (AF3)", "calibrated on 24 holo\npositives + 72 decoys and\nsite ablations; refold\nboth states with and\nwithout the ligand")
    box(80, 1, 19, 14, "Web app", "proteins, SOG x species\nmatrix, case pages", fc="#f4faf7", ec=AQUA)

    arrow(21.5, 41, 26.6, 38), arrow(21.5, 24, 26.6, 33), arrow(21.5, 8, 26.6, 10)
    arrow(37.5, 25.6, 37.5, 21.5)
    arrow(48.5, 40, 53.6, 40), arrow(48.5, 12, 53.6, 12), arrow(48.5, 30, 53.6, 20)
    arrow(75.5, 15, 79.6, 30), arrow(75.5, 8, 79.6, 8)
    save(fig, "fig1_pipeline")


# ---------------------------------------------------------------------------------------------------------------
def fig2_absence():
    """False-absence rates on known binding sites (numbers.json observation.false_A, observation.rescue;
    results/observation/<snap>/rescue_breakdown.csv for the > 10-residue row)."""
    fa = NUM["observation"]["false_A"]
    tune_k1 = next(r for r in NUM["observation"]["rescue"]["tune"] if r["k"] == 1)
    with open(R / f"results/observation/{SNAP}/rescue_breakdown.csv") as fh:
        gt10 = next(float(r["false_A_after"]) for r in csv.DictReader(fh)
                    if r["criterion"] == "A_primary" and r["half"] == "tune" and r["strata"] == "site_size_class"
                    and r["k"] == "1" and r["site_size_class"] == ">10")

    labels = ["metal ions", "small molecules", "cofactors", "all known sites",
              "all, after re-detection\nrescue", "> 10-residue sites,\nafter rescue"]
    vals = [fa["by_ligand_class"]["metal"]["false_A"], fa["by_ligand_class"]["small_molecule"]["false_A"],
            fa["by_ligand_class"]["cofactor"]["false_A"], fa["A_primary"]["false_A"],
            tune_k1["false_A_after"], gt10]
    colors = [BLUE_LIGHT] * 3 + [BLUE] + [AQUA] * 2
    fig, ax = plt.subplots(figsize=(6.2, 2.5))
    y = list(range(len(labels)))[::-1]
    ax.barh(y, vals, color=colors, height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v + 0.012, yi, f"{v:.2f}", va="center", fontsize=7.5, color=INK)
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 0.72)
    ax.set_xlabel("fraction of known binding sites called 'absent' in a confident model")
    ax.axvline(0.05, color=ORANGE, lw=1, ls="--")
    ax.text(0.055, y[0] + 0.45, "preregistered target 0.05", color=ORANGE, fontsize=7, va="bottom")
    hgrid(ax, "x")
    ax.set_title(f"How often a real pocket is missed (n = {fa['n_known']:,} experimentally known sites)",
                 loc="left", color=INK)
    save(fig, "fig2_absence")


# ---------------------------------------------------------------------------------------------------------------
def fig3_partA():
    """(a) per-class significance rates, O1 (numbers.json partA.per_class r_adj);
    (b) F_analog under two weightings and two orthologies (partA / partA_O2 F_* and strata.both_plddt_ge90)."""
    pc = NUM["partA"]["per_class"]
    cls = [("H1s", "orthologs, same site"), ("H2s", "paralogs, same site"), ("H1x", "orthologs, different site"),
           ("H3j", "shared domain family"), ("H4", "remote structural only"), ("N", "no detectable homology")]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={"width_ratios": [1.1, 1]})
    y = list(range(len(cls)))[::-1]
    for yi, (k, lab) in zip(y, cls):
        v, (lo, hi) = pc[k]["r_adj"]["value"], pc[k]["r_adj"]["ci95"]
        c = BLUE if k != "N" else ORANGE
        a.plot([lo, hi], [yi, yi], color=c, lw=1.6, solid_capstyle="round")
        a.plot(v, yi, "o", color=c, ms=5, mec="white", mew=0.8)
    a.set_xscale("log")
    a.set_yticks(y, [lab for _, lab in cls])
    a.set_xlabel("P(pocket pair is significantly similar)")
    a.set_xlim(2e-7, 0.08)
    hgrid(a, "x")
    a.set_title("Similarity rate by relationship (O1)", loc="left", color=INK)
    panel(a, "a")

    rows = [("O1, all", NUM["partA"], None), ("O2, all", NUM["partA_O2"], None),
            ("O1, pLDDT >= 90", NUM["partA"], "both_plddt_ge90"), ("O2, pLDDT >= 90", NUM["partA_O2"], "both_plddt_ge90")]
    yy = list(range(len(rows)))[::-1]
    for yi, (lab, src, stratum) in zip(yy, rows):
        for key, col, off in (("per_family_pair", ORANGE, 0.14), ("per_pair", BLUE, -0.14)):
            if stratum is None:
                d = src["F_headline_per_family_pair" if key == "per_family_pair" else "F_per_pair"]["F_analog"]
            else:
                d = src["strata"][stratum][f"{key}.F_analog"]
            v, (lo, hi) = d["value"], d["ci95"]
            b.plot([lo, hi], [yi + off] * 2, color=col, lw=1.6)
            b.plot(v, yi + off, "o", color=col, ms=5, mec="white", mew=0.8)
            b.text(hi + 0.02, yi + off, f"{v:.2f}", va="center", fontsize=7, color=INK)
    b.set_yticks(yy, [r[0] for r in rows])
    b.set_xlim(0, 1.08)
    b.set_xlabel("F_analog: share of similar pairs with no detectable homology")
    hgrid(b, "x")
    b.plot([], [], "o-", color=ORANGE, label="each family pair counts once (preregistered)")
    b.plot([], [], "o-", color=BLUE, label="each pocket pair counts once")
    b.legend(loc="upper center", bbox_to_anchor=(0.45, -0.25), fontsize=7, ncol=1)
    b.set_title("Upper bound on convergence", loc="left", color=INK)
    panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig3_partA")


# ---------------------------------------------------------------------------------------------------------------
def fig4_descent():
    """Pocket retention (M1) vs sequence identity, orthologs vs paralogs, O1 (numbers.json
    partA_descent.M1_qm90_cap_by_identity)."""
    d = NUM["partA_descent"]["M1_qm90_cap_by_identity"]
    bins = [k.split("|")[1] for k in d if k.startswith("ortholog|")]
    x = [float(b.split("-")[0]) + 0.05 for b in bins]
    fig, ax = plt.subplots(figsize=(5.2, 2.9))
    for kind, col, dx in (("ortholog", BLUE, -0.006), ("paralog", ORANGE, 0.006)):
        m = [d[f"{kind}|{b}"]["m1"] for b in bins]
        lo = [d[f"{kind}|{b}"]["ci"][0] for b in bins]
        hi = [d[f"{kind}|{b}"]["ci"][1] for b in bins]
        xs = [xi + dx for xi in x]
        ax.vlines(xs, lo, hi, color=col, lw=1.4)
        ax.plot(xs, m, "-o", color=col, lw=1.6, ms=4.5, mec="white", mew=0.8, label=f"{kind} pairs")
    ax.set_xlabel("sequence identity between the two proteins")
    ax.set_ylabel("P(site kept | testable)")
    ax.set_xticks([i / 10 for i in range(0, 11, 2)], [f"{i * 10}%" for i in range(0, 11, 2)])
    ax.set_ylim(0, 1)
    hgrid(ax)
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("Pocket retention rises with sequence identity; orthologs = paralogs", loc="left", color=INK)
    save(fig, "fig4_descent")


# ---------------------------------------------------------------------------------------------------------------
def fig5_partB():
    """(a) real vs null-input call counts (numbers.json partB.O{1,2}_funnel_and_null);
    (b) evidence funnel per category (results/cases/<snap>/case_list_v2.parquet.json counts)."""
    pb = NUM["partB"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1, 1.25]})
    groups = [("O1 gains", pb["O1_funnel_and_null"], "gain"), ("O1 losses", pb["O1_funnel_and_null"], "loss"),
              ("O2 gains", pb["O2_funnel_and_null"], "gain"), ("O2 losses", pb["O2_funnel_and_null"], "loss")]
    w = 0.26
    for i, (lab, src, t) in enumerate(groups):
        e = src["empirical_null_FDR"]
        vals = [e[f"{t}_null_curveball"]["real"], e[f"{t}_null_curveball"]["null_mean"],
                e[f"{t}_null_pseudo_clade"]["null_mean"]]
        for j, (v, col) in enumerate(zip(vals, (BLUE, MUTED, "#c3c2b7"))):
            a.bar(i + (j - 1) * w, v, w * 0.92, color=col)
    a.set_xticks(range(len(groups)), [g[0].replace(" ", "\n") for g in groups])
    a.set_ylabel("number of event calls")
    for lab, col in (("real data", BLUE), ("null: shuffled states", MUTED), ("null: random pseudo-clades", "#c3c2b7")):
        a.bar(0, 0, color=col, label=lab)
    a.legend(fontsize=7, loc="upper right")
    hgrid(a)
    a.set_title("Calls on real data vs scrambled data", loc="left", color=INK)
    panel(a, "a")

    cl = json.loads((R / f"results/cases/{SNAP}/case_list_v2.parquet.json").read_text())["counts"]
    cats = [("O1:gain", "O1 gains"), ("O1:loss", "O1 losses"), ("O2:gain", "O2 gains"), ("O2:loss", "O2 losses")]
    stages = [("n", "candidates"), ("supported_or_convincing", "persistence + 1 other filter"),
              ("artifacts_pass", "pass artifact checks"), ("convincing_headline", "convincing (all checks)")]
    cols = ["#c3c2b7", BLUE, BLUE_LIGHT, ORANGE]
    w = 0.2
    for i, (k, _) in enumerate(cats):
        for j, ((s, _), col) in enumerate(zip(stages, cols)):
            v = cl[k][s]
            b.bar(i + (j - 1.5) * w, v, w * 0.92, color=col)
            if s == "convincing_headline":
                b.text(i + (j - 1.5) * w, v + 1.5, str(v), ha="center", fontsize=7, color=ORANGE, fontweight="bold")
    for (s, lab), col in zip(stages, cols):
        b.bar(0, 0, color=col, label=lab)
    b.set_xticks(range(len(cats)), [c[1].replace(" ", "\n") for c in cats])
    b.set_ylabel("cases")
    b.set_ylim(0, 150)
    b.legend(fontsize=7, loc="upper left", ncol=2, columnspacing=1.0)
    hgrid(b)
    b.set_title("Candidates and the checks they pass", loc="left", color=INK)
    panel(b, "b")
    fig.tight_layout(w_pad=2)
    save(fig, "fig5_partB")


# ---------------------------------------------------------------------------------------------------------------
def fig6_cofold():
    """(a) calibration AUCs (numbers.json cofold.calibration.single_readout_auc);
    (b) ligand-free refold site recall, present vs absent side (numbers.json cofold.apo_pool)."""
    cal = NUM["cofold"]["calibration"]
    auc = cal["single_readout_auc"]
    names = [("lig_plddt", "ligand confidence (pLDDT)"), ("site_recall", "site recall"), ("pae_min", "min. ligand PAE"),
             ("centroid", "ligand-site distance")]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw={"width_ratios": [1, 1.25]})
    y = list(range(len(names)))[::-1]
    for yi, (k, lab) in zip(y, names):
        a.plot([auc[k]["ci_lo"], auc[k]["ci_hi"]], [yi, yi], color=BLUE, lw=1.6)
        a.plot(auc[k]["auc"], yi, "o", color=BLUE, ms=5, mec="white", mew=0.8)
        a.text(auc[k]["ci_hi"] + 0.015, yi, f"{auc[k]['auc']:.2f}", va="center", fontsize=7)
    m = cal["auc_by_ligand_class"]["metal"]["lig_plddt"]
    a.plot([m["ci_lo"], m["ci_hi"]], [-1, -1], color=ORANGE, lw=1.6)
    a.plot(m["auc"], -1, "o", color=ORANGE, ms=5, mec="white", mew=0.8)
    a.text(m["ci_hi"] + 0.015, -1, f"{m['auc']:.2f}", va="center", fontsize=7)
    a.set_yticks(y + [-1], [n[1] for n in names] + ["ligand pLDDT, metals only"])
    a.axvline(0.5, color=MUTED, ls="--", lw=1)
    a.set_xlim(0.25, 1.08)
    a.set_xlabel("AUC: true binders vs decoys / ablated sites")
    hgrid(a, "x")
    n = cal["panel"]
    a.set_title(f"Calibration ({n['n_positives']} holo + {n['n_negatives']} negatives)", loc="left", color=INK)
    panel(a, "a")

    ap = {k: v for k, v in NUM["cofold"]["apo_pool"].items() if not k.startswith("_")}
    order = sorted(ap, key=lambda k: (0 if k.startswith("gain") else 1 if k.startswith("dup") else 2, k))
    order = [k for k in order if ap[k]["absent_recall_median"] is not None]
    labels = []
    for i, k in enumerate(order):
        p, q = ap[k]["present_recall_median"], ap[k]["absent_recall_median"]
        col = BLUE if k.startswith(("gain", "dup")) else ORANGE
        b.plot([i, i], [q, p], color=col, lw=1.4, alpha=0.6)
        b.plot(i, p, "o", color=col, ms=5.5, mec="white", mew=0.8)
        b.plot(i, q, "o", mfc="white", mec=col, ms=5.5, mew=1.3)
        labels.append(k.split("_")[2] if k.startswith("control") else k.split("_")[1])
    b.set_xticks(range(len(order)), [lab.replace("HOG00", "") for lab in labels], fontsize=6, rotation=45, ha="right")
    ng = sum(k.startswith(("gain", "dup")) for k in order)
    b.text((ng - 1) / 2, -0.23, "gains", ha="center", fontsize=7, color=BLUE,
           transform=b.get_xaxis_transform())
    b.text((ng + len(order) - 1) / 2, -0.23, "false-loss controls", ha="center", fontsize=7, color=ORANGE,
           transform=b.get_xaxis_transform())
    b.set_ylim(-0.05, 1.08)
    b.set_ylabel("site recall in a fresh\nAF3 model (median of 3)")
    b.plot([], [], "o", color=MUTED, label="side where the pocket was seen")
    b.plot([], [], "o", mfc="white", mec=MUTED, mew=1.3, label="side where it was missing")
    b.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.27), ncol=2)
    hgrid(b)
    b.set_title("Ligand-free refolds of both states", loc="left", color=INK)
    panel(b, "b")
    fig.subplots_adjust(wspace=0.8)
    save(fig, "fig6_cofold")


# ---------------------------------------------------------------------------------------------------------------
def fig7_process():
    """Commits per hour over the run (record/commit_timeline.tsv), with gates from record/RESEARCH_TREE.md.
    Times are shown in CET (UTC+1), the clock the project ran on."""
    rows = list(csv.DictReader(open(R / "record/commit_timeline.tsv"), delimiter="\t"))
    t = [datetime.strptime(r["commit_time_utc"], "%Y-%m-%dT%H:%M:%SZ") + timedelta(hours=1) for r in rows]
    t0 = datetime(2026, 10, 30, 9, 0)
    hours = [(ti - t0).total_seconds() / 3600 for ti in t]
    nb = int(max(hours)) + 2
    counts = [0] * nb
    for h in hours:
        counts[int(h)] += 1
    fig, ax = plt.subplots(figsize=(7.0, 2.6))
    ax.bar([i + 0.5 for i in range(nb)], counts, width=0.85, color=BLUE)

    def h(day, hh, mm=0):
        return (datetime(2026, 10, day, hh, mm) - t0).total_seconds() / 3600

    ax.axvspan(h(30, 18, 7), h(30, 22, 56), color=ORANGE, alpha=0.12, lw=0)
    ax.axvspan(h(30, 23, 58), h(31, 20, 39), color=ORANGE, alpha=0.22, lw=0)
    ax.text((h(30, 18, 7) + h(30, 22, 56)) / 2, max(counts) * 0.8, "compute\nstalled", ha="center", fontsize=6.5,
            color=INK2, va="top")
    ax.text((h(30, 23, 58) + h(31, 20, 39)) / 2, max(counts) * 0.8,
            "lead agent idle (no wake-up scheduled)\nunattended snapshot-2 rerun 02:00-05:20 failed in 3 stages",
            ha="center", fontsize=6.5, color=INK2, va="top")
    for i, (lab, x) in enumerate((("design frozen", h(30, 10, 47)), ("first results", h(30, 12, 13)),
                                  ("cases frozen", h(30, 18, 30)), ("deadline", h(31, 9, 45)),
                                  ("report final", h(31, 22, 4)))):
        ax.axvline(x, color=INK2, lw=0.8, ls=":")
        ax.text(x, max(counts) * (1.08 if i % 2 == 0 else 1.0), lab, ha="center", fontsize=6.5, color=INK2,
                bbox=dict(fc="white", ec="none", pad=0.6))
    ticks = list(range(0, nb + 1, 3))
    ax.set_xticks(ticks, [(t0 + timedelta(hours=k)).strftime("%a %H:%M") for k in ticks], fontsize=6.5, rotation=0)
    ax.set_xlim(0, nb)
    ax.set_ylim(0, max(counts) * 1.2)
    ax.set_ylabel("commits per hour")
    hgrid(ax)
    ax.set_title(f"Activity of the agent team over the run ({len(rows)} commits)", loc="left", color=INK, pad=14)
    save(fig, "fig7_process")


if __name__ == "__main__":
    fig1_pipeline()
    fig2_absence()
    fig3_partA()
    fig4_descent()
    fig5_partB()
    fig6_cofold()
    fig7_process()
