#!/usr/bin/env python
"""Duplication-associated gains (exploratory, NOT preregistered; PI follow-up after D28).

For primary gain calls on a paralog branch (parent = duplication node) that fail orthology_split, compare the
derived paralog copy D with its sister copy S (the other children of the duplication node) species by species:
  pattern   species s with >= 1 P member in D and, in S, >= 1 A member and no P/W member (A = primary pLDDT >= 90,
            D11, after rescue)
  swapped   the same test with the paralog labels swapped (S has P, D has only A) = the label-swap null
Category "duplication-associated gain (exploratory, not preregistered)" if n_pattern >= 3 and n_pattern > n_swapped.
Scope: candidate gains only.
Writes results/partB/<snap>/<run>/duplication_gains_exploratory.parquet + DUPLICATION_GAINS.md; candidates.parquet is
not modified.

    python scripts/partB/dup_gain_test.py --sog-run O1_r3
"""
from __future__ import annotations

import argparse
import glob
import sys
import time
from pathlib import Path

import pandas as pd
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pocketevo import recon  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()
LABEL = "duplication-associated gain (exploratory, not preregistered)"


def species_states(t, lo, hi, stt, sog):
    out = {}
    for k in range(lo, hi):
        x = t.node_name[t.leaves[k]]
        out.setdefault(t.species[t.leaves[k]], []).append(stt.get((sog, x), "U"))
    return out


def main():
    ap = argparse.ArgumentParser()
    recon.add_partb_args(ap, sog_run="O1_r3", skip=("orthology",))
    ap.add_argument("--orthology", default=None, choices=["O1", "O2"])
    ap.add_argument("--min-species", type=int, default=3)
    a = ap.parse_args()
    orth = a.orthology or ("O2" if a.sog_run.startswith("O2") else "O1")
    out = ROOT / f"results/partB/{SNAP}/{a.sog_run}"
    cand = pd.read_parquet(out / "candidates.parquet")
    c = cand[cand["primary"].astype(bool) & (cand["type"] == "gain") & cand["flag_paralog_branch"].astype(bool) &
             (cand["orthology_split"] == "fail")].copy()
    ox = ROOT / ("data/orthology/ancestra_O2/accurate/hogs.orthoxml.gz" if orth == "O2"
                 else "data/orthology/ancestra_tier1_v0/accurate/hogs.orthoxml.gz")
    trees = {t.tree_id: t for t in recon.iter_orthoxml_trees(ox, keep=set(c["og_id"]))}
    stt = {}
    for tau in sorted(set(c["tau"])):
        if orth == "O2":
            st = ds.dataset(sorted(glob.glob(f"/tmp/pe_big/partB/{SNAP}/{a.sog_run}/states_resc_t{tau}/shard=*.parquet")),
                            format="parquet").to_table(columns=["sog_id", "accession", "state"]).to_pandas()
        else:
            st = pd.read_parquet(ROOT / f"results/partB/events/{SNAP}/{a.sog_run}/tau{tau}/primary/rescue/sog_states/"
                                        "bucket=all.parquet", columns=["sog_id", "accession", "state"])
        st = st[st["sog_id"].isin(set(c["sog_id"]))]
        stt[tau] = dict(zip(zip(st["sog_id"], st["accession"]), st["state"]))
    rows = []
    for r in c.itertuples():
        t = trees[r.og_id]
        node = int(r.node)
        par = int(t.parent[node])
        sd = species_states(t, t.leaf_lo[node], t.leaf_hi[node], stt[r.tau], r.sog_id)
        ss = {}
        for ch in t.children(par):
            if int(ch) == node:
                continue
            for sp, v in species_states(t, t.leaf_lo[ch], t.leaf_hi[ch], stt[r.tau], r.sog_id).items():
                ss.setdefault(sp, []).extend(v)
        common = set(sd) & set(ss)

        def pat(x, y):  # x has P, y has only A among informative (>= 1 A, no P/W)
            return sorted(s for s in common if "P" in x[s] and "A" in y[s] and not ({"P", "W"} & set(y[s])))
        pp, sw = pat(sd, ss), pat(ss, sd)
        rows.append(dict(sog_id=r.sog_id, og_id=r.og_id, tau=r.tau, node=node, level=r.level,
                         duplication_level=t.level[par], n_species_both_copies=len(common),
                         n_species_pattern=len(pp), n_species_swapped=len(sw), species_pattern=";".join(pp),
                         species_swapped=";".join(sw),
                         duplication_associated_gain=bool(len(pp) >= a.min_species and len(pp) > len(sw)),
                         label=LABEL, pass_artifacts=r.pass_artifacts, artifact_fails=r.artifact_fails))
    df = pd.DataFrame(rows)
    # background (selection-free) null: every other duplication node child of the same SOG trees, both orientations
    bg = []
    for r in c.drop_duplicates(["sog_id", "tau"]).itertuples():
        t = trees[r.og_id]
        for par in range(t.n_nodes):
            if t.kind[par] != recon.DUP:
                continue
            for ch in t.children(par):
                ch = int(ch)
                if ch == int(r.node):
                    continue
                sd = species_states(t, t.leaf_lo[ch], t.leaf_hi[ch], stt[r.tau], r.sog_id)
                ss = {}
                for o in t.children(par):
                    if int(o) != ch:
                        for sp, v in species_states(t, t.leaf_lo[o], t.leaf_hi[o], stt[r.tau], r.sog_id).items():
                            ss.setdefault(sp, []).extend(v)
                cm = set(sd) & set(ss)
                n = sum(1 for sp in cm if "P" in sd[sp] and "A" in ss[sp] and not ({"P", "W"} & set(ss[sp])))
                bg.append(n)
    bg = pd.Series(bg, dtype=float)
    f = out / "duplication_gains_exploratory.parquet"
    df.to_parquet(f, index=False)
    n = int(df["duplication_associated_gain"].sum()) if len(df) else 0
    L = [f"# Duplication-associated gains — {a.sog_run}, snapshot {SNAP} (EXPLORATORY, NOT PREREGISTERED)\n",
         f"Generated {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())} by `scripts/partB/dup_gain_test.py`; git "
         f"`{recon.git_commit()}`. Scope: primary gains on a paralog branch that fail `orthology_split` "
         f"({len(df)} of {int((cand['primary'].astype(bool) & (cand['type'] == 'gain')).sum())} primary gains).\n",
         "Test: species with the pocket (P) in the derived paralog copy and an observed absence (A only; pLDDT >= 90, "
         f"D11, after rescue) in the sister copy; category if >= {a.min_species} such species and more than under "
         "the label-swap null (sister P, derived A only). No FDR control; labelled "
         f"\"{LABEL}\".\n",
         f"**{n} of {len(df)} pass.** Label-swap null: total pattern species {int(df['n_species_pattern'].sum()) if len(df) else 0} "
         f"vs swapped {int(df['n_species_swapped'].sum()) if len(df) else 0}; SOGs with >= {a.min_species} swapped "
         f"species: {int((df['n_species_swapped'] >= a.min_species).sum()) if len(df) else 0}.\n",
         f"Background null (no selection): all other duplication-node copies of the same SOG trees, both "
         f"orientations: {len(bg)} copies, {int((bg >= a.min_species).sum()) if len(bg) else 0} "
         f"({(bg >= a.min_species).mean() if len(bg) else float('nan'):.3f}) reach >= {a.min_species} pattern species. "
         "The label-swap null is conditioned on the gain selection (derived copy P, outgroups A) and is therefore "
         "weak; read the background rate instead.\n",
         "| sog_id | level | duplication_level | n_species_both | n_pattern | n_swapped | dup-associated gain | pass_artifacts |",
         "|---|---|---|---|---|---|---|---|"]
    for r in df.sort_values("n_species_pattern", ascending=False).itertuples():
        L.append(f"| {r.sog_id} | {r.level} | {r.duplication_level} | {r.n_species_both_copies} | {r.n_species_pattern} | "
                 f"{r.n_species_swapped} | {r.duplication_associated_gain} | {r.pass_artifacts} |")
    (out / "DUPLICATION_GAINS.md").write_text("\n".join(L) + "\n")
    for x in (f, out / "DUPLICATION_GAINS.md"):
        recon.write_sidecar(x, SNAP, 0, {"script": "scripts/partB/dup_gain_test.py", "argv": sys.argv[1:]})
    print(n, len(df), "background", len(bg), int((bg >= a.min_species).sum()) if len(bg) else 0, df[["n_species_pattern", "n_species_swapped"]].sum().to_dict() if len(df) else {})


if __name__ == "__main__":
    main()
