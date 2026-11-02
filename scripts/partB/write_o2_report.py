#!/usr/bin/env python
"""O2 Part B report: results/partB/<snap>/O2_r1/{candidates.parquet, CANDIDATES.md} from run_events.py (O2 states
with the rescue overlay, prep_o2_states.py), run_nulls.py (>= 20 seeds) and o1_o2_agreement.py."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pocketevo import recon  # noqa: E402
from write_candidates import md  # noqa: E402

SNAP = recon.cli_value("--snapshot") or recon.default_snapshot()   # --snapshot / $PE_SNAPSHOT / prereg


def main():
    ap = argparse.ArgumentParser()
    recon.add_partb_args(ap, orthology="O2", sog_run="O2_r1")
    ap.add_argument("--o1-run", default="O1_r3")
    a = ap.parse_args()
    R1, R2 = a.o1_run, a.sog_run
    out = ROOT / f"results/partB/{SNAP}/{R2}"
    ev = ROOT / f"results/partB/events/{SNAP}/{R2}/tau0.2/primary"
    th = json.loads((ROOT / "results/partB/sim/measured/d11_primary_known_gt10/thresholds.json").read_text())
    e = pd.read_parquet(ev / "events.parquet")
    c = e[e["called"]].copy()
    sup = {(k, t): bool(x and x["recall"] >= 1e-4 and x["fdr"] <= th["event_fdr_target"])
           for k, v in th["chosen"].items() for t, x in v.items()}
    c["fdr_supported"] = [sup.get((f"{int(g)}:{int(l)}", t), False) for g, l, t in zip(c["gain_cost"], c["loss_cost"], c["type"])]
    c["label"] = np.where(c["fdr_supported"], "FDR-calibrated candidate", "not FDR-calibrated (audit)")
    c["primary"] = (c["gain_cost"] == 2) & (c["loss_cost"] == 1)
    c = c.rename(columns={"char": "sog_id"})
    fun = pd.read_csv(ev / "funnel.csv")
    nul = pd.read_csv(ROOT / f"results/partB/nulls/{SNAP}/{R2}/tau0.2/primary/funnel_real_vs_null.csv")
    rows = []
    for typ in ("gain", "loss"):
        real = nul[(nul["input"] == "real") & (nul["type"] == typ)]["called"].iloc[0]
        for inp in ("null_curveball", "null_pseudo_clade"):
            q = nul[(nul["input"] == inp) & (nul["type"] == typ)]
            rows.append(dict(type=typ, null=inp, real_called=int(real), null_called_mean=float(q["called"].mean()),
                             null_called_max=int(q["called"].max()), n_null_seeds=len(q),
                             empirical_null_FDR=min(1.0, float(q["called"].mean()) / real) if real else None,
                             simulated_FDR=(th["chosen"]["2:1"][typ]["fdr"] if th["chosen"]["2:1"].get(typ) and
                                            th["chosen"]["2:1"][typ]["recall"] >= 1e-4 else None)))
            c.loc[c["type"] == typ, f"{inp}_called_mean_same_cut"] = float(q["called"].mean())
    nc = pd.DataFrame(rows)
    f = out / "candidates.parquet"
    c.to_parquet(f, index=False)
    prep = json.loads((out / "o2_states_prep.json").read_text())
    agf = ROOT / f"results/partB/{SNAP}/{R1}/o1_o2_gain_agreement.json"
    agr = json.loads(agf.read_text()) if agf.exists() else None   # O1 lane of this snapshot may not have run yet
    p = c[c["primary"]]
    L = [f"# Part B candidates — O2 (O2_r1 + rescue overlay), snapshot {SNAP} (development subset, preliminary)\n",
         f"Generated {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())} by `scripts/partB/write_o2_report.py`; prereg "
         f"sha256 `{recon.prereg_sha256()}`; git `{recon.git_commit()}`. Same code, cuts (O1 measured-noise thresholds) "
         "and nulls as O1_r3. Not yet verified.\n",
         "## States\n",
         f"{R2} states τ 0.2 (D8 flags excluded, D11 applied by the SOG stage) + rescue overlays {prep.get('overlays')} "
         f"(k = 1): A in the rescue pool → W {prep['A_overlay_W']} / A {prep['A_overlay_A']}; `rescue:pending` A "
         f"outside the pool → U ({prep['A_pending_to_U']}); A without atlas-filtered pockets kept as A "
         f"({prep['A_not_pending_kept']}, not re-described). P calls are not rescued on O2 (overlay is A-only).\n",
         "## Funnel (real)\n", md(fun) + "\n",
         "## Real vs null at the same cut (≥ 20 seeds; empirical null FDR = mean null calls / real calls)\n",
         md(nc) + "\n",
         "## Calls (2:1)\n",
         f"**{int((p['type'] == 'gain').sum())} gains (FDR-calibrated cut) + {int((p['type'] == 'loss').sum())} losses "
         "(not FDR-calibrated, D21).** Yield preflight: ≥ 5 → secondary stratum not triggered.\n",
         md(p.sort_values("score_S_partial", ascending=False).head(20),
            ["type", "label", "og_id", "sog_id", "level", "clade_lca_species_tree", "n_species_P", "n_species_A",
             "out_A_lineages", "f1_single_event", "score_S_partial"]) + "\n",
         "## O1 → O2 agreement for gains (misassigned-orthology check)\n"]
    if agr is None:
        L += [f"Not computed yet: {agf.relative_to(ROOT)} is missing (run_partB.sh --orthology O2 --stages "
              f"agreement,report --o1-run {R1} after Part B O1 of this snapshot).\n"]
    else:
        L += [f"{agr['n_O1_gains']} {R1} primary gains: re-called on O2 {agr['recovered_called']}, any O2 MPR gain event "
         f"{agr['recovered_event']}, O1 derived members in one O2 SOG {agr['o2_same_sog']}. "
         f"{agr['diagnosis']['reading']}\n",
         f"Of the {agr['diagnosis']['O1_A_members_total']} outgroup A members behind the O1 gains, "
         f"{agr['diagnosis']['O1_A_members_in_another_O2_root_HOG']} are in a different O2 root HOG than the derived "
         "clade (O2 does not treat them as members of the same family).\n",
         "```json\n" + json.dumps({k: v for k, v in agr["diagnosis"].items() if k != "reading"}, indent=1) + "\n```\n"]
    L += [
         "## Not done on O2\n- Specificity: tier-3/L sites were built without fingerprints (`--no-fingerprints`).\n"
         "- Artifact cascade and P rescue on O2 members; τ 0.3 sensitivity on O2.\n"]
    (out / "CANDIDATES.md").write_text("\n".join(L))
    for x in (f, out / "CANDIDATES.md"):
        recon.write_sidecar(x, SNAP, 0, {"script": "scripts/partB/write_o2_report.py"})
    print(out / "CANDIDATES.md")


if __name__ == "__main__":
    main()
