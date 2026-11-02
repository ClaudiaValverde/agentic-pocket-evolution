#!/usr/bin/env python
"""G5 status + snapshot-2 vs snapshot-1 comparison -> results/pipeline/G5_STATUS.md.

  --status-only : stage table from /tmp/pe_big/logs/g5/<snap2>/*.jsonl (called after every stage by scripts/g5_launch.sh)
  default       : stage table + comparison (SOG counts, §8.1 F1, D20 clean-test F1, O2 rescue, Part A headline numbers,
                  O1/O2 gain/loss after rescue and after the artifact cascade with overlap, O1 reproducibility check)
  --notes FILE  : markdown appended under "Reproducibility notes" (explanations of O1 / Part A differences)

    .venv/bin/python scripts/g5_compare.py --s2 <snap2> --run-id g5 [--status-only]
"""
import argparse
import glob
import json
import subprocess
from pathlib import Path

import pandas as pd

P = Path(__file__).resolve().parents[1]
S1 = "20261030T085136Z"
LOGROOT = "/tmp/pe_big/logs/g5"   # logs stay off the home quota
BASE = {"O1": "O1_r3", "O2": "O2_r1"}          # snapshot-1 runs with the same rules (D8, D11, verifier)
BENCH1 = {"O1": "O1", "O2": "O2_r1"}           # O1_r3 has the same SOGs as O1 (benchmark lives there)


def now():
    return subprocess.run(["date"], capture_output=True, text=True).stdout.strip()


def jload(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def f1(path, split_col=None, split=None):
    try:
        R = pd.read_csv(path)
    except Exception:
        return None
    m = (R.subset == "primary") & (R.tau.round(3) == 0.2)
    m &= (R[split_col] == split) if split_col else (R.half == "test")
    r = R[m]
    return None if r.empty else (round(float(r.f1.iloc[0]), 3), round(float(r.f1_ci_lo.iloc[0]), 3),
                                 round(float(r.f1_ci_hi.iloc[0]), 3))


def candidates(snap, release):
    """Post-rescue candidate table (O1: write_candidates.py; O2: write_o2_report.py), with the artifact-cascade
    columns once run_artifacts.py ran. No fallback to the pre-rescue run_events.py table: comparing a pre-rescue
    table of one snapshot with the post-rescue table of the other is what made the first G5 status overstate the
    O1 differences (Jaccard 0.5 / 0.25)."""
    f = P / f"results/partB/{snap}/{release}/candidates.parquet"
    return pd.read_parquet(f) if f.exists() else None


def primary(c):
    """Primary calls: cost 2:1, tau 0.2 (DESIGN §5; the artifact cascade summary counts these)."""
    if c is None:
        return None
    m = (c["tau"].round(3) == 0.2) if "tau" in c else True
    if "primary" in c:
        m &= c["primary"].astype(bool)
    return c[m]


def partA_numbers(snap, rid):
    if snap == S1:
        n = (jload(P / "results/numbers.json") or {}).get("partA") or {}
        if n.get("snapshot_id") == S1:
            return n, "results/numbers.json[partA]"
        return jload(P / f"results/partA/{S1}/v2_1/numbers_partA.json") or {}, f"results/partA/{S1}/v2_1"
    return jload(P / f"results/partA/{snap}/{rid}/numbers_partA.json") or {}, f"results/partA/{snap}/{rid}"


def fmt_ci(v, ci=None, nd=4):
    if v is None:
        return "n/a"
    s_ = f"{v:.{nd}f}"
    if ci is not None and len(ci) == 2 and None not in ci:
        s_ += f" [{ci[0]:.{nd}f}, {ci[1]:.{nd}f}]"
    return s_


def partA_rows(n):
    """name -> formatted value for the Part A headline rows."""
    out = {}
    for key, lab in (("F_per_pair", "per pair"), ("F_headline_per_family_pair", "per family pair")):
        x = (n.get(key) or {}).get("F_analog") or {}
        out[f"F_analog {lab}"] = fmt_ci(x.get("value"), x.get("ci95"))
    fu = [r for r in (n.get("functional") or []) if r.get("subset") == "all" and r.get("category") == "any"
          and r.get("label_variant") == "full" and not r.get("exploratory") and r.get("exclude_train", True)]
    for lab, k in (("per family pair", "F_conv_func.per_family_pair"), ("per pair", "F_conv_func.per_pair")):
        out[f"F_conv,func {lab}"] = (fmt_ci(fu[0].get(k), [fu[0].get(k + ".ci_lo"), fu[0].get(k + ".ci_hi")])
                                     if fu else "n/a")
    d = [r for r in (n.get("descent_2d_summary") or []) if r.get("dataset") == "O1_tau0.2" and r.get("half") == "all"]
    for lab, k in (("M1 ortholog", "M1_ortholog"), ("M1 paralog", "M1_paralog"),
                   ("M1 ortholog - paralog (identity x time)", "diff_identity_x_time")):
        out[lab] = fmt_ci(d[0].get(k), [d[0].get(k + ".ci_lo"), d[0].get(k + ".ci_hi")]) if d else "n/a"
    return out


def partA_descent(snap, rid, n):
    """descent_2d_summary rows: numbers.json for snapshot 1, the run's csv otherwise (11_numbers runs before 16)."""
    n = dict(n)
    base = P / (f"results/partA/{S1}/v2_1" if snap == S1 else f"results/partA/{snap}/{rid}")
    for key, f in (("descent_2d_summary", "descent_2d_summary.csv"), ("functional", "functional_convergence.csv")):
        if not n.get(key) and (base / f).exists():
            n[key] = pd.read_csv(base / f).to_dict("records")
    return n


def tiers_identical(s2):
    """Tiers 1-2 manifest rows (unit, sizes, mtimes, row counts, method config) identical between the snapshots?"""
    try:
        a = pd.read_parquet(P / f"data/snapshot/{S1}/manifest.parquet")
        b = pd.read_parquet(P / f"data/snapshot/{s2}/manifest.parquet")
    except Exception:
        return None
    cols = [c for c in a.columns if c in b.columns]
    res = {}
    for t in (1, 2):
        x = a[a.tier == t][cols].sort_values("unit_id").reset_index(drop=True)
        y = b[b.tier == t][cols].sort_values("unit_id").reset_index(drop=True)
        res[t] = (len(x), bool(x.equals(y)))
    return res


def overlap(c1, c2):
    if c1 is None or c2 is None:
        return None
    key = [k for k in ("char", "sog_id") if k in c1.columns and k in c2.columns][:1]
    key += [k for k in ("type", "node", "clade") if k in c1.columns and k in c2.columns]
    if not key:
        return None
    a = set(map(tuple, c1[key].astype(str).to_numpy()))
    b = set(map(tuple, c2[key].astype(str).to_numpy()))
    return {"key": key, "n1": len(a), "n2": len(b), "both": len(a & b), "jaccard": round(len(a & b) / max(len(a | b), 1), 3)}


def status_table(s2, rid):
    rows = []
    for f in sorted(glob.glob(f"{LOGROOT}/{s2}/*.jsonl")):
        for line in open(f):
            r = json.loads(line)
            st = "BLOCKED" if r.get("blocked") else ("ok" if r.get("rc") == 0 else f"FAILED rc={r.get('rc')}")
            rows.append(f"| {Path(f).stem} | {r['stage']} | {st} | {r.get('wall_s', 0):.0f} | "
                        f"{(r.get('blocked') or (r.get('stderr_tail') or '')[-160:]).replace('|', '/').replace(chr(10), ' ')} |")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s2", required=True)
    ap.add_argument("--run-id", default="g5")
    ap.add_argument("--status-only", action="store_true")
    ap.add_argument("--note", default="")
    ap.add_argument("--notes", default=None, help="markdown file appended as 'Reproducibility notes'")
    ap.add_argument("--out", default=str(P / "results/pipeline/G5_STATUS.md"))
    a = ap.parse_args()
    s2, rid = a.s2, a.run_id
    out = Path(a.out)
    L = [f"# G5 status: snapshot {s2} vs {S1}", "", f"Updated {now()} by `scripts/g5_compare.py` (launcher "
         f"`scripts/g5_launch.sh`, logs `{LOGROOT}/{s2}/`). Development subset, preliminary.", ""]
    if a.note:
        L += [a.note, ""]
    m2 = jload(P / f"data/snapshot/{s2}/manifest.json") or jload(f"/tmp/pe_big/snapshot/{s2}/manifest.json")
    m1 = jload(P / f"data/snapshot/{S1}/manifest.json")
    if m2 and m1:
        t1, t2 = m1["per_tier"]["3"], m2["per_tier"]["3"]
        L += [f"Tier 3: {t1['units']:,} -> {t2['units']:,} units; {t1['pockets']:,} -> {t2['pockets']:,} pockets.", ""]
    L += ["| lane | stage | status | wall s | note |", "|---|---|---|---|---|"] + status_table(s2, rid) + [""]
    if not a.status_only:
        L += ["## Snapshot 2 vs snapshot 1", "", "| quantity | snapshot 1 | snapshot 2 |", "|---|---|---|"]
        for o in ("O1", "O2"):
            f1_ = jload(P / f"results/sogs/{S1}/{BASE[o]}/funnel.json") or {}
            f2_ = jload(P / f"results/sogs/{s2}/{o}_{rid}/funnel.json") or {}
            for k in ("sogs", "sogs_multispecies", "members_with_model", "sites_in_run_ogs"):
                L.append(f"| {o} {k} | {f1_.get(k, 'n/a')} | {f2_.get(k, 'n/a')} |")
            for k in ("P", "A", "U", "G"):
                L.append(f"| {o} state {k} | {f1_.get('states', {}).get(k, 'n/a')} | {f2_.get('states', {}).get(k, 'n/a')} |")
            L.append(f"| {o} §8.1 test F1 (τ 0.2) [CI] | {f1(P / f'results/sogs/{S1}/{BENCH1[o]}/tau_benchmark.csv')} | "
                     f"{f1(P / f'results/sogs/{s2}/{o}_{rid}/tau_benchmark.csv')} |")
        L.append(f"| O2 D20 clean-test F1 (τ 0.2) | "
                 f"{f1(P / f'results/sogs/{S1}/O2_r1_split/tau_benchmark_o2_split.csv', 'split_o2', 'test')} | "
                 f"{f1(P / f'results/sogs/{s2}/O2_{rid}_split/tau_benchmark_o2_split.csv', 'split_o2', 'test')} |")
        r1 = jload(P / f"results/sogs/{S1}/O2_r1_rescue/rescue_summary.json") or {}
        r2 = jload(P / f"results/sogs/{s2}/O2_{rid}_rescue/rescue_summary.json") or {}
        L.append(f"| O2 rescue pool rows -> W / A | {r1.get('n_selected', 'n/a')} -> {r1.get('result', 'n/a')} | "
                 f"{r2.get('n_selected', 'n/a')} -> {r2.get('result', 'n/a')} |")
        # Part A headline (O1 universe = tiers 1-2: identical inputs -> reproducibility check)
        n1, src1 = partA_numbers(S1, rid)
        n2, src2 = partA_numbers(s2, rid)
        a1, a2 = partA_rows(partA_descent(S1, rid, n1)), partA_rows(partA_descent(s2, rid, n2))
        for k in a1:
            L.append(f"| Part A {k} [95% CI] | {a1[k]} | {a2.get(k, 'n/a')} |")
        # M1 like-for-like: the reported snapshot-1 M1 used the O1 (pre-D11) states; the G5 O1 run uses the O1_r3
        # rules (D8, D11). Snapshot-1 reference on O1_r3 states: results/pipeline/G5/<s2>/s1_ref_O1_r3__descent_2d_summary.csv
        ref = P / f"results/pipeline/G5/{s2}/s1_ref_O1_r3__descent_2d_summary.csv"
        if ref.exists():
            r_ = pd.read_csv(ref)
            r_ = r_[(r_.dataset == "O1_r3_tau0.2") & (r_.half == "all")].to_dict("records")
            ar = partA_rows({"descent_2d_summary": [dict(x, dataset="O1_tau0.2") for x in r_]})
            for k in ("M1 ortholog", "M1 paralog", "M1 ortholog - paralog (identity x time)"):
                L.append(f"| Part A {k}, like-for-like (snapshot 1 on O1_r3 states) | {ar[k]} | {a2.get(k, 'n/a')} |")
        # Part B: primary (2:1, tau 0.2) calls after rescue, then after the artifact cascade; overlap by event key
        for o in ("O1", "O2"):
            c1, c2 = primary(candidates(S1, BASE[o])), primary(candidates(s2, f"{o}_{rid}"))
            for t in ("gain", "loss"):
                for stage, sel in (("after rescue", None), ("after artifact cascade", "pass_artifacts")):
                    d1 = None if c1 is None else c1[c1.type == t]
                    d2 = None if c2 is None else c2[c2.type == t]
                    if sel:
                        d1 = d1[d1[sel].astype(bool)] if d1 is not None and sel in d1 else None
                        d2 = d2[d2[sel].astype(bool)] if d2 is not None and sel in d2 else None
                    ov = overlap(d1, d2)
                    v1 = "n/a" if d1 is None else str(len(d1))
                    v2 = "n/a" if d2 is None else str(len(d2))
                    ovs = "" if ov is None else f" (in both: {ov['both']}, Jaccard {ov['jaccard']})"
                    L.append(f"| {o} {'gains' if t == 'gain' else 'losses'} {stage} (2:1, τ 0.2) | {v1} | {v2}{ovs} |")
            c1a, c2a = candidates(S1, BASE[o]), candidates(s2, f"{o}_{rid}")
            if c1a is not None and c2a is not None:
                ov = overlap(c1a, c2a)
                L.append(f"| {o} all FDR/audit candidate rows, all costs and τ (unique events) | {len(c1a)} ({ov['n1']}) | "
                         f"{len(c2a)} ({ov['n2']}; in both {ov['both']}, Jaccard {ov['jaccard']}) |")
        L += ["", f"Overlap key: {'+'.join(overlap(*(primary(candidates(S1, BASE['O1'])),) * 2)['key']) if candidates(S1, BASE['O1']) is not None else 'sog_id+type+node'}"
              f" (same event = same SOG, same type, same tree node). Snapshot-1 comparators: {BASE['O1']}, {BASE['O2']} "
              f"(same rules: D8, D11, verifier); Part A {src1} vs {src2}.", ""]
        ti = tiers_identical(s2)
        L += ["## Reproducibility of the O1 side", ""]
        if ti:
            L.append("Tiers 1-2 manifest rows (units, file sizes, mtimes, row counts, method config) between the snapshots: "
                     + ", ".join(f"tier {t}: {n_} units, {'identical' if same else 'DIFFERENT'}" for t, (n_, same) in ti.items())
                     + ". O1 SOGs, Part A (tiers 1-2 universe) and Part B O1 therefore have identical inputs: every "
                     "difference in their rows above measures reproducibility (code version, seeds, nondeterminism), "
                     "not new data. Only tier 3 (O2 members, web index) sees new data.")
            L.append("")
        if a.notes and Path(a.notes).exists():
            L += [Path(a.notes).read_text().rstrip(), ""]
        L += ["", "n/a = the stage did not run or was blocked (see the table above)."]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L) + "\n")
    print(out)


if __name__ == "__main__":
    main()
