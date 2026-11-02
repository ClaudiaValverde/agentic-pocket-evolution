#!/usr/bin/env python3
"""Ligand-free refold + BURROW for the pool cases that cannot be ligand co-folded.

Project decision (2026-10-30): Tier 2 (metal-only ligand) and Tier 3 (no cognate ligand) get NO
ligand co-folding -- metals are uninterpretable by co-folding (panel AUC 0.58) and
guessing a ligand would produce an uninformative null. Instead every such case, and every
false-loss control, gets the preregistered ligand_free_refold (protein only, 3 seeds x 5
samples) followed by `pocketcompare describe --af-mode`, asking an independent question:

    does the pocket appear in the present-side members and not in the absent-side ones?

Scored with site recall as a CONTINUOUS value (not thresholded at 0.8): the calibration
controls showed true positives reach only 0.44-1.00 because a 4.5 A ligand footprint and a
BURROW cavity are different objects.

False-loss controls: the absent-side members are the pre-D11 state-A members recovered
from results/sogs/<snap>/O1_r2/web/sog_states/, i.e. the models whose A call the quality
filters later overturned. The expected result is that the pocket REAPPEARS in the fresh
model, confirming the loss was a prediction artifact.

Constructs: proteins over --max-aa are folded as the TED domain that contains the
footprint, with footprint positions renumbered into the construct; recorded as `construct`.

Usage: run_apo_pool.py --build --run --score --gpus 0,1,2
"""
from __future__ import annotations
import argparse, glob, json, os, re, subprocess, sys, time
from pathlib import Path
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import colabfold_msa as CF

P = Path(__file__).resolve().parents[2]  # repo root
SNAP = "20261030T085136Z"
AF3_PY = os.environ.get("AF3_ENV", "/tmp/af3_env") + "/bin/python"
RUN = P / "scripts/cofold/run_af3.sh"
PC = P / ".venv/bin/pocketcompare"
BIG = Path("/tmp/pe_big/cofold/apo_pool")
RES = P / "results/cofold/apo_pool"

FALSE_LOSS = {
    "control_falseloss_HOG0027591_g2_n66_O1": "HOG0027591:g2",
    "control_falseloss_HOG0028108_g15_n235_O1": "HOG0028108:g15",
    "control_falseloss_HOG0030328_g4_n103_O1": "HOG0030328:g4",
    "control_falseloss_HOG0034353_g3_n4659_O1": "HOG0034353:g3",
    "control_falseloss_HOG0035718_g5_n12_O1": "HOG0035718:g5",
}
TIER23 = {
    "gain_HOG0056135_g10_n70_O2": "tier2_metal",      # membrane dipeptidase, Zn
    "gain_HOG0034532_g4_n8_O2": "tier2_metal",        # connexin, Ca
    "dupgain_HOG0046679_g11_n33_O2": "tier3_no_ligand",
    "gain_HOG0025964_g2_n28_O1": "tier3_no_ligand",   # IQGAP, ~1.6 kaa -> TED construct
    "gain_HOG0026737_g5_n16_O1": "tier3_no_ligand",   # VAC14
}


def log(m): print(f"[apo_pool] {m}", file=sys.stderr, flush=True)


def load_states(release: str):
    import pandas as pd
    fs = glob.glob(str(P / f"results/sogs/{SNAP}/{release}/web/sog_states/bucket=*.parquet"))
    return pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)


def ted_construct(acc: str, site: list[int], seq: str):
    """Return (sub_sequence, renumbered_site, label) for the TED domain holding the site."""
    import pandas as pd
    td = pd.read_parquet(P / "data/annotations/ted_domains.parquet",
                         columns=["accession", "chopping", "nres_domain", "plddt",
                                  "af_fragment", "seq_matches_uniprot"])
    td = td[(td.accession == acc) & (td.af_fragment == "F1") & td.seq_matches_uniprot]
    best, best_cov = None, -1
    for _, r in td.iterrows():
        segs = []
        for part in str(r.chopping).split("_"):
            if "-" in part:
                a, b = part.split("-")
                segs.append((int(a), int(b)))
        idx = [i for a, b in segs for i in range(a, b + 1)]
        cov = len(set(idx) & set(site))
        if cov > best_cov:
            best, best_cov, bidx = r, cov, idx
    if best is None or best_cov == 0:
        return None
    keep = sorted(set(bidx))
    sub = "".join(seq[i - 1] for i in keep if 0 < i <= len(seq))
    remap = {p: j + 1 for j, p in enumerate(keep)}
    return sub, sorted(remap[p] for p in site if p in remap), \
        f"TED:{best.chopping}({len(keep)}aa,cov{best_cov}/{len(site)})"


def members_for(cid: str, max_aa: int) -> list[dict]:
    """2-3 present-side and 2-3 absent-side members, with footprint and construct."""
    out = []
    doc = yaml.safe_load((P / "cases" / cid / "case.yaml").read_text())
    by_acc = {p["accession"]: p for p in doc.get("proteins") or []}
    if cid in FALSE_LOSS:
        sog = FALSE_LOSS[cid]
        r2, r3 = load_states("O1_r2"), load_states("O1_r3")
        pres = r3[(r3.sog_id == sog) & (r3.state == "P")].sort_values(
            "footprint_plddt", ascending=False).head(3)
        absn = r2[(r2.sog_id == sog) & (r2.state == "A")].sort_values(
            "footprint_plddt").head(3)
        picks = [(r.accession, "present", float(r.footprint_plddt or 0)) for _, r in pres.iterrows()]
        picks += [(r.accession, "absent_pre_d11", float(r.footprint_plddt or 0)) for _, r in absn.iterrows()]
    else:
        sides = {"P": [], "A": []}
        for p in doc.get("proteins") or []:
            if p.get("state") in sides:
                sides[p["state"]].append(p)
        picks = []
        for st, lab in (("P", "present"), ("A", "absent")):
            for p in sorted(sides[st], key=lambda x: -(x.get("footprint_plddt") or 0))[:3]:
                picks.append((p["accession"], lab, float(p.get("footprint_plddt") or 0)))
    # reference footprint: from any present-side member that lists residues
    ref = None
    for p in doc.get("proteins") or []:
        if p.get("state") == "P" and p.get("residues"):
            ref = (p["accession"], sorted({int(r.split(":")[1]) for r in p["residues"].split()
                                           if ":" in r}))
            break
    for acc, side, plddt in picks:
        try:
            seq = CF.fetch_uniprot(acc)
        except Exception as e:  # noqa: BLE001
            log(f"{cid} {acc}: UniProt fetch failed {e}")
            continue
        own = by_acc.get(acc, {}).get("residues")
        site = sorted({int(r.split(":")[1]) for r in own.split() if ":" in r}) if own else []
        site_src = "own_footprint"
        if not site and ref:
            import case_cofold as CC
            site = CC.map_positions(CF.fetch_uniprot(ref[0]), seq, ref[1])
            site_src = f"mapped_from:{ref[0]}"
        construct = "full_length"
        if len(seq) > max_aa and site:
            tc = ted_construct(acc, site, seq)
            if tc:
                seq, site, construct = tc
        out.append({"case": cid, "acc": acc, "side": side, "atlas_plddt": plddt,
                    "seq": seq, "length": len(seq), "site": site,
                    "site_source": site_src, "construct": construct})
    return out


def build(max_aa: int) -> list[dict]:
    BIG.mkdir(parents=True, exist_ok=True)
    allm = []
    for cid in list(FALSE_LOSS) + list(TIER23):
        ms = members_for(cid, max_aa)
        log(f"{cid}: {len(ms)} members "
            f"({sum(m['side'].startswith('present') for m in ms)} present / "
            f"{sum(m['side'].startswith('absent') for m in ms)} absent)")
        allm += ms
    for m in allm:
        prov = CF.fetch_msa(m["acc"] if m["construct"] == "full_length"
                            else f"{m['acc']}_{m['construct']}",
                            m["seq"], host=CF.HOST, mode="env", cache=CF.CACHE,
                            min_interval=20.0)
        name = f"{m['acc'].lower()}_apo"
        d = BIG / "inputs" / m["case"]
        d.mkdir(parents=True, exist_ok=True)
        out = d / f"{name}.json"
        CF.write_af3_json(str(out), name, m["seq"], CF.load_a3m(m["seq"], CF.CACHE),
                          [], [1, 2, 3], prov)
        meta = out.with_suffix(".meta.json")
        if meta.exists():
            (d / "meta").mkdir(exist_ok=True)
            meta.rename(d / "meta" / meta.name)
        m["input"], m["name"] = str(out), name
        m["msa_n_seqs"] = prov["n_seqs"]
    (BIG / "members.json").write_text(json.dumps(allm, indent=1))
    log(f"built {len(allm)} ligand-free jobs")
    return allm


def run(gpus: list[int]) -> None:
    ms = json.loads((BIG / "members.json").read_text())
    cases = sorted({m["case"] for m in ms})
    load = {c: sum(m["length"] for m in ms if m["case"] == c) for c in cases}
    q = {g: [] for g in gpus}
    for c in sorted(cases, key=lambda c: -load[c]):
        g = min(q, key=lambda g: sum(load[x] for x in q[g]))
        q[g].append(c)
    procs = []
    for g, cs in q.items():
        if not cs:
            continue
        script = " ; ".join(f'bash {RUN} -i {BIG/"inputs"/c} -j apo_{c} -g {g}' for c in cs)
        lf = open(BIG / f"gpu{g}.out", "w")
        procs.append((g, subprocess.Popen(["bash", "-c", script], stdout=lf,
                                          stderr=subprocess.STDOUT,
                                          env={**os.environ, "AF3_OUT_ROOT": str(BIG / "af3")}),
                      lf, time.time(), cs))
        log(f"GPU {g}: {len(cs)} cases ({sum(load[x] for x in cs)} aa)")
    for g, p, f, t0, cs in procs:
        rc = p.wait(); f.close()
        log(f"GPU {g} done rc={rc} in {time.time()-t0:.0f}s")


def score() -> None:
    import pandas as pd
    import zstandard
    RES.mkdir(parents=True, exist_ok=True)
    ms = json.loads((BIG / "members.json").read_text())
    work = BIG / "describe"; work.mkdir(exist_ok=True)
    rows = []
    for m in ms:
        jd = BIG / "af3" / f"apo_{m['case']}" / m["name"]
        top = [f for f in glob.glob(f"{jd}/*_model.cif*") if "seed-" not in Path(f).parent.name]
        if not top:
            log(f"missing model for {m['acc']}")
            continue
        raw = open(top[0], "rb").read()
        if top[0].endswith(".zst"):
            raw = zstandard.ZstdDecompressor().decompress(raw, max_output_size=2**31)
        cif = work / f"{m['acc']}.cif"; cif.write_bytes(raw)
        pq = work / f"{m['acc']}.parquet"
        r = subprocess.run([str(PC), "describe", str(cif), "-o", str(pq), "--af-mode"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            log(f"{m['acc']}: describe failed {r.stderr[-300:]}")
            continue
        df = pd.read_parquet(pq)
        site = {f"A:{p}" for p in m["site"]}
        best = None
        for _, p in df.iterrows():
            res = set(str(p["residues"]).split())
            if not res or not site:
                continue
            inter = len(res & site)
            cand = {"recall": inter / len(site), "precision": inter / len(res),
                    "jaccard": inter / len(res | site),
                    "burrow_rank": int(p["burrow_rank"]),
                    "burrow_score": round(float(p["burrow_score"]), 4),
                    "pocket_plddt": round(float(p.get("plddt_mean", float("nan"))), 1)}
            if best is None or cand["recall"] > best["recall"]:
                best = cand
        rows.append({k: m[k] for k in ("case", "acc", "side", "length", "construct",
                                       "site_source", "atlas_plddt", "msa_n_seqs")}
                    | {"n_site": len(site), "n_pockets": len(df)}
                    | ({k: round(v, 3) if isinstance(v, float) else v
                        for k, v in best.items()} if best else {"recall": 0.0}))
        log(f"{m['case'][:28]:28s} {m['acc']:12s} {m['side']:15s} "
            f"recall={rows[-1].get('recall')} rank={rows[-1].get('burrow_rank')}")
    df = pd.DataFrame(rows)
    df.to_csv(RES / "apo_burrow_pool.tsv", sep="\t", index=False)
    summ = []
    for c, g in df.groupby("case"):
        pres = g[g.side.str.startswith("present")].recall
        absn = g[g.side.str.startswith("absent")].recall
        summ.append({"case": c, "n_present": len(pres), "n_absent": len(absn),
                     "present_recall_median": round(float(pres.median()), 3) if len(pres) else None,
                     "absent_recall_median": round(float(absn.median()), 3) if len(absn) else None,
                     "delta_recall": round(float(pres.median() - absn.median()), 3)
                     if len(pres) and len(absn) else None})
    pd.DataFrame(summ).to_csv(RES / "apo_pool_summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(summ).to_string(index=False))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--gpus", default="0,1,2")
    ap.add_argument("--max-aa", type=int, default=1000)
    a = ap.parse_args()
    if a.build:
        build(a.max_aa)
    if a.run:
        run([int(g) for g in a.gpus.split(",")])
    if a.score:
        score()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
