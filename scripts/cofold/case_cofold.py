#!/usr/bin/env python3
"""Run and score a co-folding case, exactly as frozen in config/prereg.yaml (cofold).

A case is a YAML spec:

    case: cdk2_toy
    positive_control: {acc: P24941, pdb: 1hck, chain: A, ligand: ATP, seqnum: 400}
    site: [10, 11, 12, ...]              # 1-based positions of the known site (optional)
    ligands: [CCD:ATP, CCD:MG]           # cognate ligand(s), shared by all members
    decoys:  [[CCD:CFF]]                 # each entry is one decoy ligand set
    seeds:   [1, 2, 3]
    members:
      - {accession: P24941, state: derived,   label: human_CDK2}
      - {accession: P24941, state: ancestral, label: CDK2_site_ablated,
         sequence_fasta: /tmp/pe_big/ablation/cdk2_gly.fa}

What it does
  1. MSA source. `--msa auto` (default) uses the cached AF3 jackhmmer MSAs only when
     EVERY member matches a cached sequence exactly; otherwise it fetches ColabFold
     MMseqs2 MSAs for ALL members, so provenance is identical within the case
     (prereg cofold.msa, review change 11). `--msa colabfold|cached` forces it. The
     chosen source and the per-member depth are recorded in msa_source.tsv.
     `reuse_msa_for_variant` is never used across the derived/ancestral boundary.
  2. Templates. `--no-templates` always (prereg cofold.templates.primary).
  3. Blinding. Every job is named `<case>_<12-hex>` where the hex is HMAC(salt, member|
     ligandset). The readout table is written keyed by hash only; the key lives in a
     separate blind_key.json. `--unblind` joins them and applies the decision rule.
  4. Jobs. members x (cognate + each decoy), 3 seeds x 5 samples, spread over --gpus
     as balanced directory queues (one AF3 process per GPU, model built once).
  5. Decision rule (prereg cofold.decision_rule, statistic = per-seed means, must hold
     in ALL seeds): with d = mean over derived members, a = mean over ancestral members,
         D_plddt(seed) = d.lig_plddt - a.lig_plddt
         D_pae(seed)   = a.pae_min   - d.pae_min
         X_plddt(seed) = mean over members of (cognate.lig_plddt - decoy.lig_plddt)
         X_pae(seed)   = mean over members of (decoy.pae_min   - cognate.pae_min)
     the case is called `gain_in_derived` iff D_plddt > X_plddt and D_pae > X_pae in
     every seed; `loss_in_derived` iff the same holds with the signs reversed; else
     `inconclusive`.
  6. Positive-control gate (prereg cofold.positive_control_gate): the holo complex named
     in the spec is folded and its top-ranked ligand RMSD to the crystal pose measured.
     RMSD > 2.0 A -> the whole case is labelled `uninterpretable`, never dropped.

Usage
  case_cofold.py --spec cases/cofold/cdk2_toy.yaml --build --run --unblind
  case_cofold.py --spec cases/<case_id>/case.yaml --build --run --unblind   # schema case-v1
  case_cofold.py --spec ... --unblind            # re-score without re-running

A `schema: case-v1` file from scripts/cases/build_cases.py is accepted directly: members
are taken by pocket state (P -> derived, A -> ancestral, at most 2 per side, highest
footprint pLDDT first), and the site is the union of the present-side footprints. The
case author must add `cofold_ligands:` (and optionally `cofold_decoys:`).
"""
from __future__ import annotations

import argparse
import hashlib
import hmac
from concurrent.futures import ThreadPoolExecutor
import json
import os
import statistics as stats
import subprocess
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cofold_readouts import summarise_job, eval_sample  # noqa: E402
import colabfold_msa as CF  # noqa: E402

P = Path(__file__).resolve().parents[2]  # repo root
AF3_PY = os.environ.get("AF3_ENV", "/tmp/af3_env") + "/bin/python"
MK = P / "scripts/cofold/make_af3_input.py"
RUN = P / "scripts/cofold/run_af3.sh"
BIG = Path("/tmp/pe_big/cofold/cases")
PREREG = P / "config/prereg.yaml"


def log(m: str) -> None:
    print(f"[case_cofold] {m}", file=sys.stderr, flush=True)


def load_prereg() -> dict:
    c = yaml.safe_load(PREREG.read_text())["cofold"]
    return {
        "seeds": 3,
        "plddt_true_min": c["thresholds_from_benchmark"]["ligand_plddt_true_min"]["value"],
        "pae_true_max": c["thresholds_from_benchmark"]["min_pl_pae_true_max"]["value"],
        "centroid_true_max": c["thresholds_from_benchmark"]["centroid_distance_true_max_A"]["value"],
        "gate_rmsd_max": c["positive_control_gate"]["ligand_rmsd_max_A"]["value"],
        "failed_gate_label": c["positive_control_gate"]["failed_gate_label"]["value"],
        "blinding": c["blinding"]["value"],
        "primary": c["readouts"]["primary"]["value"],
    }


# ----------------------------------------------------------------- spec helpers
def member_sequence(m: dict) -> str:
    if m.get("sequence"):
        return m["sequence"].strip().upper()
    if m.get("sequence_fasta"):
        return CF.read_fasta_multi(m["sequence_fasta"])[0][1]
    return CF.fetch_uniprot(m["accession"])


def cached_exact(acc: str, seq: str) -> bool:
    sys.path.insert(0, str(P / "scripts/cofold"))
    import make_af3_input as MKI
    c = MKI.find_cache(acc, MKI.DEFAULT_CACHE_DIRS)
    if not c:
        return False
    return MKI.load_cached_protein(c)["sequence"] == seq


def map_positions(ref_seq: str, tgt_seq: str, positions: list[int]) -> list[int]:
    """Map 1-based positions from ref_seq onto tgt_seq by global pairwise alignment."""
    import gemmi
    aln = gemmi.align_string_sequences(list(ref_seq), list(tgt_seq), [])
    m, i, j = {}, 0, 0
    for a, b in zip(aln.add_gaps(ref_seq, 1), aln.add_gaps(tgt_seq, 2)):
        if a != "-" and b != "-":
            m[i + 1] = j + 1
        i += a != "-"
        j += b != "-"
    return sorted(m[p] for p in positions if p in m)


def adapt_case_v1(doc: dict, path: Path, max_per_side: int = 3) -> dict:
    """Translate a `schema: case-v1` file (scripts/cases/build_cases.py) into a cofold spec.

    Members are chosen by POCKET STATE, not by the tree role, so the same rule works for
    a gain and for a loss: state P (pocket present, footprint residues listed) becomes
    `derived`, state A (pocket absent) becomes `ancestral`. The tree role is kept in the
    label so the unblinded table still shows it. W (rescued) and U (unobserved) are never
    folded -- a case whose sides are W/U has nothing for co-folding to compare.

    The site is the footprint of the present-side members (`residues`, "A:123" form).
    Ligands must be supplied by the case author as `cofold_ligands` / `cofold_decoys`
    (or `ligands` / `decoys`); there is no way to guess them, and guessing would make the
    readouts meaningless.
    """
    prots = doc.get("proteins") or []
    sides = {"P": [], "A": []}
    for p in prots:
        st = p.get("state")
        if st in sides:
            sides[st].append(p)
    # prefer the highest-confidence models on each side, then alphabetical for determinism
    def rank(p):
        return (-(p.get("footprint_plddt") or 0), str(p.get("accession")))

    def own_site(p) -> list[int]:
        return sorted({int(r.split(":")[1]) for r in str(p.get("residues") or "").split()
                       if ":" in r})

    members = []
    for st, state_label in (("P", "derived"), ("A", "ancestral")):
        for p in sorted(sides[st], key=rank)[:max_per_side]:
            members.append({"accession": p["accession"], "state": state_label,
                            "label": f"{p['accession']}_{p.get('role', 'na')}_{st}",
                            "site": own_site(p)})
    # Residue numbers are per protein, so each member carries its OWN site. Absent-side
    # members have no footprint of their own: map the best present-side footprint onto
    # them by pairwise alignment of the two sequences (the positional-homology call the
    # OG MSA makes at scale; recorded as site_source).
    ref = next((m for m in members if m["state"] == "derived" and m["site"]), None)
    if ref is None:
        raise SystemExit(f"{path}: no present-side member carries footprint residues; "
                         "cannot define a site for the readouts.")
    ref_seq = CF.fetch_uniprot(ref["accession"])
    for m in members:
        if m["site"]:
            m["site_source"] = "own_footprint"
            continue
        m["site"] = map_positions(ref_seq, CF.fetch_uniprot(m["accession"]), ref["site"])
        m["site_source"] = f"mapped_from:{ref['accession']}"
        log(f"  {m['label']}: site mapped from {ref['accession']} "
            f"({len(m['site'])}/{len(ref['site'])} positions aligned)")
    ligands = doc.get("cofold_ligands") or doc.get("ligands")
    if not ligands:
        raise SystemExit(
            f"{path}: no cofold_ligands/ligands in the case file. Co-folding needs the "
            "cognate ligand(s) as CCD:<code> or SMILES:<smiles>; add cofold_ligands "
            "(and cofold_decoys, or let --auto-decoys pick property-matched ones).")
    missing = [s for s in ("P", "A") if len(sides[s]) == 0]
    if missing:
        raise SystemExit(f"{path}: no members in state(s) {missing}; nothing to compare. "
                         "States present: "
                         f"{ {p.get('state') for p in prots} }")
    spec = {"case": doc.get("case_id") or path.parent.name,
            "source_schema": "case-v1", "source_file": str(path),
            "members": members, "ligands": ligands,
            "decoys": doc.get("cofold_decoys") or doc.get("decoys") or [],
            "seeds": doc.get("seeds", [1, 2, 3])}
    if doc.get("positive_control"):
        spec["positive_control"] = doc["positive_control"]
    log(f"case-v1 {spec['case']}: {len(members)} members "
        f"({len(sides['P'])} P / {len(sides['A'])} A available); "
        f"site sizes {[len(m['site']) for m in members]}")
    return spec


def prefetch_msas(case_files: list[str], max_per_side: int = 2,
                  workers: int = 4) -> None:
    """Fetch ColabFold MSAs for the members any of these cases would fold.

    Runs the same member-selection rule as adapt_case_v1 but needs no ligands, so it can
    be run before the case author has chosen them -- the MSA depends only on the sequence.
    Everything is cached by sequence hash, so the later --build is then free.
    """
    want: dict[str, str] = {}
    for f in case_files:
        doc = yaml.safe_load(Path(f).read_text())
        if not str(doc.get("schema", "")).startswith("case-v1"):
            continue
        sides: dict[str, list] = {"P": [], "A": []}
        for pr in doc.get("proteins") or []:
            if pr.get("state") in sides:
                sides[pr["state"]].append(pr)
        for st in ("P", "A"):
            for pr in sorted(sides[st], key=lambda x: (-(x.get("footprint_plddt") or 0),
                                                       str(x.get("accession"))))[:max_per_side]:
                want[pr["accession"]] = doc.get("case_id", Path(f).parent.name)
    log(f"prefetch: {len(want)} distinct accessions over {len(case_files)} case files")

    def one(item):
        acc, case = item
        try:
            seq = CF.fetch_uniprot(acc)
            CF.fetch_msa(acc, seq, host=CF.HOST, mode="env", cache=CF.CACHE,
                         min_interval=20.0)
            return acc, None
        except Exception as e:  # noqa: BLE001
            return acc, f"{case}: {e}"

    # Submissions stay spaced at one per 20 s by the shared rate limiter; the workers only
    # overlap the 1-5 min server-side wait, which is what makes 48 members feasible.
    ok, fails = 0, []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for acc, err in ex.map(one, sorted(want.items())):
            if err:
                fails.append(acc)
                log(f"prefetch {acc}: FAILED {err}")
            else:
                ok += 1
    log(f"prefetch done: {ok} cached, {len(fails)} failed {fails}")


def load_spec(path: str | Path) -> dict:
    p = Path(path)
    doc = yaml.safe_load(p.read_text())
    if str(doc.get("schema", "")).startswith("case-v1"):
        return adapt_case_v1(doc, p)
    return doc


def blind_id(salt: bytes, case: str, member_label: str, ligset: str) -> str:
    return hmac.new(salt, f"{case}|{member_label}|{ligset}".encode(), hashlib.sha256).hexdigest()[:12]


def ligset_tag(ligs: list[str]) -> str:
    # Accepts "CCD:ATP"/"SMILES:..." and bare CCD codes (recovered jobs store bare codes).
    return "+".join(l.split(":", 1)[1] if ":" in l else l for l in ligs)


# ----------------------------------------------------------------- build
def build(spec: dict, msa_mode: str, workdir: Path) -> dict:
    case = spec["case"]
    seeds = spec.get("seeds", [1, 2, 3])
    members = spec["members"]
    ligsets = [("cognate", spec["ligands"])] + [
        (f"decoy{i+1}", d) for i, d in enumerate(spec.get("decoys", []))]
    inputs = workdir / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    # The blinding salt must be STABLE for a workdir. A fresh salt on every build gives
    # every job a new blind id, so a later re-score writes ids that no longer match the
    # AF3 output directories already on disk, and silently finds zero readouts.
    keyfile = workdir / "blind_key.json"
    if keyfile.exists():
        salt = bytes.fromhex(json.loads(keyfile.read_text())["salt_hex"])
        log(f"{case}: reusing the existing blinding salt ({keyfile.name})")
    else:
        salt = os.urandom(16)

    seqs = {m["label"]: member_sequence(m) for m in members}
    all_cached = all(cached_exact(m["accession"], seqs[m["label"]]) for m in members)
    source = msa_mode if msa_mode != "auto" else ("cached" if all_cached else "colabfold")
    log(f"{case}: MSA source = {source} (all members cached-exact: {all_cached})")

    msa_rows, msa_prov = [], {}
    if source == "colabfold":
        for m in members:
            lab = m["label"]
            seq = seqs[lab]
            prov = CF.fetch_msa(lab, seq, host=CF.HOST, mode="env", cache=CF.CACHE,
                                min_interval=20.0)
            msa_prov[lab] = prov
            msa_rows.append({"member": lab, "acc": m["accession"], "state": m["state"],
                             "source": "colabfold_mmseqs2", "n_seqs": prov["n_seqs"],
                             "depth_over_sqrtL": prov["neff_proxy"], "ticket": prov["ticket"]})
    else:
        for m in members:
            lab = m["label"]
            msa_rows.append({"member": lab, "acc": m["accession"], "state": m["state"],
                             "source": "cached_af3_jackhmmer", "n_seqs": None,
                             "depth_over_sqrtL": None, "ticket": None})

    jobs, key = [], {}
    for m in members:
        lab = m["label"]
        seq = seqs[lab]
        for kind, ligs in ligsets:
            bid = blind_id(salt, case, lab, kind)
            name = f"{case}_{bid}"
            out = inputs / f"{name}.json"
            if source == "colabfold":
                CF.write_af3_json(str(out), name, seq, CF.load_a3m(seq, CF.CACHE),
                                  ligs, seeds, msa_prov[lab])
            else:
                cmd = [AF3_PY, str(MK), "--acc", m["accession"], "--seq", seq,
                       "--no-uniprot", "--no-templates", "--seeds",
                       ",".join(map(str, seeds)), "--name", name, "--out", str(out)]
                for l in ligs:
                    cmd += ["--ligand", l]
                r = subprocess.run(cmd, text=True, capture_output=True)
                if r.returncode != 0:
                    raise SystemExit(f"input build failed for {name}: {r.stderr}")
            jobs.append({"blind_id": bid, "name": name, "input": str(out),
                         "length": len(seq), "ligands": ligs,
                         "site": m.get("site") or spec.get("site") or []})
            key[bid] = {"member": lab, "acc": m["accession"], "state": m["state"],
                        "ligand_kind": kind, "ligands": ligset_tag(ligs),
                        "site_source": m.get("site_source", "spec"),
                        "n_site": len(m.get("site") or spec.get("site") or [])}

    # positive control: holo complex of the family, folded under the same settings
    pc = spec.get("positive_control")
    if pc:
        pseq = CF.fetch_uniprot(pc["acc"])
        name = f"{case}_poscontrol"
        out = inputs / f"{name}.json"
        pligs = pc.get("ligands") or spec["ligands"]
        if source == "colabfold":
            prov = CF.fetch_msa(f"PC_{pc['acc']}", pseq, host=CF.HOST, mode="env",
                                cache=CF.CACHE, min_interval=20.0)
            CF.write_af3_json(str(out), name, pseq, CF.load_a3m(pseq, CF.CACHE),
                              pligs, seeds, prov)
        else:
            cmd = [AF3_PY, str(MK), "--acc", pc["acc"], "--seq", pseq, "--no-uniprot",
                   "--no-templates", "--seeds", ",".join(map(str, seeds)),
                   "--name", name, "--out", str(out)]
            for l in pligs:
                cmd += ["--ligand", l]
            subprocess.run(cmd, check=True, text=True, capture_output=True)
        jobs.append({"blind_id": "poscontrol", "name": name, "input": str(out),
                     "length": len(pseq), "ligands": pligs,
                     "site": pc.get("site") or spec.get("site") or []})

    tmpj = workdir / "jobs.json.tmp"
    tmpj.write_text(json.dumps(jobs, indent=1))
    tmpj.replace(workdir / "jobs.json")
    tmpk = workdir / "blind_key.json.tmp"
    tmpk.write_text(json.dumps({"salt_hex": salt.hex(), "map": key}, indent=1))
    tmpk.replace(workdir / "blind_key.json")
    write_tsv(workdir / "msa_source.tsv", msa_rows)
    log(f"{case}: {len(jobs)} jobs built in {inputs}")
    return {"jobs": jobs, "source": source}


# ----------------------------------------------------------------- run
def run(workdir: Path, gpus: list[int]) -> None:
    jobs = json.loads((workdir / "jobs.json").read_text())
    queues = {g: [] for g in gpus}
    for j in sorted(jobs, key=lambda j: -j["length"]):           # longest-first balance
        g = min(queues, key=lambda g: sum(x["length"] for x in queues[g]))
        queues[g].append(j)
    procs = []
    for g, qjobs in queues.items():
        if not qjobs:
            continue
        qd = workdir / f"queue_gpu{g}"
        qd.mkdir(parents=True, exist_ok=True)
        for j in qjobs:
            (qd / Path(j["input"]).name).write_text(Path(j["input"]).read_text())
        logf = open(workdir / f"gpu{g}.out", "w")
        procs.append((g, subprocess.Popen(
            ["bash", str(RUN), "-i", str(qd), "-j", f"{workdir.name}_gpu{g}", "-g", str(g)],
            stdout=logf, stderr=subprocess.STDOUT,
            env={**os.environ, "AF3_OUT_ROOT": str(workdir / "af3")}), logf, time.time()))
        log(f"GPU {g}: {len(qjobs)} jobs ({sum(x['length'] for x in qjobs)} aa)")
    for g, p, f, t0 in procs:
        rc = p.wait()
        f.close()
        log(f"GPU {g} done rc={rc} in {time.time()-t0:.0f}s")


# ----------------------------------------------------------------- score
def per_seed(rows: list[dict], key: str) -> dict[int, float]:
    out = {}
    for sd in sorted({r["seed"] for r in rows}):
        v = [r[key] for r in rows if r["seed"] == sd and r.get(key) is not None]
        if v:
            out[sd] = stats.mean(v)
    return out


def collect(workdir: Path, spec: dict) -> tuple[list[dict], dict]:
    """Blinded readouts: one row per job, keyed by blind_id only."""
    jobs = json.loads((workdir / "jobs.json").read_text())
    # Which ligand do the readouts describe? pick_ligand() otherwise takes the LARGEST
    # non-ion ligand, which for a cofactor-plus-substrate system measures the cofactor:
    # the LDH/MDH dive folds PYR/OAA together with NADH, and NADH (44 atoms) won every
    # time. The ligand of interest is the one that DIFFERS between the cognate set and the
    # decoy set(s); anything shared by both is context and is never the focus.
    sets = [{l.split(":")[-1] for l in spec["ligands"]}] + \
           [{l.split(":")[-1] for l in d} for d in (spec.get("decoys") or [])]
    shared = set.intersection(*sets) if len(sets) > 1 else set()
    blinded, perseed = [], {}
    for j in jobs:
        jd = next((p for p in (workdir / "af3").glob(f"*/{j['name']}") if p.is_dir()), None)
        if jd is None:
            log(f"missing output for {j['name']}")
            continue
        focus = None
        if shared:
            diff = [c for c in (j.get("ligands") or [])
                    if c.split(":")[-1] not in shared]
            if len(diff) == 1:
                focus = diff[0].split(":")[-1]
        rows, summ = summarise_job(jd, j.get("site") or spec.get("site") or [],
                                   want_lig=focus, pb_top=True)
        if summ is None:
            continue
        perseed[j["blind_id"]] = {
            "lig_plddt": per_seed(rows, "lig_plddt"),
            "pae_min": per_seed(rows, "prot_lig_pae_min"),
            "centroid": per_seed(rows, "site_centroid_dist"),
            "recall": per_seed(rows, "site_recall"),
        }
        blinded.append({"blind_id": j["blind_id"], "ligands": ligset_tag(j["ligands"]),
                        **{k: v for k, v in summ.items() if k != "job"}})
    write_tsv(workdir / "blinded_readouts.tsv", blinded)
    (workdir / "blinded_per_seed.json").write_text(json.dumps(perseed, indent=1))
    return blinded, perseed


def gate(workdir: Path, spec: dict, prereg: dict) -> dict:
    """Positive-control gate: top-ranked ligand RMSD to the crystal pose."""
    pc = spec.get("positive_control")
    if not pc:
        return {"status": "no_positive_control", "pass": False}
    import analyze_af3 as AA
    jd = next((p for p in (workdir / "af3").glob(f"*/{spec['case']}_poscontrol")
               if p.is_dir()), None)
    if jd is None:
        return {"status": "not_run", "pass": False}
    jobs = json.loads((workdir / "jobs.json").read_text())
    pcjob = next((j for j in jobs if j["blind_id"] == "poscontrol"), {})
    rows, summ = summarise_job(jd, pcjob.get("site") or spec.get("site") or [], pb_top=False)
    top = max(rows, key=lambda r: r["ranking"] or -1)
    ref = (pc["pdb"], pc["chain"], pc["ligand"], pc.get("seqnum"), top["lig_chain"])
    refdata = AA.ref_site(*ref[:4])
    ev = AA.eval_sample(top["path"], ref, refdata)
    rmsd = ev.get("lig_rmsd")
    ok = rmsd is not None and rmsd <= prereg["gate_rmsd_max"]
    return {"status": "run", "pass": bool(ok), "lig_rmsd": rmsd,
            "centroid_dist": ev.get("lig_centroid_dist"), "site_recall": ev.get("site_recall"),
            "lig_plddt": ev.get("lig_plddt"), "threshold_A": prereg["gate_rmsd_max"],
            "ref": f"{pc['pdb']}:{pc['ligand']}"}


def d23_binds(per_seed_entry: dict, th: dict, min_seeds: int = 2) -> tuple[bool, int]:
    """Frozen D23 composite for ONE job: ligand pLDDT >= t_plddt AND (centroid <= t_centroid
    OR site recall >= t_recall), in at least `min_seeds` of 3 seeds.

    The geometry term is not optional. Without it a job where AF3 parks the ligand 30 A
    from the footprint still scores on confidence alone, which is exactly how the MetE and
    choline-kinase cases first produced a spurious directional verdict.
    """
    n = 0
    for s, pl in (per_seed_entry.get("lig_plddt") or {}).items():
        if pl is None:
            continue
        cen = (per_seed_entry.get("centroid") or {}).get(s)
        rec = (per_seed_entry.get("recall") or {}).get(s)
        geom = (cen is not None and cen <= th["t_centroid"]) or \
               (rec is not None and rec >= th["t_recall"])
        if pl >= th["t_plddt"] and geom:
            n += 1
    return n >= min_seeds, n


def decide(workdir: Path, spec: dict, prereg: dict) -> dict:
    key = json.loads((workdir / "blind_key.json").read_text())["map"]
    perseed = json.loads((workdir / "blinded_per_seed.json").read_text())
    seeds = [str(s) for s in spec.get("seeds", [1, 2, 3])]

    def pick(state: str, kind: str) -> list[str]:
        """blind ids for one state and either the cognate ligand or any decoy."""
        return [b for b, k in key.items() if k["state"] == state
                and ((k["ligand_kind"] == "cognate") if kind == "cognate"
                     else (k["ligand_kind"] != "cognate"))]

    der_c, anc_c = pick("derived", "cognate"), pick("ancestral", "cognate")
    der_d, anc_d = pick("derived", "decoy"), pick("ancestral", "decoy")
    if not (der_c and anc_c):
        return {"verdict": "inconclusive", "reason": "missing derived or ancestral cognate runs"}

    def mean_over(bids: list[str], metric: str, seed: str):
        v = [perseed[b][metric][seed] for b in bids
             if b in perseed and seed in perseed[b][metric]]
        return stats.mean(v) if v else None

    rows, passes_gain, passes_loss = [], [], []
    for sd in seeds:
        d_p, a_p = mean_over(der_c, "lig_plddt", sd), mean_over(anc_c, "lig_plddt", sd)
        d_e, a_e = mean_over(der_c, "pae_min", sd), mean_over(anc_c, "pae_min", sd)
        # within-protein cognate-vs-decoy gap, averaged over all members that have a decoy
        xp = [mean_over(c, "lig_plddt", sd) - mean_over(d, "lig_plddt", sd)
              for c, d in ((der_c, der_d), (anc_c, anc_d))
              if d and mean_over(d, "lig_plddt", sd) is not None]
        xe = [mean_over(d, "pae_min", sd) - mean_over(c, "pae_min", sd)
              for c, d in ((der_c, der_d), (anc_c, anc_d))
              if d and mean_over(d, "pae_min", sd) is not None]
        X_p, X_e = (stats.mean(xp) if xp else None), (stats.mean(xe) if xe else None)
        D_p, D_e = (d_p - a_p if None not in (d_p, a_p) else None), \
                   (a_e - d_e if None not in (d_e, a_e) else None)
        g = (None not in (D_p, D_e, X_p, X_e)) and D_p > X_p and D_e > X_e
        l = (None not in (D_p, D_e, X_p, X_e)) and (-D_p) > X_p and (-D_e) > X_e
        passes_gain.append(g)
        passes_loss.append(l)
        rows.append({"seed": sd, "derived_lig_plddt": rnd(d_p), "ancestral_lig_plddt": rnd(a_p),
                     "derived_pae_min": rnd(d_e, 3), "ancestral_pae_min": rnd(a_e, 3),
                     "D_plddt": rnd(D_p), "X_plddt_decoy_gap": rnd(X_p),
                     "D_pae": rnd(D_e, 3), "X_pae_decoy_gap": rnd(X_e, 3),
                     "gain_rule_met": g, "loss_rule_met": l})
    write_tsv(workdir / "decision_per_seed.tsv", rows)
    # --- frozen D23 composite, per member job (the concordance requirement)
    th = json.loads((P / "results/cofold/calibration/threshold.json").read_text()
                    )["d23_composite"]["thresholds"]
    binds, geom = {}, {}
    for b, k in key.items():
        if b not in perseed:
            continue
        ok, n = d23_binds(perseed[b], th)
        binds[b] = ok
        cen = [v for v in (perseed[b].get("centroid") or {}).values() if v is not None]
        rec = [v for v in (perseed[b].get("recall") or {}).values() if v is not None]
        geom[b] = {"state": k["state"], "kind": k["ligand_kind"], "acc": k["acc"],
                   "binds": ok, "n_seeds_binding": n,
                   "median_centroid_A": round(stats.median(cen), 2) if cen else None,
                   "median_site_recall": round(stats.median(rec), 3) if rec else None}
    der_binds = [b for b in der_c if binds.get(b)]
    anc_binds = [b for b in anc_c if binds.get(b)]
    any_cognate_binds = bool(der_binds or anc_binds)
    cens = [geom[b]["median_centroid_A"] for b in (der_c + anc_c)
            if geom.get(b, {}).get("median_centroid_A") is not None]
    recs = [geom[b]["median_site_recall"] for b in (der_c + anc_c)
            if geom.get(b, {}).get("median_site_recall") is not None]
    ligand_off_footprint = (not any_cognate_binds and cens
                            and min(cens) > th["t_centroid"]
                            and (not recs or max(recs) < 0.1))

    g = gate(workdir, spec, prereg)
    if not g["pass"]:
        verdict = prereg["failed_gate_label"]
    elif ligand_off_footprint:
        # AF3 never puts the cognate ligand at the atlas footprint in ANY member, on
        # either side, so the comparison says nothing about the footprint. A difference in
        # ligand pLDDT here only describes how confidently the ligand was placed SOMEWHERE
        # ELSE, and must not be read as a binding difference.
        verdict = "uninformative_ligand_not_at_footprint"
    elif all(passes_gain) and der_binds and not anc_binds:
        verdict = "gain_in_derived"
    elif all(passes_loss) and anc_binds and not der_binds:
        verdict = "loss_in_derived"
    elif any_cognate_binds and anc_binds and not der_binds:
        verdict = "contradicts_call"
    else:
        verdict = "inconclusive"
    return {"case": spec["case"], "verdict": verdict, "positive_control_gate": g,
            "d23_composite": {"thresholds": th, "per_job": geom,
                              "derived_members_binding": [geom[b]["acc"] for b in der_binds],
                              "ancestral_members_binding": [geom[b]["acc"] for b in anc_binds],
                              "any_cognate_binds": any_cognate_binds,
                              "ligand_off_footprint": bool(ligand_off_footprint)},
            "verdict_basis": "frozen D23 composite (d23_composite) gating the per-seed rule",
            "legacy_prereg_rule": {
                "status": "SUPERSEDED: computed for the record only. It ranks on ligand pLDDT and "
                          "min PAE with no geometry term, so it can call a direction from "
                          "placements far from the footprint. The verdict uses d23_composite.",
                "rule": "all seeds must meet D > X for BOTH ligand pLDDT and min P-L PAE",
                "seeds_gain": passes_gain, "seeds_loss": passes_loss},
            "per_seed": rows,
            "primary_readout": prereg["primary"]}


def rnd(x, n: int = 1):
    return None if x is None else round(x, n)


def write_tsv(path: Path, rows: list[dict]) -> None:
    """Atomic: write to a sibling temp file and rename, so a reader never sees a partial
    table and a crash mid-write cannot leave concatenated output."""
    if not rows:
        return
    keys = list(rows[0].keys())
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join("" if r.get(k) is None else str(r[k]) for k in keys) + "\n")
    tmp.replace(path)


def unblind(workdir: Path, spec: dict, prereg: dict, results_dir: Path) -> dict:
    """Join the blinded table to the key and apply the decision rule."""
    key = json.loads((workdir / "blind_key.json").read_text())["map"]
    lines = (workdir / "blinded_readouts.tsv").read_text().strip().split("\n")
    hdr = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        r = dict(zip(hdr, line.split("\t")))
        k = key.get(r["blind_id"], {})
        rows.append({"member": k.get("member", "positive_control"),
                     "state": k.get("state", "control"),
                     "ligand_kind": k.get("ligand_kind", "cognate"), **r})
    rows.sort(key=lambda r: (r["state"], r["member"], r["ligand_kind"]))
    results_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(results_dir / "case_table.tsv", rows)
    dec = decide(workdir, spec, prereg)
    (results_dir / "decision.json").write_text(json.dumps(dec, indent=1))
    for f in ("msa_source.tsv", "decision_per_seed.tsv"):
        if (workdir / f).exists():
            (results_dir / f).write_text((workdir / f).read_text())
    log(f"verdict = {dec['verdict']}; tables in {results_dir}")
    return dec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec")
    ap.add_argument("--msa", default="auto", choices=["auto", "colabfold", "cached"])
    ap.add_argument("--gpus", default="0,1,2", help="GPU 3 is reserved for another agent")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--unblind", action="store_true")
    ap.add_argument("--workdir")
    ap.add_argument("--results-dir")
    ap.add_argument("--prefetch-msas", nargs="+", metavar="CASE_YAML",
                    help="fetch ColabFold MSAs for the members these case-v1 files would "
                         "fold, then exit (no ligands needed)")
    a = ap.parse_args()
    if a.prefetch_msas:
        prefetch_msas(a.prefetch_msas)
        return 0

    if not a.spec:
        ap.error("--spec is required unless --prefetch-msas is given")
    spec = load_spec(a.spec)
    prereg = load_prereg()
    workdir = Path(a.workdir) if a.workdir else BIG / spec["case"]
    results = Path(a.results_dir) if a.results_dir else P / "results/cofold/cases" / spec["case"]
    workdir.mkdir(parents=True, exist_ok=True)

    if a.build:
        build(spec, a.msa, workdir)
    if a.run:
        run(workdir, [int(g) for g in a.gpus.split(",")])
    if a.unblind:
        collect(workdir, spec)
        dec = unblind(workdir, spec, prereg, results)
        print(json.dumps({k: v for k, v in dec.items() if k != "per_seed"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
