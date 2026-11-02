#!/usr/bin/env python3
"""Summarise AF3 benchmark runs: runtime/memory and ligand pocket placement.

For each job dir under results/cofold/bench/<job>/ it parses run_af3.log
(wall time, per-seed inference time, peak GPU memory) and, for jobs with a
reference holo PDB, every sample model:
  * superpose predicted protein chain A onto the holo chain (CA atoms of
    sequence-aligned residues),
  * site = holo residues with any heavy atom within 5 A of the crystal ligand,
  * lig_centroid_dist: distance between predicted and crystal ligand centroids,
  * lig_rmsd: heavy-atom RMSD over atoms with matching names (CCD ligands only;
    no symmetry correction),
  * site_recall: fraction of site residues with a heavy atom within 5 A of the
    predicted ligand,
  * lig_plddt: mean pLDDT (B-factor) over predicted ligand atoms.
Decoy jobs are evaluated against the cognate ligand's site, i.e. "does AF3 put
the decoy in the known pocket?".

Usage: analyze_af3.py [--bench-dir DIR]
Writes runtime.tsv, placement_samples.tsv, placement_summary.tsv to --out-dir.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import statistics as stats
from pathlib import Path

import gemmi
import numpy as np
import zstandard

P = Path(__file__).resolve().parents[2]  # repo root
PDB = "/cluster/scratch/jjaenes/pdb"

# job-name prefix -> (pdb id, holo chain, crystal ligand resname, crystal ligand seqnum or None,
#                    predicted query-ligand chain)
REFS = {
    # MSA-source test (docs/COFOLD.md Sec.2.1): same protein+ligand, no templates,
    # three MSA provenances.
    "msatest_cdk2_cf": ("1hck", "A", "ATP", 400, "B"),            # ColabFold unpaired
    "msatest_cdk2_cacheunpaired": ("1hck", "A", "ATP", 400, "B"),  # cached unpaired only
    "cdk2_atp_mg_notmpl": ("1hck", "A", "ATP", 400, "B"),
    "cdk2_atp_mg": ("1hck", "A", "ATP", 400, "B"),
    "cdk2_decoy_azm": ("1hck", "A", "ATP", 400, "B"),
    "cdk2_decoy_cff": ("1hck", "A", "ATP", 400, "B"),
    "ca2_zn_azmsmiles": ("3hs4", "A", "AZM", None, "C"),
    "ca2_zn_azm": ("3hs4", "A", "AZM", None, "C"),
    "ca2_zn_decoy_atp": ("3hs4", "A", "AZM", None, "C"),
    "ca2_zn_decoy_cff": ("3hs4", "A", "AZM", None, "C"),
    "yADH1_nad_zn": ("4w6z", "A", "8ID", 403, "B"),
    "yADH1_decoy_cff": ("4w6z", "A", "8ID", 403, "D"),
    "pygl_plp": ("1fa9", "A", "PLP", 999, "B"),
}


def read_model(path: str) -> gemmi.Structure:
    raw = open(path, "rb").read()
    if path.endswith(".zst"):
        raw = zstandard.ZstdDecompressor().decompress(raw, max_output_size=2**31)
    doc = gemmi.cif.read_string(raw.decode())
    st = gemmi.make_structure_from_block(doc.sole_block())
    st.setup_entities()
    return st


def heavy(res):
    return [a for a in res if a.element.name != "H"]


def ref_site(pdb, chain, lig, num):
    st = gemmi.read_structure(f"{PDB}/{pdb[1:3]}/{pdb}.cif.gz")
    st.remove_hydrogens()
    ch = st[0][chain]
    cands = [r for r in ch if r.name == lig and (num is None or r.seqid.num == num)]
    if num is None:  # pick ligand copy closest to the catalytic zinc (CA2) / first copy
        zn = [r for r in ch if r.name == "ZN"]
        if zn:
            zpos = zn[0][0].pos
            cands.sort(key=lambda r: min(a.pos.dist(zpos) for a in r))
    ligres = cands[0]
    lpos = np.array([a.pos.tolist() for a in heavy(ligres)])
    poly = [r for r in ch if r.het_flag != "H" and gemmi.find_tabulated_residue(r.name)
            and gemmi.find_tabulated_residue(r.name).is_amino_acid()]
    site = []
    for i, r in enumerate(poly):
        xyz = np.array([a.pos.tolist() for a in heavy(r)])
        if np.min(np.linalg.norm(xyz[:, None] - lpos[None], axis=-1)) < 5.0:
            site.append(i)
    return poly, ligres, site


def ca_map(poly_ref, poly_pred):
    s1 = gemmi.one_letter_code([r.name for r in poly_ref])
    s2 = gemmi.one_letter_code([r.name for r in poly_pred])
    aln = gemmi.align_string_sequences(list(s1), list(s2), [])
    m, i, j = {}, 0, 0
    for a, b in zip(aln.add_gaps(s1, 1), aln.add_gaps(s2, 2)):
        if a != "-" and b != "-":
            m[i] = j
        i += a != "-"
        j += b != "-"
    return m


def eval_sample(model_path, ref, refdata):
    pdb, chain, lig, num, pred_lig_chain = ref
    poly_ref, ligres, site = refdata
    st = read_model(model_path)
    pred_poly = list(st[0]["A"])
    m = ca_map(poly_ref, pred_poly)
    pr, pp = [], []
    for i, j in m.items():
        if "CA" in [a.name for a in poly_ref[i]] and "CA" in [a.name for a in pred_poly[j]]:
            pr.append(poly_ref[i]["CA"][0].pos)
            pp.append(pred_poly[j]["CA"][0].pos)
    sup = gemmi.superpose_positions(pr, pp)  # transform pred -> ref
    st[0].transform_pos_and_adp(sup.transform)
    plig = st[0][pred_lig_chain][0]
    patoms = heavy(plig)
    ppos = np.array([a.pos.tolist() for a in patoms])
    rpos = np.array([a.pos.tolist() for a in heavy(ligres)])
    cdist = float(np.linalg.norm(ppos.mean(0) - rpos.mean(0)))
    rnames = {a.name: np.array(a.pos.tolist()) for a in heavy(ligres)}
    common = [(np.array(a.pos.tolist()), rnames[a.name]) for a in patoms if a.name in rnames]
    rmsd = (float(np.sqrt(np.mean([np.sum((x - y) ** 2) for x, y in common])))
            if plig.name == ligres.name or (lig == "8ID" and plig.name == "NAD") else None)
    recall = 0
    for i in site:
        if i not in m:
            continue
        xyz = np.array([a.pos.tolist() for a in heavy(st[0]["A"][m[i]])])
        if np.min(np.linalg.norm(xyz[:, None] - ppos[None], axis=-1)) < 5.0:
            recall += 1
    return dict(ca_rmsd=round(sup.rmsd, 2), lig_resname=plig.name, lig_centroid_dist=round(cdist, 2),
                lig_rmsd=None if rmsd is None else round(rmsd, 2), n_common_atoms=len(common),
                site_recall=round(recall / max(1, len(site)), 2), n_site=len(site),
                lig_plddt=round(float(np.mean([a.b_iso for a in patoms])), 1))


def parse_log(log):
    txt = Path(log).read_text()
    seeds = [(int(s), float(t)) for s, t in
             re.findall(r"Running model inference with seed (\d+) took ([\d.]+) seconds", txt)]
    feat = [float(t) for t in re.findall(r"Featurising data with seed \d+ took ([\d.]+) seconds", txt)]
    tail = re.search(r"# rc=(\d+) wall_s=([\d.]+) peak_gpu_mem_mib=(\d+)", txt)
    oom = "RESOURCE_EXHAUSTED" in txt or "Out of memory" in txt
    return seeds, feat, tail, oom


def n_tokens(job_dir):
    """Tokens = polymer residues + ligand heavy atoms (AF3 tokenisation), from the top model."""
    f = [x for x in glob.glob(f"{job_dir}/*/*_model.cif*") if "seed-" not in os.path.basename(os.path.dirname(x))]
    if not f:
        return None, None
    st = read_model(f[0])
    ntok, nres = 0, None
    for ch in st[0]:
        poly = ch.get_polymer()
        if len(poly):
            ntok += len(poly)
            nres = nres or len(poly)
        else:
            ntok += sum(len(heavy(r)) for r in ch)
    return ntok, nres


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench-dir", default="/tmp/pe_big/cofold/bench",
                    help="full AF3 outputs (bulky, regenerable)")
    ap.add_argument("--out-dir", default=str(P / "results/cofold/bench"), help="where TSVs go")
    a = ap.parse_args()
    bench = Path(a.bench_dir)
    outd = Path(a.out_dir); outd.mkdir(parents=True, exist_ok=True)
    rt_rows, samp_rows, summ_rows = [], [], []
    for jd in sorted(p for p in bench.iterdir() if p.is_dir()):
        job = jd.name
        log = jd / "run_af3.log"
        if not log.exists():
            continue
        seeds, feat, tail, oom = parse_log(log)
        ntok, nres = n_tokens(jd)
        first = seeds[0][1] if seeds else None
        later = [t for _, t in seeds[1:]]
        rt_rows.append(dict(job=job, residues=nres, tokens=ntok, n_seeds=len(seeds),
                            first_seed_s=first,
                            later_seed_s=round(stats.mean(later), 1) if later else None,
                            featurise_s=round(sum(feat), 1), wall_s=tail.group(2) if tail else None,
                            peak_gpu_mib=tail.group(3) if tail else None,
                            rc=tail.group(1) if tail else None, oom=oom))
        ref = next((REFS[k] for k in sorted(REFS, key=len, reverse=True) if job.startswith(k)), None)
        if ref is None:
            continue
        refdata = ref_site(*ref[:4])
        rank_csv = glob.glob(f"{jd}/*/*_ranking_scores.csv")
        ranks = {}
        if rank_csv:
            for line in open(rank_csv[0]).read().split("\n")[1:]:
                if line:
                    s, i, r = line.split(",")
                    ranks[(int(s), int(i))] = float(r)
        rows = []
        for mp in sorted(glob.glob(f"{jd}/*/seed-*_sample-*/*_model.cif*")):
            mm = re.search(r"seed-(\d+)_sample-(\d+)", mp)
            sd, si = int(mm.group(1)), int(mm.group(2))
            sc = json.load(open(glob.glob(os.path.join(os.path.dirname(mp), "*summary_confidences.json"))[0]))
            chains = sc["chain_ids"]
            li = chains.index(ref[4])
            r = dict(job=job, seed=sd, sample=si, ranking=ranks.get((sd, si)), iptm=sc["iptm"],
                     ptm=sc["ptm"], lig_chain_iptm=sc["chain_iptm"][li],
                     prot_lig_pair_iptm=sc["chain_pair_iptm"][0][li],
                     prot_lig_pae_min=sc["chain_pair_pae_min"][0][li],
                     **eval_sample(mp, ref, refdata))
            rows.append(r)
        samp_rows += rows
        if rows:
            top = max(rows, key=lambda r: r["ranking"] or -1)
            summ_rows.append(dict(job=job, ref=f"{ref[0]}:{ref[2]}", n_samples=len(rows),
                                  top_ranking=top["ranking"], top_iptm=top["iptm"],
                                  top_pair_iptm=top["prot_lig_pair_iptm"],
                                  top_pae_min=top["prot_lig_pae_min"],
                                  top_lig_plddt=top["lig_plddt"],
                                  top_centroid_dist=top["lig_centroid_dist"],
                                  top_lig_rmsd=top["lig_rmsd"], top_site_recall=top["site_recall"],
                                  frac_centroid_lt2A=round(np.mean([r["lig_centroid_dist"] < 2 for r in rows]), 2),
                                  frac_in_site=round(np.mean([r["site_recall"] >= 0.5 for r in rows]), 2),
                                  median_pair_iptm=stats.median(r["prot_lig_pair_iptm"] for r in rows)))

    def write(name, rows):
        if not rows:
            return
        keys = list(rows[0].keys())
        with open(outd / name, "w") as f:
            f.write("\t".join(keys) + "\n")
            for r in rows:
                f.write("\t".join("" if r[k] is None else str(r[k]) for k in keys) + "\n")
        print(f"== {name}")
        print(open(outd / name).read())

    write("runtime.tsv", rt_rows)
    write("placement_samples.tsv", samp_rows)
    write("placement_summary.tsv", summ_rows)


if __name__ == "__main__":
    main()
