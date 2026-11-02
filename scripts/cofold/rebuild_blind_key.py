#!/usr/bin/env python3
"""Rebuild a case workdir's blind_key.json/jobs.json from the AF3 inputs on disk.

Needed once, because build() used to mint a fresh HMAC salt on every call: a re-score
therefore produced blind ids that no longer matched the AF3 output directories already
folded, and collect() silently found zero readouts. The salt is now persisted, but the
already-folded Tier 1 cases have outputs keyed by the ORIGINAL salt, which is gone.

Recovery does not need the salt. Each job's input JSON is still on disk and contains the
protein sequence and the ligand CCD codes, so every output directory can be matched back
to (member, ligand_kind) by exact sequence identity plus ligand composition. Blinding is
preserved: the mapping still lives only in blind_key.json.

Usage: rebuild_blind_key.py --spec cases/<id>/cofold/cofold_spec.yaml --workdir /tmp/...
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_cofold as CC


def ligtags(seqs) -> set[str]:
    out = set()
    for s in seqs:
        lig = s.get("ligand")
        if lig:
            out |= set(lig.get("ccdCodes") or ([lig["smiles"]] if "smiles" in lig else []))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    ap.add_argument("--workdir", required=True)
    a = ap.parse_args()
    spec = CC.load_spec(a.spec)
    wd = Path(a.workdir)
    seq_of = {m["label"]: CC.member_sequence(m) for m in spec["members"]}
    by_seq = {v: k for k, v in seq_of.items()}
    member = {m["label"]: m for m in spec["members"]}
    ligsets = {"cognate": {l.split(":", 1)[1] for l in spec["ligands"]}}
    for i, d in enumerate(spec.get("decoys", [])):
        ligsets[f"decoy{i+1}"] = {l.split(":", 1)[1] for l in d}

    jobs, key, unmatched = [], {}, []
    for jf in sorted(wd.glob("inputs/*.json")):
        if jf.name.endswith(".meta.json"):
            continue
        doc = json.loads(jf.read_text())
        name = doc["name"]
        outs = [p for p in wd.glob(f"af3/*/{name}") if p.is_dir()]
        if not outs:
            continue                      # input written but never folded
        prot = next((s["protein"]["sequence"] for s in doc["sequences"] if "protein" in s), None)
        tags = ligtags(doc["sequences"])
        bid = name.split("_")[-1]
        if name.endswith("_poscontrol"):
            jobs.append({"blind_id": "poscontrol", "name": name, "input": str(jf),
                         "length": len(prot or ""), "ligands": sorted(tags),
                         "site": spec.get("site") or []})
            continue
        lab = by_seq.get(prot)
        kind = next((k for k, v in ligsets.items() if v == tags), None)
        if lab is None or kind is None:
            unmatched.append((name, lab, sorted(tags)))
            continue
        m = member[lab]
        jobs.append({"blind_id": bid, "name": name, "input": str(jf),
                     "length": len(prot), "ligands": sorted(tags),
                     "site": m.get("site") or spec.get("site") or []})
        key[bid] = {"member": lab, "acc": m["accession"], "state": m["state"],
                    "ligand_kind": kind, "ligands": "+".join(sorted(tags)),
                    "site_source": m.get("site_source", "spec"),
                    "n_site": len(m.get("site") or spec.get("site") or []),
                    "recovered": "matched by sequence + ligand composition"}
    tmp = wd / "jobs.json.tmp"; tmp.write_text(json.dumps(jobs, indent=1)); tmp.replace(wd / "jobs.json")
    old = json.loads((wd / "blind_key.json").read_text()) if (wd / "blind_key.json").exists() else {}
    tmp = wd / "blind_key.json.tmp"
    tmp.write_text(json.dumps({"salt_hex": old.get("salt_hex", "recovered"),
                               "map": key, "note": __doc__.strip().split("\n")[0]}, indent=1))
    tmp.replace(wd / "blind_key.json")
    print(f"recovered {len(key)} member jobs + "
          f"{sum(j['blind_id']=='poscontrol' for j in jobs)} positive control; "
          f"states={sorted({v['state'] for v in key.values()})} "
          f"kinds={sorted({v['ligand_kind'] for v in key.values()})}")
    if unmatched:
        print(f"UNMATCHED ({len(unmatched)}): {unmatched[:4]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
