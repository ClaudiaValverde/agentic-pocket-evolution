#!/usr/bin/env python
"""Extract single AFDB v6 model(s) from a local proteome tar without scanning the whole tar.

Uses the atlas tar index (<atlas>/tar_index/<tar>.parquet: name, offset_data, size) for a direct seek;
falls back to tarfile scanning if the index is missing.

    python scripts/extract_afdb_member.py --tar UP000005640_9606_HUMAN_v6.tar --acc P04035 P07900 --out /tmp/afdb
"""
import argparse
import tarfile
from pathlib import Path

import pandas as pd

AFDB_V6 = Path("/cluster/work/beltrao/shared/alphafold/v6")
ATLAS = Path("/cluster/scratch/gmagnusson/pocketcompare_afdb")


def extract(tar: Path, names: list[str], out: Path) -> list[Path]:
    out.mkdir(parents=True, exist_ok=True)
    idx_path = ATLAS / "tar_index" / f"{tar.name}.parquet"
    done = []
    if idx_path.exists():
        idx = pd.read_parquet(idx_path).set_index("name")
        with open(tar, "rb") as fh:
            for n in names:
                if n not in idx.index:
                    raise KeyError(f"{n} not in {idx_path}")
                off, size = int(idx.at[n, "offset_data"]), int(idx.at[n, "size"])
                fh.seek(off)
                (out / n).write_bytes(fh.read(size))
                done.append(out / n)
    else:
        with tarfile.open(tar) as tf:
            for n in names:
                m = tf.getmember(n)
                (out / n).write_bytes(tf.extractfile(m).read())
                done.append(out / n)
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tar", required=True, help="tar file name (in AFDB v6 dir) or path")
    ap.add_argument("--acc", nargs="+", required=True, help="UniProt accessions (fragment F1) or full member names")
    ap.add_argument("--fragment", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tar = Path(a.tar) if "/" in a.tar else AFDB_V6 / a.tar
    names = [x if x.endswith(".cif.gz") else f"AF-{x}-F{a.fragment}-model_v6.cif.gz" for x in a.acc]
    for p in extract(tar, names, Path(a.out)):
        print(p)


if __name__ == "__main__":
    main()
