#!/usr/bin/env python
"""Verify the AFDB fragment numbering rule used by pocketevo.sites against the model files.

For every structure with fragment > 1 (plus its fragment 1) in the snapshot units of the given tiers, read the
mmCIF ``_struct_ref_seq`` category from the local AFDB tar (direct seek via the atlas tar_index; only the
first part of each gzipped member is decompressed) and check
    db_align_beg - pdbx_auth_seq_align_beg == 200*(k-1)          (offset rule)
    db_align_end - db_align_beg + 1 <= 1400, == 1400 unless last  (window rule)
Writes a TSV of all checked fragments and prints a summary; exits non-zero on any mismatch.

    python scripts/verify_fragment_offsets.py --snapshot data/snapshot/20261030T085136Z --tiers 1 \
        --out results/sites/20261030T085136Z/fragment_offsets.tsv
"""
from __future__ import annotations

import argparse
import json
import sys
import zlib
from pathlib import Path

import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pocketevo.sites import FRAGMENT_STEP, FRAGMENT_WINDOW, struct_ref_seq  # noqa: E402


def read_header(fh, off: int, size: int, nbytes: int = 1 << 16) -> str:
    """Decompress just enough of a gzip member to get the mmCIF header (``_struct_ref_seq`` precedes atoms)."""
    fh.seek(off)
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    out = b""
    pos = 0
    while pos < size and b"_struct_ref_seq.seq_align_end" not in out:
        chunk = fh.read(min(nbytes, size - pos))
        pos += len(chunk)
        out += d.decompress(chunk)
    return out.decode(errors="replace")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", required=True)
    ap.add_argument("--tiers", default="1")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    snap = Path(a.snapshot)
    mj = json.loads((snap / "manifest.json").read_text())
    atlas = Path(mj["atlas_root"])
    man = pd.read_parquet(snap / "manifest.parquet")
    man = man[man.tier.isin([int(x) for x in a.tiers.split(",")])]
    wl = pd.read_parquet(atlas / "worklist.parquet", columns=["unit_id", "sources"]).set_index("unit_id")
    rows = []
    for _, u in man.iterrows():
        st = pq.read_table(atlas / u.structures_path, columns=["structure_stem", "accession", "fragment", "n_res"])
        multi = set(pc.filter(st.column("accession"), pc.greater(st.column("fragment"), 1)).to_pylist())
        if not multi:
            continue
        st = st.to_pandas()
        st = st[st.accession.isin(multi)]
        for tar in wl.at[u.unit_id, "sources"]:
            idx = pd.read_parquet(atlas / "tar_index" / f"{Path(tar).name}.parquet").set_index("name")
            with open(tar, "rb") as fh:
                for r in st.itertuples():
                    name = f"{r.structure_stem}.cif.gz"
                    if name not in idx.index:
                        continue
                    d = struct_ref_seq(read_header(fh, int(idx.at[name, "offset_data"]), int(idx.at[name, "size"])))
                    rows.append(dict(unit_id=u.unit_id, accession=r.accession, fragment=r.fragment, n_res=r.n_res,
                                     db_align_beg=int(d["db_align_beg"]), db_align_end=int(d["db_align_end"]),
                                     auth_beg=int(d["pdbx_auth_seq_align_beg"]),
                                     auth_end=int(d["pdbx_auth_seq_align_end"])))
    df = pd.DataFrame(rows)
    df["offset_mmcif"] = df.db_align_beg - df.auth_beg
    df["offset_rule"] = FRAGMENT_STEP * (df.fragment - 1)
    df["window"] = df.db_align_end - df.db_align_beg + 1
    last = df.groupby("accession").fragment.transform("max") == df.fragment
    df["ok_offset"] = df.offset_mmcif == df.offset_rule
    df["ok_window"] = (df.window <= FRAGMENT_WINDOW) & ((df.window == FRAGMENT_WINDOW) | last) & \
        (df.auth_end - df.auth_beg + 1 == df.window) & (df.n_res == df.window)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.sort_values(["accession", "fragment"]).to_csv(a.out, sep="\t", index=False)
    n_bad = int((~df.ok_offset).sum() + (~df.ok_window).sum())
    print(f"checked {len(df)} fragments of {df.accession.nunique()} accessions: "
          f"offset mismatches {int((~df.ok_offset).sum())}, window mismatches {int((~df.ok_window).sum())}")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
