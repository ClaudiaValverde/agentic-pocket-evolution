"""Per-residue sequence + pLDDT from AFDB v6 mmCIF models (tier 1/2 tars, read via the atlas tar_index offsets).

Only two mmCIF blocks are parsed (no atom records): `_entity_poly.pdbx_seq_one_letter_code_can` (model sequence)
and the `_ma_qa_metric_local` loop (per-residue pLDDT, = the CA B-factor in AFDB files), plus the
`_ma_target_ref_db_details` key/values (UniProt accession, CRC64 checksum, align begin/end).

Storage (compact cache, one parquet per tar): pLDDT is stored as `plddt_q = floor(2 * pLDDT)` in uint8 (half-pLDDT
units, 0..200). This is exact for every threshold on an integer or half-integer pLDDT (v < 70 <=> q < 140); means
are estimated as mean(q)/2 + 0.25 (max error 0.25 pLDDT). See `plddt_from_q`.

Fragments (AFDB human proteins > 2,700 aa): fragment k covers UniProt positions 200(k-1)+1 .. 200(k-1)+1400
(verified, results/sites/.../SITES.md). The combined per-position value is taken from the covering fragment in
which the position is farthest from a fragment edge (ties -> lower fragment); sequences of overlapping fragments
must agree (else flagged `fragment_seq_conflict`).
"""
from __future__ import annotations

import hashlib
import re
import zlib

import numpy as np

FRAGMENT_STEP = 200
FRAGMENT_WINDOW = 1400
Q_SCALE = 2  # plddt_q = floor(Q_SCALE * pLDDT)

_THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
          "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
          "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "UNK": "X"}
_REFDB_RE = re.compile(rb"^_ma_target_ref_db_details\.(\S+)\s+(.*?)\s*$", re.M)


def plddt_from_q(q: np.ndarray) -> np.ndarray:
    """Midpoint estimate of pLDDT from the uint8 half-unit code."""
    return q.astype(np.float32) / Q_SCALE + 0.5 / Q_SCALE


def q_threshold(plddt: float) -> int:
    """q code such that (pLDDT < plddt) <=> (q < q_threshold(plddt)) for integer/half-integer thresholds."""
    return int(round(plddt * Q_SCALE))


def parse_cif_text(txt: bytes) -> dict:
    """Return {'sequence', 'plddt_q' (uint8), 'refdb': {...}} from raw (decompressed) mmCIF bytes."""
    # sequence
    key = b"_entity_poly.pdbx_seq_one_letter_code_can"
    i = txt.find(key)
    if i < 0:
        raise ValueError("no pdbx_seq_one_letter_code_can")
    eol = txt.find(b"\n", i)
    inline = txt[i + len(key):eol].strip()
    if inline:  # short sequences are written on the key line (optionally quoted)
        seq = inline.strip(b"'\"").decode()
    else:       # multi-line text field: ';' at line start ... '\n;'
        if txt[eol + 1:eol + 2] != b";":
            raise ValueError("unexpected pdbx_seq_one_letter_code_can layout")
        k = txt.find(b"\n;", eol + 1)
        seq = re.sub(rb"\s+", b"", txt[eol + 2:k]).decode()
    # local pLDDT loop
    hdr = b"_ma_qa_metric_local."
    i = txt.find(hdr)
    if i < 0:
        raise ValueError("no _ma_qa_metric_local loop")
    cols = []
    pos = i
    while txt.startswith(hdr, pos):
        e = txt.find(b"\n", pos)
        cols.append(txt[pos + len(hdr):e].strip().decode())
        pos = e + 1
    end = txt.find(b"\n#", pos - 1)
    body = txt[pos:end].split()
    nc = len(cols)
    arr = np.array(body, dtype=object).reshape(-1, nc)
    ci_seq, ci_val, ci_comp = cols.index("label_seq_id"), cols.index("metric_value"), cols.index("label_comp_id")
    seq_id = arr[:, ci_seq].astype(np.int64)
    val = arr[:, ci_val].astype(np.float64)
    n = len(seq)
    if len(seq_id) != n or not np.array_equal(seq_id, np.arange(1, n + 1)):
        raise ValueError(f"pLDDT rows {len(seq_id)} != sequence length {n} or not contiguous")
    comp = "".join(_THREE.get(c.decode() if isinstance(c, bytes) else c, "X") for c in arr[:, ci_comp])
    if comp != seq:
        raise ValueError("pLDDT loop residue names disagree with entity sequence")
    q = np.floor(np.clip(val, 0, 100) * Q_SCALE).astype(np.uint8)
    refdb = {m.group(1).decode(): m.group(2).decode().strip('"') for m in _REFDB_RE.finditer(txt)}
    return {"sequence": seq, "plddt_q": q, "refdb": refdb}


def parse_cif_gz(blob: bytes) -> dict:
    return parse_cif_text(zlib.decompress(blob, 16 + zlib.MAX_WBITS))


def combine_fragments(frags: dict[int, tuple[str, np.ndarray]]) -> tuple[str, np.ndarray, bool]:
    """frags: {k: (seq, q)}. Returns (sequence, q, conflict) in UniProt numbering."""
    if list(frags) == [1]:
        s, q = frags[1]
        return s, q, False
    L = max(FRAGMENT_STEP * (k - 1) + len(s) for k, (s, _) in frags.items())
    seq = np.full(L, "", dtype="<U1")
    best_d = np.full(L, -1, dtype=np.int64)
    out = np.zeros(L, dtype=np.uint8)
    conflict = False
    for k in sorted(frags):
        s, q = frags[k]
        off = FRAGMENT_STEP * (k - 1)
        idx = np.arange(len(s))
        d = np.minimum(idx, len(s) - 1 - idx)
        sl = slice(off, off + len(s))
        sa = np.array(list(s))
        seen = seq[sl] != ""
        if np.any(seq[sl][seen] != sa[seen]):
            conflict = True
        better = d > best_d[sl]
        seq_sl = seq[sl]; seq_sl[~seen] = sa[~seen]; seq[sl] = seq_sl
        o = out[sl]; o[better] = q[better]; out[sl] = o
        b = best_d[sl]; b[better] = d[better]; best_d[sl] = b
    if np.any(seq == ""):
        conflict = True
    return "".join(seq.tolist()), out, conflict


def seq_md5(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()


# CRC64 (UniProt/SwissProt ISO 3309 polynomial), used to compare against `seq_db_sequence_checksum`
_CRC_TABLE = None


def crc64(seq: str) -> str:
    global _CRC_TABLE
    if _CRC_TABLE is None:
        poly = 0xD800000000000000
        t = []
        for i in range(256):
            c = i
            for _ in range(8):
                c = (c >> 1) ^ poly if c & 1 else c >> 1
            t.append(c)
        _CRC_TABLE = t
    crc = 0
    for ch in seq.encode():
        crc = _CRC_TABLE[(crc ^ ch) & 0xFF] ^ (crc >> 8)
    return f"{crc:016X}"
