"""Observation states per SOG x member: DESIGN.md §3.3, config/prereg.yaml `states`.

P  present        : a site of the member with column Jaccard(site, footprint) >= tau.
W  weak-present   : no P, but a rescue-setting / re-described filtered pocket overlaps the footprint in
                    >= k settings (`apply_rescue`, hook; the atlas stores only kept pockets, so W needs
                    re-description — produced by the rescue stage, not here).
A  absent-observed: none of the A_primary criteria fails (aligned >= 0.8, mean pLDDT >= 90, frac(<70) <= 0.2,
                    length ratio in [0.7, 1.3] x OG median, no terminal truncation, no SP/TM overlap when
                    annotated, model status ok, model sequence == source sequence).
                    `state_s70` = the same rule with the pLDDT >= 70 secondary stratum.
U  uninformative  : any criterion fails (first failing criterion -> `reason`; all of them -> `fail_bits`).
G  gene absent    : species of the OG's root clade (orthoxml TaxRange) without any member (caller adds them).

Interpretation choices (documented in results/sogs/.../SOGS.md, flagged to the PI):
* a member whose best site overlaps the footprint but with Jaccard < tau, or that has a Leiden-member site of the
  SOG that fails the Jaccard rule, is U(`partial_site_overlap`), not A: a pocket is present at the position, so
  absence is not observed.
* members with a model-sequence mismatch are U even if they have a matching site (prereg: mismatch -> U).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .afmodel import Q_SCALE, q_threshold
from .sog import MSA, footprint_jaccard

# failure bits, in the order used to pick `reason` (first failing)
FAIL_BITS = [
    ("model_sequence_mismatch", 1 << 0),
    ("model_status", 1 << 1),
    ("partial_site_overlap", 1 << 2),
    ("aligned<0.8", 1 << 3),
    ("terminal_truncation", 1 << 4),
    ("length_ratio", 1 << 5),
    ("plddt<90", 1 << 6),
    ("frac_lt70>0.2", 1 << 7),
    ("sp_tm_overlap", 1 << 8),
    ("metal_only_site", 1 << 9),
    ("likely_oligomeric_interface", 1 << 10),
    ("footprint_residues<=10", 1 << 11),   # D11 (DESIGN §12): A needs the member's footprint to span > N residues
]
BIT = dict(FAIL_BITS)


def reason_names(params: dict) -> dict:
    """Reason strings with the prereg thresholds substituted (labels only)."""
    return {
        "aligned<0.8": f"aligned<{params['A_aligned_frac']:g}",
        "plddt<90": f"plddt<{params['A_plddt']:g}",
        "frac_lt70>0.2": f"frac_lt70>{params['A_frac_lt70']:g}",
        "footprint_residues<=10": f"footprint_residues<={params.get('A_min_footprint_residues', 0):g}",
    }


def aligned_matrix(msa: MSA, values: list, fill=0, dtype=np.uint8) -> np.ndarray:
    """Per-residue arrays (UniProt order) placed into MSA columns; gaps get `fill`."""
    out = np.full(msa.mat.shape, fill, dtype=dtype)
    for r, v in enumerate(values):
        if v is None:
            continue
        v = np.asarray(v)
        out[r, msa.pos2col[r][:len(v)]] = v[:len(msa.pos2col[r])]
    return out


def segments_mask(segs: str | None, n: int) -> np.ndarray | None:
    """'1-20;45-67' -> bool mask over UniProt positions 1..n. None if not annotated (unknown)."""
    if segs is None or (isinstance(segs, float) and np.isnan(segs)):
        return None
    m = np.zeros(n, dtype=bool)
    for part in str(segs).split(";"):
        if not part:
            continue
        a, b = part.split("-")
        a, b = max(int(a), 1), min(int(b), n)
        if b >= a:
            m[a - 1:b] = True
    return m


def compute_states(msa: MSA, member: pd.DataFrame, F: sp.csr_matrix, S: sp.csr_matrix, site_member: np.ndarray,
                   site_ids: np.ndarray, labels: np.ndarray, tau: float, params: dict,
                   sog_flags: dict | None = None) -> pd.DataFrame:
    """States for every SOG (rows of F) x MSA member (rows of msa).

    member: DataFrame aligned with msa.ids with columns
      n_res (int), plddt_q (uint8 arrays, UniProt order), seq_ok (bool), status_ok (bool),
      sptm (bool mask over UniProt positions or None = not annotated), og_median_len (float, same for all rows)
    site_member: MSA row of each site (rows of S); labels: SOG of each site (Leiden), for `in_sog`.
    sog_flags: {sog_index: set of SOG-level U reasons} (e.g. metal_only_site) — hook, may be None.
    """
    n, L = msa.mat.shape
    n_sogs = F.shape[0]
    if n_sogs == 0:
        return pd.DataFrame()
    NG = msa.nongap
    Q = aligned_matrix(msa, list(member.plddt_q), fill=0)
    q90 = q_threshold(params["A_plddt"])
    q70s = q_threshold(params["A_plddt_secondary"])
    q70 = q_threshold(70.0)  # frac(<70) criterion: the residue cut-off is part of the prereg key name
    Ff = F.astype(np.float64)
    fp_size = np.asarray(F.sum(axis=1)).ravel().astype(np.float64)
    n_al = np.asarray(Ff @ NG.T.astype(np.float64))                          # (sogs, members)
    sum_q = np.asarray(Ff @ (Q.astype(np.float64) * NG).T)
    n_lt70 = np.asarray(Ff @ ((Q < q70) & NG).T.astype(np.float64))
    sptm_known = np.array([m is not None for m in member.sptm])
    SPTM = aligned_matrix(msa, [m if m is not None else None for m in member.sptm], fill=0, dtype=np.uint8)
    n_sptm = np.asarray(Ff @ SPTM.T.astype(np.float64))
    with np.errstate(invalid="ignore", divide="ignore"):
        al_frac = n_al / fp_size[:, None]
        mean_pl = np.where(n_al > 0, sum_q / np.maximum(n_al, 1) / Q_SCALE + 0.5 / Q_SCALE, np.nan)   # reported
        # conservative lower bound (stored pLDDT = floor(2 pLDDT)/2): used for the A thresholds, so a true mean
        # in [90.0, 90.25) from quantisation can never pass as >= 90
        mean_pl_lb = np.where(n_al > 0, sum_q / np.maximum(n_al, 1) / Q_SCALE, np.nan)
        frac70 = np.where(n_al > 0, n_lt70 / np.maximum(n_al, 1), np.nan)
    # terminal truncation: a footprint column before the member's first / after its last residue
    first = np.array([p[0] if len(p) else L for p in msa.pos2col])
    last = np.array([p[-1] if len(p) else -1 for p in msa.pos2col])
    Fc = F.tocsr()
    fmin = np.array([Fc.indices[Fc.indptr[s]:Fc.indptr[s + 1]].min() if Fc.indptr[s + 1] > Fc.indptr[s] else 0
                     for s in range(n_sogs)])
    fmax = np.array([Fc.indices[Fc.indptr[s]:Fc.indptr[s + 1]].max() if Fc.indptr[s + 1] > Fc.indptr[s] else -1
                     for s in range(n_sogs)])
    trunc = (fmin[:, None] < first[None, :]) | (fmax[:, None] > last[None, :])
    # best site per (sog, member)
    best_j = np.zeros((n_sogs, n), dtype=np.float64)
    best_site = np.full((n_sogs, n), -1, dtype=np.int64)
    if S.shape[0]:
        J = footprint_jaccard(F, S)
        order = np.lexsort((J.col, -J.data))  # highest J first, ties -> lowest site index
        r, c, d = J.row[order], J.col[order], J.data[order]
        m = site_member[c]
        key = r * n + m
        _, first_idx = np.unique(key, return_index=True)
        best_j[r[first_idx], m[first_idx]] = d[first_idx]
        best_site[r[first_idx], m[first_idx]] = c[first_idx]
    in_sog = np.zeros((n_sogs, n), dtype=bool)
    if len(labels):
        in_sog[labels, site_member] = True
    lr = member.n_res.to_numpy(np.float64) / float(member.og_median_len.iloc[0])
    len_bad = (lr < params["A_len_lo"]) | (lr > params["A_len_hi"])

    bits = np.zeros((n_sogs, n), dtype=np.int32)
    bits |= np.where(~member.seq_ok.to_numpy(bool), BIT["model_sequence_mismatch"], 0)[None, :]
    bits |= np.where(~member.status_ok.to_numpy(bool), BIT["model_status"], 0)[None, :]
    # a site overlapping the footprint below tau, or a Leiden-member site of this SOG: absence is not observed
    bits |= np.where(((best_j > 0) | in_sog) & (best_j < tau - 1e-12), BIT["partial_site_overlap"], 0)
    bits |= np.where(~(al_frac >= params["A_aligned_frac"] - 1e-12), BIT["aligned<0.8"], 0)
    bits |= np.where(trunc, BIT["terminal_truncation"], 0)
    if params.get("A_min_footprint_residues") is not None:   # D11: member footprint must span > N residues
        bits |= np.where(~(n_al > params["A_min_footprint_residues"]), BIT["footprint_residues<=10"], 0)
    bits |= np.where(len_bad, BIT["length_ratio"], 0)[None, :]
    sum_q_ok = n_al > 0
    pl_bad = ~(sum_q_ok & (mean_pl_lb >= params["A_plddt"] - 1e-9))
    pl_bad_s = ~(sum_q_ok & (mean_pl_lb >= params["A_plddt_secondary"] - 1e-9))
    f70_bad = ~(sum_q_ok & (frac70 <= params["A_frac_lt70"] + 1e-12))
    bits |= np.where(f70_bad, BIT["frac_lt70>0.2"], 0)
    sptm_bad = (n_sptm / np.maximum(fp_size[:, None], 1)) > params["A_sptm_overlap_max"] + 1e-12
    bits |= np.where(sptm_bad, BIT["sp_tm_overlap"], 0)
    if sog_flags:
        for s, fl in sog_flags.items():
            for f in fl:
                bits[s, :] |= BIT[f]
    bits_primary = bits | np.where(pl_bad, BIT["plddt<90"], 0)
    bits_second = bits | np.where(pl_bad_s, BIT["plddt<90"], 0)

    present = (best_j >= tau - 1e-12) & member.seq_ok.to_numpy(bool)[None, :]
    state = np.where(present, "P", np.where(bits_primary == 0, "A", "U"))
    state_s70 = np.where(present, "P", np.where(bits_second == 0, "A", "U"))
    names = reason_names(params)

    def first_reason(b):
        for nm, bit in FAIL_BITS:
            if b & bit:
                return names.get(nm, nm)
        return ""
    uniq = {b: first_reason(b) for b in np.unique(bits_primary)}
    reason = np.vectorize(uniq.get, otypes=[object])(bits_primary)
    # P rows: reason empty except a mismatch (which forces U above)
    reason = np.where(state == "P", "", reason)
    ss = np.repeat(np.arange(n_sogs), n)
    mm = np.tile(np.arange(n), n_sogs)
    bs = best_site.ravel()
    out = pd.DataFrame({
        "sog_index": ss.astype(np.int32),
        "accession": np.asarray(msa.ids, dtype=object)[mm],
        "state": state.ravel(),
        "state_s70": state_s70.ravel(),
        "reason": reason.ravel(),
        "fail_bits": np.where(state.ravel() == "P", 0, bits_primary.ravel()).astype(np.int16),
        "site_id": np.where(bs >= 0, np.asarray(site_ids, dtype=object)[np.maximum(bs, 0)] if len(site_ids) else None, None),
        "jaccard": best_j.ravel().astype(np.float32),
        "in_sog": in_sog.ravel(),
        "footprint_plddt": mean_pl.ravel().astype(np.float32),
        "aligned_frac": al_frac.ravel().astype(np.float32),
        "frac_lt70": frac70.ravel().astype(np.float32),
        "len_ratio": np.tile(lr, n_sogs).astype(np.float32),
        "sptm_checked": np.tile(sptm_known, n_sogs),
    })
    out.loc[out.jaccard == 0, "site_id"] = None
    return out


def apply_rescue(states: pd.DataFrame, rescue: pd.DataFrame, k: int) -> pd.DataFrame:
    """Hook for the rescue stage (re-description is produced by another agent; not run here).

    rescue: one row per (sog_id, accession) that was re-described, with `n_settings_hit` = number of the prereg
    rescue settings (6: `--min-score none` [= atlas-filtered pockets], no AF mode, seeds s1/s2, fast/accurate)
    in which a site overlaps the footprint with Jaccard >= tau, and `n_settings_run`.
    Symmetric (prereg states.rescue): A -> W if n_settings_hit >= k; P -> W if n_settings_hit < k
    (P supported in too few settings). Rows not re-described are unchanged (`rescue_applied` False).
    """
    out = states.merge(rescue[["sog_id", "accession", "n_settings_hit", "n_settings_run"]],
                       on=["sog_id", "accession"], how="left")
    out["rescue_applied"] = out.n_settings_run.notna()
    hit = out.n_settings_hit.fillna(0)
    a_to_w = out.rescue_applied & (out.state == "A") & (hit >= k)
    p_to_w = out.rescue_applied & (out.state == "P") & (hit < k)
    for m in (a_to_w, p_to_w):
        out.loc[m, "state"] = "W"
        out.loc[m, "reason"] = ("rescue:" + hit[m].astype(int).astype(str) + "/"
                                + out.loc[m, "n_settings_run"].astype(int).astype(str))
    return out
