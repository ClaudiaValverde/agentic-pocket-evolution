"""Toy tests with known answers for scripts/partB/crossog.py (cross-OG site matching, exploratory D31)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / "src"))
sys.path.insert(0, str(P / "scripts"))
sys.path.insert(0, str(P / "scripts/partB"))
import crossog as C  # noqa: E402


def test_cigar_matches_aln_strings():
    # q: ABCDEFG (start 1), t: XBCYEFG (start 3); alignment q B C - D E vs t B C Z - E
    qaln, taln = "BC-DE", "BCZ-E"
    a = C.map_from_aln(qaln, taln, 2, 3)
    b = C.map_from_cigar("2M1D1I1M", 2, 3)
    assert a[0].tolist() == [2, 3, 5] and a[1].tolist() == [3, 4, 6]
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


def test_row_positions():
    row = np.frombuffer(b"--AB-C", np.uint8)
    assert C.row_positions(row).tolist() == [0, 0, 1, 2, 0, 3]


def _toy():
    # OG x: rep row "-ABCDE" (cols 0..5); OG y: rep row "ABC--DE" (cols 0..6). identity alignment rep_x i -> rep_y i
    src_colpos = C.row_positions(np.frombuffer(b"-ABCDE", np.uint8))
    dst_p2c = np.flatnonzero(np.frombuffer(b"ABC--DE", np.uint8) != C.GAP)
    amap = np.zeros(7, np.int64)
    amap[1:6] = np.arange(1, 6)
    return src_colpos, dst_p2c, amap


def test_project_footprint_known_answer():
    src_colpos, dst_p2c, amap = _toy()
    # footprint cols 0 (gap in rep), 1 (res1), 4 (res4), 5 (res5) -> y cols 0, 5, 6; 1 unmapped
    cols, nu = C.project_footprint(np.array([0, 1, 4, 5]), src_colpos, amap, dst_p2c)
    assert cols.tolist() == [0, 5, 6] and nu == 1
    # unaligned residue counts as unmapped
    amap2 = amap.copy(); amap2[4] = 0
    cols, nu = C.project_footprint(np.array([1, 4, 5]), src_colpos, amap2, dst_p2c)
    assert cols.tolist() == [0, 6] and nu == 1


def test_directional_jaccard_and_link():
    # mapped {0,5,6} + 1 unmapped vs F_y {5,6}: inter 2, union 3+1 -> 0.5
    assert C.directional_jaccard(np.array([0, 5, 6]), 1, np.array([5, 6])) == 0.5
    assert C.directional_jaccard(np.array([1, 2]), 0, np.array([5, 6])) == 0.0
    # same site in two toy OGs links at tau; a site elsewhere does not
    src_colpos, dst_p2c, amap = _toy()
    c, nu = C.project_footprint(np.array([3, 4, 5]), src_colpos, amap, dst_p2c)     # res 3,4,5 -> y cols 2,5,6
    assert C.directional_jaccard(c, nu, np.array([2, 5, 6])) == 1.0
    assert C.directional_jaccard(c, nu, np.array([0, 1])) < C.TAU


def test_shuffle_destroys_correspondence():
    t = np.arange(1, 101)
    s = C.shift_targets(t, 100, 30)
    assert s[0] == 31 and s[-1] == 30 and sorted(s.tolist()) == t.tolist()
    rng = np.random.default_rng(0)
    for _ in range(50):
        o = C.null_offset(rng, 100)
        assert 25 <= o <= 75


def test_bins_merge_down_and_pvalue():
    rng = np.random.default_rng(0)
    null_id = np.concatenate([rng.uniform(0.2, 0.3, 1000), rng.uniform(0.9, 1.0, 10)])   # top bin sparse
    null_d = np.concatenate([rng.uniform(0, 1, 1000), np.zeros(10)])
    qb, nb = C.bin_index(np.array([0.25, 0.95]), null_id, min_n=500)
    # 0.95 bin (10 pairs) merges down until it reaches the populated 0.2-0.3 bin
    assert qb[0] == qb[1]
    p = C.empirical_p(np.array([1.1, -1.0]), qb, null_d, nb)
    assert p[0] == 1 / (1 + 1010) and p[1] == 1.0


def test_bh_known():
    q = C.bh(np.array([0.01, 0.04, 0.03, 0.5]))
    assert np.allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.5], atol=1e-9)


def test_call_b():
    y = pd.DataFrame({"taxid": [1, 2, 3, 3], "state": ["A", "A", "A", "U"]})
    assert C.call_b(3, y)["call_b"]
    assert not C.call_b(2, y)["call_b"]                              # present side needs >= 3 species
    y2 = pd.concat([y, pd.DataFrame({"taxid": [4], "state": ["W"]})])
    assert not C.call_b(5, y2)["call_b"]                             # any P/W on the absent side kills the call
    y3 = pd.DataFrame({"taxid": [1, 1, 2], "state": ["A", "A", "A"]})
    assert not C.call_b(5, y3)["call_b"]                             # A must span >= 3 species


def test_family_labels_and_pairs():
    mem = pd.DataFrame({"og_id": ["X", "X", "Y", "Y", "Z", "Z", "W"],
                        "accession": ["x1", "x2", "y1", "y2", "z1", "z2", "w1"]})
    ted = pd.DataFrame({"accession": ["x1", "x2", "y1", "z1"], "cath_label": ["3.40.50.720"] * 2 + ["3.40.50.720", "1.10.10.10"]})
    pf = pd.DataFrame({"accession": ["w1", "z1", "z2", "x1"], "pfam": ["PF1", "PF1", "PF1", "PF9"], "clan": [None, None, None, None]})
    labs = C.family_labels(mem, ted, pf)
    with_ted = {og for og, (t, _) in labs.items() if t}
    assert labs["X"][0] == {"3.40.50.720"} and labs["Y"][0] == {"3.40.50.720"}   # y1 is 1/2 = 50 % -> label
    assert C.same_family(labs["X"], labs["Y"], True, True).startswith("TED:")
    assert C.same_family(labs["X"], labs["Z"], True, True) is None
    # W has no TED -> compared through Pfam: Z (PF1 in 2/2) matches W
    assert "W" not in with_ted
    assert C.same_family(labs["Z"], labs["W"], "Z" in with_ted, False) == "Pfam:PF1"


def test_cos_dist():
    A = np.array([[1.0, 0.0], [0.0, 2.0]])
    d = C.cos_dist(A, A)
    assert np.allclose(d, [[0, 1], [1, 0]])
