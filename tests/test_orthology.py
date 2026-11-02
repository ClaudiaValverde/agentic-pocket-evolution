"""Unit tests for pocketevo.orthology on a toy OrthoXML with known ortholog / paralog answers.

Run: .venv/bin/python -m pytest tests/test_orthology.py -v

Toy gene tree (species tree ((HUMAN,MOUSE)Mammalia,FLY)Bilateria; YEAST outside, own root group):
  root S@Bilateria
    D@Mammalia  (duplication in the mammal lineage: copies A and B)
      S@Mammalia  hA, mA          (copy A)
      S@Mammalia  hB, mB          (copy B)
    fly1
  root S@Mammalia  hC, mC (another family)
  yeast1 : singleton (no group)
Known answers: hA-mA ortholog@Mammalia; hA-hB paralog@Mammalia; hA-mB paralog; hA-fly1 ortholog@Bilateria;
hA-hC unrelated; yeast1 unrelated to all; unknown ids -> unknown.
"""
from __future__ import annotations

import gzip
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pocketevo.orthology import HogIndex, parse_orthoxml, to_frames, write_tables  # noqa: E402

TOY = """<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" origin="toy" version="0.5">
  <species name="HUMAN" taxonId="4"><database name="t" version="1"><genes>
    <gene id="1" protId="hA"/><gene id="2" protId="hB"/><gene id="3" protId="hC"/></genes></database></species>
  <species name="MOUSE" taxonId="5"><database name="t" version="1"><genes>
    <gene id="4" protId="mA"/><gene id="5" protId="mB"/><gene id="6" protId="mC"/></genes></database></species>
  <species name="FLY" taxonId="6"><database name="t" version="1"><genes>
    <gene id="7" protId="fly1"/></genes></database></species>
  <species name="YEAST" taxonId="7"><database name="t" version="1"><genes>
    <gene id="8" protId="yeast1"/></genes></database></species>
  <taxonomy>
    <taxon id="1" name="root"><taxon id="2" name="Bilateria"><taxon id="3" name="Mammalia">
      <taxon id="4" name="HUMAN"/><taxon id="5" name="MOUSE"/></taxon><taxon id="6" name="FLY"/></taxon>
      <taxon id="7" name="YEAST"/></taxon>
  </taxonomy>
  <groups>
    <orthologGroup id="G1" taxonId="2"><property name="TaxRange" value="Bilateria"/>
      <paralogGroup id="G1.1" taxonId="3">
        <orthologGroup id="G1.1a" taxonId="3"><property name="TaxRange" value="Mammalia"/>
          <geneRef id="1"/><geneRef id="4"/></orthologGroup>
        <orthologGroup id="G1.1b" taxonId="3"><property name="TaxRange" value="Mammalia"/>
          <geneRef id="2"/><geneRef id="5"/></orthologGroup>
      </paralogGroup>
      <geneRef id="7"/>
    </orthologGroup>
    <orthologGroup id="G2" taxonId="3"><property name="TaxRange" value="Mammalia"/>
      <geneRef id="3"/><geneRef id="6"/></orthologGroup>
  </groups>
</orthoXML>
"""
LEVEL_TAXID = {"Bilateria": 33213, "Mammalia": 40674}


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    d = tmp_path_factory.mktemp("toy")
    p = d / "toy.orthoxml.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(TOY)
    return str(p)


@pytest.fixture(scope="module")
def idx(toy):
    return HogIndex.from_orthoxml(toy, LEVEL_TAXID)


@pytest.mark.parametrize("a,b,rel,level", [
    ("hA", "mA", "ortholog", "Mammalia"),
    ("hB", "mB", "ortholog", "Mammalia"),
    ("hA", "hB", "paralog", "Mammalia"),
    ("hA", "mB", "paralog", "Mammalia"),
    ("hA", "fly1", "ortholog", "Bilateria"),
    ("mB", "fly1", "ortholog", "Bilateria"),
    ("hC", "mC", "ortholog", "Mammalia"),
    ("hA", "hC", "unrelated", None),
    ("hA", "yeast1", "unrelated", None),
    ("hA", "nope", "unknown", None),
    ("hA", "hA", "same", None),
])
def test_relations(idx, a, b, rel, level):
    r = idx.relation(a, b)
    assert r["relation"] == rel
    assert r["level"] == level
    assert idx.relation(b, a)["relation"] == rel  # symmetric


def test_levels_and_taxids(idx):
    assert idx.relation("hA", "fly1")["level_taxid"] == 33213
    assert idx.relation("hA", "hB")["mrca_hog"] == "G1.1"
    hs = idx.hogs_of("hA")
    assert [h[1] for h in hs] == ["G1", "G1.1", "G1.1a"]
    assert [h[2] for h in hs] == ["S", "D", "S"]
    assert idx.hogs_of("yeast1") == []


def test_tables_roundtrip(toy, tmp_path):
    px = parse_orthoxml(toy)
    h, p = to_frames(px, LEVEL_TAXID)
    assert len(h) == 5                       # G1, G1.1, G1.1a, G1.1b, G2
    assert sorted(h[h.parent == -1].hog_name) == ["G1", "G2"]
    assert int(h.loc[h.hog_name == "G1", "n_genes"].iloc[0]) == 5
    assert set(p[p.root == -1].protein) == {"yeast1"}
    write_tables(px, str(tmp_path), LEVEL_TAXID)
    j = HogIndex.load(str(tmp_path))
    assert j.relation("hA", "mB")["relation"] == "paralog"
    assert j.relation("hA", "mA")["level_taxid"] == 40674


def test_alias(toy):
    h, p = to_frames(parse_orthoxml(toy), LEVEL_TAXID)
    j = HogIndex(h, p, alias={"hA_redundant": "hA"})
    assert j.relation("hA_redundant", "mA")["relation"] == "ortholog"


def test_matches_ancestra_pairs(toy):
    """Cross-check against Ancestra's own flattening (skipped if ancestra is not importable)."""
    orx = pytest.importorskip("ancestra.orthoxml")
    ox = orx.read_orthoxml(toy)
    orth, para = orx.hogs_to_pairs(ox, with_paralogs=True)
    idx = HogIndex.from_orthoxml(toy)
    for a, b in orth:
        assert idx.relation(a, b)["relation"] == "ortholog"
    for a, b in para:
        assert idx.relation(a, b)["relation"] == "paralog"
