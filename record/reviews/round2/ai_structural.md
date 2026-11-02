# Round-2 review — AI failure modes x structural biology / co-folding (DESIGN v1)
Read: DESIGN v1 (§12), reviews/round1/*, docs/COFOLD.md. Written: config/prereg.yaml (draft; every `status:
proposed` key needs a PI value before G1), docs/VERIFY_CHECKLIST.md.
**Verdict.** v1 absorbed round 1 honestly; §12 is an accurate ledger with three exceptions. The residual risk has
shifted from "unenforced promises" to "promises enforced on a pipeline that may have no Part B yield".

## 1. Round-1 points mis-implemented (§12 claims them as accepted)
1. **Percolation guard dropped (structural 1).** §3.1 keeps the single-linkage Jaccard-0.3 union-find and only
   *flags* sites > 60 residues. The review asked for overlap coefficient >= 0.5 AND Jaccard >= 0.3 plus a
   component-size cap with a Leiden split. Detecting percolation after the footprint and fingerprint are already
   destroyed is not a guard. Adopt the AND rule + cap (`sites.merge_overlap_coefficient`, `max_site_residues_split`).
2. **Rescue k never calibrated (structural 8).** ">= 2/6 settings" survived verbatim as a frozen constant; the
   review asked for the smallest k meeting a preregistered false-A target, and for rescuing P calls too. k is now
   `tunable` in prereg, grid {1,2,3,4}, 0.05 false-A target; P-rescue proposed.
3. **Interface/TED rules scoped to case members only.** §3.3 sets U for oligomeric interfaces "for case members";
   §5.4 demotes the PAE rule to "TED as proxy". At scale this leaves interface A calls and inter-domain footprints
   in the Part A rates and in candidate generation. Apply both at SOG scale: interface -> U everywhere; footprint
   spanning >= 2 TED domains -> U unless a fetched PAE passes.

## 2. Round-1 points still open (not in v1 at all)
4. **Blinding of co-folding scoring (AI-failure B8).** No primary readout is named among the four, and the
   derived/ancestral labels are visible to whoever scores the models. Name ligand pLDDT primary with the
   benchmark cut, hash state labels, require the per-seed difference in all 3 seeds.
5. **Tuning/validation circularity (B1).** tau is tuned on UniProt ECO / BioLiP2 ortholog sites; the same
   annotations are evidence filter 2. Exclude tune-half annotations from filter 2 for the same SOG, report overlap.
6. **Case provenance and the verifier's authority (AI-failure A4, A6).** No `cases/<id>/manual.yaml`
   (author_agent, source, retrieved_at), no lint against accession/ligand literals in `src/`, and nothing forbids
   a producer marking its own tree node. The control list in §5.3 is currently hand-typed into DESIGN itself.
7. **Smaller gaps**, all now in prereg.yaml: positive-control *gate* + "uninterpretable" label (A10); MSA-method
   flip -> U for case candidates (structural 2); ~10-decoy property-matched panel (structural 10); ASR posterior
   >= 0.8 over the footprint and the shuffled-ancestor control (structural E); homodimer and cofactor-present
   co-folds (structural 7); "unknown" pair class when TED is missing (data 2); PocketMiner; figure sidecars;
   named fallback counters + the `except:` grep test.

## 3. New flaws
8. **The tune/test key is not stable.** §0 splits on `sha1(hog_id) % 2`, but root-HOG ids are reassigned when O1
   becomes O2, so tau frozen at G2 lands in a different split at G3 and the "test half" is no longer held out.
   Split on a stable key (lexicographically smallest member accession) and force the control-list families into
   the test half. Also `sha1(x) % 2` is undefined on a digest object — specify `int(hexdigest, 16) % 2`.
9. **Part B may be arithmetically empty, and nothing checks this before G3.** A gain needs >= 2 outgroup members
   at footprint mean pLDDT >= 90 with frac(<70) <= 0.2, plus >= 3 present species at >= 80 % informative. Tier-1
   eukaryotes average pLDDT 73.7 with 38 % of lining residues below 70. Add a power preflight to the G2 5-family
   slice: count SOGs with >= 2 informative outgroup members under the primary cut, and if the yield is near zero,
   declare Part B loss-only or prokaryote-first *before* seeing candidates, not after.
10. **E < 0.01 is per query, so false positives scale with the query set.** ~3.2 M queries means ~32,000
    significant pairs from the threshold alone. r_c and F_analog are ratios, so the bias partly cancels, but the
    absolute "N pairs with analogous sites" headline does not. Report the expected FP count beside every r_c.
11. **The M1 re-detection ceiling is a scalar applied across a gradient.** Re-detection is itself pLDDT- and
    source-dependent and pLDDT falls with evolutionary distance here, so dividing a distance-dependent rate by a
    constant ceiling leaves exactly the trend M1 measures. Stratify the ceiling by pLDDT band and model source.
12. **Full-atlas scaling is still prose.** Principle 0 promises the same commands on 0.8 B pockets, but exact
    all-vs-all is infeasible (data change 9) and §9 names neither the IVF-PQ + exact-rerank design nor the fixed
    clade-stratified query set. Write it into §9 now, while the commands are being written.
13. **S has an undefined term.** `min(state confidence)` is defined nowhere and parsimony yields no posterior.
    Define it (e.g. rescue-setting agreement x footprint pLDDT percentile) or drop it — it allocates GPU time.

## 4. Keep visible
The four co-folding positive controls reproduce at 0.28-1.22 A and decoys separate only on ligand pLDDT and min
P-L PAE. That benchmark, not intuition, sets the case thresholds: freeze them from those distributions before the
first case runs, and run the poly-Gly ablation on the same four systems first (6 GPU-min).
