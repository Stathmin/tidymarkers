# Decision: what a TaqMan/5'-nuclease design may claim about mismatch

Status: **decided 2026-10-22.** Closes `taqman_status.proposal.md`. **No `R/`
model change** — the current constants stay, but their status changes from
"calibrated" to "engineering default", and the claims attached to them change.

## The decision

1. **The mismatch model is an engineering default. It is not calibrated.** No
   open source measures a matched-vs-mismatched TaqMan probe readout.
   `TAQMAN_MISMATCH_DG` (1.5 kcal/mol) and the probe-centrality rule keep their
   values, because the package needs *a* policy to place a probe base, and moving
   a default without a measurement only exchanges one guess for another.
2. **Neither default may be cited as literature-supported.** `LEGACY.md` /
   `R/marker_taqman.R` already say so (L49-63); the decision binds future prose
   and tests to it. `wheat2019` and `citrus2020` support TaqMan-MGB *working*
   (384/384 calls in `citrus2020`), not the mismatch threshold.
3. **Centrality is not a rule; it is an unmeasured preference.** In
   `citrus2020` the allele base sits at centrality 0.73-0.95 in eight assays, and
   five of eight FAM probes are 13-18 nt — below the 18:25 default range — while
   still calling all 48 cultivars. The design keeps centre placement for a *new*
   probe (no better information), but a landed assay is never graded "wrong" for
   an off-centre allele base.
4. **A gap in the probe footprint is a real, observed failure mode and stays
   modelled.** `apol1_2019`: a 6-bp deletion under the probe leaves 1.4 % of
   1,110 samples uncalled and 5.9 % mis-typed. This is the one TaqMan failure
   mode the package has evidence for; it is not a mismatch claim and is not
   withdrawn.

## Why the search is the finding

`ABSTRACT:"delta Ct" AND ABSTRACT:probe AND ABSTRACT:mismatch` returns **0**
hits in Europe PMC. The papers that state the trade-off (`yao2006` five-mismatch
tolerance, `stadhouders2010` position dependence, `lna2004` MGB/LNA) are
publisher-closed; their metadata was checked and is recorded in the proposal.
Under the evidence rule an identifier is provenance, not evidence, so none of
them contributes a number.

## Contradiction named

Two claims were both held at once: "a single mismatch under 1.5 kcal/mol is not
discriminated by the exonuclease step" and "the probe is centred so the mismatch
is maximally destabilising". The first has no supporting measurement; the second
is contradicted by every assay in `citrus2020`. They are kept as defaults,
downgraded in status, and never quoted as findings.

## Resulting `R/` change

None to the values. Comment and documentation status only (already in
`R/marker_taqman.R` L49-63); this file records why.

## Test that asserts it

`tests/testthat/test-taqman-status.R`: the probe-centre preference is a default
(an off-centre landed design is not marked "failed" by any check row); the
mismatch constant is documented as an engineering default (source text guard on
`R/marker_taqman.R`); `apol1_2019`'s two rates stay attributable to a footprint
gap, not a mismatch.

## What would overturn it

Any open source reporting (i) ΔCt or endpoint ratio between matched and
single-mismatched probe over a set of positions, (ii) the same probe with and
without MGB, or (iii) the same comparison at two probe lengths. A licensed copy
is acceptable; a paywall is never bypassed.
