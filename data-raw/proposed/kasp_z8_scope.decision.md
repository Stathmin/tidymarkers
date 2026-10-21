# Decision: what a KASP "conversion rate" is allowed to mean

Status: **decided 2026-10-21.** Closes `kasp_z8_scope.proposal.md`. **No `R/`
file changes** — this is a property of the data, and the guard that enforces it
lives in the test that reads the data.

## The decision

1. **A conversion rate is not one number; it is four different measurements.**
   `design` (a primer/probe triple exists for the SNV), `manufacture` (it was
   made and plated), `work` (on real DNA it produced a call), `usable` (the call
   was two-allele and informative on the material). Every KASP row in this family
   carries `conversion_definition` naming **which** of the four it is.
2. **Rates with different definitions are never pooled, averaged or compared.**
   This is the `Z8.3` guard, now enforced: two rows may share a value only if
   they share a definition.
3. **A per-assay-x-sample rate is a different quantity.** `wheat_platform2019`'s
   6.5 % is a fraction of 2,784 assay x sample points, not of assays; its rows
   keep `scope = "assay_x_sample"` and are never read as an assay conversion.
4. **A single-inbred-pair polymorphism rate is a property of the pair.**
   `sharplab2021` 46/50 and `maize_kasp2020` 71/100 are polymorphism rates on one
   parent pair. Two pairs both giving 92 % proves nothing about either.
5. **`fully_decomposed` means the per-SNP failure REASONS are known**, not that
   design and work were measured separately. `makhoul2020_multi_assembly_sanger`
   is the one row carrying it, for the reason column only.
6. **An in-silico yield and a cluster-quality metric are neither success rates**
   (2.42 % of 2,948,985 candidate sites; separation angle, compactness, NTC
   distance). They stay out of every conversion table.
7. **The `design_ambiguous` number is not a design rate.** In
   `wheat_platform2019`, "thirty-nine of the 50 passed assay design criteria
   across the three platforms" (Methods, an intersection) reuses TaqMan's
   per-platform 39/50 (Results). The extract records both and flags the Methods
   one.

## Why nothing is decomposed

None of the six existing `conversion_rates.csv` rows reports both stages on the
same denominator, so no row can be split into design and work with the files we
hold. The three new sources are **additive, not corrective**: new rows, never
pooled with the CAPS/raw-band rows (Liu 2016, soybean, Snap2Caps), which measure a
band on a gel. A pooled rate built on one source is not a pooled rate —
`maize_kasp2020` is the only row with a work numerator on the design denominator.

## Contradiction named

Two sources print **92 %** and mean different things (`sharplab2021`:
polymorphism between inbreds 1212 × B73; `maize_kasp2020`: calls produced on the
LM11 × CML25 parent pair). Any table putting them in one column is wrong, and the
guard is written so it cannot.

## Resulting `R/` change

None, by design. A conversion rate is not a design parameter; enforcing its
definition in `R/` would imply a model the data cannot support.

## Test that asserts it

`tests/testthat/test-z8-outcomes.R`: every KASP row of `conversion_rates.csv`
carries a non-NA `conversion_definition`; no two rows share a rate value while
their definitions differ; any `assay_x_sample` row states its scope; the two 92 %
rows are labelled differently.
