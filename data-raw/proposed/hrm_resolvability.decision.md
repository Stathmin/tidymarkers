# Decision: HRM resolvability

Status: **decided 2026-10-21.** Closes `hrm_resolvability.proposal.md`.

## The model, stated once

HRM's outcome is a **difference of melting peaks between the two homozygotes**,
not a yes/no property of a SNP, and not a per-primer score.

1. **The pass bar is `dTm >= 0.5 °C` on the class mean.** This is the only
   numeric bar any fetched source states as its own rule ("a dTm value >0.5 C
   guaranteed accurate genotyping", JSBBS 2018) — and the source itself says
   it held "using the Roche LightCycler 96 system" and that some class III/IV
   markers still genotyped at 0.2–0.3 °C. So: a **platform-scoped default for a
   class**, never a physical constant.
2. **A class mean below the bar is `marginal`, not a failure**, when a member of
   that class still cleared it. The class can work; it is not dependable. A class
   whose *only* measured member is below the bar is `below_threshold`.
3. **A bimorphism with no landed measurement is `unknown`** — never `pass`, never
   `fail`. No fetched source measures it, so nothing is claimed.
4. **A multi-locus panel is not answered per SNP.** With 3 C>T loci, HRM alone
   left 14 of 36 control genotypes in 7 indistinguishable pairs. The outcome is a
   *partition* of genotypes, which a per-SNP flag cannot express.
5. **`|ΔG3bp|` is a filter, not a guarantee.** The paper's own numbers disagree
   with its prose (prose mean 0.63 vs recomputed 0.69; prose t 7.6 vs 8.15), so
   the free-energy value is recorded and never thresholded in the gate. 1.43
   kcal/mol stays Exp. 3's choice, not a general rule.

## The measurements this rests on

`hrm_jsbbs2018_dtm_by_bimorphism.csv` (38 + 18 soybean/rice markers, one
instrument), per class min / mean dTm:

| class | n | min | mean |
|---|---|---|---|
| A/C | 4 | 0.514 | 0.589 |
| A/G | 18 | 0.245 | 0.649 |
| C/G | 1 | 0.198 | 0.198 |
| C/T | 28 | 0.392 | 0.692 |
| G/T | 5 | 0.367 | 0.519 |

The code's lookup is asserted against that extract by
`tests/testthat/test-hrm-jsbbs2018.R`, so the constant cannot drift from the file.

## Contradictions named, not resolved silently

1. **The source contradicts itself on where 0.5 °C applies**: it states the bar
   for LightCycler 96 and then reports class III/IV markers called correctly at
   0.2–0.3 °C. Both clauses are kept; the threshold is applied to the class mean
   and the sub-threshold members stay visible as `dtm_min`.
2. **The free-energy filter does not order the classes**: G/T has the second
   lowest mean dTm but is not second on mean `|ΔG3bp|` (A/G is). A free-energy
   rank is therefore not used as a resolvability rank.
3. **The widely cited transversion table (Liew 2004)** could not be fetched
   (closed access), so **nothing** in this decision cites or encodes it — not even
   its A::T↔G::C 0.8–1.4 °C figure.

## Lenviency

No HRM design is **failed** by this decision. The criterion reports `unknown` for
`marginal`/`below_threshold`/unmeasured, and carries its reason in `detail`. We
already constrain HRM by amplicon length, Tm balance and geometry; adding a hard
resolvability gate would drop designs whose class the evidence cannot condemn.

## Resulting `R/` change

* `R/marker_hrm.R`: `HRM_DTM_RESOLVABLE` (0.5), `HRM_BIMORPHISM_DTM` (the table
  above), `hrm_resolvability(ref, alt)`, and `hrm_resolvability_row(marker)`.
* `R/marker_checks.R`: `marker_design_checks()` appends that row for
  `scheme == "hrm"`, so an HRM evaluation always says whether its class has a
  curve to read.

## Test that asserts it

`tests/testthat/test-hrm-jsbbs2018.R` — the lookup equals the extract row by row;
the statuses for a resolvable class (C/T, A/C), a resolvable-but-not-assured one
(G/T, min below the bar), a below-threshold one (C/G), and `unknown` for an
unmeasured or non-substitution pair; the 3-locus partition counts.
`tests/testthat/test-hrm.R` — the evaluation carries `family_hrm_resolvable` with
the measured detail.
