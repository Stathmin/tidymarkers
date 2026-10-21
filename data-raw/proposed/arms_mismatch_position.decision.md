# Decision: the deliberate-mismatch position

Status: **decided 2026-10-21.** Applies to the shared AS-primer engine
(`design_as_primers()`) and therefore to every scheme built on it: KASP, ARMS /
tetra-primer, AS-LAMP, AS-RPA, CAST-PCR, ASQ, duplex. It does **not** change any
existing verdict at `stab_ats = NA` (a 3'-terminal-mismatch-only design), which
is what the `f494a50` discipline protects.

## What was wrong before this decision

`design_as_primers()` defaults to `stab_ats = c(NA, 2L, 3L)` and both positions
were ranked in one pool by discrimination alone. That is a silent average of the
contradiction: whichever position wins on one locus is returned, so the package
had no position policy at all — it had a coin flip with the literature's name on
it.

## The model, stated once

Position is counted **negative from the primer's 3' end** (last base = −1). Two
layers, never averaged:

1. **The terminal base (−1) carries the discrimination.** A lone terminal
   mismatch is `unknown` at the specificity policy, never a pass. The extra
   mismatch is a *rescuer* for a weak terminal class, not a substitute for it.
2. **The extra mismatch goes at −3 by default.** −2 is the **fallback, for a
   weak terminal class**: Ye 2001 arranges strong/medium/weak terminal classes and
   the extra mismatch is a rescue for the *weak* ones, so a weak terminal base
   (the package's own index: T:G / T:C / T:T, i.e. Kwok's `minimal` class, plus
   the G:T / A:C wobble `specificity_policy()` already flags) gets −2. Every
   other class, including `context_dependent`, keeps the measured −3: the
   rule of thumb is not measured, and we do not spend it where it is not needed.

## Rationale, and its limits

Two controlled scans, at two chemistries, both peak at −3:

* Liu 2012 — 1686 AS primers, AS-PCR: 25.1 / 32.9 / 29.1 % specific at the 2nd /
  3rd / 4th base from the 3' end (139/554, 189/575, 162/557 = 490/1686). U-shaped,
  peaking at the 3rd base.
* Lefever 2019 — qPCR, 6 SNPs × positions 1–5 (0 = terminal): optimum at
  position 3, i.e. −3.

−2 has **no controlled measurement**: Ye 2001 states the position and reports no
per-primer outcome we hold; Zhu 2021's −2 primer is a `primer_scan` (position
read from Table S1, ranking only in a figure).

Their **magnitudes are not comparable** (percent polymorphic vs ΔCq, different
chemistries), so nothing is pooled: the *position* is shared evidence, the *size*
is not. No ΔΔG number is attached to the −3 default, and it is carried as
directional, not calibrated (`calibrated = FALSE`).

## Contradictions named, not silently resolved

1. **Ye 2001 (−2) vs Liu 2012 + Lefever 2019 (−3).** Both rows stay in
   `z7_priors.csv`; Ye's is `rule of thumb`, `measured = FALSE`. The default is a
   policy choice, not a finding that Ye was wrong.
2. **Zhu 2021's "fourth from the 3' end" primer also carries −1** (wxt F4:
   mismatch at −1 and −4), so it is not evidence about −4 alone; only its
   wxc F4 (−2, clean) is used, and only as `primer_scan`.
3. **Liu's own curve is U-shaped** (2nd < 4th < 3rd). "The 3rd base is best" is a
   peak between two inferior neighbours — reading it as "later is better" is the
   error this file exists to prevent.

## Lenviency

Neither scan covers KASP's fluorescent readout (Kwok's data is Taq kinetics), so
a design at −2 outside the fallback case, or at −4, is not **rejected**: it is
kept, reported with its position, and the policy marks it `marginal`. Only the
ranking changes. Hard-failing it would delete assays the evidence cannot
condemn — and we already carry many constraints.

## Resulting `R/` change

* `R/marker_arms.R`: `ARMS_EXTRA_MISMATCH_POSITION <- -3L`,
  `ARMS_EXTRA_MISMATCH_FALLBACK <- -2L`, and
  `arms_extra_mismatch_position(terminal_class)` — −3 by default, −2 only for
  `minimal` (the weak terminal class), `NA` for a class with no measurement.
* `R/marker_design.R`: `design_as_primers()` ranks the chosen position first,
  the fallback second (only when that configuration's terminal base is weak),
  no-mismatch third, and any other position last, and records the rule and the
  weak-terminal flag on each row. The pool is unchanged;
  only the order is.

Both are directional: no primer is added or removed by the policy, so a design
that existed before still exists.

## Test that asserts it

`tests/testthat/test-arms-mismatch-position.R`: the constants; the −3 / −2 /
`NA` returns of `arms_extra_mismatch_position()`; `design_as_primers()` puts −3
ahead of −2 at equal discrimination and marks the row; `z7_priors.csv` keeps the
−3 row `primary literature` and the −2 row `rule of thumb`, unmeasured; and the
landed ledger (`arms_mismatch_position.csv`) still has two measurements at −3 and
none at −2.
