# Decision: the RPA 3'-anchor mismatch scheme

Status: **decided 2026-10-21.** Closes `rpa_anchor_rule.proposal.md`. The
proposal's own open questions are settled here, against the landed files — no new
paper is needed, and the 279 position-class cells already answer the S1-vs-S2
question empirically.

## The model, stated once

The anchor is the **4 nt at the primer's 3' end** (offset 0 = terminal), as the
source defines it ("pre-ante-penultimate to 3'-terminal"); that width is asserted
three ways in the proposal and is not revisited here.

`w(position, mismatch)` is a penalty in log-odds (Firth coefficient on reaction
success), negative = worse:

| offset | weight | evidence |
|---|---|---|
| 0 (terminal) | per-base, from the 12 terminal cells | measured |
| 1, 2, 3 | one class penalty per offset | measured |

**Ranking, which is what a design uses.** The package does not need the absolute
log-odds; it needs the order. Across the 279 landed position cells (99 at offset
1, 108 at 2, 72 at 3) the non-terminal offsets are **not** ordered:

| offset | n | mean | min | max |
|---|---|---|---|---|
| 1 (penultimate) | 99 | −3.82 | −7.36 | −1.80 |
| 2 (ante-penultimate) | 108 | −3.74 | −7.41 | −1.94 |
| 3 (pre-ante-penultimate) | 72 | −3.74 | −6.60 | −1.96 |

So the proposal's "3 is worse than 2, which is worse than 1" does **not** hold on
the paper's own cells: offsets 2 and 3 are indistinguishable on the mean, and
offset 1 has the best mean but the widest spread. The rank asserted by the design
rule is therefore **terminal ≫ non-terminal, with no order among the non-terminal
positions.**

**Terminal weights**: used as printed, including the two unstarred cells (A–C,
T–G = −0.66). Not weak evidence — RPA is *tolerant of exactly these two*, which
the source's own probability (0.984) states. Encoding them as "unknown" would
invent a rung the data does not show; encoding them as 0 would overstate the
tolerance. They are used as printed and labelled unstarred.

**Pair terms**: when two adjacent mismatches fall inside the anchor, the
**measured class coefficient** from `rpa_mismatch_position.csv` replaces the sum
of two singles. Summing is closer (C−C + 1n gives −6.14 against −6.77 measured),
but the class cell is the measured value.

## Contradictions named, not silently resolved

1. **The design rule's own ordering does not survive its data** (table above).
   Both are kept: the rule is quoted, the cells are the ranking.
2. **S1 and S2 print slightly different values for the same class** because they
   come from separate model fits. The decision uses **S1** (the per-base sheet,
   which is where the terminal cells live, so one fit describes the whole anchor)
   and says so; S2's `1n`/`2n`/`3n` are kept as a cross-check.
3. **The field data contradict any per-mismatch error rate**: AS-RPA missed 8 of
   54 rifampicin-resistant isolates, and sickle-cell AS-RPA reached 100 % for
   both alleles on a design this table scores at −1.71. One genotype cannot
   calibrate a weight, and a coefficient is not an error rate.

## What the scheme must NOT be used for

- Predicting a call rate (measured above: it does not follow from the log-odds).
- Ranking alleles at a transition.
- Anything outside the 4-nt anchor, or with more than two mismatches.

## Lenviency

The weights **rank**; they do not reject. `rpa_variant_gate()` keeps the
geometry it already had (the variant within the 3' window, or under the probe's
span), and the anchor score is reported beside it. Turning the score into a gate
would fail designs on a quantity the source's own field data show is not
predictive.

## Resulting `R/` change

* `R/marker_rpa.R`: `RPA_ANCHOR_NT <- 4L`, `rpa_anchor_weights()` (terminal
  cells in the extract's row order, one class penalty per non-terminal offset,
  `evidence` per entry), `rpa_anchor_weight()`, and
  `rpa_anchor_score(primer, template)` summing the weights over the mismatches
  inside the anchor. The existing gate is unchanged.

The extract prints each ordered mismatch **template-first** (`C-T` is a primer T
reading a template C); the scorer follows that. Cells are resolved positionally,
in the extract's own order, because on this R a character subscript containing a
hyphen does not resolve reliably - a defect found here, not a design choice, and
the test pins the order so it cannot drift.

## Test that asserts it

`tests/testthat/test-rpa-anchor.R`: the weight table equals the extracts (the 12
terminal cells of `rpa_mismatch_terminal.csv`, including the two unstarred ones);
`rpa_anchor_score()` is 0 for a perfect anchor and equals the printed terminal
weight for a single 3'-terminal mismatch; the 279 cells do NOT order the
non-terminal offsets (mean at 2 and 3 within 0.1 of each other) — which is the
decision's stated ranking; and the gate's behaviour is unchanged.
