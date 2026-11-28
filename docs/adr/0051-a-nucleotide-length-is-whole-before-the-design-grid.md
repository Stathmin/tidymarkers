# 0051 — A nucleotide length is whole before the design grid

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0045](0045-a-position-is-a-whole-number-before-the-cast.md)
  (a position is a whole number before the cast), [ADR 0044](0044-a-count-has-one-definition.md)
  (one `is_whole_number()` predicate backs every count), [ADR 0050](0050-a-genomic-interval-is-whole-before-the-read.md)
  (the same defect one level below)

## Context

Two exported design entries took the length of an oligonucleotide and cast it
with `as.integer()` before looking at it:

| entry | cast |
|---|---|
| `design_crispr_marker(spacer_lengths = )` | `as.integer(spacer_lengths)` before the candidate grid |
| `specificity_policy(length = )` | `as.integer(length)` before `mismatch_weights()` |

Measured on the branch before the fix, on the `test-crispr.R` fixture:

```r
design_crispr_marker(g, snp, cas = "cas12a", spacer_lengths = c(20.9, 21.5))
#> a 20 nt guide
specificity_policy(length = 20.9)$weight
#> 20 weights
specificity_policy(length = NA)$weight
#> Error in if (!n_nt) ... : missing value where TRUE/FALSE needed
```

A spacer length is not decoration: it fixes where the guide's seed and the cut
site sit, so a truncated `20.9` is a guide the caller did not ask for — the same
wrong-answer class as ADR 0045, one level below the target table. A missing
length reached a base R `if` instead of a package message.

The first version of the new gate also exposed a latent defect in the shared
predicate: `bad_position()` combined its fraction test with `&` on a
`is.numeric()` guard, but `&` does not short-circuit, so a character column
reached `round()` and raised "non-numeric argument to mathematical function"
instead of the message ADR 0045 promises. That fix is recorded as an amendment
to ADR 0045; this decision only names the length gate.

## Decision

`check_nucleotide_length(x, arg, multi = FALSE)` in `R/marker_validate.R`
admits one or more positive whole numbers of nucleotides through
`bad_position()`, and returns the integer vector. A series is accepted only
when the caller's entry takes a series:

- `design_crispr_marker(spacer_lengths = )` — `multi = TRUE`, the candidate
  series;
- `specificity_policy(length = )` — one target-site length.

An accepted whole length behaves exactly as before; no grid, guide or weight
ladder changes for it.

## Alternatives rejected

- **Round the length.** A rounded guide is still a guide the caller never
  named, the reason ADR 0042 rejects rounding a count.
- **Check inside `crispr_guide_candidates()` / `mismatch_weights()`.** Both are
  internal and are also called from already-validated paths; the boundary is
  the exported argument, and one gate at the entry is the rule that held for
  the radius and the worker count.
- **Accept a vector everywhere.** `specificity_policy(length = c(20, 21))`
  produced "the condition has length > 1" from `mismatch_weights()`; the
  scalar is the contract, so the gate refuses the series rather than letting a
  downstream `if` do it.

## Consequences

- A fractional, zero, negative, `NA`, `Inf`, non-numeric or (where a scalar is
  required) vector length is refused with the package message before any grid
  or weight ladder exists.
- Accepted whole lengths are unchanged; no guide, weight or verdict moves.
- The type-safe `bad_position()` means the same refusal now also works on a
  character coordinate column, which is the ADR 0045 amendment.

## Verification

- `tests/testthat/test-crispr.R`, "design_crispr_marker refuses a spacer length
  it would have truncated": `c(20.9, 21.5)`, `20.5`, `0`, `-1`, `NA`, `Inf` and
  `c(20, NA)` are refused with "positive whole numbers of nucleotides";
  `20:23` still designs a guide of a named length.
- `tests/testthat/test-specificity.R`, "specificity_policy refuses a target
  length it would have truncated": `20.9`, `0`, `-1`, `NA`, `Inf`, `c(20, 21)`
  and `"x"` are refused; `20L` and `"20"` still give a 20-position ladder.
- `Rscript dev/test.R "^crispr$"` — 152 pass, 0 fail; `"^specificity$"` — 324
  pass, 0 fail, 2 skip.
