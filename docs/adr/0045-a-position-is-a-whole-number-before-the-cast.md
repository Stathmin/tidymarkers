# 0045 — A position is a whole number before the cast

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)
  (a position outside the positive integer range is refused, not coerced),
  [ADR 0042](0042-a-threshold-is-a-count-not-a-number-to-truncate.md) (a
  threshold is a count, not a number to truncate)

## Context

Three normalisers cast a `pos` column before the check that claims to guard it:

| path | cast | check that follows |
|---|---|---|
| `normalise_targets()` | `as.integer(x$pos)` | `check_target_loci()`: "positive integer `pos`" |
| `target_rows()` (grammar form) | `as.integer(pos)` | `check_target_loci()` |
| `normalise_genotypes()` | `as.integer(g$pos)` | only `is.na()` |
| `normalise_haplotypes()` | `as.integer(.data$pos)` | `loci_bad`: `is.na` / `< 1` |

`check_target_loci()` documents that "a POS that is not a positive integer
would reach the designer as `NA`/`0`", but `as.integer()` had already truncated
the fraction, so the branch could only ever see `NA` and `0`. Measured on the
branch:

```r
panel <- tibble::tibble(
  line = paste0("L", 1:8), pos = 250.9, ref = "A", alt = "G",
  genotype = c(2, 2, 2, 2, 0, 0, 0, 0)
)
marker_informativeness(tg_at_250, genotypes = panel)
#> n_lines 8, alt_freq 0.5, status "ok"
```

The row names a locus at 250.9 and was scored as the target at 250: a panel
record for a position that does not exist decided the informativeness of a
different variant. `marker_targets(chrom = "Chr1", pos = 250.9, ...)` likewise
produced a target labelled `Chr1:250` with no complaint.

## Decision

A position is checked on the value the caller gave, before any cast, with one
vectorised predicate `bad_position(x, min = 1L)` in `R/marker_validate.R`:

```r
iv <- suppressWarnings(as.integer(x))
is.na(x) | is.na(iv) | iv < min | (is.numeric(x) & (!is.finite(x) | x != round(x)))
```

A character `"250"` that casts exactly is still accepted (CSV input), as
before. `check_target_pos()` keeps the existing message and holds the target
path's check in one place; the POS branch leaves `check_target_loci()`, which
now judges only the chromosome and the alleles. The genotype and haplotype
normalisers call the same predicate on their raw column.

## Alternatives rejected

- **Round the position instead of refusing it.** `250.9` would become `251`,
  still not the locus the caller named; ADR 0018 already refuses a position it
  cannot use rather than coercing it.
- **Keep the check after the cast and compare `as.integer()` with the input.**
  That is the same predicate written twice, and the comparison is easy to drop
  in a later edit - which is how this defect appeared.
- **Refuse in `read_variants_vcf()` only.** The panel path does not come from a
  VCF: a hand-built tibble reaches `marker_informativeness()` directly.

## Consequences

- A fractional, `Inf` or otherwise uncastable `pos` is refused in targets,
  grammar targets, genotype panels and haplotype tables, with the row and the
  value in the message.
- Every previously accepted whole position behaves as before; no statistic or
  verdict changes for an accepted input.
- `check_target_loci()` no longer checks POS; `check_target_pos()` does, and
  runs earlier.

## Verification

- `tests/testthat/test-plan.R`: `pos = 250.9` and `pos = Inf` in the grammar
  form and `250.9` in a data table abort with "positive integer `pos` ... row 1
  is 250.9" / "is Inf"; `0`, `NA`, `-1` keep their previous messages.
- `tests/testthat/test-informativeness.R`: a panel at 250.9 against a target at
  250 aborts with "positive whole `pos` ... row 1 is 250.9" (was `status "ok"`).
- `tests/testthat/test-haplotypes.R`: a haplotype row at 250.9 aborts with
  "positive integer".
- `Rscript dev/test.R "plan|informativeness|haplotype"` — 314 pass, 0 fail.

## Amendment (2026-11-28): the predicate is type-safe

`bad_position()` combined its fraction test with `&` on a `is.numeric()` guard:

```r
is.na(x) | is.na(iv) | iv < min |
  (is.numeric(x) & (!is.finite(x) | x != round(x)))
```

`&` evaluates both sides, so a character column still reached `round()`. A
target row with `pos = "abc"` therefore raised base R's "non-numeric argument to
mathematical function" instead of the message this ADR promises — for every
caller, not only the new one.

The fraction test is now built only when the column is numeric:

```r
frac <- if (is.numeric(x)) !is.finite(x) | x != round(x) else rep(FALSE, length(x))
is.na(x) | is.na(iv) | iv < min | frac
```

`tests/testthat/test-plan.R` adds `pos = "abc"` and expects "positive integer
`pos` ... row 1 is abc"; `Rscript dev/test.R "^plan$"` — 189 pass, 0 fail. No
accepted position changes.

## Amendment (2026-11-28): the annotation reader applies the same rule

`marker_features()` — the exported reader behind [marker_genome()] windows —
built its table with `start = as.integer(start_v)` and `end = as.integer(end_v)`
on both the NCBI feature-table and the GFF3 branch. A file coordinate of
`2801.5` was read as `2801`, `abc` became `NA` with a base R warning, and `-1`
passed as a start. A feature window could therefore move by a base or lose its
bounds without the reader saying so.

The reader now runs `bad_position(start_v) | bad_position(end_v)` on the raw
columns before the cast and aborts with
`` `path` has a start/end that is not a positive whole number (row N). ``
The accepted forms are unchanged: `"2801"` casts exactly and a whole number
column passes.

Verification: `tests/testthat/test-annotate.R` refuses `2801.5`, `-1` and `abc`
and reads `2801`/`2900` unchanged; `Rscript dev/test.R annotate` — 10 pass,
0 fail, 5 skip.

## Amendment (2026-11-28): the mask helper applies the same rule

`mask_positions()` cast its `positions` vector with `as.integer()` before the
range filter, so `2.9` masked base 2 and `"abc"` raised base R's
"NAs introduced by coercion" warning instead of a message. Its documented
contract keeps positions outside the sequence ignored, so this is not
`check_scan_radius()` (negative and zero stay silently ignored); the rule is
only "whole in its spelling", which `bad_position(positions,
min = -.Machine$integer.max)` answers.

The helper now refuses `2.9`, `Inf`, `NA` and `"abc"` with
`` `positions` must be whole numbers. `` and keeps `c(0L, 99L)` ignored as
before. Verification: `tests/testthat/test-masking.R` refuses the four values
and keeps the range case; `Rscript dev/test.R "^masking$|^io$"` — 243 pass,
0 fail, 3 skip.
