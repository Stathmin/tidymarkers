# 0050 — A genomic interval is whole before the read

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0045](0045-a-position-is-a-whole-number-before-the-cast.md)
  (a position is a whole number before the cast), [ADR 0044](0044-a-count-has-one-definition.md)
  (one `is_whole_number()` predicate backs every count)

## Context

`genome_seq()` is the single accessor every design path reads a region through
(`snp_template()`, `strand_gc_hint()`, the plot mappers, the specificity
windows). It cast both bounds before looking at them:

```r
genome_seq <- function(genome, chrom, start = 1L, end = genome_length(genome, chrom)) {
  start <- as.integer(start)
  end <- as.integer(end)
```

Measured on the branch before the fix, on the two-row fixture of
`test-genome.R`:

| request | in-memory `DNAStringSet` | indexed `marker_genome` |
|---|---|---|
| `start = 2.9, end = 6` | `"CGTAC"` — a span beginning at base 2 | `"CGTAC"` |
| `start = NA, end = 6` | `"ACGTAC"` — read from base 1 | `Error: 'start' or 'end' cannot contain NAs` |
| `start = 6, end = 2` | IRanges: "each range must have an end that is greater or equal to its start minus one" | same base R error |

Two defects in one accessor:

1. A fractional bound was truncated, so the caller received a base it did not
   name — the same wrong-answer class as ADR 0045, one level below the target
   table. A window widened by one base can include or exclude a SNP.
2. The same invalid request behaved differently per backend: `NA` read a
   sequence from an in-memory genome and failed on an indexed one, so a test
   against one source proved nothing about the other. The accessor's own
   contract says the two "behave identically".

## Decision

`check_coord(x, arg)` in `R/marker_validate.R` is the one gate for a single
coordinate: it refuses anything `bad_position()` refuses, plus a vector of
length other than one, and returns the integer. `genome_seq()` calls it on both
bounds before the cast, then refuses `start > end` with its own message naming
both numbers.

An accepted whole coordinate behaves exactly as before; the gate is placed in
the shared accessor, so every design path inherits it and no caller keeps a
second cast.

## Alternatives rejected

- **Guard each caller (`snp_template()`, the plot mappers, …).** That is how
  the family drifted before: some callers validated, the accessor did not.
- **Round the bounds.** Rounding invents a span the caller never named, the
  same reason ADR 0042 rejects rounding a count.
- **Leave the ordering to IRanges.** The base R message names neither argument
  nor value, and it is not a package condition; the accessor knows both bounds.
- **Also refuse an end past the sequence length.** Left to the engine on
  purpose: the length is backend-specific (a `BSgenome` may not be loaded), and
  the engine's out-of-range message already names the region.

## Consequences

- A fractional, negative, `NA`, `Inf`, vector or reversed interval is refused
  with the package message before any sequence is read, on both backends.
- `genome_seq()` now raises its own error where the in-memory path silently
  read from base 1; no accepted call changes its result.
- The helper is the coordinate counterpart of `bad_position()`, so a new
  reader of a single coordinate has one call to make.

## Verification

- `tests/testthat/test-genome.R`, "genome_seq refuses a coordinate it would
  have truncated": for both a `marker_genome` and a `DNAStringSet`,
  `start = 2.9`, `end = 6.9`, `0`, `NA`, `Inf` and `c(1, 2)` are refused with
  "positive whole number", `6..2` with "must not be past", and `2..6` still
  reads `"CGTAC"`.
- `Rscript dev/test.R "^genome$"` — 48 pass, 0 fail, 1 skip.
