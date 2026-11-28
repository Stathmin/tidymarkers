# ADR 0071 — A region bound is a coordinate before the cast

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0050
  (a coordinate is checked before the read), ADR 0068 (a position is whole in
  its own spelling), ADR 0069 (a GenBank location is read whole or refused)

## Context

`normalise_region()` accepted a data frame or list with `chrom`, `start` and
`end`, cast the bounds with `as.integer()` and validated the result:

```r
out <- tibble::tibble(
  chrom = as.character(x$chrom),
  start = as.integer(x$start),
  end   = as.integer(x$end)
)
if (anyNA(c(out$chrom, out$start, out$end)) || !nzchar(out$chrom) ||
    out$start < 1L || out$end < out$start) abort_msg(...)
```

A fractional bound therefore reached the comparison already shifted:

| input | before the fix |
|---|---|
| `list(chrom = "Chr1", start = 10.5, end = 20L)` | region `10..20` |
| `list(chrom = "Chr1", start = 10L, end = 20.9)` | region `10..20` |
| `list(chrom = "Chr1", start = 10.9, end = 10.1)` | inverted window becomes `10..10` |

The shifted region decides which variants fall inside the block, so the boundary
is a statistic, not a cosmetic difference. The `GRanges` branch is unaffected:
its start/end are integers by construction. `normalise_region()` is reached by
`marker_plan(region = )` and `read_genbank(region = )`.

## Decision

- The data-frame/list branch checks `bad_position(x$start)` and
  `bad_position(x$end)` **before** `as.integer()`, and refuses a fractional,
  `NA`, `Inf`, non-numeric or non-positive bound with a message naming both
  values.
- The `GRanges` branch keeps its current cast; the existing chromosome and
  `start <= end` check stays as the last rule for both branches.

## Alternatives rejected

- **Check the cast value.** `as.integer(10.5)` is `10`, which is a plausible
  coordinate; only the spelling can tell that the caller named a different one.
- **Round the bound.** A rounded region includes a base the caller did not name.
- **Validate only after the existing `start <= end` check.** The check runs on
  the truncated values and cannot recover the inverted `10.9..10.1` window.

## Consequences

- `marker_plan(region = )` and `read_genbank(region = )` read the interval the
  caller wrote, or refuse it.
- Integer and `GRanges` regions behave exactly as before; no accepted input
  changes the variants it expands to.

## Verification

- `tests/testthat/test-plan.R`: `normalise_region(list(chrom = "Chr1", start =
  10.5, end = 20L))` and the string `"20.5"` bound are refused with "whole";
  the empty-chromosome refusal and the integer region test are unchanged.
- Before the fix `Rscript dev/test.R '^plan$'` was 189 pass / 2 fail; after it
  191 pass / 0 fail.
