# ADR 0073 — A diagnostic mask position is whole in its spelling

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0068
  (a position is whole in its own spelling), ADR 0071 (a region bound is a
  coordinate before the cast)

## Context

`normalise_diagnostic_positions()` accepts the subgenome-diagnostic mask of a
KASP construct as a data frame, a list with `chrom`/`pos`, a `GRanges`, or a
named list `chrom -> positions`. The non-`GRanges` branches cast `pos` with
`as.integer()` and validated only the result:

```r
chrom <- as.character(x$chrom)
pos <- as.integer(x$pos)
out <- tibble::tibble(chrom = ath_chrom(chrom), pos = pos)
if (anyNA(out$chrom) || anyNA(out$pos) || any(out$pos < 1L)) abort_msg(...)
```

A fractional position therefore reached the mask as a plausible base:

| input | before the fix |
|---|---|
| `data.frame(chrom = "Chr1", pos = 10.9)` | mask position `10` |
| `list(chrom = "Chr1", pos = 10.9)` | mask position `10` |
| `list(Chr1 = c(10.9, 20))` | mask positions `10`, `20` |

The mask decides whether the common primer's 3' base is reported as
`unique`/`multicopy`; a shifted coordinate claims a diagnostic position the
caller never named. The `GRanges` branch is unaffected: its start is an integer
by construction.

## Decision

- The data-frame and both list branches check `bad_position(pos_raw)` **before**
  `as.integer()`, and refuse a fractional, `NA`, `Inf`, non-numeric or
  non-positive position with the mask message naming a "positive whole-number
  position".
- The `GRanges` branch keeps its current cast; the existing chromosome and
  `pos >= 1` check stays as the last rule for every branch.

## Alternatives rejected

- **Check the cast value.** `as.integer(10.9)` is `10`, a plausible base; only
  the spelling shows that the caller named another one.
- **Round the position.** A rounded mask flags a base the caller did not name.
- **Validate only after the cast.** The existing check runs on truncated values
  and cannot recover the fraction.

## Consequences

- `assemble_kasp_assay(diagnostic_positions = )` and
  `assemble_primer_assay(diagnostic_positions = )` report a mask the caller
  wrote, or refuse it.
- Integer masks and `GRanges` masks behave exactly as before.

## Verification

- `tests/testthat/test-y32-homeolog.R`: the data frame, the list, the named list
  and the character `"10.9"` are refused with "whole-number position"; the
  integer mask still returns `c(10L, 20L)`.
- `Rscript dev/test.R "y32-homeolog"` — 69 pass / 0 fail.
