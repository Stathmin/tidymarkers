# ADR 0069 — A GenBank location is read whole or refused

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0068
  (a position column is whole in its own spelling)

## Context

`read_genbank()` parses the FEATURES location of each record with
`parse_location()`. Every boundary is cast with `suppressWarnings(as.integer())`
and the result is kept only when it is not `NA`:

```r
s <- suppressWarnings(as.integer(bounds[1]))
e <- suppressWarnings(as.integer(bounds[2]))
if (!is.na(s)) { starts <- c(starts, s); ends <- c(ends, e); ... }
```

`as.integer()` truncates, so the guard can only see `NA`, never a fraction.
Measured on the branch:

| location | returned |
|---|---|
| `100.5..200` | `start 100`, `end 200` |
| `join(100..200, 300..400.9)` | `start 100`, `end 400` |
| `2801.5` | `start 2801`, `end 2801` |
| `abc` | `start NA`, `end NA` |

A feature at a coordinate the file never named is published as if it were at
the truncated one; a feature whose location is unreadable becomes a row with
`NA` bounds instead of an error. `marker_features()` already refuses the same
input on its two branches (ADR 0045 amendment, ZV93); the GenBank reader is
the remaining member of that family.

## Decision

- Every boundary token in a location goes through `bad_position()` before the
  cast. A token that is not a positive whole number refuses the record with
  `` GenBank location `<location>` is not a whole range. ``, naming the
  location as it appears in the file.
- The check sits in `parse_location()`, the one place that turns text into
  boundaries, so `join`, `order`, `complement`, `^` and a single base all
  inherit it.
- Fuzzy markers `<` and `>` keep their meaning (the number they qualify is the
  boundary); a location that names no whole boundary is refused rather than
  published with `NA` bounds.

## Alternatives rejected

- **Keep the truncated coordinate.** It is a different locus: a primer or an
  overlap check would read a window the file does not describe.
- **Keep `NA` bounds for an unreadable part.** `NA` in a `start`/`end` column
  is indistinguishable from a feature with no location; the caller cannot tell
  a malformed record from a missing one.
- **Fix only the decimal case.** `as.integer("1e3")` is 1000 as well; the
  spelling rule of ADR 0068 covers both without a second predicate.

## Consequences

- A GenBank file with a fractional, exponent-spelled or unreadable feature
  boundary is refused by name instead of producing shifted or `NA` coordinates.
- Every whole location (`1..120`, `complement(10..60)`, `join(10..30,40..60)`)
  parses exactly as before, so no accepted record changes.
- The refusal happens on read, before the `.gbidx` cache is trusted: a cached
  index of a previously accepted file cannot reintroduce the truncated row,
  because the cache stores the parsed tbl and is only reused for the same
  file size and mtime.

## Verification

- `tests/testthat/test-io.R`: a record with `CDS 100.5..120` is refused with
  the location in the message; the existing parse, cache and mixed-strand
  fixtures still pass.
- Focused run `Rscript dev/test.R "^io$"` — 232 pass, 0 fail.
