# ADR 0078 — A plotted feature interval is whole in its spelling

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0068
  (a position is whole in its own spelling), ADR 0077 (a locus row and its
  range bridge are whole in their spelling)

## Context

`normalise_features()` (`R/marker_plots_map.R`) prepares the `features` table
that `plot_primers()`, `plot_marker_map()` and the doomed-target context draw as
a gene/feature row. It cast both bounds with `as.integer(features$start)` and
`as.integer(features$end)` without judging their spelling:

| caller's feature row | before the fix |
|---|---|
| `start = 10.9, end = 50` | arrow began at 10 |
| `start = -5, end = 50` | negative arrow accepted |

A feature file read through `marker_features()` carries integer bounds, but the
`features` argument is public and also accepts any data frame or list. A
fractional bound therefore moved the annotation to a neighbour the caller never
named, and a negative bound drew an interval outside the chromosome.

## Decision

- `normalise_features()` judges `features$start` and `features$end` with
  `bad_position()` before the cast and refuses a bound that is not a positive
  whole number: "`features` needs a positive whole `start`/`end` on every row."
  Every row is judged, so one bad row among good ones is refused too.
- A missing or empty feature table keeps its current meaning (no annotation
  row), so only a supplied bound is judged.

## Rejected alternatives

- Round with `round()`: still designs an interval the caller did not name.
- Clamp to the plot window: hides the malformed input and draws a wrong arrow.
- Stay tolerant because the plot is diagnostic: a figure is evidence, and a
  shifted gene arrow is a wrong one.

## Consequences

- A malformed feature table fails where it enters the plot with a message that
  names both bounds.
- Every valid feature table — integer, `GRanges`-derived or file-derived — draws
  exactly as before.

## Verification

`tests/testthat/test-plots.R` refuses a fractional second row (`start = pos -
10.5` beside a whole first row) and a negative `start = -5L` in the feature-row
test and still renders the integer table and the gene record:
`Rscript dev/test.R "plots"` — 97 pass, 0 fail. `tests/testthat/test-alignment.R`
draws its two-row feature track — 52 pass, 0 fail.
