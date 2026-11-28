# ADR 0077 — A locus row and its range bridge are whole in their spelling

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0068
  (a position is whole in its own spelling), ADR 0076 (a panel variant table is
  checked at its seam)

## Context

Two seams still cast a position that a caller supplies, without judging its
spelling: `locus_row()` and `marker_as_granges()`.

`locus_row()` is the one canonical row behind three exported entry points —
`design_markers()`, `design_marker_batch()` and `design_kasp_markers()`. It cast
with `as.integer(loc$pos)`, so a locus table row at `10.9` was designed at `10`,
and `-5` reached the designer as a position:

| caller's row | before the fix |
|---|---|
| `Chr1:10.9 A>G` | designed at `Chr1:10` |
| `Chr1:-5 A>G` | accepted as `Chr1:-5` |

`marker_as_granges()` is exported and accepts any variant table. It cast with
`as.integer(df$pos)` the same way, so `10.9` became a `GRanges` start at 10 and
`-5` became a negative range.

Both are silent shifts: the caller names one locus and the package designs or
reports another.

## Decision

- `locus_row()` judges `loc$pos` with `bad_position()` before the cast and
  refuses a value that is not a positive whole number: "`loc` needs a positive
  whole `pos`; got 10.9."
- `marker_as_granges()` judges `df$pos` before the cast and refuses with the
  writer's message, naming the row: "Variant 1 has a POS outside the positive
  integer range (10.9)."
- A batch keeps its contract: `design_kasp_markers()` calls `locus_row()` inside
  its per-locus `tryCatch`, so a bad row becomes an `error` verdict with `NA`
  position instead of aborting the batch. `design_markers()`' summary reports
  such a row as `NA`, never as a truncated neighbour.

## Rejected alternatives

- **Validate in each exported entry point.** The three designers share one
  builder; a guard in one of them would leave the other two open.
- **Cast and compare.** `as.integer("10.9")` is `10`, a legal position; only the
  spelling shows the caller named something else (ADR 0068).
- **Drop the bad row.** The batch contract is one verdict per input row; a
  dropped row changes the row count and hides the error.

## Consequences

- `design_markers()`, `design_marker_batch()` and `design_kasp_markers()` design
  exactly the locus they were given, or report the row as an error.
- `marker_as_granges()` builds a range only from a positive whole position.
- Integer tables and `GRanges` inputs behave exactly as before.

## Verification

- `tests/testthat/test-batch.R`: `locus_row()` keeps `10L` and refuses `10.9`
  and `-5`; the full batch suite passes (63 pass / 0 fail).
- `tests/testthat/test-io.R`: `marker_as_granges()` refuses a table at `10.9`
  and `-5`; `Rscript dev/test.R "^io$"` — 247 pass / 0 fail.
