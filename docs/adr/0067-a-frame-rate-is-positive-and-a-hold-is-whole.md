# ADR 0067 — A frame rate is positive and finite, a hold is a whole count

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0044 (one whole-number predicate), ADR 0049 (a worker count is
  whole before the cast), ADR 0064 (a terminal window is whole before it is cut)

## Context

`render_design_search()` turns a search into a GIF. Two of its arguments reach
arithmetic and `rep()` raw:

| argument | before the fix |
|---|---|
| `fps =` | `assemble_gif(..., delay_cs = round(100 / fps))` |
| `hold =` | `rep(s, if (s %in% winners) hold else 0L)` while building the frame list |

Measured on R 4.6.1: `fps = 0` gives `delay_cs = Inf`, which ImageMagick
receives as `-delay Inf`; `fps = -1` gives a negative delay; `fps = Inf` gives
a zero delay; `fps = NA` or a vector reaches `round()` and fails without naming
the argument. `hold = 2.9` is truncated by `rep()` to two holds, `hold = -1`
raises base R's `invalid 'times' value`, and `Inf`/`NA`/a vector fail the same
way.

This is the same defect class as ADR 0049: a caller-supplied number crosses
into a writer before the package has said what a valid number is. It is not a
truncating cast on its own, because a frame rate need not be whole.

## Decision

- `hold` is one non-negative whole number: `check_scan_radius(hold, "hold",
  min = 0L)` runs at the entry, before the frame list is built. Zero stays
  legal — it means "do not hold the winner".
- `fps` is one positive finite number: a numeric scalar with `is.finite(fps)`
  and `fps > 0`. A fractional rate stays legal because the GIF delay is
  centiseconds (`fps = 2.5` is a 40 cs delay), so this is a separate rule, not
  `check_scan_radius()`.
- Both checks run before `ggplot2` is loaded, before any frame is drawn and
  before the output name is declared.

## Alternatives rejected

- **Make `fps` a whole number too.** It would refuse a legitimate 2.5 fps and
  confuse a rate with a count; the unit of a GIF delay is not a frame.
- **Clamp `fps` to a positive value.** A caller who asked for zero frames per
  second has a bug; silently choosing a rate hides it.
- **Leave `hold` to `rep()`.** Its error names `times`, not `hold`, and the
  fraction is truncated without any message.

## Consequences

- `fps = 0`, `-1`, `Inf`, `NA`, a character or a vector is refused with
  `` `fps` must be one positive finite number. ``
- `hold = 2.9`, `-1`, `Inf`, a character or a vector is refused with the
  package's whole-number message; `0L` and a whole count behave as before.
- No exported signature changes; `assemble_gif()` keeps its own `delay_cs`
  contract (the internal caller now passes a finite integer).

## Verification

`tests/testthat/test-search-anim.R` refuses the six bad `fps` values and the
five bad `hold` values before any frame, and the existing GIF-name test still
runs. Focused run `Rscript dev/test.R search-anim` —
`FAIL 0 | WARN 0 | SKIP 0 | PASS 30`; `dev/check.R --quick --no-render` passes.
