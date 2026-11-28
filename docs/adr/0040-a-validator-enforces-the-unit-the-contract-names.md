# 0040 — A validator enforces the unit the contract names

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0013](0013-metrics-gate-is-not-a-green-delta.md) (the
  metrics gate is not a green delta)

## Context

The metrics contract (§0.1 of the three-day plan) fixes each field's type, not
only its name:

| Field | Contract |
| ----- | -------- |
| `duration_ms` | integer |
| `retries` | integer |
| `tokens_in` / `tokens_out` | integer/null |
| `queue_depth` | integer |
| `gpu_util` | number/null, GPU load as a **fraction** `0`–`1` |

`validate_row()` in `dev/metrics.R` is the gate that decides whether a JSONL
row may enter the daily report. It checked only `>= 0` for the integer fields
and admitted `gpu_util <= 100`, with the message "percentage in [0, 100]". A
hand-written or foreign row therefore passed `--check` with a fractional
duration and a percentage GPU load:

```r
row$duration_ms <- 1.5   # passes
row$queue_depth <- 0.5   # passes
row$gpu_util <- 100      # passes
```

`metrics_row()` never produces these values (`as_int()` coerces and stops on a
fraction), so the hole was only on the validation path — exactly the path that
reads files the writer did not produce. p95 would then be computed over a
mixed-unit column, and `gpu_util` would mix `0.75` and `75` in one report with
no error.

## Decision

1. A field the contract calls an integer is validated as a whole, finite,
   non-negative number. `1.5` is not close enough.
2. `gpu_util` is a fraction in `[0, 1]`; the percentage range is rejected and
   the error message says "fraction", not "percentage".
3. The validator's message names the unit it enforces. A message that says
   "integer" is evidence of an integer check; if the check admits a fraction,
   the message and the check have drifted, and the gate is wrong.

## Alternatives rejected

- **Widen the contract to percentages** because the message already said so.
  The plan is the published contract; the code message is not a second
  authority. `0.75` and `75` in one column would still be ambiguous.
- **Coerce on read** (`round(duration_ms)`) so the report stays clean. That
  silently changes a measurement; a report must refuse an unmeasured or
  malformed duration, as `--emit` already does with no `--ms`.
- **Check integers only in `metrics_row()`.** The writer already does. The
  defect is on the validation path, which reads rows the writer did not make.

## Consequences

- `--check` refuses a row whose duration, retries, token counts or queue depth
  is fractional, and a GPU load outside `0`–`1`.
- Every row in `data/metrics/*.jsonl` (all `gpu_util` are `null`) still
  validates; no report number moves.
- The writer is unchanged; `as_int()` remains the single coercion at the
  production edge.

## Verification

- `tests/testthat/test-metrics-repair.R`, "the validator enforces the unit the
  contract names": `1.5` is refused in each integer field; `gpu_util = 0.75`
  passes and `100` is refused with the fraction message.
- A raw JSONL with `duration_ms = 1.5` and `gpu_util = 100` gives
  `dev/metrics.R --check` two problems and exit 1; before the fix it reported
  `ok`.
- `Rscript dev/test.R "^metrics-repair$"` — 20 pass, 0 fail;
  `Rscript dev/metrics.R --check 2026-11-28` — 36 rows, exit 0.
