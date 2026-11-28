# ADR 0085 — The geometry helpers judge the numbers they are built from

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast)

## Context

`primer_interval()`, `as_3p_interval()` and `as_anchor()` are the exported
geometry helpers every design path goes through. They built an
`IRanges`/coordinate arithmetic from their arguments without judging them, and
`IRanges` truncates silently. Measured:

| call | returned |
|---|---|
| `primer_interval(1000.5, "+", 20L)` | `980-1000`, **width 21** for a 20 nt primer |
| `primer_interval(1000L, "+", 20.5)` | width 21, not the 20 the caller named |
| `primer_interval(1000L, "+", 20L, gap = 3.5)` | gap measured as 3 |
| `as_3p_interval(1000.5, "+")` | an interval anchored at 1000 |
| `as_anchor(1000.5, 1L, "+")` | anchor 1000.5, truncated by every caller |

The width case is the worst: a primer of `length = 20` came back as a 21-base
interval, so the reported primer and the primer the caller asked for differ by
one base at the 3' end — the end that discriminates the allele.

## Decision

1. `primer_interval()` judges `snp_pos` and `length` with `is_whole_number()`
   (minimum 1) and `gap` with `is_whole_number(min = 0L)`, before any interval
   is built.
2. `as_3p_interval()` judges `snp_pos` and `window`; `as_anchor()` judges
   `v_start`, `ref_len` and `variant_offset` (minimum 0 for the offset).
3. The messages name the argument; whole doubles (`1000`) still pass, so the
   internal grid, which passes integers, is unaffected.
4. No cast is added: once the value is known whole, `IRanges()` stores it as an
   integer without truncating, and a cast here would only duplicate the check.

## Rejected alternatives

- Guard only the exported callers' inputs: the helpers are exported and
  documented as the single source of truth for geometry, so a direct call is a
  supported entry point.
- Round instead of refuse: `round(1000.5)` is 1000 or 1001 depending on the
  value, so a rounded primer is still a different primer than the caller named.
- Move the checks into `IRanges()` callers only: the width bug is produced by
  the constructor itself, and the caller cannot see that the interval it gets
  back is one base wider.

## Consequences

- A fractional position, length, gap, window or offset is refused by argument
  name instead of producing a shifted interval.
- Measured cost: `primer_interval()` 333 → 341 µs per call on 20k calls (2.5%),
  dominated by the S4 construction the three checks do not touch.
- `block_in_seg()` is left unguarded: it is internal and receives `ref_len` and
  `primer_len` that the design path already judged; it is not exported.

## Verification

`Rscript dev/test.R '^geometry$'` — 26 pass, 0 fail, including eight new
refusals and the whole-double/zero-offset pass-through.
`Rscript dev/test.R '^(mnp|arms|design-common-primers|design-metrics)$'` — 177
pass, 0 fail: no design path passes a fractional value.
`Rscript dev/check-casts.R` — `53 position/threshold cast(s), 22 guarded in
place, 31 registered with a reason; self-test passed`.
