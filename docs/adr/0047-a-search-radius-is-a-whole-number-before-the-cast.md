# 0047 — A search radius is a whole number before the cast

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0042](0042-a-threshold-is-a-count-not-a-number-to-truncate.md)
  (a threshold is a count, not a number to truncate), [ADR 0044](0044-a-count-has-one-definition.md)
  (one `is_whole_number()` predicate backs every count), [ADR 0046](0046-a-design-key-names-every-input-the-design-reads.md)
  (the same defect in the key that names a design)

## Context

Six specificity paths narrowed the caller's search radius with
`as.integer()` before any check:

| path | cast before the check |
|---|---|
| `in_silico_pcr()` | `as.integer(max_mismatch)` |
| `in_silico_pcr_pool()` | `as.integer(max_mismatch)` |
| `exact_hit_table()` | `as.integer(max_mismatch)`, then a check on the cast |
| `specificity_bowtie()` | `as.integer(max_mismatches)`, then the `<= 3` check |
| `specificity_seed()` | `k` and `max_mismatch` used uncast |
| `specificity_of()` | `as.integer(k)` / `as.integer(max_mismatch)` in the seed comparison |

Measured on the branch before the fix:

```r
in_silico_pcr(f, r, genome, max_mismatch = 1.9)   # also accepts -1, Inf, "1"
specificity_seed(primers, genome, k = 2.5)        # seed "CAC" (3 nt)
specificity_seed(primers, genome, max_mismatch = Inf)
#> Error in if (max_mismatch == 0L) ... : missing value where TRUE/FALSE needed
```

`1.9` became `1`: the scan used a narrower radius than the caller named.
That is not a cosmetic difference — a product one substitution away is an
off-target, so the specificity verdict changes. `-1` reached Biostrings'
`'max.mismatch' must be a non-negative integer`, `Inf` became `NA`, and
`"1"` was accepted through coercion. `k = 2.5` silently produced a 3-nt
seed; `k = 0` reached `empty patterns are not supported`.

## Decision

`check_scan_radius(x, arg, min = 0L)` in `R/marker_validate.R` is the one
gate for a radius:

```r
check_scan_radius <- function(x, arg, min = 0L) {
  if (!is_whole_number(x, min = min)) {
    abort_msg("`{arg}` must be one {if (min > 0L) \"positive\" else \"non-negative\"} whole number.")
  }
  as.integer(x)
}
```

Every path above calls it on the value the caller gave, before the cast:

- `max_mismatch` / `max_mismatches` — non-negative whole number;
- `k` — positive whole number (a zero-length seed has no 3' end);
- `specificity_bowtie()` keeps its `<= 3` message after the gate, so the
  engine's own limit is still named.

A character radius is refused, as every other count is under ADR 0044. The
gate returns the integer, so no caller keeps a second cast.

## Alternatives rejected

- **Round the radius instead of refusing it.** `1.9` would search a radius
  the caller never named; the same reason ADR 0042 rejects rounding a count.
- **Keep four separate checks.** That is how the paths drifted apart: two
  checked after the cast, two did not check at all.
- **Accept `"1"` by casting it exactly, as a position does.** A radius is a
  direct function argument, not a CSV column; accepting a string here would
  keep a second admission rule for the same kind of value.

## Consequences

- A fractional, negative, `Inf`, vector or character radius is refused with
  the package message before any engine runs.
- Every previously accepted whole radius behaves as before; no hit count or
  verdict changes for an accepted input.
- The cast and the check live in one function, so a new specificity path has
  one call to make and no rule to copy.

## Verification

- `tests/testthat/test-specificity.R`: `in_silico_pcr()`,
  `in_silico_pcr_pool()`, `exact_hit_table()`, `specificity_seed()` and
  `specificity_of()` refuse `1.9`, `-1`, `Inf`, `"1"` and `c(0, 1)` with
  "non-negative whole number"; `k = 2.5`, `0`, `-1` and `Inf` give "positive
  whole number"; the `0L` / `1L` controls still return a tibble.
- `tests/testthat/test-bowtie.R`: a fractional radius is refused before the
  index is opened, and `4L` still gives "at most 3".
- `Rscript dev/test.R "specificity|insilico|scan|bowtie"` — 494 pass, 0 fail.
