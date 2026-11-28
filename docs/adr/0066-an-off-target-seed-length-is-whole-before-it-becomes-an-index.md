# 0066 — An off-target seed length is whole before it becomes an index

- Status: accepted
- Date: 2026-11-28
- Supersedes: none

## Context

`offtarget_profile()` is the exported 3'-end-aware grading of a hit table; its
`seed_len` override reaches `mismatch_weights()` in
`R/marker_specificity_policy.R`, which treats it as a raw vector index:

```r
seed <- seq.int(max(1L, n_nt - seed_len + 1L), n_nt)
w[seed] <- pol$weight_seed
```

Measured on a 21 nt primer and a hit table from `exact_hit_table()`:

```r
mismatch_weights(21L, "pcr", 0.5)   # seeds nothing: same as 0, no refusal
mismatch_weights(21L, "pcr", 4.7)   # truncated to the 4 nt seed
mismatch_weights(21L, "pcr", Inf)   # every position seeds: length 21, all 3
mismatch_weights(21L, "pcr", -1L)   # length 23: the weight vector grew past
                                    # the primer, indexed by positions 22-23
```

The consequence is a real grading change, not a rounding: on an inner mismatch
at position 10, `offtarget_profile(..., seed_len = Inf)` reported penalty 3
where the policy (seed 5) reported 1. `NA` is the documented sentinel "use the
policy's seed" and must stay legal.

## Decision

`offtarget_profile()` names `seed_len` before any profile is computed:

```r
seed_len <- if (isTRUE(is.na(seed_len))) NA_integer_ else
  check_scan_radius(seed_len, "seed_len", min = 1L)
```

`isTRUE()` keeps the NA sentinel scalar-only: a vector containing `NA` is not
the sentinel and is refused. One is the floor (a zero-length seed is not a
seed).

## Alternatives rejected

- **Guard inside `mismatch_weights()`.** It sees the count but not the
  argument name the caller used, and `NA` means something different there (the
  policy default) than it does at the exported entry.
- **`as.integer(seed_len)`.** That is the silent truncation the defect is
  about; `Inf` still becomes `NA` and `-1` still grows the vector.
- **Refuse `NA`.** It is the documented way to ask for the scheme's seed and
  the default of both `offtarget_profile()` and `mismatch_profile()`.

## Consequences

- A fractional, zero, negative, `Inf`, vector or non-numeric `seed_len` is
  refused by name; `NA` and whole `>= 1` values behave as before.
- `mismatch_weights()` itself stays unguarded as an internal helper; only the
  exported entry is a trust boundary.

## Verification

- `tests/testthat/test-specificity-3p.R`, "offtarget_profile names its seed
  length before the weights": `0.5`, `20.9`, `Inf`, `-1`, a vector and `"5"`
  are refused by name; `seed_len = 5L` equals `5`; the default equals
  `NA_integer_`.
- Red before the fix: `Inf` changed the inner-mismatch penalty from 1 to 3 and
  `-1L` returned a 23-long weight vector for a 21 nt primer.
- Green after the fix: `dev/test.R specificity` — 411 pass, 0 fail, 2 skip.
