# 0065 — Two pool limits are whole base-pair counts before the verdict

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0047](0047-a-search-radius-is-a-whole-number-before-the-cast.md)
  (a search radius is a whole number before the cast)

## Context

`primer_set_compatibility()` is the exported "can these primers share a tube"
answer, and the base-pair rules behind `marker_tubes()` call it. It compared
its two limits raw:

```r
d$pass <- dplyr::if_else(d$kind == "self", d$run <= max_selfcomp, d$run < min_cross_run)
```

`run` is a whole base-pair count, so a fractional limit silently changes the
verdict rather than merely being imprecise:

```r
p <- c(AS1 = "GAAGGTGACCAAGTTCATGCT", AS2 = "GAAGGTCGGAGTCAACGGATT")
# a 6 bp cross run under the default 6L is a failure (6 < 6 is FALSE)
primer_set_compatibility(p, min_cross_run = 6.5)
#> the same 6 bp row now passes (6 < 6.5 is TRUE)
```

`max_selfcomp = 5.4` is indistinguishable from `5L`, `Inf`/`NA` make the verdict
all-fail or error without naming an argument, and a vector limit recycles. The
caller named a threshold the function did not honour.

## Decision

Both limits are named before the pairwise table is built:

```r
min_cross_run <- check_scan_radius(min_cross_run, "min_cross_run", min = 0L)
max_selfcomp <- check_scan_radius(max_selfcomp, "max_selfcomp", min = 0L)
```

Zero stays legal: it deliberately disables the corresponding rule (`run < 0`
never fails, `run <= 0` only fails a zero run). The rule is about whole counts,
not about a minimum useful value.

## Alternatives rejected

- **Guard in `marker_tubes()` only.** It is one caller; the exported function is
  the trust boundary and its other callers keep the defect.
- **`as.integer()` both limits.** That hides the caller's mistake and still
  turns `Inf` into `NA`.
- **Floor `min_cross_run` at 1.** A zero threshold is an explicit "do not
  screen cross dimers" switch already reachable with `0L`; refusing it would be
  a behaviour change outside this defect.

## Consequences

- A fractional, negative, `Inf`, `NA`, vector or non-numeric limit is refused by
  name before any pair is scored.
- Whole calls, including `min_cross_run = 0L` and `max_selfcomp = 0L`, return
  the same tables as before.

## Verification

- `tests/testthat/test-pool.R`, "primer_set_compatibility names its two
  base-pair limits": `6.5`, `-1`, `NA`, `Inf`, a vector and `"6"` are refused
  for both arguments; `min_cross_run = 6, max_selfcomp = 5` equals the default
  call.
- Red before the fix: `min_cross_run = 6.5` flipped a 6 bp cross row to pass;
  the error cases named no argument.
- Green after the fix: `dev/test.R pool` — 37 pass, 0 fail.
