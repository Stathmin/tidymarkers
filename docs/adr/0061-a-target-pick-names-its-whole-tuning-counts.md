# 0061 — A target pick names its whole tuning counts before the scan

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0060](0060-a-target-window-names-its-whole-locus.md)
  (a target window names its whole locus)

## Context

`pick_target_snp()` in `R/marker_target.R` ranks loci before one is designed.
Three of its tuning counts reached a dependency without a guard:

```r
IRanges::IRanges(.data$pos - clean_3p + 1L, .data$pos - 1L)   # truncates
dplyr::slice_head(short, n = max(1L, as.integer(max_tries)))  # truncates
dplyr::slice_head(ranked, n = max(1L, specificity_candidates)) # raw dplyr error
```

Measured at HEAD on the cheap target fixture (80 bp `ACGT` repeat, gene covers
it, one substitution at 20):

```r
pick_target_snp(v, g, gene, specificity_candidates = 0L, clean_3p = 2.9)
#> identical to clean_3p = 2L     # the guard band is two bases, not 2.9
pick_target_snp(v, g, gene, specificity_candidates = 0L, max_tries = 0.5)
#> returned a target             # as.integer(0.5) is 0, then max(1L, 0)
pick_target_snp(v, g, gene, specificity_candidates = 1.9)
#> Error: `n` must be a round number, not the number 1.9.
```

ADR 0042–0060 fixed the same class at other entry points; this is the ranking
pre-filter that precedes the target window.

## Decision

`pick_target_snp()` names the three counts after the contract block, before any
variant row is read:

```r
clean_3p <- check_nucleotide_length(clean_3p, "clean_3p")
specificity_candidates <- check_scan_radius(
  specificity_candidates, "specificity_candidates", min = 0L
)
max_tries <- check_scan_radius(max_tries, "max_tries", min = 1L)
```

`specificity_candidates = 0` is the documented "skip the seed scan", so its
bound is non-negative; `max_tries` is one or more attempts. `check_nucleotide_length()`
and `check_scan_radius()` already exist (ADR 0042/0058/0059).

## Alternatives rejected

- **Cast with `max(1L, as.integer(x))`.** That invents an attempt the caller did
  not ask for, and `clean_3p = 2.9` would still build the 2-base band.
- **Let `slice_head()` raise.** Its message names `n`, not the argument, and
  `clean_3p` raises nothing at all.
- **Guard inside the ranking helpers.** The three counts are the public entry's
  own contract; the helpers are shared with already-checked callers.

## Consequences

- A fractional, zero (for `clean_3p`/`max_tries`), negative, `Inf`, vector or
  non-numeric tuning count is refused by name before the scan.
- `specificity_candidates = 0` keeps its documented "skip the seed scan".
- Whole calls return the same target and the same ranking as before.

## Verification

- `tests/testthat/test-target.R`, "the target pick names its whole tuning counts
  before the scan": `clean_3p = 2.9`, `specificity_candidates = 1.9`,
  `max_tries = 0.5` and `max_tries = Inf` are refused by name; the whole call
  still returns the fixture's `id`.
- Red before the fix: `dev/test.R target` — 1 failure (`2.9` built the 2-base
  band) and 1 error (`1.9` reached dplyr).
- Green after the fix: `dev/test.R target` — 43 pass, 0 fail.
