# 0062 — A target screen names its whole window and budget before the grid

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0061](0061-a-target-pick-names-its-whole-tuning-counts.md)
  (a target pick names its whole tuning counts)

## Context

`screen_snp_targets()` in `R/marker_screen.R` is the cheap vectorised
promising/doomed screen that runs before `pick_target_snp()`. Three of its
arguments reach window arithmetic or a comparison without a guard:

```r
lo <- max(1L, min(g$pos) - lmax - clean_3p)      # fractional lo
grid <- tidyr::crossing(strand = strands, length = lengths)
w3 <- IRanges::IRanges(p - clean_3p + 1L, p + clean_3p - 1L)
hot_ok <- hot <= hotzone                         # 7.9 tolerates 7
```

Measured at HEAD on the committed `example_region` fixture:

```r
screen_snp_targets(loci, genome, lengths = 20.9)
#> Error in `check_coord()`: `start` must be one positive whole number;
#>   got 2818.1
screen_snp_targets(loci, genome, clean_3p = 4.9)
#> Error in `check_coord()`: `start` must be one positive whole number;
#>   got 2812.1
screen_snp_targets(loci, genome, hotzone = 7.9)
#> 16/16   # identical to hotzone = 7L: one variant fewer than the caller named
```

The refusal names `start`/`end` from `genome_seq()`, not the argument, so the
caller cannot see which value is wrong; `hotzone = 7.9` is silently narrowed.
ADR 0042–0061 fixed the same class at the designers, the readers and the target
pick; this is the screen that runs before all of them.

## Decision

`screen_snp_targets()` names the three values before it derives the scan window:

```r
lengths <- check_nucleotide_length(lengths, "lengths", multi = TRUE)
clean_3p <- check_scan_radius(clean_3p, "clean_3p", min = 1L)
hotzone <- check_scan_radius(hotzone, "hotzone", min = 0L)
```

`lengths` is the placement series, so a series is accepted; `clean_3p` is a
3'-terminal window and must hold at least one base; `hotzone` is a tolerated
count, so zero is meaningful. The helpers already exist (ADR 0044/0047/0051).

## Alternatives rejected

- **Let `genome_seq()` refuse.** The message names `start`/`end`, which the
  caller never passed, and `hotzone` raises nothing at all.
- **Cast with `as.integer()` first.** That is the defect: `20.9` becomes a
  20 nt placement while `n_tried` and `score` count the caller's number.
- **Guard inside `screen_one_chrom()`.** It runs once per chromosome after the
  window is derived; the wait for `hotzone` is a comparison, not a window.

## Consequences

- A fractional, zero (for `lengths`/`clean_3p`), `Inf`, vector or non-numeric
  value is refused by argument name before any window is read.
- `hotzone = 0` keeps its meaning: only a placement with no other variant.
- Whole calls return the same `n_tried`, `n_ok` and verdict as before.

## Verification

- `tests/testthat/test-screen.R`, "the screen refuses a fractional placement
  before the window is read": `lengths = 20.9`, `c(20L, Inf)`,
  `clean_3p = 4.9`/`0.5`, `hotzone = 7.9`/`-1` are refused by name; a whole
  call still returns one row.
- Red before the fix: the two fractional windows abort inside `genome_seq()`
  with `start`/`end`, and `hotzone = 7.9` returns the integer-call score.
- Green after the fix: `dev/test.R screen` — 23 pass, 0 fail.
