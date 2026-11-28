# 0060 — A target window names its whole locus before the cut

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0059](0059-an-as-primer-anchor-is-whole-before-the-window.md)
  (an AS-primer anchor is whole before the window)

## Context

`snp_template()` in `R/marker_target.R` builds the window every single-locus
designer works on. It checked `flank` but cast the locus with
`as.integer(pos)` on both sides of the window, without a guard:

```r
lo <- max(1L, as.integer(pos) - flank)
hi <- min(len, as.integer(pos) + flank)
```

Measured at HEAD on a 100 bp window, `flank = 10L`:

```r
snp_template(g, "Chr1", 50.9, 10L)
#> snp_i 11, same sequence as 50L      # window built for another SNP
snp_template(g, "Chr1", 0, 10L)
#> snp_i -9                            # variant outside the cut window
design_ssr_marker(g, snp(pos = 150.9), flank = 150L)
#> identical result to pos = 150L
```

`strand_gc_hint()` had the same cast with a worse message: a fractional `pos`
or `flank` reached `get_region()` and was refused as `start`/`end`, naming an
argument the caller never passed.

`snp_template()` is the shared seam for the 26 designer call sites, so one
guard reaches every family at once, the same shape as ADR 0056 and ADR 0058.

## Decision

`snp_template()` names the locus before it cuts a window:

```r
pos <- check_coord(pos, "pos")
```

and `strand_gc_hint()` names both of its coordinates:

```r
pos <- check_coord(pos, "pos")
flank <- check_base_pairs(flank, "flank")
```

`check_coord()` already exists (ADR 0042/0059); `check_base_pairs()` is the
guard `snp_template()` already used for `flank`.

## Alternatives rejected

- **Check `snp$pos` in each designer.** That duplicates the guard at 26 call
  sites and leaves an internal caller free to reintroduce the truncation.
- **Cast with `round()`.** The window would still be built for a SNP the table
  does not name.
- **Leave `strand_gc_hint()` to `get_region()`.** It does abort, but the
  message names `start`/`end`, so the caller cannot tell which argument to fix.

## Consequences

- A fractional, zero, negative, `Inf`, vector or non-numeric `pos` is refused
  with "`pos` must be one positive whole number" before any window is cut.
- A `pos` outside the chromosome still reaches the existing clip: the window is
  the part the sequence has, which is the documented behaviour.
- Whole calls return the same windows and the same markers as before.
- Every designer built on `snp_template()` inherits the refusal; none needs its
  own copy of the check.

## Verification

- `tests/testthat/test-genome.R`, "snp_template names whole coordinates before
  the window is cast": `pos = 20.9` and `pos = 0` are refused by name;
  `strand_gc_hint()` refuses a fractional `pos` by name and a fractional
  `flank` as base pairs; the whole case still returns a 21 bp window.
- Red before the fix: `dev/test.R genome` — 3 failures and 1 error (the two
  fractional `pos` calls and `pos = 0` did not throw; `strand_gc_hint()`
  reported `end` instead of `pos`).
- Green after the fix: `dev/test.R genome` — 174 pass, 0 fail, 1 skip.
