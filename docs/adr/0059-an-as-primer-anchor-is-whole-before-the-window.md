# 0059 — An AS-primer anchor and linked offset are whole before the window

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0058](0058-a-repeat-scan-names-whole-units-and-copies.md)
  (a repeat scan names whole units and copies), [ADR 0056](0056-the-shared-primer-grid-names-whole-lengths.md)
  (the shared primer grid names whole lengths)

## Context

`design_as_primers()` in `R/marker_design.R` passes its geometry straight to
`IRanges` and `Biostrings`. `IRanges` truncates a fractional start, and the
function's own `as.integer(linked$offset)` truncates again, so a fractional
request builds a molecule for a different coordinate than the caller named.
Measured at HEAD with `set.seed(3)`, a 200 bp random window, SNP at 100:

```r
design_as_primers(W, 100.9, ref = "A", alt = "G", strands = "+")
#> tmpl_start 81, tmpl_end 100        # exactly the 100L result
design_as_primers(W, 100L, ref = "A", alt = "G", strands = "+",
                  variant_offset = 1.9)
#> tmpl_start 82, but alt_primer "TTCATTCAACACTGGCCAGG"
#> vs 1L's "...CCCAG"                  # anchor truncated, block placed
linked <- tibble::tibble(offset = 96.9, alt = "C")
#> exactly the 96L result
linked$offset = 0.5 / NA / 201L
#> raw Biostrings/IRanges errors: "some ranges in 'at' are off-limits",
#> "'start' or 'width' cannot contain NAs"
```

The reported `tmpl_start`/`tmpl_end` and the molecule therefore disagree for
`variant_offset = 1.9`, and a bad `linked$offset` reaches a dependency error
instead of a package message. ADR 0042–0058 fixed the same class at other
entry points; this is the AS-primer geometry that those guards do not cover.

## Decision

`design_as_primers()` names the caller's geometry before it builds a template:

```r
v_start <- check_coord(v_start, "v_start")
variant_offset <- check_scan_radius(variant_offset, "variant_offset", min = 0L)
if (!is.null(linked) && nrow(linked) > 0L) {
  off <- suppressWarnings(as.integer(linked$offset))
  bad <- bad_position(linked$offset, min = 1L) | (!is.na(off) & off > length(window))
  if (any(bad)) {
    abort_msg(
      "`linked$offset` must name whole positions inside `window` (1..{length(window)}).",
      i = "Got {paste(format(linked$offset[bad]), collapse = ', ')}."
    )
  }
}
```

`check_coord()` and `check_scan_radius()` already exist (ADR 0042/0052);
`bad_position()` is the shared vectorised predicate used by the neighbouring
`stab_ats` guard. `variant_offset` is a distance from the variant block, so `0`
is legal and the bound is non-negative rather than positive.

## Alternatives rejected

- **Round the values.** The primer would still be built at a coordinate the
  caller did not name, and the reported geometry would still describe the
  rounded value.
- **Guard inside `as_anchor()`/`make_as_primer()`.** The window replacement for
  `linked` happens before either call, so the earliest truncation would remain;
  the shared builder is the wrong seam for a range check.
- **Let `IRanges` raise.** Its messages name no package argument and the
  fractional case raises nothing at all.

## Consequences

- A fractional, zero, `Inf`, vector or non-numeric `v_start` is refused by name;
  a negative, fractional, `Inf` or vector `variant_offset` is refused by name.
- A `linked$offset` that is fractional, `NA`, below 1 or past
  `length(window)` is refused with one message naming the window's range.
- Whole calls return the same primers and the same rows as before.
- `linked` with zero rows keeps the current pass-through (it is not a table of
  variants, and the template is unchanged).

## Verification

- `tests/testthat/test-as-primer.R`, "the AS-primer geometry names whole
  positions before the window": `v_start = 30.9`, `variant_offset = 1.9`/`Inf`,
  and `linked$offset` of `20.9`/`NA`/`0L`/`77L` are refused; `77L` names the
  window's 1..76 range.
- Red before the fix: `dev/test.R as-primer` — 2 failures (`30.9` and `1.9` did
  not throw) and 1 error (`Inf` reached the `IRanges` comparison).
- Green after the fix: `dev/test.R as-primer` — 44 pass, 0 fail.
- `Rscript dev/check.R --quick --no-render` passed on the preceding commit.
