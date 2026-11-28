# 0054 — A stabilising mismatch is a whole distance inside the primer

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0044](0044-a-count-has-one-definition.md) (one
  `is_whole_number()` predicate backs every count), [ADR 0051](0051-a-nucleotide-length-is-whole-before-the-design-grid.md)
  (a nucleotide length is whole before the design grid), [ADR 0045](0045-a-position-is-a-whole-number-before-the-cast.md)
  (a position is whole before the cast)

## Context

`design_as_primers()` is the exported AS-primer grid. Two of its arguments
reach `substr()` through plain arithmetic:

```r
cfgs <- tibble::crossing(strand, as_length = lengths, stab_at = stab_ats) ...
i <- length - stab_at + 1L            # 3'-relative -> 5' index
p <- paste0(substr(p, 1L, i - 1L), stab_base, substr(p, i + 1L, nchar(p)))
```

`lengths` was not checked at all; `stab_ats` was refused only when negative.
Measured on the branch before the fix (76 bp ACGT window, `v_start = 30`,
`ref = "A"`):

```r
design_as_primers(w, 30, "A", "G", lengths = 10.9, stab_ats = NA)
#> as_length column says 10.9; every ref_primer has nchar 11
design_as_primers(w, 30, "A", "G", lengths = 10L, stab_ats = 2.9)
#> the primer is identical to stab_at = 3 (GTACGTAAGA), yet the row says 2.9
#> and arms_stab_policy(2.9) reads as.integer(2.9) = 2 -> "fallback"
design_as_primers(w, 30, "A", "G", lengths = 10L, stab_ats = 0.5)
#> a 10 nt primer whose 3'-terminal base was replaced: the allele-specific
#> base is gone
design_as_primers(w, 30, "A", "G", lengths = 10L, stab_ats = 12)
#> an 11 nt primer for a 10 nt configuration: a base was prepended, not
#> replaced
```

So a fractional position was truncated by `substr()`, then truncated
differently by `as.integer()` in the policy label; a position at or past the
5' end turned a replacement into an insertion. The grid's row described a
molecule it had not built — the same class as ZV75 (CRISPR/RPA lengths) and
ZV68 (positions), in an argument those fixes did not reach.

## Decision

`design_as_primers()` checks the caller's values before building any
configuration:

- `lengths` goes through `check_nucleotide_length(lengths, "lengths", multi =
  TRUE)` (ADR 0051's helper), so it is one or more positive whole numbers of
  nucleotides; the grid uses the returned integer vector.
- each non-`NA` `stab_at` must pass `bad_position(stab_ats, min = 1L)`; the
  values are then cast once.
- each non-`NA` `stab_at` must not exceed `min(lengths)`: a stabilising
  mismatch is a base *inside* the primer, so a position the shortest requested
  primer does not have is refused rather than turned into a longer oligo.

Accepted whole positions and lengths are unchanged: `stab_ats = c(NA, 2L, 3L)`
with `lengths = 20:21` still returns the same 14 configurations, and the
primer length equals the row's `as_length`.

## Alternatives rejected

- **Round the position.** A rounded position is still not the distance the
  caller named, and the policy label would still disagree with the molecule.
- **Clamp a position past the 5' end to the 5' base.** That invents a
  configuration; the caller named a base that is not in the primer.
- **Let `substr()`'s truncation stand and document it.** The row carries
  `as_length`/`stab_at`; a table that describes a different molecule than the
  one returned is the defect, not a convention.

## Consequences

- A fractional, zero, `Inf`, vector or non-numeric `lengths` is refused with
  the package message before the grid.
- A fractional, zero, `Inf` or non-numeric stabilising position is refused
  before the grid; a position past the shortest primer is refused with its own
  message.
- Whole values keep the same primers and the same `stab_policy` labels.
  `design_arms_marker()` keeps its earlier, more specific message for a
  negative `stab_ats` (the signed policy constant).

## Verification

- `tests/testthat/test-as-primer.R`, "a fractional or out-of-primer
  stabilising mismatch is refused before the grid": `2.9`, `0.5` and `Inf`
  as `stab_ats`, `23L` against `lengths = 22L`, and `22.9`/`Inf` as `lengths`
  are refused; the refusal runs before the engine, so the test needs no
  OligoArrayAux.
- the same file's multi-ALT test asserts every `ref_primer` has
  `nchar == as_length`.
- `Rscript dev/test.R "as-primer|arms|z7-priors"` — 334 pass, 0 fail, 2 skip.
