# 0057 — A LAMP length is checked before the swallowed grid

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0056](0056-the-shared-primer-grid-names-whole-lengths.md)
  (the shared primer grid names whole lengths), [ADR 0051](0051-a-nucleotide-length-is-whole-before-the-design-grid.md)
  (a nucleotide length is whole before the design grid)

## Context

ADR 0056 put `check_nucleotide_length()` in `primer_candidates()`, so a
fractional `lengths` is refused before the grid. `design_lamp_marker()`
passes its `arm_lengths` and `outer_lengths` through `lamp_arm_sites()` and
then `safe_candidates()`, which wraps the grid in a blanket handler:

```r
out <- tryCatch(
  primer_candidates(seq, starts, lengths, strand, tm_range, gc_range, ...),
  error = function(e) NULL
)
if (is.null(out)) empty else out
```

The handler was written for the composition screener's empty-set case, but it
catches every error, including the new input refusal. Measured on the branch
before the fix (`snp_fixture()`, `flank = 150L`):

```r
design_lamp_marker(f$genome, f$snp, arm_lengths = 20.9)
#> No LAMP primer set passed the variant-placement / inner-Tm / structure /
#> cross-dimer / amplicon filters.
design_lamp_marker(f$genome, f$snp, outer_lengths = 21.9)
#> (the same "no set passed" message)
```

Neither call raised the package message that names the bad length. A caller
who typed a fraction was told the locus had no design — a false negative that
invites changing the locus rather than the argument. ADR 0056 had assumed the
grid's refusal reached the caller for LAMP; it does not, because
`safe_candidates()` sits in front of it.

## Decision

`design_lamp_marker()` checks both length arguments before it builds the
template:

```r
arm_lengths <- check_nucleotide_length(arm_lengths, "arm_lengths", multi = TRUE)
outer_lengths <- check_nucleotide_length(outer_lengths, "outer_lengths", multi = TRUE)
```

The check is at the exported entry, so the message names the argument the
caller actually passed. `safe_candidates()` keeps its blanket handler and its
empty-set contract: changing it touches every LAMP/RPA probe path and is a
separate decision.

## Alternatives rejected

- **Remove the `tryCatch()` in `safe_candidates()`.** It is the root cause,
  but the handler is the contract three call families were built against; a
  bad length would then be reported as the grid's `lengths` even when the
  caller passed `arm_lengths`, and any other grid error would start
  propagating through the LAMP ladder.
- **Check only `arm_lengths`.** `outer_lengths` reaches the same swallowed
  grid and produced the same false "no set passed".
- **Leave it and document the message.** "No set passed" is wrong for a
  malformed argument; a false negative about the locus is worse than a
  refusal about the input.

## Consequences

- A fractional, zero, `Inf`, vector or non-numeric `arm_lengths` or
  `outer_lengths` is refused with "`arm_lengths`/`outer_lengths` must be one
  or more positive whole numbers of nucleotides" before any Tm work.
- Whole values keep the same LAMP sets and the same ladder behaviour.
- `safe_candidates()` still hides other grid errors from LAMP callers; that
  limit is recorded here rather than changed in this fix.

## Verification

- `tests/testthat/test-lamp.R`, "a fractional LAMP arm or outer length is
  refused, not read as no set": `arm_lengths = 20.9` and
  `outer_lengths = 21.9` are refused with "positive whole numbers of
  nucleotides".
- `dev/test.R "lamp|primer-pairs|hrm"` — 693 pass, 0 fail, 0 skip.
- `design_lamp_marker()` with the defaults still returns a 4-primer set in
  the same test file.
