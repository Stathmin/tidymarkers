# 0055 — A common-primer search names whole lengths, gaps and anchor

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0051](0051-a-nucleotide-length-is-whole-before-the-design-grid.md)
  (a nucleotide length is whole before the design grid), [ADR 0050](0050-a-genomic-interval-is-whole-before-the-read.md)
  (a coordinate is whole before the read), [ADR 0044](0044-a-count-has-one-definition.md)
  (one `is_whole_number()` predicate backs every count)

## Context

`design_common_primers()` is the exported common-primer search used by
`design_arms_marker()`. It builds its geometry from three caller values and
reports all three back in the candidate table:

```r
start <- snp_i + grid$gap
end <- start + grid$L - 1L
seq_plus <- substring(as.character(window), cand$start, cand$end)
```

None of the three was checked. `tidyr::crossing()` keeps `20.9` as `20.9`,
base `substring()` truncates it, and the table reports the untruncated value.
Measured on the branch before the fix (160 bp random ACGT window, seed 1,
`snp_i = 80L`, `gaps = 20:40`):

```r
design_common_primers(w, 80L, lengths = 20.9, gaps = 20:40)
#> len column says 20.9, start 104 and end 123.9; nchar(com) is 20
design_common_primers(w, 80.9, lengths = 20L, gaps = 20:40)
#> start column says 100.9; substring() read from 100
```

A gap fraction shifts the reported `start`/`end` the same way. The row
described a geometry and a primer the function had not built — the same class
as ZV78 (`design_as_primers()`), in the search those AS primers are paired
with.

## Decision

`design_common_primers()` checks the caller's values before the grid:

- `snp_i` goes through `check_coord(snp_i, "snp_i")` (ADR 0050's helper), so
  it is one positive whole number.
- `lengths` and `gaps` go through `check_nucleotide_length(..., multi = TRUE)`
  (ADR 0051's helper), so each is one or more positive whole numbers of
  nucleotides; the grid uses the returned integer vectors.

Accepted whole values are unchanged: the default search and the internal
`design_arms_marker()` call still build the same candidates, now with integer
`len`/`start`/`end` columns.

## Alternatives rejected

- **Cast with `as.integer()` at the top.** That hides the fraction from the
  caller and still reports a geometry they did not name; the point is to
  refuse, not to round.
- **Validate only `lengths`.** `gaps` and `snp_i` reach the same arithmetic
  and the same reported columns.
- **Let `substring()`'s truncation stand and document it.** A row carrying
  `len` 20.9 for a 20 nt primer is the defect, not a convention.

## Consequences

- A fractional, zero, `Inf`, vector or non-numeric `lengths`, `gaps` or
  `snp_i` is refused with the package message before any window is read.
- Whole values keep the same candidates; the reported `len`, `start` and
  `end` now match the molecule and the coordinates.
- `design_arms_marker()` inherits the guard through its own `com_gaps`
  argument without a second check.

## Verification

- `tests/testthat/test-design-common-primers.R`, "the lengths, gaps and
  anchor are whole numbers before the grid": `20.9`, `c(20L, 20.5)`,
  `NA_real_` as `lengths`; `20.5` and `c(20L, Inf)` as `gaps`; `8.5` as
  `snp_i` are refused with the package messages. The 16-base window cannot
  yield a candidate, so each refusal is raised before the filter message.
- `tests/testthat/test-design-common-primers.R` and the callers in
  `test-thermo.R` still pass with whole defaults.
- `Rscript dev/test.R "design-common-primers|as-primer|thermo"` — 83 pass,
  0 fail.
