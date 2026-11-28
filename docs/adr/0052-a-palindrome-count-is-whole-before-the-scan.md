# 0052 — A palindrome's counts are whole before `findPalindromes()`

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0044](0044-a-count-has-one-definition.md) (one
  `is_whole_number()` predicate backs every count), [ADR 0047](0047-a-search-radius-is-a-whole-number-before-the-cast.md)
  (a search radius is whole before the scan), [ADR 0051](0051-a-nucleotide-length-is-whole-before-the-design-grid.md)
  (the same defect in the design grid)

## Context

`has_palindrome()` is the exported QC entry that answers "does this oligo form
an inverted repeat". Its three counts — `min_armlength`, `max_looplength`,
`max_mismatch` — went to `Biostrings::findPalindromes()` unvalidated; only the
fast path's own predicate looked at them, and it looked at them *after*
deciding whether to take that path.

Measured on the branch before the fix:

```r
has_palindrome("ACGTACGT", max_mismatch = 1.9)
#> TRUE          # silently a 1-mismatch budget, not the 1.9 asked for
has_palindrome("ACGTACGT", max_mismatch = NA)
#> Error: 'max.mismatch' must be a single integer
has_palindrome("ACGTACGT", max_mismatch = -1L)
#> Error: 'max.mismatch' must be a non-negative integer
has_palindrome("ACGTACGT", max_mismatch = c(0L, 1L))
#> Error: 'max.mismatch' must be a single integer
has_palindrome("ACGTACGT", min_armlength = 0L)
#> Error: 'min.armlength' must be >= 2
```

The fractional case is the dangerous one: a tolerance nobody named was accepted
and *loosened* the answer to TRUE, while the other malformed counts raised
Biostrings' messages instead of the package's. This is the same silent
truncation ADR 0047 removed from the specificity radius, on a path that never
received the guard.

## Decision

`has_palindrome()` validates every count with the shared `check_scan_radius()`
before the sequence is looked at:

- `max_mismatch` — non-negative whole number (0 is the strict default);
- `min_armlength` — positive whole number (the Biostrings floor of 2 stays
  where it is, so the same values reach the same engine);
- `max_looplength` — non-negative whole number.

The fast-path predicate then only has to test the mismatch tolerance
(`identical(max_mismatch, 0L)`) and the arm floor (`min_armlength >= 2L`),
because the values are already integers.

## Alternatives rejected

- **Validate inside `findPalindromes()`.** Not ours, and the fast path never
  reaches it.
- **Keep the `is.numeric()`/`floor()` conjunction and only fix the fractional
  case.** That is a second definition of a whole number, the thing ADR 0044
  removes; it also keeps the engine messages for `NA`, negative and vector
  counts.
- **Round the tolerance.** A rounded mismatch budget is still a budget the
  caller did not name.

## Consequences

- A fractional, negative, `NA`, `Inf`, vector, non-numeric, or too-small count
  is refused with the package message before any scan.
- Accepted whole counts produce the same answer as before; `max_mismatch = 0L`
  keeps the fast path.
- The three `@param` entries now state the unit and the admission rule.

## Verification

- `tests/testthat/test-palindromes.R`, "has_palindrome refuses a count it would
  have handed to findPalindromes": `1.9`, `-1`, `NA`, `Inf`, `c(0, 1)`, `"0"`
  as `max_mismatch`, `0`/`4.5` as `min_armlength`, `2.5` as `max_looplength` are
  refused with "whole number"; `ACGTACGT` is still TRUE, `ACGTTGCA` still
  FALSE, and `max_mismatch = 1L` equals `max_mismatch = 1`.
- `Rscript dev/test.R "^palindromes$"` — 51 pass, 0 fail; `"^thermo$|^mcp$"` —
  894 pass, 0 fail.
