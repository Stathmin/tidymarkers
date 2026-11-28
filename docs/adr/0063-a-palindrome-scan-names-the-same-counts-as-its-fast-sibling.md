# 0063 — A palindrome scan names the same counts as its fast sibling

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0052](0052-a-palindrome-count-is-whole-before-the-scan.md)
  (a palindrome count is whole before the scan)

## Context

`has_palindrome()` in `R/marker_thermo.R` was given the whole-number guard for
its three counts (ADR 0052). `find_palindromes()` — the exported function the
same file documents as its slow, hit-table sibling — was not:

```r
h <- Biostrings::findPalindromes(
  Biostrings::DNAString(s),
  min.armlength = min_armlength, max.looplength = max_looplength,
  max.mismatch = max_mismatch
)
```

Measured at HEAD on `c(hit = "ACGTACGT", none = "AACCGG")`:

```r
find_palindromes(seq, max_mismatch = 1.9)    # 1/0 hits: the same as 1L
find_palindromes(seq, min_armlength = 4.9)   # 1/0 hits: a 4 nt arm
find_palindromes(seq, max_looplength = 2.9)  # 1/0 hits: a 2 nt loop
find_palindromes(seq, min_armlength = 0L)
#> Error in `.find_palindromes()`: 'min.armlength' must be >= 2
find_palindromes(seq, max_mismatch = Inf)
#> Error in `if (max.mismatch < 0L)`: missing value where TRUE/FALSE needed
```

A fraction narrows the scan silently; zero and `Inf` raise messages that name
`min.armlength` or nothing, not the argument the caller passed.

## Decision

`find_palindromes()` names the three counts at its entry, with the same bounds
`has_palindrome()` uses:

```r
min_armlength <- check_scan_radius(min_armlength, "min_armlength", min = 1L)
max_looplength <- check_scan_radius(max_looplength, "max_looplength", min = 0L)
max_mismatch <- check_scan_radius(max_mismatch, "max_mismatch")
```

The sibling calls share one rule; the pair can no longer drift.

## Alternatives rejected

- **Guard inside `has_palindrome()` only and let the hit-table function
  raise.** That is the defect: the two exported entry points to one Biostrings
  scan answered a fractional argument differently.
- **Move the guard into a helper both call.** One line at each exported entry
  is smaller and keeps each function's own contract visible; a shared wrapper
  would also have to carry three argument names.

## Consequences

- A fractional, negative, `Inf`, vector or non-numeric count is refused by name
  before `Biostrings::findPalindromes()`.
- Whole calls return the same hit tables as before.
- `min_armlength` still reaches Biostrings' own `>= 2` message for `1L`, the
  same as `has_palindrome()`; that bound belongs to Biostrings, not to the
  count rule.

## Verification

- `tests/testthat/test-palindromes.R`, "find_palindromes names the same three
  counts as has_palindrome": `max_mismatch` `1.9`/`-1`/`NA`/`Inf`/vector/`"0"`,
  `min_armlength` `0`/`4.5` and `max_looplength` `2.5` are refused by name; the
  whole call still returns 1 hit for `ACGTACGT` and 0 for `AACCGG`.
- Red before the fix: `1.9` returned the same table as `1L`; `0` and `Inf`
  raised the dependency messages above.
- Green after the fix: `dev/test.R "palindromes|design-metrics"` — 88 pass,
  0 fail.
