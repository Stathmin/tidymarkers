# 0058 — A repeat scan names whole units and copies

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0052](0052-a-palindrome-count-is-whole-before-the-scan.md)
  (a palindrome count is whole before the scan), [ADR 0051](0051-a-nucleotide-length-is-whole-before-the-design-grid.md)
  (a nucleotide length is whole before the design grid)

## Context

`repeat_regex()` in `R/marker_thermo.R` turns the repeat arguments into a
back-reference pattern:

```r
sprintf("([ACGT]{%d})\\1{%d,}", as.integer(unit), as.integer(min_copies) - 1L)
```

`as.integer()` truncates both numbers before the scan. Measured at HEAD on
`"ACACAC"`:

```r
find_repeats("ACACAC", min_copies = 1.9)
#> [[1]] "AC" "ACA" "CAC" "ACAC"   # a 2-copy hit at a 1.9 request
find_repeats("ACACAC", min_copies = 2.9)
#> [[1]] "AC"                      # silently the 2-copy hit
find_repeats("ACACAC", min_copies = 0.5)
#> Error in {min,max} interval. (U_REGEX_BAD_INTERVAL,
#>   context=([ACGT]{2})\1{-1,})
find_repeats("ACACAC", units = 2.9)
#> [[1]] "AC"                      # scanned unit 2
has_repeat("ACAC", min_copies = 1.9)
#> TRUE                            # 2 copies accepted
```

The same truncation reaches `design_ssr_marker()`, because
`ssr_span_at()` calls `repeat_regex()`. With `snp_fixture()` (600 bp,
seed 2, SNP at 300, `flank = 150L`) it returned a marker whose repeat unit
and copy count were the integer call's, so a fractional request silently
produced a marker for a different repeat.

`repeat_regex()` is shared by `find_repeats()`, `has_repeat()`, `has_run()`,
`ssr_span_at()` (and so `design_ssr_marker()`), `repeat_spans()` (the plot
layer) and `marker_screen()`. One guard at the shared builder reaches all of
them, the same shape as ADR 0051 and ADR 0056.

## Decision

`repeat_regex()` names both counts before it builds the pattern:

```r
unit <- check_nucleotide_length(unit, "units", multi = TRUE)
min_copies <- check_scan_radius(min_copies, "min_copies", min = 1L)
```

`check_scan_radius()` already exists for `max_mismatches` and `max_lengths`
(ADR 0042/0052); a copy count is one positive whole number. The message names
the caller's argument (`units`, `min_copies`), which is also the public name in
`find_repeats()`, `has_repeat()` and `design_ssr_marker()`.

## Alternatives rejected

- **Round with `round()` or `as.integer()`.** The scan would still use a
  number the caller did not name, and the copy count would still change
  silently.
- **Check in each exported function.** That duplicates the guard at four or
  more entries and leaves the next internal caller free to reintroduce the
  truncation.
- **Guard only `min_copies`.** A fractional `units` scans a different repeat
  unit, which is the same class of wrong answer; both numbers form one regex.

## Consequences

- A fractional, zero, `Inf`, vector or non-numeric `units` or `min_copies` is
  refused with a package message before any regex is built; `0.5` no longer
  reaches a raw stringi interval error.
- Whole values keep the same patterns and the same hits.
- `find_repeats(units = integer(0))` now aborts instead of scanning nothing.
  `has_repeat()` already short-circuits an empty `units` to `FALSE` and keeps
  that contract, so the two functions differ on an empty set by design.
- The plot layer (`repeat_spans()`) inherits the refusal instead of drawing a
  repeat at a truncated unit.

## Verification

- `tests/testthat/test-ssr.R`, "a repeat scan names whole units and copies":
  `units = 2.9`, `min_copies = 2.9`, `min_copies = 0.5`, a two-element
  `min_copies` and `Inf` are refused with "positive whole numbers of
  nucleotides" / "positive whole number"; whole calls still return the same
  hits; `design_ssr_marker()` refuses both fractions by name.
- Red before the fix: `dev/test.R ssr` — 3 failures (the two fractional
  requests did not throw and `0.5` raised the stringi error).
- Green after the fix: `dev/test.R "ssr|thermo|palindromes"` — 167 pass,
  0 fail.
