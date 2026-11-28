# 0064 — A terminal window is whole before it is cut

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0047](0047-a-search-radius-is-a-whole-number-before-the-cast.md)
  (a search radius is a whole number before the cast)

## Context

`terminal_dg()` in `R/marker_qc.R` reports the duplex free energy of a primer's
3'-terminal `n`-mer, the "GC clamp" diagnostic. It cut the window with a raw
`n`:

```r
mer <- substr(s, nchar(s) - n + 1L, nchar(s))
```

`substr()` truncates a fractional start, so the annealed molecule was not the
window the caller named:

```r
p <- "ACGTACGTACGTACGTACGT"
terminal_dg(p, n = 5L)    # -2.812  (TACGT)
terminal_dg(p, n = 5.9)   # -3.35   (GTACGT, a 6-mer)
terminal_dg(p, n = 6.9)   # -5.178  (TGTACGT, a 7-mer)
terminal_dg(p, n = 0.5)   # 0       (a zero-width mer)
terminal_dg(p, n = Inf)   # NA + "NAs introduced by coercion to integer range"
```

The value is reported next to `tm` and `hairpin_dg` in the design diagnostics,
so a 6-mer dG labelled as a 5-mer is a wrong number in a user-facing table, not
a cosmetic difference.

## Decision

`terminal_dg()` names `n` before it cuts the window:

```r
n <- check_scan_radius(n, "n", min = 1L)
```

The helper already exists (ADR 0047); `n` is the number of terminal bases, so
one is the floor. A primer shorter than the window keeps its documented
`NA_real_`.

## Alternatives rejected

- **`as.integer(n)` before `substr()`.** That is the defect: `5.9` becomes 5,
  but `0.5` still becomes a zero-width mer and `Inf` still becomes `NA`.
- **Refuse only inside `oligo_dg_vec()`.** It sees the mer, not the argument
  name, and the zero-width mer is built before it is called.
- **Let `substr()` warn.** A warning is not a refusal, and the first call in a
  batch still returns a plausible number.

## Consequences

- A fractional, zero, negative, `Inf`, `NA`, vector or non-numeric `n` is
  refused by name before `substr()`.
- `n = 5L` and `n = 5` return the same number as before; a primer shorter than
  the window still answers `NA_real_`.

## Verification

- `tests/testthat/test-design-metrics.R`, "terminal_dg names its terminal
  window before it is cut": `5.9`, `6.9`, `0.5`, `0`, `-1`, `Inf`, `NA`, a
  vector and `"5"` are refused by name; `n = 6L` equals `n = 6`, differs from
  the 5-mer, and `"ACGT"` at `n = 5L` is still `NA`.
- Red before the fix: `5.9`/`6.9` returned the 6-/7-mer dG, `0.5` returned 0,
  and `Inf` returned `NA` with a coercion warning.
- Green after the fix: `dev/test.R "palindromes|design-metrics"` — 88 pass,
  0 fail.
