# ADR 0068 — A position column is whole in its own spelling, and a selector keeps it whole

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast), ADR 0045
  amendment (a character column must not reach `round()`), ADR 0053 (a sample
  selector names a row, it does not round)

## Context

`bad_position()` (ADR 0045) judged a numeric column on its value and a character
column only through `as.integer()`:

```r
frac <- if (is.numeric(x)) !is.finite(x) | x != round(x) else rep(FALSE, length(x))
is.na(x) | is.na(iv) | iv < min | frac
```

`as.integer()` on text truncates or re-parses, so the character branch could
never see a fraction:

```r
as.integer("250.9")  #> 250
as.integer("1e3")    #> 1000
```

Measured on the branch, three exported paths still accepted a position that is
not a whole number:

| path | input | before the fix |
|---|---|---|
| `read_variant_genotypes()` | VCF `POS = 2801.5` | row at 2801, no message |
| `as_variant_row(x, pos = 250.9)` | selector | the variant at 250 |
| `as_variant_row(x)` | row with `pos = 250.9` | returned with `pos = 250` |

The streaming genotype reader checked the cast (`is.na(pos) || pos < 1L`), so
the record passed; `as_variant_row()` cast the selector with `as.integer()`
before matching and cast the chosen row after `is.na()`, so a fraction could
only become `NA`, never be seen. `pos = Inf` in the selector raised base R's
`NAs introduced by coercion to integer range` warning instead of a package
message.

This is the residual of the family ADR 0045 opened: the helper was used as the
single rule, but the rule's character branch did not answer for a spelling that
`as.integer()` silently rewrites.

## Decision

- `bad_position()` judges a character value by its spelling, not by what
  `as.integer()` makes of it: `^ *[+-]?[0-9]+ *$`. `"250"`, `"+250"` and
  `" 250 "` stay accepted (CSV input); `"250.9"`, `"2.0"`, `"1e3"` and `""` are
  refused. A non-numeric, non-character value keeps its previous answer.
- `as_variant_row()` checks `pos =` with `bad_position()` before the selector
  matches, and checks the chosen row's own `pos` with `bad_position()` before
  `as.integer()`; the existing messages and the "exactly one variant" rule are
  unchanged.
- `read_variant_genotypes()` checks the raw `POS` field with `bad_position()`
  before the cast and keeps its record-numbered message.

## Alternatives rejected

- **Compare the cast with the parsed double.** `as.numeric("1e3") == 1000` is
  true, so the exponent spelling would pass and a later reader would disagree
  about what the text meant; the spelling is the only stable answer.
- **Trim and compare `as.character(as.integer(x))` with the input.** It refuses
  `"0250"` and `"+250"`, which are ordinary CSV spellings, for no gain.
- **Fix `as_variant_row()` only.** The VCF reader and every other caller of the
  helper would keep the hole; the root is the shared predicate.

## Consequences

- A fractional or exponent-spelled position is refused by name in targets,
  genotype panels, haplotype tables, annotation readers, `as_variant_row()` and
  the streaming genotype reader.
- Every position that is spelled as an integer behaves exactly as before, so no
  accepted input changes its statistic or verdict.
- `as_variant_row(pos = Inf)` now aborts with the package message instead of a
  base R warning followed by "0 matched".

## Verification

- `tests/testthat/test-validate.R`: `bad_position()` truth table — `"250"`,
  `"+250"`, `" 250 "` pass; `"250.9"`, `"1e3"`, `"2.0"`, `""` refused.
- `tests/testthat/test-io.R`: VCF `POS = 2801.5` refused by
  `read_variant_genotypes()`; `as_variant_row(pos = 10.9)` and `(pos = Inf)`
  refused with the selector message; a row at `250.9` refused with the
  "positive integer `pos`" message.
- Focused run `Rscript dev/test.R "validate|^io$"` — 251 pass, 0 fail;
  shared-caller run `Rscript dev/test.R "plan|informativeness|haplotype|annotate|^design|target|pool|genome"`
  — 616 pass, 0 fail.
