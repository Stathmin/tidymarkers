# 0044 — A count has one definition

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0042](0042-a-threshold-is-a-count-not-a-number-to-truncate.md)
  (a threshold is a count, not a number to truncate),
  [ADR 0040](0040-a-validator-enforces-the-unit-the-contract-names.md) (a
  validator enforces the unit the contract names)

## Context

Four argument gates ask the same question - "is this a single whole number in
range?" - and each spelled it out by hand:

| gate | file | condition |
|---|---|---|
| `max_markers` | `R/marker_one_pot.R` | `is.numeric`, length 1, `!is.na`, `>= 1`, `%% 1 == 0` |
| `spacer_len` | `R/marker_specificity_policy.R` | non-number first coerced to `NA_real_`, then length, finite, range, rounding |
| `min_lines` | `R/marker_informativeness.R` | numeric, length 1, finite, range, rounding (ADR 0042) |
| `min_evidence` | `R/marker_ingest.R` | numeric, length 1, `!is.na`, `>= 0`, `!= round(...)` |

They drifted in two ways that a single definition would not allow:

- `marker_one_pot(max_markers = Inf)`: `Inf %% 1` is `NaN`, `NaN != 0` is `NA`,
  so the condition evaluated to `NA` and `if` raised base R's
  `missing value where TRUE/FALSE needed` instead of
  `` `max_markers` must be one positive whole number ``.
- `marker_ingest_literature(min_evidence = Inf)`: `Inf >= 0` is `TRUE` and
  `Inf != round(Inf)` is `FALSE`, so the criterion passed, every row with a
  finite evidence count passed too, and the recorded `criteria$min_evidence`
  became `NA` through `as.integer()`. A bar that was never met read as "no bar".

`min_lines = Inf` and `spacer_len = Inf` already abort, but by two different
conditions; the rule had four homes and three shapes.

## Decision

One internal predicate, `is_whole_number(x, min = 1L, max = .Machine$integer.max)`
in `R/marker_validate.R`:

```r
is.numeric(x) && length(x) == 1L && is.finite(x) &&
  x >= min && x <= max && x == round(x)
```

`is.numeric()` first and `&&` short-circuit, so `NA`, `NaN`, `Inf`, a
character, a logical and a vector all answer `FALSE` and the caller reaches its
own message. Each gate keeps its wording and its own `as.integer()` after the
check; only the boolean is shared. `min = 0L` covers a non-negative count.

## Alternatives rejected

- **A helper that also aborts and casts** (`as_count(x, arg, message)`). It
  would force four different messages through one template ("positive whole
  number of nucleotides", "not `.val`"), so the shared part would shrink to the
  wording nobody agrees on.
- **Only fix `marker_one_pot` and `marker_ingest` in place.** That leaves four
  copies of the same condition; the next gate added would be a fifth.
- **Accept `Inf` as "no upper bar".** A criterion nobody can meet is not a
  criterion; a caller who wants no bar passes the largest count, and the
  recorded `criteria` then says what was applied.

## Consequences

- `max_markers = Inf` and `min_evidence = Inf` now abort with the package
  message; `min_evidence = Inf` no longer records `NA` as its criterion.
- `min_lines`, `spacer_len`, `max_markers` and `min_evidence` accept and refuse
  exactly the same values they did, except `Inf` for the last two.
- New file `R/marker_validate.R` (about 30 lines, one function).

## Verification

- `tests/testthat/test-validate.R`: the predicate's truth table - `1L`, `4`,
  `.Machine$integer.max` and `0` at `min = 0L` pass; `0`, `4.9`, a two-element
  vector, `numeric(0)`, `NA`, `NaN`, `±Inf`, `"4"`, `TRUE` and
  `.Machine$integer.max + 1` fail.
- `tests/testthat/test-one-pot.R`: `max_markers = Inf` gives "positive whole
  number" (was base R's `missing value where TRUE/FALSE needed`).
- `tests/testthat/test-ingest.R`: `min_evidence = Inf` gives "min_evidence".
- `Rscript dev/test.R "validate|one-pot|ingest|specificity|informativeness"` —
  557 pass, 0 fail, 2 skips (crisprVerse and the external cross-check, both
  pre-existing).
