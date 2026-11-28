# ADR 0086 — A frequency weighs nothing until it is a proportion

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0038 (a count is admitted only in its domain), ADR 0045 (a
  position is a whole number before the cast), ADR 0068 (a position column is
  whole in its spelling)

## Context

A `pop_freq` is the share of the population that carries a variant, and a
haplotype's `freq` is the share of the population that carries that haplotype.
Both are proportions: one finite number in `[0, 1]`. Every reader that accepts
them cast the column and used the value without judging its domain. Measured on
the tree before this decision:

| call | returned |
|---|---|
| `marker_chunk(df)` with `pop_freq = 5` | the weight `5`, passed through |
| `marker_chunk(df)` with `pop_freq = -0.2` | the weight `-0.2` |
| `apply_maf(v, maf = <pos = 5.5, maf = 0.02>)` | no override applied, no message |
| `apply_maf(v, maf = <pos = 5, maf = 5>)` | `pop_freq = 5` |
| `normalise_hap_freq(c(1.5, 0, ...), nms)` | `1.5` for the first haplotype |
| `normalise_hap_freq(tibble(haplotype, freq = "5%"), nms)` | share `0` for that haplotype, silently |
| `variant_load_in(iv, c(10, 20, 30), c(0.5, NA, 2))` | `NA` |

Three of these change the answer, not just the report:

- `variant_load_in()` sums the weights of the variants under a primer, so a
  percent-scale table counts one variant as five, and a negative weight makes a
  site *better* than a clean one — the primer ranking flips.
- An `NA` weight returned `NA` for the whole load, and the dropout path
  (`zone_dropout()` → `dropout_pairs()`) dropped the variant outright, so an
  unknown frequency read as **zero cost**. Measured: the same footprint reported
  `as_foot_dropout = 0` with `pop_freq = NA` and `0.009` with `pop_freq = 1`.
- A fractional `pos` in a MAF table matched no variant, so the override was
  silently lost while the design continued with the alignment-implied weight. An
  unparseable share became `0`, removing that haplotype from the population.

The load documentation already promised the opposite of the `NA` behaviour:
"absent weights every variant counts 1". A missing value in a bound table is the
same absence.

## Decision

1. `bad_proportion()` is the vectorised domain rule, next to `bad_position()`:
   a present value must be one finite number in `[0, 1]`. A character value is
   judged by its spelling (`"0.05"` passes, `"5%"` does not); a factor, logical,
   complex or list column answers `TRUE` (its `as.numeric()` is a code, a `0/1`
   or a lost imaginary part). `NaN` is refused; `NA` answers `FALSE` — an
   absence is resolved by the caller, not by the predicate.
2. `check_variant_weights()` applies that rule where weights enter scoring
   (`marker_variant_alleles()` builds `vweight` once per assay) and maps `NA`
   to `1`, the documented whole-variant weight. `variant_load_in()` and
   `zone_dropout()` keep the same fallback for a caller that hands them a raw
   vector, so an unknown weight can no longer read as zero cost.
3. `apply_maf()` judges `pos` with `bad_position()` before the join and the
   `maf`/`freq` column with `bad_proportion()`. `NA` in the column still means
   "no override for this row" — the documented meaning — and the column is cast
   once, after the check.
4. `normalise_hap_freq()` refuses an `NA`, non-finite or out-of-range share
   before the name lookup. A haplotype *absent from the table* still carries
   share `0`, as documented; an unparseable share for a named haplotype no
   longer disappears into that same `0`.

## Rejected alternatives

- Clamp to `[0, 1]`: `5` becomes `1` and `-0.2` becomes `0`, so a percent-scale
  table still scores and a caller never learns the column is in the wrong unit.
- Refuse `NA` everywhere: a variant table often binds rows that carry no
  `pop_freq` (the package's own `freq_fixture()` does), and there the value is
  genuinely unknown. Refusing it would make the documented optional column
  mandatory per row.
- Keep `NA` as "no evidence": for a *weight* that is not the same as zero
  evidence. A variant of unknown frequency may be common, so counting it once is
  the conservative reading; the load's own documentation says so.
- Judge only in the readers: `assemble_kasp_assay(variants =)` and the panel
  accept a table directly, so the reader is not the only trust boundary.

## Consequences

- A percent-scale `pop_freq`, a negative share and a misspelled `freq` are
  refused by column and row instead of skewing the ranking.
- A fractional `pos` in a MAF table is refused instead of silently not matching.
- An unknown weight counts the variant once, so a partially weighted table
  degrades to the plain count instead of `NA` or zero cost.
- Tables that were silently wrong now stop the run; no valid table changes its
  answer (the test suite's frequency fixtures pass unchanged).

## Verification

`Rscript dev/test.R 'alignment|frequencies'` — 74 pass, 0 fail, including the
four new refusals and the `NA`-weight load equality.
`Rscript dev/test.R 'selection|chunk|haplotypes|panel|lax|informativeness|plots'`
— 959 pass, 0 fail: no design path passes a value outside the domain.
