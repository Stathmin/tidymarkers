# 0042 — A threshold is a count, not a number to truncate

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)
  (a position outside the positive integer range is refused, not coerced),
  [ADR 0038](0038-a-depth-the-writer-refuses-is-not-a-depth-the-reader-reports.md)
  (a depth the writer refuses is not a depth the reader reports)

## Context

`marker_informativeness(genotypes =, min_lines = )` gates a genotype panel on a
minimum sample size. The argument was cast before it was checked:

```r
min_lines <- as.integer(min_lines)[1L]
if (is.na(min_lines) || min_lines < 1L) {
  abort_msg("`min_lines` must be a positive whole number.")
}
```

`as.integer()` truncates towards zero, so the check never saw the fraction. With
a panel of four lines and `min_lines = 4.9`, the call asked for at least five
lines and the gate answered for four:

```r
panel <- tibble::tibble(
  line = paste0("L", 1:4), pos = 250L, ref = "A", alt = "G",
  genotype = c(2, 0, 0, 0)
)
marker_informativeness(tg, genotypes = panel, min_lines = 4.9)
#> status "ok", n_lines 4, pic 0.1875
```

The same cast accepted `"4"` as a character threshold and turned `Inf` into
`NA`, so the error message named a positive whole number while three different
non-numbers reached the gate. The truncation is the defect: a threshold that
the caller sets to five must not pass four.

## Decision

Validate the argument in its own domain, before any cast:

- numeric, length one, finite;
- at least `1` and at most `.Machine$integer.max`;
- equal to its own rounding (`4L` and `4` pass, `4.9` and `4.5` do not).

Only then cast to integer. The same rule already guards `spacer_len`
(`R/marker_specificity_policy.R`) and `max_markers` (`R/marker_one_pot.R`); the
wording follows them.

## Alternatives rejected

- **Round instead of truncate.** `4.9` would become 5, which is not what the
  caller wrote either; a threshold is a count, so a fraction is a caller error,
  not a value to normalise.
- **Accept a fraction and use `ceiling()`.** That hides a unit mistake in the
  caller and still answers a question nobody asked; the package's other numeric
  gates refuse instead of reinterpreting (ADR 0018, 0038).
- **Keep the cast and check `min_lines != round(min_lines)` after it.** The
  fraction is already gone at that point; the check must run on the input.

## Consequences

- `min_lines = 4.9`, `"4"`, `Inf` and `1e10` are refused with one message.
- `min_lines = 4L` and `4` are unchanged; no panel verdict moves.
- No output column, artefact field or replay verdict changes.

## Verification

- `tests/testthat/test-informativeness.R`, "a sample-size threshold is a whole
  number, not a truncated one": `4.9`, `"4"` and `Inf` abort; a four-line panel
  is `"ok"` at `min_lines = 4L` and `"insufficient_sample"` at `5L`.
- Before the fix the first call returned `status = "ok"`; after it the same call
  aborts. `Rscript dev/test.R informativeness` — 48 pass, 0 fail (was 43).
