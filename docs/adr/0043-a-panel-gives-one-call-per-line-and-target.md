# 0043 — A panel gives one call per line and target

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0041](0041-a-genotype-is-its-allele-multiset.md) (a genotype
  is its allele multiset),
  [ADR 0020](0020-a-gt-that-names-a-missing-allele-is-not-a-call.md) (a GT
  that names a missing allele is not a call)

## Context

`marker_informativeness(genotypes = )` scores a target against the lines that
carry a call there. `genotypes_for_target()` filtered the panel to the target
and then collapsed duplicate lines with

```r
d[!duplicated(d$line), , drop = FALSE]
```

which keeps the **first** row in file order. When one line arrived twice with
different calls - two panels pasted together, a merge that kept both sides -
the answer was whichever row happened to come first. Measured on the branch:

```r
panel <- tibble::tibble(
  line = c("L1", "L1", "L2", "L3", "L4", "L5"), pos = 250L,
  ref = "A", alt = "G", genotype = c(2, 0, 0, 0, 0, 0)
)
marker_informativeness(tg, genotypes = panel)                 #> alt_freq 0.2, pic 0.2688, "ok"
marker_informativeness(tg, genotypes = panel[c(2, 1, 3:6), ]) #> alt_freq 0,   pic 0,      "monomorphic"
```

The same set of calls gave two verdicts. The `NA` row had the milder version of
the same problem: `NA` beside a call dropped the line when `NA` came first and
counted it when the call came first.

## Decision

A genotype panel is a set of calls keyed by `(line, target)`:

- two rows for one line at one target with **the same** dosage are one call;
- a missing call beside a present one is that line's call, not a second call;
- two **different** dosages are a caller error: abort, naming the line, the
  target and both dosages.

Duplicate collapsing keeps a present call over a missing one, so the `NA` case
no longer depends on order.

## Alternatives rejected

- **Keep the first row and document "one row per line".** The input contract
  already said that; the code still answered instead of refusing, and the
  answer changed with the row order, which is the defect.
- **Average the two calls.** There is no dosage between `0` and `2` that the
  panel measured; inventing one would fabricate data.
- **Drop the line when calls conflict.** That silently shrinks the panel and
  can flip `insufficient_sample`; the caller must decide which panel is right.

## Consequences

- A conflicting panel aborts before any statistic is computed; identical
  duplicates and `NA`/call pairs keep their previous counts.
- `n_lines`, `alt_freq`, `pic` and `distinguishable_pairs` do not depend on row
  order for any panel that is accepted.
- No column or artefact field changes; haplotype input is untouched.

## Verification

- `tests/testthat/test-informativeness.R`, "a line with two different calls at
  one target is refused, not ordered": both row orders abort, the message names
  `Line 'L1' ... (dosages 0, 2)`.
- Same file, "the same call twice is one line, with or without a missing row":
  a duplicated ALT call counts as four lines at `alt_freq = 0.25`, and the
  `NA`-beside-call panel gives one `alt_freq`/`pic` in both row orders.
- Before the fix the two orders gave `pic` 0.2688 and 0; `Rscript dev/test.R
  informativeness` — 56 pass, 0 fail (was 48).
