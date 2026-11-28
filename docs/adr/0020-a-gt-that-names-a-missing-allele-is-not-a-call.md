# 0020 — A GT that names a missing allele is not a call

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)
  (one allele and `POS` rule for all three readers),
  [ADR 0019](0019-an-absent-allele-is-not-a-base.md) (an absent allele is
  absent before it is normalised)

## Context

`read_variant_genotypes()` turns one `FORMAT` field into a dosage and a depth
through `genotype_call(gt, ad)`. The documented contract says two different
things about a `GT` the reader cannot use:

- "A genotype naming an allele other than REF/ALT is `NA`, never silently
  folded to the first ALT."
- "Dosage comes from `GT`; when `GT` is missing, `AD` decides it (rounded
  `2 * alt / (ref + alt)`)."

The second sentence is a fallback for an absent call. The code implemented it
as a fallback for any `NA` dosage, so it also answered for a `GT` that was
present and unparseable. A record split from a multi-allelic site, or written
by a caller that renumbered alleles, can carry `GT = 2/2` next to
`REF = A`, `ALT = G`. Allele `2` does not exist in that record, so the call is
corrupt — but the reader returned the `AD` estimate:

| record | GT | AD | dosage before |
| --- | --- | --- | --- |
| biallelic | `2/2` | `6,4` | `1` |
| biallelic | `0/2` | `6,4` | `1` |
| biallelic | `0/1` | `6,4` | `1` (correct control) |

The first two rows contradict the first sentence, and the third row shows why
the bug is silent: the corrupt call and the real call are the same number. A
corrupt call was therefore replaced by a depth estimate without a count, a
message or a `NA`, and `marker_informativeness()` sums that dosage
downstream.

## Decision

The `AD` fallback applies only when `GT` is absent or a no-call. A `GT` that
is present and names an allele the record does not have is `NA`.

```r
no_call <- is.na(gt) || !nzchar(gt) || gt %in% c(".", "./.", ".|.")
dosage <- if (no_call) {
  NA_real_
} else {
  alleles <- strsplit(gt, "[/|]", perl = TRUE)[[1L]]
  a <- suppressWarnings(as.integer(alleles))
  a <- a[!is.na(a)]
  if (length(a) != 2L || any(a > 1L)) NA_real_ else sum(a)
}
if (no_call && !is.na(depth) && depth > 0) {
  dosage <- round(2 * alt_depth / depth)
}
```

- Depth is still reported for such a record. The read happened; only the call
  cannot be trusted.
- A no-call (`GT = .`, `./.`, `.|.`) keeps its pre-existing behaviour: `AD`
  fills the dosage when `AD` is present. The bug was the conflation of "no
  call" with "a call the record cannot support", not the fallback itself.
- A record with no `GT` column at all is unchanged: `AD` decides.

## Alternatives rejected

- **Drop the record and count it.** A corrupt `GT` still carries a usable
  depth, and the biallelic record itself may be fine. Dropping removes a read
  to punish its annotation.
- **Fold allele `2` onto the first ALT.** That is the silent reinterpretation
  the contract forbids; `2/2` would become a false hom-alt.
- **Keep the `NA`-dosage test as the trigger and only document AD.** The help
  page already promises `NA`; documenting the contradiction leaves the
  numeric answer wrong.
- **Treat a no-call as invalid too.** That would change every `./.` record
  with an `AD` column, a wider behaviour change than this defect needs, and
  `AD` is real evidence behind a missing call.

## Consequences

- A dosage is a call or `NA`. It is never a depth estimate standing in for a
  call the record contradicts.
- The `AD` fallback has one trigger: `GT` absent or no-call. Two sentences in
  the help page and one branch in the code now say the same thing.
- `marker_informativeness()` no longer receives a dosage derived from a
  corrupt genotype.

## Verification

- `tests/testthat/test-io.R`, "a GT naming a missing allele is NA, not a
  guess from AD": `2/2:6,4` and `0/2:6,4` give `NA`; `0/1:6,4`, `.:6,4` and
  `./.:6,4` give `1`; all five keep `depth = 10`. Before the change the first
  two returned `1` (probe `/tmp/opencode/probe_gt.R`).
- `Rscript dev/test.R "^io$"` — 163 pass, 0 fail.

## Amendment (2026-11-28) — a negative allele index is not a dosage

The decision above said "an allele the record does not have". The code checked
only the upper bound (`any(a > 1L)`), so a malformed `GT` whose allele is a
negative integer was summed and returned as a *negative* dosage:

| record | GT | AD | dosage before |
| --- | --- | --- | --- |
| biallelic | `-1/0` | `6,4` | `-1` |
| biallelic | `-1/-1` | `6,4` | `-2` |

No allele index below zero exists, so the same rule applies: a genotype whose
alleles are not both `0` or `1` is `NA`, and `AD` still does not answer for it
because the call is present. One bound, not two:

```r
if (length(a) != 2L || any(a < 0L | a > 1L)) NA_real_ else sum(a)
```

`marker_informativeness()` therefore never sums a negative dosage. Verified in
`tests/testthat/test-io.R` ("a GT naming a missing allele is NA, not a guess
from AD"): `-1/0:6,4` and `-1/-1:6,4` are `NA`, `depth = 10`;
`Rscript dev/test.R "^io$"` — 213 pass, 0 fail.
