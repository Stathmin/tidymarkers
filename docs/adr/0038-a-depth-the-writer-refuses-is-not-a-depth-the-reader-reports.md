# 0038 — A depth the writer refuses is not a depth the reader reports

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0018](0018-one-allele-and-position-rule-for-three-vcf-readers.md)
  (one allele and `POS` rule for the three readers) and
  [ADR 0028](0028-a-vcf-writer-emits-only-records-its-readers-admit.md) (a VCF
  writer emits only records its readers admit)

## Context

`read_variant_genotypes()` reads `FORMAT/AD` into `depth` and, when `GT` is
absent or no-call, into the dosage. The VCF specification types `AD` as
`Number=R, Type=Integer`; `write_variants_vcf()` enforces that on the way out:
a depth that is negative, fractional or outside the integer range is refused
(ADR 0028, "depth outside the non-negative integer range").

The reader accepted any two finite numbers:

```r
if (length(d) >= 2L && all(is.finite(d[1:2]))) {
  ref_depth <- d[[1L]]
  alt_depth <- d[[2L]]
}
```

Measured on the branch before the fix:

| `GT:AD` | dosage | depth |
| --- | --- | --- |
| `.:-5,10` | `4` | `5` |
| `0/1:5,-10` | `1` | `-5` |
| `0/1:6,4.5` | `1` | `10.5` |
| `.:1e10,1e10` | `1` | `2e10` |

The first row is the sharpest: a negative read count produced a dosage of `4`,
outside the documented `0`/`1`/`2`/`NA` domain. The other three rows report a
depth the package's own writer refuses. A table read from a file and then
written back cannot round-trip.

## Decision

The `AD` fallback uses only counts the writer would accept: both values finite,
non-negative, whole, at most `.Machine$integer.max`, and their sum within that
range. Otherwise `AD` is treated as absent:

```r
counts <- if (length(d) >= 2L) d[1:2] else numeric()
valid <- length(counts) == 2L && all(is.finite(counts)) &&
  all(counts >= 0 & counts == floor(counts) & counts <= .Machine$integer.max) &&
  sum(counts) <= .Machine$integer.max
```

- An invalid `AD` never changes a valid `GT`: the call still comes from `GT`,
  and only `depth` becomes `NA`.
- A no-call `GT` with an invalid `AD` stays `NA`/`NA`; the fallback has no
  evidence to answer from.
- No message is added. `AD` is annotation, and the readers already count
  dropped records for alleles; a malformed depth column has no separate count
  in the contract.

## Alternatives rejected

- **Clamp negative counts to zero.** That invents a read that was not there and
  makes `-5,10` look like a 10-read site.
- **Round fractional counts.** VCF says `Integer`; a fractional count is not a
  measurement to repair, and the writer would still refuse the result.
- **Refuse the whole record.** The allele and `GT` are still usable; only the
  depth annotation is corrupt. Dropping it would discard a call to punish a
  count.

## Consequences

- The dosage domain is `0`/`1`/`2`/`NA` for every input, not only for
  well-formed `AD`.
- A genotype table read by `read_variant_genotypes()` is accepted by
  `write_variants_vcf()` with the same depth value.
- The rule lives in `genotype_call()` only; `read_variant_depth()` returns raw
  `DP`/`AD` columns and is a separate raw view (its contract does not compute a
  dosage). *(Amended 2026-11-28: see below — the depth reader obeys the same
  rule, because its columns feed the confidence gate, not only a printout.)*

## Verification

- `tests/testthat/test-io.R`, "an AD the writer refuses is not a depth the
  reader reports": `-5,10`, `5,-10`, `6,4.5` and `1e10,1e10` give
  `depth = NA`; a valid `GT` keeps its dosage; `6,4` still gives depth `10`.
- Before the fix the same rows returned `5`, `-5`, `10.5` and `2e10`, and the
  first returned dosage `4` (probe `/tmp/opencode/probe_gt.R`).
- `write_variants_vcf()` refuses each of the four depths by the documented
  message, so the reader and the writer now name the same admissible set.
- `Rscript dev/test.R "^io$"` — 220 pass, 0 fail.

## Amendment 2026-11-28 (night shift): the depth reader obeys the same rule

The Consequences above called `read_variant_depth()` a separate raw view. That
is wrong on the path that matters: its `depth`/`alt_depth` columns are the input
to `variant_confidence()`, and the `confident` column it produces gates the
solver through `filter_targets_by_depth()`. A malformed `DP` is therefore not a
printout; it is a gate decision.

Measured on the branch before this amendment:

| field | read back | `write_variants_vcf()` |
| --- | --- | --- |
| `GT:AD:DP = 0/1:-5,10:-5` | `depth -5`, `alt_depth 10` | refuses depth `-5` |
| `GT:AD:DP = 0/1:6,4.5:10.5` | `depth 10.5`, `alt_depth 4.5` | refuses depth `10.5` |
| `GT:AD:DP = 0/1:1,1:1e10` | `depth 1e10`, `alt_depth 1` | refuses depth `1e10` |
| sites-only `INFO/DP=-5` | `depth -5` | refuses depth `-5` |

The third row is the sharpest: `1e10` exceeds the VCF `Integer` range and was
still `confident = TRUE`, because the depth term saturates.

### Decision

`vcf_count()` admits the writer's domain — finite, non-negative, whole, at most
`.Machine$integer.max` — and otherwise reports `NA`:

- `DP` (per-sample and sites-only `INFO/DP`) is `NA` when it is not a count;
- an `AD` field is admitted only when both of its first two counts are counts
  and their sum fits `integer`; otherwise `alt_depth` is `NA`. A valid `DP` is
  independent of a corrupt `AD`, and vice versa.
- No error is raised: the reader's contract already says an absent `DP`/`AD` is
  `NA`, and a malformed count is annotation corruption, not a reason to discard
  a usable allele row.

`genotype_call()` keeps its own pair check because it also constrains the sum
and the dosage.

### Verification

- `tests/testthat/test-io.R`, "read_variant_depth reports only counts the writer
  accepts (ADR 0038)": the three bad rows give `depth = NA` and, where the `AD`
  is malformed, `alt_depth = NA`; `6,4:10` gives `10`/`4`, `0,0:0` gives
  `0`/`0`, and `write_variants_vcf()` accepts every depth the reader reports.
  Sites-only `DP=-5`/`DP=7` give `NA`/`7`.
- `Rscript dev/test.R "^io$"` — 224 pass, 0 fail, 3 skip (was 220).
