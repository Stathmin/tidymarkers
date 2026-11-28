# 0046 — A design key names every input the design reads

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0042](0042-a-threshold-is-a-count-not-a-number-to-truncate.md)
  (a threshold is a count, not a number to truncate), [ADR 0045](0045-a-position-is-a-whole-number-before-the-cast.md)
  (validate the value the caller gave, not a later cast)

## Context

`marker_solve()` de-duplicates its targets and memoises each design by
`target_key()`:

```r
target_key <- function(tg) paste(tg$chrom, tg$pos, tg$ref, tg$alt, sep = "|")
```

But `target_snp()` forwards two further target fields to the designer, and each
one changes the design:

- `strand` — `design_common_primers()` reverse-complements a minus-strand locus
  before it places an allele-specific primer.
- `alt2` — a second ALT makes the site multi-allelic and changes the ALT primer
  to an IUPAC code.

So two rows with the same `chrom/pos/ref/alt` but a different strand or a
second ALT were one key. Measured on the branch before the fix:

```r
p <- marker_plan(genome) + marker_tech(schemes = "snp") +
  marker_targets(chrom = "Chr1", pos = 20L, ref = "T", alt = "G", strand = "-") +
  marker_targets(chrom = "Chr1", pos = 20L, ref = "T", alt = "G", strand = "+")
target_key(p$targets)
#> "Chr1|20|T|G" "Chr1|20|T|G"
```

A two-row solve returned two report rows but one design, and the cache entry
for the second row was looked up under the first row's key. Single-row solves
of the same locus with `+` and `-` are not identical (`all.equal()` reports the
complemented `snp`), so the collapsed row was answered by a design for the
other orientation. `alt2 = "A"` likewise produces a different design.

## Decision

`target_key()` names every field `target_snp()` forwards:

```r
target_key <- function(tg) {
  alt2 <- if ("alt2" %in% names(tg)) tg$alt2 else rep(NA_character_, nrow(tg))
  alt2[is.na(alt2)] <- ""
  strand <- if ("strand" %in% names(tg)) tg$strand else rep("+", nrow(tg))
  paste(tg$chrom, tg$pos, tg$ref, tg$alt, alt2, strand, sep = "|")
}
```

`alt2 = NA` is the empty field and a missing `strand` is `"+"`, matching
`normalise_targets()`'s defaults, so tables that predate either column keep
their old key semantics. Two truly identical rows still collapse to one design;
the key is context, not identity.

## Alternatives rejected

- **Drop `alt2`/`strand` from `target_snp()`.** They are documented inputs of
  `marker_targets()` (`strand` orients the alleles, `alt2` names a
  multi-allelic site); dropping them would refuse designs the package already
  promises.
- **Hash the whole target row.** `label` and `type` are reporting fields; the
  same locus under two labels would then be designed twice for no gain.
- **Add only `strand`.** `alt2` changes the ALT primer just as directly, and a
  key that names one design input but not the other leaves the same trap.

## Consequences

- Rows that differ only in `strand` or `alt2` get their own design and their
  own cache entry; a repeated identical row is still solved once.
- Cache keys change shape. Nothing persisted depends on them: `marker_replay()`
  re-solves from the record with `cache = NULL`, so old artefacts still replay.
- `artefact_lookup()` uses the same key and therefore finds the design of the
  row it is asked about.

## Verification

- `tests/testthat/test-plan.R`: "a target's design key includes alt2 and strand"
  seeds the cache per row with a tagged marker and requires each report row to
  carry its own tag for `strand = c("-", "+")` and for `alt2 = "A"` vs `NA`;
  the same row twice still carries one object twice.
- `Rscript dev/test.R plan` — 215 pass, 0 fail.
