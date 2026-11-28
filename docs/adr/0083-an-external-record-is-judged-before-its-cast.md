# ADR 0083 — An external record is judged before its cast

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0045 (a position is a whole number before the cast),
  ADR 0068 (a position column is whole in its spelling), ADR 0033 (a download
  declares its path before the request)

## Context

Two getters parse a cached Ensembl REST response and cast its coordinates
straight into the package's table:

```r
# ensembl_gene()
start = as.integer(j$start), end = as.integer(j$end), strand = as.integer(j$strand)
# ensembl_variants()
pos = as.integer(.data$start)
```

`dev/casts.yml` registered both with the reason "`start`/`end` come from the
Ensembl REST record; the cast canonicalises its integer fields". That is a
schema assumption, not a check: nothing in the code fails if the assumption
stops holding. Measured with a cache file written by hand:

| record | returned |
|---|---|
| `start: 268901.5` | `start = 268901` |
| `end: 271503.5` | `end = 271503` |
| `strand: 1.5` | `strand = 1` |
| `start: "268901.9"` (variant) | `pos = 268901` |

Each is silent. A truncated gene start moves the target window and every
primer designed in it; a truncated strand flips which side the arms sit on; a
truncated variant position names another locus. The record is written by a
remote service and read back through a cache file, so it is a trust boundary:
the caller of `ensembl_gene()` cannot see what the server sent, and a schema
change is exactly the event the `need_cols()` message already anticipates.

## Decision

1. `ensembl_gene()` judges `start` and `end` with `bad_position()` and the
   strand with `is_whole_number(..., min = -1L, max = 1L)` before any cast.
   The refusal names the record id and both values.
2. `ensembl_variants()` judges `start` with `bad_position()` before the
   `transmute()` cast and names the offending variant id.
3. The two `dev/casts.yml` entries are removed: the casts are guarded in place,
   and a registry entry that describes a check which does not exist is what let
   the hole stay invisible. The cast in `ensembl_variants()` reads `df$start`
   rather than `.data$start` so the gate can see the guard it names.
4. A whole record is accepted exactly as before, including a whole character
   spelling (ADR 0068): only a non-whole, `NA`, `Inf` or non-positive value is
   refused.

## Rejected alternatives

- Trust the schema and keep the registry entry: the entry has to be re-audited
  by hand forever, and the failure it hides is a silently moved marker.
- Cast through `round()`: `268901.5` becomes `268902`, which is a different
  wrong locus than `268901`. A position that is not whole is refused, not
  rounded (ADR 0045).
- Judge the parsed object in a shared reader: `fromJSON()` returns a list for
  the gene and a data frame for the variants, and a shared reader would have to
  know both shapes. The two call sites are three lines each.

## Consequences

- A cached record with a fractional or textual fraction coordinate fails at the
  getter with a message naming the record, instead of producing a table whose
  positions look valid.
- The ledger drops from 33 to 31 entries and the guarded count rises from 19 to
  21; `dev/check-casts.R` reports the new split.
- The JSON schema is still not type-checked as a whole (e.g. a missing
  `alleles` column is caught by `need_cols()`); this ADR covers the coordinates
  the package casts, which is where a silent wrong answer lives.

## Verification

`Rscript dev/test.R '^target$'` — 48 pass, 0 fail, including four new refusals
(gene `start`, gene `end`, gene `strand`, variant `start` as a number and as
text). `Rscript dev/check-casts.R` — `52 position/threshold cast(s), 21 guarded
in place, 31 registered with a reason; self-test passed`.
