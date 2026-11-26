# 0003 — DuckDB/`duckplyr` is deferred until an aggregation is measured

- Status: accepted
- Date: 2026-11-26
- Supersedes: none

## Context

The proposal offered DuckDB with `duckplyr` as an embedded analytics layer for
cross-project aggregations ("all marker outcomes", "quality trends across
releases"). In this repository the aggregation question is real but not yet
answerable: the pooled corpus that would be aggregated (`marker_trace_schema()`
records across designs, ZY5) is externally blocked, and the population-scale
output that would justify an out-of-core engine does not exist yet (`ZX2`, the
wheat fixture). Today `marker_artefact()` and `criterion_spec()` return tables
that fit comfortably in memory; `dplyr`/`vctrs` already do the work.

The dependency policy in `ROADMAP.md` ZA4 refuses `data.table` without a second
adapter. `duckplyr` is the same kind of addition — a second query engine next to
`dplyr` — and it pulls `duckdb`.

Neither `duckdb` nor `duckplyr` is installed in the R 4.6.1 / Bioconductor 3.23
library this repository runs against; `targets` and `data.table` are.

## Decision

Deferred. Not installed now, and never in `Imports` or `Suggests` as a runtime
dependency. If adopted it is dev-only (a script under `dev/`, `.Rbuildignore`d),
the same boundary as `targets` (ADR 0002).

Adoption requires a measurement, not a hunch: the same aggregation run through
`duckplyr` and `dplyr`, with wall time and peak memory recorded in `results/`,
and the two answers compared on the same input.

## Consequences

- No new dependency; `dplyr` stays the only aggregation path.
- A future out-of-core need is not pre-optimised: the pooled corpus will be
  built without a DuckDB shim in the way.
- Because `duckdb` is not in the library, any adoption starts with an install
  decision, which is an operator action on a pinned R 4.6.1 tree — not a
  silent `install.packages()` in a pipeline.

## Alternatives

- Adopt now "for the aggregate views". Rejected: there is no table to aggregate
  yet, so the benchmark that would justify it cannot be run.
- Adopt `data.table` instead. Rejected by the same ZA4 rule: a second engine
  without a second consumer.

## Revisit trigger

Both must hold:

1. Population-scale or pooled output exists (ZX2 / ZY5), and
2. one aggregation either exceeds available RAM or is run by hand more than
   once a week.

Until then this ADR stands. Re-check `requireNamespace("duckdb")` before any
plan, and keep the engine under `dev/`.

## References

- [`ROADMAP.md`](../../ROADMAP.md) — ZA4 (dependency policy), ZX2, ZY5.
- [`R/marker_trace_schema.R`](../../R/marker_trace_schema.R) — the pooled record.
- [`R/marker_artefact.R`](../../R/marker_artefact.R) — what a design exports today.
