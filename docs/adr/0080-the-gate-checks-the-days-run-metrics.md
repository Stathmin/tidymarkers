# ADR 0080 — The gate checks the day's run metrics

- Status: accepted
- Date: 2026-11-28
- Builds on: ADR 0011 (run metrics: one JSONL per run, numbers and ids only)

## Context

`dev/metrics.R --check` validates one day's run rows against the metrics
contract: a row must carry an id-shaped `object_id` of at most 120 characters,
a known stage and status, and integer durations and counts. Nothing in
`dev/check.R` ran it. On 2026-11-28 three repair rows in
`data/metrics/2026-11-28-repair.jsonl` carried prose `object_id`s of 128–153
characters, so `--check` refused the file, while the committed report still
described 39 rows (the file had grown to 60). The invalid data and the stale
report were found by hand.

The cast gate (ADR 0079) made the same point for code: a rule that only runs
when someone remembers it is not enforced. The metrics contract has the same
shape — it is checkable in well under a second with base R and `jsonlite`.

## Decision

- `dev/check.R` gains a `run metrics` step after `position casts`. It runs
  `Rscript dev/metrics.R --check` with no argument, which defaults to the
  current date, in every mode including `--quick`.
- The gate checks the current day only. An old malformed day does not block an
  unrelated change; repairing it is its own commit, as the three ids and six
  escaped journal rows were repaired today.

## Rejected alternatives

- Check every committed day: one historical bad row would fail every branch
  until someone edited old data, which is not the change under review.
- Regenerate the report in the gate: the report also writes a file and compares
  against `baseline.json`; a statistic is not a gate verdict, and a green gate
  must not leave a diff.
- Keep the check manual: measured cost is the invalid file above — it passed
  every green gate that day.

## Consequences

- A metrics row with a bad id, stage, status or count fails the gate on the day
  it lands, before it reaches the report.
- A day's file stays as committed once the day ends; a repair commit is still
  the way to clean old rows, and `--report` names malformed rows it drops.
- The step costs about 0.2 s and adds no dependency beyond the ones
  `dev/metrics.R` already uses.

## Verification

`Rscript dev/metrics.R --check` — `ok metrics 2026-11-28: 1 file(s), 60
row(s)`. `Rscript dev/check.R --quick` runs `run metrics` green. Appending a
row whose `object_id` is 121 characters makes the step fail with the id rule.
