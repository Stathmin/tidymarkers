# 0013 — The metrics report is not a gate, and a mixed-sample p95 is not a regression

- Status: proposed
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0011](0011-release-artefact-and-mcp-contract.md) (the daily
  metrics contract) and §0.1 of
  [`docs/plans/2026-11-27-integration-3day.md`](../plans/2026-11-27-integration-3day.md)

## Context

The three-day plan makes the metrics file the daily gate, not a dashboard.
§0.1 gate 1 says `p95(stage) <= p95_baseline(stage)` for every stage with
`n >= 20`, and failure criterion 9 says a stage p95 above baseline fails the
day. §8 repeats it for the programme.

The implementation does not enforce it. On 2026-11-27:

```sh
$ Rscript dev/metrics.R --report 2026-11-27; echo $?
0
$ Rscript dev/metrics.R --check 2026-11-27; echo $?
0
```

The report it writes says:

| stage   | n  | p95 ms | delta p95 |
| ------- | -- | ------ | --------- |
| ingest  | 40 | 30     | 900.0     |
| design  | 56 | 13646  | 27.3      |
| checks  | 30 | 144    | 1340.0    |
| mcp     | 27 | 293    | 946.4     |

`--check` validates the JSONL schema and rejects content; `--report` prints the
deltas and exits 0 whatever they are. The only green/red signal in the daily
gate is therefore the schema, not the regression the plan named.

The second half of the defect is that the printed comparison is not
like-for-like. A stage row pools every route that reached that stage. The
baseline (`baseline.json`, run `baseline-15a3c41` on `main`) has `ingest`
`n = 21`, `p95 = 3 ms`; the current file has `ingest` `n = 40`, `p95 = 30 ms`,
because the scenario runs added rows that the baseline never sampled. The
`checks` stage grows the same way (21 -> 30). Neither number is a slowdown of
one call: the sample composition changed. Enforcing `p95 <= p95_baseline` on
that pair would fail the gate for the wrong reason, and a tolerance would hide
a real one.

## Decision

Split the contract into the two things it conflates.

1. **Schema and content stay a gate.** `dev/metrics.R --check <date>` keeps
   exiting 1 on an unknown key, `prompt`, `content`, `sequence`, a target
   label or a path. That part is already enforced.
2. **A p95 comparison is enforced only on a comparable sample.** Two samples
   are comparable when they cover the same run set (the plan's per-route run
   ids) and both have `n >= 20` for that `(stage, route)`. When they are not
   comparable the gate reports `not_run` and names both sample sizes; `not_run`
   is not green and is not silent.
3. **The stage table stays informational.** `metrics-<date>.md` keeps printing
   the per-stage delta, but the delta is labelled as a report statistic, not a
   verdict. The gate verdict lives in its own `## Gate` block with one row per
   `(stage, route)` and the reason for every `not_run`.
4. **A comparable regression fails.** `--report` exits 1 when a comparable
   `(stage, route)` has `p95 > p95_baseline` or a billed `cost_per_route`
   above a zero baseline. A regression cannot be downgraded to a note.

## Alternatives

- **Set a tolerance (for example 2x) and fail anyway.** Rejected: the current
  deltas are 9-14x from composition, so a tolerance large enough to pass them
  would also pass a real 2x regression.
- **Re-baseline on every run.** Rejected: the baseline then always equals the
  current sample and the gate is green by construction.
- **Drop the p95 gate and keep only cost.** Rejected: §0.1 gate 1 and failure
  criterion 9 name latency; removing it would be a scope cut, not a decision.

## Consequences

- The report stops implying a verdict it does not compute. A reader sees
  `delta p95` and, separately, `gate: not_run (n 40 vs baseline 21)`.
- The gate becomes enforceable only after a repeated run of the same suite;
  until then it is honestly `not_run`, with the sample sizes in the file.
- `baseline.json` needs the run ids per stage so "same run set" is checkable
  rather than guessed from `n`. That is the first implementation step, and it
  changes the baseline schema once, additively.
- `dev/metrics.R --report` gains an exit status. The day gate must keep
  calling it; a green schema check alone is no longer enough.

## Gates

1. `dev/metrics.R --report <date>` writes `## Gate` with one row per
   `(stage, route)`: `ok`, `fail` or `not_run` plus the two sample sizes.
2. A fixture pair with equal run sets and a 2x p95 exits 1; the same pair
   within baseline exits 0; unequal `n` exits 0 with `not_run`, never `ok`.
3. `dev/metrics.R --check <date>` still rejects unknown keys and content.
4. `docs/plans/2026-11-27-integration-3day.md` §0.1 gate 1 names the
   comparable-sample rule; README does not promise a stage-level gate.

## Revisit trigger

A repeated run of the same J1-J7 suite gives two comparable samples and the
first real `ok`/`fail` verdict. If the comparable samples still disagree on
`n`, the run ids are not stable enough and the baseline key moves to the
scenario id.
