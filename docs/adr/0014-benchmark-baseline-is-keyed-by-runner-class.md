# 0014 — A timing baseline is keyed by runner class, not by hostname

- Status: accepted
- Date: 2026-11-27
- Supersedes: none
- Builds on: [ADR 0011](0011-release-artefact-and-mcp-contract.md) (the release
  gate) and `dev/bench.R`

## Context

`dev/bench.R --gate` compares each case against `dev/bench-baseline.csv` and
fails above 2x. A timing baseline belongs to one machine: the same call differs
by an order of magnitude across CPUs. The script therefore records
`Sys.info()[["nodename"]]` and treats a foreign baseline as `not_run`, not as a
regression.

That rule is right locally and useless in CI. GitHub hands out a fresh
`fv-az…` nodename on every run, so a baseline captured on `ubuntu-latest` can
never match the next run; the gate would print `not_run` forever.

## Decision

1. `dev/bench.R` reads the host key from `TIDYMARKERS_BENCH_HOST` and falls back
   to `Sys.info()[["nodename"]]`.
2. The CI step sets `TIDYMARKERS_BENCH_HOST=github-ubuntu-latest`, a stable name
   for the runner class (image + R release), and the baseline is committed for
   that key.
3. A missing key stays `not_run` with both names in the message and exit 0. A
   cold first run never pretends a regression; the operator records the
   baseline with `--update` on that runner.

## Alternatives

- **Match on nodename.** Rejected: GitHub nodenames are ephemeral, so the gate
  never compares.
- **No host check.** Rejected: a laptop baseline compared on CI reports a
  slowdown that is only a different CPU.
- **Relative within-run timing.** Rejected: it detects instability, not a
  regression against a known cost.

## Consequences

- The key names a runner class, so two physically different `ubuntu-latest`
  machines share one baseline. The 2x tolerance absorbs the hardware spread;
  a tighter budget needs a per-image baseline key.
- A runner image or R upgrade changes the class; the baseline is refreshed
  deliberately, and the gate reports `not_run` until it is.
- `TIDYMARKERS_BENCH_HOST` also lets an operator compare two local machines
  without editing the CSV.

## Verification

```sh
TIDYMARKERS_BENCH_HOST=other Rscript dev/bench.R --gate
# bench-gate: not_run (baseline host vm140.server4, this host other; ...)
TIDYMARKERS_BENCH_HOST=vm140.server4 Rscript dev/bench.R --gate
# bench-gate: ok (every case within 2x of dev/bench-baseline.csv)
```
