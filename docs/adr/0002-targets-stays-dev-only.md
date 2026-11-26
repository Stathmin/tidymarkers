# 0002 — `targets` stays a dev-only freshness holder

- Status: accepted
- Date: 2026-11-26
- Supersedes: none

## Context

`targets` is already used here, not proposed: `dev/targets/_targets.R` defines a
short pipeline that hashes what the generated artifacts read and reruns only the
steps a change invalidated. The chain is `engine` → `figures` → `manifest`, plus
the Z8.2 harness (`panel_assays` + `panel_gates` → `z8_outcomes_file`).
`dev/targets.R` is the front end (`--make`, `--outdated`, `--check`), and
`dev/check.R` calls `--make` when `targets` is installed, falling back to
`dev/render-all.R` when it is not.

`dev/` is `.Rbuildignore`d. Every current call site is in `dev/`, and
`DESCRIPTION` has no `targets` entry — while `ROADMAP.md` ZA4 lists `targets`
under Suggests. That is a drift between prose and `DESCRIPTION`, and this ADR
fixes which side is right.

## Decision

`targets` stays a dev-only tool with no `DESCRIPTION` entry.

- No `_targets.R` at the package root and no target that reads or writes under
  `R/` or `inst/`. The engine remains the only author of a figure; `targets`
  decides only *when* a renderer runs.
- No widening to `data-raw/` regeneration or test fixtures without a
  measurement (see the trigger).
- ZA4 is corrected: `targets` is a dev tool, not a package Suggest. Adding it to
  `Suggests` would promise a dependency that nothing under `R/` calls, and
  `R CMD check` does not inspect `dev/` anyway.

## Consequences

- A bare machine still renders: the gate's fallback path stays.
- The gate has two routes to the same figure bytes, which is the point — the
  output is content-hashed, so a diff means the data changed, not the runner.
- `dev/targets/_targets` becomes a second cache store next to `results/`; both
  are `.Rbuildignore`d and untracked.
- A reader who trusts ZA4's old list would try `renv`-style installation of a
  dependency the package does not declare; the roadmap line now says so.

## Alternatives

- Put `targets` in `Suggests`. Rejected: no `R/` call site, so it would only
  satisfy `dev/` scripts that are not shipped.
- Drop it and keep `dev/render-all.R` alone. Rejected: `tar_outdated()` answers
  "would a fresh clone see this stale?" for free, and `dev/check.R` already
  uses it.

## Revisit trigger

Either of these, measured on a warm tree:

- `data-raw/` regeneration or an artifact step costs more than roughly ten
  minutes, so rerunning everything for one edit is the bottleneck; or
- more than three ordered scripts are run by hand every week.

`targets` may then hold those steps too — still under `dev/`, still without a
`DESCRIPTION` entry.

## References

- [`dev/targets/_targets.R`](../../dev/targets/_targets.R) — the pipeline.
- [`dev/targets.R`](../../dev/targets.R) — the front end and the fallback.
- [`dev/check.R`](../../dev/check.R) — the caller.
- [`ROADMAP.md`](../../ROADMAP.md) — ZA4 (dependency policy), ZA8.1/Z8.2.
