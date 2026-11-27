# 0002 — The CI gate fails when a linter, a test runner or the tests are absent

- Date: 2026-11-27
- Status: accepted
- Supersedes: the informal "a missing tool is a message" policy

## Context

`dev/check.R` already fails when `styler` or `lintr` is absent, but CI never ran
that script: the workflow called `R CMD check` and the licence gate only. A
missing linter therefore could not red a pull request, and an empty
`tests/testthat/` would still produce a green `R CMD check`. A gate that passes
because the tool did not run proves nothing.

## Decision

Strict policy, enabled now:

1. **A missing tool is a failure.** `styler`, `lintr` and `testthat` are
   required by `dev/check.R`; a missing one stops the gate with the tool name.
2. **A missing test surface is a failure.** `dev/check.R` requires at least one
   `tests/testthat/test-*.R` file and a loadable `testthat`.
3. **CI runs the same gate.** The `ubuntu-latest (release)` leg installs
   `styler`, `lintr` and `roxygen2` as extra packages and runs
   `Rscript dev/check.R --quick --no-render` before `R CMD check`. The other
   legs keep `R CMD check`; style and lint are platform-independent, so one leg
   is enough to fail the merge without quadrupling the wall clock.
4. **`R CMD check --as-cran` remains the test runner.** It still executes the
   suite on every platform, with `error-on: "warning"`.

## Consequences

- A pull request that removes a linter, removes the tests, or breaks style is
  red before the expensive matrix finishes.
- A local `Rscript dev/check.R --quick --no-render` is now byte-for-byte the
  policy CI enforces, so "green locally" and "green in CI" mean the same gate.
- The gate needs `roxygen2` in CI to compare `man/`; installing it is the only
  new dependency, and it was already a `Suggests` for maintainers.
