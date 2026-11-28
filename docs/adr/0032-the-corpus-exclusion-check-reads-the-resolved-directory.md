# 0032 — The corpus-exclusion check reads the resolved directory

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0030](0030-a-derived-path-component-is-a-declared-name.md)
  (the caller-named root may be a link) and
  [ADR 0031](0031-a-named-root-is-normalised-before-a-sibling-is-derived.md)
  (normalise the named root before deriving from it)

## Context

`marker_export_outcomes()` refuses a handover directory inside the scanned
corpus root. The stated reason is concrete: the next `marker_outcome_corpus()`
scan would read the exported crates back as more corpus, so the handover would
change the data it was taken from.

The check compared spellings: `fs::path_abs(dir)` against `fs::path_abs(root)`.
ADR 0030 allows the caller to name a link; following it is the ordinary meaning
of naming it. But naming a link whose *target* is inside the corpus is a
different path from the one the check measured, and the invariant still fails:

```r
file.symlink(fs::path(root, "handover"), dir)  # dir resolves inside root
marker_export_outcomes(corpus, dir)            # spelling check passes
```

Measured 2026-11-28 on the branch, corpus fixture of `test-outcomes-crate.R`, an
existing `<root>/handover` linked as `dir`: the export succeeded, wrote both
crates under the corpus root, and the next scan would read those two
`outcomes.csv` files as new artefacts. The new test failed (40 pass / 1 fail)
before the check was resolved; after the fix the call is refused and the test
passes (41 pass / 0 fail).

## Decision

1. **The exclusion is decided on resolved paths.** After `fs::dir_create(dir)`,
   the function resolves both `dir` and `root` with `fs::path_real()` and
   repeats the same "inside the scanned corpus root" refusal on the real
   paths.
2. **The link stays allowed, the target is judged.** A link at `dir` that
   resolves outside the root is still a correct handover (ADR 0030); only a
   target inside the root is refused.
3. **The refusal precedes every write.** The identity map and the crates are
   written after both checks, so a refused handover leaves neither.

## Consequences

- Test: `test-outcomes-crate.R`, "a handover linked into the corpus is refused"
  — a link to an existing directory inside the root is refused with the
  existing "outside the scanned corpus root" message; the crate directories are
  not created.
- The check runs after `fs::dir_create(dir)`, so the directory exists when it is
  resolved. A link whose target does not exist is not a corpus leak; it fails
  at resolution with an `fs` error instead of being silently followed.
- This closes the last spelling-versus-location gap in the handover: the
  outside-root rule (this ADR), the sibling identity map (ADR 0031), the derived
  crate components (ADR 0030) and the file names (ADR 0029).
