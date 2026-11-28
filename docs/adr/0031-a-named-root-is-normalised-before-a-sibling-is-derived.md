# 0031 — A named root is normalised before a sibling is derived from it

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0029](0029-an-output-writer-creates-its-files-not-writes-through-a-link.md)
  (point 4: a file beside the directory is still a declared name) and
  [ADR 0030](0030-a-derived-path-component-is-a-declared-name.md) (a derived
  path component is a declared name)

## Context

`marker_export_outcomes()` writes the identity map as a sibling of the
handover directory: `<dir>-identity.csv`. The name is built by string
concatenation, and the directory was used exactly as the caller spelled it.
A trailing separator changes what that concatenation means:

```r
dir <- "/tmp/handover/out/"
paste0(dir, "-identity.csv")  # "/tmp/handover/out/-identity.csv"
```

The path is now *inside* the directory the function just created, not beside
it. The guard and the write still agreed with each other, so nothing failed:
the map was written as a hidden `-identity.csv` in the handover, and `dir`
looked normal. The other checks resolve `dir` with `fs::path_abs()`, which drops
the trailing separator, so the outside-root check and the sibling derivation
disagreed about which path `dir` names.

This is the promise of ADR 0029 point 4 in the other direction: the name
`<dir>-identity.csv` is declared, but the *base* of the derivation was never
normalised. A handover whose identity map travels with the pseudonymised crates
re-links every pseudonym, which is exactly what the handover exists to prevent.

Measured 2026-11-28 on the branch, corpus fixture of `test-outcomes-crate.R`,
`dir = paste0(tempdir, "/out/")`: the map was written at `<dir>/-identity.csv`
(inside), and `<dir>-identity.csv` did not exist. Before the fix the new test
recorded two failures; after it, 40 pass / 0 fail for the file.

## Decision

1. **Normalise the named root once, before deriving anything from it.**
   `marker_export_outcomes()` resolves `dir` with `fs::path_abs()` immediately
   after validating it, and every later use (the outside-root check, the sibling
   map, the crates) reads the normalised value.
2. **A sibling name is derived from a path, not from the caller's spelling of
   it.** String concatenation is allowed only for a normalised base; without a
   trailing separator, `<dir>-identity.csv` is beside `dir` by construction.
3. **The fix is normalisation, not a second guard.** The identity map remains
   guarded by ADR 0029 point 4; this ADR removes the spelling that moved it.

## Consequences

- A trailing separator no longer changes where the map is written. Test:
  `test-outcomes-crate.R`, "a trailing slash does not put the identity map
  inside the handover" — the map exists beside the directory, the hidden
  `-identity.csv` inside it does not, and both crates are still written.
- The rule is general for this family: a path the caller names may be a link
  (ADR 0030, consequence 3), but its *spelling* is not a second name. Normalise
  it before deriving sibling names or comparing it with a scanned root.
- `fs::path_abs()` makes a relative `dir` absolute; the handover already stores
  only absolute paths in the crate metadata, so no output changes for callers
  that passed a relative path.
