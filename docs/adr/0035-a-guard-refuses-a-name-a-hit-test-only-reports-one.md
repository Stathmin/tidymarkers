# 0035 — A guard refuses a name; a hit test only reports one

- Status: accepted (amended)
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0029](0029-an-output-writer-creates-its-files-not-writes-through-a-link.md)
  (a writer creates its files, point 1: all names before the first write) and
  [ADR 0033](0033-a-cache-and-a-rendered-file-are-declared-names.md)
  (a cache entry is a declared name, ZV44: a link is not a cache hit)

## Context

`dev/check-writers.R` (ZV45) makes the writer roster a checked ledger: a
top-level function in `R/` that calls a write primitive must call a guard in the
same body, or be justified in `dev/writers.yml`. The gate is only as sound as
its list of guard calls, and that list named `cache_present(`:

```r
guard_pattern <- "refuse_links\\(|declared_path\\(|declared_index\\(|cache_present\\("
```

`cache_present(path)` answers whether the name is a package-made cache hit —
`file_exists()` and not a link. It returns a logical; it refuses nothing and
touches no file. A body such as

```r
function(path) {
  if (cache_present(path)) return(path)
  writeLines("x", path)
}
```

therefore counted as guarded while writing an undeclared name. Measured on the
branch before the fix: with the old pattern, the self-test's cache-tested
writer body matched `cache_present\(` and the gate printed `ok`. With the
pattern narrowed to the three refusals it is reported as
`undeclared writer`. No current writer depended on the false positive — all 17
guarded bodies call `refuse_links()`/`declared_path()`/`declared_index()` — so
the fix changes no production write.

## Decision

1. **A guard is a call that refuses.** `guard_pattern` names only
   `refuse_links()`, `declared_path()` and `declared_index()`. Each throws on a
   link at the declared name. A predicate that returns `TRUE`/`FALSE` is not a
   guard, however useful it is before the write.
2. **The self-test carries the class.** A third synthetic body asks
   `cache_present()` and then writes; it must be flagged like the bare writer.
   A scanner that regains the false positive reddens immediately.
3. **Hit tests stay.** `cache_present()` remains the correct test for "this
   cache entry is ours" (ADR 0033, ZV44) and is still used by readers. The
   decision is about what the gate accepts as a guard, not about the helper.

## Second hole: the call form is not the rule

The same audit found the write pattern accepted one spelling per primitive:
`readr::write_csv` but not a bare `write_csv`. In a package the imports are
attached, so the bare call is the natural form and `writeLines()` was the only
primitive matched without a prefix. A top-level writer written as

```r
function(path) write_csv(x, path)
```

escaped the scanner entirely — not because it was justified, but because the
regex named the namespace. Measured on the branch: the old pattern returns no
match for that body (`grepl()` is `FALSE`); the widened pattern returns
`write_csv`.

Decision: the namespace prefix is optional in `write_pattern`
(`(?:[A-Za-z0-9._]+::)?`), while the primitive names stay the documented list.
The self-test carries the bare imported call next to the namespaced one, so the
two spellings are checked as equal. No current writer relied on the hole: all
24 bodies already used one of the matched spellings, which is why the count did
not move.

## Consequences

- A new writer that tests a cache and writes passes the gate only after it
  declares the name; the ledger cannot be satisfied by a read-side answer.
- The gate's vocabulary matches ADR 0029's rule ("the names are declared
  before the first write"), so a helper added later must be added to
  `guard_pattern` only if it can refuse.
- The check remains a source heuristic: it does not prove the guard runs before
  the write, only that the same body contains a refusing call. That limit is
  unchanged and is why `dev/writers.yml` names the declaring caller for each
  registered helper.
- The scan reads every function definition, not only a top-level one: a closure
  inside a function is its own writer with its own guard (ADR 0036). The early
  text scan flattened nested bodies, so an outer guard covered a nested write;
  that hole is closed and is measured in ADR 0036.

## Amendment 2026-11-28 — the nested-closure limit was stated wrongly

The consequences above originally said a closure defined inside a function was
"outside" the scan. Measured, the deparsed outer body contained the nested
definition, so the nested write was scanned as part of the outer function and
any guard in the outer text covered it. A second artifact came from the same
text match: a body that only built the string `writeLines(` counted as a
writer. ADR 0036 replaces the text scan with a walk over function units and
matches calls. Nothing in this decision changes: a guard is still a refusing
call, and a hit test is still not a guard.
