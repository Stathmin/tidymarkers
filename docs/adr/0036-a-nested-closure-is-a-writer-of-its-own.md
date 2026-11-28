# 0036 — A nested closure is a writer of its own

- Status: accepted
- Date: 2026-11-28
- Supersedes: none
- Builds on: [ADR 0029](0029-an-output-writer-creates-its-files-not-writes-through-a-link.md)
  (a writer creates its files) and
  [ADR 0035](0035-a-guard-refuses-a-name-a-hit-test-only-reports-one.md)
  (the gate's vocabulary names refusals, and the call form is not the rule)

## Context

`dev/check-writers.R` (ZV45, ZV46) deparsed each top-level
`name <- function(...)` body to text and then grepped it for a write primitive
and a guard. The consequence paragraph of ADR 0035 said a closure defined inside
a function is "outside" the scan. Measured, it is the opposite: the body text
contains the nested definition, so the nested write is *inside* the outer unit
and every guard in the outer text covers it.

Two synthetic bodies, run against the branch at `036cefe`:

```r
outer_guard_launders_nested =
  "function(d) { f <- function(p) writeLines(\"x\", p); refuse_links(d, \"a.csv\"); f(d) }"
string_false_positive =
  "function(x) paste0(\"writeLines(\", x, \")\")"
```

The first was reported `guarded = TRUE`: one unrelated guard in the outer body
laundered a write to `p`, a name the nested closure never declared. The second
was reported as a writer, because the primitive appeared in a string the body
builds. That false positive is why `reprex_from_table()` sat in
`dev/writers.yml` with a "false positive" reason instead of being clean.

Neither is a production escape today: the only nested writer,
`render_design_search()$draw`, renders frames into `withr::local_tempdir()`.
Both are holes in the gate's rule, and the second one hides a real writer behind
a ledger entry that says it is not one.

## Decision

1. **The unit is a function definition, at any depth.** The checker walks the
   parse tree and collects every `function` expression. A top-level definition
   is named `outer`; a definition assigned inside it is `outer/inner`. A closure
   passed as an argument is `<anonymous>` under its enclosing unit. The unit key
   in the gate output and in `dev/writers.yml` is that path, so a nested writer
   is named, not absorbed.
2. **The guard must be in the same unit.** A write call and a refusing call are
   collected from one function body only; the walk does not descend into a
   nested definition, which is scanned as its own unit. A guard in the caller
   does not cover a write in the callee, and a guard in one closure does not
   cover a write in a sibling closure.
3. **The rule matches calls, not text.** A primitive inside a string or a
   comment writes nothing and is not a hit. A call is matched by its head,
   with or without a namespace prefix (ADR 0035, second hole); a computed head
   is not a match.
4. **The ledger follows the name.** `render_design_search` becomes
   `render_design_search/draw` with the same reason; the `reprex_from_table`
   entry is deleted, because the body no longer reports a write. A registered
   symbol that stops writing unguarded still fails as stale, so the deletion is
   proven by the gate, not assumed.

## Evidence

- Self-test grows two bodies. `laundered` must be reported as
  `laundered/f`; `nested_declared` (guard inside the closure) must not be
  reported; `strings` must not be reported. On the old text scan the first was
  `guarded = TRUE` and the third was a hit.
- Full scan of `R/`: `24 writer(s), 17 guard in place, 7 registered` before,
  `23 writer(s), 17 guard in place, 6 registered` after. The one lost hit is
  the `reprex_from_table()` string, which is exactly the false positive; no
  guarded body changed.
- No production write changed. The nested `draw()` writer was already inside a
  private tempdir; it is now registered under its qualified name.

## Consequences

- A new closure that writes a caller's path fails the gate even when its
  enclosing function declares a name. The remedy is a guard inside the closure
  (or, for a private tempdir, a ledger entry on the qualified name).
- The gate output names nested units, so a failure points at the closure that
  writes rather than at the file's top-level symbol.
- The scan is still a source heuristic: it proves a refusing call exists in the
  same function, not that it runs before the write. That limit is unchanged and
  is why the ledger names the declaring caller.
- A function defined by `assign()` or built from text at run time is still
  invisible to the scan. Neither form exists in `R/` today; a future one needs
  the scanner extended, not a ledger entry.

## References

- [`dev/check-writers.R`](../../dev/check-writers.R) — the unit walk, the call
  match and the self-test.
- [`dev/writers.yml`](../../dev/writers.yml) — the ledger by qualified name.
- [`ROADMAP.md`](../../ROADMAP.md) — `ZV47`.
