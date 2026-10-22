---
name: r-profiling
description: Profile an R pipeline with Rprof, read raw stacks, export a tidy profile (+ RDF), and check the claim it supports.
---

# R profiling — measure, then claim

Use when a call is slow, when a "fix" is about to be written without a profile,
or when a timing claim is about to enter a report. Local gate
`dev/check.R` and `AGENTS.md` win over this skill.

## The one rule

A profile measures where the **interpreter** spent samples, not where the
**algorithm** spends work. A fast function called 184 149 times and a slow one
called six times look alike in the top rows. Never write "X dominates" from
`self` alone: pair every number with a count from the code ("how many hits did
it report?"). A profiling claim is a **hypothesis** until the smallest case is
reproduced.

## Workflow

1. **Bounded run.** Never profile an unbounded call. `timeout 300 Rscript …`,
   and profile to a file:
   `Rprof(".tmp/prof/x.out", interval = 0.01)`.
2. **Read raw stacks even when the run is killed.** A partial profile is still
   a profile — the samples before the kill name the hot stack. This is how the
   KASP hang was found: variant 2 never finished, the raw stacks did.
3. **One hypothesis**, written down: "`mismatch_profile` runs once per hit and
   is called 184 149 times by `orient`".
4. **Smallest case.** Reproduce in a standalone `Rscript` file **without the
   package**: take the one call, one input. If the standalone version is fast,
   the cost is in the caller, not the callee.
5. **Fix, then prove `all.equal()`** against the saved pre-change output on
   real inputs (10 cases) and edge inputs (indels, `N`, multi-base codes).
6. **Re-profile.** Report before/after on the same input, same interval.

## Tool: `dev/profile.R`

```sh
export PATH=/opt/R-4.6.1/bin:$PATH
export UNAFOLDDAT=/root/development/tidymarkers/.tools/share/oligoarrayaux
export R_LIBS_SITE=/opt/R-4.6.1/lib/R/library

Rscript dev/profile.R run    .tmp/probe/v2.R          # source it under Rprof
Rscript dev/profile.R stacks .tmp/probe/p2.rprof      # read an existing file
```

Writes to `.tmp/prof/<name>/`:

| file | contents |
|---|---|
| `<name>.self.csv` | per function: `self_seconds`, `total_seconds`, `pct_self`, `pct_under`, sample counts |
| `<name>.frames.csv` | long table: one row per frame per sample — `sample`, `depth`, `function_name`, `caller` |
| `<name>.ttl` | the same facts as RDF 1.1 N-Triples |
| `<name>.rds` | the tables, for reuse in a script |

From R (after `source("dev/profile.R")`):

```r
x <- profile_read(".tmp/probe/p2.rprof")
profile_report(x)                     # top rows, one screen
profile_callers(x, "make_as_primer")  # who repeated it — summaryRprof cannot
```

`summaryRprof()` aggregates by function name and **hides callers**. The
`frames` table keeps them, so "who repeats this?" is a `table()`.

## Reading the rows

- **S4 churn.** A profile topped by `validObject`, `initialize`, `getClassDef`,
  `standardGeneric`, `getClassesFromCache` is **object construction** (IRanges /
  DataFrame / DNAStringSet here), not the algorithm. Build fewer objects; do not
  rewrite the maths.
- **Grid multiplication.** Count the objects a call builds (strands x lengths x
  positions x bases) before optimising any one. Multiply, do not iterate:
  vectorise the inner computation over the whole hit vector.
- **Cheap doomed path first.** When a whole grid is built and then discarded,
  answer the discard condition before the grid exists (`as_grid_can_pass()` is
  the precedent).
- **`.Call`/`.Call2`/`.XString.matchPattern` at the top** means the cost is in
  C library loops reached many times — reduce the number of calls, not the
  per-call code.
- **`compiler:::tryCmpfun` / `genCode` / `mayCallBrowser`** is byte-compiling a
  closure on first use. Pay it once, outside the loop.

## RDF export — why, and its limit

No RDF serialiser is installed here (no `rdflib`, no `RDF`/`redland`/`jsonld`),
and this repo adds no dependency for a dev tool. **N-Triples** is line-oriented,
so base R writes it: `<s> <p> <o> .` per line. Vocabulary is
`https://tidymarkers.org/profile#`:

- `ProfileRun` — `script`, `commit`, `elapsedSeconds`, `sampleInterval`, `samples`
- `Function` — `functionName`, `selfSeconds`, `totalSeconds`, `selfSamples`,
  `totalSamples`, `run`
- edges: `<run> hasFunction <fn>`, `<fn> calledBy <caller>`

Numbers are written as `xsd:decimal` (`-?digits.digits`: no bare integer, no
exponent, no leading zeros).

**Reading it back in R — no external parser needed.** The profiler ships both
readers, so an audit re-runs in R alone:

- `profile_read_tidy(path)` — `readr::read_rds()` on the `.rds`. This is the
  object to prefer: typed columns, no parsing, no ambiguity.
- `profile_ttl_read(path)` — a real N-Triples reader, not a grep. Columns
  `subject`, `predicate`, `object`, `object_iri`, `datatype`, `value`; a
  decimal comes back as a number, so a report recomputes from the file alone.
  A malformed line `stop()`s with its line number instead of being dropped.

No RDF package is installed (`rdf`, `ntriples`, `turtle`, `sparql`: none
exported by any installed package; no `rapper`/`riot` binary; no Python `rdflib`).
`redland`, `rdflib` and `jsonld` **are** on CRAN and deliberately not added —
a dev tool does not owe the package a dependency.

Two parsing facts on this build, both cost real time:

- `regexec` returns the **whole match** as element 1, then groups. An offset
  computed from `m[[1L]]` is the full match, not group 1.
- Escaped patterns written as one literal in a file can fail to match where
  the same pattern built with `paste0` succeeds. Build the pattern.

Why not an R debugging standard: R's call graph knows nothing about RDF, and a
profile is a graph of samples, not of resources. RDF is for **interchange** —
several runs across commits, compared in a store. For one run, the tidy
`.rds`/`.csv` is the better object. The export exists so a claim can be audited
later without re-running.

## Accuracy pass — ask before reporting

1. Is the metric the right one? `self` for "what is hot", `total` for "what
   subtree costs", sample counts for "how often called".
2. Is the interval short enough? At `0.01 s` a 2 s run has 200 samples — the
   top row is +/- a few percent. Say so; do not quote 3 significant figures.
3. Did the fix change behaviour? `all.equal()` the output, not just the time.
4. Is the comparison fair? Same input, same interval, same commit lineage.
5. Would the number be true on a bigger input? Say "on variant 2" — not "in
   general".

## Cache gotcha

`pkgload::load_all()` caches compiled/lazyloaded bodies: an edited function can
still behave like the old body. Set `recompile = TRUE` before believing a
profile of edited code. R 4.6.1 also byte-compiles on first use, so the first
iteration of a loop is not representative.

## Shell notes

`sh`, not bash: no `time`, no `<<<`, and multi-line `Rscript -e` breaks — write
a `.R` file under `.tmp/probe/`.
