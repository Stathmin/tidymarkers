# Graphify health (generated)

Run `Rscript dev/graphify-health.R --write` after `graphify update .`.
Do not edit by hand.

The code graph is a **limitation of orientation, not proof of correctness**: an empty graph and clean code look the same.

## Structural invariants

These must stay at zero. `dev/graphify-health.R --check` fails if one moves.

| metric | value |
| --- | ---: |
| dangling endpoints | 0 |
| exact duplicate edges | 0 |
| same-endpoint (collapsed) edges | 0 |
| import cycles | 0 |
| self-loops outside `.tools/` | 0 |

## Known, accepted

None. Every structural invariant above is zero; nothing is accepted as-is.

## Community labels

A label that echoes a filename is honest for a single-file community and misleading for a multi-file one. `--check` warns (never fails) when any label is `missing`, `duplicate` or `echo-stale`.

| status | communities |
| --- | ---: |
| `missing` | 0 |
| `duplicate` | 0 |
| `echo-stale` | 0 |
| `echo-honest` | 0 |
| `named` | 432 |

No missing, duplicate or multi-file file-echo labels.

### Relabel provenance

| date | model | requested | relabelled | kept | note |
| --- | --- | ---: | ---: | ---: | --- |
| 2026-11-16 | `Ternary-Bonsai-2-27B-tight` | 433 | 427 | 6 |  |
| 2026-11-16 | `Ternary-Bonsai-2-27B-tight` | 6 | 1 | 5 |  |
| 2026-11-16 | `manual` | 1 | 1 | 0 | cid 27 kept proposing the rejected Soybean Worked Panel (collides with cid 60, the fixture implementation); renamed by hand to say test-suite, not implementation |
| 2026-11-16 | `Ternary-Bonsai-2-27B-tight` | 46 | 45 | 1 |  |
| 2026-11-16 | `manual` | 1 | 1 | 0 | cid 343 (data-raw/soybean_worked_panel.R, 3 symbols) kept proposing the rejected Soybean Worked Panel, which collides with cid 60; named by hand for the data-raw script, not the fixture |
| 2026-11-16 | `Ternary-Bonsai-2-27B-tight` | 35 | 32 | 3 |  |
| 2026-11-16 | `manual` | 3 | 3 | 0 | strict-purpose: Title Case, 2-5 words, no file/symbol/heading echoes |

## Size

| metric | value |
| --- | ---: |
| nodes | 8 382 |
| edges | 12 196 |
| communities | 432 |
| nodes without `source_file` | 381 |

## Nodes without a `source_file`

Most are expected. `external_function` is the AST record of a call into
another package: it has no file here, and `explain` on one shows its
callers, not its definition. `package_function` is published deliberately by
`dev/graphify-handrolled.R` for a hand-roll's replacement; when the same
`pkg::fn` already exists as an AST node the injection reuses that node, so
`explain "pkg::fn"` resolves to one node and shows the `supersedes` edge.
Two kinds are accidental and unexplained by the extractor: the 13 `concept`
nodes are `import`/`source()` targets -- Python modules (`argparse`, `json`,
`urllib_request`) and R `source()` paths rendered as absolute ids
(`root_development_..._xlsx_read_r`) -- with a raw id for a label and no
file, and `external_class` has the same shape.
`explain` on any of them is vague by construction -- there is no definition
to show. The external-symbol total is a drift signal; the concept total should
be zero, and is not.

| kind | count |
| --- | ---: |
| external_function | 308 |
| package | 55 |
| concept | 14 |
| package_function | 3 |
| external_class | 1 |

## Granularity and interfaces

Graphify's R extractor emits a node per function (a nested closure included), per file, and per assignment, top-level or local. The split is where the noise is: 100 nested function nodes under `R/` (closures, a helper defined inside another function) and 0 nested assignment(s) -- against 2 188 inside `tests/` -- `test_that()` and `local()` bodies mint a node per local -- and 249 in `data-raw/`. Most such nodes are leaves reachable only through their file's `contains` edge; a leaf adds no seam, so the leaf ratio -- not the raw node count -- is the entity-multiplication signal.

| metric | value |
| --- | ---: |
| nodes: `function` | 1 612 |
| nodes: `variable` | 5 565 |
| nodes: `file` | 361 |
| degree-1 (leaf) nodes | 6 330 (76%) |
| …of which `variable` | 5 560 |
| nested (indented) nodes | 2 679 (33% of located) |
| labels shared by >1 node | 664 (3 654 nodes) |
| nodes under `R/` | 1 210 (14%) |
| nodes under `data-raw/`, `tests/`, `dev/` | 6 355 (76%) |
| exported functions with a node | 180 / 180 |
| S3 methods with a node | 91 / 92 |
| generics with a node | 5 / 12 |
| `generic`/`method` pairs joined by an edge | 1 |

Every export (175/175) and all but one S3 method have a node. The absent method is `+.marker_plan()`: it is defined as `"+.marker_plan" <- function`, and graphify's `binding()` accepts only an identifier left-hand side, so a string-named assignment is dropped. The 7 generics without a node are external (`+`, `plot`, `print` and four `vctrs::`): the package implements their methods but does not define the generic, so no package node exists. That is honest, not a gap. graphify's R extractor emits no `generic -> method` relation and no exported/internal flag (`interface_dispatch.py` covers C#, Java, Swift and other languages, not R), so the single pair counted above is a plain `calls` edge (`marker_primers.marker_lax()` calls `marker_primers()` in its body), not dispatch. "What implements `marker_chunk()`" is a label-prefix query, not a traversal. The duplicate labels are not duplicate entities: two `ext()` nodes are two files' helpers, but `explain ext` is ambiguous because R has no namespaces.

### Function / variable nodes: top-level vs nested

| kind | top-level | nested |
| --- | ---: | ---: |
| `function` | 1 419 | 193 |
| `variable` | 3 082 | 2 483 |

### Nested nodes by directory

| dir | nested nodes |
| --- | ---: |
| `tests/` | 2 188 |
| `data-raw/` | 249 |
| `dev/` | 133 |
| `inst/` | 6 |
| `R/` | 100 |

### `R/` files by size

Largest files first; `functions` is the `function` node count for the file.

| file | functions | lines |
| --- | ---: | ---: |
| `R/marker_specificity.R` | 41 | 1 949 |
| `R/marker_plots.R` | 41 | 1 544 |
| `R/marker_criteria.R` | 36 | 1 500 |
| `R/marker_assemble.R` | 22 | 1 455 |
| `R/marker_checks.R` | 15 | 1 062 |
| `R/marker_core.R` | 27 | 987 |
| `R/marker_plot_methods.R` | 47 | 903 |
| `R/marker_lamp.R` | 19 | 884 |
| `R/marker_thermo.R` | 32 | 847 |

## Code ↔ documentation

The graph cannot see `man/`, so this reads the source: for each export it compares roxygen `@param` names with `formals()`, resolving `@rdname` topics (one Rd block documents several aliases). It is a distance, not a faithfulness check: whether the prose describes the body needs a model, and that layer is not built.

| metric | value |
| --- | ---: |
| exports checked | 180 |
| exports without a roxygen block | 0 |
| undocumented parameters | 0 |
| documented-but-absent parameters | 0 |
| missing `@return` | 1 |
| missing `@examples` | 0 |

Undocumented and phantom parameters are gated by `--check`. Missing `@return`: `design_asq_marker`.

## Semantic extraction (local LLM)

`dev/graphify-semantic.sh` is the only step that calls a model. Its output is cached under `graphify-out/cache/semantic/` and keyed by `semantic_hash` in the manifest, so an unchanged file is not re-asked. A graph built by `graphify update .` is AST-only: provenance tags every node and edge it writes `_origin: ast`, and the AST pass itself emits INFERRED `references` edges. So the INFERRED total below is **not** evidence the model ran, and a `concept` label does not prove semantic origin either. Count the pass by `_origin == "semantic"`; `nodes from the model` below is that count and is the only number here that the LLM pass moves. `concept` is a trap: graphify's AST pass mints `concept` nodes for Python `import` and R `source()` targets, so a nonzero `concept` count is an extractor artifact, not model output. For this corpus the LLM pass is not the first lever anyway: the R docs are structured (roxygen, `man/*.Rd`, vignette headings and code chunks), so deterministic stubs and the `@param` check above yield doc-to-code facts with no model.

| metric | value |
| --- | ---: |
| document-like files | 46 |
| with a `semantic_hash` | 0 |
| nodes from the model | 0 |
| concept nodes (any origin) | 14 |
| …of which import/source targets (no `source_file`) | 14 |
| INFERRED edges | 1 803 |
| …from the AST pass | 1 803 |
| …from the model | 0 |

## God nodes

Defined here as the top targets by **in-degree over `calls`**: how often the
rest of the code calls this. (Graphify's own ranking mixes document nodes and
is not reproduced.)

Read this as the **third-party plumbing this package leans on**, not as the package's own interfaces. The R resolver resolves `pkg::fn()` calls to an `external_function` node but a bare internal call only inside one file, so this table is in-degree from one file at a time, and the internal helpers below (which `dev/graphify-calls.R` links across files) are under-counted here. A function called only from tests, `@examples` or a demo is not a seam at all: `ath_genome()` has 0 `R/` callers and is a TAIR10 demo fixture, not an interface.

| node | callers |
| --- | ---: |
| `tibble::tibble()` | 184 |
| `cli::cli_abort()` | 112 |
| `glue::glue()` | 59 |
| `dplyr::mutate()` | 48 |
| `dplyr::filter()` | 46 |
| `purrr::list_rbind()` | 46 |
| `Biostrings::DNAStringSet()` | 43 |
| `fs::path()` | 40 |
| `utils::read.csv()` | 40 |
| `purrr::map()` | 38 |

### Internal seams

The package's own interfaces, measured: distinct files under `R/` that call each function, from the cross-file pass. Broad and shallow is the house style (many small helpers, no hub), so the list is long-tailed rather than headed by one god node.

| function | calling files |
| --- | ---: |
| `snp_template()()` | 15 |
| `map_spec()()` | 13 |
| `map_amplicon()()` | 12 |
| `spec_row()()` | 11 |
| `gc_frac()()` | 10 |
| `scheme_conditions()()` | 10 |
| `marker_design_checks()()` | 9 |
| `primer_tm()()` | 9 |
| `engine_has()()` | 8 |
| `genome_seq()()` | 8 |
| `has_repeat()()` | 8 |
| `map_pair_rows()()` | 8 |

## Known gaps

Adopted from GWASILIY's graph practice.

| gap | status | note |
| --- | --- | --- |
| hand-roll ledger | data | `dev/handrolled.yml` holds the decisions; `dev/check-handrolled.R` reads it |
| `supersedes` graph edges | adopted | `dev/graphify-handrolled.R` injects one `supersedes` edge per entry, reusing the AST `pkg::fn()` node when it exists and creating a `package::function` node otherwise; idempotent, and `dev/check.R` runs it |
| default label prompt on rebuild | guarded | `graphify update .` re-runs the weak built-in prompt on changed communities; the strict pass and the origin ledger record and repair it |
| `dev/relabel-bonsai.R` helper | shipped | repeatable strict pass; defaults to missing/duplicate/echo-stale, `--all-echo` also renames the honest single-file echoes, `--cids=` repairs one name, and a proposal equal to the current label is kept, not counted as a rename. A rebuild that renumbers communities drops label keys; the pass restores them as `missing` and names them |
| vignettes/ and man/ coverage | partial | `.graphifyrc` maps `.Rmd` to markdown, so `vignettes/*.Rmd` **are** indexed (91 `document` nodes). `man/*.Rd` is not, and `.graphifyrc` leaves it unmapped on purpose: mapping it measures as a net loss (189 `document` nodes but communities 432 -> 600, 218 singletons, every curated label renumbered). Read `man/*.Rd` as text for doc<->code faithfulness instead |
| single-file echo labels (254) | skipped | a file-echo name is honest for a one-file community; the strict pass renames only missing/duplicate/echo-stale, `--all-echo` opts in |
| nested-local assignment leaves | accepted | `variable` nodes include function locals: `test_that()`/`local()` bodies in `tests/` mint 2 130 nested nodes ( 2 679 nested nodes in all), and most are leaves with only their file's `contains` edge. They hold no seam; graphify has no per-kind knob, so the only lever is a `.graphifyignore` path exclusion, which also drops the `tests/` map |
| R S3 dispatch edges | gap | every export and 91/92 S3 methods have nodes (`+.marker_plan()` is dropped: a string-named LHS fails graphify's identifier-only `binding()`); the 7 generics without nodes are external implementations. No `generic -> method` relation and no exported/internal flag (`interface_dispatch.py` covers C#, Java, Swift and others, not R), so "what implements `marker_chunk()`" is a label-prefix query |
| code ↔ documentation faithfulness | partial | `dev/graphify-health.R` compares each export's `@param` names with `formals()` and gates undocumented or phantom parameters; whether the prose describes the body is not checked -- that needs the local model, and it is not built |
| cross-file call resolution | shipped | The R resolver binds a bare call only inside the caller's own file, and resolves a cross-file call only through an `imports_from` edge, which it creates from an explicit `source()` call. An R package has neither: it declares dependencies in `NAMESPACE`/`Collate`, not with `source()`, so the AST alone left **0 cross-file** `calls` edges and every community a file island. `dev/graphify-calls.R` closes it with  652 INFERRED edges keyed on a package-wide unique function name (ambiguous names, shadowed locals and test fixtures are skipped, never guessed). Run it after every rebuild: a rebuild is AST-only and drops the pass. |
| semantic layer in graph.json | drift | `graphify update .` writes an AST-only graph (`_origin: ast` on every node and edge); the semantic cache under `graphify-out/cache/semantic/` holds chunk results that are **not** merged into `graph.json`, so a current graph has 0 `_origin: semantic` nodes and 0 `semantic_hash` in the manifest until `dev/graphify-semantic.sh` runs to completion. Re-running it is a decision, not maintenance: the pass does not finish a document at bonsai's 31 tok/s (see the lever table), so budget hours and expect `finish=length` |
| curated label → graph.json | adopted | `dev/relabel-bonsai.R` writes `graphify-out/.graphify_labels.json` and then propagates every name into the `community_name` cached on each node in `graph.json`, so `explain`/`query`/`path` agree with the health table. Propagation runs even when nothing was relabelled, so a rebuild without a rename still picks up names curated in an earlier pass |

## Configuration levers

What is wired into graphify here, and what is deliberately not. `open` is usable in the installed build but not adopted; `skipped` has no consumer; `rejected` was tried and measured to lose more than it gains.

| lever | status | note |
| --- | --- | --- |
| `graphify update .` | shipped | code-only re-extract, no LLM; first step of Rerun below |
| strict label pass + origin ledger | shipped | `dev/relabel-bonsai.R`; all communities named, every pass recorded |
| `supersedes` edges from the hand-roll ledger | shipped | see Known gaps; answers "what replaced this hand-roll" |
| local LLM semantic backend | shipped | `extract --backend=openai` via `dev/graphify-semantic.sh`; reads `OPENAI_BASE_URL`/`OPENAI_MODEL` (bonsai), 2 slots, `--token-budget 12000`. Measured on this box: it never finishes a document. One 18.6 KB doc took 692 s, 6449 prompt tokens and 16000 output tokens, `finish=length`, unfinished; the head is valid JSON but the node ids degenerate (`..._final_final`), and output grows with document length at 31 tok/s, so `ROADMAP.md` would need ~46 min of generation per cycle. The thinking knob is exonerated: a 76-token probe answers clean JSON in 5.5 s with and without it. No per-file semantic skip and no output cap exist (`--token-budget` is prompt-side, `GRAPHIFY_MAX_OUTPUT_TOKENS` a global ceiling). So: routine rebuilds use `graphify update .` (code-only); the semantic pass is a deliberate, rare event, not a step. Coverage and yield in 'Semantic extraction' |
| `dev/graphify-faithfulness.R` (per-export doc<->code) | shipped | one call per export, prose against body, verdict `match`/`partial`/ `contradicts`; 173 exports, median 5.6 s, 103/69/1 on the run of 2026-11-16. Small contexts are the reason this finishes where the document pass cannot: a roxygen block plus 70 body lines stays under ~1.5 k prompt tokens at 31 tok/s. `@rdname` aliases are skipped (their prose belongs to the topic owner) and the prompt forbids `contradicts` when the body delegates or is marked elided: without those two rules 6 of 7 `contradicts` were artefacts. Writes gitignored `graphify-out/faithfulness.md`; not in `dev/check.R` (needs the model) |
| `diagnose multigraph` | shipped | `graphify diagnose multigraph --json` independently confirms 0 dangling, 0 duplicate and 0 collapsed edges, and lists the missing-`source_file` warnings (the external nodes in the table above); it does not replace `--check`, which gates the same invariants plus cycles and labels |
| `export callflow-html` | open | Mermaid architecture/call-flow HTML from the same graph; a candidate for the architecture doc, but no reader is wired to it yet |
| `GRAPHIFY_QUERY_LOG` | open | opt-in log of agent `query`/`explain` calls; would show which nodes agents actually ask for, so the graph can be tuned to real questions. Off by default; enable when there is a reason to read the log |
| `export wiki` | shipped | 411 articles in gitignored `graphify-out/wiki/`, one per community plus an index; template-only, no LLM, so its quality is the curated labels |
| git hook (`graphify hook install`) | shipped | post-commit/checkout/merge/rewrite plus the graph merge driver; it runs `update` only (no LLM), so `dev/graphify-health.R --write` stays manual |
| `.graphifyignore` | shipped | `data-raw/proposed/`, `results/`, `.tmp/`, `.tools/` are design scratch or vendored tooling, not package API; a rebuild prunes them, and the file can only ever add exclusions |
| per-kind node pruning | skipped | no CLI or rc knob drops `variable` leaves; only a `.graphifyignore` path exclusion can, at the cost of the `tests/`/`data-raw/`/`dev/` map |
| `--mode deep` | skipped | richer INFERRED edges, but a full re-extract and it re-names changed communities with the weak built-in prompt; no consumer here |
| `dev/graphify-calls.R` | shipped | the AST pass resolves `calls` only inside one file, so every community was a file island; this deterministic post-pass re-reads `R/` with srcref attribution and adds  652  cross-file edges ( 652  carry `_origin = calls-pass`) |
| `language.rd=markdown` (man/ in the graph) | rejected | adds 189 `document` nodes, but takes communities 432 -> 600 with 218 singletons (one island per man page) and renumbers every curated label. A gitignored `.md` mirror cannot work either: graphify honours .gitignore, so the stubs never enter the graph. Doc<->code faithfulness reads `man/*.Rd` as text instead |
| `graphify-mcp` | shipped | 0.9.92 does ship an MCP server (`graphify-mcp`, module `graphify.serve`): stdio or Streamable HTTP over the same `graph.json`, with `--graph`, `--port` and `--stateless`. Verified by `graphify-mcp --help`; register it in the editor's `context_servers` to query the graph without a shell-out |
| `query` / `path` / `explain` / `affected` | shipped | verified from the shell: `explain "Biostrings::complement"` resolves one node and shows its `supersedes` edge, `affected "normalise_allele()"` lists its caller; the same four are reachable over MCP (see the row above) |
| `export neo4j` / `graphml` / `obsidian` / `svg` | skipped | no consumer in this repository |

## Provenance

Built at commit `24bb510bf84c93820a872a158de757f3d2079043`.
`HEAD` is `35892c9ea0630e5ab433c86f74146aeba7dd4dc5`, so the graph is **stale**; rerun `graphify update .` and this script.

## Rerun

```sh
graphify update .
Rscript dev/graphify-handrolled.R
Rscript dev/graphify-health.R --write
graphify export wiki
graphify export html
# when a label warns: Rscript dev/relabel-bonsai.R [--all-echo]
# one wrong name: Rscript dev/relabel-bonsai.R --cids=N
# semantic pass (opt-in, local LLM): dev/graphify-semantic.sh
# refresh the analysis sidecar: graphify cluster-only . --no-label
```

Do not pass `--no-cluster` to `update`: it drops the community
assignment and duplicates some edges, so the health numbers above
mislead. `cluster-only` is deterministic on an unchanged graph but can
renumber after a rebuild, so re-check the labels and rerun
`dev/relabel-bonsai.R` if it warns.
