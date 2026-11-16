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
| `named` | 428 |

No missing, duplicate or multi-file file-echo labels.

### Relabel provenance

| date | model | requested | relabelled | kept | note |
| --- | --- | ---: | ---: | ---: | --- |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 74 | 69 | 5 |  |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 5 | 5 | 0 |  |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 3 | 3 | 0 |  |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 251 | 249 | 2 |  |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 2 | 2 | 0 |  |
| 2026-11-06 | `manual` | 1 | 1 | 0 | cid 408: orphan external yaml node; the model kept proposing a rejected placeholder, fixed by hand |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 4 | 2 | 2 | cids 210/217 re-asked with --cids after the unchanged-check fix; the model repeated the same honest single-file name, so both stay |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 38 | 37 | 1 | after the .graphifyignore rebuild renumbered communities and dropped 38 label keys; 37 named here, cid 264 fixed by hand |
| 2026-11-06 | `manual` | 1 | 1 | 0 | cid 264: single-file SECURITY.md community; model proposed the member echo 'Security Policy', fixed by hand |
| 2026-11-06 | `Ternary-Bonsai-2-27B-tight` | 4 | 4 | 0 | targeted --cids=190,210,217,246; thin but honest names deepened to say what the cluster is for |
| 2026-11-10 | `Ternary-Bonsai-2-27B-tight` | 12 | 12 | 0 |  |
| 2026-11-10 | `Ternary-Bonsai-2-27B-tight` | 88 | 85 | 3 |  |
| 2026-11-10 | `Ternary-Bonsai-2-27B-tight` | 3 | 1 | 2 |  |
| 2026-11-10 | `manual` | 2 | 2 | 0 |  |
| 2026-11-10 | `Ternary-Bonsai-2-27B-tight` | 41 | 39 | 2 |  |
| 2026-11-10 | `Ternary-Bonsai-2-27B-tight` | 2 | 1 | 1 |  |
| 2026-11-10 | `manual` | 1 | 1 | 0 |  |
| 2026-11-10 | `Ternary-Bonsai-2-27B-tight` | 35 | 34 | 1 |  |
| 2026-11-10 | `manual` | 1 | 1 | 0 | cid 413: the rebuild split propagate() off the relabel script; kept distinct from "Bonsai Relabeling Script" by hand |
| 2026-11-12 | `Ternary-Bonsai-2-27B-tight` | 2 | 2 | 0 |  |
| 2026-11-12 | `Ternary-Bonsai-2-27B-tight` | 89 | 86 | 3 |  |
| 2026-11-16 | `Ternary-Bonsai-2-27B-tight` | 104 | 101 | 3 |  |
| 2026-11-16 | `Ternary-Bonsai-2-27B-tight` | 3 | 3 | 0 |  |

## Size

| metric | value |
| --- | ---: |
| nodes | 8 268 |
| edges | 11 345 |
| communities | 428 |
| nodes without `source_file` | 380 |

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
| concept | 13 |
| package_function | 3 |
| external_class | 1 |

## Granularity and interfaces

Graphify's R extractor emits a node per function (a nested closure included), per file, and per assignment, top-level or local. The split is where the noise is: 100 nested nodes under `R/` -- all closures (a helper defined inside another function), no nested assignments (0) -- against 2 170 inside `tests/` -- `test_that()` and `local()` bodies mint a node per local -- and 249 in `data-raw/`. Most such nodes are leaves reachable only through their file's `contains` edge; a leaf adds no seam, so the leaf ratio -- not the raw node count -- is the entity-multiplication signal.

| metric | value |
| --- | ---: |
| nodes: `function` | 1 587 |
| nodes: `variable` | 5 430 |
| nodes: `file` | 355 |
| degree-1 (leaf) nodes | 6 269 (76%) |
| …of which `variable` | 5 425 |
| nested (indented) nodes | 2 615 (33% of located) |
| labels shared by >1 node | 641 (3 554 nodes) |
| nodes under `R/` | 1 206 (15%) |
| nodes under `data-raw/`, `tests/`, `dev/` | 6 193 (75%) |
| exported functions with a node | 180 / 180 |
| S3 methods with a node | 91 / 92 |
| generics with a node | 5 / 12 |
| `generic`/`method` pairs joined by an edge | 1 |

Every export (175/175) and all but one S3 method have a node. The absent method is `+.marker_plan()`: it is defined as `"+.marker_plan" <- function`, and graphify's `binding()` accepts only an identifier left-hand side, so a string-named assignment is dropped. The 7 generics without a node are external (`+`, `plot`, `print` and four `vctrs::`): the package implements their methods but does not define the generic, so no package node exists. That is honest, not a gap. graphify's R extractor emits no `generic -> method` relation and no exported/internal flag (`interface_dispatch.py` covers C#, Java, Swift and other languages, not R), so the single pair counted above is a plain `calls` edge (`marker_primers.marker_lax()` calls `marker_primers()` in its body), not dispatch. "What implements `marker_chunk()`" is a label-prefix query, not a traversal. The duplicate labels are not duplicate entities: two `ext()` nodes are two files' helpers, but `explain ext` is ambiguous because R has no namespaces.

### Function / variable nodes: top-level vs nested

| kind | top-level | nested |
| --- | ---: | ---: |
| `function` | 1 400 | 187 |
| `variable` | 3 005 | 2 425 |

### Nested nodes by directory

| dir | nested nodes |
| --- | ---: |
| `tests/` | 2 170 |
| `data-raw/` | 249 |
| `dev/` | 87 |
| `inst/` | 6 |
| `R/` | 100 |

### `R/` files by size

Largest files first; `functions` is the `function` node count for the file.

| file | functions | lines |
| --- | ---: | ---: |
| `R/marker_specificity.R` | 41 | 1 949 |
| `R/marker_plots.R` | 41 | 1 544 |
| `R/marker_criteria.R` | 36 | 1 495 |
| `R/marker_assemble.R` | 22 | 1 436 |
| `R/marker_checks.R` | 15 | 1 055 |
| `R/marker_core.R` | 27 | 975 |
| `R/marker_plot_methods.R` | 47 | 902 |
| `R/marker_lamp.R` | 19 | 880 |
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
| document-like files | 48 |
| with a `semantic_hash` | 0 |
| nodes from the model | 0 |
| concept nodes (any origin) | 13 |
| …of which import/source targets (no `source_file`) | 13 |
| INFERRED edges | 1 077 |
| …from the AST pass | 1 077 |
| …from the model | 0 |

## God nodes

Defined here as the top targets by **in-degree over `calls`**: how often the
rest of the code calls this. (Graphify's own ranking mixes document nodes and
is not reproduced.)

Read this as the **third-party plumbing this package leans on**, not as the package's own interfaces. The R resolver resolves `pkg::fn()` calls to an `external_function` node but a bare internal call only inside one file (see 'cross-file call resolution' below), so no internal helper can accumulate the in-degree a `tibble::tibble()` does. The real internal seams are the R-only cross-file fan-out measured from source (comments stripped, only `R/` callers, test fixtures and roxygen examples excluded): `snp_template` (16 files), `spec_row` (12), `map_spec` (13), `scheme_conditions` (10), `gc_frac` (10), `primer_candidates` (9), `primer_tm` (9), `has_palindrome` (7), `ath_chrom` (5) -- broad and shallow, with no hub, not this table. A function called only from tests, `@examples` or a demo is not a seam at all: `ath_genome()` has 0 `R/` callers and is a TAIR10 demo fixture, not an interface. This list is hand-maintained and drifts; compute it in the script from the same stripped source scan before trusting the counts.

| node | callers |
| --- | ---: |
| `tibble::tibble()` | 184 |
| `cli::cli_abort()` | 111 |
| `glue::glue()` | 59 |
| `dplyr::mutate()` | 48 |
| `dplyr::filter()` | 46 |
| `purrr::list_rbind()` | 46 |
| `Biostrings::DNAStringSet()` | 43 |
| `utils::read.csv()` | 40 |
| `fs::path()` | 39 |
| `purrr::map()` | 38 |

## Known gaps

Adopted from GWASILIY's graph practice.

| gap | status | note |
| --- | --- | --- |
| hand-roll ledger | data | `dev/handrolled.yml` holds the decisions; `dev/check-handrolled.R` reads it |
| `supersedes` graph edges | adopted | `dev/graphify-handrolled.R` injects one `supersedes` edge per entry, reusing the AST `pkg::fn()` node when it exists and creating a `package::function` node otherwise; idempotent, and `dev/check.R` runs it |
| default label prompt on rebuild | guarded | `graphify update .` re-runs the weak built-in prompt on changed communities; the strict pass and the origin ledger record and repair it |
| `dev/relabel-bonsai.R` helper | shipped | repeatable strict pass; defaults to missing/duplicate/echo-stale, `--all-echo` also renames the honest single-file echoes, `--cids=` repairs one name, and a proposal equal to the current label is kept, not counted as a rename. A rebuild that renumbers communities drops label keys; the pass restores them as `missing` and names them |
| vignettes/ and man/ coverage | partial | `.graphifyrc` maps `.Rmd` to markdown, so `vignettes/*.Rmd` **are** indexed (91 `document` nodes). `man/*.Rd` is not, and `.graphifyrc` leaves it unmapped on purpose (174 files would dominate a semantic run); export plain `.md` stubs via `dev/graphify-docs.R` if the manuals must be graph concepts |
| single-file echo labels (254) | skipped | a file-echo name is honest for a one-file community; the strict pass renames only missing/duplicate/echo-stale, `--all-echo` opts in |
| nested-local assignment leaves | accepted | `variable` nodes include function locals: `test_that()`/`local()` bodies in `tests/` mint 2 113 nested nodes ( 2 615 nested nodes in all), and most are leaves with only their file's `contains` edge. They hold no seam; graphify has no per-kind knob, so the only lever is a `.graphifyignore` path exclusion, which also drops the `tests/` map |
| R S3 dispatch edges | gap | every export and 91/92 S3 methods have nodes (`+.marker_plan()` is dropped: a string-named LHS fails graphify's identifier-only `binding()`); the 7 generics without nodes are external implementations. No `generic -> method` relation and no exported/internal flag (`interface_dispatch.py` covers C#, Java, Swift and others, not R), so "what implements `marker_chunk()`" is a label-prefix query |
| code ↔ documentation faithfulness | partial | `dev/graphify-health.R` compares each export's `@param` names with `formals()` and gates undocumented or phantom parameters; whether the prose describes the body is not checked -- that needs the local model, and it is not built |
| cross-file call resolution | gap | The R resolver binds a bare call only inside the caller's own file, and resolves a cross-file call only through an `imports_from` edge, which it creates from an explicit `source()` call. An R package has neither: it declares dependencies in `NAMESPACE`/`Collate`, not with `source()`, so of 747 `calls` edges whose ends both have a `source_file`, **0 cross a file**. `pkg::fn()` targets do resolve cross-file because they name the package. Consequence: `explain`/`affected` under-report impact, communities are file islands, and a source-scan (not the graph) is the only way to see the real seams. Fix: a deterministic post-pass keyed on a package-wide unique function name (skip ambiguous names, shadowed locals and test fixtures), emitting INFERRED `calls` edges -- the same pattern `dev/graphify-handrolled.R` uses for `supersedes`. |
| semantic layer in graph.json | drift | `graphify update .` writes an AST-only graph (`_origin: ast` on every node and edge); the semantic cache under `graphify-out/cache/semantic/` holds chunk results that are **not** merged into `graph.json`, so a current graph has 0 `_origin: semantic` nodes and 0 `semantic_hash` in the manifest until `dev/graphify-semantic.sh` runs to completion. Re-run it when the bonsai slots are free (see the lever note on slot contention) |
| curated label → graph.json | adopted | `dev/relabel-bonsai.R` writes `graphify-out/.graphify_labels.json` and then propagates every name into the `community_name` cached on each node in `graph.json`, so `explain`/`query`/`path` agree with the health table. Propagation runs even when nothing was relabelled, so a rebuild without a rename still picks up names curated in an earlier pass |

## Configuration levers

What is wired into graphify here, and what is deliberately not. `open` is usable in the installed build but not adopted; `skipped` has no consumer.

| lever | status | note |
| --- | --- | --- |
| `graphify update .` | shipped | code-only re-extract, no LLM; first step of Rerun below |
| strict label pass + origin ledger | shipped | `dev/relabel-bonsai.R`; all communities named, every pass recorded |
| `supersedes` edges from the hand-roll ledger | shipped | see Known gaps; answers "what replaced this hand-roll" |
| local LLM semantic backend | shipped | `extract --backend=openai` via `dev/graphify-semantic.sh`; reads `OPENAI_BASE_URL`/`OPENAI_MODEL` (bonsai), 2 slots, `--token-budget 12000`. `GRAPHIFY_MAX_OUTPUT_TOKENS` must stay at 16000: a lower cap truncates, splits and retries forever. The server has `-np 2`, so a second client on the same port (a Merovingen `acceptance_run.py` did this) holds both slots with long generations and this pass blocks as the deferred request -- that, not a config fault, is why a few chunks can take minutes. Check `llamacpp:requests_deferred` and that no other job owns the slots before blaming the backend. Coverage and yield in 'Semantic extraction' |
| `diagnose multigraph` | shipped | `graphify diagnose multigraph --json` independently confirms 0 dangling, 0 duplicate and 0 collapsed edges, and lists the missing-`source_file` warnings (the external nodes in the table above); it does not replace `--check`, which gates the same invariants plus cycles and labels |
| `export callflow-html` | open | Mermaid architecture/call-flow HTML from the same graph; a candidate for the architecture doc, but no reader is wired to it yet |
| `GRAPHIFY_QUERY_LOG` | open | opt-in log of agent `query`/`explain` calls; would show which nodes agents actually ask for, so the graph can be tuned to real questions. Off by default; enable when there is a reason to read the log |
| `export wiki` | shipped | 411 articles in gitignored `graphify-out/wiki/`, one per community plus an index; template-only, no LLM, so its quality is the curated labels |
| git hook (`graphify hook install`) | shipped | post-commit/checkout/merge/rewrite plus the graph merge driver; it runs `update` only (no LLM), so `dev/graphify-health.R --write` stays manual |
| `.graphifyignore` | shipped | `data-raw/proposed/`, `results/`, `.tmp/`, `.tools/` are design scratch or vendored tooling, not package API; a rebuild prunes them, and the file can only ever add exclusions |
| per-kind node pruning | skipped | no CLI or rc knob drops `variable` leaves; only a `.graphifyignore` path exclusion can, at the cost of the `tests/`/`data-raw/`/`dev/` map |
| `--mode deep` | skipped | richer INFERRED edges, but a full re-extract and it re-names changed communities with the weak built-in prompt; no consumer here |
| `query` / `path` / `explain` / `affected` | shipped | verified from the shell: `explain "Biostrings::complement"` resolves one node and shows its `supersedes` edge, `affected "normalise_allele()"` lists its caller. 0.9.92 has no MCP server, so agents use this CLI |
| `export neo4j` / `graphml` / `obsidian` / `svg` | skipped | no consumer in this repository |

## Provenance

Built at commit `3300637e9ffeb99f11384366584c445e3688f81e`.
The graph matches `HEAD`.

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
