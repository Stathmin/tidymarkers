# Graphify health (generated)

Run `Rscript dev/graphify-health.R --write` after `graphify update .`.
Do not edit by hand.

The code graph is a **limitation of orientation, not proof of correctness**: an empty graph and clean code look the same.

## Structural invariants

These must stay at zero. `dev/graphify-health.R --check` fails if one moves.

| metric | value |
| --- | ---: |
| dangling endpoints | 0 |
| self-loops | 1 |
| exact duplicate edges | 0 |
| same-endpoint (collapsed) edges | 0 |
| import cycles | 0 |

## Size

| metric | value |
| --- | ---: |
| nodes | 7 597 |
| edges | 10 338 |
| communities | 398 |
| nodes without `source_file` | 349 |

## Nodes without a `source_file`

Symbols from other packages are recorded as nodes but have no file in this
repository. They are expected; the count is a drift signal, not a defect.

| kind | count |
| --- | ---: |
| external_function | 290 |
| package | 54 |
| concept | 4 |
| external_class | 1 |

## God nodes

Defined here as the top targets by **in-degree over `calls`**: how often the
rest of the code calls this. (Graphify's own ranking mixes document nodes and
is not reproduced.)

| node | callers |
| --- | ---: |
| `tibble::tibble()` | 171 |
| `cli::cli_abort()` | 94 |
| `glue::glue()` | 59 |
| `dplyr::filter()` | 46 |
| `dplyr::mutate()` | 45 |
| `purrr::list_rbind()` | 45 |
| `Biostrings::DNAStringSet()` | 41 |
| `utils::read.csv()` | 37 |
| `purrr::map()` | 36 |
| `Biostrings::DNAString()` | 34 |

## Provenance

Built at commit `5d2e48c5a4cb85a983abee96cbfd4cb141544188`.
The graph matches `HEAD`.

## Rerun

```sh
graphify update .
Rscript dev/graphify-health.R --write
```

Do not pass `--no-cluster`: it drops the community assignment and
duplicates some edges, so the health numbers above mislead.
