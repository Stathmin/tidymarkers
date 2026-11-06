# =============================================================================
#  dev/graphify-health.R  --  is the code graph healthy, and what does it miss?
# =============================================================================
#  `graphify-out/` is a developer tool, not a CI dependency: the gate must run
#  on a bare machine, so a missing graph is a message and this script exits 0.
#
#  The graph is a way to find a file, not a proof that the code is right.  The
#  report names that limit in a fixed sentence (checked by dev/check-docs.R),
#  because "the graph is empty" and "the code is clean" look the same.
#
#  Metrics are read straight from graph.json, so no `graphify` binary is needed:
#    * structural invariants that must stay at zero (dangling endpoints,
#      duplicate edges, import cycles) -- a change here is a real regression;
#    * descriptive counts (nodes, edges, communities, nodes without a
#      source_file, god nodes) -- drift is a warning, not a failure: editing
#      any file moves them, and a gate that fires on every edit is not a gate.
#
#  Run:  Rscript dev/graphify-health.R --write   (refresh the report)
#        Rscript dev/graphify-health.R --check   (fail on a regression)
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

GRAPH <- "graphify-out/graph.json"
LABELS <- "graphify-out/.graphify_labels.json"
ORIGIN <- "graphify-out/.graphify_labels.origin.json"
REPORT <- "docs/generated/graphify-health.md"

# stable sentence: dev/check-docs.R asserts it verbatim, so do not reword here
# without rewording there (that coupling is the point -- it is the roadmap's check)
LIMITATION <- paste(
  "The code graph is a **limitation of orientation, not proof of correctness**:",
  "an empty graph and clean code look the same."
)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--write" %in% args) "write" else if ("--check" %in% args) "check" else NA
if (is.na(mode)) {
  message("usage: Rscript dev/graphify-health.R --write | --check")
  quit(status = 2L)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  message("skip (jsonlite not installed)")
  quit(status = 0L)
}
if (!file.exists(GRAPH)) {
  message("skip (no ", GRAPH, " - run `graphify update .` first; CI has no graph)")
  quit(status = 0L)
}

g <- jsonlite::fromJSON(GRAPH, simplifyVector = FALSE)
nodes <- g$nodes
links <- g$links

node_id <- function(n) n$id %||% ""
link_rel <- function(e) e$relation %||% e$type %||% ""
link_end <- function(e, which) e[[which]] %||% ""

ids <- vapply(nodes, node_id, character(1))
labels <- vapply(nodes, function(n) n$label %||% n$id %||% "", character(1))
has_source <- vapply(nodes, function(n) nzchar(n$source_file %||% ""), logical(1))
kind <- vapply(
  nodes,
  function(n) n$metadata$kind %||% n$file_type %||% "unknown",
  character(1)
)
communities <- unique(vapply(nodes, function(n) as.character(n$community %||% NA), character(1)))
communities <- communities[!is.na(communities)]

src <- vapply(links, link_end, character(1), which = "source")
tgt <- vapply(links, link_end, character(1), which = "target")
rel <- vapply(links, link_rel, character(1))

# --- structural invariants: any non-zero value is a defect --------------------
dangling <- sum(!src %in% ids | !tgt %in% ids)
self_loops <- sum(src == tgt)
duplicate_edges <- sum(duplicated(paste(src, tgt, rel, sep = "\r")))
collapsed <- sum(duplicated(paste(pmin(src, tgt), pmax(src, tgt), sep = "\r")))

# --- import cycles (Tarjan SCC on the import subgraph) ------------------------
import_edges <- which(rel %in% c("imports", "imports_from"))
imp <- unique(c(src[import_edges], tgt[import_edges]))
adj <- stats::setNames(vector("list", length(imp)), imp)
for (i in import_edges) {
  adj[[src[[i]]]] <- unique(c(adj[[src[[i]]]], tgt[[i]]))
}

tarjan_scc <- function(adj) {
  counter <- 0L
  index <- list()
  low <- list()
  on_stack <- list()
  stack <- character()
  out <- list()
  connect <- function(v) {
    counter <<- counter + 1L
    index[[v]] <<- counter
    low[[v]] <<- counter
    stack <<- c(stack, v)
    on_stack[[v]] <<- TRUE
    for (w in adj[[v]]) {
      if (is.null(index[[w]])) {
        connect(w)
        low[[v]] <<- min(low[[v]], low[[w]])
      } else if (isTRUE(on_stack[[w]])) {
        low[[v]] <<- min(low[[v]], index[[w]])
      }
    }
    if (low[[v]] == index[[v]]) {
      comp <- character()
      repeat {
        w <- stack[[length(stack)]]
        stack <<- stack[-length(stack)]
        on_stack[[w]] <<- FALSE
        comp <- c(comp, w)
        if (identical(w, v)) break
      }
      out[[length(out) + 1L]] <<- comp
    }
  }
  for (v in names(adj)) if (is.null(index[[v]])) connect(v)
  out
}

scc <- if (length(imp)) tarjan_scc(adj) else list()
import_cycles <- Filter(function(c) {
  length(c) > 1L || any(adj[[c[[1L]]]] == c[[1L]])
}, scc)

# --- descriptive metrics ------------------------------------------------------
god <- sort(table(tgt[rel == "calls"]), decreasing = TRUE)
god <- head(god, 10L)

# --- community labels: named, file-echo, duplicate, missing ------------------
# graphify's default label is the hub node's name, which echoes a filename for
# single-file communities.  That is honest; on a multi-file community it is
# misleading, so the report tracks the drift rather than failing on it.
label_tbl <- if (file.exists(LABELS)) {
  jsonlite::fromJSON(LABELS, simplifyVector = TRUE)
} else {
  character()
}
dup_labels <- names(which(table(label_tbl[nzchar(label_tbl)]) > 1L))
is_file_echo <- function(x) {
  grepl(" + ", x, fixed = TRUE) |
    grepl("(^|\\.)(R|r|md|py|yml|yaml|js|json|csv|tsv)$", x)
}
node_comm <- vapply(nodes, function(n) as.character(n$community %||% NA), character(1))
node_file <- vapply(nodes, function(n) n$source_file %||% "", character(1))

# the one known self-loop is inside the vendored ELK layout script; any other
# self-loop is a real extraction defect, so it is gated (the count alone would
# hide a swap: one leaves, another enters, the total still reads 1)
self_loop_file <- node_file[match(src[src == tgt], ids)]
unexpected_self_loops <- sum(!grepl("\\.tools/", self_loop_file))
comm_nodes <- table(node_comm)[communities]
comm_files <- vapply(
  communities,
  function(c) {
    f <- node_file[node_comm %in% c & nzchar(node_file)]
    length(unique(f))
  },
  integer(1)
)
label_name <- vapply(communities, function(c) {
  x <- label_tbl[[c]]
  if (is.null(x) || !nzchar(trimws(x))) "" else x
}, character(1))
label_status <- ifelse(
  !nzchar(label_name), "missing",
  ifelse(
    label_name %in% dup_labels, "duplicate",
    ifelse(is_file_echo(label_name) & comm_files >= 2L, "echo-stale",
      ifelse(is_file_echo(label_name), "echo-honest", "named")
    )
  )
)
label_status <- factor(label_status, levels = c(
  "missing", "duplicate", "echo-stale", "echo-honest", "named"
))
label_counts <- table(label_status)
label_actionable <- sum(label_counts[c("missing", "duplicate", "echo-stale")])

passes <- if (file.exists(ORIGIN)) {
  jsonlite::fromJSON(ORIGIN, simplifyVector = FALSE)$passes %||% list()
} else {
  list()
}

# --- generated report ---------------------------------------------------------
fmt <- function(x) formatC(x, big.mark = " ", format = "d")
rows <- function(...) paste0("| ", paste(c(...), collapse = " | "), " |")
report <- c(
  "# Graphify health (generated)",
  "",
  "Run `Rscript dev/graphify-health.R --write` after `graphify update .`.",
  "Do not edit by hand.",
  "",
  LIMITATION,
  "",
  "## Structural invariants",
  "",
  "These must stay at zero. `dev/graphify-health.R --check` fails if one moves.",
  "",
  rows("metric", "value"),
  rows("---", "---:"),
  rows("dangling endpoints", fmt(dangling)),
  rows("exact duplicate edges", fmt(duplicate_edges)),
  rows("same-endpoint (collapsed) edges", fmt(collapsed)),
  rows("import cycles", fmt(length(import_cycles))),
  rows("self-loops outside `.tools/`", fmt(unexpected_self_loops)),
  "",
  "## Known, accepted",
  "",
  "Not zero, and not gated: the cause is understood and lives outside this package.",
  "",
  rows("metric", "value", "status"),
  rows("---", "---:", "---"),
  rows(
    "self-loops", fmt(self_loops),
    paste(
      "accepted: one recursive `collectEdges()` self-edge in the vendored",
      "`.tools/elk/elk-layout.js`; `--check` fails on any self-loop outside `.tools/`"
    )
  ),
  "",
  "## Community labels",
  "",
  paste(
    "A label that echoes a filename is honest for a single-file community and",
    "misleading for a multi-file one. `--check` warns (never fails) when any",
    "label is `missing`, `duplicate` or `echo-stale`."
  ),
  "",
  rows("status", "communities"),
  rows("---", "---:"),
  vapply(
    names(label_counts),
    function(s) rows(paste0("`", s, "`"), fmt(as.integer(label_counts[[s]]))),
    character(1)
  ),
  "",
  if (label_actionable > 0L) {
    c(
      "### Actionable",
      "",
      rows("cid", "nodes", "files", "label"),
      rows("---:", "---:", "---:", "---"),
      {
        idx <- which(label_status %in% c("missing", "duplicate", "echo-stale"))
        idx <- idx[order(label_status[idx], -as.integer(comm_nodes[idx]))]
        vapply(idx, function(i) {
          rows(
            communities[[i]], fmt(as.integer(comm_nodes[[i]])),
            fmt(comm_files[[i]]), paste0("`", label_name[[i]], "`")
          )
        }, character(1))
      }
    )
  } else {
    "No missing, duplicate or multi-file file-echo labels."
  },
  "",
  "### Relabel provenance",
  "",
  if (!length(passes)) {
    "No relabel pass recorded (`.graphify_labels.origin.json` absent)."
  } else {
    c(
      rows("date", "model", "requested", "relabelled", "kept", "note"),
      rows("---", "---", "---:", "---:", "---:", "---"),
      vapply(passes, function(p) {
        rows(
          p$date %||% "?",
          paste0("`", p$model %||% "?", "`"),
          fmt(as.integer(p$requested %||% 0L)),
          fmt(as.integer(p$relabelled %||% 0L)),
          fmt(as.integer(p$kept_previous %||% 0L)),
          p$note %||% ""
        )
      }, character(1))
    )
  },
  "",
  "## Size",
  "",
  rows("metric", "value"),
  rows("---", "---:"),
  rows("nodes", fmt(length(nodes))),
  rows("edges", fmt(length(links))),
  rows("communities", fmt(length(communities))),
  rows("nodes without `source_file`", fmt(sum(!has_source))),
  "",
  "## Nodes without a `source_file`",
  "",
  "Two kinds, both expected. `external_function` is the AST record of a call",
  "into another package: it has no file here, and `explain` on one shows its",
  "callers, not its definition. `package_function` is published deliberately by",
  "`dev/graphify-handrolled.R` for a hand-roll's replacement; when the same",
  "`pkg::fn` already exists as an AST node the injection reuses that node, so",
  "`explain \"pkg::fn\"` resolves to one node and shows the `supersedes` edge.",
  "The total is a drift signal, not a defect.",
  "",
  rows("kind", "count"),
  rows("---", "---:"),
  {
    kind_counts <- sort(table(kind[!has_source]), decreasing = TRUE)
    vapply(
      names(kind_counts),
      function(k) rows(k, fmt(as.integer(kind_counts[[k]]))),
      character(1)
    )
  },
  "",
  "## God nodes",
  "",
  "Defined here as the top targets by **in-degree over `calls`**: how often the",
  "rest of the code calls this. (Graphify's own ranking mixes document nodes and",
  "is not reproduced.)",
  "",
  rows("node", "callers"),
  rows("---", "---:"),
  vapply(
    names(god),
    function(n) rows(paste0("`", labels[match(n, ids)], "`"), fmt(as.integer(god[[n]]))),
    character(1)
  ),
  "",
  "## Known gaps",
  "",
  "Adopted from GWASILIY's graph practice.",
  "",
  rows("gap", "status", "note"),
  rows("---", "---", "---"),
  rows(
    "hand-roll ledger", "data",
    "`dev/handrolled.yml` holds the decisions; `dev/check-handrolled.R` reads it"
  ),
  rows(
    "`supersedes` graph edges", "adopted",
    paste(
      "`dev/graphify-handrolled.R` injects one `supersedes` edge per entry,",
      "reusing the AST `pkg::fn()` node when it exists and creating a",
      "`package::function` node otherwise; idempotent, and `dev/check.R` runs it"
    )
  ),
  rows(
    "default label prompt on rebuild", "guarded",
    paste(
      "`graphify update .` re-runs the weak built-in prompt on changed",
      "communities; the strict pass and the origin ledger record and repair it"
    )
  ),
  rows(
    "`dev/relabel-bonsai.R` helper", "shipped",
    paste(
      "repeatable strict pass; defaults to missing/duplicate/echo-stale,",
      "`--all-echo` also renames the honest single-file echoes, `--cids=`",
      "repairs one name, and a proposal equal to the current label is kept,",
      "not counted as a rename. A rebuild that renumbers communities drops",
      "label keys; the pass restores them as `missing` and names them"
    )
  ),
  rows(
    "vignettes/ and man/ coverage", "gap",
    paste(
      "graphify does not index `.Rmd` or `.Rd`, so `vignettes/*.Rmd` and",
      "`man/*.Rd` never enter the graph; only the prose in `docs/`, `README.md`,",
      "`NEWS.md`, `ROADMAP.md` and `VALIDATION.md` can. Export or convert them",
      "first if the manuals must be graph concepts"
    )
  ),
  rows(
    "single-file echo labels (254)", "skipped",
    paste(
      "a file-echo name is honest for a one-file community; the strict pass",
      "renames only missing/duplicate/echo-stale, `--all-echo` opts in"
    )
  ),
  "",
  "## Configuration levers",
  "",
  paste(
    "What is wired into graphify here, and what is deliberately not. `open` is",
    "usable in the installed build but not adopted; `skipped` has no consumer."
  ),
  "",
  rows("lever", "status", "note"),
  rows("---", "---", "---"),
  rows(
    "`graphify update .`", "shipped",
    "code-only re-extract, no LLM; first step of Rerun below"
  ),
  rows(
    "strict label pass + origin ledger", "shipped",
    "`dev/relabel-bonsai.R`; all communities named, every pass recorded"
  ),
  rows(
    "`supersedes` edges from the hand-roll ledger", "shipped",
    "see Known gaps; answers \"what replaced this hand-roll\""
  ),
  rows(
    "local LLM semantic backend", "open (docs)",
    paste(
      "`dev/graphify-semantic.sh` holds the safe flags; `extract --backend=openai`",
      "reads `OPENAI_BASE_URL`/`OPENAI_MODEL` (bonsai). Measured 2026-11-06:",
      "llama-server -c 65536 -np 2 gives 32768 tokens per slot and 2 slots,",
      "~32 output tokens/s, ~3000 tokens per document; `--token-budget 6000`",
      "and `--max-concurrency 2` fit. A full 54-document pass is still hours,",
      "so only `docs` is opted in; see the extraction result below"
    )
  ),
  rows(
    "`export wiki`", "shipped",
    paste(
      "411 articles in gitignored `graphify-out/wiki/`, one per community plus an",
      "index; template-only, no LLM, so its quality is the curated labels"
    )
  ),
  rows(
    "git hook (`graphify hook install`)", "shipped",
    paste(
      "post-commit/checkout/merge/rewrite plus the graph merge driver; it runs",
      "`update` only (no LLM), so `dev/graphify-health.R --write` stays manual"
    )
  ),
  rows(
    "`.graphifyignore`", "shipped",
    paste(
      "`data-raw/proposed/`, `results/`, `.tmp/`, `.tools/` are design scratch or",
      "vendored tooling, not package API; a rebuild prunes them, and the file can",
      "only ever add exclusions"
    )
  ),
  rows(
    "`--mode deep`", "skipped",
    paste(
      "richer INFERRED edges, but a full re-extract and it re-names changed",
      "communities with the weak built-in prompt; no consumer here"
    )
  ),
  rows(
    "`query` / `path` / `explain` / `affected`", "shipped",
    paste(
      "verified from the shell: `explain \"Biostrings::complement\"` resolves one",
      "node and shows its `supersedes` edge, `affected \"normalise_allele()\"`",
      "lists its caller. 0.9.92 has no MCP server, so agents use this CLI"
    )
  ),
  rows(
    "`export neo4j` / `graphml` / `obsidian` / `svg`", "skipped",
    "no consumer in this repository"
  ),
  "",
  "## Provenance",
  "",
  paste0("Built at commit `", g$built_at_commit %||% "unknown", "`.")
)

head_commit <- suppressWarnings(system2(
  "git", c("rev-parse", "HEAD"),
  stdout = TRUE, stderr = FALSE
))
if (length(head_commit) == 1L) {
  report <- c(
    report,
    if (identical(substr(head_commit, 1L, nchar(g$built_at_commit %||% "")), g$built_at_commit %||% "")) {
      "The graph matches `HEAD`."
    } else {
      paste0(
        "`HEAD` is `", head_commit, "`, so the graph is **stale**; rerun ",
        "`graphify update .` and this script."
      )
    }
  )
}

report <- c(report, "", "## Rerun", "", "```sh", "graphify update .",
            "Rscript dev/graphify-handrolled.R",
            "Rscript dev/graphify-health.R --write",
                      "graphify export wiki",
            "graphify export html",
            "# when a label warns: Rscript dev/relabel-bonsai.R [--all-echo]",
            "# one wrong name: Rscript dev/relabel-bonsai.R --cids=N",
            "# semantic pass (opt-in, local LLM): dev/graphify-semantic.sh",
            "# refresh the analysis sidecar: graphify cluster-only . --no-label",
            "```", "",
            "Do not pass `--no-cluster` to `update`: it drops the community",
            "assignment and duplicates some edges, so the health numbers above",
            "mislead. `cluster-only` is deterministic on an unchanged graph but can",
            "renumber after a rebuild, so re-check the labels and rerun",
            "`dev/relabel-bonsai.R` if it warns.")

if (identical(mode, "write")) {
  dir.create(dirname(REPORT), showWarnings = FALSE, recursive = TRUE)
  writeLines(report, REPORT, useBytes = TRUE)
  message("wrote ", REPORT, " (", length(nodes), " nodes, ", length(links), " edges)")
  quit(status = 0L)
}

# --- --check: regressions fail; drift only warns ------------------------------
regressions <- c(
  if (dangling > 0L) paste0(dangling, " dangling endpoint(s)"),
  if (duplicate_edges > 0L) paste0(duplicate_edges, " exact duplicate edge(s)"),
  if (collapsed > 0L) paste0(collapsed, " same-endpoint edge(s)"),
  if (unexpected_self_loops > 0L) {
    paste0(unexpected_self_loops, " self-loop(s) outside .tools/")
  },
  if (length(import_cycles)) paste0(length(import_cycles), " import cycle(s)")
)
if (length(regressions)) {
  message("FAIL graphify health: ", paste(regressions, collapse = "; "))
  quit(status = 1L)
}
if (label_actionable > 0L) {
  message(
    "WARN ", label_actionable, " community label(s) to curate ",
    "(missing/duplicate/echo-stale); see ", REPORT
  )
}
if (!file.exists(REPORT)) {
  message("WARN no ", REPORT, " - run `Rscript dev/graphify-health.R --write`")
} else if (!any(grepl(LIMITATION, readLines(REPORT, warn = FALSE), fixed = TRUE))) {
  message("WARN ", REPORT, " does not name the limitation; rerun --write")
}
message("ok graphify health: no dangling, duplicate, collapsed or cyclic edges")
