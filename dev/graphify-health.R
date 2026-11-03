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
  rows("self-loops", fmt(self_loops)),
  rows("exact duplicate edges", fmt(duplicate_edges)),
  rows("same-endpoint (collapsed) edges", fmt(collapsed)),
  rows("import cycles", fmt(length(import_cycles))),
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
  "Symbols from other packages are recorded as nodes but have no file in this",
  "repository. They are expected; the count is a drift signal, not a defect.",
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
            "Rscript dev/graphify-health.R --write", "```", "",
            "Do not pass `--no-cluster`: it drops the community assignment and",
            "duplicates some edges, so the health numbers above mislead.")

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
  if (length(import_cycles)) paste0(length(import_cycles), " import cycle(s)")
)
if (length(regressions)) {
  message("FAIL graphify health: ", paste(regressions, collapse = "; "))
  quit(status = 1L)
}
if (!file.exists(REPORT)) {
  message("WARN no ", REPORT, " - run `Rscript dev/graphify-health.R --write`")
} else if (!any(grepl(LIMITATION, readLines(REPORT, warn = FALSE), fixed = TRUE))) {
  message("WARN ", REPORT, " does not name the limitation; rerun --write")
}
message("ok graphify health: no dangling, duplicate, collapsed or cyclic edges")
