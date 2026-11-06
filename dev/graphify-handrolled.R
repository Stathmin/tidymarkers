#!/usr/bin/env Rscript
# =============================================================================
#  dev/graphify-handrolled.R  --  show what replaced a hand-roll
# =============================================================================
#  Reads the ledger in dev/handrolled.yml and writes two things into
#  graphify-out/graph.json for every row with a `replacement`:
#    * a `package::function` node (file_type "package", no source_file), and
#    * a `supersedes` edge from that node to the hand-roll's `symbol()` node.
#  The graph then answers "what replaced this hand-roll" and a community name
#  like "Primer Thermodynamics QC" can point at its replacement.
#
#  Idempotent: rows from a previous run carry `_origin = "handrolled"` and are
#  dropped before the new ones are added.  This mutates a gitignored artifact
#  (graphify-out/ is a developer tool), exactly as GWASILIY's graph step does.
#
#  Self-skips (exit 0) when there is no graph or ledger: CI has no graph.
#
#  Run:  Rscript dev/graphify-handrolled.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

GRAPH <- "graphify-out/graph.json"
LEDGER <- "dev/handrolled.yml"

if (!file.exists(GRAPH) || !file.exists(LEDGER)) {
  message("skip (no ", GRAPH, " or ", LEDGER, " - graphify-out/ is a developer tool)")
  quit(status = 0L)
}
for (pkg in c("jsonlite", "yaml")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("skip (", pkg, " not installed)")
    quit(status = 0L)
  }
}

registry <- yaml::read_yaml(LEDGER)$entries
g <- jsonlite::fromJSON(GRAPH, simplifyVector = FALSE)
nodes <- g$nodes %||% list()
links <- g$links %||% list()

# drop a previous injection so the operation is idempotent
keep_node <- vapply(nodes, function(n) !identical(n[["_origin"]], "handrolled"), logical(1))
keep_link <- vapply(links, function(e) !identical(e[["_origin"]], "handrolled"), logical(1))
nodes <- nodes[keep_node]
links <- links[keep_link]

ids <- vapply(nodes, function(n) n$id %||% "", character(1))
labels <- vapply(nodes, function(n) n$label %||% "", character(1))

new_nodes <- list()
new_links <- list()
skipped <- character()
seen <- character()
for (row in registry) {
  replacement <- row$replacement %||% ""
  if (!nzchar(replacement)) next
  target <- ids[match(paste0(row$symbol, "()"), labels)]
  if (!length(target) || is.na(target) || !nzchar(target)) {
    skipped <- c(skipped, row$symbol)
    next
  }
  # The AST already records a called package function as a `<pkg>::<fn>()`
  # node with no source_file.  Reuse it instead of publishing a second,
  # near-identically named node, so `explain "pkg::fn"` is unambiguous and the
  # external-node count does not grow by one per registry entry.
  existing <- ids[match(paste0(replacement, "()"), labels)]
  existing <- existing[!is.na(existing) & nzchar(existing)]
  replacement_id <- if (length(existing)) existing[[1L]] else paste0("package:", replacement)
  if (!replacement_id %in% seen && !replacement_id %in% ids) {
    new_nodes[[length(new_nodes) + 1L]] <- list(
      id = replacement_id,
      label = replacement,
      file_type = "package",
      # jsonlite turns NA into `{}` under auto_unbox and graphify's incremental
      # scanner then fails on the node; empty strings stay strings.
      source_file = "",
      source_location = "",
      metadata = list(kind = "package_function", status = row$status %||% ""),
      `_origin` = "handrolled"
    )
    seen <- c(seen, replacement_id)
  }
  new_links[[length(new_links) + 1L]] <- list(
    source = replacement_id,
    target = target,
    relation = "supersedes",
    confidence = 1,
    source_file = "",
    source_location = "",
    weight = 1,
    `_origin` = "handrolled"
  )
}

g$nodes <- c(nodes, new_nodes)
g$links <- c(links, new_links)
jsonlite::write_json(g, GRAPH, auto_unbox = TRUE, pretty = TRUE, null = "null")

cat(
  "  ok   ", length(new_links), " supersedes edge(s) injected; ",
  length(skipped), " symbol(s) without a graph node\n",
  sep = ""
)
if (length(skipped)) {
  cat("  skipped: ", paste(skipped, collapse = ", "), "\n", sep = "")
}
