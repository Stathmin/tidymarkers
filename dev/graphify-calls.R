#!/usr/bin/env Rscript
# =============================================================================
#  dev/graphify-calls.R  --  add cross-file `calls` edges, deterministically
# =============================================================================
#  The AST scanner never links a call in one file to a definition in another,
#  so every community is a file island and `explain`/`affected` under-report.
#  This script reparses each R/*.R file, finds every call site, attributes it
#  to the enclosing function, and emits a `calls` edge when that function calls
#  a *different* file's function.  No LLM; pure parse + name resolution.
#
#  Idempotent: links from a previous run carry `_origin = "calls-pass"` and are
#  dropped before the new ones are added.  This mutates a gitignored artifact
#  (graphify-out/ is a developer tool), exactly as dev/graphify-handrolled.R does.
#
#  Self-skips (exit 0) when there is no graph or jsonlite: CI has no graph.
#
#  Run:  Rscript dev/graphify-calls.R
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

if (!file.exists(GRAPH) || !requireNamespace("jsonlite", quietly = TRUE)) {
  message("skip (no ", GRAPH, " or jsonlite not installed - graphify-out/ is a developer tool)")
  quit(status = 0L)
}

g <- jsonlite::fromJSON(GRAPH, simplifyVector = FALSE)
nodes <- g$nodes %||% list()
links <- g$links %||% list()

# drop a previous injection so the operation is idempotent
keep_node <- vapply(nodes, function(n) !identical(n[["_origin"]], "calls-pass"), logical(1))
keep_link <- vapply(links, function(e) !identical(e[["_origin"]], "calls-pass"), logical(1))
nodes <- nodes[keep_node]
links <- links[keep_link]
nodes_before <- length(nodes)

ids <- vapply(nodes, function(n) n$id %||% "", character(1))
labels <- vapply(nodes, function(n) n$label %||% "", character(1))
srcs <- vapply(nodes, function(n) n$source_file %||% "", character(1))
kinds <- vapply(
  nodes,
  function(n) n$metadata$kind %||% n$file_type %||% "unknown",
  character(1)
)

# --- defs: unambiguous R/ function name -> (id, source_file) ----------------
# graphify labels a function node `name()`.  A name defined in more than one
# R/ file is ambiguous and never resolves.
is_fn <- kinds == "function" & grepl("^R/", srcs)
fn_name <- sub("\\(\\)$", "", labels[is_fn])
fn_id <- ids[is_fn]
fn_src <- srcs[is_fn]
name_tbl <- split(seq_along(fn_name), fn_name)
amb_names <- names(name_tbl)[
  vapply(name_tbl, function(idx) length(unique(fn_src[idx])) > 1L, logical(1))
]
ok <- !(fn_name %in% amb_names)
defs <- list()
for (i in which(ok)) {
  defs[[fn_name[[i]]]] <- list(id = fn_id[[i]], source_file = fn_src[[i]])
}
# (file, name) -> node id, for the *enclosing* function's node
file_slot <- vapply(seq_along(fn_name), function(i) fn_id[[i]], character(1))
names(file_slot) <- paste0(fn_src, "\t", fn_name)

# --- existing edges, so we never emit a duplicate ----------------------------
existing <- paste(
  vapply(links, function(e) e$source %||% "", character(1)),
  vapply(links, function(e) e$target %||% "", character(1)),
  vapply(links, function(e) e$relation %||% "", character(1)),
  sep = "\t"
)

new_links <- list()
target_deg <- character()
same_file <- 0L
shadowed <- 0L
unresolved_enclosing <- 0L
already_present <- 0L

is_fn_def <- function(x) {
  if (!is.call(x)) {
    return(FALSE)
  }
  h <- as.character(x[[1]])[[1]]
  if (h %in% c("<-", "=", "<<-") && length(x) >= 3L) {
    return(is.symbol(x[[2]]) && is.call(x[[3]]) &&
      identical(as.character(x[[3]][[1]])[[1]], "function"))
  }
  FALSE
}

for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  parsed <- tryCatch(parse(f, keep.source = TRUE), error = function(.) NULL)
  if (is.null(parsed)) next
  pd <- tryCatch(getParseData(parsed), error = function(.) NULL)
  if (is.null(pd) || !nrow(pd)) next
  srs <- attr(parsed, "srcref")

  # top-level function definitions with their line ranges; other top-level code
  # has no enclosing function and its calls are skipped below.
  scopes <- list()
  local_names <- character()
  for (i in seq_along(parsed)) {
    if (!is_fn_def(parsed[[i]])) next
    nm <- as.character(parsed[[i]][[2]])[[1]]
    local_names <- c(local_names, nm)
    if (is.null(srs) || length(srs) < i) next
    sr <- srs[[i]]
    scopes[[length(scopes) + 1L]] <- list(name = nm, first = sr[[1]], last = sr[[3]])
  }

  # every function call site in the file, with its line.  getParseData is the
  # parse tree in tabular form; SYMBOL_FUNCTION_CALL is exactly the called name.
  sfc <- pd[pd$token == "SYMBOL_FUNCTION_CALL", , drop = FALSE]
  if (!nrow(sfc)) next
  # `pkg::fn` records the called name under a `::` (NS_GET) parent: already in
  # the graph as an external node, so skip it here.
  ns_parents <- unique(pd$parent[pd$token == "NS_GET"])
  sfc <- sfc[!(sfc$parent %in% ns_parents), , drop = FALSE]
  # operators and infix/backtick names fail the identifier test; 1-2 char names
  # are dropped per spec.
  txt <- sfc$text
  keep <- make.names(txt) == txt & nchar(txt) > 2L
  sfc <- sfc[keep, , drop = FALSE]
  if (!nrow(sfc)) next

  for (k in seq_len(nrow(sfc))) {
    called <- sfc$text[[k]]
    call_line <- sfc$line1[[k]]

    # innermost enclosing function by line range
    hit <- vapply(scopes, function(s) s$first <= call_line && call_line <= s$last, logical(1))
    if (!any(hit)) next
    spans <- vapply(scopes[hit], function(s) s$last - s$first, numeric(1))
    enc <- scopes[[which(hit)[[which.min(spans)]]]]

    d <- defs[[called]]
    if (is.null(d)) next
    if (identical(d$source_file, f)) {
      same_file <- same_file + 1L
      next
    }
    if (called %in% local_names) {
      shadowed <- shadowed + 1L
      next
    }
    src_id <- unname(file_slot[paste0(f, "\t", enc$name)])
    if (!length(src_id) || is.na(src_id) || !nzchar(src_id)) {
      unresolved_enclosing <- unresolved_enclosing + 1L
      next
    }
    key <- paste(src_id, d$id, "calls", sep = "\t")
    if (key %in% existing) {
      already_present <- already_present + 1L
      next
    }
    new_links[[length(new_links) + 1L]] <- list(
      source = src_id,
      target = d$id,
      relation = "calls",
      confidence = 0.6,
      source_file = f,
      source_location = paste0("line ", call_line),
      weight = 1,
      `_origin` = "calls-pass"
    )
    existing <- c(existing, key)
    target_deg <- c(target_deg, d$id)
  }
}

edges_added <- length(new_links)

if (edges_added > 5000L) {
  cat(
    "  ABORT calls-pass: ", edges_added, " edges would be added (limit 5000)\n",
    "  ambiguous names: ", length(amb_names), "; shadowed: ", shadowed,
    "; same-file: ", same_file, "; unresolved enclosing: ", unresolved_enclosing,
    "\n",
    sep = ""
  )
  quit(status = 1L)
}

g$nodes <- nodes
g$links <- c(links, new_links)
jsonlite::write_json(g, GRAPH, auto_unbox = TRUE, pretty = TRUE, null = "null")

cat(
  "  ok   calls-pass: ", edges_added, " cross-file call edge(s) added\n",
  "  nodes: ", nodes_before, " -> ", length(nodes), "\n",
  "  ambiguous names skipped: ", length(amb_names), "\n",
  "  shadowed skipped: ", shadowed, "\n",
  "  same-file skipped: ", same_file, "\n",
  "  unresolved enclosing skipped: ", unresolved_enclosing, "\n",
  "  already present: ", already_present, "\n",
  sep = ""
)
if (length(target_deg)) {
  lbl_of <- setNames(labels, ids)
  tbl <- sort(table(factor(target_deg, levels = unique(target_deg))), decreasing = TRUE)
  top <- head(tbl, 10L)
  cat("  top targets by added in-degree:\n")
  for (i in seq_along(top)) {
    cat("    ", sprintf("%4d", as.integer(top[[i]])), "  ",
      lbl_of[[names(top)[[i]]]], "\n",
      sep = ""
    )
  }
}
