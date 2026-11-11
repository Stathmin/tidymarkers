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

# --- semantic (LLM) coverage: only dev/graphify-semantic.sh writes a hash ---
# The manifest records one `semantic_hash` per file; an empty string means the
# AST pass saw it but the model did not.  Code files never get one, so the
# denominator is the document-like set that graphify can send to the backend.
MANIFEST <- "graphify-out/manifest.json"
sem_done <- 0L
sem_total <- 0L
if (file.exists(MANIFEST)) {
  mf <- jsonlite::fromJSON(MANIFEST, simplifyVector = FALSE)
  sem_hash <- vapply(mf, function(x) x$semantic_hash %||% "", character(1))
  doclike <- grepl(
    "\\.(md|mdx|qmd|rmd|txt|rst|html|ya?ml|j2|jinja2?|pdf|png|jpe?g|gif|svg|webp)$",
    names(mf),
    ignore.case = TRUE
  )
  sem_total <- sum(doclike)
  sem_done <- sum(nzchar(sem_hash[doclike]))
}
# Provenance matters: the AST pass emits INFERRED `references` edges, so the
# INFERRED total is not evidence of the model. Count what the pass actually
# minted (`_origin == "semantic"`) and what the AST made separately.
origin_of <- function(x) x[["_origin"]] %||% ""
is_semantic <- function(x) identical(origin_of(x), "semantic")
concept_nodes <- sum(kind == "concept")
# A `concept` node with no source_file is not a model-made idea: it is an AST
# `import` target (Python `argparse`, `json`, ...) or a raw file-id label. Keep
# the two apart so the semantic table does not credit the LLM with them.
concept_unlocated <- sum(kind == "concept" & !has_source)
sem_nodes <- sum(vapply(nodes, is_semantic, logical(1)))
inferred_edges <- sum(vapply(
  links, function(e) identical(e$confidence %||% "", "INFERRED"), logical(1)
))
inferred_semantic <- sum(vapply(
  links,
  function(e) identical(e$confidence %||% "", "INFERRED") && is_semantic(e),
  logical(1)
))
inferred_ast <- inferred_edges - inferred_semantic

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

# --- granularity: what the AST pass makes a node of ---------------------------
# Code nodes are function- and top-level-assignment-level, never file-level for
# package code; the question is how many are leaves reachable only through their
# file's `contains` edge.  A leaf adds no seam, so the leaf ratio is the
# entity-multiplication signal, not the raw node count.
node_file <- vapply(nodes, function(n) n$source_file %||% "", character(1))
deg <- table(c(src, tgt))
node_deg <- as.integer(deg[ids])
node_deg[is.na(node_deg)] <- 0L
is_leaf <- node_deg <= 1L
kind_tbl <- table(kind)
leaf_kind <- vapply(
  c("function", "variable", "file"),
  function(k) sum(is_leaf & kind == k), integer(1)
)
dupl_counts <- table(labels)
dupl_labels <- sum(dupl_counts > 1L)
dupl_nodes <- sum(dupl_counts[dupl_counts > 1L])
pkg_share <- sum(grepl("^R/", node_file))
script_share <- sum(grepl("^(data-raw|tests|dev)/", node_file))

# --- nesting: top-level binding or function local? ---------------------------
# `source_location` is `L<n>`; a line indented from column 1 is a local (a nested
# closure, or a `test_that()`/`local()` body).  Only `R/` is written
# top-level-first, so the split is reported per directory: it shows where the
# entity multiplication actually lives.
node_loc <- vapply(nodes, function(n) n$source_location %||% "", character(1))
line_cache <- list()
indent_at <- function(file, loc) {
  ln <- suppressWarnings(as.integer(sub("^L", "", loc)))
  if (is.na(ln) || !nzchar(file) || !file.exists(file)) return(NA_integer_)
  if (is.null(line_cache[[file]])) line_cache[[file]] <<- readLines(file, warn = FALSE)
  lines <- line_cache[[file]]
  if (ln < 1L || ln > length(lines) || !nzchar(lines[[ln]])) return(0L)
  x <- lines[[ln]]
  nchar(x) - nchar(sub("^[ \t]+", "", x))
}
node_indent <- mapply(indent_at, node_file, node_loc, USE.NAMES = FALSE)
located <- sum(!is.na(node_indent))
is_nested <- !is.na(node_indent) & node_indent > 0L
top_dir <- sub("/.*$", "", node_file)
nested_total <- sum(is_nested)
fn_top <- sum(!is_nested & has_source & kind == "function")
fn_nested <- sum(is_nested & kind == "function")
var_top <- sum(!is_nested & has_source & kind == "variable")
var_nested <- sum(is_nested & kind == "variable")
var_nested_tests <- sum(is_nested & kind == "variable" & top_dir == "tests")
fn_nested_r <- sum(is_nested & kind == "function" & top_dir == "R")
var_nested_r <- sum(is_nested & kind == "variable" & top_dir == "R")
nested_by_dir <- vapply(
  c("tests", "data-raw", "dev", "inst", "R"),
  function(d) sum(is_nested & top_dir == d),
  integer(1)
)
r_files <- list.files("R", pattern = "[.]R$")
r_fns <- vapply(
  r_files,
  function(f) sum(kind == "function" & node_file == paste0("R/", f)),
  integer(1)
)
r_lines <- vapply(
  r_files,
  function(f) length(readLines(file.path("R", f), warn = FALSE)),
  integer(1)
)
r_size <- order(-r_lines)[seq_len(min(9L, length(r_files)))]

# --- interfaces: exports and S3 dispatch --------------------------------------
# graphify's R extractor emits one node per exported function and per S3 method,
# but no `generic -> method` edge and no exported/internal flag, so the graph
# answers "who calls what" but not "what implements this generic".
ns <- if (file.exists("NAMESPACE")) readLines("NAMESPACE", warn = FALSE) else character()
ns_arg <- function(line) {
  sub("^[^(]+\\((.*)\\)$", "\\1", line)
}
exp_names <- unlist(lapply(ns_arg(grep("^export\\(", ns, value = TRUE)), function(s) {
  trimws(strsplit(s, ",", fixed = TRUE)[[1L]])
}))
s3_parts <- lapply(ns_arg(grep("^S3method\\(", ns, value = TRUE)), function(s) {
  p <- trimws(gsub('"', "", strsplit(s, ",", fixed = TRUE)[[1L]]))
  # S3method(vctrs::vec_cast, ...) registers a method on a foreign generic;
  # graphify's node label omits the `pkg::` prefix, so strip it here.
  p <- sub("^[^:]+::", "", p)
  if (length(p) == 2L) p else NULL
})
s3_parts <- Filter(Negate(is.null), s3_parts)
fn_bare <- sub("\\(\\)$", "", labels[grepl("\\(\\)$", labels)])
s3_labels <- vapply(s3_parts, function(p) paste0(p[[1L]], ".", p[[2L]], "()"), character(1))
generics <- unique(vapply(s3_parts, function(p) paste0(p[[1L]], "()"), character(1)))
n_exp_nodes <- sum(exp_names %in% fn_bare)
n_s3_nodes <- sum(s3_labels %in% labels)
n_generic_nodes <- sum(generics %in% labels)
id_of <- function(x) ids[match(x, labels)]
method_ids <- id_of(s3_labels[s3_labels %in% labels])
generic_ids <- id_of(generics[generics %in% labels])
dispatch_edges <- sum(src %in% generic_ids & tgt %in% method_ids) +
  sum(tgt %in% generic_ids & src %in% method_ids)
pct <- function(n, d) sprintf("%.0f%%", 100 * n / max(d, 1L))
# --- code <-> documentation distance (deterministic) --------------------------
# The graph cannot see man/ or vignettes/, so compare the source directly: per
# export, roxygen @param names against formals().  `@rdname` merges several
# aliases into one Rd topic, so a function's docs and its documented names are
# the union over its topic; a documented name is "phantom" only when no alias in
# the topic has that formal.  This is distance, not faithfulness: whether the
# prose describes the body needs a model.
doc_forms <- list()
doc_recs <- list()
if (dir.exists("R")) {
  for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
    dlines <- readLines(f, warn = FALSE)
    ex <- parse(f, keep.source = TRUE)
    srs <- attr(ex, "srcref")
    for (i in seq_along(ex)) {
      e <- ex[[i]]
      if (!is.call(e) || !as.character(e[[1]]) %in% c("<-", "=")) next
      rhs <- e[[3]]
      if (!is.call(rhs) || !identical(rhs[[1]], as.name("function"))) next
      nm <- as.character(e[[2]])
      fa <- names(as.list(rhs[[2]]))
      doc_forms[[nm]] <- if (is.null(fa)) character(0) else setdiff(fa, "...")
      sr <- srs[[i]]
      if (is.null(sr)) next
      j <- as.integer(sr[[1]]) - 1L
      blk <- character()
      while (j >= 1L && grepl("^#'", dlines[[j]])) {
        blk <- c(dlines[[j]], blk)
        j <- j - 1L
      }
      if (!length(blk)) next
      txt <- sub("^#+'\\s?", "", blk)
      rd <- sub("^@rdname\\s+(\\S+).*$", "\\1", grep("^@rdname\\s+", txt, value = TRUE))
      raw <- sub("^@param\\s+(\\S+).*$", "\\1", grep("^@param\\s+", txt, value = TRUE))
      doc_recs[[nm]] <- list(
        topic = if (length(rd)) rd[[1L]] else nm,
        params = unlist(strsplit(raw, ",", fixed = TRUE), use.names = FALSE),
        return = any(grepl("^@return\\b", txt)),
        examples = any(grepl("^@examples\\b", txt))
      )
    }
  }
}
doc_ok <- length(doc_recs) > 0L
doc_topic <- vapply(doc_recs, function(r) r$topic, character(1))
doc_siblings <- function(nm) names(doc_recs)[doc_topic %in% doc_topic[[nm]]]
doc_merge <- function(nm, what) {
  unique(unlist(lapply(doc_siblings(nm), function(m) doc_recs[[m]][[what]])))
}
doc_any <- function(nm, what) {
  any(vapply(doc_siblings(nm), function(m) isTRUE(doc_recs[[m]][[what]]), logical(1)))
}
doc_missing <- setdiff(exp_names, names(doc_recs))
doc_undoc <- list()
doc_phantom <- list()
doc_no_return <- character()
doc_no_examples <- character()
for (nm in intersect(exp_names, names(doc_recs))) {
  p <- doc_merge(nm, "params")
  u <- setdiff(doc_forms[[nm]], p)
  if (length(u)) doc_undoc[[nm]] <- u
  tf <- unique(unlist(lapply(doc_siblings(nm), function(m) doc_forms[[m]])))
  ph <- setdiff(p, c(tf, "..."))
  if (length(ph)) doc_phantom[[nm]] <- ph
  if (!doc_any(nm, "return")) doc_no_return <- c(doc_no_return, nm)
  if (!doc_any(nm, "examples")) doc_no_examples <- c(doc_no_examples, nm)
}
n_doc_undoc <- sum(lengths(doc_undoc))
n_doc_phantom <- sum(lengths(doc_phantom))


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

# `.tools/` is ignored now, so the vendored ELK self-loop is gone; the gate
# stays on any self-loop outside `.tools/` because a swap in the total alone
# would hide it (one leaves, another enters, the count still reads the same).
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
  if (self_loops > 0L) {
    c(
      "Not zero, and not gated: the cause is understood and lives outside this package.",
      "",
      rows("metric", "value", "status"),
      rows("---", "---:", "---"),
      rows(
        "self-loops", fmt(self_loops),
        paste(
          "accepted: a recursive `collectEdges()` self-edge in the vendored",
          "`.tools/elk/elk-layout.js`; `--check` fails on any self-loop outside `.tools/`"
        )
      )
    )
  } else {
    "None. Every structural invariant above is zero; nothing is accepted as-is."
  },
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
  "Most are expected. `external_function` is the AST record of a call into",
  "another package: it has no file here, and `explain` on one shows its",
  "callers, not its definition. `package_function` is published deliberately by",
  "`dev/graphify-handrolled.R` for a hand-roll's replacement; when the same",
  "`pkg::fn` already exists as an AST node the injection reuses that node, so",
  "`explain \"pkg::fn\"` resolves to one node and shows the `supersedes` edge.",
  "Two kinds are accidental and unexplained by the extractor: the 13 `concept`",
  "nodes are `import`/`source()` targets -- Python modules (`argparse`, `json`,",
  "`urllib_request`) and R `source()` paths rendered as absolute ids",
  "(`root_development_..._xlsx_read_r`) -- with a raw id for a label and no",
  "file, and `external_class` has the same shape.",
  "`explain` on any of them is vague by construction -- there is no definition",
  "to show. The external-symbol total is a drift signal; the concept total should",
  "be zero, and is not.",
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
  "## Granularity and interfaces",
  "",
  paste(
    "Graphify's R extractor emits a node per function (a nested closure",
    "included), per file, and per assignment, top-level or local. The split is",
    "where the noise is:", fmt(fn_nested_r), "nested nodes under `R/` --",
    "all closures (a helper defined inside another function), no nested",
    paste0("assignments (", fmt(var_nested_r), ") -- against"),
    fmt(nested_by_dir[["tests"]]), "inside `tests/` -- `test_that()`",
    "and `local()` bodies mint a node per local -- and",
    fmt(nested_by_dir[["data-raw"]]), "in `data-raw/`. Most such nodes are",
    "leaves reachable only through their file's `contains` edge; a leaf adds no",
    "seam, so the leaf ratio -- not the raw node count -- is the",
    "entity-multiplication signal."
  ),
  "",
  rows("metric", "value"),
  rows("---", "---:"),
  rows("nodes: `function`", fmt(as.integer(kind_tbl[["function"]] %||% 0L))),
  rows("nodes: `variable`", fmt(as.integer(kind_tbl[["variable"]] %||% 0L))),
  rows("nodes: `file`", fmt(as.integer(kind_tbl[["file"]] %||% 0L))),
  rows(
    "degree-1 (leaf) nodes",
    paste0(fmt(sum(is_leaf)), " (", pct(sum(is_leaf), length(nodes)), ")")
  ),
  rows("…of which `variable`", fmt(as.integer(leaf_kind[["variable"]] %||% 0L))),
  rows(
    "nested (indented) nodes",
    paste0(fmt(nested_total), " (", pct(nested_total, located), " of located)")
  ),
  rows(
    "labels shared by >1 node",
    paste0(fmt(dupl_labels), " (", fmt(dupl_nodes), " nodes)")
  ),
  rows(
    "nodes under `R/`",
    paste0(fmt(pkg_share), " (", pct(pkg_share, length(nodes)), ")")
  ),
  rows(
    "nodes under `data-raw/`, `tests/`, `dev/`",
    paste0(fmt(script_share), " (", pct(script_share, length(nodes)), ")")
  ),
  rows(
    "exported functions with a node",
    paste0(fmt(n_exp_nodes), " / ", fmt(length(exp_names)))
  ),
  rows(
    "S3 methods with a node",
    paste0(fmt(n_s3_nodes), " / ", fmt(length(s3_labels)))
  ),
  rows(
    "generics with a node",
    paste0(fmt(n_generic_nodes), " / ", fmt(length(generics)))
  ),
  rows("`generic`/`method` pairs joined by an edge", fmt(dispatch_edges)),
  "",
  paste(
    "Every export (175/175) and all but one S3 method have a node. The absent",
    "method is `+.marker_plan()`: it is defined as `\"+.marker_plan\" <- function`,",
    "and graphify's `binding()` accepts only an identifier left-hand side, so a",
    "string-named assignment is dropped. The 7 generics without a node are",
    "external (`+`, `plot`, `print` and four `vctrs::`): the package implements",
    "their methods but does not define the generic, so no package node exists.",
    "That is honest, not a gap. graphify's R extractor emits no",
    "`generic -> method` relation and no exported/internal flag",
    "(`interface_dispatch.py` covers C#, Java, Swift and other languages, not R),",
    "so the single pair counted above is a plain `calls` edge",
    "(`marker_primers.marker_lax()` calls `marker_primers()` in its body), not",
    "dispatch. \"What implements `marker_chunk()`\" is a label-prefix query, not a",
    "traversal. The duplicate labels are not duplicate entities: two `ext()`",
    "nodes are two files' helpers, but `explain ext` is ambiguous because R has",
    "no namespaces."
  ),
  "",
  "### Function / variable nodes: top-level vs nested",
  "",
  rows("kind", "top-level", "nested"),
  rows("---", "---:", "---:"),
  rows("`function`", fmt(fn_top), fmt(fn_nested)),
  rows("`variable`", fmt(var_top), fmt(var_nested)),
  "",
  "### Nested nodes by directory",
  "",
  rows("dir", "nested nodes"),
  rows("---", "---:"),
  vapply(
    names(nested_by_dir),
    function(d) rows(paste0("`", d, "/`"), fmt(as.integer(nested_by_dir[[d]]))),
    character(1)
  ),
  "",
  "### `R/` files by size",
  "",
  "Largest files first; `functions` is the `function` node count for the file.",
  "",
  rows("file", "functions", "lines"),
  rows("---", "---:", "---:"),
  vapply(
    r_size,
    function(i) rows(paste0("`R/", r_files[[i]], "`"), fmt(r_fns[[i]]), fmt(r_lines[[i]])),
    character(1)
  ),
  "",
  "## Code ↔ documentation",
  "",
  paste(
    "The graph cannot see `man/`, so this reads the source: for",
    "each export it compares roxygen `@param` names with `formals()`, resolving",
    "`@rdname` topics (one Rd block documents several aliases). It is a distance,",
    "not a faithfulness check: whether the prose describes the body needs a",
    "model, and that layer is not built."
  ),
  "",
  rows("metric", "value"),
  rows("---", "---:"),
  rows("exports checked", fmt(length(exp_names))),
  rows("exports without a roxygen block", fmt(length(doc_missing))),
  rows("undocumented parameters", fmt(n_doc_undoc)),
  rows("documented-but-absent parameters", fmt(n_doc_phantom)),
  rows("missing `@return`", fmt(length(doc_no_return))),
  rows("missing `@examples`", fmt(length(doc_no_examples))),
  "",
  paste0(
    "Undocumented and phantom parameters are gated by `--check`. Missing ",
    "`@return`: ",
    if (length(doc_no_return)) paste0("`", paste(doc_no_return, collapse = "`, `"), "`") else "none",
    "."
  ),
  "",
  "## Semantic extraction (local LLM)",
  "",
  paste(
    "`dev/graphify-semantic.sh` is the only step that calls a model. Its output",
    "is cached under `graphify-out/cache/semantic/` and keyed by `semantic_hash`",
    "in the manifest, so an unchanged file is not re-asked. A graph built by",
    "`graphify update .` is AST-only: provenance tags every node and edge it",
    "writes `_origin: ast`, and the AST pass itself emits INFERRED `references`",
    "edges. So the INFERRED total below is **not** evidence the model ran, and a",
    "`concept` label does not prove semantic origin either. Count the pass by",
    "`_origin == \"semantic\"`; `nodes from the model` below is that count and is",
    "the only number here that the LLM pass moves. `concept` is a trap: graphify's",
    "AST pass mints `concept` nodes for Python `import` and R `source()` targets,",
    "so a nonzero `concept` count is an extractor artifact, not model output.",
    "For this corpus the LLM pass is not the first lever anyway: the R docs are",
    "structured",
    "(roxygen, `man/*.Rd`, vignette headings and code chunks), so deterministic",
    "stubs and the `@param` check above yield doc-to-code facts with no model."
  ),
  "",
  rows("metric", "value"),
  rows("---", "---:"),
  rows("document-like files", fmt(sem_total)),
  rows("with a `semantic_hash`", fmt(sem_done)),
  rows("nodes from the model", fmt(sem_nodes)),
  rows("concept nodes (any origin)", fmt(concept_nodes)),
  rows("…of which import/source targets (no `source_file`)", fmt(concept_unlocated)),
  rows("INFERRED edges", fmt(inferred_edges)),
  rows("…from the AST pass", fmt(inferred_ast)),
  rows("…from the model", fmt(inferred_semantic)),
  "",
  "## God nodes",
  "",
  "Defined here as the top targets by **in-degree over `calls`**: how often the",
  "rest of the code calls this. (Graphify's own ranking mixes document nodes and",
  "is not reproduced.)",
  "",
  paste(
    "Read this as the **third-party plumbing this package leans on**, not as the",
    "package's own interfaces. The R resolver resolves `pkg::fn()` calls to an",
    "`external_function` node but a bare internal call only inside one file (see",
    "'cross-file call resolution' below), so no internal helper can accumulate the",
    "in-degree a `tibble::tibble()` does. The real internal seams are the R-only",
    "cross-file fan-out measured from source (comments stripped, only `R/`",
    "callers, test fixtures and roxygen examples excluded): `snp_template`",
    "(16 files), `spec_row` (12), `map_spec` (13), `scheme_conditions` (10),",
    "`gc_frac` (10), `primer_candidates` (9), `primer_tm` (9), `has_palindrome`",
    "(7), `ath_chrom` (5) -- broad and shallow, with no hub, not this table.",
    "A function called only from tests, `@examples` or a demo is not a seam at",
    "all: `ath_genome()` has 0 `R/` callers and is a TAIR10 demo fixture, not an",
    "interface. This list is hand-maintained and drifts; compute it in the script",
    "from the same stripped source scan before trusting the counts."
  ),
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
    "vignettes/ and man/ coverage", "partial",
    paste(
      "`.graphifyrc` maps `.Rmd` to markdown, so `vignettes/*.Rmd` **are** indexed",
      "(91 `document` nodes). `man/*.Rd` is not, and `.graphifyrc` leaves it",
      "unmapped on purpose (174 files would dominate a semantic run); export plain",
      "`.md` stubs via `dev/graphify-docs.R` if the manuals must be graph concepts"
    )
  ),
  rows(
    "single-file echo labels (254)", "skipped",
    paste(
      "a file-echo name is honest for a one-file community; the strict pass",
      "renames only missing/duplicate/echo-stale, `--all-echo` opts in"
    )
  ),
  rows(
    "nested-local assignment leaves", "accepted",
    paste(
      "`variable` nodes include function locals: `test_that()`/`local()` bodies",
      "in `tests/` mint", fmt(var_nested_tests), "nested nodes (",
      fmt(nested_total), "nested nodes in all), and most are leaves with only",
      "their file's `contains` edge. They hold no seam; graphify has no per-kind",
      "knob, so the only lever is a `.graphifyignore` path exclusion, which also",
      "drops the `tests/` map"
    )
  ),
  rows(
    "R S3 dispatch edges", "gap",
    paste(
      "every export and 91/92 S3 methods have nodes (`+.marker_plan()` is",
      "dropped: a string-named LHS fails graphify's identifier-only `binding()`);",
      "the 7 generics without nodes are external implementations. No",
      "`generic -> method` relation and no exported/internal flag",
      "(`interface_dispatch.py` covers C#, Java, Swift and others, not R), so",
      "\"what implements `marker_chunk()`\" is a label-prefix query"
    )
  ),
  rows(
    "code ↔ documentation faithfulness", "partial",
    paste(
      "`dev/graphify-health.R` compares each export's `@param` names with",
      "`formals()` and gates undocumented or phantom parameters; whether the",
      "prose describes the body is not checked -- that needs the local model,",
      "and it is not built"
    )
  ),
  rows(
    "cross-file call resolution", "gap",
    paste(
      "The R resolver binds a bare call only inside the caller's own file, and",
      "resolves a cross-file call only through an `imports_from` edge, which it",
      "creates from an explicit `source()` call. An R package has neither: it",
      "declares dependencies in `NAMESPACE`/`Collate`, not with `source()`, so of",
      "747 `calls` edges whose ends both have a `source_file`, **0 cross a file**.",
      "`pkg::fn()` targets do resolve cross-file because they name the package.",
      "Consequence: `explain`/`affected` under-report impact, communities are",
      "file islands, and a source-scan (not the graph) is the only way to see the",
      "real seams. Fix: a deterministic post-pass keyed on a package-wide unique",
      "function name (skip ambiguous names, shadowed locals and test fixtures),",
      "emitting INFERRED `calls` edges -- the same pattern",
      "`dev/graphify-handrolled.R` uses for `supersedes`."
    )
  ),
  rows(
    "semantic layer in graph.json", "drift",
    paste(
      "`graphify update .` writes an AST-only graph (`_origin: ast` on every node",
      "and edge); the semantic cache under `graphify-out/cache/semantic/` holds",
      "chunk results that are **not** merged into `graph.json`, so a current graph",
      "has 0 `_origin: semantic` nodes and 0 `semantic_hash` in the manifest until",
      "`dev/graphify-semantic.sh` runs to completion. Re-run it when the bonsai",
      "slots are free (see the lever note on slot contention)"
    )
  ),
  rows(
    "curated label → graph.json", "adopted",
    paste(
      "`dev/relabel-bonsai.R` writes `graphify-out/.graphify_labels.json` and",
      "then propagates every name into the `community_name` cached on each node",
      "in `graph.json`, so `explain`/`query`/`path` agree with the health table.",
      "Propagation runs even when nothing was relabelled, so a rebuild without a",
      "rename still picks up names curated in an earlier pass"
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
    "local LLM semantic backend", "shipped",
    paste(
      "`extract --backend=openai` via `dev/graphify-semantic.sh`; reads",
      "`OPENAI_BASE_URL`/`OPENAI_MODEL` (bonsai), 2 slots, `--token-budget 12000`.",
      "`GRAPHIFY_MAX_OUTPUT_TOKENS` must stay at 16000: a lower cap truncates,",
      "splits and retries forever. The server has `-np 2`, so a second client on",
      "the same port (a Merovingen `acceptance_run.py` did this) holds both slots",
      "with long generations and this pass blocks as the deferred request -- that,",
      "not a config fault, is why a few chunks can take minutes. Check",
      "`llamacpp:requests_deferred` and that no other job owns the slots before",
      "blaming the backend. Coverage and yield in 'Semantic extraction'"
    )
  ),
  rows(
    "`diagnose multigraph`", "shipped",
    paste(
      "`graphify diagnose multigraph --json` independently confirms 0 dangling,",
      "0 duplicate and 0 collapsed edges, and lists the missing-`source_file`",
      "warnings (the external nodes in the table above); it does not replace",
      "`--check`, which gates the same invariants plus cycles and labels"
    )
  ),
  rows(
    "`export callflow-html`", "open",
    paste(
      "Mermaid architecture/call-flow HTML from the same graph; a candidate for",
      "the architecture doc, but no reader is wired to it yet"
    )
  ),
  rows(
    "`GRAPHIFY_QUERY_LOG`", "open",
    paste(
      "opt-in log of agent `query`/`explain` calls; would show which nodes agents",
      "actually ask for, so the graph can be tuned to real questions. Off by",
      "default; enable when there is a reason to read the log"
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
    "per-kind node pruning", "skipped",
    paste(
      "no CLI or rc knob drops `variable` leaves; only a `.graphifyignore` path",
      "exclusion can, at the cost of the `tests/`/`data-raw/`/`dev/` map"
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

report <- c(
  report, "", "## Rerun", "", "```sh", "graphify update .",
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
  "`dev/relabel-bonsai.R` if it warns."
)

if (identical(mode, "write")) {
  dir.create(dirname(REPORT), showWarnings = FALSE, recursive = TRUE)
  writeLines(report, REPORT, useBytes = TRUE)
  message("wrote ", REPORT, " (", length(nodes), " nodes, ", length(links), " edges)")
  quit(status = 0L)
}

# --- --check: regressions fail; drift only warns ------------------------------
regressions <- c(
  if (dangling > 0L) paste0(dangling, " dangling endpoint(s)"),
  if (doc_ok && length(doc_missing)) {
    paste0(length(doc_missing), " export(s) without a roxygen block")
  },
  if (doc_ok && n_doc_undoc > 0L) paste0(n_doc_undoc, " undocumented parameter(s)"),
  if (doc_ok && n_doc_phantom > 0L) {
    paste0(n_doc_phantom, " documented-but-absent parameter(s)")
  },
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
