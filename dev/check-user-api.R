# =============================================================================
#  dev/check-user-api.R  --  the user API surface gate
# =============================================================================
#  Checks that docs/user-api-surface.csv and docs/user-api-graph.dot still
#  describe the exported surface:
#
#    1. every NAMESPACE export has a kind = "export" row, and no row is stale
#    2. every S3method(generic, class) has a kind = "s3method" row
#    3. every row has a known role, a boolean user_facing and a class
#    4. every export is named in docs/user-api-surface.md
#    5. every confirmed edge in the MD exists, fits by class, and is in the graph
#    6. docs/user-api-graph.dot is exactly what this script generates
#
#  The graph is derived, never hand-edited.  `--write` regenerates the DOT file:
#
#      Rscript dev/check-user-api.R --write
#      Rscript dev/check-user-api.R
#
#  Exit 1 on any failure.  Base R only: the gate runs before any package load.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

failures <- character()
source(file.path(root, "dev", "check-reporter.R"))

csv_path <- file.path(root, "docs", "user-api-surface.csv")
md_path <- file.path(root, "docs", "user-api-surface.md")
dot_path <- file.path(root, "docs", "user-api-graph.dot")

ROLES <- c(
  "ingest", "design", "specificity", "assembly", "criteria",
  "decision", "artefact", "export", "service"
)
UNTYPED <- c("none", "any", "character", "numeric", "integer", "logical", "tibble", "path")
TERMINAL <- c(
  "files", "artefact", "decision", "plot", "svg", "report", "recipe",
  "explain", "readout", "risk", "replay", "pareto", "informativeness",
  "verdict", "tube_pairs", "registry", "capacity"
)

read_exports <- function(path) {
  lines <- readLines(path, warn = FALSE)
  ex <- grep("^export\\(", lines, value = TRUE)
  ex <- sub("^export\\((.*)\\)$", "\\1", ex)
  unique(unlist(strsplit(ex, "[[:space:]]+")))
}

read_s3 <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sm <- grep("^S3method\\(", lines, value = TRUE)
  inner <- sub("^S3method\\((.*)\\)$", "\\1", sm)
  parts <- strsplit(inner, ",", fixed = TRUE)
  vapply(parts, function(p) {
    generic <- gsub('"', "", trimws(p[[1L]]))
    cls <- gsub('"', "", trimws(p[[2L]]))
    paste0(generic, ".", cls)
  }, character(1))
}

if (!file.exists(csv_path)) {
  note("missing ", csv_path)
  quit(status = 1L)
}
if (!file.exists(md_path)) note("missing ", md_path)

api <- read.csv(csv_path, colClasses = "character", check.names = FALSE)
required_cols <- c(
  "function", "kind", "role", "user_facing", "required_args",
  "optional_args", "inputs", "output", "output_class", "entry_point",
  "leaf", "notes"
)
missing_cols <- setdiff(required_cols, names(api))
if (length(missing_cols)) {
  note("CSV is missing columns: ", paste(missing_cols, collapse = ", "))
  quit(status = 1L)
}

exports <- read_exports(file.path(root, "NAMESPACE"))
s3 <- read_s3(file.path(root, "NAMESPACE"))

# 1. coverage both ways -------------------------------------------------------
csv_export <- api[['function']][api$kind == "export"]
csv_s3 <- api[['function']][api$kind == "s3method"]

miss_exp <- setdiff(exports, csv_export)
stale_exp <- setdiff(csv_export, exports)
miss_s3 <- setdiff(s3, csv_s3)
stale_s3 <- setdiff(csv_s3, s3)

if (length(miss_exp)) note("export without a CSV row: ", paste(miss_exp, collapse = ", "))
if (length(stale_exp)) note("CSV export row is stale: ", paste(stale_exp, collapse = ", "))
if (length(miss_s3)) note("S3 method without a CSV row: ", paste(miss_s3, collapse = ", "))
if (length(stale_s3)) note("CSV S3 row is stale: ", paste(stale_s3, collapse = ", "))
if (anyDuplicated(api[['function']])) {
  note("duplicate function rows: ", paste(unique(api[['function']][duplicated(api[['function']])]), collapse = ", "))
}
if (!length(miss_exp) && !length(stale_exp) && !length(miss_s3) && !length(stale_s3)) {
  ok("CSV covers ", length(exports), " exports and ", length(s3), " S3 methods")
}

# 2. field contracts ---------------------------------------------------------
bad_role <- setdiff(unique(api$role), ROLES)
if (length(bad_role)) note("unknown role: ", paste(bad_role, collapse = ", "))
bad_kind <- setdiff(unique(api$kind), c("export", "s3method"))
if (length(bad_kind)) note("unknown kind: ", paste(bad_kind, collapse = ", "))
bad_bool <- unique(api$user_facing[!api$user_facing %in% c("TRUE", "FALSE")])
if (length(bad_bool)) note("user_facing is not TRUE/FALSE: ", paste(bad_bool, collapse = ", "))
bad_entry <- unique(api$entry_point[!api$entry_point %in% c("TRUE", "FALSE")])
if (length(bad_entry)) note("entry_point is not TRUE/FALSE: ", paste(bad_entry, collapse = ", "))
bad_leaf <- unique(api$leaf[!api$leaf %in% c("TRUE", "FALSE")])
if (length(bad_leaf)) note("leaf is not TRUE/FALSE: ", paste(bad_leaf, collapse = ", "))
empty_class <- api[['function']][!nzchar(api$output_class)]
if (length(empty_class)) note("empty output_class: ", paste(empty_class, collapse = ", "))
empty_inputs <- api[['function']][!nzchar(api$inputs)]
if (length(empty_inputs)) note("empty inputs (use `none`): ", paste(empty_inputs, collapse = ", "))
if (!length(bad_role) && !length(bad_kind) && !length(bad_bool) &&
  !length(bad_entry) && !length(bad_leaf) && !length(empty_class) &&
  !length(empty_inputs)) {
  ok("all rows have role, kind, booleans, inputs and an output class")
}

# 3. the MD names every export ----------------------------------------------
md <- if (file.exists(md_path)) paste(readLines(md_path, warn = FALSE), collapse = "\n") else ""
unmentioned <- exports[!vapply(
  exports, function(fn) grepl(paste0("`", fn, "()`"), md, fixed = TRUE), logical(1)
)]
if (length(unmentioned)) {
  note("export not named in the MD inventory: ", paste(unmentioned, collapse = ", "))
} else {
  ok("MD names all ", length(exports), " exports")
}

# 4. class fit and graph -----------------------------------------------------
split_inputs <- function(x) {
  if (!nzchar(x) || x == "none") character() else strsplit(x, ";", fixed = TRUE)[[1L]]
}
fits <- function(out, input) {
  any(out %in% input) || ("assay" %in% out && "marker" %in% input)
}

api_exp <- api[api$kind == "export", , drop = FALSE]
special <- data.frame(
  a = "+.marker_plan", b = "decision", c = "plan;targets;context;tech",
  d = "plan", e = "TRUE", f = "FALSE", g = "FALSE",
  stringsAsFactors = FALSE
)
names(special) <- c("function", "role", "inputs", "output_class", "user_facing", "entry_point", "leaf")
nodes <- rbind(
  api_exp[, c("function", "role", "inputs", "output_class", "user_facing", "entry_point", "leaf")],
  special
)

edge <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(i) {
  out <- nodes$output_class[[i]]
  if (out %in% UNTYPED) {
    return(NULL)
  }
  to <- vapply(seq_len(nrow(nodes)), function(j) {
    if (i == j) {
      return(FALSE)
    }
    any(fits(out, split_inputs(nodes$inputs[[j]])))
  }, logical(1))
  if (!any(to)) {
    return(NULL)
  }
  data.frame(from = nodes[['function']][[i]], to = nodes[['function']][to], class = out)
}))
if (is.null(edge)) edge <- data.frame(from = character(), to = character(), class = character())
edge <- unique(edge)
indeg <- table(factor(edge$to, levels = nodes[['function']]))
outdeg <- table(factor(edge$from, levels = nodes[['function']]))

# confirmed edges ------------------------------------------------------------
md_lines <- strsplit(md, "\n", fixed = TRUE)[[1L]]
start <- match("<!-- confirmed-edges:start -->", md_lines)
end <- match("<!-- confirmed-edges:end -->", md_lines)
if (is.na(start) || is.na(end) || end <= start) {
  note("MD has no confirmed-edges block")
  conf <- data.frame(from = character(), to = character(), source = character())
} else {
  block <- md_lines[(start + 1L):(end - 1L)]
  block <- block[nzchar(trimws(block))]
  m <- regexec("^- `([^`]+)` -> `([^`]+)` — (.*)$", block)
  hits <- regmatches(block, m)
  ok_lines <- vapply(hits, length, integer(1)) == 4L
  if (!all(ok_lines)) note("malformed confirmed edge line(s): ", paste(block[!ok_lines], collapse = " | "))
  conf <- do.call(rbind, lapply(hits[ok_lines], function(h) {
    data.frame(from = h[[2L]], to = h[[3L]], source = h[[4L]], stringsAsFactors = FALSE)
  }))
  if (is.null(conf)) conf <- data.frame(from = character(), to = character(), source = character())
}

if (nrow(conf)) {
  known <- c(nodes[['function']], "report", "panel")
  bad_end <- conf[!conf$from %in% known | !conf$to %in% known, , drop = FALSE]
  if (nrow(bad_end)) {
    note("confirmed edge names an unknown function: ",
      paste(unique(c(bad_end$from, bad_end$to)), collapse = ", "))
  }
  fit_ok <- vapply(seq_len(nrow(conf)), function(i) {
    a <- conf$from[[i]]
    b <- conf$to[[i]]
    if (a %in% c("report", "panel") || b %in% c("report", "panel")) {
      return(TRUE)
    }
    out <- nodes$output_class[match(a, nodes[['function']])]
    any(fits(out, split_inputs(nodes$inputs[match(b, nodes[['function']])])))
  }, logical(1))
  if (!all(fit_ok)) {
    note("confirmed edge does not fit by class: ",
      paste(paste0(conf$from[!fit_ok], " -> ", conf$to[!fit_ok]), collapse = ", "))
  }
  in_graph <- paste(conf$from, conf$to) %in% paste(edge$from, edge$to)
  if (!all(in_graph)) {
    note("confirmed edge missing from the generated graph: ",
      paste(paste0(conf$from[!in_graph], " -> ", conf$to[!in_graph]), collapse = ", "))
  }
  if (all(fit_ok) && all(in_graph)) ok("all ", nrow(conf), " confirmed edges fit and are in the graph")
}

# 5. DOT ---------------------------------------------------------------------
palette <- c(
  ingest = "#cfe8ff", design = "#d9f2d9", specificity = "#ffe6cc",
  assembly = "#e6d9f2", criteria = "#fff2b3", decision = "#f2d9e6",
  artefact = "#d9ecec", export = "#e8e8e8", service = "#f5f5f5"
)
shape_of <- function(i) {
  if (nodes$entry_point[[i]] == "TRUE") {
    "invhouse"
  } else if (nodes$leaf[[i]] == "TRUE") {
    "box"
  } else {
    "ellipse"
  }
}
dot_escape <- function(x) gsub('"', '\\\\"', x, fixed = TRUE)

dot <- c(
  "// Generated by dev/check-user-api.R --write. Do not edit by hand.",
  "digraph user_api {",
  '  graph [rankdir=LR, overlap=false, splines=true, fontname="Helvetica", bgcolor="white"];',
  '  node [fontname="Helvetica", fontsize=9, style="filled", color="#555555"];',
  '  edge [color="#b0b0b0", arrowsize=0.6];'
)
for (role in ROLES) {
  idx <- which(nodes$role == role)
  if (!length(idx)) next
  dot <- c(dot, sprintf('  subgraph cluster_%s {', role))
  dot <- c(dot, sprintf('    label="%s"; style="rounded,dashed"; color="#cccccc";', role))
  for (i in idx) {
    style <- if (nodes$user_facing[[i]] == "TRUE") "filled" else "filled,dashed"
    dot <- c(dot, sprintf(
      '    "%s" [label="%s", shape=%s, style="%s", fillcolor="%s"];',
      dot_escape(nodes[['function']][[i]]), dot_escape(nodes[['function']][[i]]),
      shape_of(i), style, palette[[role]]
    ))
  }
  dot <- c(dot, "  }")
}
conf_key <- paste(conf$from, conf$to)
for (i in seq_len(nrow(edge))) {
  is_conf <- paste(edge$from[[i]], edge$to[[i]]) %in% conf_key
  attrs <- if (is_conf) {
    ' [color="#222222", penwidth=2.2]'
  } else {
    sprintf(' [color="#c8c8c8", penwidth=0.6]')
  }
  dot <- c(dot, sprintf(
    '  "%s" -> "%s"%s;', dot_escape(edge$from[[i]]),
    dot_escape(edge$to[[i]]), attrs
  ))
}
dot <- c(dot, "}")

write_dot <- "--write" %in% commandArgs(trailingOnly = TRUE)
if (write_dot) {
  writeLines(dot, dot_path)
  ok("wrote ", basename(dot_path), " (", nrow(nodes), " nodes, ", nrow(edge), " edges)")
} else if (!file.exists(dot_path)) {
  note("missing ", dot_path, "; run with --write")
} else {
  committed <- readLines(dot_path, warn = FALSE)
  if (!identical(committed, dot)) {
    note("DOT is stale; run `Rscript dev/check-user-api.R --write`")
  } else {
    ok("DOT matches the generated graph (", nrow(nodes), " nodes, ", nrow(edge), " edges)")
  }
}

# 6. summary -----------------------------------------------------------------
counts <- table(factor(nodes$role, levels = ROLES))
message(
  "  api  ", nrow(nodes), " nodes, ", nrow(edge), " typed edges; roles: ",
  paste(names(counts), as.integer(counts), sep = "=", collapse = " ")
)
message(
  "  api  ", sum(nodes$entry_point == "TRUE"), " entries, ",
  sum(nodes$leaf == "TRUE"), " leaves, ",
  sum(indeg == 0 & outdeg == 0), " isolates"
)

if (length(failures)) {
  message("\nFAIL user-api: ", length(failures), " problem(s)")
  for (f in failures) message("  - ", f)
  quit(status = 1L)
}
message("ok user-api")
