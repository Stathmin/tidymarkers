# =============================================================================
#  marker_decision_solver.R  --  the seam to the external layout/routing engine
# =============================================================================
#  A decision DAG is *placed* by ELK (layered) and *routed* by libavoid
#  (obstacle-avoiding orthogonal), both vendored under `.tools/elk/`.  This file
#  is the only place that knows the engine exists:
#
#      structure  (marker_decision_graph)   what is connected to what
#      solver     (this file)               where it sits, how an edge is drawn
#      renderer   (marker_decision_render)  ink: colour, text, page
#
#  The exchange is one JSON document each way over stdin/stdout, so no R code
#  ever computes a port, a lane, a gap or a bend.  Coordinates come back in
#  points with y pointing DOWN (SVG's convention, not ggplot's).
# =============================================================================

#' Directory holding the vendored layout engine
#'
#' The engine (`elkjs` + `libavoid` + `opentype.js`) ships inside the package
#' sources under `.tools/elk/node_modules`, so a render never downloads anything.
#' When the sources are not the working copy (an installed package), or
#' `TIDYKASP_ELK_DIR` points elsewhere, that path wins.
#' @return a path, or `NA_character_` when no vendored engine is visible.
#' @keywords internal
#' @noRd
decision_solver_dir <- function() {
  cand <- c(
    Sys.getenv("TIDYKASP_ELK_DIR", ""),
    file.path(getwd(), ".tools", "elk"),
    file.path("..", "..", ".tools", "elk")
  )
  hit <- cand[nzchar(cand) & file.exists(file.path(cand, "render-layout.js"))]
  if (!length(hit)) NA_character_ else normalizePath(hit[[1L]], mustWork = TRUE)
}

#' The node executable the solver runs under, or `NA` when there is none
#' @keywords internal
#' @noRd
decision_solver_node <- function() {
  node <- Sys.which("node")
  if (nzchar(node)) node else NA_character_
}

#' Is the layout engine usable right now?
#' @return `TRUE` when the script, its modules and a `node` are all present.
#' @keywords internal
#' @noRd
decision_solver_available <- function() {
  dir <- decision_solver_dir()
  !nzchar(Sys.getenv("TIDYKASP_NO_ELK")) &&
    !is.na(dir) &&
    dir.exists(file.path(dir, "node_modules")) &&
    !is.na(decision_solver_node())
}

#' The font file the solver measures text with
#'
#' Text is measured from the *real* face: the solver parses the TTF with
#' opentype.js, so a box is exactly as wide as the label it holds.  A missing
#' file falls back to Liberation Sans, then DejaVu, then Helvetica.
#' @param bold measure with the bold face (hard-gate boxes are drawn bold).
#' @return a path to a TTF.
#' @keywords internal
#' @noRd
decision_font_file <- function(bold = FALSE) {
  face <- if (bold) "LiberationSans-Bold.ttf" else "LiberationSans-Regular.ttf"
  fallback <- if (bold) "DejaVuSans-Bold.ttf" else "DejaVuSans.ttf"
  cand <- c(
    Sys.getenv("TIDYKASP_ELK_FONT", ""),
    file.path("/usr/share/fonts/truetype/liberation", face),
    file.path("/usr/share/fonts/truetype/dejavu", fallback),
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
  )
  hit <- cand[nzchar(cand) & file.exists(cand)]
  if (!length(hit)) stop("no TrueType font found for the decision layout", call. = FALSE)
  hit[[1L]]
}

#' Solve one decision graph: boxes measured, placed and edges routed
#'
#' @param graph the structure (`marker_decision_graph()`, or any graph whose
#'   nodes carry `name`, `label`, `stage`, `kind`, `hard`).
#' @param style a [decision_style()].
#' @return a list with `width`, `height`, `nodes` (id, x, y, w, h, lines),
#'   `edges` (from, to, kind, points) and `groups` (id, x, y, w, h), in points
#'   with y down; or `NULL` when no engine is available (the caller then falls
#'   back and says so).
#' @keywords internal
#' @noRd
decision_solve <- function(graph, style = decision_style()) {
  if (!decision_solver_available()) {
    return(NULL)
  }
  n <- decision_graph_nodes(graph)
  e <- decision_edges(graph)
  req <- list(
    nodes = lapply(seq_len(nrow(n)), function(i) {
      list(
        id = n$name[i], label = n$label[i],
        group = n$group[i], kind = n$kind[i],
        status = n$status[i], bold = isTRUE(n$hard[i])
      )
    }),
    edges = lapply(seq_len(nrow(e)), function(i) {
      list(from = e$from[i], to = e$to[i], kind = e$kind[i])
    }),
    style = list(
      font_file = decision_font_file(FALSE),
      bold_file = decision_font_file(TRUE),
      font_size = style$font_size,
      line_height = style$line_height,
      max_width = style$wrap,
      pad_x = style$pad_x,
      pad_y = style$pad_y,
      margin = style$margin,
      segment_penalty = style$segment_penalty,
      crossing_penalty = style$crossing_penalty,
      shape_buffer = style$shape_buffer,
      nudge = style$nudge
    )
  )
  dir <- decision_solver_dir()
  # `stdin` in processx is a *file path* (a JSON string is parsed as a filename
  # and overruns ENAMETOOLONG), so the request goes through a temp file and the
  # reply comes back on a pipe
  req_file <- tempfile(fileext = ".json")
  on.exit(unlink(req_file), add = TRUE)
  writeLines(jsonlite::toJSON(req, auto_unbox = TRUE, null = "null", digits = NA), req_file)
  res <- tryCatch(
    processx::run(
      decision_solver_node(),
      args = c(file.path(dir, "render-layout.js"), req_file),
      stdout = "|", stderr = "|", error_on_status = FALSE, timeout = 60000
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    message("decision layout engine failed to run: ", conditionMessage(res))
    return(NULL)
  }
  if (res$status != 0L) {
    message("decision layout engine failed (status ", res$status, "): ", trimws(res$stderr))
    return(NULL)
  }
  out <- tryCatch(jsonlite::fromJSON(res$stdout, simplifyDataFrame = TRUE), error = function(e) NULL)
  if (is.null(out) || is.null(out$nodes) || !nrow(out$nodes)) {
    message("decision layout engine returned no placement")
    return(NULL)
  }
  out
}

#' The structure a solver consumes: ids, labels, groups, kinds - no positions
#'
#' The schema carries `x`/`y`/`order` because the text export walks it in
#' pipeline order; none of that reaches the solver.  `stage` becomes the group
#' (the band a node is drawn inside) and `hard` becomes `bold`.
#' @param graph a `tbl_graph` from [marker_decision_graph()].
#' @return a tibble: `name`, `label`, `group`, `kind`, `hard`, `status`.
#' @keywords internal
#' @noRd
decision_graph_nodes <- function(graph) {
  n <- tidygraph::as_tibble(graph)
  tibble::tibble(
    name = n$name,
    label = n$label,
    group = as.character(n$stage),
    kind = as.character(n$kind),
    hard = if ("hard" %in% names(n)) true_or_false(n$hard) else FALSE,
    status = if ("status" %in% names(n)) as.character(n$status) else "neutral"
  )
}

# logical vector, never NA (an NA would break `isTRUE` in the solver request)
true_or_false <- function(x) {
  !is.na(x) & x
}
