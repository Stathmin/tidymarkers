# =============================================================================
#  marker_decision_layout.R  --  where a decision graph's nodes sit
# =============================================================================
#  The decision schema (marker_decision.R) says *what* is connected to what; it
#  says nothing about where to draw it.  That is a layout problem, and a layout
#  problem with boxes, stage containers and 8-14 edges into one sink needs a
#  solver, not a hand-tuned grid.  This file is the seam: one interface,
#
#      decision_layout(g) -> tibble(name, x, y) + a stage/bbox attribute
#
#  with two adapters behind it.  `"sugiyama"` runs the layered algorithm with
#  every node's rank *pinned to its stage*, so a stage is one column; `"manual"`
#  is the pure-R grid that needs nothing and is what the package falls back to
#  (and what the tests pin).  Callers (plot_decision_tree) never see
#  which one ran except through the returned `engine` attribute.
# =============================================================================

# node geometry, in inches at 72 dpi: one place for the assumption the plot and
# the DOT writer share.  `stringr::str_wrap(width = 16)` is the line width, so a
# label is at most ~16 characters of the 9 pt face wide.
decision_geometry <- function() {
  list(
    char_width = 0.062, # inches per character at 9 pt Helvetica
    line_height = 0.155, # inches per wrapped line
    pad_x = 0.14, pad_y = 0.09, # box padding, inches (both sides)
    wrap = 16L, fontsize = 9,
    ranksep = 0.55, nodesep = 0.18
  )
}

# the graph's stage order, as the layout adapters both want it
decision_stage_order <- function(nodes) unique(as.character(nodes$stage))

# ---- shared: the boxes the solver must respect ------------------------------

#' Node boxes for an external layout solver
#'
#' DOT has no idea of the text a box will hold, so every node that carries a
#' `label` is handed a fixed size derived from the *wrapped* label.  Without
#' this, an external solver sizes the box from its own idea of a label and a
#' 9 pt box comes out a few millimetres wide while the actual label is a
#' centimetre.
#' @keywords internal
#' @noRd
decision_dot_boxes <- function(nodes, geom = decision_geometry()) {
  lines <- vapply(strsplit(nodes$label_wrapped, "\n", fixed = TRUE), length, integer(1))
  chars <- nchar(gsub("\n", " ", nodes$label_wrapped, fixed = TRUE))
  tibble::tibble(
    name = nodes$name,
    width = chars * geom$char_width + 2 * geom$pad_x,
    height = lines * geom$line_height + 2 * geom$pad_y
  )
}

# is the rank-pinned layout usable?
# `TIDYKASP_NO_DOT` is honoured as an alias: it was the switch before the
# adapter could box and route on its own, and a pinned test still sets it
#' @keywords internal
#' @noRd
decision_sugiyama_available <- function() {
  if (nzchar(Sys.getenv("TIDYKASP_NO_SUGIYAMA")) || nzchar(Sys.getenv("TIDYKASP_NO_DOT"))) {
    return(FALSE)
  }
  requireNamespace("igraph", quietly = TRUE)
}

# ---- adapter 1: the manual grid (no external dependency) --------------------

#' The hand-placed fallback layout
#'
#' One column per stage, stages wrapped `ncol` per row so the picture stays near
#' 4:3 instead of a 1:2.3 strip, nodes stacked inside the column.  Deterministic
#' and dependency-free: the tests pin it with `TIDYKASP_NO_SUGIYAMA=1`.
#' @keywords internal
#' @noRd
decision_layout_manual <- function(nodes) {
  order <- decision_stage_order(nodes)
  nodes$stage <- factor(as.character(nodes$stage), levels = order)
  ncol <- max(2L, ceiling(sqrt(length(order) * 1.6)))
  col_of <- match(as.character(nodes$stage), order) - 1L
  row_of <- col_of %/% ncol
  nodes <- nodes |>
    dplyr::group_by(.data$stage) |>
    dplyr::mutate(lane = dplyr::row_number()) |>
    dplyr::ungroup()
  nodes$x <- (as.integer(col_of) %% ncol) * 3.6 + as.integer(row_of) * 0.55
  nodes$y <- -as.integer(row_of) * (max(nodes$lane) + 2.2) - nodes$lane
  nodes$x[nodes$kind == "relax"] <- nodes$x[nodes$kind == "relax"] + 0.45
  out <- as.data.frame(nodes[order(nodes$order), c("name", "x", "y")])
  attr(out, "routes") <- list()
  out
}

# ---- adapter 2: pinned Sugiyama (stages are ranks) --------------------------

#' The rank-pinned Sugiyama adapter
#'
#' A decision graph is not a free graph: its stages *are* ranks.  An unpinned
#' solver spends a layer per edge, so every gate chain becomes its own column -
#' duplex came out 1:4.5 tall because 23 assay checks stacked as 23 ranks.
#' Pinning every node to its stage's layer makes each stage one column, which
#' is what the picture claims, and leaves the solver to order the column.
#' (`ggraph`'s Sugiyama honours a graph-wide `rank = same` for *every* node,
#' which puts all stages on a single row; a per-node pin is how to say "these
#' share a column, those share the next".)
#'
#' The solver answers *which box comes first in a column*; the column order is
#' the schema's, and the metric is ours, because only this file knows how big a
#' box will be.  Columns are spaced by the widest box, and every row sits on
#' one shared pitch, so a route between two stages leaves a row and lands on a
#' row.
#' @keywords internal
#' @noRd
decision_layout_sugiyama <- function(nodes, edges) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is needed for the rank-pinned layout")
  }
  geom <- decision_geometry()
  boxes <- decision_dot_boxes(nodes, geom)
  w <- boxes$width[match(nodes$name, boxes$name)]
  # drop duplicate (from, to) pairs: the schema repeats a couple of `doomed`
  # edges, and a multigraph buys the solver nothing
  key <- paste(edges$from, edges$to, sep = "\r")
  edges <- edges[!duplicated(key), , drop = FALSE]
  ig <- igraph::graph_from_data_frame(
    as.data.frame(edges[, c("from", "to")]),
    directed = TRUE,
    vertices = data.frame(name = nodes$name)
  )
  stages <- decision_stage_order(nodes)
  # igraph numbers layers from the bottom, so a stage's layer index is reversed
  layers <- length(stages) - match(as.character(nodes$stage), stages)
  ly <- igraph::layout_with_sugiyama(ig, layers = layers, hgap = 1, vgap = 1)
  xy <- ly$layout[match(nodes$name, igraph::V(ig)$name), , drop = FALSE]

  # one slot height for the whole graph: the tallest box plus the gap the
  # geometry already reserves, so a column's rows never touch and two stages'
  # rows line up for comparison
  slot <- max(boxes$height) + 2 * geom$pad_y + 0.16

  x <- numeric(nrow(nodes))
  y <- numeric(nrow(nodes))
  left <- 0
  for (k in seq_along(stages)) {
    idx <- which(as.character(nodes$stage) == stages[k])
    if (!length(idx)) next
    # the solver's order first; a column it left unlayered (a stage whose nodes
    # have no edge among them, e.g. the primer scheme's lone site node) falls
    # back to the schema's own order, never to the row order of a data frame
    idx <- idx[order(xy[idx, 2L], nodes$order[idx])]
    # all boxes of a column share one centre: a column reads as a column
    x[idx] <- left + max(w[idx]) / 2
    # every column shares one row pitch and starts at the same height, so an
    # arrow between two stages is a straight elbow at a row, not an S-bend from
    # a box to a gap
    y[idx] <- rev(seq_along(idx) - 1L) * -slot
    left <- left + max(w[idx]) + 0.5
  }
  data.frame(name = nodes$name, x = x, y = y, stringsAsFactors = FALSE)
}

# a routed edge, in the order the plot draws it
decision_routes <- function(nodes, edges) {
  nodes$stage <- as.character(nodes$stage)
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  route <- function(e) {
    a <- nodes[nodes$name == e$from, , drop = FALSE]
    b <- nodes[nodes$name == e$to, , drop = FALSE]

    if (nrow(a) == 0L || nrow(b) == 0L) {
      return(NULL)
    }
    if (identical(a$stage, b$stage)) {
      # inside a column: leave the bottom of one box, land on the top of the
      # next, step sideways just enough to be seen
      from_bottom <- a$y - a$h / 2
      mid <- (from_bottom + b$y + b$h / 2) / 2
      side <- a$x + a$w / 2 + 0.09
      list(
        from = e$from, to = e$to,
        x = c(a$x, a$x, side, side, b$x),
        y = c(from_bottom, mid, mid, b$y + b$h / 2, b$y + b$h / 2)
      )
    } else {
      # across columns: leave the right side, land on the left side, and keep
      # the arrow in the gap between the two stage bands
      gx <- (a$x + a$w / 2 + b$x - b$w / 2) / 2
      list(from = e$from, to = e$to, x = c(a$x + a$w / 2, gx, gx, b$x - b$w / 2), y = c(a$y, a$y, b$y, b$y))
    }
  }
  out <- lapply(seq_len(nrow(edges)), function(i) route(edges[i, , drop = FALSE]))
  Filter(Negate(is.null), out)
}

# ---- the seam ---------------------------------------------------------------

#' Where the decision graph's nodes sit
#'
#' The single interface both the picture and the tests use: node names and
#' edges in, coordinates in inches (y up, as `ggplot2` wants) out.  The
#' `"sugiyama"` adapter is used when `igraph` is installed and `TIDYKASP_NO_SUGIYAMA`
#' is unset; otherwise (or if the solver fails) the manual grid is used and a
#' message says so - a missing solver must never cost a render.  The two adapters do
#' **not** agree on coordinates (one is solved, one is a grid); they agree on
#' the node set, which is what a test may rely on.
#' @param g a `marker_decision_graph()` tidygraph, or its node/edge tables.
#' @param engine `"auto"` (default), `"sugiyama"` or `"manual"`.
#' @return a data frame of `name`, `x`, `y` (inches) carrying attributes
#'   `stage` (the scheme's stage order), `nodes` (the laid-out node table with
#'   per-node `w`/`h` in inches), `edges`, `routes` (one routed polyline per
#'   edge the solver routed, `from`/`to`/`x`/`y`) and `engine` (`"sugiyama"` /
#'   `"manual"`).
#' @keywords internal
#' @noRd
decision_layout <- function(g, engine = c("auto", "sugiyama", "manual")) {
  engine <- match.arg(engine)
  nodes <- if (inherits(g, "tbl_graph")) tidygraph::as_tibble(g) else g$nodes
  edges <- if (inherits(g, "tbl_graph")) decision_edges(g) else g$edges
  nodes <- nodes |>
    dplyr::arrange(.data$order) |>
    dplyr::mutate(label_wrapped = stringr::str_wrap(.data$label, decision_geometry()$wrap))
  stages <- decision_stage_order(nodes)
  nodes$stage <- factor(as.character(nodes$stage), levels = stages)
  edges$from_stage <- as.character(nodes$stage[match(edges$from, nodes$name)])
  edges$to_stage <- as.character(nodes$stage[match(edges$to, nodes$name)])

  solve <- switch(engine,
    sugiyama = TRUE,
    manual = FALSE,
    auto = decision_sugiyama_available()
  )
  pos <- NULL
  if (solve) {
    pos <- tryCatch(decision_layout_sugiyama(nodes, edges), error = function(e) {
      message("layered layout failed (", conditionMessage(e), "); using the manual grid")
      NULL
    })
  }
  if (is.null(pos)) {
    pos <- decision_layout_manual(nodes)
    got <- "manual"
  } else {
    got <- "sugiyama"
  }
  boxes <- decision_dot_boxes(nodes)
  nodes$w <- boxes$width[match(nodes$name, boxes$name)]
  nodes$h <- boxes$height[match(nodes$name, boxes$name)]
  nodes$x <- pos$x[match(nodes$name, pos$name)]
  nodes$y <- pos$y[match(nodes$name, pos$name)]
  routes <- decision_routes(nodes, edges)
  out <- as.data.frame(pos)
  attr(out, "stage") <- stages
  attr(out, "nodes") <- as.data.frame(nodes)
  attr(out, "edges") <- as.data.frame(edges)
  attr(out, "routes") <- routes
  attr(out, "engine") <- got
  out
}

#' The bounding box of a laid-out graph
#' @keywords internal
#' @noRd
decision_bbox <- function(layout) {
  n <- attr(layout, "nodes")
  pad <- 0.55 # room for the stage labels above each cluster
  list(
    xmin = min(n$x - n$w / 2) - pad, xmax = max(n$x + n$w / 2) + pad,
    ymin = min(n$y - n$h / 2) - 0.25, ymax = max(n$y + n$h / 2) + pad
  )
}
