# =============================================================================
#  marker_decision_render.R  --  structure + style in, one SVG out
# =============================================================================
#  This is the whole plotting seam.  It takes the *structure* of a decision
#  graph (nodes, groups, kinds) and a *style* (palette, font size, page), asks
#  the solver where everything goes, and writes the markup.  There is no
#  coordinate arithmetic in this file, by design: a position only ever comes
#  from `decision_solve()`, so the picture cannot disagree with its solver.
# =============================================================================

#' The look of a decision figure
#'
#' Style is ink and page, nothing else: colours, the text face and size, how a
#' label wraps, how much air a box gets, and the page shape.  Placement is not
#' here - the solver owns it - except for the spacing *inputs* it is tuned with.
#' @param palette a list with `kind`, `status`, `edge`, `edge_lty` and `band`
#'   (see [decision_palette()]).
#' @param font_size label size in points.
#' @param line_height line box, as a multiple of the font size.
#' @param wrap widest a label may get before it wraps (points).
#' @param pad_x,pad_y air between text and box edge (points).
#' @param margin air between the drawing and the page edge (points).
#' @param segment_penalty,target_penalty router preferences: cost of an extra
#'   bend and of a turn away from the target side.  Higher gives straighter,
#'   fewer-bend routes.
#' @param crossing_penalty cost of one edge crossing another.
#' @param shape_buffer,nudge distance kept from a box, and between parallel
#'   segments.
#' @param title,subtitle,caption the figure's text (empty drops the part).
#' @param legend draw the kind/status/line legend.
#' @return a list of class `decision_style`.
#' @family KASP checks
#' @export
decision_style <- function(palette = decision_palette(),
                           font_size = 9,
                           line_height = 1.22,
                           wrap = 110,
                           pad_x = 5,
                           pad_y = 3,
                           margin = 14,
                           segment_penalty = 24,
                           target_penalty = 30,
                           crossing_penalty = 40,
                           shape_buffer = 3,
                           nudge = 6,
                           title = NULL,
                           subtitle = NULL,
                           caption = NULL,
                           legend = TRUE) {
  structure(
    list(
      palette = palette,
      font_size = font_size, line_height = line_height, wrap = wrap,
      pad_x = pad_x, pad_y = pad_y, margin = margin,
      segment_penalty = segment_penalty, target_penalty = target_penalty,
      crossing_penalty = crossing_penalty, shape_buffer = shape_buffer,
      nudge = nudge,
      title = title, subtitle = subtitle, caption = caption, legend = legend
    ),
    class = "decision_style"
  )
}

#' Palette of a decision figure: fill by kind, border by status, line by edge
#' @return a list with `kind`, `status`, `edge`, `edge_lty`, `band`.
#' @family KASP checks
#' @export
decision_palette <- function() {
  list(
    kind = c(
      start = "#F5F5F4", check = "#F3F4F6", gate = "#FDE68A",
      decision = "#DDD6FE", relax = "#FED7AA", merge = "#D1FAE5",
      tier = "#FEF3C7", outcome = "#BBF7D0", doomed = "#FCA5A5"
    ),
    status = c(
      ok = "#15803D", violated = "#B91C1C", relaxed = "#B45309",
      neutral = "#57534E", unknown = "#A8A29E", start = "#57534E", relax = "#9A3412"
    ),
    edge = c(
      flow = "#78716C", pass = "#15803D", fail = "#B45409",
      relax = "#B45309", doomed = "#991B1B", check = "#78716C",
      rank = "#57534E"
    ),
    edge_lty = c(
      flow = NA, pass = NA, fail = "3 2",
      relax = "5 2.5", doomed = "2 2", check = NA, rank = NA
    ),
    band = "#F5F5F4",
    text = "#1C1917",
    band_text = "#44403C"
  )
}

# a colour for `key` from a named palette, or a stable neutral when absent
pal_colour <- function(pal, key, default = "#57534E") {
  hit <- pal[names(pal) == as.character(key)]
  if (length(hit) && nzchar(hit[[1L]]) && !is.na(hit[[1L]])) hit[[1L]] else default
}

#' Draw a decision graph to SVG
#'
#' Structure in, figure out.  The graph contributes nodes, their group, their
#' kind and their status; the style contributes ink and page; the solver
#' contributes every coordinate and every bend.  The result is a self-contained
#' `<svg>` string (or a file, with `out`), so it scales, its text stays
#' selectable, and nothing needs a graphics device.
#'
#' Structure only, deliberately: a caller may hand this any data frame with
#' `name`, `label`, `stage`, `kind` (see `decision_solve()`) - it is how the
#' figures of the decision graphs are drawn outside a KASP run.
#'
#' The page is sized from the solved drawing, so it keeps the solver's near
#' square aspect rather than a fixed canvas.  When the layout engine is missing
#' or fails, a warning names why and a simple top-down placement is used (boxes
#' still measured, edges drawn as rails between rows): a render never errors.
#' @param graph a `tbl_graph` from [marker_decision_graph()].
#' @param style a [decision_style()].
#' @param out optional path; when given the SVG is written there and the path
#'   returned invisibly.
#' @return an SVG string, or invisibly `out`.
#' @export
decision_render <- function(graph, style = decision_style(), out = NULL) {
  solved <- decision_solve(graph, style)
  if (is.null(solved)) {
    message("drawing the decision graph without the layout engine (fallback placement)")
    solved <- decision_fallback(graph, style)
  }
  svg <- decision_svg(solved, style)
  if (is.null(out)) {
    return(invisible(svg))
  }
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  writeLines(svg, out, useBytes = TRUE)
  invisible(out)
}

#' The figure: title block, then the solved drawing, then the legend
#' @keywords internal
#' @noRd
decision_svg <- function(solved, style) {
  # the boxes come back with everything needed to draw them - position, size,
  # the wrapped lines, kind, status - so the renderer never joins two tables to
  # find a coordinate
  nodes <- solved$nodes
  pal <- style$palette
  # the solver returns y DOWN, which is SVG's own direction: no flip anywhere
  # in this file, only the offset that times the drawing after the title block
  head <- decision_svg_head(style)
  dy <- attr(head, "height")
  body <- character(0)

  # -- groups: one band per stage, painted under everything ------------------
  gx <- solved$groups
  if (!is.null(gx) && nrow(gx)) {
    for (i in seq_len(nrow(gx))) {
      body <- c(body, svg_rect(gx$x[i], gx$y[i] + dy, gx$w[i], gx$h[i],
        fill = pal$band, stroke = "#E7E5E4", sw = 0.5, rx = 4
      ))
      # the stage name in the gap above its own band (the root layout reserves
      # the 70 pt between bands for exactly this)
      body <- c(body, svg_text(gx$x[i] + style$pad_x, gx$y[i] + dy - 9, gx$id[i],
        size = style$font_size - 0.5, line_height = style$font_size,
        fill = pal$band_text, anchor = "start", bold = TRUE
      ))
    }
  }

  # -- edges: one <path> per edge, exactly as routed ------------------------
  ed <- solved$edges
  if (!is.null(ed) && nrow(ed)) {
    for (i in seq_len(nrow(ed))) {
      pts <- ed$points[[i]]
      if (is.null(pts) || nrow(pts) < 2L) next
      kind <- if (ed$kind[i] %in% names(pal$edge)) ed$kind[i] else "flow"
      dash <- pal_colour(pal$edge_lty, ed$kind[i], NA_character_)
      body <- c(body, svg_path(
        cbind(pts[, 1], pts[, 2] + dy),
        stroke = pal_colour(pal$edge, kind),
        sw = 0.75,
        dash = if (is.na(dash)) NULL else dash,
        marker = paste0("url(#arrow-", kind, ")")
      ))
    }
  }

  # -- boxes: fill by kind, border by status -------------------------------
  for (i in seq_len(nrow(nodes))) {
    n <- nodes[i, ]
    hard <- isTRUE(n$hard[[1L]])
    body <- c(body, svg_rect(
      n$x[[1L]], n$y[[1L]] + dy, n$w[[1L]], n$h[[1L]],
      fill = pal_colour(pal$kind, n$kind[[1L]], "#F5F5F4"),
      stroke = pal_colour(pal$status, n$status[[1L]], pal$status[["neutral"]]),
      sw = if (hard) 1.4 else 0.8, rx = 2
    ))
    # left-aligned inside the box: the solver sized the box from the widest
    # line it wrapped, so `pad_x` from the left edge is the whole inset
    body <- c(body, svg_text(
      n$x[[1L]] + style$pad_x, n$y[[1L]] + n$h[[1L]] / 2 + dy,
      n$lines[[1L]],
      size = style$font_size,
      line_height = style$line_height * style$font_size,
      fill = pal$text, anchor = "start", bold = hard
    ))
  }

  # -- legend and caption, under the drawing ------------------------------
  # both use the page inset the engine left around the boxes, so nothing sits
  # on the canvas edge
  inset <- style$margin %||% 14
  legend <- decision_svg_legend(
    style, inset, solved$height + dy, solved$width - 2 * inset
  )
  legend_h <- attr(legend, "height")
  caption <- decision_svg_caption(style, solved$height + dy + legend_h, inset)
  caption_h <- attr(caption, "height")

  svg_document(
    solved$width, solved$height + dy + legend_h + caption_h,
    c(head, body, legend, caption),
    defs = svg_arrow_defs(pal$edge),
    title = style$title
  )
}

# title + subtitle above the drawing; attr(,"height") is the space they take
decision_svg_head <- function(style) {
  out <- character(0)
  h <- 0
  if (!is.null(style$title) && nzchar(style$title)) {
    out <- c(out, svg_text(0, h + 11, style$title,
      size = 12, line_height = 14,
      fill = style$palette$text, anchor = "start", bold = TRUE
    ))
    h <- h + 17
  }
  if (!is.null(style$subtitle) && nzchar(style$subtitle)) {
    out <- c(out, svg_text(0, h + 7, style$subtitle,
      size = 8, line_height = 10,
      fill = "#57534E", anchor = "start"
    ))
    h <- h + 13
  }
  if (h > 0) h <- h + 4
  structure(out, height = h)
}

# caption under the drawing; attr(,"height") is the space it takes
decision_svg_caption <- function(style, y, x = 0) {
  if (is.null(style$caption) || !nzchar(style$caption)) {
    return(structure(character(0), height = 0))
  }
  structure(
    svg_text(x, y + 10, style$caption,
      size = 7.5, line_height = 9,
      fill = "#78716C", anchor = "start"
    ),
    height = 10 + 9 + 14
  )
}

# the legend: kind swatches, then status borders, then edge samples
decision_svg_legend <- function(style, x, y, width) {
  if (!isTRUE(style$legend)) {
    return(structure(character(0), height = 0))
  }
  pal <- style$palette
  items <- list()
  for (k in c("start", "check", "gate", "decision", "relax", "merge", "tier", "outcome", "doomed")) {
    items[[length(items) + 1L]] <- list(label = k, fill = pal_colour(pal$kind, k), stroke = "#78716C")
  }
  for (k in c("ok", "violated", "relaxed", "unknown")) {
    items[[length(items) + 1L]] <- list(label = k, fill = "none", stroke = pal_colour(pal$status, k))
  }
  items[[length(items) + 1L]] <- list(label = "hard gate", fill = "none", stroke = "#1C1917")
  for (k in c("pass", "fail", "relax", "rank", "doomed")) {
    items[[length(items) + 1L]] <- list(
      label = paste0(k, " edge"), stroke = pal_colour(pal$edge, k),
      dash = !is.na(pal_colour(pal$edge_lty, k, NA_character_))
    )
  }
  svg_legend(items, x, y + 4, width)
}

#' Place a graph without the engine: groups as columns, nodes as stacked rows
#'
#' The fallback used when `node` or the vendored modules are missing.  It knows
#' only the structure (group order, the chain order the schema declares), measures
#' nothing (the widths come from the same measurement the solver uses, via
#' `decision_measure_text()`), and draws every edge as a straight rail between
#' box centres.  It exists so a render never errors, not to be pretty.
#' @keywords internal
#' @noRd
decision_fallback <- function(graph, style) {
  n <- decision_graph_nodes(graph)
  e <- decision_edges(graph)
  # a crude advance: 0.5 em per character is only used when nothing can measure
  meas <- decision_measure_text(n$label, style$font_size, style$wrap)
  n$lines <- meas$lines
  n$w <- meas$width + 2 * style$pad_x
  n$h <- vapply(n$lines, length, integer(1)) * style$line_height * style$font_size + 2 * style$pad_y
  n$x <- rep(NA_real_, nrow(n))
  n$y <- rep(NA_real_, nrow(n))

  groups <- unique(n$group)
  gap_x <- 90
  x <- 0
  group_boxes <- list()
  for (g in groups) {
    idx <- which(n$group == g)
    # keep the schema's order: the rows are the chain the text export walks
    col <- n[idx, ]
    col_h <- sum(col$h) + 14 * (length(idx) - 1L)
    y <- 0
    for (k in seq_len(nrow(col))) {
      n$x[idx[k]] <- x
      n$y[idx[k]] <- y
      y <- y + col$h[k] + 14
    }
    group_boxes[[length(group_boxes) + 1L]] <- tibble::tibble(
      id = g, x = x - 8, y = -18, w = max(col$w) + 16, h = col_h + 26
    )
    x <- x + max(col$w) + 16 + gap_x
  }
  # centre of one box, looked up by name
  pos <- function(id) {
    k <- match(id, n$name)
    c(n$x[[k]] + n$w[[k]] / 2, n$y[[k]] + n$h[[k]] / 2)
  }
  routes <- lapply(seq_len(nrow(e)), function(i) {
    a <- pos(e$from[i])
    b <- pos(e$to[i])
    # leave from the right side of `a`, enter the left side of `b`
    list(
      from = e$from[i], to = e$to[i], kind = e$kind[i],
      points = cbind(
        c(a[1] + n$w[[match(e$from[i], n$name)]] / 2, b[1] - n$w[[match(e$to[i], n$name)]] / 2),
        c(a[2], b[2])
      )
    )
  })
  # the same shape `decision_solve()` returns: a data frame of edges whose
  # `points` is a list column, so the drawing code has one shape to read
  list(
    width = x - gap_x + 20, height = max(n$y + n$h) + 20,
    nodes = n,
    edges = tibble::tibble(
      from = vapply(routes, `[[`, character(1), "from"),
      to = vapply(routes, `[[`, character(1), "to"),
      kind = vapply(routes, `[[`, character(1), "kind"),
      points = lapply(routes, `[[`, "points")
    ),
    groups = dplyr::bind_rows(group_boxes)
  )
}

#' Wrap and measure labels without the layout engine
#'
#' Uses `systemfonts` when it is installed and a character advance otherwise;
#' either way the answer is wrapped lines and their widest measured width, so a
#' fallback box is still sized from its text rather than from a constant.
#' @keywords internal
#' @noRd
decision_measure_text <- function(label, font_size, wrap) {
  wrap <- wrap %||% 96
  words <- strsplit(as.character(label), "[[:space:]]+")
  lines <- lapply(words, function(w) {
    w <- w[nzchar(w)]
    if (!length(w)) {
      return("")
    }
    out <- w[1]
    for (k in seq_along(w)[-1]) {
      cand <- paste(out[length(out)], w[k])
      wid <- if (requireNamespace("systemfonts", quietly = TRUE)) {
        systemfonts::string_width(cand, family = "sans", size = font_size, res = 72)
      } else {
        nchar(cand) * font_size * 0.5
      }
      if (wid <= wrap) out[length(out)] <- cand else out <- c(out, w[k])
    }
    out
  })
  width <- vapply(lines, function(l) {
    max(if (requireNamespace("systemfonts", quietly = TRUE)) {
      systemfonts::string_width(l, family = "sans", size = font_size, res = 72)
    } else {
      nchar(l) * font_size * 0.5
    })
  }, numeric(1))
  list(lines = lines, width = width)
}
