# =============================================================================
#  marker_decision_svg.R  --  SVG string builders
# =============================================================================
#  Every coordinate written here came out of the solver
#  (`marker_decision_solver.R`).  These helpers never place anything: they
#  format points, rectangles, text and paths into the markup the figure needs.
# =============================================================================

# a number as SVG wants it: at most three decimals, no exponent, no "-0"
svg_num <- function(x) {
  s <- formatC(as.numeric(x), format = "f", digits = 3, drop0trailing = TRUE)
  s[s == "-0"] <- "0"
  s
}

svg_attr <- function(...) {
  parts <- list(...)
  keep <- !vapply(parts, is.null, logical(1))
  parts <- parts[keep]
  if (!length(parts)) {
    return("")
  }
  paste0(" ", paste0(names(parts), '="', unlist(parts), '"', collapse = " "))
}

svg_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

#' One SVG document from `<defs>` and body fragments
#'
#' @param width,height page size in points.
#' @param body character vector of markup, one element per element.
#' @param defs character vector of `<defs>` content (arrow markers).
#' @return a single string.
#' @keywords internal
#' @noRd
svg_document <- function(width, height, body, defs = character(0), title = NULL) {
  c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    paste0(
      '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" ',
      'width="', svg_num(width), '" height="', svg_num(height), '" ',
      'viewBox="0 0 ', svg_num(width), " ", svg_num(height), '">'
    ),
    if (!is.null(title)) paste0("<title>", svg_escape(title), "</title>"),
    if (length(defs)) c("<defs>", defs, "</defs>"),
    # explicit page-size paper: a percentage rect is legal SVG but some parsers
    # (xml2 among them) refuse it, and the page size is known here anyway
    paste0('<rect x="0" y="0" width="', svg_num(width), '" height="', svg_num(height), '" fill="#FFFFFF"/>'),
    body,
    "</svg>"
  )
}

svg_rect <- function(x, y, w, h, fill = "none", stroke = "none", sw = 1,
                     rx = NULL, dash = NULL, opacity = NULL) {
  paste0(
    "<rect", svg_attr(
      x = svg_num(x), y = svg_num(y), width = svg_num(w), height = svg_num(h),
      rx = if (is.null(rx)) NULL else svg_num(rx),
      fill = fill, stroke = stroke, `stroke-width` = svg_num(sw),
      `stroke-dasharray` = dash, opacity = opacity
    ), "/>"
  )
}

svg_line <- function(x1, y1, x2, y2, stroke = "#000000", sw = 0.7, dash = NULL,
                     opacity = NULL) {
  paste0(
    "<line", svg_attr(
      x1 = svg_num(x1), y1 = svg_num(y1), x2 = svg_num(x2), y2 = svg_num(y2),
      stroke = stroke, `stroke-width` = svg_num(sw),
      `stroke-dasharray` = dash, opacity = opacity
    ), "/>"
  )
}

#' A text run, or several lines centred in a box
#'
#' `y` is the vertical CENTRE of the whole block and `x` the left edge when
#' `anchor = "start"`.  Line offsets come from the style's measured line height.
#' @keywords internal
#' @noRd
svg_text <- function(x, y, lines, size = 9, line_height = 11.7, fill = "#1C1917",
                     anchor = "middle", bold = FALSE, italic = FALSE) {
  n <- length(lines)
  # first baseline sits half a line above centre, the rest step down
  y0 <- y - (n - 1) * line_height / 2 + size * 0.35
  anchor_attr <- if (anchor == "start") "start" else "middle"
  paste0(
    "<text", svg_attr(
      x = svg_num(x), y = svg_num(y0), `text-anchor` = anchor_attr,
      `font-family` = "Liberation Sans, Arial, Helvetica, sans-serif",
      `font-size` = svg_num(size), fill = fill,
      `font-weight` = if (bold) "bold" else NULL,
      `font-style` = if (italic) "italic" else NULL
    ), ">",
    # one <tspan> per line, each with the text's own x and exactly one y
    paste0(
      "<tspan", sapply(
        svg_num(y0 + (seq_len(n) - 1) * line_height),
        function(yy) svg_attr(x = svg_num(x), y = yy)
      ), ">", svg_escape(lines), "</tspan>",
      collapse = ""
    ),
    "</text>"
  )
}

#' One routed edge as a single `<path>`
#'
#' The route arrives as points that already turn (libavoid emits orthogonal
#' polylines), so this only joins them - it never inserts a bend.  The arrow
#' head is a marker, so the path itself ends exactly where the router put it.
#' @keywords internal
#' @noRd
svg_path <- function(points, stroke = "#57534E", sw = 0.7, dash = NULL,
                     marker = NULL, opacity = NULL) {
  if (nrow(points) < 2L) {
    return(character(0))
  }
  d <- paste0(
    "M", svg_num(points[1, 1]), " ", svg_num(points[1, 2]),
    paste0(" L", svg_num(points[-1, 1]), " ", svg_num(points[-1, 2]), collapse = "")
  )
  paste0(
    "<path", svg_attr(
      d = d, fill = "none", stroke = stroke, `stroke-width` = svg_num(sw),
      `stroke-linejoin` = "miter", `stroke-linecap` = "butt",
      `stroke-dasharray` = dash, `marker-end` = marker, opacity = opacity
    ), "/>"
  )
}

#' Lay out one line of legend items left to right, wrapping at `width`
#'
#' Each `item` is a list with `label`, and either `fill`+`stroke` (a swatch) or
#' `stroke`+`dash` (a line sample).  Returns markup; `attr(, "height")` is the
#' space it took.
#' @keywords internal
#' @noRd
svg_legend <- function(items, x, y, width, size = 6.5, step = 4) {
  out <- character(0)
  cx <- x
  row <- 0
  row_h <- size * 1.9
  for (it in items) {
    w_lab <- nchar(it$label) * size * 0.52 + 16
    if (cx > x && cx + w_lab > x + width) {
      cx <- x
      row <- row + 1
    }
    yy <- y + row * row_h + row_h / 2
    if (!is.null(it$fill)) {
      out <- c(out, svg_rect(cx, yy - 4, 9, 8,
        fill = it$fill, stroke = it$stroke, sw = 0.8, rx = 1.5
      ))
    } else {
      out <- c(out, svg_line(cx, yy, cx + 10, yy,
        stroke = it$stroke, sw = 0.9,
        dash = if (isTRUE(it$dash)) "5 2.5" else NULL
      ))
    }
    out <- c(out, svg_text(cx + 13, yy, it$label,
      size = size, line_height = size * 1.2,
      anchor = "start", fill = "#44403C"
    ))
    cx <- cx + w_lab + step
  }
  # height covers the tallest row's text plus the swatch half-height and the
  # air below it, so a following block never lands on the legend
  structure(out, height = (row + 1) * row_h + size)
}

#' The arrow-head markers, one per edge colour
#'
#' @param colours named vector of colours (the edge-kind palette).
#' @return `<marker>` elements, keyed `url(#arrow-<name>)`.
#' @keywords internal
#' @noRd
svg_arrow_defs <- function(colours, size = 5) {
  ids <- names(colours)
  vapply(seq_along(colours), function(i) {
    paste0(
      '<marker id="arrow-', ids[[i]], '" viewBox="0 0 10 10" refX="9" refY="5" ',
      'markerWidth="', svg_num(size), '" markerHeight="', svg_num(size), '" ',
      'markerUnits="userSpaceOnUse" orient="auto-start-reverse">',
      '<path d="M 0 0 L 10 5 L 0 10 z" fill="', colours[[i]], '"/>',
      "</marker>"
    )
  }, character(1))
}
