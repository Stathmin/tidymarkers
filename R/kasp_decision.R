# =============================================================================
#  kasp_decision.R  --  the selection policy as a graph (data, then picture)
# =============================================================================
#  The selection engine (kasp_select.R) is: filter by hard gates, rank the
#  survivors lexicographically on tier margins, and - when nothing passes - relax
#  the gates one at a time (least grave first).  That policy is expressed *once*,
#  as data (`kasp_decision_schema()`), read straight from the same constants the
#  engine uses (KASP_GATES_*, KASP_RELAX_*, KASP_TIERS_*, KASP_GATE_LABEL).  Add a
#  gate, a tier, or a rung of the relaxation ladder and both the engine and the
#  picture change together - the graph cannot drift from the code.
#
#  `kasp_decision_graph()` turns the schema into a `tidygraph` (nodes + edges,
#  laid out top-to-bottom with the relaxation ladder to the side) and, given an
#  assay, stamps each node with its runtime status (ok / relaxed / violated).
#  `plot_decision_tree()` draws it with `ggraph`.
# =============================================================================

# The two stages of the policy, in pipeline order, each naming its gates (in
# evaluation order), its relaxation ladder (least grave first) and its tiers.
kasp_decision_stages <- function() {
  list(
    list(
      id = "as", label = "AS primer",
      gates = KASP_GATES_AS, relax = KASP_RELAX_AS, tiers = KASP_TIERS_AS
    ),
    list(
      id = "pair", label = "AS x COM pair",
      gates = KASP_GATES_PAIR, relax = KASP_RELAX_PAIR, tiers = KASP_TIERS_PAIR
    )
  )
}

#' The decision tree as a node/edge schema (no plotting dependency)
#'
#' Reads the same constants the selection engine uses, so the graph always
#' reflects the policy.  It is a DAG that flows top-to-bottom: the main path runs
#' down the gates to the tier ranking, and the relaxation branch descends
#' alongside and rejoins the main path at the outcome (the least-bad candidate) -
#' so no edge ever loops back up.
#' @return a list with `nodes` (name, stage, kind, label, x, y, order) and
#'   `edges` (from, to, kind) tibbles.
#' @keywords internal
#' @noRd
kasp_decision_schema <- function() {
  nodes <- list()
  edges <- list()
  o <- 0 # running main-path order across stages
  for (st in kasp_decision_stages()) {
    gates <- st$gates
    tiers <- st$tiers
    relax <- st$relax
    ids <- list(
      start = paste0(st$id, ":start"),
      gate = paste0(st$id, ":gate:", gates),
      decision = paste0(st$id, ":decision"),
      tier = paste0(st$id, ":tier:", names(tiers)),
      best = paste0(st$id, ":best"),
      relax = paste0(st$id, ":relax:", relax)
    )
    row <- function(name, kind, label, order, x = 0) {
      tibble::tibble(
        name = name, stage = st$id, kind = kind, label = label,
        order = order, x = x, y = -order
      )
    }
    rows <- list(row(ids$start, "start", "candidates", o <- o + 1))
    for (i in seq_along(gates)) {
      rows <- c(rows, list(row(ids$gate[i], "gate", paste0("gate: ", KASP_GATE_SHORT[gates[i]]), o <- o + 1)))
    }
    rows <- c(rows, list(row(ids$decision, "decision", "all gates pass?", o <- o + 1)))
    decision_order <- o
    # tiers and the relaxation ladder share the band below the decision, so the
    # ladder's last rung still sits above the outcome it rejoins
    for (i in seq_along(tiers)) {
      rows <- c(rows, list(row(ids$tier[i], "tier", paste0("rank ", i, ": ", unname(tiers[i])), decision_order + i)))
    }
    for (i in seq_along(relax)) {
      rows <- c(rows, list(row(
        ids$relax[i], "relax", paste0("relax ", i, ": ", KASP_GATE_SHORT[relax[i]]),
        decision_order + i,
        x = 0.85
      )))
    }
    best_order <- decision_order + max(length(tiers), length(relax)) + 1L
    rows <- c(rows, list(row(ids$best, "outcome", "best (least-bad)", best_order)))
    o <- best_order
    nodes[[length(nodes) + 1L]] <- dplyr::bind_rows(rows)

    chain <- function(from, to, kind) tibble::tibble(from = from, to = to, kind = kind)
    stage_edges <- list(
      chain(ids$start, ids$gate[1L], "flow"),
      chain(ids$gate[-length(gates)], ids$gate[-1L], "flow"),
      chain(ids$gate[length(gates)], ids$decision, "flow"),
      chain(ids$decision, ids$tier[1L], "pass"),
      chain(ids$tier[-length(tiers)], ids$tier[-1L], "pass"),
      chain(ids$tier[length(tiers)], ids$best, "pass"),
      # the refinement branch: relax the least-grave gate, retry down the ladder,
      # then rejoin the main path at the outcome
      chain(ids$decision, ids$relax[1L], "fail"),
      chain(ids$relax[-length(relax)], ids$relax[-1L], "relax"),
      chain(ids$relax[length(relax)], ids$best, "relax")
    )
    edges[[length(edges) + 1L]] <- dplyr::bind_rows(stage_edges)
  }
  # the AS pipeline feeds the pair pipeline
  edges[[length(edges) + 1L]] <- tibble::tibble(from = "as:best", to = "pair:start", kind = "flow")
  list(nodes = dplyr::bind_rows(nodes), edges = dplyr::bind_rows(edges))
}

#' The decision tree as a `tidygraph`, annotated with an assay's outcome
#'
#' @param assay optional assembled assay; when given, each node is stamped with
#'   its runtime status from `assay$selection` (`"relaxed"`, `"violated"` or
#'   `"ok"`).
#' @return a `tidygraph::tbl_graph`; nodes carry `name`, `stage`, `kind`,
#'   `label`, `x`, `y`, `status`.
#' @examples
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' kasp_decision_graph(assay)
#' }
#' @family KASP checks
#' @export
kasp_decision_graph <- function(assay = NULL) {
  rlang::check_installed("tidygraph", reason = gettext("for the decision graph"))
  s <- kasp_decision_schema()
  nodes <- s$nodes
  nodes$status <- "neutral"
  if (!is.null(assay) && !is.null(assay$selection)) {
    sel <- assay$selection
    stamp <- function(prefix, relaxed, violated, kind) {
      keep <- nodes$kind == kind & nodes$stage == prefix
      nodes$status[keep & nodes$name %in% paste0(prefix, ":", kind, ":", relaxed)] <<- "relaxed"
      nodes$status[keep & nodes$name %in% paste0(prefix, ":", kind, ":", violated)] <<- "violated"
      nodes$status[keep] <<- ifelse(
        nodes$status[keep] == "neutral", "ok", nodes$status[keep]
      )
    }
    stamp("pair", sel$relaxed %||% character(0), sel$violated %||% character(0), "gate")
    stamp("as", character(0), sel$as %||% character(0), "gate")
    nodes$status[nodes$kind %in% c("tier", "decision", "outcome")] <- "ok"
    nodes$status[nodes$kind == "start"] <- "start"
    nodes$status[nodes$kind == "relax"] <- "relax"
  }
  tidygraph::tbl_graph(nodes = nodes, edges = s$edges, directed = TRUE, node_key = "name")
}

#' Draw the selection decision tree
#'
#' Renders the graph from [kasp_decision_graph()] with `ggraph`: gates are boxes
#' down the main path, the "all gates pass?" decision splits into the tier
#' ranking (pass, green) or the relaxation ladder on the right (fail, amber/red).
#' With an assay the traversed nodes are coloured by their outcome (green =
#' enforced, amber = relaxed, red = violated), so you can see step by step what
#' the engine looked at and where it compromised.
#' @param x an assembled assay, a `kasp_decision_graph()` result, or `NULL` for
#'   the bare schema.
#' @return a `ggraph`/`ggplot` object.
#' @examples
#' \dontrun{
#' plot_decision_tree() # the schema
#' plot_decision_tree(assay) # this assay's traversal
#' }
#' @family KASP checks
#' @export
plot_decision_tree <- function(x = NULL) {
  rlang::check_installed("ggraph", reason = gettext("to draw the decision tree"))
  g <- if (is.null(x)) {
    kasp_decision_graph()
  } else if (inherits(x, "kasp_assay")) {
    x$decision %||% kasp_decision_graph(x)
  } else {
    x
  }
  status_cols <- c(
    ok = KASP_STATUS[["ok"]], violated = KASP_STATUS[["broken"]],
    relaxed = "#e6ab02", neutral = KASP_STATUS[["neutral"]],
    start = "grey72", relax = "grey90"
  )
  ggraph::ggraph(g, layout = "manual", x = .data$x, y = .data$y) +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_linetype = .data$kind, edge_colour = .data$kind),
      width = 0.45, end_cap = ggraph::circle(2.4, "mm"),
      arrow = ggplot2::arrow(length = ggplot2::unit(1.8, "mm"), type = "closed")
    ) +
    ggraph::scale_edge_linetype_manual(
      values = c(flow = "solid", pass = "solid", fail = "dotted", relax = "dashed"),
      guide = "none"
    ) +
    ggraph::scale_edge_colour_manual(
      values = c(flow = "grey45", pass = "#1b7837", fail = "#d7301f", relax = "#b8860b"),
      guide = "none"
    ) +
    ggraph::geom_node_label(
      ggplot2::aes(label = .data$label, fill = .data$status),
      size = 2.4, label.r = ggplot2::unit(1.5, "mm"),
      label.padding = ggplot2::unit(1.3, "mm")
    ) +
    ggplot2::scale_fill_manual(values = status_cols, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.45)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 1.2)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "KASP selection decision tree") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11))
}
