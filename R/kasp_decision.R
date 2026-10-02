# =============================================================================
#  kasp_decision.R  --  the whole criterion set as data, picture and text
# =============================================================================
#  The pipeline evaluates a candidate along one tree: site facts, the AS
#  selection (gates -> ranking, with a bounded refinement branch), the AS x COM
#  selection, and finally the acceptance checks (composition, structure,
#  heteroduplexes, off-targets, in-silico PCR) ending in a verdict.  Every node
#  is declared once in kasp_criteria() (kasp_criteria.R), so the graph, the text
#  export, the selector and kasp_design_checks() cannot drift apart.
#
#  The tree is a DAG that flows top-to-bottom: hard gates are never relaxed; a
#  soft-gate refinement branch descends and rejoins the main path at the ranking
#  (never bypassing it), with a DOOMED terminal when a hard gate blocks or the
#  floors cannot be met.
# =============================================================================

#' The two selection stages of the policy, in pipeline order
#' @keywords internal
#' @noRd
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
#' Generated from [kasp_criteria()]: the target *site* checks, then each
#' selection stage (gates down the main path to "all gates pass?"; the fail edge
#' descends the soft-gate ladder and rejoins the ranking; a DOOMED terminal marks
#' a hard block or an exhausted floor), then the *assay* acceptance checks and a
#' final verdict.  Hard gates carry `hard = TRUE`.
#' @return a list with `nodes` (name, stage, kind, label, x, y, order, hard,
#'   detail) and `edges` (from, to, kind) tibbles.
#' @keywords internal
#' @noRd
kasp_decision_schema <- function() {
  crit <- kasp_criteria()
  buf <- new.env(parent = emptyenv())
  buf$rows <- list()
  buf$edges <- list()
  push <- function(df) buf$rows[[length(buf$rows) + 1L]] <- df
  add_edge <- function(df) buf$edges[[length(buf$edges) + 1L]] <- df
  o <- 0
  node <- function(name, stage, kind, label, order, x = 0, hard = FALSE, detail = NA_character_) {
    tibble::tibble(
      name = name, stage = stage, kind = kind, label = label,
      order = order, x = x, y = -order, hard = hard, detail = detail
    )
  }
  chain <- function(from, to, kind) tibble::tibble(from = from, to = to, kind = kind)

  # --- target site checks (a pre-filter) -------------------------------------
  site <- crit[crit$stage == "site", ]
  push(node("site:start", "site", "start", "candidate", o <- o + 1))
  push(node(site$id, "site", "check", site$label, o + seq_len(nrow(site)), hard = site$hard, detail = site$limit))
  add_edge(chain(c("site:start", site$id[-nrow(site)]), site$id, "check"))
  o <- o + nrow(site)

  # --- the two selection stages ----------------------------------------------
  for (st in kasp_decision_stages()) {
    gates <- st$gates
    tiers <- st$tiers
    relax <- st$relax
    hard <- !(gates %in% relax)
    ids <- list(
      gate = paste0(st$id, ":gate:", gates),
      decision = paste0(st$id, ":decision"),
      merge = paste0(st$id, ":survivors"),
      tier = paste0(st$id, ":tier:", names(tiers)),
      best = paste0(st$id, ":best"),
      doomed = paste0(st$id, ":doomed"),
      relax = paste0(st$id, ":relax:", relax)
    )
    push(node(ids$gate, st$id, "gate", unname(KASP_GATE_SHORT[gates]), o + seq_along(gates),
      hard = hard, detail = unname(KASP_GATE_LIMIT[gates])
    ))
    o <- o + length(gates)
    push(node(ids$decision, st$id, "decision", "all gates pass?", o <- o + 1))
    decision_order <- o
    push(node(ids$relax, st$id, "relax", unname(KASP_GATE_SHORT[relax]), decision_order + seq_along(relax),
      x = 0.85, detail = unname(KASP_GATE_FLOOR_LABEL[relax])
    ))
    merge_order <- decision_order + length(relax) + 1L
    push(node(ids$merge, st$id, "merge", "survivors (bounds kept)", merge_order))
    push(node(ids$tier, st$id, "tier", unname(tiers), merge_order + seq_along(tiers)))
    best_order <- merge_order + length(tiers) + 1L
    push(node(ids$best, st$id, "outcome", "best", best_order))
    push(node(ids$doomed, st$id, "doomed", "DOOMED", best_order, x = 0.85))
    o <- best_order
    add_edge(dplyr::bind_rows(
      chain(ids$gate[-length(gates)], ids$gate[-1L], "flow"),
      chain(ids$decision, ids$merge, "pass"),
      chain(ids$decision, ids$relax[1L], "fail"),
      chain(ids$relax[-length(relax)], ids$relax[-1L], "relax"),
      chain(ids$relax[length(relax)], ids$merge, "relax"),
      chain(ids$relax[length(relax)], ids$doomed, "doomed"),
      chain(ids$merge, ids$tier[1L], "pass"),
      chain(ids$tier[-length(tiers)], ids$tier[-1L], "pass"),
      chain(ids$tier[length(tiers)], ids$best, "pass")
    ))
  }
  # stage entries: the site checks feed the AS stage; the AS best feeds the pair
  add_edge(chain(site$id[nrow(site)], "as:gate:gate_discrimination", "check"))
  add_edge(chain("as:best", "pair:gate:gate_discrimination", "flow"))

  # --- the assay acceptance checks + verdict ---------------------------------
  chk <- crit[crit$stage == "assay", ]
  push(node(chk$id, "assay", "check", chk$label, o + seq_len(nrow(chk)), hard = chk$hard, detail = chk$limit))
  o <- o + nrow(chk)
  push(node("assay:verdict", "assay", "outcome", "verdict", o <- o + 1))
  add_edge(chain(c("pair:best", chk$id[-nrow(chk)]), chk$id, "check"))
  add_edge(chain(chk$id[nrow(chk)], "assay:verdict", "check"))

  list(nodes = dplyr::bind_rows(buf$rows), edges = dplyr::bind_rows(buf$edges))
}

#' The decision tree as a `tidygraph`, annotated with an assay's outcome
#'
#' @param assay optional assembled assay; when given, gate nodes are stamped from
#'   `assay$selection` (`"relaxed"`/`"violated"`/`"ok"`) and check nodes from
#'   `assay$evaluation` (`"ok"`/`"violated"`).
#' @return a `tidygraph::tbl_graph`; nodes carry `name`, `stage`, `kind`,
#'   `label`, `x`, `y`, `hard`, `status`.
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
    status <- nodes$status
    gate <- nodes$kind == "gate"
    status[gate] <- "ok"
    is_pair <- gate & nodes$stage == "pair"
    is_as <- gate & nodes$stage == "as"
    status[is_pair & nodes$name %in% paste0("pair:gate:", sel$relaxed %||% character(0))] <- "relaxed"
    status[is_pair & nodes$name %in% paste0("pair:gate:", sel$violated %||% character(0))] <- "violated"
    status[is_as & nodes$name %in% paste0("as:gate:", sel$as %||% character(0))] <- "violated"
    status[nodes$kind %in% c("tier", "decision", "merge", "outcome")] <- "ok"
    status[nodes$kind == "start"] <- "start"
    status[nodes$kind == "relax"] <- "relax"
    ev <- assay$evaluation
    if (!is.null(ev)) {
      bad <- ev$criterion[!ev$pass]
      status[nodes$kind == "check" & nodes$label %in% bad] <- "violated"
      status[nodes$kind == "check" & !nodes$label %in% bad] <- "ok"
    }
    nodes$status <- status
  }
  tidygraph::tbl_graph(nodes = nodes, edges = s$edges, directed = TRUE, node_key = "name")
}

#' Draw the selection decision tree
#'
#' Renders the graph from [kasp_decision_graph()] with `ggraph`: the target-site
#' checks, the two selection stages (hard gates never relaxed, the refinement
#' ladder rejoining the ranking, a DOOMED terminal) and the assay acceptance
#' checks ending in a verdict.  Hard limits are drawn with a **bold border**;
#' with an assay the traversed nodes are coloured by outcome (green enforced,
#' amber relaxed, red violated).
#' @param x an assembled assay, a `kasp_decision_graph()` result, or `NULL`.
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
      width = 0.4, end_cap = ggraph::circle(2.2, "mm"),
      arrow = ggplot2::arrow(length = ggplot2::unit(1.6, "mm"), type = "closed")
    ) +
    ggraph::scale_edge_linetype_manual(
      values = c(
        flow = "solid", pass = "solid", fail = "dotted",
        relax = "dashed", doomed = "solid", check = "solid"
      ),
      guide = "none"
    ) +
    ggraph::scale_edge_colour_manual(
      values = c(
        flow = "grey45", pass = "#1b7837", fail = "#d7301f",
        relax = "#b8860b", doomed = "#d7301f", check = "grey55"
      ),
      guide = "none"
    ) +
    ggraph::geom_node_label(
      ggplot2::aes(label = .data$label, fill = .data$status),
      size = 2.2, label.r = ggplot2::unit(1.2, "mm"),
      label.padding = ggplot2::unit(1.1, "mm"), linewidth = 0.2
    ) +
    ggraph::geom_node_label(
      data = function(d) d[d$hard, , drop = FALSE],
      ggplot2::aes(label = .data$label, fill = .data$status),
      size = 2.2, fontface = "bold", colour = "black", linewidth = 1.1,
      label.r = ggplot2::unit(1.2, "mm"), label.padding = ggplot2::unit(1.1, "mm")
    ) +
    ggplot2::scale_fill_manual(values = status_cols, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 1.2)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "KASP decision tree (bold = hard limit)") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11))
}

#' The decision tree as text
#'
#' A plain-text, greppable export of the same registry the graph is drawn from:
#' the target-site checks, each selection stage's gates (with hardness and
#' floors), ladder, tiers and edges, the assay acceptance checks, and the verdict.
#' @return a character vector, one element per line.
#' @examples
#' cat(kasp_decision_text(), sep = "\n")
#' @family KASP checks
#' @export
kasp_decision_text <- function() {
  crit <- kasp_criteria()
  s <- kasp_decision_schema()
  n <- s$nodes
  e <- s$edges
  lab <- stats::setNames(n$label, n$name)
  out <- c("KASP decision tree", "")
  section <- function(stage, title, kind) {
    cc <- crit[crit$stage == stage, ]
    res <- paste0("== ", title, " ==")
    if (kind == "check") {
      for (i in seq_len(nrow(cc))) {
        mark <- if (cc$hard[i]) "! HARD" else "check"
        res <- c(res, sprintf("  %-7s %-34s %s", mark, cc$label[i], cc$limit[i]))
      }
    }
    res
  }
  out <- c(out, section("site", "target site checks", "check"), "")
  for (st in kasp_decision_stages()) {
    cc <- crit[crit$stage == st$id, ]
    gate <- cc[cc$role == "gate", ]
    out <- c(out, paste0("== stage ", st$id, " (", st$label, ") =="))
    out <- c(out, "gates (hard = never relaxed):")
    for (i in seq_len(nrow(gate))) {
      mark <- if (gate$hard[i]) "! HARD" else "- soft"
      out <- c(out, sprintf("  %-7s %-22s %s", mark, gate$label[i], gate$limit[i]))
    }
    out <- c(out, "relaxation ladder (least grave first):")
    for (i in seq_along(st$relax)) {
      gl <- st$relax[i]
      out <- c(out, sprintf("  %d. %-22s floor: %s", i, KASP_GATE_SHORT[[gl]], KASP_GATE_FLOOR_LABEL[[gl]]))
    }
    out <- c(out, "tiers (lexicographic, per-tier bottleneck min):")
    for (i in seq_along(st$tiers)) out <- c(out, sprintf("  %d. %s", i, unname(st$tiers[i])))
    out <- c(out, "edges:")
    ids <- n$name[n$stage == st$id]
    ee <- e[e$from %in% ids, ]
    for (i in seq_len(nrow(ee))) {
      from <- lab[[ee$from[i]]] %||% ee$from[i]
      to <- lab[[ee$to[i]]] %||% ee$to[i]
      out <- c(out, sprintf("  %-24s -%-6s-> %s", from, ee$kind[i], to))
    }
    out <- c(out, "")
  }
  out <- c(out, section("assay", "assay acceptance checks", "check"))
  c(
    out, "",
    "DOOMED: a hard gate blocked, or every soft floor is exhausted.",
    "verdict: all checks pass -> usable assay; else see the failing node"
  )
}
