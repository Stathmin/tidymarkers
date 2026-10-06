# =============================================================================
#  kasp_decision.R  --  the whole criterion set as data, picture and text
# =============================================================================
#  The pipeline evaluates a candidate along one tree per marker chemistry
#  (see kasp_schemes.R): site facts, the AS selection (gates -> ranking, with a
#  bounded refinement branch), the AS x COM selection, an optional duplex
#  compatibility stage, and finally the acceptance checks (composition, structure,
#  heteroduplexes, off-targets, in-silico PCR) ending in a verdict.  Every node
#  is declared once in kasp_criteria() (kasp_criteria.R), so the graph, the text
#  export, the selector and kasp_design_checks() cannot drift apart.
#
#  The tree is a DAG that flows top-to-bottom: hard gates are never relaxed; a
#  soft-gate refinement branch descends and rejoins the main path at the ranking
#  (never bypassing it), with a DOOMED terminal when a hard gate blocks or the
#  floors cannot be met.
# =============================================================================

#' The selection stages of a scheme, in pipeline order
#'
#' The allele-specific stages (AS primer, then the AS x COM pair) are shared by
#' every chemistry in the KASP family; a duplex chemistry appends its
#' compatibility stage, ASQ its cassette stage, and the marker-agnostic
#' `"primer"` scheme has none (there is no allele to discriminate).
#' @param scheme a scheme id (see [kasp_schemes()]).
#' @keywords internal
#' @noRd
kasp_decision_stages <- function(scheme = "snp") {
  info <- kasp_scheme(scheme)
  stages <- list()
  if (info$family %in% c("kasp", "multiplex", "asq")) {
    stages$as <- list(
      id = "as", label = "AS primer",
      gates = KASP_GATES_AS, relax = KASP_RELAX_AS, tiers = KASP_TIERS_AS
    )
    stages$pair <- list(
      id = "pair", label = "AS x COM pair",
      gates = KASP_GATES_PAIR, relax = KASP_RELAX_PAIR, tiers = KASP_TIERS_PAIR
    )
  }
  if (info$family == "multiplex") {
    stages$duplex <- list(
      id = "duplex", label = "duplex compatibility",
      gates = KASP_GATES_DUPLEX, relax = KASP_RELAX_DUPLEX, tiers = KASP_TIERS_DUPLEX
    )
  }
  if (info$family == "asq") {
    stages$asq <- list(
      id = "asq", label = "ASQ cassette",
      gates = KASP_GATES_ASQ, relax = KASP_RELAX_ASQ, tiers = KASP_TIERS_ASQ
    )
  }
  if (info$family == "caps") {
    stages$caps <- list(
      id = "caps", label = "CAPS digest",
      gates = KASP_GATES_CAPS, relax = KASP_RELAX_CAPS, tiers = KASP_TIERS_CAPS
    )
  }
  if (info$family == "arms") {
    stages$arms <- list(
      id = "arms", label = "ARMS bands",
      gates = KASP_GATES_ARMS, relax = KASP_RELAX_ARMS, tiers = KASP_TIERS_ARMS
    )
  }
  if (info$family == "hrm") {
    stages$hrm <- list(
      id = "hrm", label = "HRM amplicon",
      gates = KASP_GATES_HRM, relax = KASP_RELAX_HRM, tiers = KASP_TIERS_HRM
    )
  }
  stages
}

#' The decision tree of a scheme as a node/edge schema (no plotting dependency)
#'
#' Generated from [kasp_criteria()] filtered to `scheme`: the target *site*
#' checks (when the scheme has any), then each selection stage (gates down the
#' main path to "all gates pass?"; the fail edge descends the soft-gate ladder
#' and rejoins the ranking; a DOOMED terminal marks a hard block or an exhausted
#' floor), then the *assay* acceptance checks and a final verdict.  Hard gates
#' carry `hard = TRUE`.  A stage whose gates are all hard routes its fail edge
#' straight to DOOMED (no refinement branch).
#' @param scheme a scheme id (see [kasp_schemes()]).
#' @return a list with `nodes` (name, stage, kind, label, x, y, order, hard,
#'   detail) and `edges` (from, to, kind) tibbles.
#' @keywords internal
#' @noRd
kasp_decision_schema <- function(scheme = "snp") {
  info <- kasp_scheme(scheme)
  crit <- kasp_criteria(scheme = scheme)
  stages <- kasp_decision_stages(scheme)
  buf <- new.env(parent = emptyenv())
  buf$rows <- list()
  buf$edges <- list()
  push <- function(df) buf$rows[[length(buf$rows) + 1L]] <- df
  add_edge <- function(df) buf$edges[[length(buf$edges) + 1L]] <- df
  o <- 0
  refine_ids <- character(0)
  doomed_sources <- character(0)
  node <- function(name, stage, kind, label, order, x = 0, hard = FALSE, detail = NA_character_) {
    tibble::tibble(
      name = name, stage = stage, kind = kind, label = label,
      order = order, x = x, y = -order, hard = hard, detail = detail
    )
  }
  chain <- function(from, to, kind) tibble::tibble(from = from, to = to, kind = kind)

  # --- entry + target site checks (a pre-filter; absent for the bare primer) ---
  push(node("site:start", "site", "start", "candidate", o <- o + 1))
  entry <- "site:start"
  site <- crit[crit$stage == "site", ]
  if (nrow(site)) {
    push(node(site$id, "site", "check", site$label, o + seq_len(nrow(site)), hard = site$hard, detail = site$limit))
    add_edge(chain(entry, site$id[1L], "check"))
    if (nrow(site) > 1L) add_edge(chain(site$id[-nrow(site)], site$id[-1L], "check"))
    o <- o + nrow(site)
    entry <- site$id[nrow(site)]
  }

  # --- the scheme's selection stages -----------------------------------------
  for (st in stages) {
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
      refine = paste0(st$id, ":refine")
    )
    push(node(ids$gate, st$id, "gate", unname(KASP_GATE_SHORT[gates]), o + seq_along(gates),
      hard = hard, detail = unname(KASP_GATE_LIMIT[gates])
    ))
    if (length(gates) > 1L) add_edge(chain(ids$gate[-length(gates)], ids$gate[-1L], "flow"))
    add_edge(chain(entry, ids$gate[1L], "flow"))
    # the last gate flows into the decision node: without this the gate chain floats
    # and the "all gates pass?" node is disconnected from the gates it summarises
    add_edge(chain(ids$gate[length(gates)], ids$decision, "flow"))
    o <- o + length(gates)
    push(node(ids$decision, st$id, "decision", "all gates pass?", o <- o + 1))
    decision_order <- o
    soft <- length(relax) > 0L
    if (soft) {
      # one refinement node per stage: the soft-gate ladder in order (not one node
      # per rung - the ladder relaxes the *same* criteria, so a rung per gate would
      # just duplicate the gate nodes)
      push(node(
        ids$refine, st$id, "relax",
        paste0("refine: ", paste(KASP_GATE_SHORT[relax], collapse = " -> ")),
        decision_order + 1L,
        x = 0.85,
        detail = paste0(KASP_GATE_SHORT[relax], " [", KASP_GATE_FLOOR_LABEL[relax], "]", collapse = "; ")
      ))
      refine_ids <- c(refine_ids, ids$refine)
      merge_order <- decision_order + 2L
    } else {
      doomed_sources <- c(doomed_sources, ids$decision)
      merge_order <- decision_order + 1L
    }
    push(node(ids$merge, st$id, "merge", "survivors (bounds kept)", merge_order))
    has_tiers <- length(tiers) > 0L
    if (has_tiers) {
      push(node(ids$tier, st$id, "tier", unname(tiers), merge_order + seq_along(tiers)))
      best_order <- merge_order + length(tiers) + 1L
    } else {
      best_order <- merge_order + 1L
    }
    push(node(ids$best, st$id, "outcome", "best", best_order))
    o <- best_order
    add_edge(chain(ids$decision, ids$merge, "pass"))
    if (soft) {
      add_edge(chain(ids$decision, ids$refine, "fail"))
      add_edge(chain(ids$refine, ids$merge, "relax"))
    }
    if (has_tiers) {
      add_edge(chain(ids$merge, ids$tier[1L], "pass"))
      if (length(tiers) > 1L) add_edge(chain(ids$tier[-length(tiers)], ids$tier[-1L], "pass"))
      add_edge(chain(ids$tier[length(tiers)], ids$best, "pass"))
    } else {
      add_edge(chain(ids$merge, ids$best, "pass"))
    }
    entry <- ids$best
  }
  # one shared DOOMED terminal (a hard gate blocked, or every soft floor is spent)
  if (length(stages)) {
    push(node("doomed", "doomed", "doomed", "DOOMED", o <- o + 1, x = 0.85))
    add_edge(chain(c(refine_ids, doomed_sources), "doomed", "doomed"))
  }

  # --- the assay acceptance checks + verdict ---------------------------------
  chk <- crit[crit$stage == "assay", ]
  if (nrow(chk)) {
    push(node(chk$id, "assay", "check", chk$label, o + seq_len(nrow(chk)), hard = chk$hard, detail = chk$limit))
    o <- o + nrow(chk)
    add_edge(chain(entry, chk$id[1L], "check"))
    if (nrow(chk) > 1L) add_edge(chain(chk$id[-nrow(chk)], chk$id[-1L], "check"))
    entry <- chk$id[nrow(chk)]
  }
  push(node("assay:verdict", "assay", "outcome", info$verdict, o <- o + 1))
  add_edge(chain(entry, "assay:verdict", "check"))

  list(nodes = dplyr::bind_rows(buf$rows), edges = dplyr::bind_rows(buf$edges))
}

#' The decision tree as a `tidygraph`, annotated with an assay's outcome
#'
#' @param assay optional assembled assay; when given, gate nodes are stamped from
#'   `assay$selection` (`"relaxed"`/`"violated"`/`"ok"`) and check nodes from
#'   `assay$evaluation` (`"ok"`/`"violated"`).
#' @param scheme a scheme id (see [kasp_schemes()]); ignored when `assay` carries
#'   its own `$decision` graph.
#' @return a `tidygraph::tbl_graph`; nodes carry `name`, `stage`, `kind`,
#'   `label`, `x`, `y`, `hard`, `status`.
#' @examples
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' kasp_decision_graph(assay)
#' }
#' @family KASP checks
#' @export
kasp_decision_graph <- function(assay = NULL, scheme = "snp") {
  rlang::check_installed("tidygraph", reason = gettext("for the decision graph"))
  s <- kasp_decision_schema(scheme)
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
#' checks, each selection stage (hard gates never relaxed, the refinement ladder
#' rejoining the ranking, a DOOMED terminal) and the assay acceptance checks
#' ending in a verdict.  Hard limits are drawn with a **bold border**; with an
#' assay the traversed nodes are coloured by outcome (green enforced, amber
#' relaxed, red violated).  One marker chemistry is drawn at a time
#' ([kasp_schemes()]).
#' @param x an assembled assay, a `kasp_decision_graph()` result, or `NULL`.
#' @param scheme a scheme id (see [kasp_schemes()]) for the bare schema.
#' @return a `ggraph`/`ggplot` object.
#' @examples
#' \dontrun{
#' plot_decision_tree() # the schema
#' plot_decision_tree(scheme = "duplex") # the multiplex policy
#' plot_decision_tree(assay) # this assay's traversal
#' }
#' @family KASP checks
#' @export
plot_decision_tree <- function(x = NULL, scheme = "snp") {
  rlang::check_installed("ggraph", reason = gettext("to draw the decision tree"))
  g <- if (is.null(x)) {
    kasp_decision_graph(scheme = scheme)
  } else if (inherits(x, "kasp_assay")) {
    x$decision %||% kasp_decision_graph(x, scheme = scheme)
  } else {
    x
  }
  title <- if (is.null(x)) {
    paste0(kasp_scheme(scheme)$label, " - decision tree (bold = hard limit)")
  } else {
    "KASP decision tree (bold = hard limit)"
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
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11))
}

#' The decision tree of a scheme as text
#'
#' A plain-text, greppable export of the same registry the graph is drawn from:
#' the scheme's facts, the target-site checks, each selection stage's gates (with
#' hardness and floors), ladder, tiers and edges, the assay acceptance checks, and
#' the verdict.
#' @param scheme a scheme id (see [kasp_schemes()]).
#' @return a character vector, one element per line.
#' @examples
#' cat(kasp_decision_text(), sep = "\n")
#' cat(kasp_decision_text("duplex"), sep = "\n")
#' @family KASP checks
#' @export
kasp_decision_text <- function(scheme = "snp") {
  info <- kasp_scheme(scheme)
  crit <- kasp_criteria(scheme = scheme)
  s <- kasp_decision_schema(scheme)
  n <- s$nodes
  e <- s$edges
  lab <- rlang::set_names(n$label, n$name)

  # one formatted line per check: a hard marker (or `soft`) then label and limit
  check_lines <- function(cc, width, soft = "check") {
    purrr::pmap_chr(
      list(cc$hard, cc$label, cc$limit),
      function(hard, label, limit) {
        mark <- if (hard) "! HARD" else soft
        glue::glue(
          "  {stringr::str_pad(mark, 7, side = 'right')} {stringr::str_pad(label, width, side = 'right')} {limit}"
        )
      }
    )
  }
  section <- function(stage, title) {
    cc <- crit[crit$stage == stage, ]
    c(glue::glue("== {title} =="), check_lines(cc, 34L))
  }
  stage_block <- function(st) {
    cc <- crit[crit$stage == st$id, ]
    gate <- cc[cc$role == "gate", ]
    res <- c(
      glue::glue("== stage {st$id} ({st$label}) =="),
      "gates (hard = never relaxed):",
      check_lines(gate, 22L, soft = "- soft")
    )
    if (length(st$relax)) {
      res <- c(res, "relaxation ladder (least grave first):", purrr::map_chr(
        seq_along(st$relax),
        function(i) {
          gl <- st$relax[i]
          glue::glue(
            "  {i}. {stringr::str_pad(KASP_GATE_SHORT[[gl]], 22, side = 'right')} floor: {KASP_GATE_FLOOR_LABEL[[gl]]}"
          )
        }
      ))
    }
    if (length(st$tiers)) {
      res <- c(res, "tiers (lexicographic, per-tier bottleneck min):", purrr::map_chr(
        seq_along(st$tiers), function(i) glue::glue("  {i}. {unname(st$tiers[i])}")
      ))
    }
    ids <- n$name[n$stage == st$id]
    ee <- e[e$from %in% ids, ]
    edges <- purrr::pmap_chr(list(ee$from, ee$to, ee$kind), function(from, to, kind) {
      f <- lab[[from]] %||% from
      t <- lab[[to]] %||% to
      glue::glue(
        "  {stringr::str_pad(f, 24, side = 'right')} -{stringr::str_pad(kind, 6, side = 'right')}-> {t}"
      )
    })
    c(res, "edges:", edges, "")
  }

  out <- c(
    glue::glue("{info$label} - decision tree"), "",
    glue::glue(
      "discrimination: {info$discrimination} | readout: {info$readout}",
      " | zygosity: {info$zygosity} | loci: {info$loci} | amplicon: {info$amplicon}"
    ), ""
  )
  if (any(crit$stage == "site")) {
    out <- c(out, section("site", "target site checks"), "")
  }
  out <- c(out, purrr::map(kasp_decision_stages(scheme), stage_block) |> purrr::list_c())
  out <- c(
    out, section("assay", "assay acceptance checks"), "",
    "DOOMED: a hard gate blocked, or every soft floor is exhausted.",
    glue::glue("verdict: all checks pass -> {info$verdict}; else see the failing node")
  )
  as.character(out)
}
