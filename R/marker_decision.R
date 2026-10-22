# =============================================================================
#  marker_decision.R  --  the whole criterion set as data, picture and text
# =============================================================================
#  The pipeline evaluates a candidate along one tree per marker chemistry
#  (see marker_schemes.R): site facts, the AS selection (gates -> ranking, with a
#  bounded refinement branch), the AS x COM selection, an optional duplex
#  compatibility stage, and finally the acceptance checks (composition, structure,
#  heteroduplexes, off-targets, in-silico PCR) ending in a verdict.  Every node
#  is declared once in marker_criteria() (marker_criteria.R), so the graph, the text
#  export, the selector and marker_design_checks() cannot drift apart.
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
#' compatibility stage, ASQ a documentation-only geometry stage, and the
#' marker-agnostic `"primer"` scheme has none (there is no allele to discriminate).
#'
#' Each stage declares *two* orders (C6), which need not agree: `gates` is the
#' order the method applies the gates (grave / cheap first - a fatal flaw is
#' tested before an expensive one), and `relax` is the ladder's order (least
#' grave first - the mildest trade-off is tried before a worse one).  Both are
#' labelled per stage in [marker_decision_text()].
#' @param scheme a scheme id (see [marker_schemes()]).
#' @keywords internal
#' @noRd
marker_decision_stages <- function(scheme = "snp") {
  info <- marker_scheme(scheme)
  # which stages the scheme runs is declared once, in marker_builtin_stages()
  # (R/marker_stages.R): the `as`/`pair` stages are shared by the KASP family
  # (`kasp`, `multiplex`, `asq`), every other stage belongs to its own family.
  builtins <- marker_builtin_stages()
  runs <- vapply(builtins, function(st) info$family %in% st$schemes, logical(1))
  stages <- lapply(builtins[runs], function(st) {
    list(
      id = st$scope, label = st$label,
      gates = marker_stage_constants(st, "gates"),
      relax = marker_stage_constants(st, "relax"),
      tiers = marker_stage_constants(st, "tiers")
    )
  })
  # the ids are the stage names, not the scopes (the shared stages carry the
  # scope "kasp")
  names(stages) <- names(builtins)[runs]
  for (nm in names(stages)) {
    stages[[nm]]$id <- nm
  }
  # discovered schemes (R/marker_<scheme>.R) contribute their own stage
  stages <- c(stages, extra_decision_stages(scheme))
  stages
}

#' The decision tree of a scheme as a node/edge schema (no plotting dependency)
#'
#' Generated from [marker_criteria()] filtered to `scheme`: the target *site*
#' checks (when the scheme has any), then each selection stage (gates down the
#' main path to "all gates pass?"; the fail edge descends the soft-gate ladder
#' and rejoins the ranking; a DOOMED terminal marks a hard block or an exhausted
#' floor), then the *assay* acceptance checks and a final verdict.  Hard gates
#' carry `hard = TRUE`.  A stage whose gates are all hard routes its fail edge
#' straight to DOOMED (no refinement branch); in a *mixed* stage each hard gate
#' routes to DOOMED itself, so a failure of it never looks relaxable (C4).  A
#' stage with no gates (a documentation-only stage, e.g. ASQ) contributes no gate
#' / decision / DOOMED node - only its declared ranking tail (C7/C13).
#' @param scheme a scheme id (see [marker_schemes()]).
#' @return a list with `nodes` (name, stage, kind, label, x, y, order, hard,
#'   detail) and `edges` (from, to, kind) tibbles.
#' @keywords internal
#' @noRd
marker_decision_schema <- function(scheme = "snp") {
  info <- marker_scheme(scheme)
  crit <- marker_criteria(scheme = scheme)
  stages <- marker_decision_stages(scheme)
  buf <- new.env(parent = emptyenv())
  buf$rows <- list()
  buf$edges <- list()
  push <- function(df) buf$rows[[length(buf$rows) + 1L]] <- df
  add_edge <- function(df) buf$edges[[length(buf$edges) + 1L]] <- df
  o <- 0
  refine_ids <- character(0)
  doomed_sources <- character(0)
  doomed_owners <- character(0)
  stage_sink_ids <- character(0)
  stage_sink_stage <- character(0)
  stage_sink_order <- numeric(0)
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
    if (length(gates)) {
      push(node(ids$gate, st$id, "gate", scheme_label("short", gates), o + seq_along(gates),
        hard = hard, detail = scheme_label("limit", gates)
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
      # every stage that can kill a candidate gets its own DOOMED sink, so the
      # graph says *which* stage killed it; the shared rail below the stages only
      # groups them for display (the `doomed` kind link from each stage).
      stage_doomed <- paste0(st$id, ":doomed")
      if (soft) {
        # one node per rung, in the declared order (least grave first) and each
        # labelled with the gate it protects and that gate's floor: the ladder is
        # the *structure* the engine walks (rung 1, then 2, ...), so a single
        # "refine" node would hide the tier a candidate actually stops at.
        rung_ids <- paste0(st$id, ":relax:", relax)
        push(node(
          rung_ids, st$id, "relax",
          paste0(seq_along(relax), ". ", scheme_label("short", relax)),
          decision_order + seq_along(relax),
          x = 0.85,
          detail = paste0(scheme_label("short", relax), " [", scheme_label("floor", relax), "]")
        ))
        refine_ids <- c(refine_ids, rung_ids)
        merge_order <- decision_order + length(relax) + 1L
        if (length(rung_ids) > 1L) add_edge(chain(rung_ids[-length(rung_ids)], rung_ids[-1L], "relax"))
        # the rungs are tried in ladder order from the failing decision
        add_edge(chain(ids$decision, rung_ids[1L], "fail"))
        add_edge(chain(rung_ids[length(rung_ids)], ids$merge, "relax"))
        # a rung whose ladder is fully spent has no rung left to try: the last rung
        # dooms (a failure at an earlier rung falls through to the next)
        doomed_sources <- c(doomed_sources, rung_ids[length(rung_ids)])
        doomed_owners <- c(doomed_owners, stage_doomed)
        # a *hard* gate in a mixed stage dooms through the decision: a candidate
        # that fails one can never be relaxed away, so it must not be drawn as
        # merely "relaxable" (C4).  The hard gates of a stage share that path.
        if (any(hard)) {
          doomed_sources <- c(doomed_sources, ids$decision)
          doomed_owners <- c(doomed_owners, stage_doomed)
        }
      } else {
        doomed_sources <- c(doomed_sources, ids$decision)
        doomed_owners <- c(doomed_owners, stage_doomed)
        merge_order <- decision_order + 1L
      }
      add_edge(chain(ids$decision, ids$merge, "pass"))
      stage_sinks <- doomed_sources[doomed_owners == stage_doomed]
      if (length(stage_sinks)) {
        stage_sink_ids <- c(stage_sink_ids, stage_doomed)
        stage_sink_stage <- c(stage_sink_stage, st$id)
        # between the decision it answers and the survivors it feeds, so the sink
        # sits on its own rank and every edge still flows downhill
        stage_sink_order <- c(stage_sink_order, decision_order + 0.5)
      }
      # a rung whose ladder is fully spent has no rung left to try: the last rung
      # dooms (a failure at an earlier rung falls through to the next)
      if (length(rung_ids)) {
        doomed_sources <- c(doomed_sources, rung_ids[length(rung_ids)])
        doomed_owners <- c(doomed_owners, stage_doomed)
        stage_sink_order[stage_sink_stage == st$id] <- merge_order + 0.5
      }
    } else {
      # a documentation-only stage (`gates = character(0)`, e.g. ASQ): it has no
      # computable gate, so it emits no gate / decision / DOOMED node - only its
      # declared ranking tail, hung directly off the previous stage (C7/C13)
      add_edge(chain(entry, ids$merge, "flow"))
      merge_order <- o + 1L
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
    if (has_tiers) {
      add_edge(chain(ids$merge, ids$tier[1L], "pass"))
      if (length(tiers) > 1L) add_edge(chain(ids$tier[-length(tiers)], ids$tier[-1L], "pass"))
      add_edge(chain(ids$tier[length(tiers)], ids$best, "pass"))
    } else {
      add_edge(chain(ids$merge, ids$best, "pass"))
    }
    entry <- ids$best
  }
  # one DOOMED sink per killing stage (a hard gate blocked, or every soft floor
  # is spent for that stage), so the structure names the killer instead of one
  # undifferentiated rail
  for (k in seq_along(stage_sink_ids)) {
    push(node(stage_sink_ids[k], stage_sink_stage[k], "doomed", "DOOMED", stage_sink_order[k]))
  }
  if (length(stage_sink_ids)) {
    # the last rung of a spent ladder dooms (an earlier rung failure simply falls
    # through to the next rung); every hard gate / unrelaxable decision too
    add_edge(tibble::tibble(from = doomed_sources, to = doomed_owners, kind = "doomed"))
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
#'   `assay$evaluation` (`"ok"`/`"violated"`), and the assay's own `$scheme`
#'   picks the tree.  The gate stamping is per stage id, so every chemistry's own
#'   stage (CAPS, ARMS, Sanger, ...) is coloured, not only the shared `as`/`pair`
#'   ones (C13).
#' @param scheme a scheme id (see [marker_schemes()]); `NULL` (the default) uses
#'   `assay$scheme`, falling back to `"snp"`.
#' @return a `tidygraph::tbl_graph`; nodes carry `name`, `stage`, `kind`,
#'   `label`, `x`, `y`, `hard`, `status`.
#' @examples
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' marker_decision_graph(assay)
#' }
#' @family KASP checks
#' @export
marker_decision_graph <- function(assay = NULL, scheme = NULL) {
  rlang::check_installed("tidygraph", reason = gettext("for the decision graph"))
  # the assay carries the chemistry it was assembled as, so the runtime graph is
  # *that* scheme's tree, never another's (C12)
  scheme <- scheme %||% (if (!is.null(assay)) assay$scheme else NULL) %||% "snp"
  s <- marker_decision_schema(scheme)
  nodes <- s$nodes
  if (!is.null(assay) && !is.null(assay$population)) {
    has_pop_node <- any(grepl(
      "frequency|haplo|ALT freq|hap covered", nodes$label,
      ignore.case = TRUE
    ))
    if (!has_pop_node) {
      msg <- paste0(
        "Population/haplotype data are present, but this decision schema has no ",
        "population/haplotype node; those gaps are not enforced by the DAG."
      )
      cli::cli_warn(gettext(msg))
    }
  }
  nodes$status <- "neutral"
  if (!is.null(assay) && !is.null(assay$selection)) {
    sel <- assay$selection
    status <- nodes$status
    gate <- nodes$kind == "gate"
    status[gate] <- "ok"
    # a scheme records its relaxed / violated / AS-stage gates by gate id; the node
    # name is that id prefixed with its stage, so build the candidate names for
    # *every* stage the schema shows (C13) - a discovered chemistry's own gate
    # nodes are coloured too, not only the shared as / pair ones.
    stage_ids <- unique(nodes$stage[gate])
    gate_nodes <- function(gs) {
      unlist(lapply(stage_ids, function(id) paste0(id, ":gate:", gs %||% character(0))), use.names = FALSE)
    }
    status[gate & nodes$name %in% gate_nodes(sel$relaxed)] <- "relaxed"
    status[gate & nodes$name %in% gate_nodes(sel$violated)] <- "violated"
    status[gate & nodes$name %in% gate_nodes(sel$as)] <- "violated"
    status[nodes$kind %in% c("tier", "decision", "merge", "outcome")] <- "ok"
    status[nodes$kind == "start"] <- "start"
    status[nodes$kind == "relax"] <- "relax"
    ev <- assay$evaluation
    if (is.null(ev) && !inherits(assay, "marker_assay")) {
      ev <- tryCatch(
        marker_design_checks(assay, scheme = scheme),
        error = function(e) NULL
      )
    }
    if (!is.null(ev)) {
      # Prefer the explicit status; fall back to the legacy logical for an
      # `$evaluation` written before statuses existed (NA was `unknown`).
      st <- if ("status" %in% names(ev)) {
        as.character(marker_check_status(ev$status))
      } else {
        ifelse(is.na(ev$pass), "unknown", ifelse(ev$pass, "pass", "fail"))
      }
      bad <- ev$criterion[st == "fail"]
      unknown <- ev$criterion[st == "unknown"]
      not_app <- ev$criterion[st == "not_applicable"]
      passed <- ev$criterion[st == "pass"]
      # unknown is its own colour (grey), never `ok`; a check that does not
      # apply is neutral; only an evaluated pass turns a node green.
      status[nodes$kind == "check" & nodes$label %in% not_app] <- "neutral"
      status[nodes$kind == "check" & nodes$label %in% bad] <- "violated"
      status[nodes$kind == "check" & nodes$label %in% unknown] <- "unknown"
      status[nodes$kind == "check" & nodes$label %in% passed] <- "ok"
    }
    nodes$status <- status
  }
  tidygraph::tbl_graph(nodes = nodes, edges = s$edges, directed = TRUE, node_key = "name")
}

# `tidygraph::as_tibble(g, "edges")` stores `from`/`to` as node INDICES.
# Every renderer needs names (DOT would otherwise invent nodes "1", "2", ...).
decision_edges <- function(g) {
  n <- tidygraph::as_tibble(g)
  e <- tidygraph::as_tibble(g, "edges")
  e$from <- n$name[as.integer(e$from)]
  e$to <- n$name[as.integer(e$to)]
  e
}

# One palette for the ggplot view and the Graphviz renderer.
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
      flow = "#78716C", pass = "#15803D", fail = "#B91C1C",
      relax = "#B45309", doomed = "#991B1B", check = "#78716C"
    )
  )
}

#' Draw the selection decision tree
#'
#' Renders the graph from [marker_decision_graph()] as a boxed, left-to-right
#' flow: every stage is its own cluster, every node a box drawn to the size of
#' its wrapped label, and edges are *routed* by a layout solver so they leave
#' a box, pass between boxes and land on the next one.  Positions come from
#' [decision_layout()]: a rank-pinned layered layout (one column per stage) when
#' `igraph` is installed (set `TIDYKASP_NO_SUGIYAMA=1` to force the manual grid).  Node **fill** is the node
#' kind (start/check/gate/decision/relax/merge/tier/outcome/doomed, a light
#' palette), while text and border carry the outcome (ok / violated / relaxed /
#' neutral).  Hard limits keep a **bold border**.  One marker chemistry is drawn
#' at a time ([marker_schemes()]); [marker_decision_text()] is the greppable
#' double.
#' @param x an assembled assay, a `marker_decision_graph()` result, or `NULL`.
#' @param scheme a scheme id (see [marker_schemes()]) for the bare schema; `NULL`
#'   (the default) uses the assay's own `$scheme` for an assay, else `"snp"`.

#' @return a `ggraph`/`ggplot` object.
#' @examples
#' \dontrun{
#' plot_decision_tree() # the schema
#' plot_decision_tree(scheme = "duplex") # the multiplex policy
#' plot_decision_tree(assay) # this assay's traversal
#' }
#' @family KASP checks
#' @export
plot_decision_tree <- function(x = NULL, scheme = NULL) {
  rlang::check_installed("ggraph", reason = gettext("to draw the decision tree"))
  scheme_id <- scheme %||%
    (if (inherits(x, "marker_assay")) x$scheme else NULL) %||%
    "snp"
  info <- marker_scheme(scheme_id)
  g <- if (is.null(x)) {
    marker_decision_graph(scheme = scheme_id)
  } else if (inherits(x, "marker_assay")) {
    # the assay's stored graph already *is* its scheme's tree (C12)
    x$decision %||% marker_decision_graph(x, scheme = scheme_id)
  } else {
    x
  }
  # positions come from the layout seam; drawing never places a node itself
  lay <- decision_layout(g)
  nodes <- attr(lay, "nodes")
  edges <- attr(lay, "edges")
  stage_levels <- attr(lay, "stage")
  nodes$stage <- factor(as.character(nodes$stage), levels = stage_levels)
  # fill = what the node *is* (light, modern, no blue); text/border = outcome
  pal <- decision_palette()
  kind_cols <- pal$kind
  status_cols <- pal$status
  edge_cols <- pal$edge
  # ---- boxes: one rect per node, at the size the solver reserved for it ------
  # This is what makes the picture readable: a solved layout knows where a box
  # is but not how big; a drawn box that disagrees with the reserved space is
  # how a boxed flow turns back into labels on a field.
  nodes$xmin <- nodes$x - nodes$w / 2
  nodes$xmax <- nodes$x + nodes$w / 2
  nodes$ymin <- nodes$y - nodes$h / 2
  nodes$ymax <- nodes$y + nodes$h / 2
  bbox <- decision_bbox(lay)
  # ---- edges: a routed polyline where the solver gave one, else a link -------
  # An edge whose route the solver did not report (an edge that stays inside a
  # stage, on the manual grid every edge) falls back to the straight link, but
  # it stops at the *box*, not at the anchor point.
  stop_at_box <- function(p, q2, which = c("head", "tail")) {
    which <- match.arg(which)
    dx <- q2$x - p$x
    dy <- q2$y - p$y
    if (!is.finite(dx) || !is.finite(dy) || (dx == 0 && dy == 0)) {
      return(list(x = p$x, y = p$y))
    }
    hw <- p$w / 2
    hh <- p$h / 2
    t <- min(if (dx != 0) hw / abs(dx) else Inf, if (dy != 0) hh / abs(dy) else Inf)
    list(x = p$x + dx * t, y = p$y + dy * t)
  }
  link <- tibble::tibble(
    from = edges$from, to = edges$to, kind = edges$kind,
    x = nodes$x[match(edges$from, nodes$name)],
    y = nodes$y[match(edges$from, nodes$name)],
    xend = nodes$x[match(edges$to, nodes$name)],
    yend = nodes$y[match(edges$to, nodes$name)]
  )
  box_of <- function(nm) {
    k <- match(nm, nodes$name)
    list(x = nodes$x[k], y = nodes$y[k], w = nodes$w[k], h = nodes$h[k])
  }
  for (i in seq_len(nrow(link))) {
    tail <- stop_at_box(box_of(link$from[i]), list(x = link$xend[i], y = link$yend[i]))
    head <- stop_at_box(box_of(link$to[i]), list(x = link$x[i], y = link$y[i]))
    link$x[i] <- tail$x
    link$y[i] <- tail$y
    link$xend[i] <- head$x
    link$yend[i] <- head$y
  }
  # ---- routed polylines, when the solver supplied one ------------------------
  # The route is a path, not a segment: it leaves a box's side, turns in the gap
  # between two stages and ends on a box edge, so the arrow head is never
  # buried under the box it points at.
  routes <- attr(lay, "routes")
  routed <- list()
  for (r in routes) {
    kind <- link$kind[match(paste(r$from, r$to), paste(link$from, link$to))]
    pid <- length(routed) + 1L
    # `interaction(group = …)` keeps one polyline per edge: ggplot groups paths
    # by the discrete aesthetics alone, which merged every same-kind edge into
    # a single zigzag
    pth <- tibble::tibble(
      from = r$from, to = r$to,
      kind = kind %||% "flow",
      x = r$x, y = r$y, edge = pid
    )
    routed[[length(routed) + 1L]] <- pth
  }
  # no routes (the manual grid) still draws an edge: a straight link, trimmed
  # to the boxes it joins
  if (!length(routed)) {
    routed <- tibble::tibble(
      from = c(rbind(link$from, link$from)), to = c(rbind(link$to, link$to)),
      kind = c(rbind(link$kind, link$kind)),
      x = c(rbind(link$x, link$xend)), y = c(rbind(link$y, link$yend)),
      edge = rep(seq_len(nrow(link)), each = 2L)
    )
  } else {
    routed <- dplyr::bind_rows(routed)
  }
  bands <- nodes |>
    dplyr::group_by(.data$stage, .drop = FALSE) |>
    dplyr::summarise(
      xmin = min(.data$xmin) - 0.12, xmax = max(.data$xmax) + 0.12,
      ymin = min(.data$ymin) - 0.12, ymax = max(.data$ymax) + 0.12, .groups = "drop"
    )
  bands <- bands[is.finite(bands$xmin), , drop = FALSE]
  bands$fill <- rep(c("#FAFAF9", "#F5F5F4"), length.out = nrow(bands))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = bands, inherit.aes = FALSE,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      fill = bands$fill, colour = "#E7E5E4", linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = bands, inherit.aes = FALSE,
      ggplot2::aes(x = .data$xmin + 0.08, y = .data$ymax + 0.18, label = .data$stage),
      hjust = 0, vjust = 0, size = 3.0, colour = "#44403C", fontface = "bold"
    ) +
    ggplot2::geom_path(
      data = routed,
      ggplot2::aes(
        x = .data$x, y = .data$y,
        colour = .data$kind, linetype = .data$kind,
        group = .data$edge
      ),
      linewidth = 0.35, alpha = 0.8,
      arrow = ggplot2::arrow(length = ggplot2::unit(1.5, "mm"), type = "closed"),
      linejoin = "mitre", lineend = "butt"
    ) +
    ggplot2::geom_rect(
      data = nodes, inherit.aes = FALSE,
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax,
        fill = .data$kind, colour = .data$status, linewidth = .data$hard
      )
    ) +
    ggplot2::geom_text(
      data = nodes, inherit.aes = FALSE,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label_wrapped),
      size = 9 / ggplot2::.pt, colour = "#111827", lineheight = 0.95
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        flow = "solid", pass = "solid", fail = "dotted",
        relax = "dashed", doomed = "dotted", check = "solid"
      ),
      limits = names(edge_cols), drop = FALSE, name = "edge"
    ) +
    ggplot2::scale_fill_manual(values = kind_cols, name = "node") +
    ggplot2::scale_colour_manual(
      values = c(status_cols, edge_cols), breaks = names(status_cols),
      limits = c(names(status_cols), names(edge_cols)), drop = FALSE, name = "status / edge"
    ) +
    ggplot2::scale_linewidth_manual(values = c(`TRUE` = 1.2, `FALSE` = 0.25), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(bbox$xmin, bbox$xmax), expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(limits = c(bbox$ymin, bbox$ymax), expand = ggplot2::expansion(0)) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(
      title = paste0(info$label, " \u2014 decision graph"),
      subtitle = paste0(info$discrimination, " | ", info$readout),
      caption = "stages run left to right; node fill = kind, border/text = outcome; hard limit = bold border"
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "#44403C", hjust = 0),
      plot.caption = ggplot2::element_text(size = 7.5, colour = "#78716C", hjust = 0),
      legend.position = "bottom", legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7.5)
    )
}
# nolint end

#' The decision tree of a scheme as text
#'
#' A plain-text, greppable export of the same registry the graph is drawn from:
#' the scheme's facts, the target-site checks, each selection stage's gates (with
#' hardness and floors), ladder, tiers and edges, the assay acceptance checks, and
#' the verdict.
#' @param scheme a scheme id (see [marker_schemes()]).
#' @return a character vector, one element per line.
#' @examples
#' cat(marker_decision_text(), sep = "\n")
#' cat(marker_decision_text("duplex"), sep = "\n")
#' @family KASP checks
#' @export
marker_decision_text <- function(scheme = "snp") {
  info <- marker_scheme(scheme)
  crit <- marker_criteria(scheme = scheme)
  s <- marker_decision_schema(scheme)
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
      "gates (hard = never relaxed), in application order (grave / cheap first):",
      check_lines(gate, 22L, soft = "- soft")
    )
    if (length(st$relax)) {
      res <- c(res, "relaxation ladder (least grave first):", purrr::map_chr(
        seq_along(st$relax),
        function(i) {
          gl <- st$relax[i]
          pad <- stringr::str_pad(scheme_label("short", gl), 22L, side = "right")
          glue::glue("  {i}. {pad} floor: {scheme_label('floor', gl)}")
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
    glue::glue("{info$label} - decision graph"), "",
    glue::glue(
      "discrimination: {info$discrimination} | readout: {info$readout}",
      " | zygosity: {info$zygosity} | loci: {info$loci} | amplicon: {info$amplicon}"
    ), ""
  )
  if (any(crit$stage == "site")) {
    out <- c(out, section("site", "target site checks"), "")
  }
  out <- c(out, purrr::map(marker_decision_stages(scheme), stage_block) |> purrr::list_c())
  out <- c(
    out, section("assay", "assay acceptance checks"), "",
    "DOOMED: a hard gate blocked, or every soft floor is exhausted.",
    glue::glue("verdict: all checks pass -> {info$verdict}; else see the failing node")
  )
  as.character(out)
}
