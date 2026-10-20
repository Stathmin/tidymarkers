# =============================================================================
#  marker_scheme_registry.R  --  discovery layer for per-scheme modules
# =============================================================================
#  A marker chemistry lives in one file (R/marker_<scheme>.R): its designer, its
#  class, and a `scheme_def_<scheme>()` spec (the scheme row, its selection
#  stage's gates / tiers / relaxation, and the labels for those gates).  This
#  file *discovers* those specs, so adding a chemistry is adding one file - the
#  scheme registry, the criterion registry and the decision tree all pick it up
#  automatically, with no central edit and no drift.  The existing built-in
#  schemes keep their constants; the discovered ones are appended.
# =============================================================================

#' Every discovered per-scheme spec (one per `scheme_def_*()` function)
#' @return a named list of scheme specs (name = scheme id).
#' @keywords internal
#' @noRd
extra_scheme_defs <- function() {
  ns <- asNamespace("tidymarkers")
  nms <- sort(grep("^scheme_def_", ls(ns, all.names = TRUE), value = TRUE))
  file_defs <- stats::setNames(lapply(nms, function(nm) get(nm, envir = ns)()), sub("^scheme_def_", "", nms))
  # technologies composed at runtime (marker_technology()) are schemes too
  dyn <- as.list(.marker_tech_store)
  c(file_defs, dyn)
}

#' A gate-label lookup merged over the built-in vectors and every scheme spec
#' @param field one of `short`, `limit`, `floor`, `hard`, `label`.
#' @param ids gate ids.
#' @return a character vector, one per id.
#' @keywords internal
#' @noRd
scheme_label <- function(field, ids) {
  base <- switch(field,
    short = KASP_GATE_SHORT,
    limit = KASP_GATE_LIMIT,
    floor = KASP_GATE_FLOOR_LABEL,
    hard = KASP_GATE_HARD_LABEL,
    label = KASP_GATE_LABEL
  )
  key <- switch(field,
    short = "gate_short",
    limit = "gate_limit",
    floor = "gate_floor_label",
    hard = "gate_hard_label",
    label = "gate_label"
  )
  out <- unname(base[ids])
  for (d in extra_scheme_defs()) {
    m <- d[[key]]
    if (is.null(m)) next
    miss <- is.na(out) & ids %in% names(m)
    out[miss] <- unname(m[ids[miss]])
  }
  out
}

#' Discovered schemes as `marker_schemes()` rows
#' @keywords internal
#' @noRd
extra_scheme_rows <- function() {
  defs <- extra_scheme_defs()
  if (!length(defs)) {
    return(NULL)
  }
  purrr::map(defs, function(d) {
    tibble::tibble(
      scheme = d$scheme, family = d$family, label = d$label,
      discrimination = d$discrimination, readout = d$readout, zygosity = d$zygosity,
      amplicon = d$amplicon, loci = as.integer(d$loci %||% 1L), verdict = d$verdict
    )
  }) |>
    purrr::list_rbind()
}

#' Discovered schemes' selection stages as `marker_criteria()` rows
#' @keywords internal
#' @noRd
extra_stage_rows <- function() {
  defs <- extra_scheme_defs()
  # a composed technology may carry no own gate/tier (it just re-parameterises a
  # base) - it adds no stage rows; the base's stages already cover it
  defs <- defs[vapply(defs, function(d) length(d$gates) || length(d$tiers), logical(1))]
  if (!length(defs)) {
    return(NULL)
  }
  purrr::map(defs, function(d) {
    gates <- d$gates
    tiers <- d$tiers
    dplyr::bind_rows(
      tibble::tibble(
        stage = d$scheme, role = "gate", id = paste0(d$scheme, ":gate:", gates),
        label = unname(d$gate_short[gates]), level = "alignment",
        hard = !(gates %in% d$relax), limit = unname(d$gate_limit[gates]),
        order = seq_along(gates), schemes = rep(list(d$scheme), length(gates))
      ),
      tibble::tibble(
        stage = d$scheme, role = "tier", id = paste0(d$scheme, ":tier:", names(tiers)),
        label = unname(tiers), level = "alignment", hard = FALSE,
        limit = "maximise (per-tier bottleneck min)", order = seq_along(tiers),
        schemes = rep(list(d$scheme), length(tiers))
      )
    )
  }) |>
    purrr::list_rbind()
}

#' Discovered schemes as `marker_decision_stages()` entries (one for `scheme`)
#' @keywords internal
#' @noRd
extra_decision_stages <- function(scheme = NULL) {
  defs <- extra_scheme_defs()
  defs <- defs[vapply(defs, function(d) length(d$gates) || length(d$tiers), logical(1))]
  if (!is.null(scheme)) {
    defs <- defs[names(defs) == scheme]
  }
  purrr::map(defs, function(d) {
    list(id = d$scheme, label = d$stage_label, gates = d$gates, relax = d$relax, tiers = d$tiers)
  })
}

#' The designer for a discovered scheme (a `design_<scheme>_marker()` function)
#' @keywords internal
#' @noRd
extra_scheme_designer <- function(scheme) {
  ns <- asNamespace("tidymarkers")
  nm <- paste0("design_", scheme, "_marker")
  if (exists(nm, envir = ns, inherits = FALSE)) {
    return(get(nm, envir = ns))
  }
  # a composed technology has no designer of its own: dispatch to its base's
  # designer, passing the technology's parameters (and, for the KASP core, its
  # scheme id, so the assay is assembled as this technology)
  d <- marker_technology_def(scheme)
  if (is.null(d)) {
    return(NULL)
  }
  params <- d$params %||% list()
  if (d$base %in% c("snp", "indel", "multi")) {
    function(genome, snp, variants = NULL, BPPARAM = NULL, ...) {
      rlang::exec(assemble_kasp_assay, genome, snp,
        variants = variants, scheme = scheme, BPPARAM = BPPARAM, !!!params, ...
      )
    }
  } else {
    inner <- extra_scheme_designer(d$base)
    if (is.null(inner)) {
      return(NULL)
    }
    function(genome, snp, ...) rlang::exec(inner, genome, snp, !!!params, ...)
  }
}

#' Discovered scheme ids whose `family` matches (e.g. the extra KASP-family ones)
#' @keywords internal
#' @noRd
extra_family_schemes <- function(family) {
  defs <- extra_scheme_defs()
  names(defs)[vapply(defs, function(d) identical(d$family, family), logical(1))]
}
