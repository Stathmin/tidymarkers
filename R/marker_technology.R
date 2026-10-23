# =============================================================================
#  marker_technology.R  --  compose a new technology from a base chemistry
# =============================================================================
#  The roadmap's Y1: a "technology" is a named composition of a base chemistry
#  with designer *parameters* and a specificity *policy*, assembled from blocks
#  rather than written from scratch.  Example: a "giga-KASP" = the SNP KASP with
#  longer AS primers (`params = list(as_lengths = 24:30)`) and a strict
#  off-target policy.
#
#  A technology is stored in a runtime registry that the scheme discovery layer
#  (`extra_scheme_defs()`) merges with the file-based `scheme_def_*()` specs, so a
#  composed technology is first-class everywhere a scheme is - `marker_schemes()`,
#  `marker_criteria()`, the decision DAG - with no core edit.  Its design is
#  dispatched to the base designer with the technology's parameters (see
#  `extra_scheme_designer()`).
# =============================================================================

# The runtime technology registry (scheme id -> spec).  Populated by
# marker_technology(); read by extra_scheme_defs().
.marker_tech_store <- new.env(parent = emptyenv())

#' The blocks a technology is composed from
#'
#' A catalogue of the composable axes `marker_technology()` honours: the base
#' chemistry (which brings its designer and its DAG stage), the specificity
#' policy (roadmap Y3) and the designer parameters.  It names the space so a
#' user can see what can be composed without reading the source.
#' @return a tibble with `block`, `kind`, `values`, `note`.
#' @examples
#' marker_filter_blocks()
#' @family KASP technologies
#' @export
marker_filter_blocks <- function() {
  tibble::tribble(
    ~block, ~kind, ~values, ~note,
    "base", "chemistry", "any scheme id (see marker_schemes())",
    "the designer + DAG stage the technology extends",
    "specificity", "policy", "unique | homeolog-aware | multicopy | off",
    "how a genome-level off-target is treated (Y3)",
    "params", "designer", "any named argument of the base designer",
    "e.g. as_lengths for longer primers, tm_range, gc_range"
  )
}

#' A scheme id -> its technology spec, or NULL (internal)
#' @keywords internal
#' @noRd
marker_technology_def <- function(name) {
  if (exists(name, envir = .marker_tech_store, inherits = FALSE)) {
    get(name, envir = .marker_tech_store)
  } else {
    NULL
  }
}

#' Compose a technology from a base chemistry
#'
#' Registers a new marker technology = a base scheme (its designer and DAG stage)
#' plus a specificity policy and designer parameters.  It is discovered like any
#' other scheme (see `marker_schemes()`), so `marker_criteria()`, the decision
#' tree and the designers pick it up automatically; `marker_design_at()` designs
#' it by calling the base designer with `params`.  A runtime registry, so it does
#' not survive a session restart (serialise `name`/`base`/`params`/`specificity`
#' and re-run `marker_technology()` to restore it - a declarative share format is
#' the next Y1 slice).
#' @param name the new scheme id (must not already exist unless `overwrite`).
#' @param base the chemistry to extend (see [marker_schemes()]).
#' @param params named list passed to the base designer (e.g.
#'   `list(as_lengths = 24:30)`); validated only for being a named list.
#' @param specificity optional specificity policy (see `KASP_SPECIFICITY_MODES`);
#'   defaults to the base scheme's policy.
#' @param label,verdict optional human-readable overrides.
#' @param overwrite replace an existing technology of the same `name`.
#' @return the technology spec, invisibly.
#' @examples
#' marker_technology("giga_snp", base = "snp",
#'                   params = list(as_lengths = 24:30), overwrite = TRUE)
#' marker_schemes()
#' @family KASP technologies
#' @export
marker_technology <- function(name, base = "snp", params = list(),
                              specificity = NULL, label = NULL, verdict = NULL,
                              overwrite = FALSE) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    cli::cli_abort(gettext("`name` must be a single non-empty string."))
  }
  if (!overwrite && name %in% marker_schemes()$scheme) {
    cli::cli_abort(gettext("Scheme {name} already exists; pass overwrite = TRUE to replace it."))
  }
  base_row <- marker_scheme(base) # aborts on an unknown base
  if (!is.list(params) || (length(params) && is.null(names(params)))) {
    cli::cli_abort(gettext("`params` must be a (possibly empty) named list."))
  }
  spec <- list(
    scheme = name,
    family = base_row$family,
    label = label %||% paste0(base_row$label, " [", name, "]"),
    discrimination = base_row$discrimination,
    readout = base_row$readout,
    zygosity = base_row$zygosity,
    amplicon = base_row$amplicon,
    loci = base_row$loci,
    verdict = verdict %||% base_row$verdict,
    gates = character(0),
    relax = character(0),
    tiers = character(0),
    base = base,
    params = params,
    specificity = specificity %||% marker_specificity(base)
  )
  if (!is.null(specificity) && !specificity %in% KASP_SPECIFICITY_MODES) {
    cli::cli_abort(gettext("Unknown specificity mode {specificity}."))
  }
  assign(name, spec, envir = .marker_tech_store)
  invisible(spec)
}

#' List the composed technologies of this session
#'
#' @return a tibble with `scheme`, `base`, `family`, `specificity`, `params`;
#'   zero rows when none were composed.
#' @examples
#' marker_technology("giga_snp", base = "snp",
#'                   params = list(as_lengths = 24:30), overwrite = TRUE)
#' marker_technologies()
#' @family KASP technologies
#' @export
marker_technologies <- function() {
  ids <- sort(ls(.marker_tech_store, all.names = FALSE))
  if (!length(ids)) {
    return(tibble::tibble(
      scheme = character(0), base = character(0), family = character(0),
      specificity = character(0), params = list()
    ))
  }
  purrr::map(ids, function(id) {
    d <- marker_technology_def(id)
    tibble::tibble(
      scheme = id, base = d$base, family = d$family,
      specificity = d$specificity, params = list(d$params)
    )
  }) |>
    purrr::list_rbind()
}
