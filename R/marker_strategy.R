# =============================================================================
#  marker_strategy.R  --  the rescue policy as data, and a batch designer
# =============================================================================
#  Which chemistry to try for a locus is a *policy*, not a hard-coded vector:
#  `rescue_strategy()` is that policy as a tidy object (one row per step, in the
#  order to try them), `design_marker()` walks it, and `design_markers()` applies
#  one chemistry (or the whole ladder) over a table of loci.  The policy can be
#  inspected, filtered, reordered or restricted like any other tibble.
# =============================================================================

#' Scheme ids a single locus can actually be taken by
#'
#' The default rescue ladder walks one locus at a time, so it can only contain
#' schemes that have a per-locus designer.  That excludes a duplex chemistry
#' (`loci = 2`, it needs two loci) and any registered scheme without a
#' `design_<scheme>_marker()` engine (`asq`, the marker-agnostic `primer` row).
#' @return a character vector of scheme ids, in registry order.
#' @keywords internal
#' @noRd
per_locus_schemes <- function() {
  base <- c("snp", "indel", "multi", "caps", "arms", "hrm")
  disc <- names(extra_scheme_defs())
  disc <- disc[!vapply(disc, function(s) is.null(extra_scheme_designer(s)), logical(1))]
  intersect(marker_schemes()$scheme, c(base, disc))
}

#' The rescue policy: which chemistries to try, in order
#'
#' A marker is not one chemistry.  When KASP fails a locus another may take it
#' (see [marker_schemes()]), so the "try in turn" order is an object here rather
#' than a vector buried in [design_marker()].  The default ladder holds only the
#' schemes that can be tried *one locus at a time* (see [marker_schemes()]'s `loci`;
#' a duplex needs two).  Steps can be dropped or reordered:
#' `rescue_strategy(c("hrm", "caps"))`, or filtered by family.
#' @param schemes scheme ids in order (see [marker_schemes()]).  Default: KASP
#'   (`"snp"`) first, then every other per-locus scheme in registry order.
#' @return a one-row-per-step tibble: `order`, `tried`, then the [marker_schemes()]
#'   columns (`scheme`, `family`, `label`, `discrimination`, `readout`,
#'   `zygosity`, `amplicon`, `loci`, `verdict`).
#' @examples
#' rescue_strategy()
#' rescue_strategy(c("hrm", "caps", "snp"))
#' @family KASP design
#' @export
rescue_strategy <- function(schemes = NULL) {
  ladder <- per_locus_schemes()
  if (is.null(schemes)) {
    schemes <- c(intersect("snp", ladder), setdiff(ladder, "snp"))
  }
  all <- marker_schemes()$scheme
  i <- match(schemes, all)
  unknown <- schemes[is.na(i)]
  if (length(unknown)) {
    cli::cli_abort(glue::glue(
      gettext("Unknown marker scheme(s): {paste(unknown, collapse = ', ')}; see marker_schemes().")
    ))
  }
  out <- marker_schemes()[i, , drop = FALSE]
  out$order <- seq_along(schemes)
  out$tried <- TRUE
  dplyr::relocate(out, "order", "tried")
}

#' Design markers for a table of loci with one chemistry (or the rescue ladder)
#'
#' The batch twin of [design_marker()]: one row per locus, saying which chemistry
#' took it and whether it designed.  With `scheme` given it is one chemistry for
#' every locus; with `scheme = NULL` it walks [rescue_strategy()], so a locus KASP
#' dooms falls through to CAPS / ARMS / HRM / the discovered chemistries.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param loci tibble with `chrom`, `pos`, `ref`, `alt` (an optional `alt2` is
#'   honoured by the KASP scheme).
#' @param scheme a scheme id (see [marker_schemes()]), or `NULL` for the full rescue
#'   ladder.
#' @param ... passed to the per-scheme designer (e.g. `flank`).
#' @return a one-row-per-locus tibble: `chrom`, `pos`, `ref`, `alt`, `scheme` (the
#'   winner, or `NA`) and `ok`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_markers(genome, marker_example()[1:3, ], scheme = "hrm")
#' design_markers(genome, marker_example()[1:3, ], scheme = NULL) # the ladder
#' }
#' @family KASP design
#' @export
design_markers <- function(genome, loci, scheme = "snp", ...) {
  need <- c("chrom", "pos", "ref", "alt")
  missing <- setdiff(need, names(loci))
  if (length(missing)) {
    cli::cli_abort(glue::glue(
      gettext("`loci` is missing column(s): {paste(missing, collapse = ', ')}")
    ))
  }
  purrr::map(seq_len(nrow(loci)), function(i) {
    r <- loci[i, , drop = FALSE]
    snp <- tibble::tibble(
      chrom = as.character(r$chrom), pos = as.integer(r$pos),
      ref = as.character(r$ref), alt = as.character(r$alt),
      alt2 = if ("alt2" %in% names(r)) as.character(r$alt2) else NA_character_
    )
    if (is.null(scheme)) {
      d <- design_marker(genome, snp, ...)
      chosen <- d$scheme
      ok <- is_target_designable(d$marker)
    } else {
      m <- tryCatch(suppressWarnings(marker_design_at(genome, snp, scheme, ...)), error = function(e) NULL)
      chosen <- scheme
      ok <- is_target_designable(m)
    }
    tibble::tibble(
      chrom = snp$chrom, pos = snp$pos, ref = snp$ref, alt = snp$alt,
      scheme = if (is.na(chosen)) NA_character_ else chosen, ok = ok
    )
  }) |>
    purrr::list_rbind()
}
