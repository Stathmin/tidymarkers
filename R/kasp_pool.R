# =============================================================================
#  kasp_pool.R  --  many primers in one tube: a pairwise compatibility report
# =============================================================================
#  The duplex engine (R/kasp_duplex.R) checks two loci; a multiplex needs the
#  same check over EVERY primer of the pool, because the number of unwanted
#  interactions grows quadratically with primer count and is the usual cause of
#  multiplex failure.  This is the general object: give it any set of oligos and
#  it reports, per pair, the longest Watson-Crick run (a cross-dimer) and the
#  self-complementarity (a hairpin), plus the thermodynamic dimer dG37 when the
#  OligoArrayAux engine is available.  One vocabulary for "would these primers
#  misbehave together", shared by the duplex stage and any multiplex set.
# =============================================================================

#' Pairwise compatibility of a primer pool
#'
#' Takes any set of oligos (e.g. the primers of several markers) and returns one
#' tidy row per unordered pair: a `"self"` row is that oligo's hairpin, a
#' `"cross"` row an inter-primer interaction.  `run` is the longest complementary
#' stretch (see [complementary_run()]), `dg` the dimer dG37 (OligoArrayAux, `NA`
#' without it) and `pass` the verdict.  This is the general "many primers" check
#' behind the duplex stage ([duplex_compatibility()]).
#' @param primers a named character vector, or a (named) list of such vectors
#'   whose elements are pooled; names are the pool labels.
#' @param min_cross_run,max_selfcomp the base-pair rules (see
#'   [duplex_compatibility()]).
#' @return a tidy tibble: `a`, `b`, `kind` (`"self"`/`"cross"`), `run` (bp), `dg`
#'   (kcal/mol), `pass`.  All rows pass (`all(pass)`) for a clean pool.
#' @examples
#' primer_set_compatibility(c(
#'   AS1 = "GAAGGTGACCAAGTTCATGCT", AS2 = "GAAGGTCGGAGTCAACGGATT"
#' ))
#' @family KASP specificity
#' @export
primer_set_compatibility <- function(primers, min_cross_run = 6L, max_selfcomp = 5L) {
  primers <- unlist(primers, use.names = TRUE)
  if (is.null(names(primers)) || any(!nzchar(names(primers)))) {
    names(primers) <- paste0("p", seq_along(primers))
  }
  d <- tidyr::expand_grid(i = seq_along(primers), j = seq_along(primers)) |>
    dplyr::filter(.data$i <= .data$j) |>
    dplyr::mutate(
      a = names(primers)[.data$i], b = names(primers)[.data$j],
      kind = dplyr::if_else(.data$i == .data$j, "self", "cross"),
      run = purrr::map2_int(.data$i, .data$j, function(i, j) {
        complementary_run(primers[[i]], primers[[j]])
      })
    )
  dg <- tryCatch(dimer_thermo(primers), error = function(e) NULL)
  d$dg <- if (is.null(dg) || !nrow(dg)) {
    NA_real_
  } else {
    purrr::map2_dbl(d$a, d$b, function(x, y) {
      hit <- dg$dg[(dg$a == x & dg$b == y) | (dg$a == y & dg$b == x)]
      if (length(hit)) hit[[1L]] else NA_real_
    })
  }
  d$pass <- dplyr::if_else(d$kind == "self", d$run <= max_selfcomp, d$run < min_cross_run)
  dplyr::select(d, "a", "b", "kind", "run", "dg", "pass")
}
