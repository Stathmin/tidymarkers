# =============================================================================
#  marker_priors.R  --  where every threshold comes from, and a thermodynamic check
# =============================================================================
#  Thresholds used to be folk numbers scattered across the code.  This file is
#  their provenance: one table naming each parameter, its value and the source
#  it follows (LGC/KASP practice, EasyKASP 2025, primerblast-oss), so a tuning
#  decision is auditable.  It also holds the independent thermodynamic check
#  that keeps those numbers honest: our Tm engine is the TmCalculator
#  nearest-neighbour family (SantaLucia tables).  The default salt correction is
#  the design coordinate the windows are stated on (Schildkraut 1965); the
#  NN-consistent reference model is SantaLucia 1998-2, the one Primer3/oligotm
#  applies, and it is what every external Tm comparison is stated in (ZW3.4).
# =============================================================================

#' Provenance of every KASP threshold
#'
#' A tidy table of the parameters in `KASP_LIMITS` (with the design window and
#' the thermodynamic model) mapped to the source each value follows.  It is the
#' audit trail for [marker_design_checks()]: a number with no source is a bug.
#' With `outcomes`, it answers a different question instead: how often each
#' chemistry worked or failed in the laboratory.  The counts are ordinal
#' evidence, not a success probability, and `calibrated` stays `FALSE`; the
#' solver only uses them to try a working chemistry earlier, never to change a
#' verdict.  A `result` the writer refuses is an error here rather than a row
#' silently left out of the tally.
#' @param outcomes a table from [marker_outcomes()] or
#'   [marker_outcome_corpus()], or `NULL` for the threshold provenance.
#' @return `NULL`: tibble with `parameter`, `value` and `source`.  Outcomes:
#'   one row per scheme with `works`, `fails`, `unknown`, `source` and
#'   `calibrated`.
#' @examples
#' marker_priors()
#' @family checks
#' @export
marker_priors <- function(outcomes = NULL) {
  if (!is.null(outcomes)) {
    return(prior_counts(outcomes))
  }
  structure(tibble::tribble(
    ~parameter, ~value, ~source,
    "as_tm (C)", paste(KASP_LIMITS$as_tm, collapse = "-"), "EasyKASP 2025; LGC KASP manual",
    "com_tm (C)", paste(KASP_LIMITS$com_tm, collapse = "-"), "LGC (COM Tm >= AS Tm, ~60 C)",
    "gc", paste(KASP_LIMITS$gc, collapse = "-"), "EasyKASP 2025 (~45% at the variant)",
    "amplicon (bp)", paste(KASP_LIMITS$amplicon, collapse = "-"), "EasyKASP 2025 (< 100 bp)",
    "max_run", as.character(KASP_LIMITS$max_run), "common primer-design practice (<= 4)",
    "hairpin_dg", as.character(KASP_LIMITS$hairpin_dg), "KASP practice; primerblast-oss 3'-dG gate",
    "duplex_dg", as.character(KASP_LIMITS$duplex_dg), "primerblast-oss dimer warn (dG -8)",
    "homodimer_dg", as.character(KASP_LIMITS$homodimer_dg), "primerblast-oss dimer warn (dG -8)",
    "discrimination", as.character(KASP_LIMITS$discrimination), "DECIPHER Taq-aware efficiency model",
    "discrimination margin", as.character(KASP_LIMITS$discrimination_margin), "local: min_disc near-miss band",
    "footprint_dropout", as.character(KASP_LIMITS$footprint_dropout), "local policy (frequency-weighted)",
    "clean_3p (nt)", "5", "primerblast-oss 3' window",
    "as_tm target (C)", as.character(AS_TM_TARGET), "LGC/EasyKASP (AS below the COM)",
    "com_tm target (C)", as.character(COM_TM_TARGET), "LGC/EasyKASP",
    "Tm model", NN_TABLE, "TmCalculator (SantaLucia nearest-neighbour tables)",
    "Tm salt correction", TM_SALT_METHOD,
    "design coordinate (windows on it); reference = SantaLucia 1998-2 (Primer3/oligotm)",
    "dG model", "OligoArrayAux (hybrid-min / hybrid-ss-min)", "DECIPHER's engine"
  ), class = c("priors", "tbl_df", "tbl", "data.frame"))
}

#' Count laboratory outcomes per chemistry
#'
#' An outcome row names a target and a scheme; the count is a plain tally, so a
#' scheme with no recorded run is absent rather than assumed.  Nothing is
#' pooled across artefact ids: the same `source` (the artefact ids that
#' contributed) is carried so a caller can see where the count came from.
#' @keywords internal
#' @noRd
prior_counts <- function(outcomes) {
  if (!is.data.frame(outcomes) || !all(c("result", "scheme") %in% names(outcomes))) {
    abort_msg(
      "`outcomes` must be a table from marker_outcomes() with `scheme` and `result`."
    )
  }
  o <- outcomes[!is.na(outcomes$scheme) & nzchar(outcomes$scheme), , drop = FALSE]
  if (!nrow(o)) {
    return(structure(tibble::tibble(
      scheme = character(0), works = integer(0), fails = integer(0),
      unknown = integer(0), source = character(0), calibrated = logical(0)
    ), class = c("priors", "tbl_df", "tbl", "data.frame")))
  }
  schemes <- sort(unique(o$scheme))
  bad <- o$result[is.na(o$result) | !(o$result %in% c("works", "fails", "unknown"))]
  if (length(bad)) {
    abort_msg(
      "`result` must be 'works', 'fails' or 'unknown'; got {paste(utils::head(unique(bad), 10L), collapse = ', ')}."
    )
  }
  structure(tibble::tibble(
    scheme = schemes,
    works = unname(vapply(schemes, function(s) sum(o$scheme == s & o$result == "works"), integer(1))),
    fails = unname(vapply(schemes, function(s) sum(o$scheme == s & o$result == "fails"), integer(1))),
    unknown = unname(vapply(schemes, function(s) sum(o$scheme == s & o$result == "unknown"), integer(1))),
    source = vapply(schemes, function(s) {
      ids <- if ("id" %in% names(o)) unique(o$id[o$scheme == s]) else character(0)
      ids <- ids[!is.na(ids)]
      if (!length(ids)) "outcomes" else paste(ids, collapse = ",")
    }, character(1)),
    calibrated = FALSE
  ), class = c("priors", "tbl_df", "tbl", "data.frame"))
}

#' Independent thermodynamic cross-check of a window's Tm across NN tables
#'
#' Our default Tm is one nearest-neighbour parameter set (SantaLucia 2004, via
#' TmCalculator).  The window built on it (`KASP_LIMITS$as_tm`) is only as firm as
#' the model: this runs the same oligos through several published DNA NN tables
#' and reports the spread.  A spread far below the window width means the
#' threshold is model-robust; a large spread means the window is model-dependent
#' and should be widened or the table pinned.
#' @param seq oligo(s) 5'->3'.
#' @param nn_tables DNA NN tables to compare (see [primer_tm()]); `NULL` uses
#'   SantaLucia (2004), Allawi (1997), Sugimoto (1996) and Breslauer (1986).
#' @param ... passed to [primer_tm()] (e.g. `DMSO`).
#' @return tidy tibble: `seq`, one column of Tm per table, `tm_min`, `tm_max` and
#'   `spread` (max - min) in degrees C.
#' @examples
#' thermo_crosscheck(c("GAAGGTGACCAAGTTCATGCT", "TTTTGGCGCGCCGATATATAT"))
#' @family thermodynamics
#' @export
thermo_crosscheck <- function(seq, nn_tables = NULL, ...) {
  if (is.null(nn_tables)) {
    nn_tables <- c(
      "DNA_NN_SantaLucia_2004", NN_TABLE_ALLAWI,
      "DNA_NN_Sugimoto_1996", "DNA_NN_Breslauer_1986"
    )
  }
  cols <- purrr::map(nn_tables, function(tb) primer_tm(seq, nn_table = tb, ...))
  out <- tibble::tibble(seq = as.character(seq), !!!rlang::set_names(cols, nn_tables))
  out$tm_min <- rlang::exec(pmin, !!!cols, na.rm = TRUE)
  out$tm_max <- rlang::exec(pmax, !!!cols, na.rm = TRUE)
  out$spread <- out$tm_max - out$tm_min
  out
}
