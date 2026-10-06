# =============================================================================
#  kasp_priors.R  --  where every threshold comes from, and a thermodynamic check
# =============================================================================
#  Thresholds used to be folk numbers scattered across the code.  This file is
#  their provenance: one table naming each parameter, its value and the source
#  it follows (LGC/KASP practice, EasyKASP 2025, primerblast-oss), so a tuning
#  decision is auditable.  It also holds the independent thermodynamic check
#  that keeps those numbers honest: our Tm engine is the Primer3/oligotm
#  nearest-neighbour family (TmCalculator, SantaLucia tables); running the same
#  oligos through several published NN parameter sets shows how much the Tm -
#  and therefore the window built on it - actually wobbles.
# =============================================================================

#' Provenance of every KASP threshold
#'
#' A tidy table of the parameters in `KASP_LIMITS` (with the design window and
#' the thermodynamic model) mapped to the source each value follows.  It is the
#' audit trail for [kasp_design_checks()]: a number with no source is a bug.
#' @return tibble with `parameter`, `value` and `source`.
#' @examples
#' kasp_priors()
#' @family KASP checks
#' @export
kasp_priors <- function() {
  tibble::tribble(
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
    "footprint_dropout", as.character(KASP_LIMITS$footprint_dropout), "local policy (frequency-weighted)",
    "clean_3p (nt)", "5", "primerblast-oss 3' window",
    "as_tm target (C)", as.character(AS_TM_TARGET), "LGC/EasyKASP (AS below the COM)",
    "com_tm target (C)", as.character(COM_TM_TARGET), "LGC/EasyKASP",
    "Tm model", NN_TABLE, "TmCalculator (Primer3/oligotm NN family, SantaLucia)",
    "dG model", "OligoArrayAux (hybrid-min / hybrid-ss-min)", "DECIPHER's engine"
  )
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
#' @param nn_tables DNA NN tables to compare (see [primer_tm()]).
#' @param ... passed to [primer_tm()] (e.g. `DMSO`).
#' @return tidy tibble: `seq`, one column of Tm per table, `tm_min`, `tm_max` and
#'   `spread` (max - min) in degrees C.
#' @examples
#' thermo_crosscheck(c("GAAGGTGACCAAGTTCATGCT", "TTTTGGCGCGCCGATATATAT"))
#' @family KASP thermodynamics
#' @export
thermo_crosscheck <- function(seq,
                              nn_tables = c(
                                "DNA_NN_SantaLucia_2004", "DNA_NN_Allawi_1998",
                                "DNA_NN_Sugimoto_1996", "DNA_NN_Breslauer_1986"
                              ),
                              ...) {
  out <- tibble::tibble(seq = as.character(seq))
  for (tb in nn_tables) out[[tb]] <- primer_tm(seq, nn_table = tb, ...)
  cols <- unname(out[nn_tables])
  out$tm_min <- do.call(pmin, c(cols, list(na.rm = TRUE)))
  out$tm_max <- do.call(pmax, c(cols, list(na.rm = TRUE)))
  out$spread <- out$tm_max - out$tm_min
  out
}
