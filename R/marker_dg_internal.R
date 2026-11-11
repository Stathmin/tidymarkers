# =============================================================================
#  marker_dg_internal.R  --  pure-R duplex free energy (OligoArrayAux opt-out)
# =============================================================================
#  The second adapter at the `dg` seam.  `oligo_dg()` runs a UNAFold
#  subprocess; this one sums the SantaLucia 2004 nearest-neighbour table in R.
#
#  It is DUPLEX ONLY on purpose.  A stem-only nearest-neighbour screen cannot
#  reproduce UNAFold's partition function or coaxial stacking, and the hairpin
#  gates were calibrated on UNAFold, so `hairpin_dg` / `fold` stay OligoArrayAux
#  (see marker_engines()).  Reporting a number there would be a lie; returning
#  NA lets the gate stay "unknown".
#
#  The parameters mirror TmCalculator::tm_nn (extracted from the installed
#  namespace, not memory): stacks + `init` once + `init_A/T` per A/T terminal,
#  salt as 0.368 * (L - 1) * ln(mon).  The result is a STANDARD-STATE dG (1 M
#  strands), which is what `hybrid-min -q` prints and what the gate thresholds
#  were calibrated on -- the `R * ln(CT)` term of a Tm is deliberately not
#  added, or every duplex would read ~11 kcal/mol weaker than OligoArrayAux.
#
#  Salts: `na` is the monovalent concentration in M, the same units as
#  `oligo_dg(na = )`.  Mg and dNTPs are not modelled here because `oligo_dg()`
#  does not take them either (it passes only `-N`).
# =============================================================================

# SantaLucia (2004) unified DNA nearest-neighbour parameters: dH in kcal/mol,
# dS in cal/K/mol.  Extracted from TmCalculator:::.TM_CONSTANTS.
SANTA_LUCIA_2004 <- local({
  m <- matrix(
    c(
      0.2, -5.7, 2.2, 6.9, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, -1.4,
      -7.6, -21.3, -7.2, -20.4, -7.2, -21.3, -8.5, -22.7, -8.4, -22.4,
      -7.8, -21.0, -8.2, -22.2, -10.6, -27.2, -9.8, -24.4, -8.0, -19.9,
      -7.6, -21.3, -8.4, -22.4, -7.8, -21.0, -8.2, -22.2, -8.5, -22.7,
      -8.0, -19.9
    ),
    ncol = 2, byrow = TRUE
  )
  rownames(m) <- c(
    "init", "init_A/T", "init_G/C", "init_oneG/C", "init_allA/T",
    "init_5T/A", "sym", "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA",
    "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC", "TT/AA", "AC/TG",
    "AG/TC", "TC/AG", "TG/AC", "CC/GG"
  )
  colnames(m) <- c("dh", "ds")
  m
})

#' Nearest-neighbour duplex dG without OligoArrayAux (pure R)
#'
#' The opt-out adapter for [oligo_dg()] modes that do not need UNAFold.  It sums
#' the SantaLucia (2004) nearest-neighbour table and needs no subprocess, so it
#' works on a machine where OligoArrayAux cannot be installed.  It reports a
#' standard-state dG (1 M strands), the same convention as `hybrid-min -q`, so
#' the calibrated duplex gates read the same number.  It models a perfect
#' duplex only: `seq2` must be the reverse complement of `seq1`.
#'
#' **Duplex only.**  `seq2 = NULL` (the hairpin / folding mode of
#' [oligo_dg()]) is not modelled: it returns `NA` with a warning.  Hairpin and
#' amplicon-folding numbers stay on OligoArrayAux, whose ensemble is what the
#' gates were calibrated on (see [marker_engines()]).
#'
#' The model is monovalent-salt only (`na`, M).  It does not add Mg or dNTPs,
#' for the same reason [oligo_dg()] does not: neither takes them.
#'
#' @param seq1,seq2 oligos 5'->3'; `seq2` must be the reverse complement of
#'   `seq1`.  Both are recycled over their common length.
#' @param na sodium concentration (M); zero or negative gives `NA`.
#' @param temp temperature (C).
#' @return a numeric vector of dG (kcal/mol, more negative = more stable).
#'   `NA` for a non-duplex call, an ambiguous base, `na <= 0`, or a `seq2` that
#'   is not the reverse complement of `seq1` - never a guessed number.
#' @examples
#' internal_duplex_dg("GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC")
#' @seealso [oligo_dg()], [marker_engines()], [primer_tm()]
#' @family thermodynamics
#' @export
internal_duplex_dg <- function(seq1, seq2, na = 0.05, temp = 37) {
  seq1 <- as.character(seq1)
  if (missing(seq2) || is.null(seq2)) {
    cli::cli_warn(gettext(
      "The internal engine has no hairpin/fold model; use OligoArrayAux (`oligo_dg(seq1)`)."
    ))
    return(rep(NA_real_, length(seq1)))
  }
  seq2 <- as.character(seq2)
  n <- max(length(seq1), length(seq2))
  if (length(seq1) != length(seq2)) {
    if (n %% length(seq1) != 0L || n %% length(seq2) != 0L) {
      cli::cli_abort(gettext("`seq1` and `seq2` must be the same length (or recycle)."))
    }
  }
  seq1 <- rep_len(seq1, n)
  seq2 <- rep_len(seq2, n)
  na <- rep_len(na, n)
  temp <- rep_len(temp, n)
  vapply(
    seq_len(n),
    function(i) internal_duplex_dg_one(seq1[[i]], seq2[[i]], na[[i]], temp[[i]]),
    numeric(1)
  )
}

#' @keywords internal
#' @noRd
internal_duplex_dg_one <- function(seq1, seq2, na, temp) {
  s <- toupper(gsub("[[:space:]]", "", seq1))
  partner <- toupper(gsub("[[:space:]]", "", seq2))
  bases <- strsplit(s, "", fixed = TRUE)[[1L]]
  if (length(bases) < 2L || anyNA(match(bases, c("A", "C", "G", "T")))) {
    return(NA_real_)
  }
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
  if (!identical(partner, rc)) {
    return(NA_real_)
  }
  if (!is.numeric(na) || length(na) != 1L || !is.finite(na) || na <= 0) {
    return(NA_real_)
  }
  L <- nchar(s)
  # `substring` is the vectorised slicer here; base `substr` truncates to the
  # first element on this build (see AGENTS.md).
  dinuc <- substring(s, 1:(L - 1L), 2:L)
  keys <- paste0(dinuc, "/", chartr("ACGT", "TGCA", dinuc))
  if (anyNA(match(keys, rownames(SANTA_LUCIA_2004)))) {
    return(NA_real_)
  }
  dh <- sum(SANTA_LUCIA_2004[keys, "dh"])
  ds <- sum(SANTA_LUCIA_2004[keys, "ds"])
  # initiation: `init` once, plus the terminal term for each end
  dh <- dh + SANTA_LUCIA_2004[["init", "dh"]]
  ds <- ds + SANTA_LUCIA_2004[["init", "ds"]]
  for (end in c(bases[[1L]], bases[[L]])) {
    key <- if (end %in% c("A", "T")) "init_A/T" else "init_G/C"
    dh <- dh + SANTA_LUCIA_2004[[key, "dh"]]
    ds <- ds + SANTA_LUCIA_2004[[key, "ds"]]
  }
  ds <- ds + 0.368 * (L - 1L) * log(na)
  dh - (temp + 273.15) * ds / 1000
}

#' Duplex dG through the configured engine
#'
#' One place the `dg` seam decides between OligoArrayAux and the pure-R
#' adapter, so call sites do not each re-read the option.  `seq2` is required:
#' hairpin (`seq2 = NULL`) is always [oligo_dg()], because the internal adapter
#' has no folding model.
#' @keywords internal
#' @noRd
thermo_duplex_dg <- function(seq1, seq2, na = 0.05, temp = 37) {
  if (identical(getOption("tidymarkers.thermo_engine", "auto"), "internal")) {
    internal_duplex_dg(seq1, seq2, na = na, temp = temp)
  } else {
    oligo_dg(seq1, seq2, na = na, temp = temp)
  }
}
