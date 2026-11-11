# =============================================================================
#  marker_dg_internal.R  --  pure-R duplex free energy (OligoArrayAux opt-out)
# =============================================================================
#  The second adapter at the `dg` seam.  `oligo_dg()` runs a UNAFold
#  subprocess; this one sums the SantaLucia 2004 nearest-neighbour table in R.
#
#  It models a perfect duplex and a single hairpin (stem + loop penalty).  A
#  stem-only nearest-neighbour screen cannot reproduce UNAFold's partition
#  function or coaxial stacking, so `fold` (structure prediction) stays
#  OligoArrayAux (see marker_engines()).  Reporting a structure there would be a
#  lie; [internal_hairpin_dg()] is a screen for the hairpin gate, nothing more.
#  It returns one dG per strand and never a base-pairing map.
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
#' [oligo_dg()]) returns `NA` with a warning: self-folding goes through
#' [internal_hairpin_dg()], and amplicon folding stays on OligoArrayAux, whose
#' ensemble is what the gates were calibrated on (see [marker_engines()]).
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
  if (anyNA(match(paste0(dinuc, "/", chartr("ACGT", "TGCA", dinuc)), rownames(SANTA_LUCIA_2004)))) {
    return(NA_real_)
  }
  terms <- internal_duplex_terms(s, na)
  terms[["dh"]] - (temp + 273.15) * terms[["ds"]] / 1000
}

#' dH and dS of a perfect duplex (kcal/mol, cal/K/mol) after the salt
#' term, i.e. exactly what [internal_duplex_dg_one()] integrates over `temp`.
#' @keywords internal
#' @noRd
internal_duplex_terms <- function(s, na) {
  L <- nchar(s)
  # `substring` is the vectorised slicer here; base `substr` truncates to the
  # first element on this build (see AGENTS.md).
  dinuc <- substring(s, 1:(L - 1L), 2:L)
  keys <- paste0(dinuc, "/", chartr("ACGT", "TGCA", dinuc))
  dh <- sum(SANTA_LUCIA_2004[keys, "dh"])
  ds <- sum(SANTA_LUCIA_2004[keys, "ds"])
  # initiation: `init` once, plus the terminal term for each end
  dh <- dh + SANTA_LUCIA_2004[["init", "dh"]]
  ds <- ds + SANTA_LUCIA_2004[["init", "ds"]]
  for (end in c(substring(s, 1L, 1L), substring(s, L, L))) {
    key <- if (end %in% c("A", "T")) "init_A/T" else "init_G/C"
    dh <- dh + SANTA_LUCIA_2004[[key, "dh"]]
    ds <- ds + SANTA_LUCIA_2004[[key, "ds"]]
  }
  ds <- ds + 0.368 * (L - 1L) * log(na)
  c(dh = dh, ds = ds)
}

# UNAFold's DNA hairpin loop parameters.  NOT from memory: every constant was
# recovered from OligoArrayAux output at two temperatures (37 C and -10 C) on
# single-stem fixtures, after subtracting the model's own stem terms, so the
# pair below is the loop contribution UNAFold actually adds
# (`tests/testthat/test-dg-internal.R` holds the fixtures and the gate).
# Loops of 4, 5 and 6 bases all fit the same constant (dH 3.1, dS +4.76), which
# is the known behaviour of short hairpin loops; the longer entries grow
# slightly and are interpolated.  The 12 is the polyline's anchor, not
# UNAFold's "longest specified loop" (that is 30): do not conflate the two.
HAIRPIN_LOOP_DH <- c(`3` = 3.1, `4` = 3.1, `5` = 3.1, `6` = 3.1, `8` = 3.3, `10` = 3.5, `12` = 3.7)
HAIRPIN_LOOP_DS <- c(`3` = 4.76, `4` = 4.76, `5` = 4.76, `6` = 4.76, `8` = 5.0, `10` = 5.2, `12` = 5.4)
HAIRPIN_LOOP_MIN <- 3L # shorter loops do not exist
HAIRPIN_LOOP_MAX <- 30L # UNAFold's own default `maxloop`
HAIRPIN_LOOP_ANCHOR <- max(as.integer(names(HAIRPIN_LOOP_DH)))
HAIRPIN_LOOP_LONG_DH <- 0.1 # per extra base past the anchor, up to the cap
HAIRPIN_LOOP_LONG_DS <- 0.1

#' dH and dS of a hairpin loop of `length` bases (see [internal_hairpin_dg()])
#' @keywords internal
#' @noRd
internal_hairpin_loop <- function(n) {
  n <- pmin(n, HAIRPIN_LOOP_MAX)
  dh <- stats::approx(HAIRPIN_LOOP_DH, xout = n, rule = 2)$y +
    HAIRPIN_LOOP_LONG_DH * pmax(0, n - HAIRPIN_LOOP_ANCHOR)
  ds <- stats::approx(HAIRPIN_LOOP_DS, xout = n, rule = 2)$y +
    HAIRPIN_LOOP_LONG_DS * pmax(0, n - HAIRPIN_LOOP_ANCHOR)
  cbind(dh = dh, ds = ds)
}

#' Every perfect stem in a sequence: arms `i .. i + k - 1` and
#' `j - k + 1 .. j` complementary, so base `i + d` pairs with `j - d` and the
#' loop holds `j - k - i + 1` bases (3..30 modelled).
#'
#' Mismatched or bulged stems are not modelled (see [internal_hairpin_dg()]).
#' Every stem width on the same arms is its own candidate, because a shorter
#' stem is a wider loop.
#' @return a `[dH, dS]` matrix (sum of stem stacks + one loop penalty), one row
#'   per candidate; zero rows when no stem exists.
#' @keywords internal
#' @noRd
hairpin_candidates <- function(s, na) {
  n <- nchar(s)
  k_min <- max(HAIRPIN_LOOP_MIN + 2L, 4L) # the gate screens >= 4 bp stems
  # a 12-mer can only hold a 4 bp stem, so the gate floor wins over "fits"
  k_max <- max((n - HAIRPIN_LOOP_MIN) %/% 2L, k_min)
  if (n < 2L * k_min) {
    return(matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("dh", "ds"))))
  }
  strand <- strsplit(s, "", fixed = TRUE)[[1L]]
  comp <- chartr("ACGT", "TGCA", strand)
  # `i` is the first base of the 5' arm, `j` the last base of the 3' arm; pair
  # `d` is `i + d` against `j - d`, so widening walks into the stem and stops
  # at the first mismatch.  The loop then has `j - k - i + 1` bases.
  i <- rep(seq_len(n), each = n)
  j <- rep(seq_len(n), times = n)
  ok <- j >= i + 2L * k_min - 1L
  if (!any(ok)) {
    return(matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("dh", "ds"))))
  }
  i <- i[ok]
  j <- j[ok]
  m <- integer(length(i))
  live <- seq_along(i)
  for (d in 0:(k_max - 1L)) {
    inside <- i[live] + d < j[live] - d &
      comp[i[live] + d] == strand[j[live] - d]
    m[live] <- d + 1L
    live <- live[inside]
    if (length(live) == 0L) {
      break
    }
  }
  keep <- m >= k_min
  i <- i[keep]
  j <- j[keep]
  m <- m[keep]
  # every stem length k_min..m from the same arms: a shorter stem is a wider
  # loop and a different candidate
  arm <- rep(seq_along(i), m - k_min + 1L)
  width <- sequence(m - k_min + 1L) + k_min - 1L
  # width `w` pairs `i .. i + w - 1` with `j - w + 1 .. j`, loop in between
  loop <- internal_hairpin_loop(j[arm] - width - i[arm] + 1L)
  if (length(arm) == 0L) {
    return(matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("dh", "ds"))))
  }
  stem <- t(vapply(
    substring(s, i[arm], i[arm] + width - 1L),
    function(a) internal_duplex_terms(a, na),
    numeric(2)
  ))
  stem + loop
}

#' Hairpin dG of a single strand without OligoArrayAux (pure R screen)
#'
#' The hairpin counterpart of [internal_duplex_dg()]: the most stable fold of
#' `seq1` alone, built from a perfect stem (reverse-complementary arms, scored
#' with the same SantaLucia (2004) nearest-neighbour table) plus one destabilising
#' loop penalty that grows with the loop length.  Same standard-state
#' convention, salt term and units as [oligo_dg()], so the calibrated hairpin
#' gate (dG >= -3 kcal/mol) reads the same number.
#'
#' **This is a screen, not UNAFold.**  It has no partition function, no
#' multi-state ensemble, no coaxial stacking, no interior or bulge loops and no
#' mismatched or bulged stems; the loop penalty is a polyline recovered from
#' OligoArrayAux output and interpolated, then capped at 30 bases.  It returns
#' one dG per sequence and **never a structure**: `fold`
#' ([fold_amplicon()]) stays on OligoArrayAux.  What it is calibrated for is the
#' hairpin gate and the ordering of primer candidates; a near-threshold number
#' should be confirmed with OligoArrayAux.
#'
#' Only stems of at least 4 bp are considered, so a sequence with no such stem
#' returns `0` (no stabilisation relative to the unfolded state) - the same
#' mapping [oligo_dg()] applies to UNAFold's "no structure" sentinel.
#'
#' @param seq1 oligo(s) 5'->3'; a vector is answered element-wise.
#' @param na sodium concentration (M); zero or negative gives `NA`.
#' @param temp temperature (C).
#' @return a numeric vector of dG (kcal/mol, more negative = a more stable
#'   fold), one per sequence.  `0` when no stem of at least 4 bp is stable, and
#'   `NA` for an ambiguous base or `na <= 0` - never a guessed number.
#' @examples
#' internal_hairpin_dg(c("GGGGGAAAACCCCC", "ACGTACGTACGT"))
#' @seealso [oligo_dg()], [internal_duplex_dg()], [hairpin_dg_thermo()],
#'   [marker_engines()]
#' @family thermodynamics
#' @export
internal_hairpin_dg <- function(seq1, na = 0.05, temp = 37) {
  seq1 <- as.character(seq1)
  n <- length(seq1)
  na <- rep_len(na, n)
  temp <- rep_len(temp, n)
  vapply(
    seq_len(n),
    function(i) internal_hairpin_dg_one(seq1[[i]], na[[i]], temp[[i]]),
    numeric(1)
  )
}

#' @keywords internal
#' @noRd
internal_hairpin_dg_one <- function(seq1, na, temp) {
  s <- toupper(gsub("[[:space:]]", "", seq1))
  bases <- strsplit(s, "", fixed = TRUE)[[1L]]
  if (length(bases) < 2L || anyNA(match(bases, c("A", "C", "G", "T")))) {
    return(NA_real_)
  }
  if (!is.numeric(na) || length(na) != 1L || !is.finite(na) || na <= 0) {
    return(NA_real_)
  }
  terms <- hairpin_candidates(s, na)
  if (nrow(terms) == 0L) {
    return(0)
  }
  # UNAFold drops a candidate whose dG is positive from the ensemble, so a
  # negative-minimum-only floor is the honest mapping (see oligo_dg())
  min(0, min(terms[, "dh"] - (temp + 273.15) * terms[, "ds"] / 1000))
}

#' Duplex dG through the configured engine
#'
#' One place the `dg` seam decides between OligoArrayAux and the pure-R
#' adapter, so call sites do not each re-read the option.  Hairpin (`seq2 =
#' NULL`) is not routed here: its adapter is [internal_hairpin_dg()].
#' @keywords internal
#' @noRd
thermo_duplex_dg <- function(seq1, seq2, na = 0.05, temp = 37) {
  if (identical(getOption("tidymarkers.thermo_engine", "auto"), "internal")) {
    internal_duplex_dg(seq1, seq2, na = na, temp = temp)
  } else {
    oligo_dg(seq1, seq2, na = na, temp = temp)
  }
}
