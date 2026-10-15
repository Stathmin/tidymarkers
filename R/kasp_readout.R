# =============================================================================
#  kasp_readout.R  --  what a marker's signal model does (and does not) predict
# =============================================================================
#  Designing a construct is not predicting its signal.  A KASP melt curve, a
#  droplet partition, a Cas collateral-cleavage burst - most chemistries here are
#  NOT simulated, and the honest thing is to say so rather than let a design read
#  as a promise of kinetics.  `kasp_readout()` is the one entry point: a scheme
#  that *does* model its readout attaches a `$readout` table to the marker
#  (currently the hydrolysis-probe discrimination, see `probe_discrimination()`);
#  every other scheme gets an explicit "not simulated" note.
#
#  The assumption notes and the model name are structured metadata (English, like
#  the scheme gate labels), not translated messages.
# =============================================================================

#' Per-allele discrimination of hydrolysis probes (DECIPHER) - a readout model
#'
#' The probe-level analogue of the KASP primer discrimination: each allele-specific
#' probe is scored for on-/off-target hybridisation+cleavage efficiency against its
#' own and the other allele's binding site ([pcr_efficiency()], DECIPHER).  The
#' "on" target is the perfect complement of the probe; the "off" target differs at
#' the variant base.  `start`/`end` are TEMPLATE coordinates of the probe window.
#' Needs OligoArrayAux.
#' @param seq the plus-strand template (`DNAString`).
#' @param snp_i template coordinate of the variant.
#' @param probes a tibble (`allele`, `probe`, `start`, `end`, `strand`), template coords.
#' @param alt the ALT allele for a single-base site.
#' @return a tibble (`allele`, `probe`, `on`, `off`, `d_eff`), one row per probe.
#' @keywords internal
#' @noRd
probe_discrimination <- function(seq, snp_i, probes, alt) {
  s <- as.character(seq)
  ref_base <- substr(s, snp_i, snp_i)
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
  comp <- function(x) as.character(Biostrings::complement(Biostrings::DNAString(x)))
  purrr::pmap(list(probes$allele, probes$probe, probes$start, probes$end, probes$strand), function(al, pr, a, b, st) {
    other <- if (identical(al, "REF")) alt else ref_base
    # the on target is the perfect complement; the off target differs at the
    # variant base (index `j` in the target, counting from its 5' end)
    t_on <- rc(pr)
    vp <- if (identical(st, "+")) snp_i - a + 1L else b - snp_i + 1L
    j <- nchar(pr) - vp + 1L
    t_off <- t_on
    stringr::str_sub(t_off, j, j) <- if (identical(st, "+")) comp(other) else other
    on <- pcr_efficiency(pr, t_on)
    off <- pcr_efficiency(pr, t_off)
    tibble::tibble(allele = al, probe = pr, on = on, off = off, d_eff = on - off)
  }) |>
    purrr::list_rbind()
}

#' The readout assumption for a scheme that does not simulate its signal
#' @keywords internal
#' @noRd
readout_assumption <- function(scheme) {
  notes <- c(
    hrm = "melting curve not modelled (needs a uMELT/PolyMelt-class model)",
    dpcr = "droplet partitioning / Poisson counting not modelled",
    crispr = "Cas collateral-cleavage kinetics not modelled",
    lamp = "isothermal amplification kinetics not modelled",
    rpa = "isothermal amplification kinetics not modelled",
    caps = "gel band sizes reported; digestion not simulated",
    arms = "gel band sizes reported; the reaction is not simulated",
    indel_size = "band sizes reported; the gel is not simulated",
    ssr = "fragment sizes reported; the gel is not simulated",
    sanger = "read / peak model not simulated",
    snapshot = "peak heights not simulated",
    castpcr = "real-time Delta-Ct not simulated (the probe is locus-specific)",
    snp = "melt / FRET profile not simulated",
    indel = "melt / FRET profile not simulated",
    multi = "melt / FRET profile not simulated",
    lax = "melt / FRET profile not simulated",
    duplex = "the two-locus signal not simulated",
    asq = "extension signal not simulated",
    primer = "no readout (primers only)"
  )
  if (scheme %in% names(notes)) notes[[scheme]] else "the signal is not simulated for this chemistry (constructs only)"
}

#' The readout a marker's scheme models (or says it does not)
#'
#' A design carries constructs, not a signal: this returns the readout model's
#' output when the scheme simulated one (a `table`, currently the hydrolysis-probe
#' discrimination) and otherwise an explicit assumption of what is *not* modelled,
#' so a design never reads as a promise of kinetics.
#' @param marker a designed marker (a `kasp_*` object).
#' @return a list: `scheme`, `simulated` (logical), `model`, `note`, `table`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' kasp_readout(design_taqman_marker(genome, kasp_example("SNP")[1, ]))
#' }
#' @family KASP design
#' @export
kasp_readout <- function(marker) {
  scheme <- sub("^kasp_", "", class(marker)[1L])
  if (!is.null(marker$readout)) {
    return(list(
      scheme = scheme, simulated = TRUE,
      model = "hydrolysis-probe discrimination (DECIPHER)",
      note = NA_character_, table = marker$readout
    ))
  }
  list(
    scheme = scheme, simulated = FALSE, model = NA_character_,
    note = readout_assumption(scheme), table = NULL
  )
}
