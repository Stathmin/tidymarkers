# =============================================================================
#  kasp_qc.R  --  repeats, palindromes, structural QC
# =============================================================================
# =============================================================================
#  5. Structural QC : Tm, hairpins, self- and cross-complementarity
# =============================================================================

#' Full structural profile of one primer
#' @param primer oligo(s) 5'->3'.
#' @param BPPARAM optional `BiocParallel` backend; the hairpin and terminal-dG
#'   maps are independent OligoArrayAux subprocesses.
#' @return a tibble, one row per primer (`primer`, `len`, `tm`, `gc`,
#'   `hairpin_dg`, `terminal_dg`).
#' @examples
#' \dontrun{
#' primer_thermo("GAAGGTGACCAAGTTCATGCT")
#' }
#' @family KASP thermodynamics
#' @export
primer_thermo <- function(primer, BPPARAM = NULL) {
  s <- as.character(primer)
  tibble::tibble(
    primer       = s,
    len          = nchar(s),
    tm           = primer_tm(s),
    gc           = gc_frac(s),
    hairpin_dg   = oligo_dg_vec(s, BPPARAM = BPPARAM), # hybrid-ss-min, kcal/mol
    terminal_dg  = terminal_dg(s, BPPARAM = BPPARAM)
  ) # 3'-terminal n-mer duplex, kcal/mol
}

#' Duplex free energy of a primer's 3'-terminal n-mer (3'-end stability)
#'
#' The last `n` bases are annealed to their reverse complement and scored with
#' OligoArrayAux (kcal/mol, 37 C); more negative = a more stable 3' anchor, the
#' "GC clamp" a polymerase needs to initiate.  Used as a design diagnostic.
#' @param primer oligo(s) 5'->3'.
#' @param n terminal window length.
#' @param BPPARAM optional `BiocParallel` backend; the distinct terminal n-mers
#'   are annealed in parallel.
#' @return a numeric vector of dG (kcal/mol), one per primer (`NA_real_` when
#'   OligoArrayAux is unavailable).
#' @examples
#' \dontrun{
#' terminal_dg("GAAGGTGACCAAGTTCATGCT")
#' }
#' @family KASP thermodynamics
#' @export
terminal_dg <- function(primer, n = 5L, BPPARAM = NULL) {
  s <- as.character(primer)
  mer <- substr(s, nchar(s) - n + 1L, nchar(s))
  mer[nchar(s) < n] <- NA_character_
  um <- unique(mer[!is.na(mer)])
  dg <- rlang::set_names(rep(NA_real_, length(um)), um)
  if (length(um)) {
    rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(um)))
    dg <- rlang::set_names(oligo_dg_vec(um, rc, BPPARAM = BPPARAM), um)
  }
  unname(dg[mer])
}

#' Is the 3'-terminal allelic mismatch a weak "wobble" pair (G:T / A:C)?
#'
#' `on_allele`/`off_allele` are the plus-strand alleles.  The primer carries the
#' "on" base and faces the complement of the "off" base, so the mismatch set is
#' `{on, complement(off)}`; G:T and A:C are the least destabilising and give the
#' poorest allelic discrimination (LGC/KBioscience note).  Reported as a flag.
#' @param on_allele,off_allele plus-strand alleles (the primer carries `on`, the
#'   other allele is `off`).
#' @return a logical vector, one per allele pair.
#' @examples
#' wobble_mismatch("A", "G") # TRUE (A:C mismatch)
#' wobble_mismatch("A", "C") # FALSE
#' @family KASP thermodynamics
#' @export
wobble_mismatch <- function(on_allele, off_allele) {
  faced <- as.character(Biostrings::complement(Biostrings::DNAStringSet(toupper(off_allele))))
  pair <- purrr::map2_chr(
    toupper(on_allele), faced,
    ~ paste(sort(c(.x, .y)), collapse = "")
  )
  pair %in% c("GT", "AC")
}

#' Duplex (dimer) free energy for every primer pair (OligoArrayAux, DNA)
#'
#' `dg` is the hybridisation free energy of the two strands (kcal/mol, 37 C);
#' more negative = stronger dimer.  The diagonal is the homodimer of each primer.
#' Returns the FULL matrix (both orders); symmetry is imposed here, once.  Each
#' pair is one OligoArrayAux subprocess, so for large primer sets pass a
#' `BPPARAM` to run the n(n+1)/2 pairs in parallel.
#' @param primers named character vector of oligos.
#' @param BPPARAM optional `BiocParallel` backend; pairs are split across workers.
#' @return a tibble of every ordered primer pair (`a`, `b`, `dg`; kcal/mol).
#' @examples
#' \dontrun{
#' dimer_thermo(c(AS = "GAAGGTGACCAAGTTCATGCT", COM = "AGCATGAACTTGGTCACCTTC"))
#' }
#' @family KASP thermodynamics
#' @export
dimer_thermo <- function(primers, BPPARAM = NULL) {
  nm <- names(primers)
  # upper triangle (i <= j) is all DISTINCT pairs; the full grid is then a join
  grid <- tidyr::crossing(a = nm, b = nm) |>
    dplyr::mutate(
      lo = pmin(match(.data$a, nm), match(.data$b, nm)),
      hi = pmax(match(.data$a, nm), match(.data$b, nm))
    )
  keys <- dplyr::distinct(grid, .data$lo, .data$hi)
  one <- function(k) oligo_dg(primers[[keys$lo[k]]], primers[[keys$hi[k]]])
  keys$dg <- if (is.null(BPPARAM)) {
    purrr::map_dbl(seq_len(nrow(keys)), one)
  } else {
    purrr::list_c(bp_lapply(seq_len(nrow(keys)), one, BPPARAM))
  }
  grid |>
    dplyr::left_join(keys, by = c("lo", "hi")) |>
    dplyr::transmute(a = .data$a, b = .data$b, dg = .data$dg)
}
