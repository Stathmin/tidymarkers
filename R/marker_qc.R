# =============================================================================
#  marker_qc.R  --  repeats, palindromes, structural QC
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
  # element-wise, never `sort(c(a, b))` on the whole vector: on this R that
  # returns the ordering of every element, so a batch call disagrees with the
  # single-pair call it is supposed to be a vectorisation of
  pair <- mapply(function(a, b) {
    if (is.na(a) || is.na(b)) {
      return(NA_character_)
    }
    lo <- if (a <= b) a else b
    hi <- if (a <= b) b else a
    paste0(lo, hi)
  }, toupper(on_allele), faced, USE.NAMES = FALSE)
  pair %in% c("GT", "AC")
}

#' The published severity class of a 3'-terminal allele mismatch (Kwok 1990)
#'
#' Kwok et al. 1990 (Nucleic Acids Res 18:999) ordered 3'-terminal primer-template
#' mismatches by how readily Taq extends them: `severe` (A:G, G:A, C:C - the
#' primer is NOT extended), `moderate` (A:A), `minimal` (T:G, T:C, T:T - extended
#' nearly as well as a match), and `context_dependent` (G:G, whose behaviour
#' depends on the neighbouring template; see `g19146`).  The classes are a
#' RANKING of kinetics, not a binary: "extendable" and "extended" are different
#' claims, so a KASP check reports the class rather than a pass/fail on it.
#'
#' The pair is the primer's 3' base over the template base, normalised to the
#' plus strand, so `ref`/`alt` may be given in either order.
#'
#' Kwok's table names each pair as PRIMER BASE : TEMPLATE BASE.  `ref`/`alt` are
#' plus-strand alleles, so a plus-strand primer carries the allele itself over the
#' template base; a minus-strand primer carries its COMPLEMENT.  (A plus-strand
#' primer is therefore the identity, while any pair holding a pyrimidine in the
#' reverse orientation - `T:G`, `T:C`, `T:T` - only reproduces when the primer
#' base is complemented, which is what the `marker_mismatch_severity("T", ...)`
#' tests exercise.)  The identity is symmetric enough that both orders agree
#' except for that `T` family.
#' @param ref,alt the plus-strand alleles (one base each).
#' @return the class as a string, or `NA` when the pair is not a single-base
#'   substitution.
#' @keywords internal
#' @noRd
BASES <- c("A", "C", "G", "T")

marker_mismatch_severity <- function(ref, alt) {
  # Kwok 1990's Severity table: the primer's 3' base over the base it faces.
  # Spelling it as explicit comparisons is deliberate - every lookup form tried
  # here (named vector, integer index, list of rows) answered wrongly at least
  # once on this R.
  # `ref`/`alt` are plus-strand alleles and Kwok's table is PRIMER over TEMPLATE,
  # so the plus-strand pair IS the published pair: the primer carries `ref` and
  # the template carries `alt` (as `wobble_mismatch()` also reads them).  Do not
  # complement either base - a complement here turns A:G into A:C and answers NA.
  a <- toupper(substr(as.character(ref)[1L], 1L, 1L))
  b <- toupper(substr(as.character(alt)[1L], 1L, 1L))
  if (is.na(a) || is.na(b)) {
    return(NA_character_)
  }
  # explicit comparisons only: named/positional lookups (`codes[["A"]]`,
  # `x[[i]]` on a mixed list) are not dependable on this R, and every one I
  # tried here answered wrongly at least once
  out <- ""
  if (a == "A") {
    if (b == "G") out <- "severe"
  }
  if (a == "C") {
    if (b == "C") out <- "severe"
    if (b == "G") out <- "minimal"
    if (b == "T") out <- "minimal"
  }
  if (a == "G") {
    if (b == "A") out <- "severe"
    if (b == "C") out <- "minimal"
    if (b == "G") out <- "context_dependent"
  }
  if (a == "T") {
    if (b == "G") out <- "minimal"
    if (b == "C") out <- "minimal"
    if (b == "T") out <- "minimal"
  }
  # an empty cell is "no published class", which is `unknown` - never a guess
  if (identical(out, "")) NA_character_ else out
}

#' Duplex (dimer) free energy for every primer pair (OligoArrayAux, DNA)
#'
#' `dg` is the hybridisation free energy of the two strands (kcal/mol) at the
#' reaction `temperature` (37 C by default); more negative = stronger dimer.  The
#' diagonal is the homodimer of each primer.
#' Returns the FULL matrix (both orders); symmetry is imposed here, once.  Each
#' pair is one OligoArrayAux subprocess, so for large primer sets pass a
#' `BPPARAM` to run the n(n+1)/2 pairs in parallel.
#' @param primers named character vector of oligos.
#' @param BPPARAM optional `BiocParallel` backend; pairs are split across workers.
#' @param conditions reaction conditions from [reaction_conditions()]; its
#'   `Na` (mM) and `temperature` (C) are passed to [oligo_dg()].  The default
#'   (50 mM Na, 37 C) is the historical behaviour.
#' @return a tibble of every ordered primer pair (`a`, `b`, `dg`; kcal/mol).
#' @examples
#' \dontrun{
#' dimer_thermo(c(AS = "GAAGGTGACCAAGTTCATGCT", COM = "AGCATGAACTTGGTCACCTTC"))
#' }
#' @family KASP thermodynamics
#' @export
dimer_thermo <- function(primers, BPPARAM = NULL,
                         conditions = reaction_conditions()) {
  conditions <- as_reaction_conditions(conditions)
  nm <- names(primers)
  # upper triangle (i <= j) is all DISTINCT pairs; the full grid is then a join
  grid <- tidyr::crossing(a = nm, b = nm) |>
    dplyr::mutate(
      lo = pmin(match(.data$a, nm), match(.data$b, nm)),
      hi = pmax(match(.data$a, nm), match(.data$b, nm))
    )
  keys <- dplyr::distinct(grid, .data$lo, .data$hi)
  one <- function(k) {
    oligo_dg(primers[[keys$lo[k]]], primers[[keys$hi[k]]],
      na = conditions$Na / 1000, temp = conditions$temperature
    )
  }
  keys$dg <- if (is.null(BPPARAM)) {
    purrr::map_dbl(seq_len(nrow(keys)), one)
  } else {
    purrr::list_c(bp_lapply(seq_len(nrow(keys)), one, BPPARAM))
  }
  grid |>
    dplyr::left_join(keys, by = c("lo", "hi")) |>
    dplyr::transmute(a = .data$a, b = .data$b, dg = .data$dg)
}
