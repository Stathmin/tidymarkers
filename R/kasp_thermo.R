# =============================================================================
#  kasp_thermo.R  --  Tm, dG and PCR efficiency (packages)
# =============================================================================
# =============================================================================
#  2. Thermodynamics via packages
# =============================================================================

#' GC fraction and homopolymer-run test
#' @keywords internal
#' @noRd
gc_frac <- function(seq) unname(Biostrings::letterFrequency(Biostrings::DNAStringSet(seq), "GC", as.prob = TRUE)[, 1])
has_run <- function(seq, n = 4L) stringr::str_detect(seq, glue::glue("([ACGT])\\1{{{n - 1L},}}"))

#' Make a parallel worker able to dispatch on Biostrings objects
#'
#' A `SnowParam`/SOCK worker is a fresh R process.  It resolves the package's
#' closures, but unless `Biostrings` is loaded there, S4 dispatch on a serialized
#' `DNAStringSet` (e.g. `genome[[i]]`, `subseq()`) falls back to base and fails
#' with "this S4 class is not subsettable".  Loading the namespace in the worker
#' fixes the dispatch; fork-based backends ([BiocParallel::MulticoreParam()])
#' inherit loaded namespaces and are unaffected.  Called at the top of every
#' worker closure that touches the genome.
#' @keywords internal
#' @noRd
worker_biostrings <- function() requireNamespace("Biostrings", quietly = TRUE)

#' Parallel map with an interactive progress bar
#'
#' A thin wrapper so every parallel map shows a progress bar only when
#' interactive.  Reproducible RNG is the backend's job: set `RNGseed` on the
#' `BPPARAM` (`MulticoreParam(RNGseed = 1)`); the design kernels here are
#' deterministic, so nothing needs seeding (`bpoptions(RNGseed = )` is not
#' supported by BiocParallel).
#' @keywords internal
#' @noRd
bp_lapply <- function(X, FUN, BPPARAM) {
  BiocParallel::bplapply(
    X, FUN,
    BPPARAM = BPPARAM,
    BPOPTIONS = BiocParallel::bpoptions(progressbar = interactive())
  )
}

#' Tandem repeats (simple sequence repeats) in a sequence
#'
#' EasyKASP screens candidate primers for SSR repeat units, not just
#' homopolymers.  Returns, for each input, the repeated units of length `units`
#' occurring at least `min_copies` times in tandem.
#' @param seq character vector of DNA sequences.
#' @param units repeat-unit lengths to test.
#' @param min_copies tandem copies required for a hit.
#' @return a list, one character vector per input, of the repeated units (unique,
#'   possibly empty).
#' @examples
#' find_repeats(c("ACACACAC", "ACGTACGT"))
#' @family KASP QC
#' @export
find_repeats <- function(seq, units = 2:4, min_copies = 3L) {
  purrr::map(seq, function(s) {
    hits <- character(0)
    for (u in units) {
      pat <- glue::glue("([ACGT]{{{u}}})\\1{{{min_copies - 1L},}}")
      m <- stringr::str_match_all(s, pat)[[1L]]
      if (nrow(m)) hits <- c(hits, m[, 2L])
    }
    unique(hits)
  })
}

#' Does a sequence contain an SSR? (see [find_repeats()])
#' @param seq character vector of DNA sequences.
#' @param units repeat-unit lengths to test.
#' @param min_copies tandem copies required for a hit.
#' @return a logical vector, one per input sequence.
#' @examples
#' has_repeat(c("ACACACAC", "ACGT"))
#' @family KASP QC
#' @export
has_repeat <- function(seq, units = 2:4, min_copies = 3L) {
  lengths(find_repeats(seq, units, min_copies)) > 0L
}

#' Inverted repeats / palindromes in a primer (library screen)
#'
#' `Biostrings::findPalindromes()` finds reverse-complementary arms - the
#' structural signature of a hairpin.  This complements the thermodynamic
#' hairpin screen in [primer_thermo()] with a sequence-level flag that does not
#' need OligoArrayAux.  The defaults are deliberately strict (arms >= 4 bp, loop
#' <= 6 bp, no mismatches): a well-designed primer returns no hits, whereas a
#' repetitive or self-symmetric one returns several.
#'
#' @param seq character vector of DNA sequences.
#' @param min_armlength minimum arm length of a palindromic hit.
#' @param max_looplength maximum loop (spacer) length.
#' @param max_mismatch mismatches allowed in the arms.
#' @return a list, one tibble per input, with `start`, `end`, `width` of each
#'   palindromic view (empty tibble when none).
#' @examples
#' find_palindromes(c(hit = "ACGTACGT", none = "AACCGG"))
#' @family KASP QC
#' @export
find_palindromes <- function(seq, min_armlength = 4L, max_looplength = 6L, max_mismatch = 0L) {
  nm <- names(seq)
  hits <- purrr::map(as.character(seq), function(s) {
    h <- Biostrings::findPalindromes(Biostrings::DNAString(s),
      min.armlength = min_armlength,
      max.looplength = max_looplength,
      max.mismatch = max_mismatch
    )
    tibble::tibble(
      start = as.integer(IRanges::start(h)),
      end = as.integer(IRanges::end(h)),
      width = as.integer(IRanges::width(h))
    )
  })
  stats::setNames(hits, nm)
}

#' Does a primer contain an inverted repeat? (see [find_palindromes()])
#' @param seq character vector of DNA sequences.
#' @param min_armlength minimum arm length of a palindromic hit.
#' @param max_looplength maximum loop (spacer) length.
#' @param max_mismatch mismatches allowed in the arms.
#' @return a logical vector, one per input sequence.
#' @examples
#' has_palindrome(c("ACGTACGT", "AACCGG"))
#' @family KASP QC
#' @export
has_palindrome <- function(seq, min_armlength = 4L, max_looplength = 6L, max_mismatch = 0L) {
  purrr::map_int(find_palindromes(seq, min_armlength, max_looplength, max_mismatch), nrow) > 0L
}

#' Nearest-neighbour Tm of oligos (perfect duplex), salt-corrected (TmCalculator)
#'
#' `res$gr` uses TmCalculator's `$` accessor for the result GRanges; `Tm` is a
#' metadata column on it.
#' @param seq oligo(s) 5'->3'.
#' @param DMSO DMSO concentration (%), as in a KASP master mix; it lowers Tm.
#' @return a numeric vector of Tm (C), one per oligo.
#' @examples
#' primer_tm("GAAGGTGACCAAGTTCATGCT")
#' primer_tm(c("ACGTACGTACGTACGTACGT", "TTTTGGCGCGCCGATATATAT"))
#' @family KASP thermodynamics
#' @export
primer_tm <- function(seq, DMSO = PCR_COND$DMSO) {
  gr <- TmCalculator::to_genomic_ranges(as.character(seq))
  res <- TmCalculator::tm_calculate(gr,
    method = "tm_nn", nn_table = NN_TABLE,
    Na = PCR_COND$Na, Mg = PCR_COND$Mg, dNTPs = PCR_COND$dNTPs, DMSO = DMSO
  )
  as.numeric(S4Vectors::mcols(res$gr)$Tm)
}

#' Is OligoArrayAux (`hybrid-min`) available?
#'
#' DECIPHER's thermodynamics runs through OligoArrayAux.  It is looked up as
#' `options(tidykasp.oligoarrayaux)`, a project `.tools/`, the user cache (see
#' [kasp_install_tools()]) and finally the `PATH`.
#'
#' @return a single logical.
#' @examples
#' has_oligoarrayaux()
#' @family KASP utilities
#' @export
has_oligoarrayaux <- function() !is.null(oligoarrayaux_paths())

#' Point OligoArrayAux at its data tables for the duration of `code`
#'
#' The thermodynamics engine (DECIPHER via OligoArrayAux) needs `UNAFOLDDAT` to
#' find its parameter tables.  This temporarily sets it, so you need not fiddle
#' with the environment by hand.
#'
#' @param path directory with the OligoArrayAux data tables.
#' @param code expression to evaluate (e.g. `assemble_kasp_assay(...)`).
#' @return the value of `code`, evaluated with `UNAFOLDDAT` set to `path`.
#' @examples
#' \dontrun{
#' with_unafloddat("/path/to/oligoarrayaux", assemble_kasp_assay(genome, snp))
#' }
#' @family KASP utilities
#' @export
with_unafloddat <- function(path, code) {
  if (!fs::dir_exists(path)) rlang::abort(glue::glue(gettext("Directory not found: {path}")))
  withr::with_envvar(c(UNAFOLDDAT = path), code)
}

#' Turn a core count into a `BiocParallel` backend (or `NULL` for one core)
#'
#' The single place `cores` becomes a `BPPARAM`: fork on Unix, a worker pool
#' elsewhere, and `NULL` for a single core - the CRAN-friendly default.  Pass a
#' `BPPARAM` explicitly to override, which is also how a caller avoids nesting
#' (an outer locus-level backend should force the inner one to `NULL`).
#' @keywords internal
#' @noRd
kasp_bpparam <- function(cores = getOption("tidykasp.cores", 1L)) {
  n <- max(1L, as.integer(cores))
  if (n <= 1L) {
    return(NULL)
  }
  if (.Platform$OS.type == "unix") BiocParallel::MulticoreParam(n) else BiocParallel::SnowParam(n)
}

#' Worker count of a `BPPARAM` (1 when `NULL`)
#' @keywords internal
#' @noRd
kasp_ncores <- function(BPPARAM = NULL) {
  if (is.null(BPPARAM)) 1L else max(1L, BiocParallel::bpnworkers(BPPARAM))
}

#' OligoArrayAux dG for many oligos, optionally in parallel
#'
#' A vectorised `oligo_dg()` (duplex when `seq2` is given, else hairpin folding).
#' Each pair is its own `hybrid-min` subprocess, so with a `BPPARAM` the batch
#' runs across cores - the COM candidates and the QC maps are the per-SNP hot
#' spot.
#' @keywords internal
#' @noRd
oligo_dg_vec <- function(seq1, seq2 = NULL, BPPARAM = NULL) {
  n <- length(seq1)
  if (n == 0L) {
    return(numeric(0))
  }
  one <- function(i) if (is.null(seq2)) oligo_dg(seq1[[i]]) else oligo_dg(seq1[[i]], seq2[[i]])
  if (is.null(BPPARAM)) {
    vapply(seq_len(n), one, numeric(1))
  } else {
    unlist(bp_lapply(seq_len(n), one, BPPARAM), use.names = FALSE)
  }
}

#' PCR amplification efficiency of primer/target pairs (DECIPHER)
#'
#' Models hybridisation (incl. side reactions) times Taq elongation efficiency,
#' so a 3'-terminal mismatch is strongly penalised - the physical basis of KASP.
#'
#' @param primers character vector, 5'->3'.
#' @param targets character vector (same length) with the reverse-complement of
#'   the binding site (as in `?DECIPHER::CalculateEfficiencyPCR`).
#' @param temp annealing temperature (C).
#' @param P primer concentration (M).
#' @param ions monovalent + divalent salt (M).
#' @param cores processors for the vectorised DECIPHER call (default
#'   `getOption("tidykasp.cores", 1L)`; the many OligoArrayAux subprocesses are
#'   the hot spot, so >1 helps a lot, and 1 keeps CRAN happy).
#' @return a numeric vector of efficiencies (0-1), one per primer/target pair.
#' @examples
#' \dontrun{
#' pcr_efficiency("GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC")
#' }
#' @family KASP thermodynamics
#' @export
pcr_efficiency <- function(primers, targets, temp = PCR_EFF$temp,
                           P = PCR_EFF$P, ions = PCR_EFF$ions,
                           cores = getOption("tidykasp.cores", 1L)) {
  rlang::check_installed("DECIPHER")
  if (!has_oligoarrayaux()) {
    rlang::abort(c(gettext("OligoArrayAux ('hybrid-min') is required for PCR efficiencies."),
      i = gettext("run kasp_install_tools(), or build it with data-raw/install-oligoarrayaux.R")
    ))
  }
  with_oligoarrayaux(as.numeric(DECIPHER::CalculateEfficiencyPCR(
    primer = Biostrings::DNAStringSet(as.character(primers)),
    target = Biostrings::DNAStringSet(as.character(targets)),
    temp = temp, P = P, ions = ions, taqEfficiency = TRUE,
    processors = max(1L, as.integer(cores))
  )))
}
