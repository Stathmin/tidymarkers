# =============================================================================
#  kasp_fm.R  --  an FM-index genome scan (Rsubread), a BLAST-index alternative
# =============================================================================
#  `specificity_seed()` counts 3'-seed matches with Biostrings and is perfect up
#  to the genome sizes you can hold in RAM; `specificity_blast()` queries an NCBI
#  BLAST database (a heuristic index).  A third, fully index-based option is an
#  FM-index: build it once with `Rsubread::buildindex()` and map primers to it
#  with `Rsubread::align()`, reporting every genomic location a primer matches.
#  That answers a *whole-primer* question ("where else would this primer bind?")
#  rather than the 3'-seed question, and it streams the genome off disk; the
#  mapping itself is multithreaded (`nthreads`), so large genomes are tractable.
# =============================================================================

# Mismatches tolerated when mapping the whole primer with the FM-index.
FM_MISMATCH <- 2L

#' Build an FM-index for a genome (`Rsubread::buildindex`)
#'
#' Writes the (uncompressed) genome FASTA once and builds a Subread FM-index next
#' to it.  Idempotent: if `basename`'s index files already exist they are reused,
#' so calling this on every run only pays the build cost once.  Attach the
#' returned basename to a `kasp_genome` with `kasp_genome(..., fm_index = <basename>)`
#' to let the assembly use it, or pass it straight to [specificity_fm()].
#'
#' @param genome a `kasp_genome` or a `DNAStringSet`.
#' @param basename index basename to create (no extension).
#' @return the index basename (invisibly used by [specificity_fm()]).
#' @examples
#' \dontrun{
#' idx <- kasp_fm_index(kasp_genome("genome.fa"))
#' }
#' @family KASP specificity
#' @export
kasp_fm_index <- function(genome, basename = fs::path(fs::path_temp(), "kasp_fm")) {
  rlang::check_installed("Rsubread", reason = gettext("for the FM-index genome scan"))
  basename <- as.character(basename) # Rsubread wants a plain string
  if (fs::file_exists(paste0(basename, ".files"))) {
    cli::cli_alert_info(gettext("Reusing the existing FM-index at {.path {basename}}."))
    return(basename)
  }
  fasta <- if (inherits(genome, "kasp_genome")) {
    genome$path
  } else {
    p <- paste0(basename, ".fa")
    fs::dir_create(fs::path_dir(p))
    Biostrings::writeXStringSet(genome, p)
    p
  }
  cli::cli_alert_info(gettext("Building a Subread FM-index at {.path {basename}} (one-off)."))
  # buildindex reports progress on stdout with no quiet switch
  invisible(utils::capture.output(Rsubread::buildindex(basename = basename, reference = fasta)))
  basename
}

#' Count genome hits of whole primers with an FM-index (`Rsubread`)
#'
#' Maps each primer to the genome with `Rsubread::align()` (`unique = FALSE`, so
#' every location is reported) and counts the alignments per primer.  This is the
#' index-based counterpart of [blast_primers()]: a whole-primer scan, not the
#' 3'-seed scan of [specificity_seed()].  The mapping is multithreaded via
#' `nthreads` and never loads the genome into R.
#'
#' @param primers named character vector of oligos (5'->3').
#' @param index FM-index basename from [kasp_fm_index()].
#' @param max_mismatches mismatches tolerated over the whole primer.
#' @param nthreads threads for `Rsubread::align()`.
#' @return tibble `primer`, `seed` (the primer itself) and `hits` (genome-wide
#'   alignment count; more than one = a potential off-target).
#' @examples
#' \dontrun{
#' idx <- kasp_fm_index(kasp_genome("genome.fa"))
#' specificity_fm(c(AS = "ACGTACGTACGTACGT"), idx)
#' }
#' @family KASP specificity
#' @export
specificity_fm <- function(primers, index, max_mismatches = FM_MISMATCH, nthreads = 1L) {
  rlang::check_installed("Rsubread", reason = gettext("for the FM-index genome scan"))
  rlang::check_installed("Rsamtools", reason = gettext("to read the alignments"))
  index <- as.character(index) # Rsubread wants a plain string
  td <- as.character(withr::local_tempdir())
  reads <- file.path(td, "reads.fa")
  bam <- file.path(td, "hits.bam")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(stats::setNames(as.character(primers), names(primers))), reads
  )
  # Rsubread prints a progress report to stdout with no quiet switch; the result
  # is the alignment file, so the chit-chat is discarded
  invisible(utils::capture.output(
    Rsubread::align(
      index = index, readfile1 = reads, output_file = bam, type = 1,
      unique = FALSE, nBestLocations = 10000L,
      maxMismatches = max_mismatches, nthreads = nthreads
    )
  ))
  q <- as.character(Rsamtools::scanBam(bam)[[1L]]$qname)
  cnt <- table(q)
  hits <- as.integer(cnt[names(primers)])
  hits[is.na(hits)] <- 0L
  tibble::tibble(primer = names(primers), seed = as.character(primers), hits = hits)
}
