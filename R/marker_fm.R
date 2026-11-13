# =============================================================================
#  marker_fm.R  --  an FM-index genome scan (Rsubread), a BLAST-index alternative
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
#' returned basename to a `marker_genome` with `marker_genome(..., fm_index = <basename>)`
#' to let the assembly use it, or pass it straight to [specificity_fm()].
#'
#' Subread needs a non-empty defline on every sequence; unnamed sequences of an
#' in-memory `DNAStringSet` are named `seq1`, `seq2`, ... on the way out.  A
#' BSgenome-backed `marker_genome` has no FASTA and gets a named error; pass a
#' `BSgenome` itself and its FASTA is written once to a content-keyed temp file.
#'
#' Subread also silently rewrites any non-ACGT base to `A` when it builds the
#' index, so hit counts over N-rich regions are conservative; the notice is
#' replayed as a warning instead of being swallowed by the progress capture.
#'
#' The index is **not** used by the design pipeline: [specificity_of()] refuses
#' an FM-only genome and points at [marker_bowtie_index()], because the Subread
#' index costs ~24x more per base and [specificity_fm()] is a read mapper, not a
#' seed counter.  Building one is still supported for that use, and says so.
#'
#' @param genome a `marker_genome`, a `DNAStringSet` or a `BSgenome`.
#' @param basename index basename to create (no extension).
#' @return the index basename (invisibly used by [specificity_fm()]).
#' @examples
#' \dontrun{
#' idx <- marker_fm_index(marker_genome("genome.fa"))
#' }
#' @family specificity
#' @export
marker_fm_index <- function(genome, basename = fs::path(fs::path_temp(), "marker_fm")) {
  # the genome check comes first, so without Rsubread a BSgenome genome still says so
  src <- genome_fasta_source(genome, "marker_fm_index")
  genome <- src$genome
  rlang::check_installed("Rsubread", reason = gettext("for the FM-index genome scan"))
  basename <- as.character(basename) # Rsubread wants a plain string
  # a buildable index that the design pipeline never queries is a trap (D5)
  cli::cli_alert_warning(gettext(
    "The FM index is not used by the design pipeline: specificity_fm() is a read mapper.
     For the design scan use marker_bowtie_index() / specificity_bowtie()."
  ))
  if (fs::file_exists(paste0(basename, ".files"))) {
    cli::cli_alert_info(gettext("Reusing the existing FM-index at {.path {basename}}."))
    return(basename)
  }
  fasta <- if (inherits(genome, "marker_genome")) {
    genome$path
  } else {
    p <- paste0(basename, ".fa")
    fs::dir_create(fs::path_dir(p))
    gs <- genome
    nm <- names(gs)
    if (is.null(nm)) nm <- rep("", length(gs))
    blank <- is.na(nm) | !nzchar(nm)
    if (any(blank)) {
      # an empty defline makes the BAM header unusable downstream
      nm[blank] <- paste0("seq", which(blank))
      names(gs) <- nm
    }
    Biostrings::writeXStringSet(gs, p)
    p
  }
  cli::cli_alert_info(gettext("Building a Subread FM-index at {.path {basename}} (one-off)."))
  # buildindex reports progress on stdout with no quiet switch
  invisible(utils::capture.output(Rsubread::buildindex(basename = basename, reference = fasta)))
  log <- paste0(basename, ".log")
  rewritten <- fs::file_exists(log) &&
    any(grepl("non-ACGT", readLines(log, warn = FALSE), fixed = TRUE))
  if (rewritten) {
    cli::cli_alert_warning(gettext(
      "Subread rewrote the non-ACGT bases of this genome to 'A' when building the index;
       counts over ambiguous regions are conservative."
    ))
  }
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
#' @param index FM-index basename from [marker_fm_index()].
#' @param max_mismatches mismatches tolerated over the whole primer.
#' @param nthreads threads for `Rsubread::align()`.
#' @return tibble `primer`, `seed` (the primer itself) and `hits` (genome-wide
#'   alignment count; more than one = a potential off-target).
#' @examples
#' \dontrun{
#' idx <- marker_fm_index(marker_genome("genome.fa"))
#' specificity_fm(c(AS = "ACGTACGTACGTACGT"), idx)
#' }
#' @family specificity
#' @export
specificity_fm <- function(primers, index, max_mismatches = FM_MISMATCH, nthreads = 1L) {
  rlang::check_installed("Rsubread", reason = gettext("for the FM-index genome scan"))
  rlang::check_installed("Rsamtools", reason = gettext("to read the alignments"))
  index <- as.character(index) # Rsubread wants a plain string
  td <- as.character(withr::local_tempdir())
  reads <- file.path(td, "reads.fa")
  bam <- file.path(td, "hits.bam")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(rlang::set_names(as.character(primers), names(primers))), reads
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
  b <- Rsamtools::scanBam(bam)[[1L]]
  # an unaligned read still appears in the BAM with FLAG 4; counting it would
  # report an absent primer as one hit
  q <- as.character(b$qname[bitwAnd(b$flag, 4L) == 0L])
  cnt <- table(q)
  hits <- as.integer(cnt[names(primers)])
  hits[is.na(hits)] <- 0L
  tibble::tibble(primer = names(primers), seed = as.character(primers), hits = hits)
}
