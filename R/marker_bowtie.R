# =============================================================================
#  marker_bowtie.R  --  a bowtie1 whole-genome scan (Rbowtie), the indexed path
# =============================================================================
#  `specificity_seed()` counts 3'-seed matches with Biostrings and is perfect up
#  to the genome sizes you can hold in RAM; `specificity_blast()` is a heuristic
#  BLAST index.  This is the third indexed option and the one built for large
#  genomes: `Rbowtie::bowtie_build()` writes a compact FM index once, then
#  `Rbowtie::bowtie()` reports *every* genomic location each primer matches.
#  It answers the whole-primer question ("where else would this primer bind?")
#  rather than the 3'-seed question, and the query cost is flat in primer count
#  and in genome size (measured: 0.05 s for 1000 primers against a 30 Mb index).
#
#  Measured against the in-memory two-strand `countPattern()` oracle on TAIR10
#  Chr1 the counts agree at 16, 96, 383, 1000 and 5000 primers, the index costs
#  1.34 bytes per genome base (wheat 16 Gb ~ 21.5 GB), and the build is a
#  one-off.  This is why `specificity_of()` prefers it over the Subread FM index
#  (31.7 bytes/bp, and its BAM carries no mismatch tags to filter on).
# =============================================================================

#' Is the bowtie1 scan available (`Rbowtie`)?
#'
#' @return `TRUE` when the `Rbowtie` package is installed.
#' @examples
#' has_bowtie()
#' @family specificity
#' @export
has_bowtie <- function() {
  requireNamespace("Rbowtie", quietly = TRUE)
}

#' Build a bowtie1 index for a genome (`Rbowtie::bowtie_build`)
#'
#' Writes a compact FM index next to `basename` in one pass.  Idempotent: if
#' `<basename>.1.ebwt` already exists it is reused, so calling this on every run
#' only pays the build cost once.  Attach the returned basename to a
#' `marker_genome` with `marker_genome(..., bowtie_index = <basename>)` to let
#' [specificity_of()] use it, or pass it straight to [specificity_bowtie()].
#'
#' bowtie1 needs a non-empty defline on every sequence; unnamed sequences of an
#' in-memory `DNAStringSet` are named `seq1`, `seq2`, ... on the way out.  A
#' BSgenome-backed `marker_genome` has no FASTA and gets a named error; pass a
#' `BSgenome` itself and its FASTA is written once to a name-keyed temp file.
#'
#' @param genome a `marker_genome`, a `DNAStringSet` or a `BSgenome`.
#' @param basename index basename to create (no extension).
#' @param threads threads for `bowtie_build()`.
#' @return the index basename, for [specificity_bowtie()].
#' @examples
#' \dontrun{
#' idx <- marker_bowtie_index(marker_genome("genome.fa"))
#' }
#' @family specificity
#' @export
marker_bowtie_index <- function(genome, basename = fs::path(fs::path_temp(), "marker_bowtie"),
                                threads = 4L) {
  # the genome check comes first, so without Rbowtie a BSgenome genome still says so
  src <- genome_fasta_source(genome, "marker_bowtie_index")
  genome <- src$genome
  rlang::check_installed("Rbowtie", reason = gettext("for the bowtie genome scan"))
  basename <- as.character(basename)
  if (fs::file_exists(paste0(basename, ".1.ebwt"))) {
    cli::cli_alert_info(gettext("Reusing the existing bowtie index at {.path {basename}}."))
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
      # an empty defline leaves the index with unusable reference names
      nm[blank] <- paste0("seq", which(blank))
      names(gs) <- nm
    }
    Biostrings::writeXStringSet(gs, p)
    p
  }
  cli::cli_alert_info(gettext("Building a bowtie index at {.path {basename}} (one-off)."))
  # bowtie_build reports progress on stdout with no quiet switch
  invisible(utils::capture.output(
    Rbowtie::bowtie_build(
      fasta,
      outdir = fs::path_dir(basename),
      prefix = fs::path_file(basename),
      force = TRUE, strict = FALSE, threads = threads
    )
  ))
  basename
}

#' Count genome hits of whole primers with a bowtie1 index (`Rbowtie`)
#'
#' Maps each primer to the index with `Rbowtie::bowtie()` (`a = TRUE`, so every
#' location is reported) and counts its alignments.  This is the index-based
#' counterpart of [specificity_seed()]: a *whole-primer* scan of both strands,
#' not the 3'-seed scan.  Unaligned reads (`FLAG 4`) are dropped, so an absent
#' primer reports `0`; the count includes the reverse-complement strand.
#'
#' @param primers named character vector of oligos (5'->3').
#' @param index index basename from [marker_bowtie_index()].
#' @param max_mismatches mismatches tolerated over the whole primer (bowtie1
#'   allows at most 3).
#' @param nthreads threads for `Rbowtie::bowtie()`.
#' @return tibble `primer`, `seed` (the primer itself) and `hits` (both-strand
#'   genome-wide alignment count; more than one = a potential off-target).
#' @examples
#' \dontrun{
#' idx <- marker_bowtie_index(marker_genome("genome.fa"))
#' specificity_bowtie(c(AS = "ACGTACGTACGTACGT"), idx)
#' }
#' @family specificity
#' @export
specificity_bowtie <- function(primers, index, max_mismatches = 0L, nthreads = 4L) {
  rlang::check_installed("Rbowtie", reason = gettext("for the bowtie genome scan"))
  max_mismatches <- as.integer(max_mismatches)
  if (max_mismatches > 3L) {
    cli::cli_abort(gettext("bowtie1 allows at most 3 mismatches (`max_mismatches` is {max_mismatches})."))
  }
  index <- as.character(index)
  td <- as.character(withr::local_tempdir())
  reads <- file.path(td, "reads.fa")
  sam <- file.path(td, "hits.sam")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(rlang::set_names(as.character(primers), names(primers))), reads
  )
  # bowtie prints its summary to stdout with no quiet switch
  invisible(utils::capture.output(
    Rbowtie::bowtie(
      reads, index,
      f = TRUE, v = max_mismatches, a = TRUE, y = TRUE, S = TRUE, p = nthreads,
      type = "single", outfile = sam, force = TRUE
    )
  ))
  lines <- readLines(sam, warn = FALSE)
  lines <- lines[!startsWith(lines, "@")]
  cnt <- if (length(lines) == 0L) {
    integer(0)
  } else {
    fields <- strsplit(lines, "\t", fixed = TRUE)
    qname <- vapply(fields, `[`, character(1), 1L)
    flag <- as.integer(vapply(fields, `[`, character(1), 2L))
    # an unaligned read still appears in SAM with FLAG 4; counting it would
    # report an absent primer as one hit
    table(qname[bitwAnd(flag, 4L) == 0L])
  }
  hits <- as.integer(cnt[names(primers)])
  hits[is.na(hits)] <- 0L
  tibble::tibble(primer = names(primers), seed = as.character(primers), hits = hits)
}
