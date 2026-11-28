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
#
#  A monolithic index is the lever's limit: wheat is 16 Gb, so one index is about
#  21.5 GB and no single build fits a small host.  `marker_bowtie_index(chroms =)`
#  slices the genome into one index per chromosome, and
#  `specificity_bowtie()` sums the hits over the set.  The wheat build-and-query
#  measurement is still owed: no 16 Gb fixpured reference on this host.
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
#' `BSgenome` itself and its FASTA is written once to a content-keyed temp file.
#'
#' @section Index location:
#' The default `basename` is under `tempdir()`, so a built index is never left
#' on disk by accident: it dies with the session.  That is cheap for a small
#' genome and wrong for wheat, where the index is about 1.4x the FASTA (about
#' 21.5 GB for 16 Gb).  Build once into the cache instead and clean up:
#' `marker_bowtie_index(genome, basename = file.path(marker_cache_dir(), "wheat"))`,
#' then `unlink(Sys.glob(file.path(marker_cache_dir(), "wheat.*")))` when done.
#'
#' @param genome a `marker_genome`, a `DNAStringSet` or a `BSgenome`.
#' @param basename index basename to create (no extension); defaults to a name
#'   under [tempdir()] (see *Index location*).
#' @param threads threads for `bowtie_build()`.
#' @param chroms optional character vector of sequence names (see
#'   [genome_chroms()]).  When given, one index per chromosome is built at
#'   `paste0(basename, "_", chrom)` and a named `marker_bowtie_index_set` is
#'   returned, which [specificity_bowtie()] sums hits over; an existing
#'   per-chromosome index is reused like the whole-genome one.  `NULL` (the
#'   default) builds the one whole-genome index.
#' @return the index basename, for [specificity_bowtie()].  With `chroms`, a
#'   named `marker_bowtie_index_set` of basenames, one per chromosome; subset it
#'   (`idx[c("chr1")]`) to query fewer chromosomes.
#' @examples
#' \dontrun{
#' idx <- marker_bowtie_index(marker_genome("genome.fa"))
#' }
#' @family specificity
#' @export
marker_bowtie_index <- function(genome, basename = fs::path(fs::path_temp(), "marker_bowtie"),
                                threads = 4L, chroms = NULL) {
  threads <- check_worker_count(threads, "threads")
  # the genome check comes first, so without Rbowtie a BSgenome genome still says so
  src <- genome_fasta_source(genome, "marker_bowtie_index")
  genome <- src$genome
  rlang::check_installed("Rbowtie", reason = gettext("for the bowtie genome scan"))
  basename <- as.character(basename)
  if (is.null(chroms)) {
    return(bowtie_index_build(genome, basename, threads))
  }
  if (length(chroms) == 0L) {
    abort_msg(
      "`chroms` must name at least one sequence, or be `NULL` for one whole-genome index."
    )
  }
  bowtie_index_set(genome, basename, as.character(chroms), threads)
}

# bowtie1 needs a non-empty defline on every sequence; unnamed sequences of an
# in-memory set are named `seq1`, `seq2`, ... on the way out
bowtie_name_seqs <- function(gs) {
  nm <- names(gs)
  if (is.null(nm)) nm <- rep("", length(gs))
  blank <- is.na(nm) | !nzchar(nm)
  if (any(blank)) {
    # an empty defline leaves the index with unusable reference names
    nm[blank] <- paste0("seq", which(blank))
    names(gs) <- nm
  }
  gs
}

# one index at `basename`, reusing an existing build (the one-off cost)
bowtie_index_build <- function(genome, basename, threads) {
  if (cache_present(paste0(basename, ".1.ebwt"))) {
    cli::cli_alert_info(gettext("Reusing the existing bowtie index at {.path {basename}}."))
    return(basename)
  }
  fasta <- if (inherits(genome, "marker_genome")) {
    genome$path
  } else {
    p <- paste0(basename, ".fa")
    fs::dir_create(fs::path_dir(p))
    declared_path(p)
    Biostrings::writeXStringSet(bowtie_name_seqs(genome), p)
    p
  }
  bowtie_build_fasta(fasta, basename, threads)
}

# the build half, separated so a per-chromosome shard can index the one-sequence
# FASTA it just wrote instead of writing the sequence a second time
bowtie_build_fasta <- function(fasta, basename, threads) {
  declared_index(basename, c(
    ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"
  ))
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

# one index per chromosome, so no single build has to hold the whole genome.
# Each chromosome is read on its own (`index_chrom_seq()`), written as a
# one-sequence FASTA and removed after the build: peak disc is one chromosome
# plus its index, peak memory one chromosome.
bowtie_index_set <- function(genome, basename, chroms, threads) {
  # `seqnames` maps the requested name to the FASTA defline (`chrom_map`)
  defline <- if (inherits(genome, "marker_genome")) genome$seqnames else rlang::set_names(names(genome))
  unknown <- chroms[is.na(match(chroms, names(defline)))]
  if (length(unknown) > 0L) {
    abort_msg(
      "{.arg chroms} asks for sequence{?s} {.val {unknown}}, which {.arg genome} does not have.",
      i = "Available: {.val {names(defline)}}."
    )
  }
  bases <- rlang::set_names(paste0(basename, "_", chroms), chroms)
  for (chr in chroms) {
    if (cache_present(paste0(bases[[chr]], ".1.ebwt"))) {
      cli::cli_alert_info(gettext("Reusing the existing bowtie index at {.path {bases[[chr]]}}."))
      next
    }
    fa <- paste0(bases[[chr]], ".fa")
    index_write_chrom_fasta(index_chrom_seq(genome, chr), chr, fa)
    bowtie_build_fasta(fa, bases[[chr]], threads)
    unlink(fa)
  }
  structure(bases, class = c("marker_bowtie_index_set", "character"), chroms = chroms)
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
#' @param index index basename from [marker_bowtie_index()], or the
#'   per-chromosome set it returns (any vector of two or more basenames).  With a
#'   set, every index is queried and the hits are summed per primer, so a
#'   chromosome that lacks a primer contributes `0`.
#' @param max_mismatches mismatches tolerated over the whole primer (bowtie1
#'   allows at most 3).
#' @param nthreads threads for `Rbowtie::bowtie()`.
#' @return tibble `primer`, `seed` (the primer itself) and `hits` (both-strand
#'   genome-wide alignment count; more than one = a potential off-target), plus
#'   `engine` (`"Rbowtie"`) and an `engine_stamp` attribute.  Hits are summed
#'   over every index in `index`, so a per-chromosome set gives the whole-genome
#'   count in one row per primer.
#' @examples
#' \dontrun{
#' idx <- marker_bowtie_index(marker_genome("genome.fa"))
#' specificity_bowtie(c(AS = "ACGTACGTACGTACGT"), idx)
#' }
#' @family specificity
#' @export
specificity_bowtie <- function(primers, index, max_mismatches = 0L, nthreads = 4L) {
  nthreads <- check_worker_count(nthreads, "nthreads")
  rlang::check_installed("Rbowtie", reason = gettext("for the bowtie genome scan"))
  max_mismatches <- check_scan_radius(max_mismatches, "max_mismatches")
  if (max_mismatches > 3L) {
    abort_msg("bowtie1 allows at most 3 mismatches (`max_mismatches` is {max_mismatches}).")
  }
  # `as.character()` drops the set's names but keeps its length, which is the
  # only thing a set is told by: two or more basenames are summed
  index <- as.character(index)
  hits <- if (length(index) > 1L) {
    rows <- lapply(unname(index), function(one) {
      bowtie_primer_hits(primers, one, max_mismatches, nthreads)
    })
    Reduce(`+`, rows)
  } else {
    bowtie_primer_hits(primers, index, max_mismatches, nthreads)
  }
  marker_with_engine(
    tibble::tibble(primer = names(primers), seed = as.character(primers), hits = hits),
    "Rbowtie"
  )
}

# one index queried, as integer hits in `primers` order
bowtie_primer_hits <- function(primers, index, max_mismatches, nthreads) {
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
  hits
}
