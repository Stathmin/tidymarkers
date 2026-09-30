# =============================================================================
#  kasp_genome.R  --  an indexed genome source (no whole-genome DNA in RAM)
# =============================================================================
#  The design window (a few hundred bp around a variant) is all the *designer*
#  needs; the genome-wide verification steps are the only ones that must see the
#  whole genome, and those have their own indexed path (see kasp_scan.R: BLAST).
#
#  A `kasp_genome` is an Rsamtools-indexed FASTA you read regions from on demand
#  (`genome_seq()`), so `snp_template()`, `strand_gc_hint()` and the reference
#  lookup in `pick_target_snp()` work without loading chromosomes into memory.
#  Sequence names and lengths come from the library index (`seqinfo()`), not from
#  hand-parsing the `.fai`.  Every accessor also accepts a plain in-memory
#  `DNAStringSet`, so the two representations are interchangeable.
# =============================================================================

#' Indexed genome source backed by an Rsamtools FASTA index
#'
#' Wraps an (optionally compressed) FASTA so regions can be read on demand with
#' [Rsamtools::scanFa()], instead of holding whole chromosomes as a `DNAStringSet`
#' in memory.  The accessors [genome_chroms()], [genome_length()] and
#' [genome_seq()] accept both a `kasp_genome` and a plain `DNAStringSet`, so the
#' design code does not care which one it is given.
#'
#' @param path FASTA path (gzipped is decompressed once, see [fasta_index()]).
#' @param chrom_map optional named character vector `defline -> name`
#'   (e.g. `ATH_CHROM_MAP`); when given, only mapped sequences are kept.
#' @param blast_db optional prebuilt BLAST database path for the indexed
#'   genome-wide scan (see [kasp_blast_index()] / [specificity_blast()]).
#' @return a `kasp_genome` object.
#' @family KASP input
#' @export
kasp_genome <- function(path, chrom_map = NULL, blast_db = NULL) {
  plain <- fasta_index(path)
  fa    <- Rsamtools::FaFile(plain)
  si    <- GenomeInfoDb::seqinfo(fa)                 # names + lengths from the .fai
  raw   <- as.character(GenomeInfoDb::seqnames(si))
  len   <- stats::setNames(as.integer(GenomeInfoDb::seqlengths(si)), raw)
  norm  <- if (!is.null(chrom_map)) unname(chrom_map[raw]) else raw
  keep  <- if (!is.null(chrom_map)) !is.na(norm) else rep(TRUE, length(raw))
  seqnames <- stats::setNames(raw[keep], (if (!is.null(chrom_map)) norm else raw)[keep])
  structure(list(path = plain, fa = fa, seqnames = seqnames,
                 lengths = stats::setNames(unname(len[raw[keep]]), names(seqnames)),
                 blast_db = blast_db),
            class = "kasp_genome")
}

#' Print a KASP genome source
#'
#' Compact summary: FASTA file, sequences and whether a BLAST db is attached.
#'
#' @param x a `kasp_genome`.
#' @param ... unused.
#' @method print kasp_genome
#' @family KASP input
#' @export
print.kasp_genome <- function(x, ...) {
  cat("<kasp_genome>", fs::path_file(x$path), "\n")
  cat("  sequences:", length(x$seqnames),
      paste0("(", paste(utils::head(names(x$seqnames), 6L), collapse = ", "), ")"), "\n")
  cat("  BLAST db:", if (is.null(x$blast_db)) "none (see kasp_blast_index())" else x$blast_db, "\n")
  invisible(x)
}

#' Sequence names of a genome (indexed source or DNAStringSet)
#' @param genome a `kasp_genome` or a `DNAStringSet`.
#' @family KASP input
#' @export
genome_chroms <- function(genome) {
  if (inherits(genome, "kasp_genome")) names(genome$seqnames) else names(genome)
}

#' Length of one sequence of a genome
#' @param genome a `kasp_genome` or a `DNAStringSet`.
#' @param chrom sequence name.
#' @family KASP input
#' @export
genome_length <- function(genome, chrom) {
  if (!inherits(genome, "kasp_genome")) return(length(genome[[chrom]]))
  len <- genome$lengths[[chrom]]
  if (is.null(len)) rlang::abort(glue::glue("Unknown sequence: {chrom}"))
  len
}

#' Read a region as a DNAString, from an indexed source or a DNAStringSet
#'
#' This is the single accessor the design code uses, so an indexed genome and an
#' in-memory one behave identically.
#' @param genome a `kasp_genome` or a `DNAStringSet`.
#' @param chrom sequence name (as returned by [genome_chroms()]).
#' @param start,end 1-based inclusive coordinates; default to the whole sequence.
#' @family KASP input
#' @export
genome_seq <- function(genome, chrom, start = 1L, end = genome_length(genome, chrom)) {
  start <- as.integer(start); end <- as.integer(end)
  if (!inherits(genome, "kasp_genome")) {
    return(Biostrings::subseq(genome[[chrom]], start = start, end = end))
  }
  raw <- genome$seqnames[[chrom]]
  if (is.null(raw)) rlang::abort(glue::glue("Unknown sequence: {chrom}"))
  gr <- GenomicRanges::GRanges(raw, IRanges::IRanges(start, end))
  Rsamtools::scanFa(genome$fa, gr)[[1L]]
}
