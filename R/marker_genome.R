# =============================================================================
#  marker_genome.R  --  an indexed genome source (no whole-genome DNA in RAM)
# =============================================================================
#  The design window (a few hundred bp around a variant) is all the *designer*
#  needs; the genome-wide verification steps are the only ones that must see the
#  whole genome, and those have their own indexed path (see marker_scan.R: BLAST).
#
#  A `marker_genome` is an Rsamtools-indexed FASTA you read regions from on demand
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
#' [genome_seq()] accept both a `marker_genome` and a plain `DNAStringSet`, so the
#' design code does not care which one it is given.
#'
#' A `marker_genome` can also be built from a Bioconductor `BSgenome` object
#' (`bsgenome =`), in which case regions are read with `BSgenome::getSeq()` and no
#' FASTA needs to live on disk.
#'
#' @param path FASTA path (gzipped is decompressed once, see [fasta_index()]).
#'   Ignored when `bsgenome` is given.
#' @param chrom_map optional named character vector `defline -> name`
#'   (e.g. `ATH_CHROM_MAP`); when given, only mapped sequences are kept.
#' @param blast_db optional prebuilt BLAST database path for the indexed
#'   genome-wide scan (see [marker_blast_index()] / [specificity_blast()]).
#' @param fm_index optional prebuilt FM-index basename for the indexed
#'   genome-wide scan (see [marker_fm_index()] / [specificity_fm()]).
#' @param bsgenome a `BSgenome` object (or a package name) to read from instead
#'   of a FASTA (e.g. `"BSgenome.Athaliana.TAIR.TAIR9"`); `path` is then unused.
#' @return a `marker_genome` object.
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' marker_genome(fa, chrom_map = c("1" = "Chr1"))
#' @family KASP input
#' @export
marker_genome <- function(path = NULL, chrom_map = NULL, blast_db = NULL,
                          bsgenome = NULL, fm_index = NULL) {
  if (!is.null(bsgenome)) {
    return(bsgenome_kasp_genome(bsgenome, chrom_map, blast_db, fm_index))
  }
  plain <- fasta_index(path)
  fa <- Rsamtools::FaFile(plain)
  si <- GenomeInfoDb::seqinfo(fa) # names + lengths from the .fai
  raw <- as.character(GenomeInfoDb::seqnames(si))
  len <- rlang::set_names(as.integer(GenomeInfoDb::seqlengths(si)), raw)
  norm <- if (!is.null(chrom_map)) unname(chrom_map[raw]) else raw
  keep <- if (!is.null(chrom_map)) !is.na(norm) else rep(TRUE, length(raw))
  seqnames <- rlang::set_names(raw[keep], (if (!is.null(chrom_map)) norm else raw)[keep])
  structure(
    list(
      path = plain, fa = fa, seqnames = seqnames,
      lengths = rlang::set_names(unname(len[raw[keep]]), names(seqnames)),
      blast_db = blast_db, fm_index = fm_index
    ),
    class = "marker_genome"
  )
}

#' Build a `marker_genome` from a `BSgenome` object
#'
#' `getSeq()` on a `BSgenome` is the Bioconductor-native way to read reference
#' sequence, so a `marker_genome` can be backed by one instead of a FASTA file.
#' @keywords internal
#' @noRd
bsgenome_kasp_genome <- function(bsgenome, chrom_map, blast_db, fm_index = NULL) {
  load_namespace_quietly("BSgenome", reason = gettext("for a BSgenome genome source"))
  if (is.character(bsgenome)) {
    load_namespace_quietly(bsgenome, reason = gettext("for the BSgenome genome source"))
    bsgenome <- getExportedValue(bsgenome, bsgenome)
  }
  si <- GenomeInfoDb::seqinfo(bsgenome) # names + lengths, no DNA
  raw <- as.character(GenomeInfoDb::seqnames(si))
  len <- rlang::set_names(as.integer(GenomeInfoDb::seqlengths(si)), raw)
  norm <- if (!is.null(chrom_map)) unname(chrom_map[raw]) else raw
  keep <- if (!is.null(chrom_map)) !is.na(norm) else rep(TRUE, length(raw))
  seqnames <- rlang::set_names(raw[keep], (if (!is.null(chrom_map)) norm else raw)[keep])
  structure(
    list(
      path = NULL, fa = NULL, bsgenome = bsgenome, seqnames = seqnames,
      lengths = rlang::set_names(unname(len[raw[keep]]), names(seqnames)),
      blast_db = blast_db, fm_index = fm_index
    ),
    class = "marker_genome"
  )
}

#' Print a KASP genome source
#'
#' Compact summary: FASTA file, sequences and whether a BLAST db is attached.
#'
#' @param x a `marker_genome`.
#' @param ... unused.
#' @return `x`, invisibly (called for its side effect of printing a summary).
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' print(marker_genome(fa, chrom_map = c("1" = "Chr1")))
#' @method print marker_genome
#' @family KASP input
#' @export
print.marker_genome <- function(x, ...) {
  source <- if (!is.null(x$bsgenome)) paste0("<BSgenome> ", class(x$bsgenome)[1L]) else fs::path_file(x$path)
  cat("<marker_genome>", source, "\n")
  cat(
    "  sequences:", length(x$seqnames),
    paste0("(", paste(utils::head(names(x$seqnames), 6L), collapse = ", "), ")"), "\n"
  )
  cat("  BLAST db:", if (is.null(x$blast_db)) "none (see marker_blast_index())" else x$blast_db, "\n")
  cat("  FM index:", if (is.null(x$fm_index)) "none (see marker_fm_index())" else x$fm_index, "\n")
  invisible(x)
}

#' Sequence names of a genome (indexed source or DNAStringSet)
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @return a character vector of sequence names.
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' genome_chroms(marker_genome(fa, chrom_map = c("1" = "Chr1")))
#' @family KASP input
#' @export
genome_chroms <- function(genome) {
  if (inherits(genome, "marker_genome")) names(genome$seqnames) else names(genome)
}

#' Length of one sequence of a genome
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param chrom sequence name.
#' @return the sequence length as a single integer.
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' genome_length(marker_genome(fa, chrom_map = c("1" = "Chr1")), "Chr1")
#' @family KASP input
#' @export
genome_length <- function(genome, chrom) {
  if (!inherits(genome, "marker_genome")) {
    return(length(genome[[chrom]]))
  }
  len <- genome$lengths[[chrom]]
  if (is.null(len)) cli::cli_abort(gettext("Unknown sequence: {chrom}"))
  len
}

#' Read a region as a DNAString, from an indexed source or a DNAStringSet
#'
#' This is the single accessor the design code uses, so an indexed genome and an
#' in-memory one behave identically.
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param chrom sequence name (as returned by [genome_chroms()]).
#' @param start,end 1-based inclusive coordinates; default to the whole sequence.
#' @return a `DNAString` for the requested interval.
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGTACGT"), fa)
#' genome_seq(marker_genome(fa, chrom_map = c("1" = "Chr1")), "Chr1", 2, 6)
#' @family KASP input
#' @export
genome_seq <- function(genome, chrom, start = 1L, end = genome_length(genome, chrom)) {
  start <- as.integer(start)
  end <- as.integer(end)
  if (!inherits(genome, "marker_genome")) {
    # the whole sequence needs no copy: `subseq()` of a full range memcpys the
    # chromosome, and a genome-scale scan asks for it once per primer pair
    if (start == 1L && end == length(genome[[chrom]])) {
      return(genome[[chrom]])
    }
    return(Biostrings::subseq(genome[[chrom]], start = start, end = end))
  }
  raw <- genome$seqnames[[chrom]]
  if (is.null(raw)) cli::cli_abort(gettext("Unknown sequence: {chrom}"))
  gr <- GenomicRanges::GRanges(raw, IRanges::IRanges(start, end))
  if (!is.null(genome$bsgenome)) {
    return(BSgenome::getSeq(genome$bsgenome, gr)[[1L]])
  }
  Rsamtools::scanFa(genome$fa, gr)[[1L]]
}
