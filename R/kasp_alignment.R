# =============================================================================
#  kasp_alignment.R  --  a gene's transcripts -> an alignment -> haplotypes
# =============================================================================
#  A multi-sequence FASTA (e.g. an NCBI *seqdump* of one gene across accessions
#  or cultivars) is a ready-made **haplotype set**: align the transcripts, take
#  one as the reference, and every other sequence is a haplotype carrying its own
#  variants.  That is exactly the `haplotypes=` input [assemble_kasp_assay()] and
#  [kasp_context()] understand - so a public gene family becomes an *offline*
#  design fixture, with no genome download.
#
#  The alignment is a DECIPHER multiple alignment (already a dependency).  Variants
#  are called against the reference row and mapped to reference coordinates, with
#  substitutions and indels both reported (VCF-anchored: an indel's `pos` is the
#  base immediately before it, `ref`/`alt` the anchored allele strings).
# =============================================================================

#' Read a transcript FASTA into a haplotype set
#'
#' Reads a FASTA of one gene's transcripts (e.g. a NCBI *seqdump* across
#' accessions or cultivars) and returns one sequence per accession.  A seqdump
#' often splits an accession into coordinate-range records (`:257-1327` and
#' `:1395-1452`); with `merge_accessions = TRUE` those parts are joined in file
#' order to reconstruct the transcript (the region between the ranges is absent
#' from the dump, so the join is taken as given).
#' @param path a FASTA file.
#' @param merge_accessions join records that share an accession (the header up to
#'   `:`)?  `FALSE` keeps every record separate.
#' @return a named `DNAStringSet`, one sequence per accession (or per record).
#' @examples
#' f <- system.file("extdata", "sat2_transcripts.fa", package = "tidykasp")
#' read_transcripts(f)
#' @family KASP haplotypes
#' @export
read_transcripts <- function(path, merge_accessions = TRUE) {
  if (!fs::file_exists(path)) {
    cli::cli_abort(gettext("FASTA not found: {path}"))
  }
  s <- Biostrings::readDNAStringSet(path)
  if (!merge_accessions) {
    return(s)
  }
  acc <- sub(":.*$", "", names(s))
  parts <- split(as.character(s), factor(acc, levels = unique(acc)))
  Biostrings::DNAStringSet(vapply(parts, paste0, collapse = "", character(1)))
}

#' Call one haplotype's variants against the (gapped) reference row
#'
#' Walks the alignment columns, tracking the reference position, and emits a
#' substitution, an insertion or a deletion row per event (reference coordinates).
#' @keywords internal
#' @noRd
call_haplotype_variants <- function(ref_aln, seq_aln, ref_seq, chrom, hap) {
  r <- strsplit(ref_aln, "", fixed = TRUE)[[1L]]
  q <- strsplit(seq_aln, "", fixed = TRUE)[[1L]]
  n <- length(r)
  rp <- 0L
  i <- 1L
  pos <- integer(0)
  ref <- character(0)
  alt <- character(0)
  while (i <= n) {
    rc <- r[i]
    qc <- q[i]
    if (rc != "-" && qc != "-") {
      rp <- rp + 1L
      if (rc != qc) {
        pos <- c(pos, rp)
        ref <- c(ref, rc)
        alt <- c(alt, qc)
      }
      i <- i + 1L
    } else if (rc == "-" && qc != "-") {
      ins <- character(0)
      while (i <= n && r[i] == "-" && q[i] != "-") {
        ins <- c(ins, q[i])
        i <- i + 1L
      }
      if (rp < 1L) {
        cli::cli_abort(gettext("Haplotype {hap} has an insertion before the reference start; cannot anchor it."))
      }
      a <- substr(ref_seq, rp, rp)
      pos <- c(pos, rp)
      ref <- c(ref, a)
      alt <- c(alt, paste0(a, paste0(ins, collapse = "")))
    } else if (rc != "-" && qc == "-") {
      del <- character(0)
      while (i <= n && r[i] != "-" && q[i] == "-") {
        del <- c(del, r[i])
        rp <- rp + 1L
        i <- i + 1L
      }
      anchor <- rp - length(del)
      if (anchor < 1L) {
        cli::cli_abort(gettext("Haplotype {hap} has a deletion at the reference start; cannot anchor it."))
      }
      a <- substr(ref_seq, anchor, anchor)
      pos <- c(pos, anchor)
      ref <- c(ref, paste0(a, paste0(del, collapse = "")))
      alt <- c(alt, a)
    } else {
      i <- i + 1L
    }
  }
  tibble::tibble(haplotype = hap, chrom = chrom, pos = pos, ref = ref, alt = alt)
}

#' Align transcripts and call haplotypes against a reference
#'
#' Aligns a set of (near-identical) transcripts with [DECIPHER::AlignSeqs()], takes
#' one as the reference, and returns the reference as a one-sequence genome plus a
#' tidy **haplotype table** (`haplotype`, `freq`, `chrom`, `pos`, `ref`, `alt`) in
#' reference coordinates - the `haplotypes=` input of [assemble_kasp_assay()] /
#' [kasp_context()].  Substitutions and indels are both reported, so a real gene
#' family becomes a designable haplotype set with no genome download.
#' @param seqs a `DNAStringSet` (or a character vector) of transcripts.
#' @param reference index or name of the reference sequence (default: the first).
#' @param chrom the chromosome name to give the reference (the haplotype table's
#'   `chrom`, and the genome's name).
#' @param freq optional haplotype frequencies (a named vector, a single number, or
#'   `NULL` for equal frequencies) - the `pop_freq` weights of the design.
#' @param ... passed to [DECIPHER::AlignSeqs()].
#' @return a list: `reference` (a `DNAString`), `genome` (a one-sequence
#'   `DNAStringSet` named `chrom`), `haplotypes` (the tidy table), `alignment` (a
#'   `DNAMultipleAlignment`) and `reference_name`.
#' @examples
#' seqs <- c(
#'   ref = "ACGTACGTACGTACGTACGT", hap1 = "ACGTACGAACGTACGTACGT",
#'   hap2 = "ACGTACGTACGTACGTACGA"
#' )
#' hp <- haplotypes_from_alignment(seqs, chrom = "GENE")
#' hp$haplotypes
#' @family KASP haplotypes
#' @export
haplotypes_from_alignment <- function(seqs, reference = 1L, chrom = "SAT2",
                                      freq = NULL, ...) {
  if (is.character(seqs)) seqs <- Biostrings::DNAStringSet(seqs)
  if (!inherits(seqs, "DNAStringSet")) {
    cli::cli_abort(gettext("`seqs` must be a DNAStringSet or a character vector."))
  }
  if (length(seqs) < 2L) {
    cli::cli_abort(gettext("Need at least two sequences to call haplotypes."))
  }
  rlang::check_installed("DECIPHER", reason = gettext("to align the haplotypes"))
  aln <- DECIPHER::AlignSeqs(seqs, verbose = FALSE, ...)
  a <- stats::setNames(as.character(aln), names(aln))
  ref_name <- if (is.numeric(reference)) {
    if (reference < 1L || reference > length(a)) {
      cli::cli_abort(gettext("`reference` is out of range (1..{length(a)})."))
    }
    names(a)[[reference]]
  } else {
    if (!reference %in% names(a)) {
      cli::cli_abort(gettext("`reference` = {reference} is not a sequence name."))
    }
    reference
  }
  ref_seq <- gsub("-", "", a[[ref_name]], fixed = TRUE)
  others <- setdiff(names(a), ref_name)
  variants <- purrr::map(others, function(h) {
    call_haplotype_variants(a[[ref_name]], a[[h]], ref_seq, chrom, h)
  }) |>
    purrr::list_rbind()
  if (is.null(freq)) {
    f <- rep(1 / length(a), length(a))
  } else if (length(freq) == 1L) {
    f <- rep(freq, length(a))
  } else {
    f <- unname(freq[names(a)])
  }
  hap_tbl <- tibble::tibble(
    haplotype = variants$haplotype,
    freq = f[match(variants$haplotype, names(a))],
    chrom = variants$chrom, pos = variants$pos, ref = variants$ref, alt = variants$alt
  )
  list(
    reference = Biostrings::DNAString(ref_seq),
    genome = Biostrings::DNAStringSet(stats::setNames(ref_seq, chrom)),
    haplotypes = hap_tbl,
    alignment = Biostrings::DNAMultipleAlignment(aln),
    reference_name = ref_name
  )
}
