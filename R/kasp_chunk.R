# =============================================================================
#  kasp_chunk.R  --  a local chunk of a haplotype set, in every representation
# =============================================================================
#  A target sits in a *local chunk* - a window of a gene - and the same chunk is
#  naturally held in three interchangeable shapes:
#
#    alignment   a DNAMultipleAlignment (where the haplotypes come from)
#    variants    a VCF-like table  chrom, pos, ref, alt [, pop_freq]
#    haplotypes  a weighted haplotype table  haplotype, freq, chrom, pos, ref, alt
#
#  `kasp_chunk()` takes whichever you have plus a window and returns a chunk you
#  can convert cheaply with `chunk_variants()` / `chunk_haplotypes()` /
#  `chunk_alignment()`, and weight with `chunk_weight()`.  Coordinates in a chunk
#  are *window-local* (1..n), and `chunk_genome()` gives a one-sequence genome, so
#  a chunk is a drop-in for assemble_kasp_assay()/kasp_context().
#
#  Honest asymmetry: an alignment yields both variants and haplotypes; a compound
#  haplotype table yields its variants; but a plain variants table does NOT yield
#  haplotypes - that needs phasing, which the table does not carry.
# =============================================================================

#' A local chunk of a haplotype set, in every representation
#'
#' Accepts a [haplotypes_from_alignment()] result, a `DNAMultipleAlignment`, a
#' VCF-like variants table, or a weighted haplotype table, and (optionally) keeps
#' only a window `from`..`to` of the reference.  A chunk is *window-local*: its
#' positions start at 1, and `chunk_genome()` is a one-sequence `DNAStringSet` you
#' can hand to [assemble_kasp_assay()].
#' @param x the source (a chunk itself, a [haplotypes_from_alignment()] list, an
#'   alignment, or a table).
#' @param ... passed to the method (e.g. `from`/`to`/`chrom`, or `reference`).
#' @return a `kasp_chunk`.
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "sat2_transcripts.fa", package = "tidykasp")
#' hp <- haplotypes_from_alignment(read_transcripts(f), chrom = "SAT2")
#' ch <- kasp_chunk(hp, from = 350, to = 550) # a 200 bp local chunk
#' chunk_variants(ch)
#' }
#' @family KASP haplotypes
#' @export
kasp_chunk <- function(x, ...) UseMethod("kasp_chunk")

#' @export
kasp_chunk.kasp_chunk <- function(x, ...) x

#' @export
kasp_chunk.list <- function(x, from = NULL, to = NULL, chrom = NULL, ...) {
  if (is.null(x$reference) || is.null(x$alignment) || is.null(x$haplotypes)) {
    cli::cli_abort(gettext("A list source must come from haplotypes_from_alignment()."))
  }
  chunk_slice(x, from, to, chrom)
}

#' @export
kasp_chunk.data.frame <- function(x, reference = NULL, chrom = "chunk", ...) {
  if (is.null(reference)) {
    cli::cli_abort(gettext("A chunk from a table needs a `reference` sequence."))
  }
  if (all(c("haplotype", "freq") %in% names(x))) {
    x <- tibble::as_tibble(x)
    if (!"chrom" %in% names(x)) x$chrom <- chrom
    haps <- normalise_haplotypes(x)
    kasp_chunk_new(reference, variants_from_haplotypes(haps), haps, NULL, chrom)
  } else if (all(c("pos", "ref", "alt") %in% names(x))) {
    kasp_chunk_new(reference, as_variant_table(x, chrom), NULL, NULL, chrom)
  } else {
    cli::cli_abort(gettext("A table source needs `pos`/`ref`/`alt`, or `haplotype`/`freq`."))
  }
}

#' @export
kasp_chunk.DNAMultipleAlignment <- function(x, reference = NULL, chrom = "chunk", ...) {
  a <- as.character(x)
  if (is.null(reference)) reference <- names(a)[1L]
  if (!reference %in% names(a)) {
    cli::cli_abort(gettext("`reference` must be one of the alignment's sequence names."))
  }
  ref <- gsub("-", "", a[[reference]], fixed = TRUE)
  hap <- setdiff(names(a), reference)
  vars <- purrr::map(hap, function(h) {
    call_haplotype_variants(a[[reference]], a[[h]], ref, chrom, h)
  }) |>
    purrr::list_rbind()
  # every sequence is one haplotype: equal frequencies unless weighted later
  vars$freq <- 1 / length(a)
  kasp_chunk_new(ref, variants_from_haplotypes(vars), vars, x, chrom)
}

#' Build a chunk record (internal)
#' @keywords internal
#' @noRd
kasp_chunk_new <- function(reference, variants, haplotypes, alignment, chrom) {
  ref <- as.character(reference)
  structure(
    list(
      chrom = chrom,
      reference = Biostrings::DNAString(ref),
      reference_length = nchar(ref),
      genome = Biostrings::DNAStringSet(stats::setNames(ref, chrom)),
      variants = variants,
      haplotypes = haplotypes,
      alignment = alignment
    ),
    class = "kasp_chunk"
  )
}

#' Coerce a table to the variants schema (internal)
#' @keywords internal
#' @noRd
as_variant_table <- function(x, chrom = "chunk") {
  x <- tibble::as_tibble(x)
  out <- tibble::tibble(
    chrom = if ("chrom" %in% names(x)) as.character(x$chrom) else chrom,
    pos = as.integer(x$pos), ref = as.character(x$ref), alt = as.character(x$alt),
    pop_freq = if ("pop_freq" %in% names(x)) as.numeric(x$pop_freq) else 1,
    consequence = if ("consequence" %in% names(x)) as.character(x$consequence) else NA_character_
  )
  out
}

#' Keep only a window of a haplotypes_from_alignment() result (internal)
#' @keywords internal
#' @noRd
chunk_slice <- function(hp, from, to, chrom) {
  ref <- as.character(hp$reference)
  if (is.null(from)) from <- 1L
  if (is.null(to)) to <- nchar(ref)
  if (from < 1L || to > nchar(ref) || to < from) {
    cli::cli_abort(gettext("The window {from}..{to} is outside the reference (1..{nchar(ref)})."))
  }
  shift <- function(d) {
    if (is.null(d)) {
      return(NULL)
    }
    d <- d[d$pos >= from & d$pos <= to, , drop = FALSE]
    d$pos <- d$pos - from + 1L
    d
  }
  a <- as.character(hp$alignment)
  m <- do.call(rbind, strsplit(a, "", fixed = TRUE))
  rownames(m) <- names(a)
  cols_of_ref <- which(m[hp$reference_name, ] != "-")
  sub <- m[, cols_of_ref[from:to], drop = FALSE]
  aln <- Biostrings::DNAMultipleAlignment(
    Biostrings::DNAStringSet(apply(sub, 1L, paste0, collapse = ""))
  )
  kasp_chunk_new(
    substr(ref, from, to), shift(hp$variants), shift(hp$haplotypes), aln,
    chrom %||% hp$chrom %||% "chunk"
  )
}

#' @export
print.kasp_chunk <- function(x, ...) {
  nh <- if (is.null(x$haplotypes)) {
    "none (unphased variants)"
  } else {
    length(unique(x$haplotypes$haplotype))
  }
  na <- if (is.null(x$alignment)) "no" else "yes"
  cat(glue::glue("<kasp_chunk> {x$chrom} 1..{x$reference_length}"), "\n", sep = "")
  cat(glue::glue("  {nrow(x$variants)} variant(s); haplotypes: {nh}; alignment: {na}"), "\n", sep = "")
  invisible(x)
}

#' A chunk's parts
#'
#' The three representations (and the reference/genome) of a [kasp_chunk()].
#' @param chunk a `kasp_chunk`.
#' @return `chunk_reference()` a `DNAString`; `chunk_genome()` a one-sequence
#'   `DNAStringSet`; `chunk_variants()`/`chunk_haplotypes()` tidy tables (the
#'   latter `NULL` when the source was unphased); `chunk_alignment()` a
#'   `DNAMultipleAlignment` (or `NULL`).
#' @examples
#' \dontrun{
#' chunk_variants(kasp_chunk(hp, from = 350, to = 550))
#' }
#' @family KASP haplotypes
#' @export
chunk_reference <- function(chunk) {
  check_chunk(chunk)
  chunk$reference
}

#' @rdname chunk_reference
#' @export
chunk_genome <- function(chunk) {
  check_chunk(chunk)
  chunk$genome
}

#' @rdname chunk_reference
#' @export
chunk_variants <- function(chunk) {
  check_chunk(chunk)
  chunk$variants
}

#' @rdname chunk_reference
#' @export
chunk_haplotypes <- function(chunk) {
  check_chunk(chunk)
  chunk$haplotypes
}

#' @rdname chunk_reference
#' @export
chunk_alignment <- function(chunk) {
  check_chunk(chunk)
  chunk$alignment
}

#' @keywords internal
#' @noRd
check_chunk <- function(chunk) {
  if (!inherits(chunk, "kasp_chunk")) {
    cli::cli_abort(gettext("`chunk` must be a kasp_chunk (see kasp_chunk())."))
  }
}

#' Weight a chunk's variants by minor-allele frequencies
#'
#' Overrides the derived `pop_freq` of the matching variants (MAF-style), so a
#' chunk's weights can come from a real panel rather than the haplotype counts.
#' @param chunk a `kasp_chunk`.
#' @param maf a per-variant table (`pos`, `ref`, `alt`, and `maf`/`freq`) in the
#'   chunk's window-local coordinates.
#' @return the chunk with weighted `variants`.
#' @examples
#' \dontrun{
#' chunk_weight(ch, tibble::tibble(pos = 81, ref = "C", alt = "T", maf = 0.02))
#' }
#' @family KASP haplotypes
#' @export
chunk_weight <- function(chunk, maf) {
  check_chunk(chunk)
  chunk$variants <- if (is.null(chunk$haplotypes)) {
    apply_maf(chunk$variants, maf)
  } else {
    variants_from_haplotypes(chunk$haplotypes, maf = maf)
  }
  chunk
}
