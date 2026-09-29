# =============================================================================
#  kasp_io.R  --  input readers (VCF variants, arbitrary genome FASTA)
# =============================================================================
#  Kept separate from the design code so the design engine stays pure.
#
#  Big-file strategy:
#    * VCF   -> `vroom` (memory-mapped + multithreaded), reading only the four
#               columns we need; the header line is parsed once to learn the
#               (sample-dependent) column count.
#    * FASTA -> `Biostrings::readDNAStringSet`; for genome-scale random access
#               one should index (`Rsamtools::indexFa`) - noted in the docs.
#
#  Variants are returned with REF/ALT as STRINGS, so SNP, MNP and InDel all
#  survive; multi-allelic ALT fields are split into one row per allele and
#  symbolic alleles (`<DEL>`, `*`, `.`) are dropped.
# =============================================================================

#' Normalise an allele string to plain `ACGTN`
#'
#' Upper-cases the string and folds any IUPAC ambiguity code (R, Y, S, W, K,
#' M, B, D, H, V) - or anything else non-ACGTN - to `N`, so downstream design
#' can treat "unknown base" uniformly instead of dropping the whole record.
#' @keywords internal
#' @export
normalise_allele <- function(x) {
  x <- toupper(x)
  x <- chartr("RYSWKMBDHV", "NNNNNNNNNN", x)   # IUPAC ambiguity -> N
  gsub("[^ACGTN]", "N", x)
}

#' Read variants from a VCF/BCF file (streaming, memory-mapped)
#'
#' @param path VCF path; plain (`vcf`) or bgzipped (`vcf.gz`) are both read by
#'   `vroom`, so no external decompression is needed (pair with `tabix` for a
#'   random-access index on very large files).
#' @param chrom,start,end optional region filter (chromosome is normalised with
#'   [ath_chrom()], so `"4"` and `"Chr4"` both work).
#' @export
read_variants_vcf <- function(path, chrom = NULL, start = NULL, end = NULL) {
  rlang::check_installed("vroom", reason = "to read large VCF files")
  if (!fs::file_exists(path)) rlang::abort(glue::glue("VCF not found: {path}"))

  hdr <- readr::read_lines(path, n_max = 1000L)
  hdr <- hdr[startsWith(hdr, "#CHROM")][1L]
  if (is.na(hdr)) rlang::abort(glue::glue("{path} has no '#CHROM' header line."))
  nms <- strsplit(sub("^#", "", hdr), "\t", fixed = TRUE)[[1L]]

  ct <- vroom::cols(.default = vroom::col_skip(),
                    CHROM = "c", POS = "i", REF = "c", ALT = "c")
  v <- vroom::vroom(path, comment = "#", col_names = nms, col_types = ct,
                    col_select = c("CHROM", "POS", "REF", "ALT"), progress = FALSE)

  v <- v |>
    dplyr::transmute(chrom = ath_chrom(.data$CHROM), pos = as.integer(.data$POS),
                     ref = .data$REF, alt = .data$ALT)
  if (!is.null(chrom)) v <- dplyr::filter(v, .data$chrom == ath_chrom(chrom))
  if (!is.null(start)) v <- dplyr::filter(v, .data$pos >= start)
  if (!is.null(end))   v <- dplyr::filter(v, .data$pos <= end)

  v |>
    tidyr::separate_longer_delim("alt", delim = ",") |>          # multi-allelic -> rows
    # drop symbolic/absent alleles (`<DEL>`, `*`, `.`) before normalisation,
    # so they are not mangled into a fake `N` allele
    dplyr::filter(!stringr::str_detect(.data$ref, "^<|[*.]"),
                  !stringr::str_detect(.data$alt, "^<|[*.]")) |>
    dplyr::mutate(ref = normalise_allele(.data$ref),
                  alt = normalise_allele(.data$alt)) |>
    dplyr::filter(nzchar(.data$ref), nzchar(.data$alt))
}

#' Read an arbitrary genome FASTA, optionally normalising sequence names
#'
#' @param path FASTA path (gzipped is fine).
#' @param chrom_map optional named character vector `defline -> name`
#'   (e.g. `ATH_CHROM_MAP`); when given, only mapped sequences are kept.
#' @param seqids optional character vector of sequence names to keep after
#'   reading (a cheap subset when only a few chromosomes are needed).
#' @export
read_genome_fasta <- function(path, chrom_map = NULL, seqids = NULL) {
  if (!fs::file_exists(path)) rlang::abort(glue::glue("FASTA not found: {path}"))
  g <- Biostrings::readDNAStringSet(path)
  nm <- sub("\\s.*$", "", names(g))
  if (!is.null(chrom_map)) {
    keep <- nm %in% names(chrom_map)
    g <- g[keep]
    names(g) <- unname(chrom_map[nm[keep]])
  } else {
    names(g) <- nm
  }
  if (!is.null(seqids)) g <- g[intersect(seqids, names(g))]
  g
}
