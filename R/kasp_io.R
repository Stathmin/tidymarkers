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
#' @param x character vector of allele strings.
#' @family KASP input
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
#'   `ath_chrom()`, so `"4"` and `"Chr4"` both work).
#' @family KASP input
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

#' Ensure a FASTA has a `.fai` index, decompressing `.gz` once if needed
#'
#' `Rsamtools::indexFa()` needs a plain (or bgzipped) FASTA, so a `.gz` file is
#' first materialised next to it (once, cached) before indexing.  Returns the
#' path to the plain, indexed FASTA.
#' @param path FASTA path (gzipped is decompressed once).
#' @family KASP input
#' @export
fasta_index <- function(path) {
  rlang::check_installed("Rsamtools", reason = "for indexed FASTA random access")
  if (!fs::file_exists(path)) rlang::abort(glue::glue("FASTA not found: {path}"))
  plain <- if (grepl("\\.gz$", path)) sub("\\.gz$", "", path) else path
  if (!fs::file_exists(plain)) {
    cli::cli_alert_info("Decompressing {.file {fs::path_file(path)}} for indexing (one-off).")
    Biostrings::writeXStringSet(Biostrings::readDNAStringSet(path), plain)
  }
  fai <- paste0(plain, ".fai")
  if (!fs::file_exists(fai)) Rsamtools::indexFa(plain)
  plain
}

#' Mask positions in a DNA sequence by replacing them with `N`
#'
#' Used to blank out *other* variants flanking a target so no ambiguous base
#' leaks into a designed or exported primer.  Positions outside the sequence are
#' ignored.
#' @param seq a DNAString (or string).
#' @param positions 1-based positions to mask.
#' @family KASP input
#' @export
mask_positions <- function(seq, positions) {
  s <- Biostrings::DNAString(as.character(seq))
  positions <- as.integer(positions)
  positions <- positions[!is.na(positions) & positions >= 1L & positions <= length(s)]
  if (!length(positions)) return(s)
  Biostrings::replaceAt(s, IRanges::IRanges(positions, positions), "N")
}

#' Bundled example target variants (FRI / AT4G00650)
#'
#' A small table of real 1001-Genomes variants in *FRIGIDA*, handy for the
#' `variant_id` argument of [pick_target_snp()] and for trying the pipeline
#' without pulling data from Ensembl.
#'
#' @return tibble with `id`, `chrom`, `pos`, `ref`, `alt`, `consequence`,
#'   `transition`.
#' @family KASP input
#' @export
kasp_example <- function() {
  f <- system.file("extdata", "example_targets.csv", package = "kaspexample")
  if (!nzchar(f)) rlang::abort("Bundled example data not found.")
  readr::read_csv(f, show_col_types = FALSE)
}

#' Read an arbitrary genome FASTA, optionally normalising sequence names
#'
#' @param path FASTA path (gzipped is fine).
#' @param chrom_map optional named character vector `defline -> name`
#'   (e.g. `ATH_CHROM_MAP`); when given, only mapped sequences are kept.
#' @param seqids optional character vector of sequence names to keep after
#'   reading (a cheap subset when only a few chromosomes are needed).
#' @param use_index when `TRUE` and `seqids` is given, read only those sequences
#'   through the `Rsamtools` FASTA index (random access, no whole-genome load);
#'   otherwise the whole file is read with `Biostrings`.
#' @family KASP input
#' @export
read_genome_fasta <- function(path, chrom_map = NULL, seqids = NULL, use_index = TRUE) {
  if (!fs::file_exists(path)) rlang::abort(glue::glue("FASTA not found: {path}"))

  normalise <- function(g) {
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

  if (isTRUE(use_index) && !is.null(seqids)) {
    plain <- fasta_index(path)
    fai   <- readr::read_tsv(paste0(plain, ".fai"), show_col_types = FALSE,
      col_names = c("name", "length", "offset", "linebases", "linebytes"))
    tok    <- sub("\\s.*$", "", fai$name)
    mapped <- if (!is.null(chrom_map)) unname(chrom_map[tok]) else tok
    keep   <- (tok %in% seqids) | (!is.na(mapped) & mapped %in% seqids)
    sel    <- fai[keep, , drop = FALSE]
    if (nrow(sel) == 0L) return(Biostrings::DNAStringSet())
    gr <- GenomicRanges::GRanges(sel$name, IRanges::IRanges(1L, sel$length))
    g  <- Rsamtools::scanFa(Rsamtools::FaFile(plain), gr)
    return(normalise(g))
  }
  normalise(Biostrings::readDNAStringSet(path))
}
