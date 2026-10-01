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
  x <- chartr("RYSWKMBDHV", "NNNNNNNNNN", x) # IUPAC ambiguity -> N
  gsub("[^ACGTN]", "N", x)
}

#' Read variants from a VCF/BCF file (streaming, memory-mapped)
#'
#' @param path VCF path; plain (`vcf`) or bgzipped (`vcf.gz`) are both read by
#'   `vroom`, so no external decompression is needed (pair with `tabix` for a
#'   random-access index on very large files).
#' @param chrom,start,end optional region filter (chromosome is normalised with
#'   `ath_chrom()`, so `"4"` and `"Chr4"` both work).
#' @param backend `"vroom"` (default) streams the four columns we need and is
#'   best for large files; `"VariantAnnotation"` parses with the Bioconductor
#'   canonical reader (`readVcf`) and is the safer choice for complex records
#'   (INFO/GT, phased, breakends, missing ALT).  When a region and a bgzipped +
#'   tabix-indexed file are given, `readVcf` streams just that region
#'   (`ScanVcfParam(which=)`); both backends return the same tibble.
#' @family KASP input
#' @export
read_variants_vcf <- function(path, chrom = NULL, start = NULL, end = NULL,
                              backend = c("vroom", "VariantAnnotation")) {
  backend <- match.arg(backend)
  if (!fs::file_exists(path)) rlang::abort(glue::glue(gettext("VCF not found: {path}")))

  v <- switch(backend,
    vroom = read_variants_vcf_vroom(path),
    VariantAnnotation = read_variants_vcf_va(path, chrom, start, end)
  )

  if (!is.null(chrom)) v <- dplyr::filter(v, .data$chrom == ath_chrom(chrom))
  if (!is.null(start)) v <- dplyr::filter(v, .data$pos >= start)
  if (!is.null(end)) v <- dplyr::filter(v, .data$pos <= end)

  v |>
    # drop symbolic/absent alleles (`<DEL>`, `*`, `.`) before normalisation,
    # so they are not mangled into a fake `N` allele
    dplyr::filter(
      !stringr::str_detect(.data$ref, "^<|[*.]"),
      !stringr::str_detect(.data$alt, "^<|[*.]")
    ) |>
    dplyr::mutate(
      ref = normalise_allele(.data$ref),
      alt = normalise_allele(.data$alt)
    ) |>
    dplyr::filter(nzchar(.data$ref), nzchar(.data$alt))
}

#' Streaming VCF reader (`vroom`): CHROM/POS/REF/ALT, multi-allelic split
#' @keywords internal
#' @noRd
read_variants_vcf_vroom <- function(path) {
  rlang::check_installed("vroom", reason = gettext("to read large VCF files"))
  hdr <- readr::read_lines(path, n_max = 1000L)
  hdr <- hdr[startsWith(hdr, "#CHROM")][1L]
  if (is.na(hdr)) rlang::abort(glue::glue(gettext("{path} has no '#CHROM' header line.")))
  nms <- strsplit(sub("^#", "", hdr), "\t", fixed = TRUE)[[1L]]

  ct <- vroom::cols(
    .default = vroom::col_skip(),
    CHROM = "c", POS = "i", REF = "c", ALT = "c"
  )
  vroom::vroom(path,
    comment = "#", col_names = nms, col_types = ct,
    col_select = c("CHROM", "POS", "REF", "ALT"), progress = FALSE
  ) |>
    dplyr::transmute(
      chrom = ath_chrom(.data$CHROM), pos = as.integer(.data$POS),
      ref = .data$REF, alt = .data$ALT
    ) |>
    tidyr::separate_longer_delim("alt", delim = ",") # one row per ALT allele
}

#' Canonical VCF reader (`VariantAnnotation::readVcf`)
#'
#' With a region and an indexed (bgzip+tabix) file, only that region is read
#' (`ScanVcfParam(which=)`); the header's contig names are matched to `chrom` with
#' the same `ath_chrom()` normalisation used everywhere else.
#' @keywords internal
#' @noRd
read_variants_vcf_va <- function(path, chrom = NULL, start = NULL, end = NULL) {
  rlang::check_installed("VariantAnnotation", reason = gettext("for the readVcf backend"))
  param <- VariantAnnotation::ScanVcfParam()
  if (!is.null(chrom) && !is.null(start) && !is.null(end)) {
    seqs <- tryCatch(
      as.character(GenomeInfoDb::seqnames(VariantAnnotation::scanVcfHeader(path))),
      error = function(e) character(0)
    )
    hit <- seqs[ath_chrom(seqs) == ath_chrom(chrom)][1L]
    if (length(hit) && !is.na(hit)) {
      param <- VariantAnnotation::ScanVcfParam(
        which = GenomicRanges::GRanges(hit, IRanges::IRanges(start, end))
      )
    }
  }
  vcf <- VariantAnnotation::readVcf(path, param = param)
  rng <- SummarizedExperiment::rowRanges(vcf)
  tibble::tibble(
    chrom = ath_chrom(as.character(GenomicRanges::seqnames(rng))),
    pos   = as.integer(GenomicRanges::start(rng)),
    ref   = as.character(VariantAnnotation::ref(vcf)),
    alt   = purrr::map(VariantAnnotation::alt(vcf), as.character)
  ) |>
    tidyr::unnest_longer("alt") # one row per ALT allele
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
  rlang::check_installed("Rsamtools", reason = gettext("for indexed FASTA random access"))
  if (!fs::file_exists(path)) rlang::abort(glue::glue(gettext("FASTA not found: {path}")))
  plain <- if (grepl("\\.gz$", path)) sub("\\.gz$", "", path) else path
  if (!fs::file_exists(plain)) {
    cli::cli_alert_info(gettext("Decompressing {.file {fs::path_file(path)}} for indexing (one-off)."))
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
  if (!length(positions)) {
    return(s)
  }
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
  f <- system.file("extdata", "example_targets.csv", package = "tidykasp")
  if (!nzchar(f)) rlang::abort(gettext("Bundled example data not found."))
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
  if (!fs::file_exists(path)) rlang::abort(glue::glue(gettext("FASTA not found: {path}")))

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
    fai <- readr::read_tsv(paste0(plain, ".fai"),
      show_col_types = FALSE,
      col_names = c("name", "length", "offset", "linebases", "linebytes")
    )
    tok <- sub("\\s.*$", "", fai$name)
    mapped <- if (!is.null(chrom_map)) unname(chrom_map[tok]) else tok
    keep <- (tok %in% seqids) | (!is.na(mapped) & mapped %in% seqids)
    sel <- fai[keep, , drop = FALSE]
    if (nrow(sel) == 0L) {
      return(Biostrings::DNAStringSet())
    }
    gr <- GenomicRanges::GRanges(sel$name, IRanges::IRanges(1L, sel$length))
    g <- Rsamtools::scanFa(Rsamtools::FaFile(plain), gr)
    return(normalise(g))
  }
  normalise(Biostrings::readDNAStringSet(path))
}
