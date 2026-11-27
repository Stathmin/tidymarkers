# =============================================================================
#  marker_io.R  --  input readers (VCF variants, arbitrary genome FASTA)
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
#  symbolic alleles (`<DEL>`, `*`, `.`) are dropped - and counted, because a
#  GWAS hitlist silently losing its CNV records is how INDEL/CNV get dropped
#  by hand downstream.
# =============================================================================

#' A reader's input file, checked before anything opens it
#'
#' `file_exists()` on `NULL` or a vector is a base R error, and a directory
#' passes the existence test while every reader wants a file.  One guard names
#' the argument and the path for all three readers.
#' @param path the caller's `path`.
#' @param what the record kind, for the message: `"VCF"` or `"FASTA"`.
#' @return `path`, unchanged, when it is one existing file.
#' @keywords internal
#' @noRd
reader_path <- function(path, what) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    abort_msg("`path` must be one {what} file path.")
  }
  if (!fs::file_exists(path)) abort_msg("{what} not found: {path}")
  if (fs::is_dir(path)) {
    abort_msg("`path` names a directory, not a {what} file: {path}")
  }
  path
}

#' Normalise an allele string to plain `ACGTN`
#'
#' Upper-cases the string and folds any IUPAC ambiguity code (R, Y, S, W, K,
#' M, B, D, H, V) - or anything else non-ACGTN - to `N`, so downstream design
#' can treat "unknown base" uniformly instead of dropping the whole record.
#' @keywords internal
#' @param x character vector of allele strings.
#' @return a character vector the same length as `x`, upper-cased and with any
#'   non-`ACGTN` character (including IUPAC ambiguity codes) replaced by `N`.
#' @examples
#' normalise_allele(c("a", "R", "g"))
#' @family input
#' @export
normalise_allele <- function(x) {
  x <- toupper(x)
  x <- chartr("RYSWKMBDHV", "NNNNNNNNNN", x) # IUPAC ambiguity -> N
  gsub("[^ACGTN]", "N", x)
}

# The one admission rule both VCF readers apply: a symbolic REF/ALT (`<DEL>`)
# or an absent allele (`*`, `.`) is unrepresentable, because
# `normalise_allele()` would turn it into a fake `N`.
vcf_symbolic_allele <- function(x) stringr::str_detect(x, "^<|[*.]")

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
#' @return a tibble with columns `chrom`, `pos`, `ref`, `alt`, one row per ALT
#'   allele; symbolic/absent ALT alleles (e.g. `<DEL>`, `*`, `.`) are dropped
#'   with a message naming how many, and the kept alleles are normalised with
#'   [normalise_allele()].  A short oligo cannot be designed across a
#'   structural allele, so those records are unrepresentable rather than lost
#'   in silence.
#' @examples
#' vcf <- withr::local_tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
#'   "1\t250\t.\tA\tG\t.\t.\t."
#' ), vcf)
#' read_variants_vcf(vcf)
#' @family input
#' @export
read_variants_vcf <- function(path, chrom = NULL, start = NULL, end = NULL,
                              backend = c("vroom", "VariantAnnotation")) {
  backend <- match.arg(backend)
  path <- reader_path(path, "VCF")

  v <- switch(backend,
    vroom = read_variants_vcf_vroom(path),
    VariantAnnotation = read_variants_vcf_va(path, chrom, start, end)
  )

  bad_pos <- is.na(v$pos) | v$pos < 1L
  if (any(bad_pos)) {
    abort_msg(gettextf(
      "%s VCF record(s) have a POS outside the positive integer range (first on %s).",
      sum(bad_pos), v$chrom[[which(bad_pos)[[1L]]]]
    ))
  }

  if (!is.null(chrom)) v <- dplyr::filter(v, .data$chrom == ath_chrom(chrom))
  if (!is.null(start)) v <- dplyr::filter(v, .data$pos >= start)
  if (!is.null(end)) v <- dplyr::filter(v, .data$pos <= end)

  # drop symbolic/absent alleles (`<DEL>`, `*`, `.`) before normalisation, so
  # they are not mangled into a fake `N` allele - and say how many went
  drop_sym <- vcf_symbolic_allele(v$ref) | vcf_symbolic_allele(v$alt)
  n_dropped <- sum(drop_sym)
  if (n_dropped > 0L) {
    cli::cli_inform(gettext(paste0(
      "{n_dropped} VCF record{?s} with a symbolic or absent ALT ",
      "(e.g. `<DEL>`, `*`) dropped: short oligos cannot be designed ",
      "across structural alleles."
    )))
  }

  v |>
    dplyr::filter(!drop_sym) |>
    dplyr::mutate(
      ref = normalise_allele(.data$ref),
      alt = normalise_allele(.data$alt)
    ) |>
    dplyr::filter(nzchar(.data$ref), nzchar(.data$alt))
}

#' Read per-line genotypes from a VCF (GT/AD)
#'
#' A variant table says *what* varies; a genotype table says *who* carries it.
#' `read_variant_genotypes()` streams the sample columns of a VCF into long
#' `line`/`chrom`/`pos`/`ref`/`alt` rows, with `genotype` as ALT dosage
#' (`0`, `1`, `2`, `NA`) and `depth` from `AD` when present.  The result feeds
#' [marker_informativeness()] directly.
#'
#' The record is biallelic: a multi-allelic site is dropped and counted, because
#' a dosage against "the ALT" is ambiguous once a record carries two of them.
#' A symbolic or absent ALT (`<DEL>`, `*`, `.`) is dropped and counted for the
#' same reason [read_variants_vcf()] drops it: no short oligo exists for it, and
#' `normalise_allele()` would otherwise turn it into a fake `N` allele.
#' Dosage comes from `GT`; when `GT` is missing, `AD` decides it (rounded
#' `2 * alt / (ref + alt)`).  A genotype naming an allele other than REF/ALT
#' is `NA`, never silently folded to the first ALT.
#' @param path VCF path, plain or bgzipped.
#' @param samples optional sample names to keep; `NULL` keeps every sample.
#' @return a tibble with columns `line`, `chrom`, `pos`, `ref`, `alt`,
#'   `genotype`, `depth`, class `genotypes`; one row per sample and record.
#' @examples
#' vcf <- withr::local_tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   paste0(
#'     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
#'     "L1\tL2"
#'   ),
#'   "1\t250\t.\tA\tG\t.\t.\t.\tGT:AD\t0/0:10,0\t0/1:6,4"
#' ), vcf)
#' read_variant_genotypes(vcf)
#' @family input
#' @export
read_variant_genotypes <- function(path, samples = NULL) {
  path <- reader_path(path, "VCF")
  all_samples <- vcf_sample_names(path)
  if (!length(all_samples)) {
    abort_msg("{path} has no sample columns; there is no genotype to read.")
  }
  if (is.null(samples)) {
    samples <- all_samples
  } else {
    samples <- as.character(samples)
    missing <- setdiff(samples, all_samples)
    if (length(missing)) {
      abort_msg("{path} has no sample{?s} {.val {missing}}.")
    }
  }
  idx <- match(samples, all_samples)

  con <- file(path, "r")
  on.exit(close(con), add = TRUE)
  rows <- list()
  n_multi <- 0L
  n_sym <- 0L
  n_bad <- 0L
  n_bad_pos <- 0L
  repeat {
    lines <- readLines(con, n = 1000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines, "#") & nzchar(lines)]
    for (line in lines) {
      f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(f) < 9L + length(all_samples)) {
        n_bad <- n_bad + 1L
        next
      }
      if (grepl(",", f[[5L]], fixed = TRUE)) {
        n_multi <- n_multi + 1L
        next
      }
      if (vcf_symbolic_allele(f[[4L]]) || vcf_symbolic_allele(f[[5L]])) {
        n_sym <- n_sym + 1L
        next
      }
      pos <- suppressWarnings(as.integer(f[[2L]]))
      if (is.na(pos) || pos < 1L) {
        n_bad_pos <- n_bad_pos + 1L
        next
      }
      fmt <- strsplit(f[[9L]], ":", fixed = TRUE)[[1L]]
      gi <- match("GT", fmt)
      ai <- match("AD", fmt)
      if (is.na(gi) && is.na(ai)) {
        n_bad <- n_bad + 1L
        next
      }
      calls <- vapply(f[9L + idx], function(field) {
        v <- strsplit(field, ":", fixed = TRUE)[[1L]]
        gt <- if (is.na(gi) || gi > length(v)) NA_character_ else v[[gi]]
        ad <- if (is.na(ai) || ai > length(v)) NA_character_ else v[[ai]]
        genotype_call(gt, ad)
      }, c(genotype = NA_real_, depth = NA_real_))
      rows[[length(rows) + 1L]] <- tibble::tibble(
        line = samples,
        chrom = ath_chrom(f[[1L]]),
        pos = pos,
        ref = normalise_allele(f[[4L]]),
        alt = normalise_allele(f[[5L]]),
        genotype = unname(calls["genotype", ]),
        depth = unname(calls["depth", ])
      )
    }
  }
  if (n_multi > 0L) {
    cli::cli_inform(gettext(paste0(
      "{n_multi} multi-allelic VCF record{?s} dropped: ",
      "ALT dosage is ambiguous with more than one ALT."
    )))
  }
  if (n_sym > 0L) {
    cli::cli_inform(gettext(paste0(
      "{n_sym} VCF record{?s} with a symbolic or absent ALT ",
      "(e.g. `<DEL>`, `*`) dropped: short oligos cannot be designed ",
      "across structural alleles."
    )))
  }
  if (n_bad_pos > 0L) {
    cli::cli_inform(gettext(
      "{n_bad_pos} VCF record{?s} dropped: POS is not a positive integer."
    ))
  }
  if (n_bad > 0L) {
    cli::cli_inform(gettext(
      "{n_bad} VCF record{?s} dropped: no FORMAT/GT or FORMAT/AD."
    ))
  }
  out <- if (length(rows)) {
    tibble::as_tibble(do.call(rbind, rows))
  } else {
    tibble::tibble(
      line = character(), chrom = character(), pos = integer(),
      ref = character(), alt = character(), genotype = numeric(),
      depth = numeric()
    )
  }
  class(out) <- c("genotypes", class(out))
  out
}

#' ALT dosage and depth from one FORMAT field (internal)
#' @keywords internal
#' @noRd
genotype_call <- function(gt, ad) {
  ref_depth <- alt_depth <- NA_real_
  if (!is.na(ad) && nzchar(ad) && ad != ".") {
    d <- suppressWarnings(as.numeric(strsplit(ad, ",", fixed = TRUE)[[1L]]))
    if (length(d) >= 2L && all(is.finite(d[1:2]))) {
      ref_depth <- d[[1L]]
      alt_depth <- d[[2L]]
    }
  }
  depth <- if (is.na(ref_depth)) NA_real_ else ref_depth + alt_depth
  dosage <- if (is.na(gt) || gt %in% c(".", "./.", ".|.")) {
    NA_real_
  } else {
    alleles <- strsplit(gt, "[/|]", perl = TRUE)[[1L]]
    a <- suppressWarnings(as.integer(alleles))
    a <- a[!is.na(a)]
    if (length(a) != 2L || any(a > 1L)) NA_real_ else sum(a)
  }
  if (is.na(dosage) && !is.na(depth) && depth > 0) {
    dosage <- round(2 * alt_depth / depth)
  }
  c(genotype = dosage, depth = depth)
}

#' The sample names of a VCF, scanned without a line cap (internal)
#' @keywords internal
#' @noRd
vcf_sample_names <- function(path) {
  con <- file(path, "r")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) {
      abort_msg("{path} has no '#CHROM' header line.")
    }
    if (startsWith(line, "##") || !nzchar(line)) next
    if (!startsWith(line, "#CHROM")) {
      abort_msg("{path}: the first non-'##' line does not start with '#CHROM'.")
    }
    nms <- strsplit(sub("^#", "", line), "\t", fixed = TRUE)[[1L]]
    return(if (length(nms) > 9L) nms[-(1:9)] else character())
  }
}

#' One variant row, as a typed input
#'
#' A design function needs one locus, and a subset of a variants table is a
#' tibble that still claims to be many.  `as_variant_row()` selects exactly one
#' row and names the contract: columns `chrom`, `pos`, `ref`, `alt`, class
#' `variant`.  Zero or several matching rows is an error, not a silent
#' first-row choice.
#' @param x a variants tibble (`chrom`, `pos`, `ref`, `alt`).
#' @param id,chrom,pos optional selectors; give at most what names one row.
#' @return a one-row tibble of class `variant`.
#' @examples
#' as_variant_row(marker_example("SNP"), chrom = "Chr1")
#' @family input
#' @export
as_variant_row <- function(x, id = NULL, chrom = NULL, pos = NULL) {
  if (!is.data.frame(x)) abort_msg("`x` must be a variants table.")
  need_cols(x, c("chrom", "pos", "ref", "alt"))
  if (!is.null(id)) {
    need_cols(x, "id", i = gettext("`id` selects on the `id` column."))
    x <- x[!is.na(x$id) & as.character(x$id) %in% as.character(id), , drop = FALSE]
  }
  if (!is.null(chrom)) {
    x <- x[as.character(x$chrom) == ath_chrom(chrom), , drop = FALSE]
  }
  if (!is.null(pos)) {
    x <- x[as.integer(x$pos) %in% as.integer(pos), , drop = FALSE]
  }
  if (nrow(x) != 1L) {
    abort_msg(gettextf(
      "`as_variant_row()` needs exactly one variant; %d matched. Give `id`, `chrom` or `pos`.",
      nrow(x)
    ))
  }
  out <- tibble::as_tibble(x[, c("chrom", "pos", "ref", "alt"), drop = FALSE]) |>
    dplyr::mutate(
      chrom = ath_chrom(.data$chrom),
      pos = as.integer(.data$pos),
      ref = normalise_allele(.data$ref),
      alt = normalise_allele(.data$alt)
    )
  if (anyNA(out$pos) || any(!nzchar(out$ref)) || any(!nzchar(out$alt))) {
    abort_msg("The variant row needs a non-missing `pos` and non-empty alleles.")
  }
  class(out) <- c("variant", class(out))
  out
}

#' Streaming VCF reader (`vroom`): CHROM/POS/REF/ALT, multi-allelic split
#' @keywords internal
#' @noRd
read_variants_vcf_vroom <- function(path) {
  rlang::check_installed("vroom", reason = gettext("to read large VCF files"))
  nms <- vcf_header_names(path)

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

#' The `#CHROM` header names of a VCF, scanned without a line cap
#'
#' A VCF may carry any number of `##` metadata lines (some build thousands);
#' the header is the last such line and starts with a single `#CHROM`.  Reading
#' to find it is bounded by the metadata block, not by the file size.
#' @keywords internal
#' @noRd
vcf_header_names <- function(path) {
  con <- file(path, "r")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) {
      abort_msg("{path} has no '#CHROM' header line.")
    }
    if (startsWith(line, "##") || !nzchar(line)) next
    if (!startsWith(line, "#CHROM")) {
      abort_msg(
        "{path}: the first non-'##' line does not start with '#CHROM'."
      )
    }
    nms <- strsplit(sub("^#", "", line), "\t", fixed = TRUE)[[1L]]
    missing <- setdiff(c("CHROM", "POS", "REF", "ALT"), nms)
    if (length(missing)) {
      abort_msg(
        "{path}: the '#CHROM' header is missing column {.val {missing}}.",
        i = "A VCF must carry CHROM, POS, REF and ALT for read_variants_vcf()."
      )
    }
    return(nms)
  }
}

#' Canonical VCF reader (`VariantAnnotation::readVcf`)
#'
#' With a region and an indexed (bgzip+tabix) file, only that region is read
#' (`ScanVcfParam(which=)`); the header's contig names are matched to `chrom` with
#' the same `ath_chrom()` normalisation used everywhere else.
#' @keywords internal
#' @noRd
read_variants_vcf_va <- function(path, chrom = NULL, start = NULL, end = NULL) {
  load_namespace_quietly("VariantAnnotation", reason = gettext("for the readVcf backend"))
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
#' @return the path to the plain, `.fai`-indexed FASTA (a character scalar).
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">Chr1", "ACGTACGTACGT"), fa)
#' fasta_index(fa)
#' @family input
#' @export
fasta_index <- function(path) {
  rlang::check_installed("Rsamtools", reason = gettext("for indexed FASTA random access"))
  if (!fs::file_exists(path)) abort_msg("FASTA not found: {path}")
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
#' @return a `DNAString` with the selected positions replaced by `N`.
#' @examples
#' mask_positions("ACGTACGT", c(2, 5))
#' @family input
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

#' Bundled example target variants (one per marker type)
#'
#' A small, typed catalogue of real 1001-Genomes variants, one target per marker
#' type and spread across chromosomes and genes: biallelic `"SNP"` (transition
#' and transversion), a multi-allelic site with two ALTs (non-`NA` `alt2`) that
#' one degenerate primer covers ([marker_site_class()]), `"MNP"` (a contiguous
#' substitution block), and VCF-anchored `"insertion"`/`"deletion"`.  Each locus
#' was chosen for the lowest local variant load in a real gene rather than dumped
#' from one hyper-variable region, and the SNP rows are verified to design (see
#' `dev/render-targets.R`).  Handy for trying the pipeline without pulling data
#' from Ensembl; `ref`/`alt` are the genomic alleles, so the table feeds
#' [design_kasp_markers()] directly.  The rows are *Arabidopsis thaliana*
#' (TAIR10/Atlas) demo targets for use with [ath_genome()], not a species
#' default: pair them with any genome whose sequence names match `chrom`.
#'
#' @param type optional `type` to keep (`"SNP"`, `"MNP"`, `"insertion"`,
#'   `"deletion"`); `NULL` returns every example.
#' @return tibble with `id`, `chrom`, `pos`, `ref`, `alt`, `alt2` (second ALT of a
#'   multi-allelic site, `NA` otherwise), `type`, `consequence`, `transition`,
#'   `gene` (the AGI locus the target sits in).
#' @examples
#' marker_example()
#' marker_example("deletion")
#' @family input
#' @export
marker_example <- function(type = NULL) {
  f <- system.file("extdata", "example_targets.csv", package = "tidymarkers")
  if (!nzchar(f)) abort_msg("Bundled example data not found.")
  ex <- readr::read_csv(f, show_col_types = FALSE)
  if (!is.null(type)) ex <- ex[ex$type %in% type, , drop = FALSE]
  ex
}

read_fasta_strict <- function(path, read) {
  # Biostrings warns and silently drops a letter it does not know; every base
  # after it moves, so a shortened sequence is worse than no sequence
  tryCatch(
    withCallingHandlers(
      read(),
      warning = function(w) {
        if (grepl("invalid one-letter sequence codes", conditionMessage(w))) {
          abort_msg(
            "FASTA {path} has bases outside the DNA/IUPAC alphabet; refusing to shift coordinates silently."
          )
        }
      }
    ),
    error = function(e) {
      if (grepl("invalid DNA letters", conditionMessage(e))) {
        abort_msg(
          "FASTA {path} has bases outside the DNA/IUPAC alphabet; refusing to shift coordinates silently."
        )
      }
      stop(e)
    }
  )
}

# Duplicate names are not cosmetic: `g[["Chr1"]]` would pick one record while
# the indexer picks another, so the sequence behind a coordinate depends on the
# backend.  The header scan catches the raw file; the object check catches a
# `chrom_map` that sends two records to one name.
fasta_duplicate_name <- function(path) {
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  seen <- character(0)
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) {
      break
    }
    heads <- sub("^>", "", lines[startsWith(lines, ">")])
    if (length(heads)) {
      nm <- sub("\\s.*$", "", heads)
      dup <- anyDuplicated(c(seen, nm))
      if (dup) {
        return(nm[[dup - length(seen)]])
      }
      seen <- c(seen, nm)
    }
  }
  NULL
}

check_unique_fasta_names <- function(path) {
  dup <- fasta_duplicate_name(path)
  if (!is.null(dup)) {
    abort_msg(
      gettextf(
        "FASTA %s names sequence %s more than once; chromosome names must be unique.",
        path, dup
      ),
      class = "tidymarkers_duplicate_sequence_error"
    )
  }
  invisible(path)
}

check_unique_sequences <- function(g, path) {
  dup <- anyDuplicated(names(g))
  if (dup) {
    abort_msg(
      gettextf(
        "FASTA %s names sequence %s more than once after mapping; chromosome names must be unique.",
        path, names(g)[[dup]]
      ),
      class = "tidymarkers_duplicate_sequence_error"
    )
  }
  invisible(g)
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
#' @param build optional release id; stored as `genome_build` metadata so a
#'   [genome_index()] made from the result carries it.
#' @return a `DNAStringSet`; sequence names are stripped of their deflines and,
#'   when `chrom_map` is given, replaced by the mapped names (only mapped
#'   sequences are kept).
#' @examples
#' fa <- withr::local_tempfile(fileext = ".fa")
#' writeLines(c(">1", "ACGTACGTACGT"), fa)
#' read_genome_fasta(fa, chrom_map = c("1" = "Chr1"))
#' @family input
#' @export
read_genome_fasta <- function(path, chrom_map = NULL, seqids = NULL,
                              use_index = TRUE, build = NULL) {
  path <- reader_path(path, "FASTA")
  build <- check_genome_build(build)

  normalise <- function(g) {
    nm <- sub("\\s.*$", "", names(g))
    dup <- anyDuplicated(nm)
    if (dup) {
      abort_msg(
        gettextf(
          "FASTA %s names sequence %s more than once; chromosome names must be unique.",
          path, nm[[dup]]
        ),
        class = "tidymarkers_duplicate_sequence_error"
      )
    }
    if (!is.null(chrom_map)) {
      keep <- nm %in% names(chrom_map)
      g <- g[keep]
      names(g) <- unname(chrom_map[nm[keep]])
    } else {
      names(g) <- nm
    }
    # a name lookup would pick the first duplicate; filtering is positional
    if (!is.null(seqids)) g <- g[names(g) %in% seqids]
    if (!is.null(build)) {
      md <- S4Vectors::metadata(g) %||% list()
      md$genome_build <- build
      S4Vectors::metadata(g) <- md
    }
    check_unique_sequences(g, path)
    g
  }

  if (isTRUE(use_index) && !is.null(seqids)) {
    check_unique_fasta_names(path)
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
      return(normalise(Biostrings::DNAStringSet()))
    }
    gr <- GenomicRanges::GRanges(sel$name, IRanges::IRanges(1L, sel$length))
    g <- read_fasta_strict(path, function() {
      Rsamtools::scanFa(Rsamtools::FaFile(plain), gr)
    })
    return(normalise(g))
  }
  normalise(read_fasta_strict(path, function() {
    Biostrings::readDNAStringSet(path)
  }))
}

# =============================================================================
#  Interchange: variant tibble <-> GRanges / VCF (roadmap Z10.3)
# =============================================================================
#  One row of a variants tibble is the VCF `1-based POS + REF + one ALT` form,
#  so the GRanges anchor is that base POS (not POS + 1, the UCSC/VEP convention)
#  and an InDel REF keeps its anchor base.  `depth` (Y5.3) travels as an mcol.
#  Writing a VCF keeps one ALT per record, so a multi-allelic site becomes one
#  record per allele - which is exactly how [read_variants_vcf()] reads it back.
# =============================================================================

#' Variants tibble <-> `GRanges`
#'
#' Adapter between the package's variant tibble (`chrom`, `pos`, `ref`, `alt`, one
#' row per ALT allele, alleles as strings so SNP/MNP/InDel all survive) and a
#' Bioconductor `GRanges`.  Every other column becomes an mcol, so `depth` (and
#' `alt_depth`, `genotype`, `pop_freq`, ...) round-trips.  The range anchor is
#' the 1-based VCF `POS` (an InDel REF keeps its anchor base), so an InDel range
#' is as wide as its REF - filter with `nearest()` before comparing.
#' @param x a variants tibble (`chrom`, `pos`, `ref`, `alt`).
#' @param extra optional extra mcols to carry (e.g. `read_variant_depth()` rows);
#'   a column already in `x` is overwritten by `extra`.
#' @return [marker_as_granges()]: a `GRanges` with mcols `ref`, `alt` and every
#'   other column of `x`; `NULL` (with a `cli` message) when `GenomicRanges` is
#'   not installed.  [marker_from_granges()]: the variant tibble again, with
#'   `chrom`/`pos`/`ref`/`alt` first and the remaining mcols after.
#' @examples
#' v <- tibble::tibble(chrom = "Chr1", pos = c(10L, 20L), ref = c("A", "AT"), alt = c("G", "A"))
#' gr <- marker_as_granges(v, extra = tibble::tibble(depth = c(30, 12)))
#' if (!is.null(gr)) marker_from_granges(gr)
#' @family input
#' @export
marker_as_granges <- function(x, extra = NULL) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    cli::cli_inform(gettext("GenomicRanges is not installed; returning `NULL`."))
    return(NULL)
  }
  need <- c("chrom", "pos", "ref", "alt")
  need_cols(x, need)
  df <- tibble::as_tibble(x)
  if (!is.null(extra)) {
    extra <- tibble::as_tibble(extra)
    df[names(extra)] <- extra
  }
  ref <- as.character(df$ref)
  if (any(is.na(ref) | !nzchar(ref))) {
    abort_msg("`ref` must be a non-empty allele string on every row.")
  }
  gr <- GenomicRanges::GRanges(
    seqnames = as.character(df$chrom),
    ranges = IRanges::IRanges(as.integer(df$pos), width = nchar(ref))
  )
  S4Vectors::mcols(gr) <- df[c("ref", "alt", setdiff(names(df), need))]
  gr
}

#' @rdname marker_as_granges
#' @export
marker_from_granges <- function(x) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    cli::cli_inform(gettext("GenomicRanges is not installed; returning `NULL`."))
    return(NULL)
  }
  if (!methods::is(x, "GRanges")) abort_msg("`x` must be a GRanges.")
  m <- S4Vectors::mcols(x)
  need_cols(m, c("ref", "alt"), i = gettext("`mcols()` is the metadata-column table."))
  tibble::as_tibble(m) |>
    dplyr::mutate(
      chrom = ath_chrom(as.character(GenomicRanges::seqnames(x))),
      pos = as.integer(GenomicRanges::start(x)),
      .before = 1L
    ) |>
    dplyr::relocate("ref", "alt", .after = "pos") |>
    dplyr::mutate(dplyr::across(c("ref", "alt"), as.character))
}

#' Write a variants tibble as a VCF
#'
#' Writes exactly the four columns [read_variants_vcf()] consumes, one record per
#' ALT allele (a multi-allelic site is written once per allele, as it reads back)
#' and `depth` as `INFO/DP` when the column is present.  The round trip is
#' lossless for `chrom` (via `ath_chrom()` on read), `pos`, `ref`, `alt` and
#' `depth`; other columns (e.g. `pop_freq`) are not written.
#' @param x a variants tibble (`chrom`, `pos`, `ref`, `alt`, optional `depth`).
#' @param path output VCF path; `NULL` writes a temporary file you own.
#' @return the path written (invisibly when `path` is given).
#' @examples
#' v <- tibble::tibble(chrom = "4", pos = 10L, ref = "A", alt = "G", depth = 30L)
#' f <- write_variants_vcf(v)
#' read_variants_vcf(f)
#' @family input
#' @export
write_variants_vcf <- function(x, path = NULL) {
  need_cols(x, c("chrom", "pos", "ref", "alt"))
  x <- tibble::as_tibble(x)
  if (!"depth" %in% names(x)) x$depth <- NA_real_
  if (any(is.na(x$chrom) | is.na(x$pos) | is.na(x$ref) | is.na(x$alt))) {
    abort_msg("`x` must have no `NA` in `chrom`/`pos`/`ref`/`alt`.")
  }
  records <- paste(
    as.character(x$chrom), as.integer(x$pos), ".", as.character(x$ref),
    as.character(x$alt), ".", ".",
    vapply(x$depth, function(d) {
      if (is.na(d)) "." else paste0("DP=", as.integer(d))
    }, character(1)),
    sep = "\t"
  )
  body <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    records
  )
  if (is.null(path)) {
    path <- tempfile(fileext = ".vcf")
    readr::write_lines(body, path)
    return(path)
  }
  readr::write_lines(body, path)
  invisible(path)
}
