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
#' @param freq haplotype frequencies: `NULL` (equal), a single number, a named
#'   vector, or a two-column table (`haplotype`, `freq`) - the per-haplotype
#'   weights whose carrier sums give each variant's `pop_freq`.
#' @param maf optional per-variant minor-allele-frequency table (`chrom`?, `pos`,
#'   `ref`, `alt`, and `maf`/`freq`) used to override the haplotype-derived
#'   `pop_freq` of the matching variants in `$variants`.
#' @param ... passed to [DECIPHER::AlignSeqs()].
#' @return a list: `reference` (a `DNAString`), `genome` (a one-sequence
#'   `DNAStringSet` named `chrom`), `reference_length`, `haplotypes` (the tidy
#'   table), `variants` (`chrom`, `pos`, `ref`, `alt`, `pop_freq`, `consequence` -
#'   the table [assemble_kasp_assay()] weights variants by), `alignment` (a
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
                                      freq = NULL, maf = NULL, ...) {
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
  f <- normalise_hap_freq(freq, names(a))
  hap_tbl <- tibble::tibble(
    haplotype = variants$haplotype,
    freq = unname(f[variants$haplotype]),
    chrom = variants$chrom, pos = variants$pos, ref = variants$ref, alt = variants$alt
  )
  list(
    reference = Biostrings::DNAString(ref_seq),
    reference_length = nchar(ref_seq),
    chrom = chrom,
    genome = Biostrings::DNAStringSet(stats::setNames(ref_seq, chrom)),
    haplotypes = hap_tbl,
    variants = variants_from_haplotypes(hap_tbl, maf = maf),
    alignment = Biostrings::DNAMultipleAlignment(aln),
    reference_name = ref_name
  )
}

#' Normalise a haplotype-frequency input to one value per sequence
#' @keywords internal
#' @noRd
normalise_hap_freq <- function(freq, nms) {
  n <- length(nms)
  if (is.null(freq)) {
    return(stats::setNames(rep(1 / n, n), nms))
  }
  if (is.data.frame(freq)) {
    if (!all(c("haplotype", "freq") %in% names(freq))) {
      cli::cli_abort(gettext("A `freq` table needs `haplotype` and `freq` columns."))
    }
    v <- stats::setNames(as.numeric(freq$freq), as.character(freq$haplotype))
    out <- v[nms]
    out[is.na(out)] <- 0
    return(out)
  }
  if (length(freq) == 1L) {
    return(stats::setNames(rep(as.numeric(freq), n), nms))
  }
  if (is.null(names(freq))) {
    cli::cli_abort(gettext("A `freq` vector must be named by haplotype."))
  }
  v <- stats::setNames(as.numeric(freq), names(freq))
  out <- v[nms]
  out[is.na(out)] <- 0
  out
}

#' Each variant's pop_freq from a haplotype table (with optional MAF overrides)
#'
#' A variant's weight is the summed frequency of the haplotypes that carry it
#' (capped at 1); an optional `maf` table replaces it for the matching `(chrom,
#' pos, ref, alt)` rows - so real minor-allele frequencies win over the count
#' implied by the alignment.
#' @keywords internal
#' @noRd
variants_from_haplotypes <- function(hap_tbl, maf = NULL) {
  v <- hap_tbl |>
    dplyr::group_by(.data$chrom, .data$pos, .data$ref, .data$alt) |>
    dplyr::summarise(pop_freq = min(1, sum(.data$freq)), .groups = "drop")
  v$consequence <- NA_character_
  apply_maf(v, maf)[, c("chrom", "pos", "ref", "alt", "pop_freq", "consequence")]
}

#' Override a variants table's pop_freq from a MAF table
#' @keywords internal
#' @noRd
apply_maf <- function(v, maf = NULL) {
  if (is.null(maf)) {
    return(v)
  }
  m <- tibble::as_tibble(maf)
  fcol <- intersect(c("maf", "freq"), names(m))[1L]
  if (!all(c("pos", "ref", "alt") %in% names(m)) || is.na(fcol)) {
    cli::cli_abort(gettext("A `maf` table needs `pos`, `ref`, `alt` and `maf`/`freq`."))
  }
  m <- m[, c(intersect("chrom", names(m)), "pos", "ref", "alt", fcol)]
  names(m)[names(m) == fcol] <- "maf"
  by <- intersect(c("chrom", "pos", "ref", "alt"), names(m))
  v <- dplyr::left_join(v, m, by = by, suffix = c("", "_maf"))
  v$pop_freq <- dplyr::coalesce(v$maf, v$pop_freq)
  v$maf <- NULL
  v
}

# nucleotide fill for the alignment view
.aln_base_colours <- c(A = "#4daf4a", C = "#377eb8", G = "#ff7f00", T = "#e41a1c", "-" = "grey88", N = "grey70")

#' Plot an alignment of the haplotypes (a base x sequence matrix)
#'
#' A compact view of [haplotypes_from_alignment()]: one row per haplotype (the
#' reference first), one column per site, filled by the base - so SNPs, MNPs and
#' indels read at a glance.  `type = "variants"` (default) keeps only the columns
#' that differ from the reference; `type = "all"` plots the whole alignment.
#' @param hp a list from [haplotypes_from_alignment()].
#' @param type `"variants"` (only the variant columns) or `"all"`.
#' @param features optional feature table (`start`, `end` in reference
#'   coordinates; optional `name`/`symbol`, `strand`) drawn as a `feature` track
#'   under the alignment (e.g. the GenBank `gene`/`CDS` of the reference).
#' @param ... reserved for future options.
#' @return a `ggplot`.
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "sat2_transcripts.fa", package = "tidykasp")
#' plot_alignment(haplotypes_from_alignment(read_transcripts(f), chrom = "SAT2"))
#' }
#' @family KASP haplotypes
#' @export
plot_alignment <- function(hp, type = c("variants", "all"), features = NULL, ...) {
  rlang::check_installed("ggplot2", reason = gettext("to plot an alignment"))
  if (!is.list(hp) || is.null(hp$alignment)) {
    cli::cli_abort(gettext("`hp` must come from haplotypes_from_alignment()."))
  }
  type <- match.arg(type)
  A <- as.character(hp$alignment)
  m <- do.call(rbind, strsplit(A, "", fixed = TRUE))
  rownames(m) <- names(A)
  ref_row <- m[hp$reference_name, ]
  ref_mat <- matrix(ref_row, nrow = nrow(m), ncol = ncol(m), byrow = TRUE)
  variant_cols <- which(colSums(m != ref_mat) > 0L)
  cols <- if (identical(type, "variants")) variant_cols else seq_len(ncol(m))
  lv <- c(hp$reference_name, setdiff(rownames(m), hp$reference_name))
  # the feature track sits below the sequences
  feat <- normalise_features(features, c(1L, hp$reference_length))
  if (!is.null(feat)) lv <- c(lv, "feature")
  d <- tibble::tibble(
    haplotype = factor(rep(rownames(m), times = length(cols)), levels = lv),
    col = rep(cols, each = nrow(m)),
    base = as.vector(m[, cols, drop = FALSE])
  )
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$col, y = .data$haplotype, fill = .data$base)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = .aln_base_colours, na.value = "white", name = "base") +
    ggplot2::labs(
      x = "alignment column", y = NULL,
      title = glue::glue("Haplotype alignment ({length(variant_cols)} variant columns)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (!is.null(feat)) {
    # the alignment axis is *columns*, so map each feature's reference span onto
    # the reference row's non-gap columns, then clamp to the drawn columns
    cols_of_ref <- which(ref_row != "-")
    nref <- length(cols_of_ref)
    feat$c1 <- cols_of_ref[pmax(pmin(feat$start, nref), 1L)]
    feat$c2 <- cols_of_ref[pmax(pmin(feat$end, nref), 1L)]
    keep <- feat$c2 >= min(cols) & feat$c1 <= max(cols)
    feat <- feat[keep, , drop = FALSE]
    feat$c1 <- pmax(feat$c1, min(cols))
    feat$c2 <- pmin(feat$c2, max(cols))
    if (nrow(feat)) {
      yf <- length(lv)
      p <- p + ggplot2::geom_rect(
        data = feat,
        ggplot2::aes(
          xmin = .data$c1 - 0.5, xmax = .data$c2 + 0.5,
          ymin = yf - 0.42, ymax = yf + 0.42
        ),
        inherit.aes = FALSE, fill = "grey72", colour = "grey40"
      )
      if (any(!is.na(feat$name))) {
        p <- p + ggplot2::geom_text(
          data = feat[!is.na(feat$name), ],
          ggplot2::aes(x = (.data$c1 + .data$c2) / 2, y = yf, label = .data$name),
          inherit.aes = FALSE, size = 2.1, fontface = "bold"
        )
      }
    }
  }
  p
}
