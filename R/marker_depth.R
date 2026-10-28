# =============================================================================
#  marker_depth.R  --  sequencing depth as a confidence weight (roadmap Y5.3)
# =============================================================================
#  A genotype is only as good as the reads behind it.  Low depth hurts EVERY
#  call, not just heterozygotes (a hexaploid/paralog mapping can flip even a
#  homozygote), so the depth term 1 - exp(-depth / min_depth) weights every
#  genotype, and a heterozygote balance term down-weights only an imbalanced het
#  (its alt fraction is far from 1/2).  `read_variant_depth()` pulls DP/AD/GT
#  from a VCF, `variant_confidence()` turns them into `weight`/`confident`, and
#  `filter_targets_by_depth()` gates the solver on the confident rows.
# =============================================================================

#' Read per-sample depth (DP/AD/GT) from a VCF
#'
#' Streams the DP/AD/GT fields of one sample out of a VCF into a tidy table, one
#' row per ALT allele (a multi-allelic `ALT` is split on `,`; its `AD` is
#' `ref,alt1,alt2,...`, so the ALT's depth is `AD[alt_index + 1]`).  `chrom` is
#' normalised with `ath_chrom()`.  When the VCF has no sample columns, the reader
#' falls back to `INFO/DP` and leaves `alt_depth`/`genotype` as `NA`.  Missing
#' values stay `NA`; an absent DP/AD/GT never errors.
#' @param path VCF path; plain (`vcf`) or bgzipped (`vcf.gz`) are both read by
#'   `vroom`.
#' @param sample which sample column to read: an index (default `1L`, the first)
#'   or a column name.
#' @return a tibble with columns `chrom`, `pos`, `ref`, `alt`, `depth` (DP),
#'   `alt_depth` (AD of that ALT) and `genotype` (GT string), one row per ALT
#'   allele.
#' @examples
#' vcf <- withr::local_tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
#'   "1\t250\t.\tA\tG\t.\t.\t.\tGT:AD:DP\t0/1:5,5:10"
#' ), vcf)
#' read_variant_depth(vcf)
#' @family KASP input
#' @export
read_variant_depth <- function(path, sample = 1L) {
  rlang::check_installed("vroom", reason = gettext("to read VCF files"))
  if (!fs::file_exists(path)) cli::cli_abort(gettext("VCF not found: {path}"))
  hdr <- readr::read_lines(path, n_max = 1000L)
  hdr <- hdr[startsWith(hdr, "#CHROM")][1L]
  if (is.na(hdr)) cli::cli_abort(gettext("{path} has no '#CHROM' header line."))
  nms <- strsplit(sub("^#", "", hdr), "\t", fixed = TRUE)[[1L]]

  fmt_i <- match("FORMAT", nms)
  sample_names <- if (!is.na(fmt_i)) nms[seq.int(fmt_i + 1L, length(nms))] else character(0)
  if (!length(sample_names)) {
    return(read_variant_depth_sites(path, nms))
  }

  sample_col <- resolve_depth_sample(sample, sample_names)
  cols <- list(
    .default = vroom::col_skip(),
    CHROM = vroom::col_character(),
    POS = vroom::col_integer(),
    REF = vroom::col_character(),
    ALT = vroom::col_character(),
    FORMAT = vroom::col_character()
  )
  cols[[sample_col]] <- vroom::col_character()
  d <- vroom::vroom(
    path,
    comment = "#", col_names = nms,
    col_types = do.call(vroom::cols, cols),
    progress = FALSE
  )
  depth_from_samples(d, sample_col)
}

#' Resolve the `sample` selector to a sample-column name
#' @keywords internal
#' @noRd
resolve_depth_sample <- function(sample, sample_names) {
  if (is.numeric(sample)) {
    i <- as.integer(sample)
    if (length(i) != 1L || is.na(i) || i < 1L || i > length(sample_names)) {
      cli::cli_abort(gettext(
        "`sample` must be a sample index (1..N) or a column name in the VCF."
      ))
    }
    return(sample_names[[i]])
  }
  s <- as.character(sample)
  if (length(s) != 1L || is.na(s) || !s %in% sample_names) {
    cli::cli_abort(glue::glue(gettext(
      "Sample '{s}' is not in the VCF (samples: {paste(sample_names, collapse = ', ')})."
    )))
  }
  s
}

#' DP/AD/GT of one sample, one row per ALT allele
#' @keywords internal
#' @noRd
depth_from_samples <- function(d, sample_col) {
  alts <- stringr::str_split(d$ALT, ",")
  n_alt <- lengths(alts)
  idx <- rep(seq_len(nrow(d)), n_alt)
  alt_index <- sequence(n_alt)
  fmt <- stringr::str_split(d$FORMAT[idx], ":")
  smp <- stringr::str_split(d[[sample_col]][idx], ":")

  pick <- function(key) {
    vapply(seq_along(smp), function(i) {
      j <- match(key, fmt[[i]])
      if (is.na(j) || j > length(smp[[i]])) NA_character_ else smp[[i]][[j]]
    }, character(1))
  }

  depth <- clean_vcf_value(pick("DP"))
  ad <- clean_vcf_value(pick("AD"))
  genotype <- clean_vcf_value(pick("GT"))

  ad_split <- strsplit(ad, ",", fixed = TRUE)
  alt_depth <- vapply(seq_along(ad_split), function(i) {
    a <- ad_split[[i]]
    k <- alt_index[[i]]
    if (is.na(a[[1L]]) || k + 1L > length(a)) {
      return(NA_character_)
    }
    a[[k + 1L]]
  }, character(1))

  tibble::tibble(
    chrom = ath_chrom(d$CHROM[idx]),
    pos = as.integer(d$POS[idx]),
    ref = d$REF[idx],
    alt = unlist(alts, use.names = FALSE),
    depth = as.numeric(depth),
    alt_depth = as.numeric(alt_depth),
    genotype = genotype
  )
}

#' Sites-only VCF: `INFO/DP`, no per-sample fields
#' @keywords internal
#' @noRd
read_variant_depth_sites <- function(path, nms) {
  cols <- vroom::cols(
    .default = vroom::col_skip(),
    CHROM = "c", POS = "i", REF = "c", ALT = "c", INFO = "c"
  )
  d <- vroom::vroom(
    path,
    comment = "#", col_names = nms, col_types = cols,
    progress = FALSE
  )
  alts <- stringr::str_split(d$ALT, ",")
  n_alt <- lengths(alts)
  idx <- rep(seq_len(nrow(d)), n_alt)
  tibble::tibble(
    chrom = ath_chrom(d$CHROM[idx]),
    pos = as.integer(d$POS[idx]),
    ref = d$REF[idx],
    alt = unlist(alts, use.names = FALSE),
    depth = as.numeric(info_dp(d$INFO[idx])),
    alt_depth = NA_real_,
    genotype = NA_character_
  )
}

#' `DP=` out of an INFO field (`NA` when absent)
#' @keywords internal
#' @noRd
info_dp <- function(info) {
  vapply(strsplit(info, ";", fixed = TRUE), function(kv) {
    dp <- kv[startsWith(kv, "DP=")]
    if (!length(dp) || is.na(dp[[1L]])) {
      return(NA_character_)
    }
    sub("^DP=", "", dp[[1L]])
  }, character(1))
}

#' Turn a VCF `.`/empty field into `NA`, leave everything else
#' @keywords internal
#' @noRd
clean_vcf_value <- function(x) {
  x[!is.na(x) & x %in% c(".", "")] <- NA_character_
  x
}

#' Turn DP/AD/GT into a confidence weight
#'
#' Adds, to a table of `chrom`/`pos`/`depth` (and optional `alt_depth`/`genotype`):
#'
#' - `alt_fraction` = `alt_depth / depth` (`NA` when either is `NA` or `depth` is 0);
#' - `genotype_class` from GT when present (`"0/1"`/`"0|1"` -> `"het"`, `"0/0"` ->
#'   `"hom_ref"`, `"1/1"` -> `"hom_alt"`, otherwise `"unknown"`); when GT is absent,
#'   derived from `alt_fraction` (`< 0.2` hom_ref, `> 0.8` hom_alt, else het);
#' - `weight` in `[0, 1]`: the depth term `1 - exp(-depth / min_depth)` applies to
#'   **every** genotype (so the weight rises with depth; the penalty falls).  A
#'   `"het"` with an `alt_fraction` is multiplied by the balance term
#'   `pmax(0, 1 - abs(alt_fraction - 0.5) / 0.5)`; homozygotes are never penalised
#'   (their fraction is legitimately 0 or 1).  `NA` depth -> `NA` weight.
#' - `confident` = `!is.na(weight) & depth > 0 & weight >= min_weight`.
#'
#' The depth term is the point of roadmap Y5.3: low coverage degrades every
#' genotype, not just heterozygotes, so a low-depth homozygote is not confident.
#' @param depth a table with `chrom`, `pos`, `depth` (DP) and optional `alt_depth`
#'   (AD of the ALT) and `genotype` (GT), e.g. from [read_variant_depth()].
#' @param min_depth the depth scale of the exponential term (default `3L`): at
#'   this depth the weight is `1 - exp(-1)` ~ 0.632.
#' @param min_weight the `weight` a genotype must reach to be `confident`
#'   (default `0.5`).
#' @return the input with columns `alt_fraction`, `genotype_class`, `weight` and
#'   `confident` appended.
#' @examples
#' variant_confidence(tibble::tibble(
#'   chrom = "Chr1", pos = 1L, depth = 10L, alt_depth = 0L, genotype = "0/0"
#' ))
#' @family KASP QC
#' @export
variant_confidence <- function(depth, min_depth = 3L, min_weight = 0.5) {
  if (!is.numeric(min_depth) || length(min_depth) != 1L || is.na(min_depth) || min_depth <= 0) {
    cli::cli_abort(gettext("`min_depth` must be a single positive number."))
  }
  x <- tibble::as_tibble(depth)
  missing <- setdiff(c("chrom", "pos", "depth"), names(x))
  if (length(missing)) {
    cli::cli_abort(glue::glue(
      gettext("`depth` is missing column(s): {paste(missing, collapse = ', ')}.")
    ))
  }
  x$chrom <- ath_chrom(x$chrom)

  d <- as.numeric(x$depth)
  ad <- if ("alt_depth" %in% names(x)) as.numeric(x$alt_depth) else rep(NA_real_, nrow(x))
  gt <- if ("genotype" %in% names(x)) as.character(x$genotype) else rep(NA_character_, nrow(x))

  alt_fraction <- ifelse(is.na(d) | is.na(ad) | d == 0, NA_real_, ad / d)
  genotype_class <- genotype_class(gt, alt_fraction)

  depth_weight <- 1 - exp(-d / min_depth)
  balance <- pmax(0, 1 - abs(alt_fraction - 0.5) / 0.5)
  balance[is.na(balance)] <- 1 # no alt_fraction -> no balance penalty
  weight <- depth_weight * ifelse(genotype_class == "het", balance, 1)
  weight[is.na(d)] <- NA_real_

  out <- x
  out$alt_fraction <- alt_fraction
  out$genotype_class <- genotype_class
  out$weight <- weight
  out$confident <- !is.na(weight) & d > 0 & weight >= min_weight
  out
}

#' Classify a genotype, deriving from `alt_fraction` when GT is absent
#' @keywords internal
#' @noRd
genotype_class <- function(gt, alt_fraction) {
  gt <- as.character(gt)
  frac <- alt_fraction
  out <- rep(NA_character_, length(gt))
  out[!is.na(gt) & (gt == "0/1" | gt == "0|1")] <- "het"
  out[!is.na(gt) & gt == "0/0"] <- "hom_ref"
  out[!is.na(gt) & gt == "1/1"] <- "hom_alt"
  out[!is.na(gt) & is.na(out)] <- "unknown" # GT present but not one of the four

  absent <- is.na(gt)
  out[absent & !is.na(frac) & frac < 0.2] <- "hom_ref"
  out[absent & !is.na(frac) & frac > 0.8] <- "hom_alt"
  out[absent & !is.na(frac) & frac >= 0.2 & frac <= 0.8] <- "het"
  out[absent & is.na(frac)] <- "unknown"
  out
}

#' Drop solver targets that a depth table marks not confident (internal)
#'
#' A target row is dropped when the depth table has a row at the same
#' `(chrom, pos)` with `confident == FALSE`.  A target with **no** depth row is
#' kept: absence of DP is "unknown", not "low confidence".  `NULL` depth leaves
#' the targets unchanged.  If nothing survives, this aborts naming the depth gate
#' and its knobs.
#' @keywords internal
#' @noRd
filter_targets_by_depth <- function(targets, depth) {
  if (is.null(depth)) {
    return(targets)
  }
  tg <- tibble::as_tibble(targets)
  if (!nrow(tg)) {
    return(targets)
  }
  d <- tibble::as_tibble(depth)
  missing <- setdiff(c("chrom", "pos", "confident"), names(d))
  if (length(missing)) {
    cli::cli_abort(glue::glue(
      gettext("The depth table is missing column(s): {paste(missing, collapse = ', ')}.")
    ))
  }
  low <- paste(ath_chrom(d$chrom), d$pos, sep = "|")[!is.na(d$confident) & !d$confident]
  key <- paste(ath_chrom(tg$chrom), tg$pos, sep = "|")
  out <- tg[!key %in% low, , drop = FALSE]
  if (!nrow(out)) {
    msg <- glue::glue(
      "Every target was dropped by the depth gate ",
      "(min_depth = {attr(depth, 'min_depth') %||% 3L}, ",
      "min_weight = {attr(depth, 'min_weight') %||% 0.5})."
    )
    cli::cli_abort(msg)
  }
  out
}
