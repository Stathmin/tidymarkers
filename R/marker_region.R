# =============================================================================
#  marker_region.R  --  a region as a target, and the haplotype -> signal table
# =============================================================================
#  Roadmap Y2.1: a plan may name a *region* instead of a single site.  The
#  region is not a hard-coded list of sites: at solve time it expands to every
#  variant the context places inside it, so the solver still walks sites x
#  chemistries and picks one.  A variant outside the block never becomes a
#  target - that is the whole point of the boundary.
#
#  Linkage is only trustworthy inside a phased block, so the result carries a
#  haplotype -> signal table with an explicit `within_phased_block` flag rather
#  than implying that a marker explains everything in the region.
# =============================================================================

#' Normalise a region to one `chrom`/`start`/`end` row (internal)
#'
#' Accepts a `GRanges` (length one), or a data frame / list with `chrom`,
#' `start` and `end` (1-based, closed, as everywhere in the package).
#' @keywords internal
#' @noRd
normalise_region <- function(region) {
  if (inherits(region, "GRanges")) {
    if (length(region) != 1L) {
      cli::cli_abort(gettext("`region` must be a single GRanges interval."))
    }
    out <- tibble::tibble(
      chrom = as.character(GenomicRanges::seqnames(region)),
      start = as.integer(GenomicRanges::start(region)),
      end = as.integer(GenomicRanges::end(region))
    )
  } else if (is.data.frame(region) || is.list(region)) {
    x <- tibble::as_tibble(region)
    missing <- setdiff(c("chrom", "start", "end"), names(x))
    if (length(missing)) {
      cli::cli_abort(glue::glue(
        gettext("`region` is missing column(s): {paste(missing, collapse = ', ')}.")
      ))
    }
    if (nrow(x) != 1L) {
      cli::cli_abort(gettext("`region` must name exactly one interval."))
    }
    out <- tibble::tibble(
      chrom = as.character(x$chrom),
      start = as.integer(x$start),
      end = as.integer(x$end)
    )
  } else {
    cli::cli_abort(gettext(
      "`region` must be a GRanges or a data frame with `chrom`, `start`, `end`."
    ))
  }
  if (anyNA(c(out$chrom, out$start, out$end)) || out$start < 1L || out$end < out$start) {
    cli::cli_abort(gettext("`region` needs a chromosome and a valid 1-based start <= end."))
  }
  out
}

#' The variant context a region expands against (internal)
#' @keywords internal
#' @noRd
region_variants <- function(plan) {
  v <- plan$variants
  if (is.null(v) && !is.null(plan$haplotypes)) {
    v <- variants_from_haplotypes(normalise_haplotypes(plan$haplotypes))
  }
  if (is.null(v)) {
    cli::cli_abort(gettext(
      "A `region` needs variant context; add marker_context(variants = ) or marker_context(haplotypes = )."
    ))
  }
  v <- tibble::as_tibble(v)
  missing <- setdiff(c("chrom", "pos", "ref", "alt"), names(v))
  if (length(missing)) {
    cli::cli_abort(glue::glue(
      gettext("The variant context is missing column(s): {paste(missing, collapse = ', ')}.")
    ))
  }
  v[, c("chrom", "pos", "ref", "alt")]
}

#' Expand a plan's region to the targets inside it (internal)
#'
#' One target per distinct variant in the (1-based, closed) interval, labelled
#' and typed as a region target.  An empty block is an error, not an empty plan:
#' a region the context cannot fill cannot be designed.
#' @keywords internal
#' @noRd
region_targets <- function(plan) {
  reg <- plan$region
  v <- region_variants(plan)
  v <- v[
    as.character(v$chrom) == reg$chrom &
      as.integer(v$pos) >= reg$start & as.integer(v$pos) <= reg$end, ,
    drop = FALSE
  ]
  v <- dplyr::distinct(v, .data$chrom, .data$pos, .data$ref, .data$alt)
  if (!nrow(v)) {
    cli::cli_abort(glue::glue(
      gettext("No variant lies in the region {reg$chrom}:{reg$start}-{reg$end}.")
    ))
  }
  target_rows(
    chrom = v$chrom, pos = v$pos, ref = v$ref, alt = v$alt,
    data = NULL, type = "region", label = NULL
  )
}

#' The haplotype -> signal table of a solved region (internal)
#'
#' For every solved target, the allele each known haplotype carries and the
#' primer it is expected to bind; `signal` is the readout that follows.  A
#' marker that does not model haplotypes falls back to the target allele alone.
#' Every row says `within_phased_block = TRUE`: the context is a phased table,
#' so the statement holds inside the block and nowhere else (roadmap Y2.1).
#' @keywords internal
#' @noRd
region_signals <- function(plan, report) {
  h <- plan$haplotypes
  if (is.null(h) || !nrow(report)) {
    return(NULL)
  }
  h <- normalise_haplotypes(h)
  freqs <- h |>
    dplyr::distinct(.data$haplotype, .data$freq) |>
    tibble::deframe()
  ok <- report[report$ok, , drop = FALSE]
  if (!nrow(ok)) {
    return(NULL)
  }
  rows <- purrr::map(seq_len(nrow(ok)), function(i) {
    hits <- ok$marker[[i]]$haplotype_hits
    if (is.list(hits) && nrow(hits)) {
      return(tibble::tibble(
        target = as.character(ok$label[[i]]),
        haplotype = as.character(hits$haplotype),
        freq = as.numeric(hits$freq),
        target_allele = as.character(hits$target_allele),
        primer = as.character(hits$primer),
        binds = as.logical(hits$binds),
        signal = ifelse(hits$binds, as.character(hits$primer), "no_bind"),
        within_phased_block = TRUE
      ))
    }
    # a chemistry without a haplotype model: the target allele is the signal
    pos <- as.integer(ok$pos[[i]])
    alt <- as.character(ok$alt[[i]])
    carries <- vapply(
      names(freqs),
      function(hp) any(h$haplotype == hp & h$pos == pos & h$alt == alt),
      logical(1)
    )
    tibble::tibble(
      target = as.character(ok$label[[i]]),
      haplotype = names(freqs),
      freq = as.numeric(freqs),
      target_allele = ifelse(carries, "ALT", "REF"),
      primer = NA_character_,
      binds = NA,
      signal = ifelse(carries, "ALT", "REF"),
      within_phased_block = TRUE
    )
  })
  purrr::list_rbind(rows)
}
