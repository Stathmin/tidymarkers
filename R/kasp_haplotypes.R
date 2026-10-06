# =============================================================================
#  kasp_haplotypes.R  --  haplotype-aware targets (frequencies + co-segregation)
# =============================================================================
#  A tidy variant table says "there is a variant at this position".  A haplotype
#  table says "this *combination* of alleles travels together, and at this
#  frequency" - which is what actually decides whether a primer binds one
#  haplotype and not another.
#
#  Input shape (`haplotypes=` of assemble_kasp_assay()):
#    tibble(haplotype, freq, chrom, pos, ref, alt), one row per variant a
#    haplotype carries; a row at the target position marks an ALT carrier, and a
#    haplotype with no such row carries REF there.  `freq` is the population
#    frequency of that haplotype.
#
#  From it we derive (a) per-position weighted frequencies (the `pop_freq`
#  column the footprint load already understands), (b) co-segregating positions
#  (baked into the ALT primer via `linked=`), and (c) which primer each haplotype
#  is expected to bind, for the output.
# =============================================================================

#' Normalise and validate a haplotype table
#' @keywords internal
#' @noRd
normalise_haplotypes <- function(haplotypes) {
  need <- c("haplotype", "freq", "chrom", "pos", "ref", "alt")
  missing <- setdiff(need, names(haplotypes))
  if (length(missing)) {
    rlang::abort(glue::glue(gettext(
      "`haplotypes` is missing column(s): {paste(missing, collapse = ', ')}."
    )))
  }
  h <- tibble::as_tibble(haplotypes) |>
    dplyr::transmute(
      haplotype = as.character(.data$haplotype),
      freq = as.numeric(.data$freq),
      chrom = as.character(.data$chrom),
      pos = as.integer(.data$pos),
      ref = toupper(as.character(.data$ref)),
      alt = toupper(as.character(.data$alt))
    )
  bad <- h |>
    dplyr::summarise(lo = min(.data$freq), hi = max(.data$freq), .by = "haplotype") |>
    dplyr::filter(.data$lo != .data$hi | .data$lo < 0 | .data$hi > 1)
  if (nrow(bad)) {
    rlang::abort(glue::glue(gettext("`freq` must be one value in [0, 1] per haplotype.")))
  }
  h
}

#' Derive per-position frequencies and co-segregation from haplotypes
#'
#' Also reports which haplotypes carry the ALT allele at the target.
#' @keywords internal
#' @noRd
haplotype_variants <- function(haplotypes, snp) {
  h <- normalise_haplotypes(haplotypes)
  freqs <- h |>
    dplyr::distinct(.data$haplotype, .data$freq) |>
    tibble::deframe()
  all_hap <- names(freqs)
  alt_hap <- unique(h$haplotype[h$pos == snp$pos & h$alt == snp$alt])

  v <- h[h$pos != snp$pos, , drop = FALSE]
  v$w <- unname(freqs[v$haplotype])
  agg <- v |>
    dplyr::group_by(.data$pos, .data$ref, .data$alt) |>
    dplyr::summarise(
      pop_freq = min(1, sum(.data$w)),
      carriers = list(unique(.data$haplotype)),
      .groups = "drop"
    )
  variants <- tibble::tibble(
    chrom = snp$chrom, pos = agg$pos, ref = agg$ref, alt = agg$alt,
    pop_freq = agg$pop_freq, consequence = NA_character_
  )
  # a position co-segregates with the ALT target allele when exactly the ALT
  # carriers (and no REF carrier) have it: it rides the same haplotype
  linked <- NULL
  if (length(alt_hap) && nrow(agg)) {
    ride <- vapply(agg$carriers, function(c) setequal(c, alt_hap), logical(1))
    if (any(ride)) {
      linked <- tibble::tibble(target = snp$pos, linked = agg$pos[ride])
    }
  }
  list(
    variants = variants, linked = linked, haplotypes = h,
    alt = alt_hap, ref = setdiff(all_hap, alt_hap), freqs = freqs
  )
}

#' Which primer each haplotype is expected to bind
#'
#' For every haplotype: its allele at the target (REF/ALT), the AS primer that
#' matches it, and how many of the haplotype's own variants sit under that
#' primer's footprint (which would disrupt binding).  `binds` is `TRUE` when the
#' matching primer's footprint carries none of the haplotype's variants.
#' @keywords internal
#' @noRd
haplotype_hits <- function(haplotypes, snp, as_start, as_end, alt_hap = character(0)) {
  h <- normalise_haplotypes(haplotypes)
  freqs <- h |>
    dplyr::distinct(.data$haplotype, .data$freq) |>
    tibble::deframe()
  footprint <- IRanges::IRanges(as_start, as_end)
  purrr::imap_dfr(freqs, function(f, nm) {
    carries_alt <- nm %in% alt_hap
    own <- h$pos[h$haplotype == nm & h$pos != snp$pos]
    n_foot <- count_variants_in(footprint, own)
    tibble::tibble(
      haplotype = nm, freq = f,
      target_allele = if (carries_alt) snp$alt else snp$ref,
      primer = if (carries_alt) "AS_alt" else "AS_ref",
      footprint_variants = n_foot,
      binds = n_foot == 0L
    )
  })
}
