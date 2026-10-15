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
    cli::cli_abort(gettext(
      "`haplotypes` is missing column(s): {paste(missing, collapse = ', ')}."
    ))
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
    cli::cli_abort(gettext("`freq` must be one value in [0, 1] per haplotype."))
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
    dplyr::summarise(pop_freq = min(1, sum(.data$w)), .groups = "drop")
  variants <- tibble::tibble(
    chrom = snp$chrom, pos = agg$pos, ref = agg$ref, alt = agg$alt,
    pop_freq = agg$pop_freq, consequence = NA_character_
  )
  # each ALT haplotype's own co-segregating positions; `shared` marks the ones
  # carried by *every* ALT haplotype (the set a single ALT primer can bake)
  linked <- NULL
  if (length(alt_hap) && nrow(v)) {
    linked <- purrr::list_rbind(lapply(alt_hap, function(hp) {
      own <- v$pos[v$haplotype == hp]
      if (!length(own)) {
        return(NULL)
      }
      tibble::tibble(target = snp$pos, haplotype = hp, linked = own)
    }))
    if (!is.null(linked)) {
      carrier_sets <- lapply(alt_hap, function(hp) v$pos[v$haplotype == hp])
      linked$shared <- vapply(linked$linked, function(pos) {
        all(vapply(carrier_sets, function(cs) pos %in% cs, logical(1)))
      }, logical(1))
    }
  }
  list(
    variants = variants, linked = linked, haplotypes = h,
    alt = alt_hap, ref = setdiff(all_hap, alt_hap), freqs = freqs
  )
}

#' Group haplotypes by the alleles within a window around the target
#'
#' Two haplotypes that carry the *same* alleles around the target are the same
#' primer template: one ALT primer can bind both.  This groups a tidy haplotype
#' table by each haplotype's own variants within `window` bp of the target (either
#' side), so a designer designs one primer per *group* rather than one per
#' haplotype.  Chemistry-agnostic - it is the primer-template equivalence any
#' scheme (KASP, lax, indel, ...) can key on.
#' @param haplotypes a tidy haplotype table (`haplotype`, `freq`, `chrom`, `pos`,
#'   `ref`, `alt`).
#' @param snp one-row target (`pos`).
#' @param window half-width (bp) of the primer window around the target; `NULL`
#'   uses the longest KASP AS primer ([KASP_AS_LENGTHS]).
#' @return a tibble `haplotype`, `group` (the group's allele signature).
#' @examples
#' \dontrun{
#' haplotype_groups(haplotypes, snp)
#' }
#' @family KASP haplotypes
#' @export
haplotype_groups <- function(haplotypes, snp, window = NULL) {
  h <- normalise_haplotypes(haplotypes)
  window <- window %||% max(KASP_AS_LENGTHS)
  hp <- unique(h$haplotype)
  sig <- vapply(hp, function(x) {
    hap_signature(h, x, snp, snp$pos - window, snp$pos + window)
  }, character(1))
  tibble::tibble(haplotype = hp, group = sig)
}

#' A haplotype's own variants (`pos:ref:alt`) within a reference span
#'
#' The primer-template signature: two haplotypes with the same signature over a
#' primer's footprint are bound by the same primer.  The target base itself is
#' excluded (every ALT haplotype carries it).
#' @keywords internal
#' @noRd
hap_signature <- function(h, haplotype, snp, a, b) {
  d <- h[h$haplotype == haplotype & h$pos != snp$pos & h$pos >= a & h$pos <= b, , drop = FALSE]
  paste(sort(paste(d$pos, d$ref, d$alt, sep = ":")), collapse = ";")
}

#' Does a haplotype match a primer's template over a footprint?
#'
#' The one definition of "this haplotype is bound": its own variants over the
#' footprint equal the template's ([hap_signature()]), so a reference-template
#' primer (`template_sig = ""`) binds a haplotype only if it carries nothing
#' there.
#' @keywords internal
#' @noRd
binds_footprint <- function(h, haplotype, snp, a, b, template_sig) {
  identical(hap_signature(h, haplotype, snp, a, b), template_sig)
}

#' The template signature a baked linked set produces over a footprint
#'
#' The counterpart of [hap_signature()] for a primer that is *not* reference-based:
#' a primer carrying baked linked alleles matches only haplotypes whose footprint
#' signature equals theirs.  Coordinates are the haplotype table's own.
#' @keywords internal
#' @noRd
linked_template_sig <- function(h, positions, a, b) {
  p <- as.integer(positions)
  p <- p[!is.na(p) & p >= a & p <= b]
  if (!length(p)) {
    return("")
  }
  d <- h[h$pos %in% p, c("pos", "ref", "alt"), drop = FALSE]
  d <- d[!duplicated(d), , drop = FALSE]
  paste(sort(paste(d$pos, d$ref, d$alt, sep = ":")), collapse = ";")
}

#' The primer footprint of a designed marker (`NULL` when it has none)
#'
#' The span a marker's allele-specific primer covers, so the haplotype cover can
#' be merged by footprint.  Extend by adding methods for other chemistries.
#' @param marker a designed marker.
#' @param ... unused.
#' @return an integer vector `c(start, end)` in reference coordinates, or `NULL`.
#' @keywords internal
#' @noRd
marker_as_footprint <- function(marker, ...) UseMethod("marker_as_footprint")

#' @export
marker_as_footprint.default <- function(marker, ...) NULL

#' @export
marker_as_footprint.kasp_assay <- function(marker, ...) {
  aa <- marker$assay[1L, ]
  as.integer(c(aa$as_start, aa$as_end))
}

#' Design one ALT primer per primer-window group, merging by footprint
#'
#' The shared two-pass engine for **per-haplotype** design (used by
#' [design_lax_marker()], and available to any scheme whose ALT primer binds a
#' footprint).  First it groups the ALT haplotypes by the alleles in the primer
#' window ([haplotype_groups()]) and designs one marker per group with `design`;
#' then it merges by *footprint* - a design covers every haplotype whose variants
#' under its footprint match [binds_footprint()] - and greedily keeps the fewest
#' designs that cover every ALT haplotype.  So haplotypes that differ only
#' *outside* a primer share it.
#' @param haplotypes a tidy haplotype table.
#' @param snp one-row target.
#' @param design `function(haplotype)` returning a designed marker (or `NULL`);
#'   `haplotype` is the group representative's name.
#' @param window half-width (bp) of the primer window; `NULL` = the longest KASP
#'   AS primer (see [haplotype_groups()]).
#' @return a list: `groups` (tibble `haplotype`, `group`, `design` - one row per
#'   ALT haplotype, `design` `NULL` when nothing covers it), `designs` (the kept
#'   designs) and `uncovered` (haplotypes no design covers).
#' @examples
#' \dontrun{
#' haplotype_cover(haplotypes, snp, function(hp) design_lax_one(genome, snp, hp))
#' }
#' @family KASP haplotypes
#' @export
haplotype_cover <- function(haplotypes, snp, design, window = NULL) {
  h <- normalise_haplotypes(haplotypes)
  alt_hap <- unique(h$haplotype[h$pos == snp$pos & h$alt == snp$alt])
  empty <- tibble::tibble(haplotype = character(), group = character(), design = list())
  if (!length(alt_hap)) {
    return(list(groups = empty, designs = list(), uncovered = character()))
  }
  hg <- haplotype_groups(h, snp, window = window)
  hg <- hg[hg$haplotype %in% alt_hap, , drop = FALSE]
  reps <- split(hg$haplotype, hg$group)
  cand <- purrr::map(unname(reps), function(x) {
    d <- design(x[1L])
    fp <- marker_as_footprint(d)
    if (is.null(d) || is.null(fp)) NULL else list(rep = x[1L], design = d, fp = fp)
  })
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (!length(cand)) {
    return(list(
      groups = tibble::tibble(
        haplotype = alt_hap,
        group = NA_character_,
        design = vector("list", length(alt_hap))
      ),
      designs = list(), uncovered = alt_hap
    ))
  }
  cov <- purrr::map(cand, function(cc) {
    alt_hap[vapply(alt_hap, function(hp) {
      binds_footprint(h, hp, snp, cc$fp[1L], cc$fp[2L], hap_signature(h, cc$rep, snp, cc$fp[1L], cc$fp[2L]))
    }, logical(1))]
  })
  # greedy set cover: keep the fewest designs that cover every ALT haplotype
  left <- alt_hap
  pick <- integer(0)
  while (length(left)) {
    gain <- vapply(cov, function(x) length(intersect(x, left)), integer(1))
    if (all(gain == 0L)) break
    k <- which.max(gain)
    pick <- c(pick, k)
    left <- setdiff(left, cov[[k]])
  }
  sig <- vapply(cand, function(cc) hap_signature(h, cc$rep, snp, cc$fp[1L], cc$fp[2L]), character(1))
  which_k <- vapply(alt_hap, function(hp) {
    for (k in pick) {
      if (hp %in% cov[[k]]) {
        return(k)
      }
    }
    NA_integer_
  }, integer(1))
  list(
    groups = tibble::tibble(
      haplotype = alt_hap,
      group = dplyr::if_else(is.na(which_k), NA_character_, sig[which_k]),
      design = lapply(which_k, function(k) if (is.na(k)) NULL else cand[[k]]$design)
    ),
    designs = cand[pick],
    uncovered = alt_hap[is.na(which_k)]
  )
}

#' Which primer each haplotype is expected to bind
#'
#' For every haplotype: its allele at the target (REF/ALT), the AS primer that
#' matches it, and how many of the haplotype's own variants sit under that
#' primer's footprint *and are not baked into it* (which would disrupt binding).
#' `binds` is `TRUE` when the haplotype's footprint signature equals the primer's
#' template ([binds_footprint()]) - the one definition used everywhere: the REF
#' primer's template is the reference, and the ALT primer's is the reference plus
#' the `alt_baked` linked set.
#' @keywords internal
#' @noRd
haplotype_hits <- function(haplotypes, snp, as_start, as_end, alt_hap = character(0),
                           alt_start = NULL, alt_end = NULL, alt_baked = integer(0)) {
  h <- normalise_haplotypes(haplotypes)
  alt_start <- alt_start %||% as_start
  alt_end <- alt_end %||% as_end
  alt_sig <- linked_template_sig(h, alt_baked, alt_start, alt_end)
  freqs <- h |>
    dplyr::distinct(.data$haplotype, .data$freq) |>
    tibble::deframe()
  purrr::imap(freqs, function(f, nm) {
    carries_alt <- nm %in% alt_hap
    a <- if (carries_alt) alt_start else as_start
    b <- if (carries_alt) alt_end else as_end
    own <- h$pos[h$haplotype == nm & h$pos != snp$pos]
    if (carries_alt && length(alt_baked)) own <- setdiff(own, alt_baked)
    tibble::tibble(
      haplotype = nm, freq = f,
      target_allele = if (carries_alt) snp$alt else snp$ref,
      primer = if (carries_alt) "AS_alt" else "AS_ref",
      n_variants = length(own),
      footprint_variants = count_variants_in(IRanges::IRanges(a, b), own),
      binds = binds_footprint(h, nm, snp, a, b, if (carries_alt) alt_sig else "")
    )
  }) |>
    purrr::list_rbind()
}
