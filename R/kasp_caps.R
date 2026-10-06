# =============================================================================
#  kasp_caps.R  --  CAPS / dCAPS markers (the PCR family beyond KASP)
# =============================================================================
#  Where KASP needs the variant at an allele-specific primer's 3' end, a CAPS
#  marker anchors the variant *inside the amplicon* and reads it with a
#  restriction digest: a biallelic site whose recognition sequence is created or
#  destroyed by the variant (or introduced by a primer, dCAPS) yields a
#  different set of fragment lengths for each allele - a binary, gel-scorable
#  call.  That is why CAPS tolerates the very things that doom a KASP: a variant
#  in the 3' window, an AT-rich junction, or marginal discrimination (see
#  ROADMAP section U).  The restriction-reaction recipe (enzyme + expected
#  fragments) is stored on the object, like the KASP mix recipe is.
# =============================================================================

#' Bundled restriction-enzyme catalogue (name, site, cut)
#'
#' The recognition site (IUPAC) and the top-strand cut position (bases from the
#' 5' end of the site) for **palindromic** Type II restriction enzymes, sourced
#' from the REBASE *bionet* catalogue and regenerated deterministically by
#' `dev/fetch-rebase.R` (never typed by hand).  Palindromic sites are read on the
#' plus strand alone, so one match per site is enough.  Non-palindromic and
#' Type IIS enzymes (which cut outside/asymmetrically) are deliberately excluded -
#' the fragment model here takes a single top-strand offset.
#' @return tibble with `enzyme`, `site`, `cut`.
#' @examples
#' head(restriction_enzymes())
#' @family KASP CAPS
#' @export
restriction_enzymes <- function() {
  f <- system.file("extdata", "restriction_enzymes.csv", package = "tidykasp")
  if (!nzchar(f)) cli::cli_abort(gettext("Bundled restriction-enzyme catalogue not found."))
  readr::read_csv(f, show_col_types = FALSE)
}

#' Digest fragment lengths of one sequence, per enzyme
#'
#' For each enzyme, finds the palindromic recognition sites (IUPAC-aware) and
#' returns the top-strand cut boundaries and the resulting linear fragment
#' lengths.  An enzyme with no site returns a single band (no cut).
#' @param seq a DNAString or character sequence (the amplicon, plus strand).
#' @param enzymes enzyme table (see [restriction_enzymes()]).
#' @return tibble with `enzyme`, `n_sites`, `cuts` (list of 0-based boundaries)
#'   and `fragments` (list of integer band lengths).
#' @examples
#' caps_fragments("GAATTCTTAAGGAATTC")
#' @family KASP CAPS
#' @export
caps_fragments <- function(seq, enzymes = restriction_enzymes()) {
  s <- Biostrings::DNAString(as.character(seq))
  n <- length(s)
  purrr::pmap(enzymes, function(enzyme, site, cut) {
    hits <- Biostrings::matchPattern(site, s, fixed = FALSE)
    starts <- IRanges::start(hits) - 1L
    cut_at <- sort(unique(starts + as.integer(cut)))
    cut_at <- cut_at[cut_at > 0L & cut_at < n]
    frag <- as.integer(diff(c(0L, cut_at, n)))
    tibble::tibble(
      enzyme = enzyme, n_sites = length(starts),
      cuts = list(cut_at), fragments = list(frag)
    )
  }) |>
    purrr::list_rbind()
}

#' The fragment-length gap between two digests (how resolvable on a gel)
#'
#' The smallest absolute difference between a REF band and an ALT band (ignoring
#' equal pairs) - the closest two bands that a gel would have to separate.  `Inf`
#' when the two digest patterns are identical.
#' @keywords internal
#' @noRd
digest_gap <- function(a, b) {
  if (length(a) == length(b) && setequal(a, b)) {
    return(Inf)
  }
  d <- abs(outer(a, b, "-"))
  d <- d[d > 0]
  if (length(d)) min(d) else Inf
}

#' Enzymes whose digest distinguishes the REF and ALT amplicons (CAPS scan)
#'
#' Copies the target into a window, builds the two allele amplicons (REF vs ALT),
#' digests each with every enzyme and keeps the enzymes whose fragment-length
#' pattern differs - i.e. a natural CAPS marker.  A pattern is reported as
#' *resolvable* when the closest differing bands are at least `min_gap` apart
#' (`min_gap = 0` keeps everything that differs at all).
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank window half-width (the prospective amplicon).
#' @param enzymes enzyme table (see [restriction_enzymes()]).
#' @param min_gap minimal band separation to call the digest gel-resolvable.
#' @return tidy tibble, one row per enzyme with `n_sites_ref`, `n_sites_alt`,
#'   `fragments_ref`, `fragments_alt`, `site_change` (a recognition site is
#'   gained/lost, not just a length shift), `gap` and `resolvable`, sorted by
#'   decreasing `resolvable` then `gap`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' snp <- kasp_example("SNP")[1, ]
#' find_caps_enzymes(genome, snp, flank = 250)
#' }
#' @family KASP CAPS
#' @export
find_caps_enzymes <- function(genome, snp, flank = 250L, enzymes = restriction_enzymes(),
                              min_gap = 50L) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  ref_amp <- tmpl$seq
  ref_len <- nchar(as.character(snp$ref)[1L])
  alt_amp <- Biostrings::replaceAt(
    ref_amp, IRanges::IRanges(tmpl$snp_i, width = ref_len), as.character(snp$alt)[1L]
  )
  fr <- caps_fragments(ref_amp, enzymes)
  fa <- caps_fragments(alt_amp, enzymes)
  out <- dplyr::inner_join(
    dplyr::rename(fr, n_sites_ref = "n_sites", cuts_ref = "cuts", fragments_ref = "fragments"),
    dplyr::rename(fa, n_sites_alt = "n_sites", cuts_alt = "cuts", fragments_alt = "fragments"),
    by = "enzyme"
  )
  # a CAPS call is a *site* change (gain/loss of a recognition site), not the
  # trivial length difference an indel leaves in every enzyme's digest
  out$site_change <- out$n_sites_ref != out$n_sites_alt
  out$gap <- purrr::map2_dbl(out$fragments_ref, out$fragments_alt, digest_gap)
  out$resolvable <- out$site_change & is.finite(out$gap) & out$gap >= min_gap
  out <- dplyr::filter(out, is.finite(.data$gap))
  dplyr::arrange(out, dplyr::desc(.data$resolvable), dplyr::desc(.data$gap))
}

#' IUPAC-aware mismatch offsets between a sequence window and an enzyme site
#'
#' A position is a mismatch when the base is not one of the site's IUPAC code
#' members.  Returns 0-based offsets.
#' @keywords internal
#' @noRd
site_mismatches <- function(window, site) {
  w <- strsplit(window, "")[[1L]]
  s <- strsplit(site, "")[[1L]]
  k <- min(length(w), length(s))
  w <- w[seq_len(k)]
  s <- s[seq_len(k)]
  bad <- which(!vapply(seq_len(k), function(i) {
    grepl(w[i], Biostrings::IUPAC_CODE_MAP[[s[i]]], fixed = TRUE)
  }, logical(1)))
  bad - 1L
}

#' dCAPS: an enzyme the *primer* can turn on for one allele
#'
#' A natural CAPS site is absent at most variants, so dCAPS introduces one: a
#' primer carries a single mismatch that completes the recognition site next to
#' the variant, and the variant base itself is inside the site - so the site forms
#' for one allele only.  For every enzyme this scans the windows that contain the
#' variant, keeps those that are one (non-variant) mismatch short of the site for
#' one allele while the other allele breaks it, applies the primer's base to both
#' amplicons and digests.  Single-base substitutions only (a length-changing
#' allele cannot be completed by a substitution primer).
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank window half-width.
#' @param enzymes enzyme table (see [restriction_enzymes()]).
#' @param min_gap minimal band separation to call it gel-resolvable.
#' @return tidy tibble, one row per candidate: `enzyme`, `allele` (which allele the
#'   primer completes the site for), `site_start` (template), `site_len`,
#'   `mm_offset` (0-based, the base the primer introduces), `fix_base`, the two
#'   digest `fragments_*`, `gap` and `resolvable`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' find_dcaps_enzymes(genome, kasp_example("SNP")[1, ], flank = 250)
#' }
#' @family KASP CAPS
#' @export
find_dcaps_enzymes <- function(genome, snp, flank = 250L, enzymes = restriction_enzymes(),
                               min_gap = 50L) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    cli::cli_inform(gettext("dCAPS scanning is implemented for single-base substitutions only."))
    return(tibble::tibble())
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  ref_chr <- as.character(tmpl$seq)
  snp_i <- tmpl$snp_i
  n <- nchar(ref_chr)
  alt_chr <- as.character(Biostrings::replaceAt(
    tmpl$seq, IRanges::IRanges(snp_i, width = 1L), alt
  ))
  purrr::pmap(enzymes, function(enzyme, site, cut) {
    L <- nchar(site)
    starts <- seq.int(max(1L, snp_i - L + 1L), min(n - L + 1L, snp_i))
    if (!length(starts)) {
      return(NULL)
    }
    purrr::map(starts, function(st) {
      rw <- substr(ref_chr, st, st + L - 1L)
      aw <- substr(alt_chr, st, st + L - 1L)
      voff <- snp_i - st
      mr <- site_mismatches(rw, site)
      ma <- site_mismatches(aw, site)
      allele <- if (length(mr) == 1L && mr != voff && length(ma) >= 2L) {
        "ref"
      } else if (length(ma) == 1L && ma != voff && length(mr) >= 2L) {
        "alt"
      } else {
        NA_character_
      }
      if (is.na(allele)) {
        return(NULL)
      }
      mm <- if (allele == "ref") mr else ma
      fix <- substr(site, mm + 1L, mm + 1L)
      pos <- st + mm # 1-based template position the primer corrects
      rf <- ref_chr
      af <- alt_chr
      substr(rf, pos, pos) <- fix
      substr(af, pos, pos) <- fix
      one <- enzymes[enzymes$enzyme == enzyme, ]
      fr <- caps_fragments(rf, one)$fragments[[1L]]
      fa <- caps_fragments(af, one)$fragments[[1L]]
      tibble::tibble(
        enzyme = enzyme, allele = allele, site_start = st, site_len = L,
        mm_offset = mm, fix_base = fix,
        fragments_ref = list(fr), fragments_alt = list(fa), gap = digest_gap(fr, fa)
      )
    }) |>
      purrr::list_rbind()
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(resolvable = .data$gap >= min_gap) |>
    dplyr::arrange(dplyr::desc(.data$gap))
}

#' Design a CAPS / dCAPS marker at a variant
#'
#' The CAPS analogue of [assemble_kasp_assay()]: find the enzyme whose digest
#' separates the alleles (see [find_caps_enzymes()]), then design the flanking
#' primer pair around the discriminating site.  When no natural site exists it
#' falls back to dCAPS ([find_dcaps_enzymes()]): a primer introduces the site.
#' Returns a `kasp_caps` object with the enzyme, both digest patterns and the
#' primer pair - the restriction recipe travels with the marker, as the KASP mix
#' recipe does.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank window half-width for the scan and primer design.
#' @param enzymes enzyme table (see [restriction_enzymes()]).
#' @param min_gap minimal band separation to require (see [find_caps_enzymes()]).
#' @return a `kasp_caps` list (`snp`, `enzyme`, `fragments_ref`, `fragments_alt`,
#'   `gap`, `primers`, `amplicon`), or `NULL` (with a message) when no enzyme
#'   separates the alleles.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_caps_marker(genome, kasp_example("SNP")[1, ], flank = 250)
#' }
#' @family KASP CAPS
#' @export
design_caps_marker <- function(genome, snp, flank = 250L, enzymes = restriction_enzymes(),
                               min_gap = 50L) {
  hits <- find_caps_enzymes(genome, snp, flank, enzymes, min_gap)
  hits <- dplyr::filter(hits, .data$resolvable)
  if (!nrow(hits)) {
    return(design_dcaps_marker(genome, snp, flank, enzymes, min_gap))
  }
  best <- hits[1L, ]
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  # the discriminating cut(s): where the two digest patterns first diverge
  cuts <- setdiff(best$cuts_ref[[1L]], best$cuts_alt[[1L]])
  if (!length(cuts)) cuts <- setdiff(best$cuts_alt[[1L]], best$cuts_ref[[1L]])
  lo <- min(cuts) + 1L
  hi <- max(cuts) + 1L
  primers <- design_site_pair(tmpl$seq, lo, hi)
  if (is.null(primers)) {
    cli::cli_inform(gettext("A discriminating enzyme was found, but no flanking primer pair passed the Tm/GC filters."))
    return(NULL)
  }
  out <- list(
    snp = snp, enzyme = best$enzyme,
    fragments_ref = best$fragments_ref[[1L]], fragments_alt = best$fragments_alt[[1L]],
    gap = best$gap, primers = primers,
    amplicon = tibble::tibble(
      start = snp$pos + primers$start - tmpl$snp_i,
      end = snp$pos + primers$end - tmpl$snp_i
    )
  )
  structure(out, class = "kasp_caps")
}

#' Build a dCAPS marker from a primer-introduced site
#'
#' Takes the best [find_dcaps_enzymes()] candidate, applies the primer's base to
#' the REF amplicon, designs the site-carrying primer + a flanking mate, and
#' packages the restricted recipe as a `kasp_caps` object (`dcaps = TRUE`).
#' @keywords internal
#' @noRd
design_dcaps_marker <- function(genome, snp, flank, enzymes, min_gap) {
  d <- find_dcaps_enzymes(genome, snp, flank, enzymes, min_gap)
  d <- dplyr::filter(d, .data$resolvable)
  if (!nrow(d)) {
    cli::cli_inform(gettext("No natural or primer-introduced (dCAPS) restriction marker here."))
    return(NULL)
  }
  best <- d[1L, ]
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  ref_chr <- as.character(tmpl$seq)
  pos <- best$site_start + best$mm_offset # 1-based; the base the primer introduces
  fixed <- ref_chr
  substr(fixed, pos, pos) <- best$fix_base
  lo <- best$site_start
  hi <- best$site_start + best$site_len - 1L
  primers <- design_site_pair(fixed, lo, hi)
  if (is.null(primers)) {
    cli::cli_inform(gettext("A dCAPS site was found, but no primer pair passed the Tm/GC filters."))
    return(NULL)
  }
  out <- list(
    snp = snp, enzyme = best$enzyme, dcaps = TRUE,
    fragments_ref = best$fragments_ref[[1L]], fragments_alt = best$fragments_alt[[1L]],
    gap = best$gap, primers = primers,
    amplicon = tibble::tibble(
      start = snp$pos + primers$start - tmpl$snp_i,
      end = snp$pos + primers$end - tmpl$snp_i
    )
  )
  structure(out, class = "kasp_caps")
}

#' @export
print.kasp_caps <- function(x, ...) {
  cat("<kasp_caps>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, " via ", x$enzyme,
    if (isTRUE(x$dcaps)) " (dCAPS)" else "", "\n",
    sep = ""
  )
  cat("  REF digest:", paste(x$fragments_ref, collapse = " + "), "bp\n")
  cat("  ALT digest:", paste(x$fragments_alt, collapse = " + "), "bp\n")
  cat("  amplicon:", x$primers$product, "bp (gap", x$gap, "bp)\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  invisible(x)
}

#' Draw the expected CAPS / dCAPS gel
#'
#' A two-lane schematic gel (REF and ALT) with a band for every restriction
#' fragment length, the diagnostic bands highlighted - so the readout is visible
#' without running the digest.  Bands shared by both alleles are grey;
#' allele-diagnostic bands take the allele colour.
#' @param marker a `kasp_caps` object (see [design_caps_marker()]).
#' @return a `ggplot`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' m <- design_caps_marker(genome, kasp_example("SNP")[1, ], flank = 250)
#' plot_caps_gel(m)
#' }
#' @family KASP CAPS
#' @export
plot_caps_gel <- function(marker) {
  rlang::check_installed("ggplot2", reason = gettext("to draw the CAPS gel"))
  if (!inherits(marker, "kasp_caps")) {
    cli::cli_abort(gettext("`marker` must be a kasp_caps object (see design_caps_marker())."))
  }
  shared <- intersect(marker$fragments_ref, marker$fragments_alt)
  bands <- tibble::tibble(
    allele = rep(c("REF", "ALT"), times = c(length(marker$fragments_ref), length(marker$fragments_alt))),
    size = c(marker$fragments_ref, marker$fragments_alt)
  ) |>
    dplyr::mutate(
      diagnostic = !.data$size %in% shared,
      x = as.numeric(factor(.data$allele, levels = c("REF", "ALT"))),
      col = dplyr::case_when(
        !.data$diagnostic ~ "grey65",
        .data$allele == "REF" ~ KASP_STATUS[["neutral"]],
        TRUE ~ KASP_STATUS[["broken"]]
      )
    )
  ref_txt <- paste(marker$fragments_ref, collapse = "+")
  alt_txt <- paste(marker$fragments_alt, collapse = "+")
  gel_plot(
    bands, c("REF", "ALT"), "fragment length (bp)",
    glue::glue("CAPS gel: {marker$snp$chrom}:{marker$snp$pos} via {marker$enzyme}"),
    glue::glue("REF: {ref_txt}  |  ALT: {alt_txt} bp")
  )
}
