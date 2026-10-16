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
#
#  Design rules (researched, not invented):
#   * The assay is a co-dominant *restriction* readout - the variant creates or
#     destroys a recognition site, so the two alleles give different fragment
#     sets (Konieczny & Ausubel 1993, Plant J 4:403-410,
#     <https://doi.org/10.1046/j.1365-313X.1993.04020403.x>, the original CAPS
#     paper).  No discriminating enzyme = no marker, so that gate is hard.
#   * Most SNPs create no natural site, so a *primer* introduces one (dCAPS; Neff
#     et al. 1998, Plant J 14:387-392, <https://doi.org/10.1046/j.1365-313X.1998.00124.x>)
#     by carrying one mismatch that completes the site next to the variant; the
#     variant base itself sits inside the site, so it forms for one allele only.
#   * The amplicon must hold the site and leave two bands a plain gel resolves:
#     the widely used rule is a >= 50 bp difference for mid-size fragments on
#     1-2 % agarose (smaller differences need high-percentage agarose or PAGE -
#     hence the relaxed floor).  That makes the CAPS amplicon longer (120-1200
#     bp) than the 50-100 bp KASP optimum; see `caps_limits()`.
#   * A failed digest mimics the uncut ("REF") allele, so the readout is only
#     unambiguous against an undigested control lane or a second, constant cut -
#     a wet-lab step the recipe (`plot_caps_gel()`) makes visible.
# =============================================================================

# ---- the CAPS selection stage -------------------------------------------------
# CAPS / dCAPS: the variant is read by a restriction digest (the designer below),
# so selection is "is there an enzyme that separates the alleles, and do the
# resulting bands resolve on a gel" - none of the KASP 3'-terminal machinery
# applies.  The two declared orders (C6): gates applied in the order below
# (grave / cheap first), then the ladder (least grave first: the amplicon window
# is a preference, the Tm balance a practical constraint, and the band gap the
# very point of the assay - still bounded below).
KASP_GATES_CAPS <- c("gate_caps_site", "gate_caps_resolvable", "gate_caps_amplicon", "gate_caps_tm")
KASP_RELAX_CAPS <- c("gate_caps_amplicon", "gate_caps_tm", "gate_caps_resolvable")
KASP_TIERS_CAPS <- c(
  tier1 = "band gap (gel resolvability)", tier2 = "Tm balance", tier3 = "amplicon length"
)

KASP_CAPS_GAP_MIN <- 20L # bp: the accepted closest-band difference (a resolved gel)
KASP_CAPS_GAP_FLOOR <- 10L # bp: the band gap may be relaxed to a PAGE-resolvable 10 bp

# The window the CAPS designer selects against (and the gate limits display).  A
# CAPS amplicon is read by *fragment length*, so - unlike a KASP amplicon
# (`KASP_LIMITS$amplicon`, the 50-100 bp EasyKASP optimum) - it must hold the
# discriminating site and leave two gel-resolvable fragments: 120-1200 bp.  A
# function, not a constant, because R sources R/kasp_core.R (where KASP_LIMITS
# lives) *after* this file ("caps" < "core"), like `caps_gate_limit()` below.
#
# `band_gap` is the accepted closest-band difference measured on the *product*.
# It is 20 bp, not the 50 bp the scanning window would suggest, because a dCAPS
# marker's cut sits inside the primer that introduces the site, so its
# REF-against-ALT difference is bounded by the primer length (~20-28 bp) - still
# a scorable 2-4 % agarose / PAGE gel.  Natural sites, whose pair can centre the
# cut, clear the bar by a wide margin and win the tier.
#
# NB (C5, requested serial slot): the acceptance check `amplicon length (bp)`
# carries the KASP window (KASP_LIMITS$amplicon, 50-100 bp) for *every* scheme, so
# a wet-lab-viable CAPS amplicon would be reported as an acceptance failure.  The
# check has to become scheme-aware in kasp_criteria()/kasp_checks.R - files this
# agent does not own - so the reconciliation is reported as a requested serial
# slot.  Until it lands, the CAPS window is the gate below, not KASP_LIMITS.
caps_limits <- function() {
  list(
    amplicon = c(120, 1200), # bp
    band_gap = KASP_CAPS_GAP_MIN, # bp: closest differing bands on the product
    tm = KASP_LIMITS$as_tm # C: both flanking primers (design_site_pair's window)
  )
}

# How far each CAPS gate may be lowered (margin units; 0 = hard) - the ladder's
# floors (C8).  `gate_caps_site` is hard (no discriminating enzyme = no CAPS
# marker).  The band gap is soft but *bounded*: it may be relaxed only to a 10 bp
# difference (`KASP_CAPS_GAP_FLOOR`), still resolvable on a high-percentage gel /
# PAGE, so a relaxed winner is never an unresolvable digest - the risk the audit
# flagged.  The amplicon window widens to [100, 1220] and the Tm window to
# [53, 67] C.
KASP_GATE_FLOOR_CAPS <- c(
  gate_caps_site = 0,
  gate_caps_resolvable = KASP_CAPS_GAP_FLOOR - KASP_CAPS_GAP_MIN, # 20 -> 10 bp
  gate_caps_amplicon = -20, # [120, 1200] -> [100, 1220]
  gate_caps_tm = -2 # [55, 65] -> [53, 67]
)

# This scheme's entries in the shared gate-label lookups (assembled in
# R/kasp_select.R), kept beside the gates they name.  The limit labels are built
# by a function, not a constant, because the windows come from `caps_limits()`
# and `KASP_LIMITS` (defined *after* this file).  Real values, never a `[lo, hi]`
# placeholder (C9).
KASP_GATE_LABEL_CAPS <- c(
  gate_caps_site = "no enzyme separates the REF and ALT amplicons",
  gate_caps_resolvable = "the diagnostic bands are too close on a gel",
  gate_caps_amplicon = "the CAPS amplicon is outside the accepted window",
  gate_caps_tm = "a flanking primer Tm is outside the window"
)
KASP_GATE_SHORT_CAPS <- c(
  gate_caps_site = "restriction site", gate_caps_resolvable = "band gap",
  gate_caps_amplicon = "amplicon length (CAPS)", gate_caps_tm = "primer Tm (CAPS)"
)
caps_gate_limit <- function() {
  lim <- caps_limits()
  c(
    gate_caps_site = "an enzyme's cut separates the REF and ALT amplicons",
    gate_caps_resolvable = glue::glue("closest product bands >= {lim$band_gap} bp apart"),
    gate_caps_amplicon = glue::glue("amplicon within [{lim$amplicon[1L]}, {lim$amplicon[2L]}] bp"),
    gate_caps_tm = glue::glue(
      "both primers Tm within [{lim$tm[1L]}, {lim$tm[2L]}] C"
    )
  )
}
KASP_GATE_FLOOR_LABEL_CAPS <- c(
  gate_caps_amplicon = "widen the amplicon window",
  gate_caps_resolvable = "allow a 10 bp band gap (PAGE)",
  gate_caps_tm = "widen the Tm window"
)
KASP_GATE_HARD_LABEL_CAPS <- c(
  gate_caps_site = "no discriminating enzyme = no CAPS marker"
)

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
  found <- purrr::pmap(enzymes, function(enzyme, site, cut) {
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
    purrr::list_rbind()
  if (is.null(found) || !nrow(found)) {
    # `enzymes` carries plain columns (site/cut, not list columns), so pmap hands
    # over values - only the *result* can be empty, and then there is nothing to
    # mark resolvable
    return(tibble::tibble())
  }
  found |>
    dplyr::mutate(resolvable = .data$gap >= min_gap) |>
    dplyr::arrange(dplyr::desc(.data$gap))
}

# ---- the candidates the CAPS designer selects over ----------------------------
# `find_caps_enzymes()`/`find_dcaps_enzymes()` enumerate the *sites*; turning a
# site into a flanking pair and choosing among them is the selector's job (C7).
# These helpers keep ONE candidate schema for both site kinds (natural and
# primer-introduced), so every gate reads the same columns and the designer stays
# a thin wrapper around `select_with_ladder()`.

#' One candidate row-set (natural and dCAPS share this schema)
#' @keywords internal
#' @noRd
caps_candidates <- function(dcaps, enzyme, allele, site_change, gap, reg_lo, reg_hi,
                            cuts_ref, cuts_alt, fix_base, fix_at) {
  tibble::tibble(
    dcaps = dcaps, enzyme = enzyme, allele = allele, site_change = site_change,
    gap = gap, reg_lo = reg_lo, reg_hi = reg_hi,
    cuts_ref = cuts_ref, cuts_alt = cuts_alt,
    fix_base = fix_base, fix_at = fix_at
  )
}

#' The discriminating region a CAPS pair must span, from the two digest cut sets
#' @keywords internal
#' @noRd
caps_site_region <- function(cuts_ref, cuts_alt) {
  cuts <- setdiff(cuts_ref, cuts_alt)
  if (!length(cuts)) cuts <- setdiff(cuts_alt, cuts_ref)
  if (!length(cuts)) {
    return(NULL)
  }
  c(min(cuts) + 1L, max(cuts) + 1L)
}

#' Natural CAPS sites as candidates, down to the ladder's floor
#'
#' Not filtered by `resolvable`: the site-change gate, not a pre-filter, decides,
#' so it has something to reject (C7).  Ranked by site change then band gap, both
#' criteria the selector prizes, and truncated so the primer design runs on the
#' handful that could win.  `gap` here is the *window* gap, used only to
#' shortlist; the real gap is recomputed on the designed amplicon (see
#' `caps_grade()`).
#' @keywords internal
#' @noRd
caps_natural_candidates <- function(genome, snp, flank, enzymes, min_gap, n_max = 6L) {
  hits <- find_caps_enzymes(genome, snp, flank, enzymes, min_gap)
  hits <- dplyr::filter(hits, is.finite(.data$gap), .data$gap >= min_gap)
  hits <- utils::head(
    dplyr::arrange(hits, dplyr::desc(.data$site_change), dplyr::desc(.data$gap)), n_max
  )
  caps_candidates(
    dcaps = rep(FALSE, nrow(hits)), enzyme = hits$enzyme,
    allele = rep(NA_character_, nrow(hits)), site_change = hits$site_change,
    gap = hits$gap,
    reg_lo = rep(NA_integer_, nrow(hits)), reg_hi = rep(NA_integer_, nrow(hits)),
    cuts_ref = hits$cuts_ref, cuts_alt = hits$cuts_alt,
    fix_base = rep(NA_character_, nrow(hits)), fix_at = rep(NA_integer_, nrow(hits))
  )
}

#' dCAPS sites (a primer completes the site) as candidates
#' @keywords internal
#' @noRd
caps_dcaps_candidates <- function(genome, snp, flank, enzymes, min_gap, n_max = 6L) {
  d <- find_dcaps_enzymes(genome, snp, flank, enzymes, min_gap)
  if (!nrow(d)) {
    return(caps_candidates(
      dcaps = logical(0), enzyme = character(0), allele = character(0),
      site_change = logical(0), gap = numeric(0),
      reg_lo = integer(0), reg_hi = integer(0), cuts_ref = list(), cuts_alt = list(),
      fix_base = character(0), fix_at = integer(0)
    ))
  }
  d <- utils::head(dplyr::arrange(d, dplyr::desc(.data$gap)), n_max)
  caps_candidates(
    dcaps = rep(TRUE, nrow(d)), enzyme = d$enzyme, allele = d$allele,
    site_change = rep(TRUE, nrow(d)), gap = d$gap,
    reg_lo = d$site_start, reg_hi = d$site_start + d$site_len - 1L,
    cuts_ref = rep(list(integer(0)), nrow(d)), cuts_alt = rep(list(integer(0)), nrow(d)),
    fix_base = d$fix_base, fix_at = d$site_start + d$mm_offset
  )
}

#' Grade candidate sites: the flanking pair plus the gate margins
#'
#' Turns one row per candidate *site* into the pair that would read it and the
#' gate margins `select_with_ladder()` needs.  The accepted window is searched
#' first, so a pairing that fits it is never drawn as a relaxation; only when no
#' pair fits is the window widened to the gate's floor, which is what makes the
#' amplicon and Tm gates real decisions (C7) rather than pre-filters.  The hard
#' site gate rejects a length-only candidate.
#'
#' The pair is placed the way the method calls for: a *natural* site gets a plain
#' flanking pair, which centres the cut and so leaves a well-separated band pair;
#' a *dCAPS* site needs the primer that introduces it (`design_site_pair()`), so
#' its difference is bounded by the primer length - exactly why the gate's bar is
#' 20 bp.  The band gap is measured on the *designed product*, not the scanning
#' window: REF and ALT amplicons carry the (dCAPS) primer's base and differ only
#' at the variant, and only their digest is the gel a wet-lab run would show.
#' @keywords internal
#' @noRd
caps_grade <- function(cands, snp, tmpl, enzymes, limits, floors) {
  if (!nrow(cands)) {
    return(cands)
  }
  amp_win <- limits$amplicon
  tm_win <- limits$tm
  amp_relax <- amp_win + c(-1, 1) * (-floors[["gate_caps_amplicon"]])
  tm_relax <- tm_win + c(-1, 1) * (-floors[["gate_caps_tm"]])
  design_pair <- function(seq, lo, hi, dcaps) {
    f <- if (dcaps) design_site_pair else design_flanking_pair
    f(seq, lo, hi, min_product = amp_win[1L], max_product = amp_win[2L], tm_range = tm_win) %||%
      f(seq, lo, hi, min_product = amp_relax[1L], max_product = amp_relax[2L], tm_range = tm_relax)
  }
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  delta <- nchar(alt) - nchar(ref)
  ref_chr <- as.character(tmpl$seq)
  variant_i <- tmpl$snp_i
  # by row, not pmap: the cut sets are list-columns, and `cand$cuts_ref[[1L]]` is
  # unambiguous where pmap would hand over a one-element list
  purrr::map(seq_len(nrow(cands)), function(k) {
    cand <- cands[k, , drop = FALSE]
    seq <- ref_chr
    if (cand$dcaps) {
      substr(seq, cand$fix_at, cand$fix_at) <- cand$fix_base
    }
    reg <- if (cand$dcaps) c(cand$reg_lo, cand$reg_hi) else caps_site_region(cand$cuts_ref[[1L]], cand$cuts_alt[[1L]])
    if (is.null(reg)) {
      return(NULL)
    }
    if (!cand$dcaps) {
      # a natural site's pair must *centre* the discriminating cut(s): expand the
      # region by the accepted gap on each side before flanking it, so the
      # product's two bands come out at least that far apart.  A dCAPS site
      # cannot - its cut sits inside the very primer that introduces it, which is
      # why its difference is bounded by the primer length (see the gate's bar).
      m <- limits$band_gap
      reg <- c(max(2L, reg[1L] - m), min(nchar(seq) - 1L, reg[2L] + m))
      if (reg[1L] >= reg[2L]) {
        return(NULL)
      }
    }
    pair <- design_pair(seq, reg[1L], reg[2L], cand$dcaps)
    if (is.null(pair)) {
      return(NULL)
    }
    alt_seq <- paste0(substr(seq, 1L, variant_i - 1L), alt, substr(seq, variant_i + nchar(ref), nchar(seq)))
    one <- enzymes[enzymes$enzyme == cand$enzyme, ]
    fr <- caps_fragments(substr(seq, pair$start, pair$end), one)$fragments[[1L]]
    fa <- caps_fragments(substr(alt_seq, pair$start, pair$end + delta), one)$fragments[[1L]]
    gap <- digest_gap(fr, fa)
    if (!is.finite(gap)) {
      return(NULL) # the shortlisted cut sits outside the product: not a CAPS call
    }
    m_tm <- min(pair$tm_f - tm_win[1L], tm_win[2L] - pair$tm_f) |>
      min(pair$tm_r - tm_win[1L], tm_win[2L] - pair$tm_r)
    m_amp <- min(pair$product - amp_win[1L], amp_win[2L] - pair$product)
    tibble::tibble(
      dcaps = cand$dcaps, enzyme = cand$enzyme, allele = cand$allele, gap = gap,
      fragments_ref = list(fr), fragments_alt = list(fa), primers = list(pair),
      m_site = if (isTRUE(cand$site_change)) 1 else -1,
      m_resolvable = gap - limits$band_gap,
      m_amplicon = m_amp, m_tm = m_tm,
      tier1 = gap - limits$band_gap, tier2 = m_tm, tier3 = m_amp
    )
  }) |>
    purrr::list_rbind()
}

#' The best CAPS/dCAPS marker the gates admit, or `NULL`
#'
#' The designer's whole selection step: grade the candidates, run the ladder
#' (hard `gate_caps_site`, soft bounded amplicon / Tm / band-gap), and package the
#' winner as a `kasp_caps` object carrying the restriction recipe, the pair and
#' the gate outcome (`$selection`) so the DAG is evaluated, not decorative (C7).
#' @keywords internal
#' @noRd
caps_marker <- function(cands, snp, tmpl, enzymes, limits, floors) {
  pick <- select_with_ladder(
    caps_grade(cands, snp, tmpl, enzymes, limits, floors),
    gates = KASP_GATES_CAPS, relax_order = KASP_RELAX_CAPS,
    tiers = names(KASP_TIERS_CAPS), ties = character(0),
    gate_margin = c(
      gate_caps_site = "m_site", gate_caps_resolvable = "m_resolvable",
      gate_caps_amplicon = "m_amplicon", gate_caps_tm = "m_tm"
    ),
    gate_floor = floors
  )
  best <- pick$best
  if (is.null(best)) {
    return(NULL)
  }
  primers <- best$primers[[1L]]
  warn_relaxed(pick$violated)
  structure(
    list(
      snp = snp, enzyme = best$enzyme, dcaps = isTRUE(best$dcaps), allele = best$allele,
      fragments_ref = best$fragments_ref[[1L]], fragments_alt = best$fragments_alt[[1L]],
      gap = best$gap, primers = primers,
      amplicon = tibble::tibble(
        start = snp$pos + primers$start - tmpl$snp_i,
        end = snp$pos + primers$end - tmpl$snp_i
      ),
      selection = list(relaxed = pick$relaxed, violated = pick$violated)
    ),
    class = "kasp_caps"
  )
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
  # `min_gap` is the cheap *scanning* pre-filter (see find_caps_enzymes()); the
  # gate that decides is the band gap measured on the real product (`caps_limits()`)
  limits <- caps_limits()
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  marker <- caps_marker(
    caps_natural_candidates(genome, snp, flank, enzymes, min_gap),
    snp, tmpl, enzymes, limits, KASP_GATE_FLOOR_CAPS
  )
  if (!is.null(marker)) {
    return(marker)
  }
  # no natural site separates the alleles - let a primer introduce one (dCAPS)
  design_dcaps_marker(genome, snp, flank, enzymes, min_gap)
}

#' Build a dCAPS marker from a primer-introduced site
#'
#' Takes the best [find_dcaps_enzymes()] candidate, applies the primer's base to
#' the REF amplicon, designs the site-carrying primer + a flanking mate, and
#' packages the restricted recipe as a `kasp_caps` object (`dcaps = TRUE`).
#' @keywords internal
#' @noRd
design_dcaps_marker <- function(genome, snp, flank, enzymes, min_gap = 50L) {
  limits <- caps_limits()
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  marker <- caps_marker(
    caps_dcaps_candidates(genome, snp, flank, enzymes, min_gap),
    snp, tmpl, enzymes, limits, KASP_GATE_FLOOR_CAPS
  )
  if (is.null(marker)) {
    cli::cli_inform(gettext("No natural or primer-introduced (dCAPS) restriction marker here."))
  }
  marker
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
