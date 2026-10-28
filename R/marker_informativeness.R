# =============================================================================
#  marker_informativeness.R  --  how much a marker tells a panel apart
# =============================================================================
#  Roadmap Y4.1.  A marker that separates two parents can be worthless on the
#  population those parents came from: the cross asks "is P1 different from P2",
#  the population asks "how many lines can this call split".  This module puts a
#  number on the second question from an available panel - either line
#  genotypes (`line`, `pos`, `genotype`) or the package's haplotype table - and
#  refuses to answer when the sample is too small rather than returning a
#  confident-looking number built on two individuals.
#
#  Reported quantities (all biallelic, since a target is one ref/alt pair):
#    - pic                     Botstein's polymorphism information content,
#                              1 - p^2 - q^2 - 2 p^2 q^2 (0.375 at p = q = 0.5)
#    - expected_heterozygosity 2 p q (Nei's gene diversity)
#    - distinguishable_pairs   fraction of line pairs the marker separates
#                              (genotype panels); with haplotypes it is the
#                              frequency-weighted probability that two draws
#                              differ in the target allele
#  `informativeness` for the Y4.2 quality axis is `pic`.
# =============================================================================

#' How informative a marker is on an available panel
#'
#' A hard gate asks whether a marker can be designed; this asks whether it is
#' worth using.  A marker that cleanly separates two parents may still be
#' monomorphic across the population they belong to, so the answer depends on
#' the panel, never on the cross alone.  The function accepts either a tidy
#' genotype table (`line`, `pos`, `genotype`, optionally `chrom`/`ref`/`alt`)
#' or the package's `haplotypes` table, and returns one row per target.
#'
#' Genotypes are alt dosages: `0`/`"REF"` (ref homozygote), `1`/`"HET"`, or
#' `2`/`"ALT"`; `NA` calls are dropped.  With genotypes the sample has to be
#' large enough to say anything - below `min_lines` the estimates are `NA` and
#' `reason` says why, instead of reporting the dichotomy of two individuals.
#' A haplotype table carries population frequencies itself, so it is measured
#' as a frequency model and no sample-size gate applies.
#'
#' All three statistics are biallelic: a target is one `ref`/`alt` pair.
#' - `pic` is Botstein's polymorphism information content,
#'   `1 - p^2 - q^2 - 2 p^2 q^2` (0.375 at `p = q = 0.5`).
#' - `expected_heterozygosity` is Nei's gene diversity `2 p q`.
#' - `distinguishable_pairs` is the fraction of line pairs the marker separates
#'   by genotype; it is `2 p q` under the haplotype model.
#'
#' `status` is `"ok"`, `"monomorphic"` (the panel carries one allele - a real
#' zero, not missing data), `"insufficient_sample"` or `"no_data"`.
#'
#' @param targets a target table (`chrom`, `pos`, `ref`, `alt`, optional
#'   `label`), a [marker_plan()], or a [marker_solve()] result whose plan is
#'   read.
#' @param genotypes a tidy table with `line`, `pos` and `genotype` (plus
#'   optional `chrom`, `ref`, `alt` to disambiguate records at the same
#'   position).  One row per line and variant.
#' @param haplotypes a haplotype table (`haplotype`, `freq`, `chrom`, `pos`,
#'   `ref`, `alt`, see [haplotype_groups()]).
#' @param min_lines minimum number of genotyped lines before a genotype panel
#'   is scored (default 4); ignored for the haplotype model.
#' @return a tibble, one row per target: `target`, `chrom`, `pos`, `ref`, `alt`,
#'   `n_lines`, `alt_freq`, `expected_heterozygosity`, `pic`,
#'   `distinguishable_pairs`, `status`, `reason`.
#' @examples
#' panel <- tibble::tibble(
#'   line = paste0("L", 1:8),
#'   pos = 250L,
#'   ref = "A",
#'   alt = "G",
#'   genotype = c(2, 2, 2, 2, 0, 0, 0, 0)
#' )
#' marker_informativeness(
#'   tibble::tibble(chrom = "Chr1", pos = 250L, ref = "A", alt = "G"),
#'   genotypes = panel
#' )
#' @family KASP informativeness
#' @export
marker_informativeness <- function(targets, genotypes = NULL, haplotypes = NULL,
                                   min_lines = 4L) {
  if (!is.null(genotypes) && !is.null(haplotypes)) {
    cli::cli_abort(gettext("Give exactly one of `genotypes` and `haplotypes`, not both."))
  }
  if (is.null(genotypes) && is.null(haplotypes)) {
    cli::cli_abort(gettext("Give `genotypes` (line calls) or `haplotypes` (a frequency model)."))
  }
  min_lines <- as.integer(min_lines)[1L]
  if (is.na(min_lines) || min_lines < 1L) {
    cli::cli_abort(gettext("`min_lines` must be a positive whole number."))
  }

  tg <- informativeness_targets(targets)
  if (!nrow(tg)) {
    return(informativeness_result(
      tg,
      n_lines = integer(), alt_freq = numeric(),
      het = numeric(), pic = numeric(), pairs = numeric(),
      status = character(), reason = character()
    ))
  }
  rows <- if (!is.null(genotypes)) {
    informativeness_genotypes(tg, normalise_genotypes(genotypes), min_lines)
  } else {
    informativeness_haplotypes(tg, normalise_haplotypes(haplotypes))
  }
  purrr::list_rbind(rows)
}

#' The targets a marker_informativeness() call scores (internal)
#' @keywords internal
#' @noRd
informativeness_targets <- function(x) {
  if (inherits(x, "marker_plan")) {
    x <- x$targets
  } else if (inherits(x, "marker_solved")) {
    x <- x$plan$targets
  }
  normalise_targets(x)[, c("label", "chrom", "pos", "ref", "alt")]
}

#' Normalise a tidy genotype table (internal)
#' @keywords internal
#' @noRd
normalise_genotypes <- function(g) {
  g <- tibble::as_tibble(g)
  need <- c("line", "pos", "genotype")
  missing <- setdiff(need, names(g))
  if (length(missing)) {
    cli::cli_abort(glue::glue(
      gettext("`genotypes` is missing column(s): {paste(missing, collapse = ', ')}.")
    ))
  }
  col <- function(nm, f) if (nm %in% names(g)) f(g[[nm]]) else rep(NA_character_, nrow(g))
  out <- tibble::tibble(
    line = as.character(g$line),
    chrom = col("chrom", as.character),
    pos = as.integer(g$pos),
    ref = col("ref", function(x) toupper(as.character(x))),
    alt = col("alt", function(x) toupper(as.character(x))),
    genotype = genotype_dosage(g$genotype)
  )
  if (anyNA(out$line) || anyNA(out$pos)) {
    cli::cli_abort(gettext("`genotypes` needs a non-missing `line` and `pos` on every row."))
  }
  out
}

#' Alt dosage from a genotype call (internal)
#'
#' Accepts dosages `0`/`1`/`2` or the calls `REF`/`HET`/`ALT` (and VCF-style
#' `0/0`, `0/1`, `1/1`, with either separator).  `NA` stays `NA`.
#' @keywords internal
#' @noRd
genotype_dosage <- function(x) {
  if (is.numeric(x)) {
    bad <- !is.na(x) & !(x %in% c(0, 1, 2))
    if (any(bad)) {
      cli::cli_abort(gettext("`genotype` dosages must be 0, 1 or 2."))
    }
    return(as.numeric(x))
  }
  calls <- toupper(as.character(x))
  map <- c(
    REF = 0, HET = 1, ALT = 2,
    "0/0" = 0, "0/1" = 1, "1/1" = 2,
    "0|0" = 0, "0|1" = 1, "1|1" = 2
  )
  out <- unname(map[calls])
  bad <- !is.na(calls) & is.na(out)
  if (any(bad)) {
    cli::cli_abort(gettext(
      "`genotype` calls must be REF/HET/ALT (or 0/0, 0/1, 1/1)."
    ))
  }
  out
}

#' The rows of a genotype panel that belong to one target (internal)
#' @keywords internal
#' @noRd
genotypes_for_target <- function(panel, t) {
  keep <- panel$pos == t$pos
  if (!all(is.na(panel$chrom))) {
    keep <- keep & (is.na(panel$chrom) | panel$chrom == t$chrom)
  }
  if (!all(is.na(panel$ref))) {
    keep <- keep & (is.na(panel$ref) | panel$ref == t$ref)
  }
  if (!all(is.na(panel$alt))) {
    keep <- keep & (is.na(panel$alt) | panel$alt == t$alt)
  }
  d <- panel[keep, , drop = FALSE]
  d[!duplicated(d$line), , drop = FALSE]
}

#' Score every target against a genotype panel (internal)
#' @keywords internal
#' @noRd
informativeness_genotypes <- function(tg, panel, min_lines) {
  purrr::map(seq_len(nrow(tg)), function(i) {
    t <- tg[i, ]
    d <- genotypes_for_target(panel, t)
    dos <- d$genotype[!is.na(d$genotype)]
    n <- length(dos)
    if (n == 0L) {
      return(informativeness_result(
        t,
        n_lines = 0L, alt_freq = NA_real_, het = NA_real_, pic = NA_real_,
        pairs = NA_real_, status = "no_data",
        reason = "No line has a genotype at this target."
      ))
    }
    if (n < min_lines) {
      return(informativeness_result(
        t,
        n_lines = n, alt_freq = NA_real_, het = NA_real_, pic = NA_real_,
        pairs = NA_real_, status = "insufficient_sample",
        reason = glue::glue(
          "At least {min_lines} genotyped lines are needed; the panel has {n}."
        )
      ))
    }
    p <- sum(dos) / (2 * n)
    informativeness_result(
      t,
      n_lines = n, alt_freq = p,
      het = 2 * p * (1 - p),
      pic = informativeness_pic(p),
      pairs = genotype_pair_fraction(dos),
      status = if (p == 0 || p == 1) "monomorphic" else "ok",
      reason = if (p == 0 || p == 1) "The panel carries only one allele at this target." else NA_character_
    )
  })
}

#' Score every target against the haplotype frequency model (internal)
#' @keywords internal
#' @noRd
informativeness_haplotypes <- function(tg, h) {
  freqs <- h |>
    dplyr::distinct(.data$haplotype, .data$freq) |>
    tibble::deframe()
  total <- sum(freqs)
  purrr::map(seq_len(nrow(tg)), function(i) {
    t <- tg[i, ]
    carriers <- unique(h$haplotype[h$pos == t$pos & h$alt == t$alt])
    p <- if (total > 0) sum(freqs[carriers]) / total else NA_real_
    status <- if (is.na(p)) {
      "no_data"
    } else if (p == 0 || p == 1) {
      "monomorphic"
    } else {
      "ok"
    }
    informativeness_result(
      t,
      n_lines = length(freqs), alt_freq = p,
      het = 2 * p * (1 - p),
      pic = informativeness_pic(p),
      # two independent draws differ in the target allele with probability 2pq
      pairs = 2 * p * (1 - p),
      status = status,
      reason = switch(status,
        no_data = "The haplotype table carries no frequency mass.",
        monomorphic = "The haplotypes carry only one allele at this target.",
        NA_character_
      )
    )
  })
}

#' Botstein's PIC for a biallelic marker (internal)
#' @keywords internal
#' @noRd
informativeness_pic <- function(p) {
  q <- 1 - p
  1 - (p^2 + q^2) - 2 * p^2 * q^2
}

#' Fraction of line pairs with different genotype calls (internal)
#' @keywords internal
#' @noRd
genotype_pair_fraction <- function(dos) {
  n <- length(dos)
  pairs <- choose(n, 2)
  if (pairs == 0) {
    return(NA_real_)
  }
  same <- sum(choose(as.numeric(table(dos)), 2))
  (pairs - same) / pairs
}

#' One result row per target (internal)
#' @keywords internal
#' @noRd
informativeness_result <- function(t, n_lines, alt_freq, het, pic, pairs,
                                   status, reason) {
  tibble::tibble(
    target = t$label,
    chrom = t$chrom,
    pos = t$pos,
    ref = t$ref,
    alt = t$alt,
    n_lines = n_lines,
    alt_freq = alt_freq,
    expected_heterozygosity = het,
    pic = pic,
    distinguishable_pairs = pairs,
    status = status,
    reason = reason
  )
}
