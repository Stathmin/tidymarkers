# =============================================================================
#  marker_pareto.R  --  the front over cost, quality and linkage risk
# =============================================================================
#  Roadmap Y4.2.  A single "best" marker hides the trade-off: the cheapest
#  chemistry is rarely the most informative one, and a site at the edge of a
#  region is a worse bet than one near its middle.  With three axes the answer
#  is a Pareto front plus one recommended point, not a lexicographic minimum.
#
#  Axes (all numeric, oriented so that smaller is better except quality):
#    cost            relative chemistry cost (marker_costs())
#    quality         [0, 1]: assay quality times informativeness
#    linkage_risk    [0, 1]: Haldane recombination fraction between the marker
#                    and an anchor, from a declared physical rate
#
#  `linkage_risk` is a rule of thumb, not a calibrated probability: the rate is
#  a declared Morgans-per-bp figure, and no recombination map is read (roadmap
#  Y4.2's blocker: Z8 calibration is not available, so the axis says so).
# =============================================================================

#' A Pareto front over marker candidates
#'
#' Ranks candidates on three axes instead of one: **cost** (relative chemistry
#' cost, see [marker_costs()]), **quality** (assay quality times informativeness
#' on the available panel, see [marker_informativeness()]) and **linkage risk**
#' (the chance the marker has recombined away from the anchor).  A candidate is
#' *dominated* when another is at least as cheap, at least as good and at least
#' as low-risk, and strictly better on one axis; dominated candidates are
#' excluded from the front.  Exactly one front point is marked `recommended`:
#' the one closest to the ideal point after min-max normalisation of the three
#' axes.
#'
#' `x` is either a candidate table or a [marker_solve()] result.  A candidate
#' table needs a non-negative `cost` and a `quality` in `[0, 1]`; an optional
#' `linkage_risk` in `[0, 1]` defaults to 0.  For a solved object the cost comes
#' from its `report`; `quality` defaults to a graded [marker_risk()] verdict
#' (`low` / `medium` / `high` mapped to 1 / 0.5 / 0) and is multiplied by
#' `pic / 0.375` when `genotypes` or `haplotypes` is supplied.
#'
#' **Linkage risk is a rule of thumb.** With an `anchor` (or, for a solved
#' region, the region midpoint) the distance in bp is converted to Morgans with
#' `recombination` (Morgans per bp, default `1e-8` - 1 cM/Mb) and read through
#' Haldane's map function, `(1 - exp(-2d)) / 2`.  No recombination map is
#' consulted and the result is not a calibrated probability; it is a monotone
#' ordering of "nearer is safer" that the front can use.  With no anchor and no
#' region every point sits at the marker itself, so the risk is 0.
#'
#' @param x a candidate table (`cost`, `quality`, optional `linkage_risk` and
#'   any identifying columns) or a [marker_solve()] result.
#' @param quality optional quality override; a number per candidate (a scalar is
#'   recycled for a solved object).
#' @param linkage_risk optional linkage-risk override, one number per candidate.
#' @param anchor numeric position the marker should stay linked to, or `NULL`.
#'   Ignored when `linkage_risk` is supplied.
#' @param recombination Morgans per bp for the Haldane conversion (default
#'   `1e-8`).
#' @param genotypes,haplotypes optional panel passed to
#'   [marker_informativeness()] to weight the quality axis.
#' @param min_lines passed to [marker_informativeness()].
#' @return `x` as a tibble with `on_front` and `recommended` added; for a solved
#'   object the columns are `label`, `chrom`, `pos`, `cost`, `quality`,
#'   `linkage_risk`, `on_front`, `recommended`.
#' @examples
#' candidates <- tibble::tibble(
#'   label = c("cheap", "careful"),
#'   cost = c(1, 4),
#'   quality = c(0.5, 0.9),
#'   linkage_risk = c(0.05, 0.01)
#' )
#' marker_pareto(candidates)
#' @family KASP ranking
#' @export
marker_pareto <- function(x, quality = NULL, linkage_risk = NULL, anchor = NULL,
                          recombination = 1e-8, genotypes = NULL,
                          haplotypes = NULL, min_lines = 4L) {
  cand <- pareto_candidates(
    x, quality, linkage_risk, anchor, recombination,
    genotypes, haplotypes, min_lines
  )
  pareto_front(cand)
}

#' Build the candidate table a Pareto front is computed over (internal)
#' @keywords internal
#' @noRd
pareto_candidates <- function(x, quality, linkage_risk, anchor, recombination,
                              genotypes, haplotypes, min_lines) {
  if (inherits(x, "marker_solved")) {
    base <- tibble::tibble(
      label = x$report$label,
      chrom = x$report$chrom,
      pos = x$report$pos,
      cost = x$report$cost
    )
    if (is.null(quality)) {
      quality <- vapply(x$report$marker, pareto_risk_quality, numeric(1))
    }
    if (!is.null(genotypes) || !is.null(haplotypes)) {
      quality <- quality * pareto_informativeness(
        x, genotypes, haplotypes, min_lines
      )
    }
    if (is.null(linkage_risk)) {
      linkage_risk <- pareto_linkage_risk(base, x$plan, anchor, recombination)
    }
    base$quality <- quality
    base$linkage_risk <- linkage_risk
  } else if (is.data.frame(x)) {
    base <- tibble::as_tibble(x)
    if (!"cost" %in% names(base)) {
      cli::cli_abort(gettext("A candidate table needs a `cost` column."))
    }
    if (!is.null(quality)) {
      pareto_override_length(quality, nrow(base), "quality")
      base$quality <- quality
    }
    if (!"quality" %in% names(base)) {
      cli::cli_abort(gettext(
        "A candidate table needs a `quality` column (or an explicit `quality`)."
      ))
    }
    if (!is.null(linkage_risk)) {
      pareto_override_length(linkage_risk, nrow(base), "linkage_risk")
      base$linkage_risk <- linkage_risk
    }
    if (!"linkage_risk" %in% names(base)) {
      base$linkage_risk <- 0
    }
  } else {
    cli::cli_abort(gettext(
      "`x` must be a candidate table or a marker_solve() result."
    ))
  }
  pareto_validate(base)
  base
}

#' Check an explicit axis override has one value per candidate (internal)
#' @keywords internal
#' @noRd
pareto_override_length <- function(v, n, nm) {
  if (!is.numeric(v) || length(v) != n) {
    cli::cli_abort(gettextf("`%s` must be one number per candidate.", nm))
  }
  invisible(v)
}

#' A graded quality from a solved marker's risk verdict (internal)
#'
#' `low` / `medium` / `high` map to 1 / 0.5 / 0.  This is an ordering, not a
#' probability; it exists so a front can use the existing verdict.
#' @keywords internal
#' @noRd
pareto_risk_quality <- function(marker) {
  if (is.null(marker)) {
    return(NA_real_)
  }
  r <- tryCatch(marker_risk(marker)$risk, error = function(e) NA_character_)
  r <- as.character(r)[1L]
  if (is.na(r)) {
    return(NA_real_)
  }
  switch(r,
    low = 1,
    medium = 0.5,
    high = 0,
    NA_real_
  )
}

#' The informativeness factor for the solved candidates (internal)
#'
#' Information is a multiplier on quality: a marker that cannot tell the
#' available lines apart is worth less for this purpose.  `pic` is normalised by
#' its biallelic maximum (0.375); an estimate that is `NA` (small or missing
#' panel) is neutral, not zero.
#' @keywords internal
#' @noRd
pareto_informativeness <- function(solved, genotypes, haplotypes, min_lines) {
  inf <- marker_informativeness(
    solved,
    genotypes = genotypes, haplotypes = haplotypes, min_lines = min_lines
  )
  ref <- solved$report$label
  pic <- inf$pic[match(ref, inf$target)]
  factor <- pic / informativeness_pic(0.5)
  factor[is.na(factor)] <- 1
  factor
}

#' Haldane linkage risk between candidates and an anchor (internal)
#' @keywords internal
#' @noRd
pareto_linkage_risk <- function(base, plan, anchor, recombination) {
  if (is.null(anchor)) {
    reg <- plan$region %||% NULL
    if (is.null(reg)) {
      return(rep(0, nrow(base)))
    }
    anchor <- (as.numeric(reg$start) + as.numeric(reg$end)) / 2
  }
  d <- abs(as.numeric(base$pos) - as.numeric(anchor))
  haldane_fraction(d, recombination)
}

#' Haldane's recombination fraction from a physical distance (internal)
#' @keywords internal
#' @noRd
haldane_fraction <- function(bp, per_bp = 1e-8) {
  d <- as.numeric(bp) * per_bp
  (1 - exp(-2 * d)) / 2
}

#' Reject a malformed candidate table (internal)
#' @keywords internal
#' @noRd
pareto_validate <- function(d) {
  limits <- list(cost = c(0, Inf), quality = c(0, 1), linkage_risk = c(0, 1))
  for (nm in names(limits)) {
    v <- d[[nm]]
    if (!is.numeric(v) || anyNA(v) || any(!is.finite(v))) {
      cli::cli_abort(gettextf("`%s` must be a finite number on every candidate.", nm))
    }
    lo <- limits[[nm]][1L]
    hi <- limits[[nm]][2L]
    if (any(v < lo) || any(v > hi)) {
      cli::cli_abort(gettextf("`%s` must lie in [%s, %s].", nm, lo, hi))
    }
  }
  d
}

#' Mark the front and the recommended point (internal)
#' @keywords internal
#' @noRd
pareto_front <- function(d) {
  n <- nrow(d)
  dominated <- logical(n)
  if (n > 1L) {
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j && pareto_dominates(d[j, ], d[i, ])) {
          dominated[i] <- TRUE
          break
        }
      }
    }
  }
  d$on_front <- !dominated
  d$recommended <- FALSE
  if (any(d$on_front)) {
    score <- pareto_score(d)
    idx <- which(d$on_front)
    d$recommended[idx[which.min(score[idx])]] <- TRUE
  }
  d
}

#' Does `a` dominate `b`? (internal)
#' @keywords internal
#' @noRd
pareto_dominates <- function(a, b) {
  no_worse <- a$cost <= b$cost &&
    a$quality >= b$quality &&
    a$linkage_risk <= b$linkage_risk
  better <- a$cost < b$cost ||
    a$quality > b$quality ||
    a$linkage_risk < b$linkage_risk
  no_worse && better
}

#' Distance-to-ideal score over the three axes (internal)
#'
#' Each axis is min-max normalised over the candidate set, with cost and risk
#' oriented small-is-good and quality large-is-good.  Ties are broken by the
#' caller's row order.
#' @keywords internal
#' @noRd
pareto_score <- function(d) {
  norm <- function(v, invert = FALSE) {
    rng <- range(v)
    if (rng[1L] == rng[2L]) {
      return(rep(0, length(v)))
    }
    out <- (v - rng[1L]) / (rng[2L] - rng[1L])
    if (invert) 1 - out else out
  }
  norm(d$cost) + norm(d$linkage_risk) + norm(d$quality, invert = TRUE)
}
