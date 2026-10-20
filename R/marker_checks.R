# =============================================================================
#  marker_checks.R  --  acceptance checks (driven by the criterion registry)
# =============================================================================

#' Longest homopolymer run in a sequence
#' @keywords internal
#' @noRd
longest_run <- function(seq) {
  max(nchar(stringr::str_extract_all(seq, "([ACGT])\\1*")[[1L]]))
}

# ---------------------------------------------------------------------------
#  check statuses
# ---------------------------------------------------------------------------
# One vocabulary for what an acceptance check says about a marker:
#   pass           - evaluated, inside the bar
#   fail           - evaluated, outside the bar
#   unknown        - could not be evaluated (missing data, engine unavailable);
#                    never a pass
#   not_applicable - the check does not apply to this marker class (e.g. a
#                    CAPS-only check on a KASP assay); neither a pass nor a fail
MARKER_CHECK_STATUS_LEVELS <- c("pass", "fail", "unknown", "not_applicable")

#' Construct a check-status vector
#' @param status character; `NA` (or anything unrecognised) becomes `"unknown"`.
#' @return a `forcats` factor with levels `pass`, `fail`, `unknown`,
#'   `not_applicable`.
#' @keywords internal
#' @noRd
marker_check_status <- function(status = character()) {
  status <- as.character(status)
  status[is.na(status) | !status %in% MARKER_CHECK_STATUS_LEVELS] <- "unknown"
  forcats::fct(status, levels = MARKER_CHECK_STATUS_LEVELS)
}

#' The legacy logical `pass` view of a status
#'
#' Only an evaluated `pass` is `TRUE`; `unknown` and `not_applicable` are
#' deliberately not successes.
#' @keywords internal
#' @noRd
marker_check_pass <- function(status) {
  as.character(marker_check_status(status)) == "pass"
}

#' Check an assembled assay against the KASP acceptance criteria
#'
#' Turns the prose design manual into a tidy pass/fail table so "good" and
#' "bad" assays are graded the same way every time (thresholds: `KASP_LIMITS`).
#' The criterion list, levels and limits come from [marker_criteria()] (its
#' `role == "check"` rows), so the checks and the decision tree share one
#' vocabulary.  Reads only the assay object produced by [assemble_kasp_assay()],
#' so it is pure and cheap to test.
#'
#' Only the checks the scheme's DAG displays are evaluated: `scheme` (or the
#' assay's own `$scheme`) filters [marker_criteria()] to that chemistry, so a
#' non-allele-specific scheme is never graded on AS/COM/REF-ALT facts it does not
#' have (C3).  The values are computed lazily, so a scheme that lacks a check
#' never reads the field that check needs.  The `amplicon length (bp)` check is
#' graded against the *scheme's own* amplicon window (`marker_amplicon_window()`),
#' so a CAPS band, Sanger read, SSR / InDel-size fragment or isothermal product is
#' never failed for the 50-100 bp KASP optimum (C2/C5).
#'
#' A criterion the selection ladder lowered for this winner (see
#' `assay$selection$relaxed`/`$violated`/`$as`) is graded against the *bar the
#' selection applied*, not the ideal: the gate is soft precisely because the method
#' tolerates the deviation, so "relaxed" must not also read as "failed" (C5).  Its
#' `limit` then carries the relaxed bar and a `(relaxed)` marker.
#'
#' @param assay result of [assemble_kasp_assay()].
#' @param limits list of thresholds; defaults to `KASP_LIMITS`.
#' @param scheme the marker chemistry the assay is graded as (see
#'   [marker_schemes()]); `NULL` (the default) uses `assay$scheme`, falling back to
#'   `"snp"`.
#' @return tibble with `stage`, `criterion`, `observed`, `limit`, `pass`,
#'   `status`, `detail` and `level`.  `status` is a factor (`pass`, `fail`,
#'   `unknown`, `not_applicable`); `pass` is its logical view
#'   (`status == "pass"`), kept for existing callers.  Criteria cover the
#'   target *site*, each *feature* oligo, the *alignment*, and *off-target*
#'   specificity - as far as the scheme has them.
#' @examples
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' marker_design_checks(assay)
#' }
#' @family KASP checks
#' @export
marker_design_checks <- function(assay, limits = KASP_LIMITS, scheme = NULL) {
  scheme <- scheme %||% assay$scheme %||% if (inherits(assay, "marker_assay")) {
    "snp"
  } else {
    cls <- class(assay)[1L]
    if (startsWith(cls, "marker_")) sub("^marker_", "", cls) else "snp"
  }
  # Non-KASP markers use the same family criterion set and the generic
  # amplicon/primer/structure/specificity evaluator.  The KASP path below reads
  # the richer assay object only when that object is present.
  if (!inherits(assay, "marker_assay") && grepl("^marker_", class(assay)[1L])) {
    return(marker_generic_checks(assay, scheme, limits = limits))
  }
  crit <- marker_criteria(limits, scheme = scheme)
  crit <- crit[crit$role == "check", ]
  amp_win <- marker_amplicon_window(marker_scheme(scheme)$family, limits)

  tol <- limits$footprint_dropout %||% 0.05
  aa <- assay$assay[1L, ]
  th <- assay$thermo
  d <- assay$complementarity

  # Each ingredient is a one-shot memo: computed at most once, and only if a
  # check the scheme owns asks for it, so a chemistry without the AS primer never
  # reads th$tm / assay$amplicons / assay$ideogram.
  memo <- function(f) {
    val <- NULL
    done <- FALSE
    function() {
      if (!done) {
        val <<- f()
        done <<- TRUE
      }
      val
    }
  }
  site <- memo(function() assay$site)
  site_class <- memo(function() {
    s <- site()
    if (is.null(s)) "biallelic" else as.character(s$class[1L])
  })
  site_disc <- memo(function() {
    s <- site()
    if (is.null(s)) TRUE else isTRUE(s$disc[1L])
  })
  n_alt <- memo(function() {
    s <- site()
    if (is.null(s)) 1L else as.integer(s$n_alt[1L])
  })
  linked_skipped <- memo(function() as.integer(assay$selection$linked_skipped %||% 0L))
  # how many of the site's ALTs the ALT primer actually covers (the degenerate /
  # linked engine surfaces it); a plain biallelic site covers the one ALT
  alt_covered <- memo(function() {
    s <- site()
    if (is.null(s) || !"alt_covered" %in% names(s)) n_alt() else as.integer(s$alt_covered[1L])
  })
  n_uncovered <- memo(function() max(0L, n_alt() - alt_covered()))
  linked_dependent <- memo(function() isTRUE(assay$selection$linked_dependent))
  as_tm <- memo(function() mean(th$tm[th$name %in% c("AS_ref", "AS_alt")]))
  com_tm <- memo(function() aa$com_tm)
  hair <- memo(function() min(th$hairpin_dg))
  cross <- memo(function() min(d$dg[d$a != d$b])) # worst hetero-dimer
  homo <- memo(function() min(d$dg[d$a == d$b])) # worst self-dimer
  run <- memo(function() max(longest_run(aa$as_ref), longest_run(aa$com)))
  amp <- memo(function() assay$amplicons)
  # in-silico PCR (only if computed): REF primer -> 1 product, ALT primer -> 0
  amp_obs <- memo(function() {
    a <- amp()
    if (is.null(a)) NA_real_ else max(a$n_products)
  })
  amp_pass <- memo(function() {
    a <- amp()
    if (is.null(a)) {
      FALSE
    } else {
      ref <- a$n_products[a$primer == "AS_ref"]
      alt <- a$n_products[a$primer == "AS_alt"]
      length(ref) == 1L && length(alt) == 1L && ref == 1L && alt == 0L
    }
  })
  # off-target: exact genome matches of each primer (from the ideogram, computed
  # once).  This is an *exact* hit count - a 1-mismatch 3'-seed count is noise in
  # a 100+ Mb genome and would not mean a real off-target amplicon.
  ot_obs <- memo(function() {
    ih <- assay$ideogram$hits
    if (is.null(ih) || !nrow(ih)) NA_real_ else as.numeric(max(table(ih$primer)))
  })
  ot_pass <- memo(function() {
    o <- ot_obs()
    is.na(o) || o <= 1
  })

  # one entry per registered check, keyed by label; `detail` (only the site class
  # carries one) falls back to NA
  checks <- list(
    "target site class" = function() {
      list(observed = NA_real_, pass = site_disc(), detail = site_class())
    },
    "ALT alleles (n)" = function() {
      list(observed = as.numeric(n_alt()), pass = n_alt() >= 1L)
    },
    "linked variants skipped (n)" = function() {
      list(observed = as.numeric(linked_skipped()), pass = linked_skipped() == 0L)
    },
    "ALT primer covers all ALTs" = function() {
      list(observed = as.numeric(n_uncovered()), pass = n_uncovered() == 0L)
    },
    "design relies on linkage" = function() {
      list(observed = as.numeric(linked_dependent()), pass = !linked_dependent())
    },
    "amplicon length (bp)" = function() {
      list(observed = aa$amplicon, pass = dplyr::between(aa$amplicon, amp_win[1L], amp_win[2L]))
    },
    "AS Tm (C)" = function() {
      list(observed = as_tm(), pass = dplyr::between(as_tm(), limits$as_tm[1], limits$as_tm[2]))
    },
    "COM Tm (C)" = function() {
      list(observed = com_tm(), pass = dplyr::between(com_tm(), limits$com_tm[1], limits$com_tm[2]))
    },
    "COM Tm >= AS Tm" = function() {
      list(observed = com_tm() - as_tm(), pass = com_tm() >= as_tm())
    },
    "AS GC" = function() {
      list(observed = aa$as_gc, pass = dplyr::between(aa$as_gc, limits$gc[1], limits$gc[2]))
    },
    "COM GC" = function() {
      list(observed = aa$com_gc, pass = dplyr::between(aa$com_gc, limits$gc[1], limits$gc[2]))
    },
    "longest homopolymer (bp)" = function() {
      list(observed = run(), pass = run() < limits$max_run)
    },
    "SSR (2-4 nt, >=3 copies)" = function() {
      list(
        observed = sum(lengths(find_repeats(c(aa$as_ref, aa$com)))),
        pass = !has_repeat(aa$as_ref) && !has_repeat(aa$com)
      )
    },
    "inverted repeat (arms>=4, loop<=6)" = function() {
      list(
        observed = sum(purrr::map_int(find_palindromes(c(aa$as_ref, aa$as_alt, aa$com)), nrow)),
        pass = !any(has_palindrome(c(aa$as_ref, aa$as_alt, aa$com)))
      )
    },
    "worst hairpin dG37" = function() {
      list(observed = hair(), pass = hair() >= limits$hairpin_dg)
    },
    "worst cross-dimer dG37" = function() {
      list(observed = cross(), pass = cross() >= limits$duplex_dg)
    },
    "worst homodimer dG37" = function() {
      list(observed = homo(), pass = homo() >= limits$homodimer_dg)
    },
    "AS 3' window variants (n)" = function() {
      list(observed = aa$as_n_3p, pass = !is.na(aa$as_dropout) | aa$as_n_3p == 0L)
    },
    "AS footprint variants (n)" = function() {
      list(observed = aa$as_n_var, pass = !is.na(aa$as_foot_dropout) | aa$as_n_var == 0L)
    },
    "COM footprint variants (n)" = function() {
      list(observed = aa$n_var, pass = !is.na(aa$com_dropout) | aa$n_var == 0L)
    },
    "AS 3' window dropout" = function() {
      list(observed = aa$as_dropout, pass = is.na(aa$as_dropout) | aa$as_dropout <= tol)
    },
    "AS footprint dropout" = function() {
      list(observed = aa$as_foot_dropout, pass = is.na(aa$as_foot_dropout) | aa$as_foot_dropout <= tol)
    },
    "COM footprint dropout" = function() {
      list(observed = aa$com_dropout, pass = is.na(aa$com_dropout) | aa$com_dropout <= tol)
    },
    "discrimination" = function() {
      disc <- assay$efficiency$discrimination[1L]
      list(observed = disc, pass = disc >= limits$discrimination)
    },
    "in-silico PCR (REF 1 / ALT 0)" = function() {
      list(observed = amp_obs(), pass = amp_pass())
    },
    "off-target exact primer hits (n)" = function() {
      list(observed = ot_obs(), pass = ot_pass())
    }
  )

  rows <- lapply(crit$label, function(l) {
    if (l %in% names(checks)) checks[[l]]() else list(observed = NA_real_, pass = NA)
  })
  observed <- vapply(rows, function(r) as.numeric(r$observed), numeric(1))
  limit <- crit$limit
  pass <- vapply(rows, function(r) {
    p <- r$pass
    if (!length(p) || is.na(p[1L])) NA else isTRUE(p[1L])
  }, logical(1))
  # A criterion the ladder lowered for *this* winner is judged against the bar the
  # selection actually applied, not the ideal, so an explicitly accepted
  # (relaxed, warned) result is not also reported as an acceptance failure (C5).
  # The floors are the method's tolerances for a soft gate, so the relaxed bar is
  # still an accepted value; `limit` then shows the bar that was in force.
  relaxed <- unique(c(assay$selection$relaxed, assay$selection$violated, assay$selection$as))
  for (g in intersect(relaxed, names(KASP_GATE_CHECK))) {
    for (lab in KASP_GATE_CHECK[[g]]) {
      i <- match(lab, crit$label)
      if (is.na(i)) next
      bar <- marker_gate_relaxed_bar(lab, limits)
      limit[i] <- bar$limit
      pass[i] <- bar$test(observed[i])
    }
  }
  status <- ifelse(is.na(pass), "unknown", ifelse(pass, "pass", "fail"))
  # The KASP path is reached by a `marker_assay` (or an unclassed mock): family
  # checks that belong to another marker class do not apply here rather than
  # being silently unknown.
  native <- identical(class(assay)[1L], paste0("marker_", scheme))
  fam_labs <- MARKER_FAMILY_CHECKS[[marker_scheme(scheme)$family]] %||% character(0)
  status[crit$label %in% fam_labs & !native] <- "not_applicable"
  # Missing data is `unknown`, never a pass.  Several checks used to coerce an
  # `NA` to TRUE (a dropout with no value, an uncomputed off-target count, an
  # empty dimer table); this one place makes the rule uniform.  The one check
  # whose value is deliberately textual (`target site class`) is exempt - its
  # `observed` is NA by construction, not because data are missing.
  missing_obs <- !is.finite(observed) & !crit$label %in% "target site class"
  status[missing_obs & status != "not_applicable"] <- "unknown"
  status <- marker_check_status(status)
  tibble::tibble(
    stage = crit$stage,
    criterion = crit$label,
    observed = observed,
    limit = limit,
    pass = marker_check_pass(status),
    status = status,
    detail = vapply(rows, function(r) r$detail %||% NA_character_, character(1)),
    level = crit$level
  )
}

#' Is an evaluation usable?
#'
#' Collapses an acceptance-check table into one verdict: `"usable"` only when
#' every *applicable* check is `pass`.  A `fail` (or a legacy `pass = FALSE`)
#' is `"doomed"`; an `unknown` row keeps the verdict `"unknown"` instead of
#' silently passing.  `not_applicable` rows are ignored - a check that does not
#' apply can neither certify nor doom a marker.
#'
#' @param x a checks tibble (see [marker_design_checks()]), or a marker/assay
#'   whose `$evaluation` is read, falling back to [marker_design_checks()].
#' @param limits thresholds; see [marker_design_checks()].
#' @return a single character: `"usable"`, `"doomed"` or `"unknown"`.
#' @keywords internal
#' @noRd
evaluation_verdict <- function(x, limits = KASP_LIMITS) {
  ev <- if (is.data.frame(x)) {
    x
  } else {
    x$evaluation %||% tryCatch(
      marker_design_checks(x, limits = limits),
      error = function(e) NULL
    )
  }
  if (is.null(ev) || !nrow(ev)) {
    return("unknown")
  }
  status <- if ("status" %in% names(ev)) {
    as.character(marker_check_status(ev$status))
  } else if ("pass" %in% names(ev)) {
    ifelse(is.na(ev$pass), "unknown", ifelse(ev$pass, "pass", "fail"))
  } else {
    rep("unknown", nrow(ev))
  }
  if (any(status == "fail")) {
    return("doomed")
  }
  if (any(status == "unknown")) {
    return("unknown")
  }
  if (!any(status == "pass")) {
    return("unknown")
  }
  "usable"
}
