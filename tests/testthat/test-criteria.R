# The criterion registry is the single vocabulary for selection + acceptance,
# and check membership is family-based (KASP_FAMILY_CHECKS: marker-agnostic,
# KASP-selector *site*, and allele-specific).  These tests iterate over
# marker_schemes() and read each scheme's spec, so a new chemistry is covered with
# no edit here.

crit_mock <- function(scheme = "snp") {
  p <- kasp_primers()
  # the amplicon acceptance window is per-family, so a "clean" mock must sit inside
  # the scheme's own window (a CAPS band is not 90 bp)
  amp <- as.integer(mean(marker_amplicon_window(marker_scheme(scheme)$family)))
  list(
    assay = tibble::tibble(
      amplicon = amp, as_ref = p[[1]], as_alt = p[[2]], com = p[[3]], com_tm = 60,
      as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L, as_n_var = 0L,
      as_dropout = 0, as_foot_dropout = 0, com_dropout = 0
    ),
    thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), hairpin_dg = c(-1, -1, -1)),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3), b = rep(names(p), times = 3),
      dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
    ),
    amplicons = tibble::tibble(primer = c("AS_ref", "AS_alt"), n_products = c(1L, 0L)),
    efficiency = tibble::tibble(discrimination = 0.9)
  )
}

# the check a scheme's family should carry (used by the membership tests below)
family_checks <- function(scheme) {
  fam <- marker_scheme(scheme)$family
  classes <- KASP_FAMILY_CHECKS[[fam]] %||% character(0)
  c(
    KASP_CHECK_PRIMER,
    if ("site" %in% classes) KASP_CHECK_SITE,
    if ("as" %in% classes) KASP_CHECK_AS,
    # the family constants hold criterion ids; the registry labels are their names
    unname(names(MARKER_FAMILY_CHECKS[[fam]])) %||% character(0)
  )
}

test_that("marker_criteria declares gates, tiers and checks with metadata", {
  crit <- marker_criteria()
  expect_setequal(unique(crit$role), c("gate", "tier", "check"))
  expect_true(all(c("stage", "id", "label", "level", "hard", "limit", "order", "schemes") %in% names(crit)))
  # ZW4: the spec carries the evidence rung (the public view does not expose it yet),
  # right after the deprecated `calibrated` column so no existing column moves
  spec <- tidymarkers:::criterion_spec(KASP_LIMITS)
  expect_true(all(c("source", "calibrated", "evidence") %in% names(spec)))
  expect_equal(
    match("evidence", names(spec)), match("calibrated", names(spec)) + 1L
  )
  expect_true(all(crit$stage %in% c(
    "site", "as", "pair", "duplex", "asq", "caps", "arms", "hrm", "assay",
    marker_schemes()$scheme
  )))
  expect_true(any(crit$hard))
  # every criterion names at least one scheme, and all of them are registered
  expect_true(all(lengths(crit$schemes) >= 1L))
  expect_true(all(unlist(crit$schemes) %in% marker_schemes()$scheme))
  # the gates are the constants, and the hard ones are exactly those off the ladder
  pgate <- crit$id[crit$role == "gate" & crit$stage == "pair"]
  expect_setequal(pgate, paste0("pair:gate:", KASP_GATES_PAIR))
  expect_setequal(
    crit$id[crit$role == "gate" & crit$hard & crit$stage %in% c("as", "pair")],
    c(
      sprintf("as:gate:%s", setdiff(KASP_GATES_AS, KASP_RELAX_AS)),
      sprintf("pair:gate:%s", setdiff(KASP_GATES_PAIR, KASP_RELAX_PAIR))
    )
  )
})

test_that("the check classes partition the registered checks", {
  labels <- marker_criteria()$label[marker_criteria()$role == "check"]
  family_labels <- unique(unlist(lapply(MARKER_FAMILY_CHECKS, names), use.names = FALSE))
  expect_setequal(
    c(KASP_CHECK_PRIMER, KASP_CHECK_SITE, KASP_CHECK_AS, family_labels), labels
  )
  # 8 marker-agnostic + 5 site + 14 allele-specific + 29 family-specific = 56
  # (the 14th is `as_mismatch_type`, Z0.4: the terminal mismatch's Kwok class)
  expect_length(KASP_CHECK_PRIMER, 8L)
  expect_length(KASP_CHECK_SITE, 5L)
  expect_length(KASP_CHECK_AS, 14L)
  expect_length(family_labels, 29L)
  expect_length(labels, 56L)
  # distinct classes: a base check is never also family-specific
  expect_length(intersect(c(KASP_CHECK_PRIMER, KASP_CHECK_SITE, KASP_CHECK_AS), family_labels), 0L)
  expect_length(intersect(KASP_CHECK_PRIMER, KASP_CHECK_SITE), 0L)
  expect_length(intersect(KASP_CHECK_PRIMER, KASP_CHECK_AS), 0L)
  expect_length(intersect(KASP_CHECK_SITE, KASP_CHECK_AS), 0L)
})

test_that("every scheme carries exactly its family's check classes", {
  for (sc in marker_schemes()$scheme) {
    chk <- marker_criteria(scheme = sc)
    chk <- chk$label[chk$role == "check"]
    # C2: a scheme carries only the checks its chemistry obeys
    expect_true(all(chk %in% family_checks(sc)), info = sc)
    # ... and every marker-agnostic check is there (it holds for any primer pair)
    expect_true(all(KASP_CHECK_PRIMER %in% chk), info = sc)
    # the AS-only and site checks track the family, never the scheme
    info <- marker_scheme(sc)
    as_family <- "as" %in% (KASP_FAMILY_CHECKS[[info$family]] %||% character(0))
    expect_equal(any(KASP_CHECK_AS %in% chk), as_family, info = sc)
  }
  # the multiplex rules are duplex-only
  expect_true(all(paste0("duplex:gate:", KASP_GATES_DUPLEX) %in% marker_criteria(scheme = "duplex")$id))
  expect_false(any(grepl("duplex", marker_criteria(scheme = "snp")$id)))
})

test_that("the AS/COM oligo checks reach only families with an AS primer (C2)", {
  # an AS/COM Tm or GC window is meaningless where those oligos do not exist, so
  # the moved labels follow the `as` class and never leak to a flanking pair.
  as_com <- c("AS Tm (C)", "COM Tm (C)", "AS GC", "COM GC")
  for (sc in marker_schemes()$scheme) {
    crit <- marker_criteria(scheme = sc)
    chk <- crit$label[crit$role == "check"]
    classes <- KASP_FAMILY_CHECKS[[marker_scheme(sc)$family]] %||% character(0)
    expect_equal(all(as_com %in% chk), "as" %in% classes, info = sc)
  }
})

test_that("marker_design_checks returns exactly the registry's checks, in order", {
  for (sc in marker_schemes()$scheme) {
    crit <- marker_criteria(scheme = sc)
    checks <- crit[crit$role == "check", ]
    ch <- marker_design_checks(crit_mock(sc), scheme = sc)
    expect_equal(ch$criterion, checks$label, info = sc)
    expect_equal(ch$limit, checks$limit, info = sc)
    # a KASP mock cannot supply CAPS/LAMP/RPA/... fields; those checks are
    # unknown/not_applicable and must never be reported as a failure.  Every
    # check that is actually evaluable must pass, and `pass` is exactly the
    # logical view of `status`.
    expect_false(any(ch$status == "fail"), info = sc)
    expect_true(all(ch$pass[ch$status == "pass"]), info = sc)
    expect_equal(ch$pass, ch$status == "pass", info = sc)
  }
})

test_that("a criterion's threshold and its rendered limit cannot diverge (Z2)", {
  spec <- tidymarkers:::criterion_spec()
  crit <- marker_criteria()
  # the public view is a projection of the spec, row for row
  expect_equal(nrow(crit), nrow(spec))
  expect_equal(crit$id, spec$id)
  expect_equal(crit$limit, spec$limit)

  numeric_bar <- vapply(
    spec$bounds, function(b) is.numeric(b) && length(b) >= 2L, logical(1)
  )
  # a KASP bar is rendered from its own bound, by the one formatter
  kasp <- numeric_bar & spec$scope == "kasp"
  expect_gt(sum(kasp), 0L)
  expect_equal(
    crit$limit[kasp],
    vapply(spec$bounds[kasp], tidymarkers:::criterion_bound_text, character(1))
  )
  # a family bar keeps its published prose, but its numbers must be in that prose
  fam <- numeric_bar & spec$scope == "family"
  expect_gt(sum(fam), 0L)
  for (i in which(fam)) {
    for (v in spec$bounds[[i]][is.finite(spec$bounds[[i]])]) {
      expect_true(
        grepl(as.character(v), crit$limit[[i]], fixed = TRUE),
        info = paste(spec$id[[i]], "bound", v, "absent from", crit$limit[[i]])
      )
    }
  }
})

test_that("the scheme filter moves the amplicon bound and its limit together (Z2)", {
  for (sc in c("caps", "lamp", "ssr", "taqman")) {
    win <- tidymarkers:::marker_amplicon_window(marker_scheme(sc)$family)
    lim <- KASP_LIMITS
    lim$amplicon <- win
    spec <- tidymarkers:::criterion_spec(lim)
    i <- match("amplicon", spec$id)
    # one object: the bound the selector grades against and the text the results
    # show are the same numbers
    expect_equal(spec$bounds[[i]], win, info = sc)
    crit <- marker_criteria(scheme = sc)
    expect_equal(
      crit$limit[match("amplicon length (bp)", crit$label)],
      spec$limit[[i]],
      info = sc
    )
  }
})
