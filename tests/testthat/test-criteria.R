# The criterion registry is the single vocabulary for selection + acceptance,
# and check membership is family-based (KASP_FAMILY_CHECKS: marker-agnostic,
# KASP-selector *site*, and allele-specific).  These tests iterate over
# kasp_schemes() and read each scheme's spec, so a new chemistry is covered with
# no edit here.

crit_mock <- function(scheme = "snp") {
  p <- c(AS_ref = "CATGGACCTAGCATTGACGTA", AS_alt = "CATGGACCTAGCATTGACGTT", COM = "GCTTAGCCATGGTCAATGCTA")
  # the amplicon acceptance window is per-family, so a "clean" mock must sit inside
  # the scheme's own window (a CAPS band is not 90 bp)
  amp <- as.integer(mean(kasp_amplicon_window(kasp_scheme(scheme)$family)))
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
  classes <- KASP_FAMILY_CHECKS[[kasp_scheme(scheme)$family]] %||% character(0)
  c(
    KASP_CHECK_PRIMER,
    if ("site" %in% classes) KASP_CHECK_SITE,
    if ("as" %in% classes) KASP_CHECK_AS
  )
}

test_that("kasp_criteria declares gates, tiers and checks with metadata", {
  crit <- kasp_criteria()
  expect_setequal(unique(crit$role), c("gate", "tier", "check"))
  expect_true(all(c("stage", "id", "label", "level", "hard", "limit", "order", "schemes") %in% names(crit)))
  expect_true(all(crit$stage %in% c(
    "site", "as", "pair", "duplex", "asq", "caps", "arms", "hrm", "assay",
    kasp_schemes()$scheme
  )))
  expect_true(any(crit$hard))
  # every criterion names at least one scheme, and all of them are registered
  expect_true(all(lengths(crit$schemes) >= 1L))
  expect_true(all(unlist(crit$schemes) %in% kasp_schemes()$scheme))
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

test_that("the three check classes partition the registered checks", {
  labels <- kasp_criteria()$label[kasp_criteria()$role == "check"]
  expect_setequal(c(KASP_CHECK_PRIMER, KASP_CHECK_SITE, KASP_CHECK_AS), labels)
  # distinct classes: no criterion is both marker-agnostic and allele-specific
  expect_length(intersect(KASP_CHECK_PRIMER, KASP_CHECK_SITE), 0L)
  expect_length(intersect(KASP_CHECK_PRIMER, KASP_CHECK_AS), 0L)
  expect_length(intersect(KASP_CHECK_SITE, KASP_CHECK_AS), 0L)
})

test_that("every scheme carries exactly its family's check classes", {
  for (sc in kasp_schemes()$scheme) {
    chk <- kasp_criteria(scheme = sc)
    chk <- chk$label[chk$role == "check"]
    # C2: a scheme carries only the checks its chemistry obeys
    expect_true(all(chk %in% family_checks(sc)), info = sc)
    # ... and every marker-agnostic check is there (it holds for any primer pair)
    expect_true(all(KASP_CHECK_PRIMER %in% chk), info = sc)
    # the AS-only and site checks track the family, never the scheme
    info <- kasp_scheme(sc)
    as_family <- "as" %in% (KASP_FAMILY_CHECKS[[info$family]] %||% character(0))
    expect_equal(any(KASP_CHECK_AS %in% chk), as_family, info = sc)
  }
  # the multiplex rules are duplex-only
  expect_true(all(paste0("duplex:gate:", KASP_GATES_DUPLEX) %in% kasp_criteria(scheme = "duplex")$id))
  expect_false(any(grepl("duplex", kasp_criteria(scheme = "snp")$id)))
})

test_that("kasp_design_checks returns exactly the registry's checks, in order", {
  for (sc in kasp_schemes()$scheme) {
    crit <- kasp_criteria(scheme = sc)
    checks <- crit[crit$role == "check", ]
    ch <- kasp_design_checks(crit_mock(sc), scheme = sc)
    expect_equal(ch$criterion, checks$label, info = sc)
    expect_equal(ch$limit, checks$limit, info = sc)
    expect_true(all(ch$pass), info = sc) # a clean mock passes every check it owns
  }
})
