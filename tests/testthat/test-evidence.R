# ZW4: the evidence ladder.  A threshold's provenance is an ordered rung, not a
# boolean, and the rung travels to the verdict.  No scientific source was ever
# claimed for these numbers, so every row on this base is `rule_of_thumb`; the
# tests below pin that (no invented provenance) as well as the ladder itself.

# a minimal one-row registry entry, so `criterion_tbl()` is tested without the
# whole KASP registry
ev_one <- function(calibrated) {
  tidymarkers:::criterion_tbl(
    id = "probe", scope = "kasp", stage = "assay", role = "check",
    label = "probe", level = "alignment", kind = "numeric",
    bounds = list(c(1, 2)), units = "x", source = "package assumption",
    calibrated = calibrated, evaluator = "probe"
  )
}

# An assay that passes every acceptance check.  `good_assay()` lives in
# test-kasp-checks.R, which is KNOWN RED on this base and must not be edited, so
# the fixture is repeated here rather than forked from it.
ev_good_assay <- function() {
  p <- kasp_primers()
  list(
    assay = tibble::tibble(
      amplicon = 90L, as_ref = p[["AS_ref"]], as_alt = p[["AS_alt"]], com = p[["COM"]],
      com_tm = 60, as_gc = 0.5, com_gc = 0.5, n_var = 0L, as_n_3p = 0L,
      as_n_var = 0L, as_dropout = 0, as_foot_dropout = 0, com_dropout = 0
    ),
    thermo = tibble::tibble(
      name = names(p), tm = c(57, 57, 60), gc = c(0.5, 0.5, 0.5),
      hairpin_dg = c(-1, -1, -1)
    ),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3), b = rep(names(p), times = 3),
      dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)
    ),
    amplicons = tibble::tibble(
      primer = c("AS_ref", "AS_alt"),
      n_products = c(1L, 0L), width = c(90L, NA_integer_)
    ),
    ideogram = list(hits = tibble::tibble(
      primer = c("AS_ref", "AS_alt", "COM")
    )),
    efficiency = tibble::tibble(discrimination = 0.9)
  )
}

test_that("the evidence ladder is strictly ordered", {
  ladder <- tidymarkers:::CRITERION_EVIDENCE
  expect_equal(
    ladder,
    c("rule_of_thumb", "vendor_guideline", "primary_literature", "fitted_here")
  )
  ev <- tidymarkers:::criterion_evidence(ladder)
  expect_s3_class(ev, "ordered")
  expect_equal(levels(ev), ladder)
  # an ordered factor compares by rung, not by string
  expect_true(
    tidymarkers:::criterion_evidence("fitted_here") >
      tidymarkers:::criterion_evidence("rule_of_thumb")
  )
  expect_false(
    tidymarkers:::criterion_evidence("rule_of_thumb") >
      tidymarkers:::criterion_evidence("fitted_here")
  )
  expect_true(
    tidymarkers:::criterion_evidence("vendor_guideline") <
      tidymarkers:::criterion_evidence("primary_literature")
  )
  # unset is unset, never silently the weakest rung
  expect_true(is.na(tidymarkers:::criterion_evidence(NULL)))
  expect_true(is.na(tidymarkers:::criterion_evidence(NA)))
})

test_that("an unknown evidence value aborts, naming it and the allowed set", {
  err <- rlang::catch_cnd(tidymarkers:::criterion_evidence("peer_reviewed"))
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_match(msg, "peer_reviewed", fixed = TRUE)
  for (rung in tidymarkers:::CRITERION_EVIDENCE) {
    expect_match(msg, rung, fixed = TRUE)
  }
})

test_that("the deprecated `calibrated` flag maps onto the ladder", {
  expect_equal(as.character(ev_one(TRUE)$evidence), "fitted_here")
  expect_equal(as.character(ev_one(FALSE)$evidence), "rule_of_thumb")
  expect_true(is.na(ev_one(NA)$evidence))
  # the old column is still emitted and still says what it said
  expect_true(ev_one(TRUE)$calibrated)
  expect_false(ev_one(FALSE)$calibrated)
  # an explicit rung wins over the flag
  x <- tidymarkers:::criterion_tbl(
    id = "probe", scope = "kasp", stage = "assay", role = "check",
    label = "probe", level = "alignment", kind = "numeric",
    bounds = list(c(1, 2)), units = "x", source = "package assumption",
    calibrated = FALSE, evaluator = "probe", evidence = "vendor_guideline"
  )
  expect_equal(as.character(x$evidence), "vendor_guideline")
  expect_s3_class(x$evidence, "ordered")
  expect_equal(levels(x$evidence), tidymarkers:::CRITERION_EVIDENCE)
  # an unknown rung is a bug, not a silent fallback
  expect_error(
    tidymarkers:::criterion_tbl(
      id = "probe", scope = "kasp", stage = "assay", role = "check",
      label = "probe", level = "alignment", kind = "numeric",
      bounds = list(c(1, 2)), units = "x", source = "package assumption",
      evaluator = "probe", evidence = "guesswork"
    ),
    "guesswork"
  )
})

test_that("no criterion row claims provenance it does not have", {
  spec <- tidymarkers:::criterion_spec(KASP_LIMITS)
  expect_true("evidence" %in% names(spec))
  expect_s3_class(spec$evidence, "ordered")
  expect_equal(levels(spec$evidence), tidymarkers:::CRITERION_EVIDENCE)
  expect_false(anyNA(spec$evidence))
  expect_setequal(as.character(spec$evidence), "rule_of_thumb")
  # and nothing is promoted to a rung no source backs
  expect_false(any(as.character(spec$evidence) %in% c(
    "vendor_guideline", "primary_literature", "fitted_here"
  )))
  # the deprecated column is still there, unchanged (all FALSE here)
  expect_true("calibrated" %in% names(spec))
  expect_false(anyNA(spec$calibrated))
  expect_false(any(spec$calibrated))
})

test_that("the evidence rung reaches the verdict, without a new verdict string", {
  ch <- marker_design_checks(ev_good_assay())
  expect_equal(evaluation_verdict(ch), "usable")
  expect_equal(evaluation_evidence(ch), "rule_of_thumb")

  # a marker that fails every check has nothing passing, hence no rung
  doomed <- tibble::tibble(criterion = ch$criterion, pass = FALSE, status = "fail")
  expect_equal(evaluation_verdict(doomed), "doomed")
  expect_true(is.na(evaluation_evidence(doomed)))
  # no rows at all is likewise no evidence
  expect_true(is.na(evaluation_evidence(tibble::tibble())))

  # the verdict vocabulary itself is untouched: every branch of
  # evaluation_verdict() - all pass, a fail, an unevaluated check, nothing that
  # actually passed, and no rows at all - still answers from exactly this closed
  # set, so the evidence rung is a second axis and not a fifth verdict string
  vocabulary <- c("usable", "doomed", "unknown")
  no_status <- ch
  no_status$status <- NA_character_
  # passes on paper only: neither not_applicable nor unknown is a real pass
  nothing_passed <- ch
  nothing_passed$status <- "not_applicable"
  verdicts <- vapply(
    list(ch, doomed, no_status, nothing_passed, tibble::tibble()),
    evaluation_verdict, character(1)
  )
  expect_equal(
    setdiff(verdicts, vocabulary), character(),
    info = paste("verdicts:", paste(verdicts, collapse = ", "))
  )
  # ...and the branches above really do reach more than one of them, so the
  # subset check above is not passing on a single degenerate answer
  expect_gte(length(unique(verdicts)), 2L)
})

test_that("evaluation_evidence takes the weakest rung among passing checks only", {
  crit <- tidymarkers:::criterion_spec(KASP_LIMITS)
  # a fabricated registry where one check is strong and another weak
  lab <- crit[crit$role == "check", ]
  strong <- lab$id[[which(lab$level == "site")[1L]]]
  weak <- lab$id[[which(lab$level == "alignment")[1L]]]
  spec <- crit
  spec$evidence[match(strong, spec$id)] <- "primary_literature"
  spec$evidence <- tidymarkers:::criterion_evidence(as.character(spec$evidence))
  ch <- tibble::tibble(
    criterion = c(strong, weak, strong),
    status = c("pass", "pass", "fail")
  )
  # the pass-only row is the strong one...
  expect_equal(
    evaluation_evidence(ch[1L, ], criteria = spec),
    "primary_literature"
  )
  # ...but a weak check that also passes drags it down to the minimum
  expect_equal(evaluation_evidence(ch, criteria = spec), "rule_of_thumb")
  # a failing strong row vouches for nothing
  only_fail <- ch[3L, ]
  expect_true(is.na(evaluation_evidence(only_fail, criteria = spec)))
  # not_applicable / unknown are ignored rather than counted as passes
  skipped <- tibble::tibble(
    criterion = c(strong, weak), status = c("not_applicable", "unknown")
  )
  expect_true(is.na(evaluation_evidence(skipped, criteria = spec)))
})
