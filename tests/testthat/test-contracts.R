# ZY12: design by contract - a broken invariant fails at the call site, with a
# classed error naming the function, instead of surfacing downstream as an NA.

expect_contract_error <- function(expr, fn) {
  err <- tryCatch(expr, error = identity)
  expect_s3_class(err, "tidymarkers_contract_error")
  expect_match(conditionMessage(err), fn, fixed = TRUE)
  err
}

no_contract_error <- function(expr) {
  !inherits(tryCatch(expr, error = identity), "tidymarkers_contract_error")
}

# The cheap valid target from test-target.R: no seed scan, no tool.
target_fixture <- function() {
  seq <- strrep("ACGT", 20L)
  list(
    genome = Biostrings::DNAStringSet(c(Chr1 = seq)),
    gene = list(chrom = "Chr1", symbol = "X", start = 1L, end = nchar(seq), strand = 1L),
    variants = tibble::tibble(
      id = "v", chrom = "Chr1", pos = 20L,
      ref_genome = list(c("T", "A")), consequence = "missense_variant"
    )
  )
}

# The minimal artefact from test-z8-outcomes.R: manifest + targets, no design.
outcome_fixture <- function(dir) {
  writeLines(
    '{"schema": "tidymarkers/artefact@1", "id": "0123456789abcdef"}',
    fs::path(dir, "artefact.json")
  )
  readr::write_csv(
    tibble::tibble(target = c("t1", "t2")),
    fs::path(dir, "targets.csv")
  )
  dir
}

# `R CMD check` runs this file against an installed build, where the dev-mode
# default is off (testthat does not put `pkgload` on `loadedNamespaces()`), so
# nothing would be checked.  The option is the documented escape hatch: each
# block that needs the contracts opens the gate for its own body.
contracts_on <- function() {
  # `.local_envir` is the caller's frame: the default would restore the option
  # as soon as this helper returns.
  withr::local_options(tidymarkers.contracts = TRUE, .local_envir = parent.frame())
}

test_that("the contracts option overrides the dev-mode default both ways", {
  contracts_on()
  expect_true(contracts_enabled())
  expect_identical(
    withr::with_options(list(tidymarkers.contracts = NULL), contracts_enabled()),
    in_dev_mode()
  )
  expect_false(withr::with_options(list(tidymarkers.contracts = FALSE), contracts_enabled()))
  expect_true(withr::with_options(list(tidymarkers.contracts = TRUE), contracts_enabled()))
})

test_that("marker_contract treats NA and a vector condition as violations", {
  expect_contract_error(marker_contract(NA, "Contract violated in f: unset"), "f")
  expect_contract_error(marker_contract(c(TRUE, FALSE), "Contract violated in f: vector"), "f")
  expect_invisible(marker_contract(TRUE, "not used"))
  # invisibility is kept by `identical()`, so this pins the returned value too
  expect_identical(marker_contract(TRUE, "ok"), TRUE)
})

test_that("marker_contract names its caller, not itself, in the error call", {
  # a direct caller whose own call is `wrapper()`: under a `sys.call(0)` default
  # the reported call would be `marker_contract(FALSE, "boom")` instead
  wrapper <- function() marker_contract(FALSE, "boom")
  err <- tryCatch(wrapper(), error = identity)

  expect_s3_class(err, "tidymarkers_contract_error")
  expect_identical(conditionCall(err), quote(wrapper()))
})

test_that("pick_target_snp fails its own precondition, not downstream", {
  contracts_on()
  f <- target_fixture()
  broken <- f$variants[, c("id", "chrom", "pos")]
  err <- expect_contract_error(
    pick_target_snp(broken, f$genome, f$gene, specificity_candidates = 0L),
    "pick_target_snp"
  )
  expect_match(conditionMessage(err), "variants", fixed = TRUE)
  # the call shown to the user is the function's own, not the helper's
  expect_match(paste(deparse(conditionCall(err)), collapse = " "), "pick_target_snp", fixed = TRUE)
  # the gate is off in production: the same input fails later (or not at all),
  # but never with the contract class
  expect_true(no_contract_error(withr::with_options(
    list(tidymarkers.contracts = FALSE),
    pick_target_snp(broken, f$genome, f$gene, specificity_candidates = 0L)
  )))
})

test_that("assemble_kasp_assay fails its own precondition, not downstream", {
  contracts_on()
  f <- snp_fixture()
  bad <- f$snp
  bad$pos <- 0L
  err <- expect_contract_error(assemble_kasp_assay(f$genome, bad, flank = 150L), "assemble_kasp_assay")
  expect_match(conditionMessage(err), "position", fixed = TRUE)
  expect_match(paste(deparse(conditionCall(err)), collapse = " "), "assemble_kasp_assay", fixed = TRUE)
  expect_true(no_contract_error(withr::with_options(
    list(tidymarkers.contracts = FALSE),
    assemble_kasp_assay(f$genome, bad, flank = 150L)
  )))
})

test_that("marker_record_outcome fails its own precondition, not downstream", {
  contracts_on()
  dir <- withr::local_tempdir()
  writeLines('{"schema": "tidymarkers/artefact@1"}', fs::path(dir, "artefact.json"))
  err <- expect_contract_error(marker_record_outcome(dir, "works"), "marker_record_outcome")
  expect_match(conditionMessage(err), "artefact id", fixed = TRUE)
  expect_match(paste(deparse(conditionCall(err)), collapse = " "), "marker_record_outcome", fixed = TRUE)
  # with the gate off the old silent behaviour is back: an id-less row
  withr::with_options(list(tidymarkers.contracts = FALSE), {
    out <- marker_record_outcome(dir, "works")
    expect_true(is.na(out$id[[1L]]))
  })
})

test_that("the contracts do not fire on a scheme that legitimately designs nothing", {
  contracts_on()
  withr::local_envvar(LANGUAGE = "en")
  f <- target_fixture()
  expect_warning(
    tgt <- pick_target_snp(f$variants, f$genome, f$gene,
      scheme = "primer", specificity_candidates = 0L
    ),
    "designable for scheme primer"
  )
  expect_identical(tgt$id, "v")
})

test_that("a valid pick_target_snp keeps its return shape", {
  contracts_on()
  f <- target_fixture()
  tgt <- pick_target_snp(f$variants, f$genome, f$gene, specificity_candidates = 0L)
  expect_s3_class(tgt, "tbl_df")
  expect_equal(nrow(tgt), 1L)
  expect_identical(tgt$id, "v")
  expect_setequal(
    names(tgt),
    c("id", "chrom", "pos", "consequence", "ref_len", "ref", "alt", "transition", "n_3p", "seed_hits")
  )
  expect_null(attr(tgt, "marker"))
})

test_that("a valid assay keeps the completeness fields the contract checks", {
  contracts_on()
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  a <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp, flank = 150L))
  expect_s3_class(a, "marker_assay")
  expect_identical(a$snp$ref, f$ref)
  expect_true(all(c("scheme", "assay", "snp", "tails") %in% names(a)))
  expect_true(nzchar(a$assay$as_ref[[1L]]))
  expect_setequal(names(a$tails), c("AS_ref", "AS_alt", "COM"))
  expect_true(is.finite(a$assay$amplicon[[1L]]))
})

test_that("a recorded outcome is readable back through marker_outcomes()", {
  contracts_on()
  dir <- withr::local_tempdir()
  dir <- outcome_fixture(dir)
  out <- marker_record_outcome(dir, "works")
  expect_equal(nrow(out), 1L)
  expect_identical(out$id, "0123456789abcdef")
  expect_false(any(out$calibrated))
  expect_identical(marker_outcomes(dir), out)
})
