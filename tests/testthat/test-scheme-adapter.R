# ZA6 step 3 / ZA7: every scheme is driven through the same five-step seam.
# A scheme is behaviour, not a row: the adapter gives design/evaluate/tracks/
# decision/readout the same names for every chemistry, delegating to the
# implementation the pipeline already uses so the two cannot drift.

test_that("every scheme with a designer resolves through the adapter", {
  schemes <- per_locus_schemes()
  expect_gt(length(schemes), 0L)
  for (s in schemes) {
    expect_true(is.function(scheme_designer(s)), info = s)
    a <- marker_scheme_adapter(s)
    expect_s3_class(a, "marker_scheme_adapter")
    expect_identical(a$scheme, s)
    expect_setequal(
      names(a),
      c("scheme", "design", "evaluate", "tracks", "decision", "readout")
    )
    for (step in c("design", "evaluate", "tracks", "decision", "readout")) {
      expect_true(is.function(a[[step]]), info = paste(s, step))
    }
  }
})

test_that("the adapter refuses an unknown scheme", {
  expect_error(marker_scheme_adapter("not_a_scheme"), "Unknown marker scheme")
})

test_that("a scheme without a designer still adapts, and designs nothing", {
  # `primer` is the marker-agnostic pair subset (C13): it is registered and its
  # criteria are evaluable, but it has no locus-level designer, so design()
  # returns NULL rather than inventing an allele assay.
  expect_true("primer" %in% marker_schemes()$scheme)
  expect_null(scheme_designer("primer"))
  a <- marker_scheme_adapter("primer")
  expect_null(a$design(NULL, NULL))
})

test_that("the five steps compose on a real design and NA never reads as pass", {
  skip_if_not_installed("tidygraph")
  f <- snp_fixture()
  a <- marker_scheme_adapter("snp")
  m <- suppressWarnings(suppressMessages(
    a$design(f$genome, f$snp, flank = 150L, cores = 1L)
  ))
  skip_if(is.null(m), "no KASP design at this fixture")

  ev <- a$evaluate(m)
  expect_s3_class(ev, "tbl_df")
  expect_true(all(c("criterion", "observed", "status", "pass") %in% names(ev)))
  expect_true(is.factor(ev$status))
  # the ZA7 rule: a criterion with no observation is `unknown`, never `pass`.
  # Two rows are deliberately textual (`target site class`, `AS 3' mismatch base
  # pair`): their observation is a word in `detail`, so a missing numeric
  # `observed` is not missing data - but a pass on such a row must carry it.
  textual <- c("target site class", "AS 3' mismatch base pair")
  unobserved <- is.na(ev$observed) & !ev$criterion %in% textual
  expect_false(any(ev$status == "pass" & unobserved))
  txt <- ev[ev$criterion %in% textual, , drop = FALSE]
  expect_true(all(!is.na(txt$detail[txt$status == "pass"])))

  tr <- a$tracks(m)
  expect_s3_class(tr, "tbl_df")
  expect_gt(nrow(tr), 0L)
  expect_true(all(c("feature", "start", "end", "kind") %in% names(tr)))

  dg <- a$decision(m)
  expect_s3_class(dg, "tbl_graph")

  ro <- a$readout(m)
  expect_type(ro, "list")
  expect_identical(ro$scheme, "snp")
  expect_type(ro$simulated, "logical")
})
