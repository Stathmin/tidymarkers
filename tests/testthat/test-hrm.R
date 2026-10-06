# HRM as a marker scheme: one flanking pair over a short amplicon.

test_that("design_hrm_marker returns a short flanking pair with the variant inside", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_hrm_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "kasp_hrm")
  # the variant is inside the amplicon and the amplicon is short (HRM wants a
  # short melt domain), with the primers flanking it
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  expect_lte(m$primers$product, 150L)
  expect_true(grepl("^[ACGT]+$", m$primers$forward) && grepl("^[ACGT]+$", m$primers$reverse))
  expect_output(print(m), "kasp_hrm")
})

test_that("the hrm scheme is registered with its own decision stage and is designable", {
  expect_true("hrm" %in% kasp_schemes()$scheme)
  s <- kasp_decision_schema("hrm")
  expect_true("hrm" %in% s$nodes$stage)
  expect_true(is_target_designable(structure(list(), class = "kasp_hrm")))
})

test_that("plot_hrm_amplicon draws the marker and rejects other objects", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_hrm_marker(genome, snp, flank = 150L)
  expect_s3_class(plot_hrm_amplicon(m), "ggplot")
  expect_error(plot_hrm_amplicon(list()))
})
