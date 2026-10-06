# Sanger as a marker scheme: one flanking pair around a well-placed variant.

test_that("design_sanger_marker returns a unique amplicon with the variant inside", {
  W <- bseq(1200L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 600L, 600L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 600L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_sanger_marker(genome, snp, flank = 500L)
  expect_s3_class(m, "kasp_sanger")
  # the variant is inside the amplicon and clear of both read starts
  expect_gt(m$amplicon$start, 0L)
  expect_lt(m$amplicon$start, snp$pos)
  expect_gt(m$amplicon$end, snp$pos)
  expect_gte(snp$pos - m$amplicon$start, 50L)
  expect_gte(m$amplicon$end - snp$pos, 50L)
  # the amplicon sits in the Sanger window and the primers are plain ACGT
  expect_gte(m$primers$product, 400L)
  expect_lte(m$primers$product, 800L)
  expect_true(grepl("^[ACGT]+$", m$primers$forward) && grepl("^[ACGT]+$", m$primers$reverse))
  expect_equal(nrow(in_silico_pcr(m$primers$forward, m$primers$reverse, genome)), 1L)
})

test_that("print.kasp_sanger summarises the marker", {
  m <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 10L, ref = "A", alt = "G"),
      primers = list(forward = "ACGTACGTACGTACGTAC", reverse = "TGCATGCATGCATGCATG", product = 600L),
      amplicon = tibble::tibble(start = 5L, end = 605L)
    ),
    class = "kasp_sanger"
  )
  expect_output(print(m), "kasp_sanger")
  expect_output(print(m), "600 bp")
})

test_that("design_sanger_marker declines when no pair keeps the variant off the read start", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(300L, 3L)))
  ref <- substr(as.character(genome$Chr1), 25L, 25L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 25L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  expect_message(design_sanger_marker(genome, snp, flank = 20L), "read starts")
})

test_that("scheme_def_sanger describes the chemistry and its gates", {
  d <- scheme_def_sanger()
  expect_identical(d$scheme, "sanger")
  expect_identical(d$family, "sanger")
  expect_identical(d$loci, 1L)
  expect_setequal(d$gates, names(d$gate_short))
  expect_setequal(d$gates, names(d$gate_limit))
  expect_true(all(d$relax %in% d$gates))
  expect_true(all(d$gates %in% c(names(d$gate_floor_label), names(d$gate_hard_label))))
  expect_length(d$tiers, 3L)
})
