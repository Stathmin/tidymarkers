# LAMP as a marker scheme: F3/B3 outside, FIP/BIP inside, variant in the loop.

test_that("design_lamp_marker returns a 4-primer LAMP set with the variant in the loop", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_lamp_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "kasp_lamp")
  expect_named(m$primers, c("F3", "B3", "FIP", "BIP"))
  # inner (two-arm fusion) primers are long, the outer F3/B3 short (LAMP rule)
  expect_gte(nchar(m$primers$FIP), 40L)
  expect_lte(nchar(m$primers$FIP), 50L)
  expect_gte(nchar(m$primers$BIP), 40L)
  expect_lte(nchar(m$primers$BIP), 50L)
  expect_gte(nchar(m$primers$F3), 21L)
  expect_lte(nchar(m$primers$F3), 24L)
  # the four primers are screened jointly, not each one alone
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  expect_true(all(grepl("^[ACGT]+$", unlist(m$primers))))
  # amplicon within the LAMP window and straddling the variant
  product <- m$amplicon$end - m$amplicon$start + 1L
  expect_gte(product, 120L)
  expect_lte(product, 300L)
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  # the variant lies between F2 and B2 (the loop the inner primers read)
  expect_lte(m$loop$start, 300L)
  expect_gte(m$loop$end, 300L)
  # the inner primers are two-arm fusions, longer than the outer arms
  expect_gt(nchar(m$primers$FIP), nchar(m$primers$F3))
  expect_gt(nchar(m$primers$BIP), nchar(m$primers$B3))
  expect_output(print.kasp_lamp(m), "kasp_lamp")
  expect_output(print.kasp_lamp(m), "FIP:")
})

test_that("print.kasp_lamp names the four primers and the loop", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_lamp_marker(genome, snp, flank = 150L)
  expect_output(print.kasp_lamp(m), "variant loop:")
  expect_output(print.kasp_lamp(m), "B3:")
  expect_output(print.kasp_lamp(m), "BIP:")
})

test_that("scheme_def_lamp declares the LAMP stage consistently", {
  s <- scheme_def_lamp()
  expect_equal(s$scheme, "lamp")
  expect_equal(s$family, "lamp")
  expect_equal(s$loci, 1L)
  expect_true(all(s$gates %in% names(s$gate_short)))
  expect_true(all(s$gates %in% names(s$gate_limit)))
  expect_true(all(s$relax %in% s$gates))
  expect_true(all(s$gates %in% names(s$gate_floor_label) | !(s$gates %in% s$relax)))
  expect_identical(unname(s$tiers), c("loop-region specificity", "Tm set", "amplicon length"))
})
