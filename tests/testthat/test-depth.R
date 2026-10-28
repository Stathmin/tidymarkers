# Depth weighting (roadmap Y5.3): DP/AD -> weight/confident, and the solver gate.

test_that("a low-depth homozygote is not confident, a deep one is", {
  low <- variant_confidence(tibble::tibble(
    chrom = "Chr1", pos = 1L, depth = 1L, alt_depth = 0L, genotype = "0/0"
  ))
  expect_false(low$confident)
  expect_lt(low$weight, 0.5)

  deep <- variant_confidence(tibble::tibble(
    chrom = "Chr1", pos = 1L, depth = 10L, alt_depth = 0L, genotype = "0/0"
  ))
  expect_true(deep$confident)
})

test_that("weight is strictly increasing in depth for a homozygote", {
  x <- variant_confidence(tibble::tibble(
    chrom = "Chr1", pos = 1L, depth = c(1, 2, 3, 5, 10, 20), alt_depth = 0L
  ))
  expect_true(all(diff(x$weight) > 0))
})

test_that("the heterozygote balance term down-weights an imbalanced call", {
  balanced <- variant_confidence(tibble::tibble(
    chrom = "Chr1", pos = 1L, depth = 20L, alt_depth = 10L, genotype = "0/1"
  ))
  imbalanced <- variant_confidence(tibble::tibble(
    chrom = "Chr1", pos = 1L, depth = 20L, alt_depth = 18L, genotype = "0/1"
  ))
  expect_true(balanced$confident)
  expect_lt(imbalanced$weight, balanced$weight)
  expect_false(imbalanced$confident)
})

test_that("read_variant_depth() reads DP/AD/GT and falls back to INFO/DP", {
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:DP\t0/0:10,0:10",
    "1\t200\t.\tC\tT\t.\t.\t.\tGT:AD:DP\t0/1:5,5:10"
  ), vcf)
  d <- read_variant_depth(vcf)
  expect_equal(nrow(d), 2L)
  expect_equal(d$chrom, c("Chr1", "Chr1"))
  expect_equal(d$pos, c(100L, 200L))
  expect_equal(d$depth, c(10, 10))
  expect_equal(d$alt_depth, c(0, 5))
  expect_equal(d$genotype, c("0/0", "0/1"))

  sites <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t300\t.\tA\tG\t.\t.\tDP=7"
  ), sites)
  s <- read_variant_depth(sites)
  expect_equal(s$depth, 7)
  expect_true(is.na(s$alt_depth))
  expect_true(is.na(s$genotype))
})

test_that("the depth gate keeps confident and unknown targets, drops low ones", {
  fx <- snp_fixture()
  region <- tibble::tibble(chrom = "Chr1", start = 1L, end = 300L)
  v <- tibble::tibble(
    chrom = "Chr1", pos = c(100L, 200L, 250L), ref = "A", alt = "G"
  )
  d <- tibble::tibble(
    chrom = "Chr1", pos = c(100L, 200L),
    depth = c(1L, 20L), alt_depth = c(0L, 0L), genotype = c("0/0", "0/0")
  )
  p <- marker_plan(fx$genome, region = region) + marker_context(variants = v, depth = d)
  tg <- region_targets(p)
  expect_true(all(c(100L, 200L, 250L) %in% tg$pos))
  kept <- filter_targets_by_depth(tg, p$depth)
  expect_equal(sort(kept$pos), c(200L, 250L))
})

test_that("the solver rejects a plan whose every target is low-confidence", {
  fx <- snp_fixture()
  region <- tibble::tibble(chrom = "Chr1", start = 1L, end = 300L)
  v <- tibble::tibble(chrom = "Chr1", pos = 100L, ref = "A", alt = "G")
  d <- tibble::tibble(
    chrom = "Chr1", pos = 100L, depth = 1L, alt_depth = 0L, genotype = "0/0"
  )
  p <- marker_plan(fx$genome, region = region) + marker_context(variants = v, depth = d)
  expect_error(marker_solve(p, cores = 1L, progress = FALSE), "depth gate")
})

test_that("variant_confidence() normalises chrom so a '1' depth table gates Chr1 targets", {
  d1 <- variant_confidence(tibble::tibble(
    chrom = "1", pos = 100L, depth = 1L, alt_depth = 0L, genotype = "0/0"
  ))
  d2 <- variant_confidence(tibble::tibble(
    chrom = "Chr1", pos = 100L, depth = 1L, alt_depth = 0L, genotype = "0/0"
  ))
  expect_equal(d1$chrom, "Chr1")
  targets <- tibble::tibble(chrom = "Chr1", pos = 100L, ref = "A", alt = "G")
  expect_error(filter_targets_by_depth(targets, d1), "depth gate")
  expect_error(filter_targets_by_depth(targets, d2), "depth gate")
})

test_that("filter_targets_by_depth() normalises the target chrom against a normalised depth table", {
  d <- variant_confidence(tibble::tibble(
    chrom = "1", pos = 100L, depth = 1L, alt_depth = 0L, genotype = "0/0"
  ))
  targets <- tibble::tibble(chrom = "1", pos = 100L, ref = "A", alt = "G")
  expect_error(filter_targets_by_depth(targets, d), "depth gate")
})

test_that("an empty plan aborts with 'no targets', not the depth gate", {
  fx <- snp_fixture()
  d <- tibble::tibble(
    chrom = "Chr1", pos = 100L, depth = 1L, alt_depth = 0L, genotype = "0/0"
  )
  p <- marker_plan(fx$genome) + marker_context(depth = d)
  expect_error(marker_solve(p, cores = 1L, progress = FALSE), "no targets")
})
