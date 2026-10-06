# Geometry invariants: the shared interval helpers are the single source of
# truth for where primers sit, so these lock the strand conventions down.

test_that("primer_interval: the AS primer's 3' end is the SNP on both strands", {
  pos <- 1000L
  L <- 25L
  plus <- primer_interval(pos, "+", L, role = "AS")
  minus <- primer_interval(pos, "-", L, role = "AS")
  expect_equal(IRanges::end(plus), pos) # "+" extends downstream from the SNP
  expect_equal(IRanges::start(plus), pos - L + 1L)
  expect_equal(IRanges::start(minus), pos) # "-" extends upstream from the SNP
  expect_equal(IRanges::end(minus), pos + L - 1L)
})

test_that("primer_interval: the COM primer sits on the opposite side of the SNP", {
  pos <- 1000L
  L <- 22L
  gap <- 30L
  plus <- primer_interval(pos, "+", L, gap = gap, role = "COM")
  minus <- primer_interval(pos, "-", L, gap = gap, role = "COM")
  expect_gt(IRanges::start(plus), pos) # both fully downstream
  expect_lt(IRanges::end(minus), pos) # or fully upstream
  # the gap is measured from the SNP to the near edge
  expect_equal(IRanges::start(plus), pos + gap)
  expect_equal(IRanges::end(minus), pos - gap)
})

test_that("as_3p_interval excludes the SNP and spans window-1 bases", {
  expect_equal(c(
    IRanges::start(as_3p_interval(1000L, "+", 5L)),
    IRanges::end(as_3p_interval(1000L, "+", 5L))
  ), c(996L, 999L))
  expect_equal(c(
    IRanges::start(as_3p_interval(1000L, "-", 5L)),
    IRanges::end(as_3p_interval(1000L, "-", 5L))
  ), c(1001L, 1004L))
})

test_that("count_variants_in counts points and drops the excluded SNP", {
  iv <- IRanges::IRanges(996L, 999L)
  expect_equal(count_variants_in(iv, c(996L, 998L, 1000L, 1005L), exclude = 1000L), 2L)
  expect_equal(count_variants_in(iv, NULL), NA_integer_)
})

test_that("a split multi-allelic site counts (and weighs) as one variant", {
  iv <- IRanges::IRanges(996L, 999L)
  # a VCF multi-allelic row arrives as one row per ALT: same position, two ALTs
  expect_equal(count_variants_in(iv, c(997L, 997L, 998L)), 2L)
  # the site's weights sum, so its total non-reference frequency is counted once
  expect_equal(variant_load_in(iv, c(997L, 997L, 998L), c(0.1, 0.2, 1)), 1.3)
})

test_that("strand_gc_hint reports the GC of each flank", {
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(strrep("A", 25L), "T", strrep("G", 25L))))
  h <- strand_gc_hint(genome, "chr1", 26L, flank = 25L)
  expect_equal(h$gc, c(0, 1))
  expect_lt(
    h$distance_to_45[h$side == "upstream"],
    h$distance_to_45[h$side == "downstream"]
  )
})
