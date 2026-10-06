# Haplotype model: per-position frequencies and co-segregation from haplotypes,
# and which AS primer each haplotype is expected to bind.

# target at 250 on a fixed genome; H1 carries REF, H2/H3 carry ALT at 250 and
# both share a variant at 260 (so 260 co-segregates with the ALT allele).
hap_fixture <- function() {
  W <- withr::with_seed(4L, paste(sample(c("A", "C", "G", "T"), 500L, TRUE), collapse = ""))
  ref <- substring(W, 250L, 250L)
  snp <- tibble::tibble(chrom = "Chr1", pos = 250L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
  b240 <- substring(W, 240L, 240L)
  b260 <- substring(W, 260L, 260L)
  haps <- tibble::tibble(
    haplotype = c("H1", "H2", "H2", "H3", "H3"),
    freq = c(0.5, 0.3, 0.3, 0.2, 0.2),
    chrom = "Chr1", pos = c(240L, 250L, 260L, 250L, 260L),
    ref = c(b240, ref, b260, ref, b260),
    alt = c("A", snp$alt, "A", snp$alt, "A")
  )
  list(genome = Biostrings::DNAStringSet(c(Chr1 = W)), snp = snp, haps = haps)
}

test_that("haplotypes give per-position frequencies and co-segregation", {
  m <- hap_fixture()
  hv <- haplotype_variants(m$haps, m$snp)
  # the position at 260 is carried by H2+H3 -> freq 0.3 + 0.2
  expect_equal(hv$variants$pop_freq[hv$variants$pos == 260L], 0.5, tolerance = 1e-9)
  expect_equal(hv$variants$pop_freq[hv$variants$pos == 240L], 0.5, tolerance = 1e-9)
  # 260 rides the ALT haplotypes exactly -> linked to the target; 240 does not
  expect_equal(hv$linked$target, 250L)
  expect_equal(hv$linked$linked, 260L)
  expect_setequal(hv$alt, c("H2", "H3"))
  expect_equal(hv$ref, "H1")
})

test_that("haplotype_hits says which AS primer each haplotype should bind", {
  m <- hap_fixture()
  hv <- haplotype_variants(m$haps, m$snp)
  hits <- haplotype_hits(m$haps, m$snp, as_start = 235L, as_end = 255L, alt_hap = hv$alt)
  expect_equal(nrow(hits), 3L)
  expect_equal(hits$primer[hits$haplotype == "H1"], "AS_ref")
  expect_equal(hits$primer[hits$haplotype %in% c("H2", "H3")], c("AS_alt", "AS_alt"))
  # H2/H3 carry the 260 variant, which sits OUTSIDE this footprint -> clean bind
  expect_true(all(hits$binds[hits$haplotype %in% c("H2", "H3")]))
})

test_that("assemble_kasp_assay accepts haplotypes and reports the hits", {
  skip_if_no_oligoarrayaux()
  m <- hap_fixture()
  a <- suppressWarnings(assemble_kasp_assay(m$genome, m$snp, haplotypes = m$haps, flank = 150L, cores = 1L))
  if (inherits(a, "kasp_doomed")) skip("the synthetic haplotype did not assemble")
  expect_true(!is.null(a$haplotypes))
  expect_equal(nrow(a$haplotype_hits), 3L)
  expect_true(all(c("haplotype", "freq", "target_allele", "primer", "binds") %in% names(a$haplotype_hits)))
  # the co-segregating variant was baked into the ALT haplotype
  expect_true(!is.null(a$selection$linked) && nrow(a$selection$linked) >= 1L)
  # the layered haplotype track renders, and the full plot still assembles
  expect_s3_class(haplotype_layers(a, m$snp, span = c(200, 300)), "ggplot")
  v <- haplotype_variants(m$haps, m$snp)$variants
  expect_s3_class(plot_primers(a, v, m$snp, genome = m$genome), "patchwork")
})
