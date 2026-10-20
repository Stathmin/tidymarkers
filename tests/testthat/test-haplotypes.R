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
  # 260 rides the ALT haplotypes exactly -> linked to the target; 240 (H1) is not
  # an ALT ride, and is reported per haplotype
  expect_equal(unique(hv$linked$target), 250L)
  expect_equal(unique(hv$linked$linked), 260L)
  expect_setequal(hv$linked$haplotype, c("H2", "H3"))
  expect_true(all(hv$linked$shared))
  expect_setequal(hv$alt, c("H2", "H3"))
  expect_equal(hv$ref, "H1")
})

test_that("only the ALT haplotypes' shared variants are baked", {
  m <- hap_fixture()
  b210 <- substring(m$genome[["Chr1"]], 210L, 210L)
  haps <- dplyr::bind_rows(
    m$haps,
    tibble::tibble(haplotype = "H2", freq = 0.3, chrom = "Chr1", pos = 210L, ref = b210, alt = "A")
  )
  hv <- haplotype_variants(haps, m$snp)
  # 260 is carried by both ALT haplotypes -> shared; 210 by H2 alone -> not
  expect_true(all(hv$linked$shared[hv$linked$linked == 260L]))
  expect_false(any(hv$linked$shared[hv$linked$linked == 210L]))
})

test_that("design_as_primers bakes a linked indel length-aware", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(paste(rep("ACGT", 40L), collapse = ""))
  plain <- design_as_primers(W, 40L, "T", "A", strands = "+", lengths = 22L, stab_ats = NA)
  b <- as.character(Biostrings::subseq(W, 30L, 30L))
  linked <- tibble::tibble(offset = 30L, ref = b, alt = "")
  baked <- design_as_primers(W, 40L, "T", "A", strands = "+", lengths = 22L, stab_ats = NA, linked = linked)
  expect_equal(nchar(baked$alt_primer), 22L)
  expect_false(identical(baked$alt_primer, plain$alt_primer)) # the deletion shifted it
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

test_that("any marker can carry population haplotypes and report ALT frequency", {
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "G")
  h <- tibble::tibble(
    haplotype = c("H1", "H1", "H2", "H0"),
    freq = c(0.55, 0.55, 0.25, 0.20),
    chrom = "Chr1", pos = c(300L, 310L, 300L, 320L),
    ref = c("A", "C", "A", "G"), alt = c("G", "T", "G", "A")
  )
  m <- marker_attach_haplotypes(list(snp = snp, amplicon = tibble::tibble(start = 250L, end = 350L)), h, snp)
  expect_equal(m$population$alt_freq, 0.8)
  expect_equal(m$population$n_haplotypes, 3L)
  expect_equal(m$population$n_alt_haplotypes, 2L)
  ch <- marker_population_checks(m)
  expect_equal(ch$value[ch$name == "ALT freq"], "80%")
  expect_equal(ch$status[ch$name == "ALT freq"], "ok")
})

test_that("assemble_kasp_assay accepts haplotypes and reports the hits", {
  skip_if_no_oligoarrayaux()
  m <- hap_fixture()
  a <- suppressWarnings(assemble_kasp_assay(m$genome, m$snp, haplotypes = m$haps, flank = 150L, cores = 1L))
  if (inherits(a, "marker_doomed")) skip("the synthetic haplotype did not assemble")
  expect_true(!is.null(a$haplotypes))
  expect_equal(nrow(a$haplotype_hits), 3L)
  expect_true(all(c("haplotype", "freq", "target_allele", "primer", "binds") %in% names(a$haplotype_hits)))
  a <- marker_attach_haplotypes(a, m$haps, m$snp)
  pop <- marker_population_checks(a)
  expect_true(all(c("ALT freq", "haplotypes", "ALT haplotypes", "hap covered") %in% pop$name))
  # the co-segregating variant was baked into the ALT haplotype
  expect_true(!is.null(a$selection$linked) && nrow(a$selection$linked) >= 1L)
  # the layered haplotype track renders, and the full plot still assembles
  expect_s3_class(haplotype_layers(a, m$snp, span = c(200, 300)), "ggplot")
  v <- haplotype_variants(m$haps, m$snp)$variants
  expect_s3_class(plot_primers(a, v, m$snp, genome = m$genome), "patchwork")
})

test_that("marker_run_context stores variants and haplotypes unchanged", {
  g <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  v <- tibble::tibble(chrom = "Chr1", pos = 10L, alt = "A")
  h <- tibble::tibble(
    haplotype = c("H1", "H2"), freq = c(0.4, 0.6), chrom = "Chr1",
    pos = c(10L, 20L), ref = c("C", "G"), alt = c("A", "T")
  )
  ctx <- marker_run_context(g, variants = v, population = h)
  expect_s3_class(ctx, "marker_run_context")
  expect_identical(ctx$genome, g)
  expect_identical(ctx$variants, v)
  expect_identical(ctx$population, h)
  expect_identical(ctx$haplotypes, h)
})

test_that("a context with no population leaves the frequency unknown, not 1", {
  g <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  ctx <- marker_run_context(g)
  expect_null(ctx$population)
  expect_true(is.na(marker_run_context_frequency(ctx)))
  expect_false(isTRUE(marker_run_context_frequency(ctx) == 1))
})

test_that("a context with known haplotypes attaches $population on design", {
  f <- hap_fixture()
  ctx <- marker_run_context(f$genome, population = f$haps)
  args <- marker_run_context_args(ctx)
  expect_identical(args$haplotypes, f$haps)
  # the existing attach path computes the frequency from the haplotypes
  marker <- marker_attach_haplotypes(list(snp = f$snp), args$haplotypes, f$snp)
  expect_equal(marker$population$alt_freq, 0.5)
})
