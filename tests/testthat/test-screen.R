# The cheap target screen: promising / doomed, and the linked-variant forgiveness.

clean_genome <- function(seed = 1L, n = 400L) {
  withr::with_seed(seed, {
    w <- paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
    Biostrings::DNAStringSet(c(Chr1 = w))
  })
}

test_that("screen_snp_targets returns the canonical placement score", {
  genome <- clean_genome()
  s <- as.character(genome[["Chr1"]])
  ref <- substr(s, 200L, 200L)
  vt <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
  sc <- screen_snp_targets(vt, genome, lengths = 20:27, strands = c("+", "-"))
  expect_equal(nrow(sc), 1L)
  expect_equal(sc$n_tried, 16L)
  expect_match(sc$score, "^[0-9]+/16$")
  expect_true(sc$verdict %in% c("promising", "doomed"))
  expect_s3_class(plot_snp_screen(sc), "ggplot")
})

test_that("a GC desert (AT-only) is doomed", {
  s <- paste(rep(c("A", "T"), 200L), collapse = "")
  genome <- Biostrings::DNAStringSet(c(Chr1 = s))
  vt <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = substr(s, 200L, 200L), alt = "C")
  sc <- screen_snp_targets(vt, genome)
  expect_equal(sc$verdict, "doomed")
  expect_equal(sc$n_ok, 0L)
})

test_that("a linked 3' variant is forgiven, an unlinked one is not", {
  genome <- clean_genome(seed = 5L)
  s <- as.character(genome[["Chr1"]])
  at <- function(p) {
    b <- substr(s, p, p)
    setdiff(c("A", "C", "G", "T"), b)[1L]
  }
  vt <- tibble::tibble(
    chrom = "Chr1", pos = c(200L, 199L, 201L),
    ref = c(substr(s, 200L, 200L), substr(s, 199L, 199L), substr(s, 201L, 201L)),
    alt = c(at(200L), at(199L), at(201L))
  )
  hard <- screen_snp_targets(vt, genome)
  hard_t <- hard[hard$pos == 200L, ]
  # neighbours sit in both strands' 3' windows -> some placements fail on var_3p
  expect_match(hard_t$reasons, "var_3p")

  linked <- tibble::tibble(target = 200L, linked = c(199L, 201L))
  soft <- screen_snp_targets(vt, genome, linked = linked)
  soft_t <- soft[soft$pos == 200L, ]
  expect_false(grepl("var_3p", soft_t$reasons))
  expect_gte(soft_t$n_ok, hard_t$n_ok)
})

test_that("complementary alleles have no usable placement", {
  genome <- clean_genome(seed = 9L)
  s <- as.character(genome[["Chr1"]])
  ref <- substr(s, 100L, 100L)
  comp <- as.character(Biostrings::complement(Biostrings::DNAString(ref))) # A<->T, C<->G
  vt <- tibble::tibble(chrom = "Chr1", pos = 100L, ref = ref, alt = comp)
  sc <- screen_snp_targets(vt, genome)
  expect_equal(sc$class, "complementary")
  expect_equal(sc$score, "0/0")
  expect_equal(sc$verdict, "doomed")
  expect_equal(sc$reasons, "not_discriminable")
})

test_that("a multi-allelic site is promising when one degenerate primer covers it", {
  genome <- clean_genome(seed = 9L)
  s <- as.character(genome[["Chr1"]])
  # a REF A site whose first ALT (T) is complementary, but a second ALT (G) lets
  # a degenerate 3'-terminal base (K = G/T) exclude REF - the Ensembl allele shape
  i <- which(strsplit(s, "")[[1L]] == "A")[1L]
  vt <- tibble::tibble(chrom = "Chr1", pos = i, ref_genome = list(c("A", "T", "G")))
  sc <- screen_snp_targets(vt, genome)
  expect_equal(sc$class, "degenerate")
  expect_gt(sc$n_ok, 0L)
  expect_equal(sc$verdict, "promising")
})
