# lax KASP: one ALT primer per ALT haplotype (linked variants allowed).

lax_fixture <- function() {
  W <- withr::with_seed(2L, paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  v <- 300L
  ref <- substr(W, v, v)
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  bases <- c(substr(W, v - 4L, v - 4L), substr(W, v + 4L, v + 4L))
  other <- function(b) setdiff(c("A", "C", "G", "T"), b)[1L]
  snp <- tibble::tibble(chrom = "Chr1", pos = v, ref = ref, alt = alt)
  # H1 and H2 both carry the ALT at the target but differ in a linked SNP; H3 is REF
  haps <- tibble::tibble(
    haplotype = c("H1", "H1", "H2", "H2", "H3"),
    freq = c(0.3, 0.3, 0.2, 0.2, 0.5),
    chrom = "Chr1",
    pos = c(v, v + 4L, v, v - 4L, v + 4L),
    ref = c(ref, bases[2L], ref, bases[1L], bases[2L]),
    alt = c(alt, other(bases[2L]), alt, other(bases[1L]), other(bases[2L]))
  )
  list(
    genome = genome, snp = snp, haps = haps,
    linked_h1 = v + 4L, linked_h2 = v - 4L
  )
}

test_that("design_lax_marker designs one ALT primer per ALT haplotype", {
  skip_if_no_oligoarrayaux()
  m <- lax_fixture()
  x <- suppressWarnings(design_lax_marker(m$genome, m$snp, m$haps, flank = 150L, cores = 1L))
  expect_s3_class(x, "kasp_lax")
  # two ALT haplotypes, each with its own assembled assay
  expect_setequal(x$alt, c("H1", "H2"))
  expect_length(x$assays, 2L)
  expect_true(all(x$designed %in% x$alt))
  expect_output(print(x), "kasp_lax")
  # the joint pool of every primer is checked for dimers/hairpins, not just each
  # haplotype's own triplet
  expect_true(is.data.frame(x$compatibility))
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(x$compatibility)))
  expect_true(length(x$pool) >= length(x$designed))
  # each assay baked *its own* linked variant into the ALT primer, so the two ALT
  # primers are not the same sequence
  if (length(x$designed) == 2L) {
    p1 <- x$assays[["H1"]]$assay$as_alt[1L]
    p2 <- x$assays[["H2"]]$assay$as_alt[1L]
    expect_false(identical(p1, p2))
  }
})

test_that("design_lax_marker degenerates to pure KASP without a haplotype table", {
  skip_if_no_oligoarrayaux()
  m <- lax_fixture()
  x <- suppressWarnings(design_lax_marker(m$genome, m$snp, flank = 150L, cores = 1L))
  expect_s3_class(x, "kasp_lax")
  expect_equal(x$alt, "ALT")
  expect_length(x$assays, 1L)
})

test_that("the lax scheme is registered and has its own decision stage", {
  expect_true("lax" %in% kasp_schemes()$scheme)
  s <- kasp_decision_schema("lax")
  expect_true("lax" %in% s$nodes$stage)
  # it is a KASP-family scheme, so it also carries the as/pair stages
  expect_true(all(c("as", "pair") %in% s$nodes$stage))
})
