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
  # two ALT haplotypes, one ALT primer each, on the shared backbone
  expect_setequal(x$alt, c("H1", "H2"))
  expect_length(x$alt_primers, 2L)
  expect_true(any(!is.na(x$alt_primers)))
  expect_output(print(x), "kasp_lax")
  # one shared REF primer and one shared COM - the haplotypes vary only the ALT
  expect_equal(sum(names(x$pool) == "AS_ref"), 1L)
  expect_equal(sum(names(x$pool) == "COM"), 1L)
  expect_lte(length(x$pool), 4L) # REF + COM + up to two ALT primers
  # the joint pool of every primer is checked for dimers/hairpins
  expect_true(is.data.frame(x$compatibility))
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(x$compatibility)))
  # each designed ALT primer carries *its own* linked variant, so they differ
  if (sum(!is.na(x$alt_primers)) == 2L) {
    expect_false(identical(x$alt_primers[["H1"]], x$alt_primers[["H2"]]))
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
