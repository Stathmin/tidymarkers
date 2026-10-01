# The two-level plot vocabulary: feature level (one oligo) vs alignment level.

test_that("check_level splits criteria into feature and alignment", {
  expect_equal(check_level("AS Tm (C)"), "feature")
  expect_equal(check_level("worst hairpin dG37"), "feature")
  expect_equal(check_level("amplicon length (bp)"), "alignment")
  expect_equal(check_level("in-silico PCR (REF 1 / ALT 0)"), "alignment")
})

test_that("clamp_window follows the primer's 3' end on either strand", {
  expect_equal(clamp_window(10L, 19L, forward = TRUE, n = 5L), c(15L, 19L))
  expect_equal(clamp_window(10L, 19L, forward = FALSE, n = 5L), c(10L, 14L))
})

# A tiny real assay on a random genome (skips without OligoArrayAux).
make_plot_assay <- function(seed = 4L, n = 500L, pos = 250L) {
  for (s in seed:(seed + 19L)) {
    withr::with_seed(s, {
      w <- paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
      genome <- Biostrings::DNAStringSet(c(Chr1 = w))
    })
    ref <- as.character(Biostrings::subseq(genome[["Chr1"]], pos, pos))
    snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
    vs <- tibble::tibble(chrom = "Chr1", pos = pos, consequence = "missense_variant")
    a <- try(suppressWarnings(assemble_kasp_assay(genome, snp, variants = vs, flank = 150L)), silent = TRUE)
    if (!inherits(a, "try-error")) {
      return(list(assay = a, snp = snp, variants = vs))
    }
  }
  testthat::skip("no random genome assembled into an assay")
}

test_that("the assay plots render and carry the two-level vocabulary", {
  skip_if_not(nzchar(Sys.which("hybrid-min")), "OligoArrayAux not on PATH")
  m <- make_plot_assay()

  feats <- primer_features(m$assay)
  expect_setequal(feats$role, c("AS-ref", "AS-alt", "COM"))
  expect_true(all(feats$status %in% c("ok", "broken")))
  expect_true(all(feats$clamp %in% c("ok", "broken")))

  # styled letters: clamp bases are bold, motif bases are red, rest neutral
  let <- primer_letters(m$assay)
  expect_setequal(unique(let$face), c("plain", "bold"))
  expect_true("black" %in% let$colour)

  expect_s3_class(plot_checks(m$assay), "ggplot")
  expect_s3_class(plot_primers(m$assay, m$variants, m$snp), "patchwork")
  expect_s3_class(plot_alignment(m$assay), "ggplot")
})
