# FM-index genome scan (Rsubread): an index-based alternative to BLAST.

test_that("kasp_fm_index + specificity_fm count genome loci", {
  skip_if_not_installed("Rsubread")
  skip_if_not_installed("Rsamtools")
  withr::local_seed(1)
  motif <- paste(sample(c("A", "C", "G", "T"), 30L, TRUE), collapse = "")
  spacer <- paste(sample(c("A", "C", "G", "T"), 400L, TRUE), collapse = "")
  uniq <- paste(sample(c("A", "C", "G", "T"), 400L, TRUE), collapse = "")
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(spacer, motif, uniq, motif, uniq)))
  idx <- kasp_fm_index(genome, basename = fs::path(fs::path_temp(), "kasp_fm_test"))
  res <- specificity_fm(c(TWO = motif, ONE = substr(spacer, 1L, 30L)), idx)
  expect_equal(res$hits[res$primer == "TWO"], 2L)
  expect_equal(res$hits[res$primer == "ONE"], 1L)
  expect_true(all(c("primer", "seed", "hits") %in% names(res)))
})

test_that("a kasp_genome carries an attached FM index", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">1", "ACGTACGTACGTACGTACGT"), fa)
  gs <- kasp_genome(fa, fm_index = "/tmp/some/index")
  expect_equal(gs$fm_index, "/tmp/some/index")
  expect_s3_class(gs, "kasp_genome")
})
