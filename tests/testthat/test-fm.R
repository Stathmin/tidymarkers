# FM-index genome scan (Rsubread): an index-based alternative to BLAST.

test_that("specificity_fm returns finite whole-primer counts on a small index", {
  # one specificity_fm() call costs ~52 s regardless of genome size or primer
  # count, so this stays out of the fast suite (measured on 1.8 kb)
  skip_if_not_full()
  skip_if_not_installed("Rsubread")
  skip_if_not_installed("Rsamtools")
  withr::local_seed(3)
  motif <- paste(sample(c("A", "C", "G", "T"), 28L, TRUE), collapse = "")
  filler <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = "")
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(filler, motif, filler, motif, filler)))
  idx <- marker_fm_index(genome, basename = fs::path(withr::local_tempdir(), "fast"))
  res <- specificity_fm(
    c(TWO = motif, ONE = substr(filler, 1L, 28L), ABSENT = strrep("ACGT", 7L)), idx
  )
  expect_equal(res$primer, c("TWO", "ONE", "ABSENT"))
  expect_true(all(is.finite(res$hits)))
  expect_true(all(res$hits >= 0L))
  expect_equal(res$hits[res$primer == "TWO"], 2L)
  expect_equal(unique(res$engine), "Rsubread")
})

test_that("marker_fm_index + specificity_fm count genome loci", {
  skip_if_not_full()
  skip_if_not_installed("Rsubread")
  skip_if_not_installed("Rsamtools")
  withr::local_seed(1)
  motif <- paste(sample(c("A", "C", "G", "T"), 30L, TRUE), collapse = "")
  spacer <- paste(sample(c("A", "C", "G", "T"), 400L, TRUE), collapse = "")
  uniq <- paste(sample(c("A", "C", "G", "T"), 400L, TRUE), collapse = "")
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(spacer, motif, uniq, motif, uniq)))
  idx <- marker_fm_index(genome, basename = fs::path(fs::path_temp(), "marker_fm_test"))
  res <- specificity_fm(c(TWO = motif, ONE = substr(spacer, 1L, 30L)), idx)
  expect_equal(res$hits[res$primer == "TWO"], 2L)
  expect_equal(res$hits[res$primer == "ONE"], 1L)
  expect_true(all(c("primer", "seed", "hits", "engine") %in% names(res)))
  expect_equal(unique(res$engine), "Rsubread")
  expect_equal(marker_evaluation_stamp(res), marker_engine_stamp())
})

test_that("marker_fm_index says the design pipeline does not query it (D5)", {
  skip_if_not_full()
  skip_if_not_installed("Rsubread")
  skip_if_not_installed("Rsamtools")
  genome <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGTACGTACGTACGT"))
  base <- fs::path(withr::local_tempdir(), "warn")
  expect_message(idx <- marker_fm_index(genome, basename = base), "read mapper")
  # the reuse path says it too: the trap is the index, not the build
  expect_message(marker_fm_index(genome, basename = base), "read mapper")
})

test_that("a marker_genome carries an attached FM index", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">1", "ACGTACGTACGTACGTACGT"), fa)
  gs <- marker_genome(fa, fm_index = "/tmp/some/index")
  expect_equal(gs$fm_index, "/tmp/some/index")
  expect_s3_class(gs, "marker_genome")
})

test_that("unnamed sequences are named before the index is built", {
  skip_if_not_full()
  skip_if_not_installed("Rsubread")
  skip_if_not_installed("Rsamtools")
  withr::local_seed(2)
  motif <- paste(sample(c("A", "C", "G", "T"), 30L, TRUE), collapse = "")
  filler <- paste(sample(c("A", "C", "G", "T"), 400L, TRUE), collapse = "")
  # no names: the BAM header used to come out unusable downstream
  genome <- Biostrings::DNAStringSet(paste0(filler, motif, filler, motif, filler))
  base <- fs::path(withr::local_tempdir(), "unnamed")
  idx <- marker_fm_index(genome, basename = base)
  res <- specificity_fm(c(TWO = motif), idx)
  expect_equal(res$hits, 2L)
})

test_that("specificity_of refuses an FM-only genome and names the bowtie path", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">chr1", "ACGTACGTACGTACGTACGT"), fa)
  genome <- marker_genome(fa, fm_index = "/tmp/some/index")
  # the Subread index is 24x heavier per base than bowtie1 and flat-cost per
  # query, so it is no longer the whole-genome scan; say so, do not dispatch
  expect_error(
    specificity_of(c(AS = "ACGTACGT"), genome),
    "no longer dispatches"
  )
  expect_error(
    specificity_of(c(AS = "ACGTACGT"), genome),
    "marker_bowtie_index"
  )
})
