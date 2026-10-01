# Indexed genome-wide scan (BLAST) and the dispatcher.

blast_amp <- "ACGTACGTACGTACGTACGT"

blast_fa <- function() {
  fa <- tempfile(fileext = ".fa")
  writeLines(c(">chr1", paste0(
    strrep("A", 30L), blast_amp, strrep("T", 30L),
    strrep("A", 30L), blast_amp, strrep("T", 30L)
  )), fa)
  fa
}

test_that("kasp_blast_index builds a database and specificity_blast counts hits", {
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  skip_if_not_installed("rBLAST")
  db <- kasp_blast_index(kasp_genome(blast_fa()))
  res <- specificity_blast(c(P = blast_amp, Q = "GATTACAGATTACAGATTAC"), db, k = NA_integer_)
  expect_equal(res$hits[res$primer == "P"], 4L) # 2 copies on both strands
  expect_equal(res$hits[res$primer == "Q"], 0L)
})

test_that("specificity_of routes an indexed genome to BLAST, and needs a db", {
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  skip_if_not_installed("rBLAST")
  fa <- blast_fa()
  gs <- kasp_genome(fa, blast_db = kasp_blast_index(kasp_genome(fa)))
  res <- specificity_of(c(P = blast_amp), gs)
  expect_gte(res$hits[res$primer == "P"], 1L)
  # an indexed genome without a database is a clear error, not a silent zero
  expect_error(specificity_of(c(P = blast_amp), kasp_genome(fa)), "BLAST")
})

test_that("specificity_of uses the in-memory path for a DNAStringSet", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  res <- specificity_of(c(P = "TTTTTTAAACCCGGG"), genome, k = 9L, max_mismatch = 0L)
  expect_equal(res$hits, 4L)
})
