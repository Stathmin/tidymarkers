# Indexed genome-wide scan (BLAST) and the dispatcher.

blast_amp <- "ACGTACGTACGTACGTACGT"

blast_fa <- function(fa) {
  writeLines(c(">chr1", paste0(
    strrep("A", 30L), blast_amp, strrep("T", 30L),
    strrep("A", 30L), blast_amp, strrep("T", 30L)
  )), fa)
  fa
}

test_that("marker_blast_index builds a database and specificity_blast counts hits", {
  skip_if_not_full()
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  skip_if_not_installed("rBLAST")
  db <- marker_blast_index(marker_genome(blast_fa(withr::local_tempfile(fileext = ".fa"))))
  res <- specificity_blast(c(P = blast_amp, Q = "GATTACAGATTACAGATTAC"), db, k = NA_integer_)
  expect_equal(res$hits[res$primer == "P"], 4L) # 2 copies on both strands
  expect_equal(res$hits[res$primer == "Q"], 0L)
  expect_equal(res$engine, rep("BLAST+", nrow(res)))
})

test_that("marker_blast_index returns the FASTA prefix and is idempotent", {
  skip_if_not_full()
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  skip_if_not_installed("rBLAST")
  g <- marker_genome(blast_fa(withr::local_tempfile(fileext = ".fa")))
  db <- marker_blast_index(g)
  # the return is the FASTA path (= the BLAST db prefix), not a bare .n* path
  expect_identical(db, g$path)
  idx <- paste0(db, ".nsq")
  expect_true(file.exists(idx))
  before <- file.info(idx)$mtime
  Sys.sleep(1.1) # mtime granularity
  expect_identical(marker_blast_index(g), db)
  expect_identical(file.info(idx)$mtime, before) # an existing index is not rebuilt
})

test_that("specificity_of routes an indexed genome to BLAST, and needs a db", {
  skip_if_not_full()
  withr::local_envvar(LANGUAGE = "en")
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  skip_if_not_installed("rBLAST")
  fa <- blast_fa(withr::local_tempfile(fileext = ".fa"))
  gs <- marker_genome(fa, blast_db = marker_blast_index(marker_genome(fa)))
  res <- specificity_of(c(P = blast_amp), gs)
  expect_gte(res$hits[res$primer == "P"], 1L)
  expect_equal(attr(res, "scan"), "3'-seed")
  expect_equal(res$engine, "BLAST+")
  # an indexed genome without a database is a clear error, not a silent zero
  expect_error(specificity_of(c(P = blast_amp), marker_genome(fa)), "BLAST")
})

test_that("specificity_of uses the in-memory path for a DNAStringSet", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  res <- specificity_of(c(P = "TTTTTTAAACCCGGG"), genome, k = 9L, max_mismatch = 0L)
  expect_equal(res$hits, 4L)
  # the counts answer the seed question; the bowtie path answers the whole-primer
  expect_equal(attr(res, "scan"), "3'-seed")
  expect_equal(res$engine, "Biostrings")
})

test_that("an assembled assay tags both genome scans (D3)", {
  skip_if_not_full()
  skip_if_no_oligoarrayaux()
  skip_if_not_installed("Rbowtie")
  f <- snp_fixture()
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(f$genome, fa)
  genome <- marker_genome(fa)
  genome$bowtie_index <- marker_bowtie_index(genome,
    basename = fs::path(withr::local_tempdir(), "bw")
  )
  out <- suppressWarnings(
    assemble_kasp_assay(genome, f$snp, flank = 150L, cores = 1L, run_blast = TRUE)
  )
  expect_false(inherits(out, "marker_doomed"))
  # the two tables answer different questions: `blast` scans the AS primers with
  # their deliberate mismatch, `specificity` the template-matched allele
  expect_equal(attr(out$blast, "scan"), "whole-primer")
  expect_equal(attr(out$specificity, "scan"), "whole-primer")
  # COM carries no mismatch, so both scans must agree on it
  expect_equal(
    out$blast$hits[out$blast$primer == "COM"],
    out$specificity$hits[out$specificity$primer == "COM"]
  )
})

test_that("specificity_seed asks the indexed genome for the 3'-seed (D1)", {
  fa <- blast_fa(withr::local_tempfile(fileext = ".fa"))
  # no index at all: name the missing index, do not pretend the seed was counted
  expect_error(
    specificity_seed(c(P = blast_amp), marker_genome(fa)),
    "neither a bowtie index nor a BLAST database"
  )
  skip_if_not_full()
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  skip_if_not_installed("rBLAST")
  gs <- marker_genome(fa, blast_db = marker_blast_index(marker_genome(fa)))
  res <- specificity_seed(c(P = blast_amp, Q = "GATTACAGATTACAGATTAC"), gs, k = 15L)
  # the query is the 3'-terminal 15-mer, not the whole primer
  expect_equal(res$seed[res$primer == "P"], substr(blast_amp, 6L, 20L))
  expect_gte(res$hits[res$primer == "P"], 1L)
  expect_equal(res$hits[res$primer == "Q"], 0L)
  # a delegating branch inherits the delegate's engine tag, never its own
  expect_equal(res$engine, rep("BLAST+", nrow(res)))
})

test_that("marker_blast_index keys an in-memory genome by content (D4)", {
  skip_if_not_full()
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  skip_if_not_installed("rBLAST")
  g1 <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGT"))
  g2 <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGA"))
  db1 <- marker_blast_index(g1)
  db2 <- marker_blast_index(g2)
  expect_match(basename(db1), "^genome-[0-9a-f]{12}\\.fa$")
  expect_false(identical(db1, db2)) # two genomes never share a database
  expect_false(grepl("tair10", db1))
  # the alias shares the build, and a second call reuses the database
  idx <- paste0(db1, ".nsq")
  before <- file.info(idx)$mtime
  Sys.sleep(1.1) # mtime granularity
  expect_identical(build_blast_db(g1), db1)
  expect_identical(file.info(idx)$mtime, before)
})

test_that("specificity_of answers a named genome set one genome at a time", {
  genomes <- list(
    A = Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT")),
    B = Biostrings::DNAStringSet(c(chr1 = strrep("G", 24L)))
  )
  primers <- c(P = "TTTTTTAAACCCGGG", Q = "GGGGGGGGGGGGGGGG")
  res <- specificity_of(primers, genomes, k = 9L, max_mismatch = 0L)
  expect_equal(nrow(res), 2L * length(primers))
  expect_equal(names(res)[1], "genome")
  expect_equal(res$genome, c("A", "A", "B", "B"))
  expect_equal(attr(res, "scan"), c(A = "3'-seed", B = "3'-seed"))
  # A carries P twice on both strands and no Q run; B is the mirror image
  hits <- stats::setNames(res$hits, paste(res$genome, res$primer))
  expect_equal(unname(hits[c("A P", "B P", "A Q", "B Q")]), c(4L, 0L, 0L, 16L))
})

test_that("a genome set must be named", {
  genomes <- list(
    Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGT")),
    Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGA"))
  )
  expect_error(specificity_of(c(AS = "ACGTACGT"), genomes), "Name every genome")
  names(genomes) <- c("A", "")
  expect_error(specificity_of(c(AS = "ACGTACGT"), genomes), "Name every genome")
})

test_that("a single DNAStringSet keeps the one-row-per-primer shape", {
  genome <- Biostrings::DNAStringSet(c(chr1 = "AAACCCGGGTTTAAACCCGGGTTT"))
  res <- specificity_of(c(P = "TTTTTTAAACCCGGG"), genome, k = 9L, max_mismatch = 0L)
  expect_equal(nrow(res), 1L)
  expect_false("genome" %in% names(res))
  expect_equal(attr(res, "scan"), "3'-seed")
})
