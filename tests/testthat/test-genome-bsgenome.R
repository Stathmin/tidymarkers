# A bare BSgenome gets a real FASTA on disk, written once and then reused (№10).

test_that("a BSgenome genome is materialised to one FASTA, reused on the next call", {
  skip_if_not_installed("BSgenome.Athaliana.TAIR.TAIR9")
  skip_if_not_installed("BSgenome")
  skip_if_not_installed("Biostrings")
  bs <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9

  path <- expect_no_error(genome_fasta_path(bs, "marker_blast_index"))
  expect_true(file.exists(path))
  expect_match(fs::path_file(path), "^genome-[0-9a-f]{12}\\.fa$")

  # the file is a readable FASTA of the whole genome, records and widths included
  Rsamtools::indexFa(path)
  fa <- Rsamtools::FaFile(path)
  si <- GenomeInfoDb::seqinfo(fa)
  bs_si <- GenomeInfoDb::seqinfo(bs)
  expect_equal(length(si), length(bs_si))
  expect_setequal(GenomeInfoDb::seqnames(si), GenomeInfoDb::seqnames(bs_si))
  expect_equal(
    unname(GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqnames(bs_si)]),
    unname(GenomeInfoDb::seqlengths(bs_si))
  )

  # a second call is the same file, untouched: no rewrite
  mtime <- file.mtime(path)
  Sys.setFileTime(path, mtime - 3600)
  expect_identical(genome_fasta_path(bs, "marker_bowtie_index"), path)
  expect_equal(file.mtime(path), mtime - 3600)
})

test_that("an index builder accepts a BSgenome and builds from the materialised FASTA", {
  skip_if_not_installed("BSgenome.Athaliana.TAIR.TAIR9")
  skip_if_not_installed("BSgenome")
  skip_if_not_installed("rBLAST")
  skip_if_not(has_blast(), "NCBI BLAST+ not installed")
  bs <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9
  db <- suppressMessages(marker_blast_index(bs))
  expect_true(fs::file_exists(paste0(db, ".nsq")))
  expect_identical(db, genome_fasta_path(bs, "marker_blast_index"))
})
