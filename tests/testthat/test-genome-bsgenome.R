# A bare BSgenome gets a real FASTA on disk, written once and then reused (№10).

test_that("a BSgenome genome is materialised to one FASTA, reused on the next call", {
  skip_if_not_installed("BSgenome.Athaliana.TAIR.TAIR9")
  skip_if_not_installed("BSgenome")
  skip_if_not_installed("Biostrings")
  bs <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9

  path <- expect_no_error(genome_fasta_path(bs, "marker_blast_index"))
  expect_true(file.exists(path))
  expect_match(fs::path_file(path), "^genome-[0-9a-f]{12}\\.fa$")

  # a second call is the same file, untouched: no rewrite.  This runs before the
  # `FaFile` handle below is opened: that handle stays open for the session, and
  # on Windows it reverts the backwards `Sys.setFileTime`, so the assertion would
  # read our own artefact rather than a rewrite.
  mtime <- file.mtime(path)
  Sys.setFileTime(path, mtime - 3600)
  expect_identical(genome_fasta_path(bs, "marker_bowtie_index"), path)
  expect_equal(file.mtime(path), mtime - 3600)

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

test_that("a same-name, same-length BSgenome does not reuse another genome's FASTA", {
  skip_if_not_installed("BSgenome.Athaliana.TAIR.TAIR9")
  skip_if_not_installed("BSgenome")
  skip_if_not_installed("Biostrings")
  bs <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9
  # A rebuilt or patched package looks exactly like this: same class, package
  # name and seqinfo, one different base.  `getSeq()` reads a per-sequence cache
  # slot before the package rda, so a fresh slot overrides one chromosome.
  skip_if_not(
    ".seqs_cache" %in% slotNames(bs),
    "BSgenome build has no .seqs_cache slot to override"
  )
  bs2 <- bs
  bs2@.seqs_cache <- new.env(parent = emptyenv())
  bs2@.link_counts <- new.env(parent = emptyenv())
  chrc_name <- "ChrC"
  chrc <- BSgenome::getSeq(bs, chrc_name)
  pos <- 100L
  Biostrings::subseq(chrc, pos, pos) <- Biostrings::DNAString("A")
  assign(chrc_name, chrc, envir = bs2@.seqs_cache)
  # a direct `assign()` records no link, so BSgenome's finalizer would evict the
  # override on the next GC; pin the link count to keep the mock stable
  assign(chrc_name, 1e6L, envir = bs2@.link_counts)
  expect_equal(GenomeInfoDb::seqinfo(bs), GenomeInfoDb::seqinfo(bs2))

  # names and lengths cannot tell them apart, so the FASTA must be keyed by
  # content; `genome_signature()` never sees a `BSgenome` for exactly that reason
  path <- genome_fasta_path(bs, "marker_blast_index")
  path2 <- genome_fasta_path(bs2, "marker_blast_index")
  expect_false(identical(path, path2))

  Rsamtools::indexFa(path2)
  gr <- GenomicRanges::GRanges("ChrC", IRanges::IRanges(pos, pos))
  fa <- Rsamtools::scanFa(Rsamtools::FaFile(path2), gr)[[1L]]
  expect_identical(as.character(fa), "A")
})
