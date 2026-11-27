# ADR 0005: the genome declares its strategy and the scan names its question.

test_that("genome_index declares the capability matrix per strategy", {
  gs <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGT"))
  expected <- list(
    memory = list(
      seed_scan = TRUE, whole_primer_scan = FALSE,
      mismatch_tolerant = TRUE, max_mismatch = Inf, engine = "Biostrings"
    ),
    bowtie = list(
      seed_scan = TRUE, whole_primer_scan = TRUE,
      mismatch_tolerant = TRUE, max_mismatch = 3, engine = "Rbowtie"
    ),
    blast = list(
      seed_scan = TRUE, whole_primer_scan = TRUE,
      mismatch_tolerant = FALSE, max_mismatch = 0, engine = "BLAST+"
    ),
    fm = list(
      seed_scan = FALSE, whole_primer_scan = TRUE,
      mismatch_tolerant = TRUE, max_mismatch = 2, engine = "Rsubread"
    )
  )
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)
  for (strategy in names(expected)) {
    genome <- if (strategy == "memory") gs else marker_genome(fa)
    idx <- genome_index(genome, strategy)
    expect_identical(idx$strategy, strategy)
    expect_identical(capabilities(idx), expected[[strategy]])
    expect_identical(capabilities(idx, "engine"), expected[[strategy]]$engine)
  }
  expect_error(genome_index("not a genome", "memory"), "DNAStringSet")
})

test_that("memory and bowtie index give the same verdict on a small fixture", {
  skip_if_not_installed("Rbowtie")
  withr::local_seed(7)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  unique <- rand(25L)
  repeated <- rand(20L)
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
    rand(300L), unique, rand(300L), repeated, rand(300L), repeated, rand(300L)
  )))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(genome, fa)
  primers <- c(UNIQUE = unique, TWICE = repeated)

  mem <- specificity_of(primers, genome_index(genome, "memory"),
    k = 15L, max_mismatch = 0L
  )
  idx <- specificity_of(
    primers,
    genome_index(marker_genome(fa), "bowtie", dir = withr::local_tempdir()),
    k = 15L, max_mismatch = 0L
  )

  # the two paths answer different questions, but the locus verdict is the same
  verdict <- function(res) ifelse(res$hits == 1L, "unique", "off-target")
  expect_identical(mem$primer, idx$primer)
  expect_identical(mem$hits, idx$hits)
  expect_identical(verdict(mem), verdict(idx))
  expect_identical(verdict(mem), c("unique", "off-target"))

  # every path declares attr(x, "scan"), and the strategy names it
  expect_identical(attr(mem, "scan"), "3'-seed")
  expect_identical(attr(idx, "scan"), "whole-primer")
  expect_identical(mem$engine, rep("Biostrings", 2L))
  expect_identical(idx$engine, rep("Rbowtie", 2L))
})

test_that("the manifest names every shard and hashes its artefact set", {
  skip_if_not_installed("Rbowtie")
  withr::local_seed(11)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  gs <- Biostrings::DNAStringSet(c(
    chr1 = rand(400L), chr2 = rand(350L), chr3 = rand(300L)
  ))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)
  dir <- withr::local_tempdir()

  idx <- genome_index(marker_genome(fa), "bowtie", chroms = c("chr1", "chr2"), dir = dir)
  expect_identical(names(idx$manifest), INDEX_MANIFEST_COLUMNS)
  expect_identical(idx$manifest$chrom, c("chr1", "chr2"))
  expect_identical(idx$manifest$tool, rep("Rbowtie", 2L))
  expect_false(is.na(idx$manifest$tool_version[[1L]]))
  expect_true(all(is.na(idx$manifest$sha256))) # nothing built by the constructor

  # the object holds paths, never sequence
  expect_false(any(vapply(idx, function(x) inherits(x, "DNAStringSet"), logical(1))))
  # a lazy per-chromosome build keeps no whole chromosome in the query cache
  genome <- marker_genome(fa)
  idx <- genome_index(genome, "bowtie", chroms = c("chr1", "chr2"), dir = dir)
  primer <- as.character(Biostrings::subseq(gs[["chr1"]], 1L, 20L))
  specificity_of(c(AS = primer), idx, k = 15L, max_mismatch = 0L)
  expect_identical(length(ls(genome$cache)), 0L)

  man <- index_manifest(idx)
  expect_true(all(!is.na(man$sha256)))
  # one row per requested chromosome; the measured hash covers the shard files
  files <- sort(Sys.glob(paste0(man$path[[1L]], ".*")))
  parts <- vapply(files, function(f) {
    paste0(basename(f), ":", digest::digest(file = f, algo = "sha256"))
  }, character(1))
  expect_identical(man$sha256[[1L]], digest::digest(paste(parts, collapse = "|"), algo = "sha256"))
  expect_identical(man$bytes[[1L]], sum(as.numeric(fs::file_size(files))))
  expect_false(is.na(man$created[[1L]]))
  expect_identical(man$source_sha256[[1L]], digest::digest(file = fa, algo = "sha256"))
})

test_that("a valid shard is reused without touching its mtime", {
  skip_if_not_installed("Rbowtie")
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  gs <- Biostrings::DNAStringSet(c(chr1 = rand(200L), chr2 = rand(200L)))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)
  dir <- withr::local_tempdir()
  primer <- c(AS = as.character(Biostrings::subseq(gs[["chr1"]], 1L, 20L)))
  build <- function() {
    idx <- genome_index(marker_genome(fa), "bowtie", chroms = c("chr1", "chr2"), dir = dir)
    specificity_of(primer, idx, k = 15L, max_mismatch = 0L)
    idx
  }
  first <- build()
  shard <- paste0(first$manifest$path[[1L]], ".1.ebwt")
  before <- fs::file_info(shard)$modification_time
  build()
  expect_identical(fs::file_info(shard)$modification_time, before)
  # a fresh object over the same directory is the manifest after a restart
  restarted <- genome_index(marker_genome(fa), "bowtie",
    chroms = c("chr1", "chr2"), dir = dir
  )
  expect_identical(restarted$manifest$path, first$manifest$path)
  expect_no_error(index_verify(restarted))
})

test_that("a missing shard names the chromosome, a changed one names the shard", {
  skip_if_not_installed("Rbowtie")
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  gs <- Biostrings::DNAStringSet(c(chr1 = rand(200L), chr2 = rand(200L)))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)
  dir <- withr::local_tempdir()
  primer <- c(AS = as.character(Biostrings::subseq(gs[["chr1"]], 1L, 20L)))
  build <- function() {
    idx <- genome_index(marker_genome(fa), "bowtie", chroms = c("chr1", "chr2"), dir = dir)
    specificity_of(primer, idx, k = 15L, max_mismatch = 0L)
    # record the built shards, as a second run would
    idx$manifest <- index_manifest(idx)
    idx
  }

  tampered_idx <- build()
  cat("x", file = paste0(tampered_idx$manifest$path[[1L]], ".1.ebwt"), append = TRUE)
  expect_error(
    index_verify(tampered_idx),
    class = "tidymarkers_index_integrity_error"
  )

  unlink(Sys.glob(paste0(tampered_idx$manifest$path[2L], ".*")))
  idx <- genome_index(marker_genome(fa), "bowtie", chroms = c("chr1", "chr2"), dir = dir)
  expect_error(
    index_verify(idx),
    regexp = "chromosome chr2",
    class = "tidymarkers_index_missing_error"
  )
})

test_that("scan_question names the question and refuses an absent one", {
  gs <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGTACGTACGTACGT"))
  mem <- specificity_of(c(AS = "ACGTACGTACGTACGT"), gs, k = 12L, max_mismatch = 0L)
  expect_identical(scan_question(mem), "3'-seed")
  expect_s3_class(mem, "specificity_scan")
  expect_error(scan_question(tibble::tibble(hits = 1L)), class = "tidymarkers_scan_error")
  bad <- mem
  attr(bad, "scan") <- "seed-ish"
  expect_error(scan_question(bad), class = "tidymarkers_scan_error")
})

test_that("a capability the strategy lacks is a named error, not a warning", {
  gs <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGT"))
  expect_error(
    specificity_of(c(AS = "ACGTACGTACGTACGT"), genome_index(gs, "bowtie"), k = 12L),
    class = "tidymarkers_capability_error"
  )
  expect_error(
    specificity_of(c(AS = "ACGTACGTACGTACGT"), genome_index(gs, "blast"), max_mismatch = 1L),
    class = "tidymarkers_capability_error"
  )
})

test_that("a table cannot mix two questions for one genome", {
  base <- tibble::tibble(hits = 1L)
  mixed <- base
  attr(mixed, "scan") <- c("3'-seed", "whole-primer")
  expect_error(
    specificity_verdict(mixed, "kasp"),
    class = "tidymarkers_scan_mixed_error"
  )
  # one question per genome is allowed
  named <- empty_pcr_products()
  attr(named, "scan") <- c(A = "3'-seed", B = "whole-primer")
  expect_no_error(specificity_verdict(named, "kasp"))
})

test_that("memory and bowtie agree through specificity_verdict", {
  skip_if_not_installed("Rbowtie")
  withr::local_seed(3)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  site <- rand(24L)
  genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
    rand(200L), site, rand(200L), site, rand(200L)
  )))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(genome, fa)
  hits <- in_silico_pcr(site, site, genome)

  mem <- specificity_verdict(hits, "kasp", genome = genome_index(genome, "memory"))
  idx <- specificity_verdict(
    hits, "kasp",
    genome = genome_index(marker_genome(fa), "bowtie", dir = withr::local_tempdir())
  )
  expect_identical(mem$verdict, idx$verdict)
  expect_identical(mem$n_tolerated, idx$n_tolerated)
  expect_identical(mem$products$violations, idx$products$violations)
  expect_identical(idx$tier, "index")
})

test_that("a pangenome is a named list of per-chromosome shards", {
  skip_if_not_installed("Rbowtie")
  withr::local_seed(21)
  rand <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  site <- rand(24L)
  make_assembly <- function(with_site) {
    filler <- rand(300L)
    Biostrings::DNAStringSet(c(
      chr1 = if (with_site) paste0(filler, site, filler) else paste0(filler, filler),
      chr2 = rand(500L)
    ))
  }
  fa <- vapply(c(A = TRUE, B = FALSE, D = FALSE), function(site_in) {
    p <- tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(make_assembly(site_in), p)
    p
  }, character(1))
  dir <- withr::local_tempdir()
  pan <- lapply(names(fa), function(nm) {
    genome_index(marker_genome(fa[[nm]], build = nm), "bowtie",
      chroms = c("chr1", "chr2"), dir = file.path(dir, nm)
    )
  })
  names(pan) <- names(fa)

  # one shard per (genome_build, chrom), never one monolithic build
  for (nm in names(pan)) {
    idx <- pan[[nm]]
    expect_identical(idx$manifest$chrom, c("chr1", "chr2"))
    expect_identical(nrow(idx$manifest), 2L)
    expect_identical(idx$build, nm)
    expect_true(all(idx$manifest$genome_build == nm))
  }
  # a build is in the path, so two builds cannot share a shard file
  expect_length(unique(unlist(lapply(pan, function(x) x$manifest$path))), 6L)

  res <- specificity_of(c(AS = site), pan, k = 15L, max_mismatch = 0L)
  expect_identical(res$genome, c("A", "B", "D"))
  expect_identical(res$hits, c(1L, 0L, 0L))
  expect_identical(res$build, c("A", "B", "D"))
  expect_identical(scan_question(res), c(A = "whole-primer", B = "whole-primer", D = "whole-primer"))
  # a table spanning builds has no single verdict
  expect_error(
    specificity_verdict(res, "pcr"),
    class = "tidymarkers_mixed_build_error"
  )
})

test_that("a build is explicit and survives read_genome_fasta()", {
  gs <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGT"))
  fa <- withr::local_tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(gs, fa)

  genome <- read_genome_fasta(fa, build = "A")
  expect_identical(S4Vectors::metadata(genome)$genome_build, "A")
  idx <- genome_index(genome, "memory")
  expect_identical(idx$build, "A")
  expect_identical(idx$manifest$genome_build, "A")

  expect_error(marker_genome(fa, build = c("A", "B")), "one short release id")
  expect_error(genome_index(gs, "memory", build = 1), "one short release id")
  expect_error(
    genome_index(list(A = gs), "memory"),
    "one index per build"
  )

  # an unnamed component cannot join a pangenome
  unnamed <- list(
    A = genome_index(read_genome_fasta(fa, build = "A"), "memory"),
    B = genome_index(gs, "memory")
  )
  expect_error(
    specificity_of(c(AS = "ACGTACGTACGT"), unnamed),
    class = "tidymarkers_mixed_build_error"
  )
})

test_that("a build or chromosome name cannot leave the index directory", {
  gs <- Biostrings::DNAStringSet(c(chr1 = "ACGTACGTACGTACGTACGT"))
  expect_error(marker_genome(gs, build = "../escape"), "one short release id")
  expect_error(
    tidymarkers:::index_shards(gs, "bowtie", "a/b", tempdir()),
    "path separator"
  )
  expect_error(
    tidymarkers:::index_shards(gs, "bowtie", "..", tempdir()),
    "path separator"
  )
})

test_that("the TAIR10 gold set gives the same verdict on memory and bowtie", {
  skip_if_not_installed("Rbowtie")
  fa <- testthat::test_path("fixtures", "tair10", "genome.fa")
  expected <- utils::read.csv(
    testthat::test_path("fixtures", "tair10", "expected.csv"),
    stringsAsFactors = FALSE
  )
  primers <- c(
    Chr1_unique = "GGTTATAAATTATTAGAAGTTTTC",
    Chr2_unique = "ATTCGTCGACCAGGACGGCGGAAT",
    Chr3_unique = "CGTCGTATTCGGAAGTTTGCTCGA",
    Chr1_twice = "GTCATTTGGAGATACGAAATCAAA"
  )
  mem <- specificity_of(
    primers, genome_index(Biostrings::readDNAStringSet(fa), "memory"),
    k = SEED_K, max_mismatch = 0L
  )
  idx <- specificity_of(
    primers,
    genome_index(marker_genome(fa), "bowtie",
      chroms = c("Chr1", "Chr2", "Chr3"), dir = withr::local_tempdir()
    ),
    k = SEED_K, max_mismatch = 0L
  )
  verdict <- function(hits) ifelse(hits == 1L, "unique", "off-target")
  expect_identical(mem$primer, expected$primer)
  expect_identical(mem$hits, expected$hits)
  expect_identical(idx$hits, expected$hits)
  expect_identical(verdict(mem$hits), verdict(expected$hits))
  expect_identical(verdict(idx$hits), verdict(expected$hits))
})
