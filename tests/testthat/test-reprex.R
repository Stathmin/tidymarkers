# Emits the script to a temp file and sources it in a fresh environment.
reprex_env <- function(path) {
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  env
}

# `marker_reprex()` is not exported until roxygen runs; `:::` keeps lintr quiet.
reprex_write <- function(...) {
  withCallingHandlers(
    tidymarkers:::marker_reprex(...),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

reprex_vcf_records <- function(env) {
  strsplit(env$vcf_rows, "\t", fixed = TRUE)
}

reprex_gc <- function(x) {
  vapply(strsplit(x, "", fixed = TRUE), function(b) sum(b == "G" | b == "C"), integer(1))
}

reprex_kmers <- function(x, k = 20L) {
  n <- nchar(x)
  vapply(seq_len(n - k + 1L), function(i) substring(x, i, i + k - 1L), character(1))
}

test_that("a variant tibble keeps its gaps and chromosome order, and drops columns", {
  v <- tibble::tibble(
    chrom = c("Chr4", "Chr4", "Chr5", "Chr4"),
    pos = c(1710L, 2340L, 900L, 2500L),
    ref = c("A", "C", "G", "T"),
    alt = c("G", "T", "A", "C"),
    sample = c("S1", "S2", "S3", "S4"),
    gt = c("0/1", "1/1", "0/0", "0/1")
  )
  path <- file.path(withr::local_tempdir(), "reprex.R")
  expect_warning(
    tidymarkers:::marker_reprex(v, path, preserve = c("variant_type", "linkage", "n_variants")),
    "dropped"
  )
  expect_false(any(grepl("sample|gt|S1", readLines(path)[-1L])))

  env <- reprex_env(path)
  rec <- reprex_vcf_records(env)
  chrom <- vapply(rec, `[[`, character(1), 1L)
  pos <- as.integer(vapply(rec, `[[`, character(1), 2L))
  # Chr4 appears first and is chr1; Chr5 is chr2, in order of first appearance
  expect_identical(chrom, c("chr1", "chr1", "chr2", "chr1"))
  expect_false(any(pos %in% v$pos))
  # within Chr4: gaps 2340-1710 = 630 and 2500-2340 = 160 survive the shift
  expect_identical(pos, c(1L, 631L, 1L, 791L))
  expect_identical(diff(pos[c(1L, 2L, 4L)]), c(630L, 160L))
  # REF/ALT are replayed verbatim, so the variant type is unchanged
  expect_identical(vapply(rec, `[[`, character(1), 4L), v$ref)
  expect_identical(vapply(rec, `[[`, character(1), 5L), v$alt)
  # no sample column reaches the script
  expect_identical(lengths(rec), rep(8L, 4L))
})

test_that("a DNAStringSet becomes a synthetic sequence of the same length and GC", {
  x <- "ACGTACGTTGCAAGGCTTAACGTACGTTGCAAGGCTTAACGTACGTTGCAAGGCTTAA"
  genome <- Biostrings::DNAStringSet(c(chr1 = x))
  path <- file.path(withr::local_tempdir(), "reprex.R")
  reprex_write(genome, path)

  env <- reprex_env(path)
  expect_identical(nchar(env$syn), nchar(x))
  expect_identical(reprex_gc(env$syn), reprex_gc(x))
  # preserved to the base, not to the sequence
  expect_false(identical(unname(env$syn), x))
  expect_identical(env$nm, "s1")
  # no 20-mer of the original survives anywhere in the script
  expect_false(any(reprex_kmers(x) %in% reprex_kmers(env$syn)))
  expect_false(any(reprex_kmers(x) %in% readLines(path)))
})

test_that("a repeat_structure request keeps the dinucleotide class", {
  x <- "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
  path <- file.path(withr::local_tempdir(), "reprex.R")
  reprex_write(Biostrings::DNAStringSet(c(s1 = x)), path, preserve = "repeat_structure")
  env <- reprex_env(path)
  expect_identical(nchar(env$syn), nchar(x))
  expect_identical(reprex_gc(env$syn), reprex_gc(x))
})

test_that("an unsupported preserve value aborts and names the allowed set", {
  x <- Biostrings::DNAStringSet(c(s1 = "ACGTACGTACGTACGT"))
  err <- expect_error(reprex_write(x, preserve = "sample_id"), "does not implement")
  allowed <- c("gc_content", "variant_type", "n_variants", "linkage", "repeat_structure")
  for (p in allowed) {
    expect_match(conditionMessage(err), p, fixed = TRUE)
  }
})

test_that("the written script parses and carries the header", {
  x <- Biostrings::DNAStringSet(c(s1 = "ACGTACGTACGTACGT"))
  path <- file.path(withr::local_tempdir(), "reprex.R")
  reprex_write(x, path)
  expect_silent(parse(path))
  expect_identical(
    readLines(path)[[1L]],
    paste0("##synthesized_by=tidymarkers_reprex/", utils::packageVersion("tidymarkers"))
  )
  expect_false(any(grepl("sessionInfo", readLines(path))))
})

test_that("the path is returned invisibly, exists and is classed", {
  x <- Biostrings::DNAStringSet(c(s1 = "ACGTACGTACGTACGT"))
  dir <- withr::local_tempdir()
  out <- file.path(dir, "named.R")
  res <- withVisible(reprex_write(x, path = out))
  expect_false(res$visible)
  expect_s3_class(res$value, "marker_reprex")
  expect_identical(unclass(res$value), out)
  expect_true(file.exists(out))

  # a NULL path still writes a tempfile and returns it
  tmp <- withVisible(reprex_write(x))
  expect_false(tmp$visible)
  expect_true(file.exists(unclass(tmp$value)))
})

test_that("an unsupported input aborts", {
  txt <- file.path(withr::local_tempdir(), "x.txt")
  writeLines("x", txt)
  expect_error(reprex_write(1:3), "not a variant table")
  expect_error(reprex_write(list(a = 1)), "not a variant table")
  expect_error(reprex_write(txt), "not a VCF, a BED or a FASTA")
})

test_that("a FASTA path is read by hand and anonymised", {
  fasta <- file.path(withr::local_tempdir(), "g.fa")
  writeLines(c(">chrA", "ACGTACGTTGCAAGGCTTAA"), fasta)
  path <- file.path(withr::local_tempdir(), "reprex.R")
  reprex_write(fasta, path)
  env <- reprex_env(path)
  expect_identical(nchar(env$syn), 20L)
  expect_identical(env$nm, "s1")
})

test_that("a CRLF FASTA is read without the carriage return", {
  fasta <- file.path(withr::local_tempdir(), "g.fa")
  writeBin(charToRaw(">chrA\r\nACGTACGTTGCAAGGCTTAA\r\n"), fasta)
  path <- file.path(withr::local_tempdir(), "reprex.R")
  reprex_write(fasta, path)
  env <- reprex_env(path)
  expect_identical(nchar(env$syn), 20L)
  expect_identical(env$nm, "s1")
})

test_that("a marker_genome is anonymised through its FASTA path", {
  fasta <- file.path(withr::local_tempdir(), "g.fa")
  writeLines(c(">chrA", "ACGTACGTTGCAAGGCTTAACGTACGTTGCAAGGCTTAACGTACGTACGT"), fasta)
  genome <- structure(list(path = fasta), class = "marker_genome")
  path <- file.path(withr::local_tempdir(), "reprex.R")
  reprex_write(genome, path)
  env <- reprex_env(path)
  expect_identical(nchar(env$syn), 50L)
  expect_identical(env$nm, "s1")
  # a BSgenome-backed genome has nothing to read
  bs <- structure(list(path = NULL), class = "marker_genome")
  expect_error(reprex_write(bs), "no file to read")
})

test_that("marker_reprex() refuses to write through a link", {
  dir <- withr::local_tempdir()
  outside <- file.path(dir, "outside.R")
  writeLines("keep", outside)
  link <- file.path(dir, "reprex.R")
  skip_if(!file.symlink(outside, link), "symlinks are not supported")
  expect_error(
    reprex_write("ACGTACGTACGTACGTACGT", link),
    "which is a link"
  )
  expect_identical(readLines(outside), "keep")
})
