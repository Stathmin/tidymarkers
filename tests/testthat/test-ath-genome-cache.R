# =============================================================================
#  ath_genome() cannot be poisoned by a truncated cache
# =============================================================================
#  A half-written FASTA makes `file_exists()` stay TRUE forever, so the old
#  code failed on every later call with no way out but a manual delete.  The
#  refetch-once branch is the fix.  Both tests mock the download: nothing here
#  touches the network or the real ~35 MB file.

tair10_fake <- function(path, chromosomes = c("1", "2", "3", "4", "5", "Mt", "Pt")) {
  lines <- unlist(lapply(
    chromosomes,
    function(ch) c(paste0(">", ch, " dna:chromosome"), "ACGTACGTACGTACGT")
  ))
  con <- gzfile(path, "w")
  on.exit(close(con))
  writeLines(lines, con)
  path
}

tair10_cache_path <- function(cache) {
  fs::path(cache, "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
}

test_that("a truncated cache is deleted, refetched once and the genome returned", {
  cache <- withr::local_tempdir()
  writeBin(charToRaw("not a gzip stream at all"), tair10_cache_path(cache))

  calls <- 0L
  local_mocked_bindings(
    tair10_download = function(file, release) {
      calls <<- calls + 1L
      tair10_fake(file)
    },
    .package = "tidymarkers"
  )

  g <- ath_genome(cache = cache)
  expect_equal(calls, 1L)
  expect_s4_class(g, "DNAStringSet")
  expect_setequal(names(g), c(paste0("Chr", 1:5), "ChrM", "ChrC"))
})

test_that("a readable cache is used and the download is never called", {
  cache <- withr::local_tempdir()
  tair10_fake(tair10_cache_path(cache))

  local_mocked_bindings(
    tair10_download = function(...) stop("cache hit must not download"),
    .package = "tidymarkers"
  )

  g <- ath_genome(cache = cache)
  expect_equal(length(g), 7L)
})

test_that("a cache with a base outside the alphabet is refused, not refetched", {
  cache <- withr::local_tempdir()
  con <- gzfile(tair10_cache_path(cache), "w")
  writeLines(c(">1 dna:chromosome", "ACGTXCGT"), con)
  close(con)

  local_mocked_bindings(
    tair10_download = function(...) stop("a bad base must not trigger a download"),
    .package = "tidymarkers"
  )

  expect_error(
    ath_genome(cache = cache),
    "outside the DNA/IUPAC alphabet"
  )
})
