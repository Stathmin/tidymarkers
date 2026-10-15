# The pairwise primer-pool compatibility report and its tube-level wrapper.

test_that("primer_set_compatibility reports one row per unordered pair", {
  cc <- primer_set_compatibility(c(AS1 = "GAAGGTGACCAAGTTCATGCT", AS2 = "GAAGGTCGGAGTCAACGGATT"))
  expect_true(all(c("a", "b", "kind", "run", "dg", "pass") %in% names(cc)))
  expect_equal(sum(cc$kind == "self"), 2L)
  expect_equal(sum(cc$kind == "cross"), 1L)
  # a primer and its own reverse complement form a full-length cross-dimer
  p <- "ACGTACGTACGTACGTACGT"
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(p)))
  bad <- primer_set_compatibility(c(x = p, y = rc))
  expect_false(bad$pass[bad$kind == "cross"])
})

test_that("tube_compatibility pools oligos, drops duplicates and tolerates modified bases", {
  expect_null(tube_compatibility(character(0)))
  expect_null(tube_compatibility(c(only = "ACGTACGTACGT")))
  expect_null(tube_compatibility(c(a = "ACGTACGTACGT", b = "ACGTACGTACGT"))) # same oligo
  cc <- tube_compatibility(c(
    fwd = "GAAGGTGACCAAGTTCATGCT", rev = "GAAGGTCGGAGTCAACGGATT",
    probe = "FAM-ACGTXACGT" # a modified probe (non-ACGT marks) must not error
  ))
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(cc)))
  expect_equal(nrow(cc), 6L) # 3 self + 3 cross
})
