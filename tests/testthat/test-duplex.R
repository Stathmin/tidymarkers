# The duplex engine: two loci in one tube + the compatibility stage.

test_that("duplex_compatibility flags inter-locus complementarity and checks the budget", {
  mk <- function(as_ref, as_alt, com, amplicon) {
    list(assay = tibble::tibble(as_ref = as_ref, as_alt = as_alt, com = com, amplicon = amplicon))
  }
  a1 <- mk("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGA", "TTTTGCGCGCGCGCGCGCAA", 90L)
  # locus 2's `com` is the reverse complement of locus 1's -> a long complementary run
  a2 <- mk("GGGGCCCCGGGGCCCCGGGG", "GGGGCCCCGGGGCCCCGGGA", "TTGCGCGCGCGCGCGCAAAA", 95L)
  cc <- duplex_compatibility(list(a1, a2), readout = "co-dominant", homologous = TRUE)
  expect_setequal(cc$gate, c(
    "gate_duplex_targets", "gate_duplex_tm", "gate_duplex_complement",
    "gate_duplex_selfcomp", "gate_duplex_amplicon"
  ))
  expect_true(cc$pass[cc$gate == "gate_duplex_targets"]) # co-dominant needs homologous loci
  expect_false(cc$pass[cc$gate == "gate_duplex_complement"]) # the revcomp pair
  expect_true(cc$pass[cc$gate == "gate_duplex_amplicon"])
  # a co-dominant call on non-homologous loci is a targets-gate failure
  bad <- duplex_compatibility(list(a1, a2), readout = "co-dominant", homologous = FALSE)
  expect_false(bad$pass[bad$gate == "gate_duplex_targets"])
  # presence/absence accepts non-homologous loci
  pa <- duplex_compatibility(list(a1, a2), readout = "presence-absence", homologous = FALSE)
  expect_true(pa$pass[pa$gate == "gate_duplex_targets"])
})

test_that("design_duplex assembles two loci and reports a verdict", {
  skip_if_no_oligoarrayaux()
  W <- bseq(800L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  bases <- c(substr(W, 250L, 250L), substr(W, 550L, 550L))
  loci <- tibble::tibble(
    chrom = "Chr1", pos = c(250L, 550L), ref = bases,
    alt = vapply(bases, function(b) setdiff(c("A", "C", "G", "T"), b)[1L], character(1))
  )
  d <- suppressWarnings(design_duplex(genome, loci, flank = 150L, cores = 1L))
  expect_s3_class(d, "kasp_duplex")
  expect_true(d$verdict %in% c("pass", "warn", "fail"))
  expect_length(d$assays, 2L)
  expect_output(print(d), "kasp_duplex")
})

test_that("design_duplex rejects anything other than two loci", {
  expect_error(design_duplex(NULL, tibble::tibble(chrom = "1", pos = 1L, ref = "A", alt = "G")))
})
