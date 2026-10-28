# Informativeness of a marker on an available panel: pairwise discrimination,
# PIC and expected heterozygosity, computed from lines (genotypes) or from the
# package's own haplotype table.

informativeness_target <- function() {
  tibble::tibble(chrom = "Chr1", pos = 250L, ref = "A", alt = "G")
}

# a panel where ALT is at frequency 0.5: four ALT homozygotes, four REF
informativeness_panel <- function(genotype = c(2, 2, 2, 2, 0, 0, 0, 0)) {
  tibble::tibble(
    line = paste0("L", seq_along(genotype)),
    chrom = "Chr1", pos = 250L, ref = "A", alt = "G",
    genotype = as.numeric(genotype)
  )
}

test_that("a biallelic panel reports PIC, heterozygosity and pair discrimination", {
  r <- marker_informativeness(
    informativeness_target(),
    genotypes = informativeness_panel()
  )
  expect_s3_class(r, "tbl_df")
  expect_equal(nrow(r), 1L)
  expect_equal(r$n_lines, 8L)
  expect_equal(r$alt_freq, 0.5)
  # Botstein's biallelic PIC at p = q = 0.5 is 0.375; He is 0.5
  expect_equal(r$pic, 0.375)
  expect_equal(r$expected_heterozygosity, 0.5)
  # all cross-genotype pairs differ: 4 * 4 of choose(8, 2)
  expect_equal(r$distinguishable_pairs, 16 / 28)
  expect_equal(r$status, "ok")
  expect_true(is.na(r$reason))
})

test_that("a panel too small to estimate anything returns NA and says why", {
  # the two parents of a cross: they differ, but two lines say nothing about
  # a population
  r <- marker_informativeness(
    informativeness_target(),
    genotypes = informativeness_panel(c(2, 0))
  )
  expect_equal(r$n_lines, 2L)
  expect_equal(r$status, "insufficient_sample")
  expect_true(is.na(r$pic))
  expect_true(is.na(r$expected_heterozygosity))
  expect_match(r$reason, "At least 4")
})

test_that("a marker that separates only the parents is weak on the population", {
  # one ALT homozygote among nine REF: the cross would look informative, the
  # panel does not
  panel <- informativeness_panel(c(2, rep(0, 9)))
  r <- marker_informativeness(informativeness_target(), genotypes = panel)
  expect_equal(r$alt_freq, 0.1)
  expect_lt(r$pic, 0.2)
  expect_equal(r$distinguishable_pairs, 9 / 45)
  expect_equal(r$status, "ok")
})

test_that("a monomorphic panel is reported as such, not as missing", {
  r <- marker_informativeness(
    informativeness_target(),
    genotypes = informativeness_panel(rep(0, 6))
  )
  expect_equal(r$status, "monomorphic")
  expect_equal(r$pic, 0)
  expect_equal(r$expected_heterozygosity, 0)
  expect_equal(r$distinguishable_pairs, 0)
  expect_match(r$reason, "one allele")
})

test_that("missing genotype calls do not count towards the line number", {
  panel <- informativeness_panel(c(2, 2, 0, 0, NA))
  r <- marker_informativeness(informativeness_target(), genotypes = panel)
  expect_equal(r$n_lines, 4L)
  expect_equal(r$alt_freq, 0.5)
})

test_that("haplotypes give the frequency-weighted answer, no sample-size gate", {
  # H1 carries nothing at the target (so it is REF there); only H2 has a row
  h <- tibble::tibble(
    haplotype = c("H1", "H2"),
    freq = c(0.4, 0.6),
    chrom = "Chr1", pos = c(240L, 250L),
    ref = "A", alt = c("T", "G")
  )
  r <- marker_informativeness(informativeness_target(), haplotypes = h)
  expect_equal(r$n_lines, 2L)
  expect_equal(r$alt_freq, 0.6)
  expect_equal(r$expected_heterozygosity, 0.48)
  expect_equal(r$pic, 0.3648)
  expect_equal(r$distinguishable_pairs, 0.48)
  expect_equal(r$status, "ok")
})

test_that("a plan supplies the targets", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">Chr1", paste(rep("A", 400L), collapse = "")), fa)
  genome <- marker_genome(fa)
  tg <- tibble::tibble(
    chrom = c("Chr1", "Chr1"), pos = c(120L, 250L),
    ref = "A", alt = "G"
  )
  panel <- dplyr::bind_rows(
    informativeness_panel(),
    dplyr::mutate(informativeness_panel(), pos = 120L)
  )
  plan <- marker_plan(genome, targets = tg)
  r <- marker_informativeness(plan, genotypes = panel)
  expect_equal(nrow(r), 2L)
  expect_true(all(r$status == "ok"))
  expect_equal(r$pic, rep(0.375, 2L))
})

test_that("exactly one panel is required", {
  expect_error(marker_informativeness(informativeness_target()), "genotypes")
  expect_error(
    marker_informativeness(
      informativeness_target(),
      genotypes = informativeness_panel(),
      haplotypes = tibble::tibble()
    ),
    "one of"
  )
})
