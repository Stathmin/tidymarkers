# Allele sets and the multi-allelic (degenerate-primer) KASP criterion.

test_that("iupac_code inverts the Biostrings alphabet", {
  expect_equal(iupac_code("A"), "A")
  expect_equal(iupac_code(c("C", "G")), "S")
  expect_equal(iupac_code(c("C", "G", "T")), "B")
  expect_equal(iupac_code(c("A", "C", "G", "T")), "N")
  expect_equal(iupac_code(character()), NA_character_)
  expect_equal(iupac_code(c("A", "N")), NA_character_)
})

test_that("kasp_site_class names the designable and the hopeless sites", {
  b <- kasp_site_class("A", "G")
  expect_equal(b$class, "biallelic")
  expect_true(b$disc)

  # two SNP ALTs whose consensus excludes REF: one degenerate primer covers both
  d <- kasp_site_class("A", c("T", "G"))
  expect_equal(d$class, "degenerate")
  expect_true(d$disc)
  expect_equal(d$code, "K")

  # a lone SNP ALT is always separable, including the A/T and C/G transversions
  # (a strong same-base mismatch on the wrong template)
  c1 <- kasp_site_class("A", "T")
  expect_equal(c1$class, "biallelic")
  expect_true(c1$disc)
  expect_equal(kasp_site_class("C", "G")$class, "biallelic")

  # non-SNP ALTs are the indel path's business; a lone REF is nothing
  expect_equal(kasp_site_class("A", "AT")$class, "indel")
  expect_equal(kasp_site_class("A", "A")$class, "none")
})

test_that("kasp_snp_alts resolves an ALT roster into the set the primer must cover", {
  # a biallelic target: one ALT, nothing degenerate
  one <- kasp_snp_alts(tibble::tibble(chrom = "Chr1", pos = 1L, ref = "A", alt = "G"))
  expect_false(one$degenerate)
  expect_equal(one$alt_pass, "G")
  # an all-SNP multi-allelic target: one degenerate primer must cover BOTH ALTs
  two <- kasp_snp_alts(tibble::tibble(chrom = "Chr1", pos = 1L, ref = "A", alt = "G", alt2 = "C"))
  expect_true(two$degenerate)
  expect_setequal(two$alt_pass, c("G", "C"))
  # an "alt2" equal to the REF base is not an ALT: it must not enter the
  # degenerate 3' code, or the ALT primer would also extend the REF template
  ref_alt <- kasp_snp_alts(tibble::tibble(chrom = "Chr1", pos = 1L, ref = "A", alt = "G", alt2 = "A"))
  expect_false(ref_alt$degenerate)
  expect_false("A" %in% ref_alt$alt_pass)
})
