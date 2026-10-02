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

  # a lone Watson-Crick complement cannot be split by a 3'-terminal base
  c1 <- kasp_site_class("A", "T")
  expect_equal(c1$class, "complementary")
  expect_false(c1$disc)

  # non-SNP ALTs are the indel path's business; a lone REF is nothing
  expect_equal(kasp_site_class("A", "AT")$class, "indel")
  expect_equal(kasp_site_class("A", "A")$class, "none")
})
