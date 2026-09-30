# Package helpers: class predicate, all-pass predicate, env helper.

mk_assay <- function() {
  p <- c(AS_ref = "TGACTGCAAGTCCGATCAGTA",
         AS_alt = "TGACTGCAAGTCCGATCAGTT",
         COM    = "GCTTAGCCATGGTCAATGCTA")
  structure(list(
    snp = tibble::tibble(chrom = "Chr4", pos = 1L, ref = "G", alt = "A"),
    strand = "+",
    assay = tibble::tibble(amplicon = 90L, as_ref = p[["AS_ref"]], as_alt = p[["AS_alt"]],
                           com = p[["COM"]], com_tm = 60, as_gc = 0.5, com_gc = 0.5),
    thermo = tibble::tibble(name = names(p), tm = c(57, 57, 60), gc = 0.5, hairpin_dg = -1),
    complementarity = tibble::tibble(a = rep(names(p), each = 3),
                                     b = rep(names(p), times = 3),
                                     dg = c(-2, -3, -4, -3, -2, -4, -4, -4, -3)),
    amplicons = tibble::tibble(primer = c("AS_ref", "AS_alt"),
                               n_products = c(1L, 0L), width = c(90L, NA_integer_)),
    efficiency = tibble::tibble(discrimination = 0.9),
    specificity = tibble::tibble(primer = names(p), seed = "ACGT", hits = 1L)
  ), class = "kasp_assay")
}

test_that("is_kasp_assay recognises the class", {
  expect_true(is_kasp_assay(mk_assay()))
  expect_false(is_kasp_assay(list()))
  expect_false(is_kasp_assay(NULL))
})

test_that("passes_all reflects the acceptance checks", {
  a <- mk_assay()
  expect_true(passes_all(a))
  a$assay$amplicon <- 500L          # break the amplicon criterion
  expect_false(passes_all(a))
})

test_that("with_unafloddat sets and restores UNAFOLDDAT", {
  d <- tempfile("dat"); dir.create(d)
  old <- Sys.getenv("UNAFOLDDAT", unset = NA)
  seen <- with_unafloddat(d, Sys.getenv("UNAFOLDDAT"))
  expect_equal(seen, d)
  expect_equal(Sys.getenv("UNAFOLDDAT", unset = NA), old)   # restored afterwards
  expect_error(with_unafloddat(file.path(d, "nope"), 1), "not found")
})
