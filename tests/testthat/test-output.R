# Consolidated results writer: a few grouped files, order formats as subsets.

test_that("write_kasp_results writes grouped files and order subsets", {
  p <- c(
    AS_ref = "TGACTGCAAGTCCGATCAGTA",
    AS_alt = "TGACTGCAAGTCCGATCAGTT",
    COM = "GCTTAGCCATGGTCAATGCTA"
  )
  assay <- list(
    snp = tibble::tibble(chrom = "Chr4", pos = 269090L, ref = "G", alt = "A"),
    strand = "+",
    assay = tibble::tibble(
      rank = 1L, strand = "+", as_ref = p[["AS_ref"]], as_alt = p[["AS_alt"]],
      com = p[["COM"]], discrimination = 0.9, as_gc = 0.5, amplicon = 90L,
      as_len = 21L, as_n_var = 0L, as_n_3p = 0L, n_var = 0L, as_foot_dropout = 0
    ),
    thermo = tibble::tibble(
      name = names(p), primer = unname(p), len = 21L,
      tm = c(58, 58, 60), gc = 0.5,
      hairpin_dg = c(-1, -1, -1), terminal_dg = c(-7, -7, -8)
    ),
    specificity = tibble::tibble(primer = names(p), seed = "ACCTAGCTA", hits = 1L),
    amplicons = tibble::tibble(
      primer = c("AS_ref", "AS_alt"),
      n_products = c(1L, 0L), width = c(90L, NA_integer_)
    ),
    complementarity = tibble::tibble(
      a = rep(names(p), each = 3),
      b = rep(names(p), times = 3), dg = -3
    ),
    blast = NULL
  )
  dir <- withr::local_tempdir(pattern = "res")
  paths <- write_kasp_results(assay, list(symbol = "FRI"), dir = dir)
  for (f in c(
    "kasp_design.csv", "kasp_primer_qc.csv", "kasp_dimer_dg.csv",
    "kasp_FRI_Chr4_269090_tidy.csv", "kasp_FRI_Chr4_269090_idt.csv",
    "kasp_FRI_Chr4_269090_lgc.csv"
  )) {
    expect_true(file.exists(file.path(dir, f)), info = f)
  }
  # primer QC joins thermo + specificity + in-silico PCR into one row per primer
  qc <- readr::read_csv(file.path(dir, "kasp_primer_qc.csv"), show_col_types = FALSE)
  expect_equal(nrow(qc), 3L)
  expect_true(all(c("tm", "terminal_dg", "hits", "n_products") %in% names(qc)))
})
