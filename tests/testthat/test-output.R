# Consolidated results writer: a few grouped files, order formats as subsets.

output_fixture <- function() {
  p <- c(
    AS_ref = "TGACTGCAAGTCCGATCAGTA",
    AS_alt = "TGACTGCAAGTCCGATCAGTT",
    COM = "GCTTAGCCATGGTCAATGCTA"
  )
  list(
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
}

test_that("write_marker_results writes grouped files and order subsets", {
  assay <- output_fixture()
  dir <- withr::local_tempdir(pattern = "res")
  paths <- write_marker_results(assay, list(symbol = "FRI"), dir = dir)
  for (f in c(
    "marker_design.csv", "marker_primer_qc.csv", "marker_dimer_dg.csv",
    "marker_FRI_Chr4_269090_tidy.csv", "marker_FRI_Chr4_269090_idt.csv",
    "marker_FRI_Chr4_269090_lgc.csv"
  )) {
    expect_true(file.exists(file.path(dir, f)), info = f)
  }
  # primer QC joins thermo + specificity + in-silico PCR into one row per primer
  qc <- readr::read_csv(file.path(dir, "marker_primer_qc.csv"), show_col_types = FALSE)
  expect_equal(nrow(qc), 3L)
  expect_true(all(c("tm", "terminal_dg", "hits", "n_products") %in% names(qc)))
  # the old name is a deprecation shim over the new one
  expect_warning(
    write_kasp_results(assay, list(symbol = "FRI"), dir = dir),
    "deprecated"
  )

  # a path-like gene symbol cannot write outside `dir`
  escaped <- write_marker_results(assay, list(symbol = "../../escape"), dir = dir)
  expect_true(all(startsWith(normalizePath(escaped), normalizePath(dir))))
  expect_false(file.exists(file.path(dirname(dir), "escape_Chr4_269090_tidy.csv")))
})

test_that("the output writers refuse to write through a link", {
  assay <- output_fixture()
  gene <- list(symbol = "FRI")
  members <- c(
    "marker_design.csv", "marker_primer_qc.csv", "marker_dimer_dg.csv",
    "marker_FRI_Chr4_269090_tidy.csv", "marker_FRI_Chr4_269090_idt.csv",
    "marker_FRI_Chr4_269090_lgc.csv"
  )
  for (member in members) {
    dir <- withr::local_tempdir()
    outside <- fs::path(withr::local_tempdir(), member)
    writeLines("keep", outside)
    skip_if_not(
      file.symlink(outside, fs::path(dir, member)),
      "this filesystem cannot make a symlink"
    )
    expect_error(write_marker_results(assay, gene, dir = dir), "which is a link", info = member)
    expect_identical(readLines(outside), "keep", info = member)
  }

  # export_primers() is exported on its own, so it guards its own three names
  dir <- withr::local_tempdir()
  outside <- fs::path(withr::local_tempdir(), "marker_FRI_Chr4_269090_idt.csv")
  writeLines("keep", outside)
  skip_if_not(
    file.symlink(outside, fs::path(dir, "marker_FRI_Chr4_269090_idt.csv")),
    "this filesystem cannot make a symlink"
  )
  expect_error(export_primers(assay, gene, dir = dir), "which is a link")
  expect_identical(readLines(outside), "keep")
})

# A cache and a rendered file are declared names too (ADR 0033): the guard the
# writers share refuses a link and returns an ordinary spelling unchanged.
test_that("declared_path refuses a link and keeps an ordinary path", {
  dir <- withr::local_tempdir()
  outside <- fs::path(withr::local_tempdir(), "outside.json")
  writeLines("keep", outside)
  link <- fs::path(dir, "cache.json")
  skip_if(!file.symlink(outside, link), "this filesystem cannot make a symlink")
  expect_error(declared_path(link), "which is a link")
  expect_identical(readLines(outside), "keep")
  plain <- fs::path(dir, "plain.json")
  expect_identical(declared_path(plain), plain)
})
