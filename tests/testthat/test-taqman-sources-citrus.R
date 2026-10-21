# TaqMan / 5'-nuclease science extracts (citrus2020). Data-only, from
# inst/extdata; the generator data-raw/citrus2020.R re-derives every number
# from the fetched JATS before writing.

tq_read <- function(name) {
  utils::read.csv(testthat::test_path("..", "..", "inst", "extdata", name),
    stringsAsFactors = FALSE
  )
}

test_that("citrus2020 is a complete 48 x 8 call matrix", {
  q <- tq_read("taqman_citrus2020_calls.csv")
  expect_identical(nrow(q), 384L)
  expect_identical(length(unique(q$marker)), 8L)
  expect_identical(length(unique(q$cultivar_no)), 48L)
  expect_false(any(is.na(q$genotype) | !nzchar(q$genotype)))
  expect_true(all(q$genotype %in% c("AA", "AG", "GG", "AT", "TT", "CC", "CT", "CG")))
  # each marker's calls use exactly the two bases its two probes underline
  pr <- tq_read("taqman_citrus2020_probes.csv")
  for (i in seq_len(nrow(pr))) {
    g <- q$genotype[q$marker == pr$marker[i]]
    expect_identical(
      sort(unique(strsplit(paste(g, collapse = ""), "")[[1]])),
      sort(c(pr$fam_marked[i], pr$hex_marked[i]))
    )
  }
  # (the paper's separate "Polymorphism" sub-table lists the same pairs in a
  # DIFFERENT marker order and disagrees on three rows -- SI306 A/G vs probes
  # C/G, SI210 A/G vs A/C, SI381 T/C vs the printed T in both probes -- so the
  # calls, not that sub-table, are what the probes are checked against)
})

test_that("citrus2020's probes are the printed sequences, not a centred default", {
  pr <- tq_read("taqman_citrus2020_probes.csv")
  expect_identical(nrow(pr), 8L)
  expect_true(all(pr$probe_len_fam < 21L))
  # 5 of 8 FAM probes are shorter than the package's 18-nt window floor
  expect_identical(sum(pr$probe_len_fam < 18L), 5L)
  # no family places the allele base centrally
  expect_true(all(stats::na.omit(pr$centrality) < 0.99))
  expect_identical(length(stats::na.omit(pr$centrality)), 6L)
  # exactly one FAM/HEX pair is printed with a one-base shift
  expect_identical(sum(pr$probe_len_fam == pr$probe_len_hex), 7L)
  expect_identical(pr$marker[pr$probe_len_fam != pr$probe_len_hex], "SI381")
  # the underlined allele base is one of the two the probes carry
  expect_true(all(pr$fam_marked %in% c("A", "C", "G", "T")))
  expect_true(all(pr$hex_marked %in% c("A", "C", "G", "T")))
})

test_that("citrus2020's amplicons are the paper's 60-99 bp window", {
  a <- tq_read("taqman_citrus2020_amplicons.csv")
  expect_identical(nrow(a), 8L)
  expect_true(all(a$product_bp >= 60L & a$product_bp <= 99L))
  expect_identical(sum(a$product_bp), 642L)
})
