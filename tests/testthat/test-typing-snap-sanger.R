# =============================================================================
#  SNaPshot vs Sanger: discordance counts on 272 mosquitoes (PMC13020337)
# =============================================================================
#  These tests hold the arithmetic the prose states, and hold the fact that the
#  per-sample table the paper cites is NOT in the supplement.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("SNaPshot vs Sanger: three genes, two populations each", {
  d <- ext("snap_sanger2026_concordance.csv")
  expect_equal(nrow(d), 6L)
  expect_setequal(unique(d$gene), c("vgsc", "Ace-1", "RDL"))
  expect_true(all(c("HN-PY", "AH-SC") %in% d$population))
  expect_true(all(d$n_field_samples_total == 272L))
})

test_that("SNaPshot vs Sanger: the counts are the paper's counts", {
  d <- ext("snap_sanger2026_concordance.csv")
  expect_equal(sum(d$n_discordant[d$gene == "vgsc"]), 5L)
  expect_equal(sum(d$n_discordant[d$gene == "Ace-1"]), 8L)
  expect_equal(sum(d$n_discordant[d$gene == "RDL"]), 0L)
  # and the population split the paper prints
  v <- d[d$gene == "vgsc", ]
  expect_equal(v$n_discordant[v$population == "HN-PY"], 3L)
  expect_equal(v$n_discordant[v$population == "AH-SC"], 2L)
  a <- d[d$gene == "Ace-1", ]
  expect_equal(a$n_discordant[a$population == "HN-PY"], 5L)
  expect_equal(a$n_discordant[a$population == "AH-SC"], 3L)
})

test_that("SNaPshot vs Sanger: no count can exceed the cohort", {
  d <- ext("snap_sanger2026_concordance.csv")
  expect_true(all(d$n_discordant >= 0L))
  expect_true(all(d$n_discordant <= d$n_field_samples_total))
  # "RDL consistent" is a zero, not a missing value
  expect_false(anyNA(d$n_discordant))
})

test_that("SNaPshot vs Sanger: the organism is stated as non-plant", {
  d <- ext("snap_sanger2026_concordance.csv")
  expect_true(all(grepl("Anopheles sinensis", d$organism, fixed = TRUE)))
  expect_true(all(grepl("not a plant", d$organism, fixed = TRUE)))
})
