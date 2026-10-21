test_that("crisprBase registry table is complete and has no seed/PFS slot", {
  path <- test_path("../../data-raw/proposed/crisprbase_registry.csv")
  skip_if_not(file.exists(path), "proposed registry table not present")
  reg <- utils::read.csv(path, stringsAsFactors = FALSE)

  ## the enumeration the report makes: 8 nucleases, Cas13d only, no LwaCas13a
  expect_setequal(
    reg$key,
    c("AsCas12a", "CasRx", "Csm", "enAsCas12a", "MAD7", "SaCas9", "SpCas9", "SpG")
  )
  expect_false("LwaCas13a" %in% reg$key)

  ## slot names parsed from the class sources: no seed and no PFS slot anywhere
  slots <- strsplit(reg$class_slots[1], ",", fixed = TRUE)[[1]]
  expect_false(any(grepl("seed|pfs", slots, ignore.case = TRUE)))
  expect_true(all(c("motifs", "pam_side", "spacer_length") %in% slots))

  ## every row says the same thing about the seed model, because it is a
  ## property of the CLASS, not of the enzyme
  expect_true(all(reg$seed_model == "none (class has no seed slot)"))
  expect_equal(nrow(reg), 8L)

  ## the two RNA-targeting nucleases are the only ones with a PFS column at all,
  ## and both store the catch-all N: a PFS field exists, a PFS MODEL does not
  rna <- reg$target == "RNA"
  expect_setequal(reg$key[rna], c("CasRx", "Csm"))
  expect_true(all(reg$pfs_model[rna] == "none (catch-all N)"))
  expect_true(all(grepl("not applicable", reg$pfs_model[!rna])))

  ## a spacer length is set for every enzyme, and motifs are non-empty
  expect_true(all(reg$spacer_len > 0))
  expect_true(all(reg$n_motifs >= 1L))
})
