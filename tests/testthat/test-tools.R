# External-tool resolution (options -> project .tools -> cache -> PATH).

test_that("OligoArrayAux is resolved from options(tidykasp.oligoarrayaux)", {
  prefix <- withr::local_tempdir()
  fs::dir_create(fs::path(prefix, "bin"))
  fs::dir_create(fs::path(prefix, "share", "oligoarrayaux"))
  file.create(fs::path(prefix, "bin", "hybrid-min"))
  withr::local_options(list(tidykasp.oligoarrayaux = prefix))

  p <- tidykasp:::oligoarrayaux_paths()
  expect_equal(fs::path_file(p$exe), "hybrid-min")
  expect_equal(as.character(p$bin), as.character(fs::path(prefix, "bin")))
  expect_equal(as.character(p$dat), as.character(fs::path(prefix, "share", "oligoarrayaux")))
  expect_true(has_oligoarrayaux())
})

test_that("a lone binary without data tables still resolves (dat = NA)", {
  prefix <- withr::local_tempdir()
  fs::dir_create(fs::path(prefix, "bin"))
  file.create(fs::path(prefix, "bin", "hybrid-min"))
  withr::local_options(list(tidykasp.oligoarrayaux = prefix))

  p <- tidykasp:::oligoarrayaux_paths()
  expect_true(is.na(p$dat))
})

test_that("BLAST+ needs both binaries in the same directory", {
  d <- withr::local_tempdir()
  file.create(fs::path(d, "blastn"))
  withr::local_options(list(tidykasp.blast = d))
  before <- tidykasp:::blast_bin() # only blastn -> not a usable pair
  file.create(fs::path(d, "makeblastdb"))
  after <- tidykasp:::blast_bin()
  if (!is.null(before)) {
    expect_false(identical(as.character(before), as.character(fs::path(d))))
  }
  expect_equal(as.character(after), as.character(fs::path(d)))
})
