# External-tool resolution (options -> project .tools -> cache -> PATH).

test_that("OligoArrayAux is resolved from options(tidymarkers.oligoarrayaux)", {
  prefix <- withr::local_tempdir()
  fs::dir_create(fs::path(prefix, "bin"))
  fs::dir_create(fs::path(prefix, "share", "oligoarrayaux"))
  file.create(fs::path(prefix, "bin", "hybrid-min"))
  withr::local_options(list(tidymarkers.oligoarrayaux = prefix))

  p <- tidymarkers:::oligoarrayaux_paths()
  expect_equal(fs::path_file(p$exe), "hybrid-min")
  expect_equal(as.character(p$bin), as.character(fs::path(prefix, "bin")))
  expect_equal(as.character(p$dat), as.character(fs::path(prefix, "share", "oligoarrayaux")))
  expect_true(has_oligoarrayaux())
})

test_that("a lone binary without data tables still resolves (dat = NA)", {
  prefix <- withr::local_tempdir()
  fs::dir_create(fs::path(prefix, "bin"))
  file.create(fs::path(prefix, "bin", "hybrid-min"))
  withr::local_options(list(tidymarkers.oligoarrayaux = prefix))

  p <- tidymarkers:::oligoarrayaux_paths()
  expect_true(is.na(p$dat))
})

test_that("BLAST+ needs both binaries in the same directory", {
  d <- withr::local_tempdir()
  file.create(fs::path(d, "blastn"))
  withr::local_options(list(tidymarkers.blast = d))
  before <- tidymarkers:::blast_bin() # only blastn -> not a usable pair
  file.create(fs::path(d, "makeblastdb"))
  # resolution is memoised for the session, so a binary that appears afterwards
  # is only seen once the memo is dropped (`marker_install_tools()` does that)
  tidymarkers:::session_memo_clear()
  after <- tidymarkers:::blast_bin()
  if (!is.null(before)) {
    expect_false(identical(as.character(before), as.character(fs::path(d))))
  }
  expect_equal(as.character(after), as.character(fs::path(d)))
})

test_that("tool resolution is memoised and keyed on the option", {
  a <- withr::local_tempdir()
  fs::dir_create(fs::path(a, "bin"))
  file.create(fs::path(a, "bin", "hybrid-min"))
  b <- withr::local_tempdir()
  fs::dir_create(fs::path(b, "bin"))
  file.create(fs::path(b, "bin", "hybrid-min"))
  withr::local_options(list(tidymarkers.oligoarrayaux = a))
  first <- tidymarkers:::oligoarrayaux_paths()

  # the probe is paid once: deleting the binary does not change the answer ...
  fs::file_delete(fs::path(a, "bin", "hybrid-min"))
  expect_identical(tidymarkers:::oligoarrayaux_paths(), first)

  # ... until the memo is dropped
  tidymarkers:::session_memo_clear()
  expect_false(identical(tidymarkers:::oligoarrayaux_paths(), first))

  # a different option value is a different probe, with no clear needed
  withr::local_options(list(tidymarkers.oligoarrayaux = b))
  expect_equal(
    as.character(tidymarkers:::oligoarrayaux_paths()$exe),
    as.character(fs::path(b, "bin", "hybrid-min"))
  )
})

test_that("marker_install_tools() drops the session memo", {
  # a build puts a tool where a memoised "not found" would still say otherwise
  tidymarkers:::session_memo("sentinel", TRUE)
  testthat::local_mocked_bindings(
    oligoarrayaux_works = function(bin, dat) TRUE,
    .package = "tidymarkers"
  )
  tidymarkers::marker_install_tools(dest = withr::local_tempdir(), ask = FALSE)
  expect_length(ls(tidymarkers:::.marker_state$memo, all.names = TRUE), 0L)
})

test_that("marker_install_tools refuses a link where the download and install go", {
  dest <- withr::local_tempdir()
  outside <- fs::path(dest, "outside.txt")
  writeLines("keep", outside)
  for (member in c(
    "oligoarrayaux-3.8.1.tar.gz", "oligoarrayaux-3.8.1", "bin", "share"
  )) {
    link <- fs::path(dest, member)
    skip_if(!file.symlink(outside, link), "this filesystem cannot make a symlink")
    expect_error(
      tidymarkers::marker_install_tools(dest, ask = FALSE, force = TRUE),
      "which is a link",
      class = "rlang_error"
    )
    expect_equal(readLines(outside), "keep")
    file.remove(link)
  }
})

test_that("a missing OligoArrayAux is reported once per session, not per oligo", {
  # a 100-marker design asks oligo_dg() thousands of times; the machine's
  # missing tool is one fact, so the hint is told once
  hints <- tidymarkers:::.marker_state$hints
  rm(list = ls(hints, all.names = TRUE), envir = hints)
  testthat::local_mocked_bindings(
    oligoarrayaux_paths = function() NULL,
    .package = "tidymarkers"
  )
  withr::local_envvar(PATH = "")
  expect_warning(oligo_dg("ACGTACGTACGTTGCA"), "was not found")
  expect_silent(expect_identical(oligo_dg("ACGTACGTACGTTGCA"), NA_real_))
})
