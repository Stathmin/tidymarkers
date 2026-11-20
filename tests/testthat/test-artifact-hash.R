# The artifact gate must hash the `code` field (the package source an artifact is
# drawn from), not only the generator and its inputs: otherwise an edit to the
# plot/check code marks no figure stale.  `dev/artifacts.R` is pure enough to
# source (`main()` runs only at `sys.nframe() == 0`); `dev/` is not shipped in
# the build, so the file skips there like `test-dev-reporter.R` does.

artifacts_src <- test_path("..", "..", "dev", "artifacts.R")
skip_if_not(file.exists(artifacts_src), "dev/artifacts.R is not shipped in the build")
source(artifacts_src)

test_that("the input hash covers the code sources, not just generator + inputs", {
  dir <- withr::local_tempdir()
  writeLines("generator", file.path(dir, "gen.R"))
  writeLines("data", file.path(dir, "data.txt"))
  writeLines("code v1", file.path(dir, "code.R"))
  writeLines("artifact", file.path(dir, "figure.svg"))
  entry <- list(
    path = "figure.svg", generator = "gen.R", inputs = "data.txt", code = "code.R"
  )

  withr::with_dir(dir, {
    base <- artifact_input_hash(entry)

    # stable when nothing changes
    expect_identical(artifact_input_hash(entry), base)

    # the code source is part of the hash, not decoration
    no_code <- modifyList(entry, list(code = character()))
    expect_false(identical(artifact_input_hash(no_code), base))

    # and editing it moves the hash
    writeLines("code v2", file.path(dir, "code.R"))
    expect_false(identical(artifact_input_hash(entry), base))

    # restoring it restores the hash
    writeLines("code v1", file.path(dir, "code.R"))
    expect_identical(artifact_input_hash(entry), base)
  })
})
