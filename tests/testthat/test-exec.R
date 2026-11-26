# `run_external()` is the one seam for cold external binaries (ADR 0007 §5):
# a bad exit is data, a hang is a typed error.

test_that("run_external() returns the exit status and output", {
  res <- run_external(
    file.path(R.home("bin"), "Rscript"),
    c("-e", "cat('hi'); quit(status = 3)"),
    timeout = 60
  )
  expect_identical(as.integer(res$status), 3L)
  expect_match(res$stdout, "hi")
})

test_that("run_external() kills a hung binary and names the timeout", {
  t0 <- Sys.time()
  err <- tryCatch(
    run_external(
      file.path(R.home("bin"), "Rscript"),
      c("-e", "Sys.sleep(60)"),
      timeout = 2
    ),
    error = identity
  )
  expect_s3_class(err, "tidymarkers_timeout")
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 30)
  expect_match(conditionMessage(err), "timed out after 2 s", fixed = TRUE)
})
