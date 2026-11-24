# The validation benchmark is a dev script, not shipped code, but its contract is
# what ZY21 rests on: a fixed task count, a verdict from a closed vocabulary, and
# no `not_run`/`gap` without a stated reason.  `dev/` is not in the build, so the
# file skips there like `test-artifact-hash.R` does.

benchmark_src <- test_path("..", "..", "dev", "validation-benchmark.R")
skip_if_not(file.exists(benchmark_src), "dev/validation-benchmark.R is not shipped in the build")
source(benchmark_src)

# One benchmark run per test file: the set designs a panel and reads the shipped
# tables, so repeating it per test triples the file's wall time.
benchmark_once <- local({
  value <- NULL
  function() {
    if (is.null(value)) value <<- benchmark_tasks()
    value
  }
})

test_that("the benchmark reproduces its own task set", {
  first <- benchmark_once()
  second <- benchmark_tasks()

  expect_equal(nrow(first), 25L)
  expect_identical(first$verdict, second$verdict)
  expect_identical(first$task, second$task)
  expect_true(all(first$verdict %in% c("pass", "fail", "not_run", "gap")))
  expect_equal(sum(first$family == "soybean"), 13L)
})

test_that("every unestablished or unrunnable task states why", {
  x <- benchmark_once()
  silent <- x$verdict %in% c("not_run", "gap") & is.na(x$note)
  expect_false(any(silent))
})

test_that("the recorded CSV round-trips", {
  x <- benchmark_once()
  path <- benchmark_write(x, tempfile(fileext = ".csv"))
  back <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "")
  expect_equal(back$task, x$task)
  expect_equal(back$verdict, x$verdict)
  expect_equal(back$observed, x$observed)
})
