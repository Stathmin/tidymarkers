# `abort_msg()` is the single chokepoint for translated abort messages.  These
# tests pin its contract: one literal in, the same rendered condition out, with
# cli bullet names, condition class, call and parent chain preserved.

test_that("abort_msg keeps a single message", {
  err <- rlang::catch_cnd(abort_msg("plain message"))
  expect_s3_class(err, "rlang_error")
  expect_identical(conditionMessage(err), "plain message")
})

test_that("abort_msg interpolates from the caller frame", {
  raise <- function(x) abort_msg("bad {x}", i = "hint about {x}")
  err <- rlang::catch_cnd(raise("value"))
  expect_match(conditionMessage(err), "bad value", fixed = TRUE)
  expect_identical(err$body, c(i = "hint about value"))
})

test_that("abort_msg preserves the condition class", {
  expect_error(abort_msg("nope", class = "tidymarkers_test_error"),
    class = "tidymarkers_test_error"
  )
  err <- rlang::catch_cnd(abort_msg("nope", class = c("a_error", "b_error")))
  expect_s3_class(err, "a_error")
  expect_s3_class(err, "b_error")
})

test_that("abort_msg chains the parent condition", {
  parent <- rlang::catch_cnd(rlang::abort("root cause"))
  err <- rlang::catch_cnd(abort_msg("wrapper", parent = parent))
  expect_match(conditionMessage(err), "wrapper", fixed = TRUE)
  expect_identical(err$parent, parent)
})

test_that("abort_msg preserves unknown bullet names", {
  err <- rlang::catch_cnd(abort_msg("top", weird = "bullet text"))
  expect_identical(names(err$body), "weird")
  expect_identical(unname(err$body), "bullet text")
})

test_that("abort_msg looks the message up in the compiled Russian catalogue", {
  english <- "Bundled restriction-enzyme catalogue not found."
  mo <- system.file("po", "ru", "LC_MESSAGES", "R-tidymarkers.mo",
    package = "tidymarkers"
  )
  skip_if(!nzchar(mo), "no compiled Russian catalogue")
  withr::local_envvar(LANGUAGE = "ru")
  ru <- gettext(english, domain = "R-tidymarkers")
  skip_if(identical(ru, english), "Russian catalogue not installed in this run")
  err <- rlang::catch_cnd(abort_msg(english))
  expect_identical(conditionMessage(err), ru)
})
