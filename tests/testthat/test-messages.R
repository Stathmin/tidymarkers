# User-facing messages are part of the API surface: snapshot them so a wording
# change shows up as a conscious diff.  Pinned to English here (testthat already
# sets LANGUAGE; the guard documents the intent).

test_that("abort messages are stable", {
  withr::local_envvar(LANGUAGE = "en")
  expect_snapshot(error = TRUE, with_unafloddat("/definitely/not/here", 1))
  expect_snapshot(
    error = TRUE,
    kasp_mix_recipe(n_reactions = 1L, mix_uL = 1, stock_uM = 1)
  )
})
