# `hash_object()`/`hash_file()` are the one seam that names the digest algorithm
# (ADR 0007 §5).  A path and a value are different inputs: `hash_file()` reads
# the bytes, `hash_object()` serialises the R object.

test_that("hash_object() hashes the object, not a same-named path", {
  # `digest`'s default serialisation is what every existing manifest used, so
  # the object hash must stay `digest::digest()`'s answer, not raw bytes.
  expect_identical(
    hash_object("abc"),
    digest::digest("abc", algo = "sha256")
  )
  f <- withr::local_tempfile()
  writeLines("abc", f)
  expect_false(identical(hash_object(f), hash_file(f)))
})

test_that("hash_file() hashes the bytes and returns NA when absent", {
  f <- withr::local_tempfile()
  writeBin(charToRaw("abc"), f)
  expect_equal(
    hash_file(f),
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
  )
  expect_true(is.na(hash_file(file.path(tempdir(), "tidymarkers-no-such-file"))))
  expect_true(is.na(hash_file(NA_character_)))
})
