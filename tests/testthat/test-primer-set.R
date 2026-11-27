# The tube layer takes a typed primer set, not the candidate table a design
# step returns.  These are the two ends of that seam: the accessor's class and
# the converter that claims a plain named vector.

test_that("as_primer_set classes a named vector and is idempotent", {
  p <- as_primer_set(kasp_primers())
  expect_s3_class(p, "primer_set")
  expect_true(is.character(p))
  expect_identical(names(p), names(kasp_primers()))
  expect_identical(as.character(p), unname(as.character(kasp_primers())))

  again <- as_primer_set(p)
  expect_identical(again, p)
  expect_s3_class(again, "primer_set")

  # an empty set is still the typed object: a marker with no oligos is not an
  # error, it is an empty primer_set
  expect_s3_class(as_primer_set(character(0)), "primer_set")
  expect_error(as_primer_set(list(a = "ACGT")), "named character")
  expect_error(as_primer_set(1), "named character")
})

test_that("marker_primers returns the typed set from every method", {
  m <- structure(
    list(primers = list(forward = "ACGT", reverse = "TGCA")),
    class = "marker_unknown"
  )
  p <- marker_primers(m)
  expect_s3_class(p, "primer_set")
  expect_identical(as.character(p), c("ACGT", "TGCA"))
  expect_identical(names(p), c("forward", "reverse"))

  # a marker exposing no oligos is an empty typed set, not NULL
  empty <- marker_primers(structure(list(), class = "marker_unknown"))
  expect_s3_class(empty, "primer_set")
  expect_length(empty, 0L)
})
