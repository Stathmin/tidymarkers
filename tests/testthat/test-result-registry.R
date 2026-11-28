# Result generics dispatch through a class -> function registry, so the package
# carries one S3 method per generic instead of one per result class.

test_that("the plot registry is a named list of layouts", {
  reg <- marker_plot_registry()
  expect_true(all(startsWith(names(reg), "marker_")))
  expect_true(all(vapply(reg, is.list, logical(1))))
  expect_identical(anyDuplicated(names(reg)), 0L)
})

test_that("plot() is the only plot S3 method and lookup walks the class vector", {
  expect_null(getS3method("plot", "marker_assay", optional = TRUE))
  expect_true(is.function(getS3method("plot", "marker", optional = TRUE)))
  wrapped <- structure(list(), class = c("marker_explain", "marker_assay", "marker"))
  expect_identical(
    marker_registry_lookup(marker_plot_registry(), wrapped, "plot"),
    marker_plot_registry()$marker_assay
  )
})

test_that("result objects carry marker as their last class", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  a <- suppressWarnings(design_asq_marker(f$genome, f$snp, cores = 1L))
  expect_identical(class(a)[1L], "marker_assay")
  expect_identical(class(a)[length(class(a))], "marker")
})

test_that("an unknown marker class aborts with the internal class", {
  x <- structure(list(), class = c("marker_unknown", "marker"))
  expect_error(plot(x), class = "tidymarkers_registry_error")
})

test_that("print and map_spec registries mirror the plot contract", {
  pr <- marker_print_registry()
  ms <- marker_map_spec_registry()
  expect_identical(anyDuplicated(names(pr)), 0L)
  expect_identical(anyDuplicated(names(ms)), 0L)
  expect_true(all(vapply(pr, is.function, logical(1))))
  expect_true(all(vapply(ms, is.function, logical(1))))
  expect_null(getS3method("print", "marker_taqman", optional = TRUE))
  expect_true(is.function(getS3method("print", "marker", optional = TRUE)))
  expect_true(is.function(getS3method("marker_map_spec", "marker", optional = TRUE)))
})

test_that("a marker class without a print entry still reaches the registry error", {
  x <- structure(list(), class = c("marker_unknown", "marker"))
  expect_error(print(x), class = "tidymarkers_registry_error")
})
