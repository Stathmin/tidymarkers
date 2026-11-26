# Result generics dispatch through a class -> function registry, so the package
# carries one S3 method per generic instead of one per result class.

test_that("the plot registry is a named list of functions", {
  reg <- marker_plot_registry()
  expect_true(all(startsWith(names(reg), "marker_")))
  expect_true(all(vapply(reg, is.function, logical(1))))
  expect_identical(anyDuplicated(names(reg)), 0L)
})

test_that("plot() is the only plot S3 method and lookup walks the class vector", {
  expect_null(getS3method("plot", "marker_assay", optional = TRUE))
  expect_true(is.function(getS3method("plot", "marker", optional = TRUE)))
  wrapped <- structure(list(), class = c("marker_explain", "marker_assay", "marker"))
  expect_identical(
    marker_registry_lookup(marker_plot_registry(), wrapped, "plot"),
    marker_plot_assay
  )
})

test_that("result objects carry marker as their last class", {
  f <- snp_fixture()
  a <- suppressWarnings(design_asq_marker(f$genome, f$snp, cores = 1L))
  expect_identical(class(a)[1L], "marker_assay")
  expect_identical(class(a)[length(class(a))], "marker")
})

test_that("an unknown marker class aborts with the internal class", {
  x <- structure(list(), class = c("marker_unknown", "marker"))
  expect_error(plot(x), class = "tidymarkers_registry_error")
})
