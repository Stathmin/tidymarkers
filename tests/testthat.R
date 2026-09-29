# Test runner for the (non-package) demo:  Rscript tests/testthat.R
# The design functions are sourced directly instead of building a package.
library(testthat)
library(here)
for (f in c("kasp_functions.R", "kasp_io.R")) source(here::here("R", f))
test_dir(here::here("tests", "testthat"))
