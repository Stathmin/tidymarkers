# Test runner for the (non-package) demo:  Rscript tests/testthat.R
# The design functions are sourced directly instead of building a package.
library(testthat)
library(here)
source(here::here("R", "kasp_functions.R"))
test_dir(here::here("tests", "testthat"))
