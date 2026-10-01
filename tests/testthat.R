# Test runner for the package:  Rscript tests/testthat.R
# The package must be installed (R CMD INSTALL) so test_check() can load it;
# testthat runs the tests in an environment that inherits the package namespace,
# so internal (non-exported) helpers are reachable without `:::`.
library(testthat)
library(tidykasp)
test_check("tidykasp")
