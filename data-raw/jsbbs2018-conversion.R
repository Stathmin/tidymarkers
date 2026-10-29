# =============================================================================
#  data-raw/jsbbs2018-conversion.R  --  Z8.1: the Exp. 3 failures, decomposed
# =============================================================================
#  Yamagata et al. 2018 (Breed Sci, 10.1270/jsbbs.18048, PMC6198901) reports
#  24 NNNs-HRM primer pairs in Exp. 3 and names BOTH of its failure classes:
#
#    * 3 showed no polymorphism between the parents
#    * 3 showed a dominant segregation pattern
#    * 18 (75 %) segregated co-dominantly
#
#  Those two classes add up (3 + 3 = 6 = 24 - 18), which is what makes this the
#  second row in conversion_rates.csv whose failures decompose.  This script
#  ASSERTS the paper's printed totals against the committed row and against
#  hrm_jsbbs2018_funnel.csv, which data-raw/jsbbs2018.R generates from the
#  supplement.  It writes nothing: a mismatch is a failure, not an average.
#
#  Run:  Rscript data-raw/jsbbs2018-conversion.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

ext <- function(name) {
  readr::read_csv(fs::path(root, "inst", "extdata", name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

c <- ext("conversion_rates.csv")
y <- c[c$study == "jsbbs2018_exp3", ]
stopifnot(nrow(y) == 1L)

## the paper's printed numbers
stopifnot(
  y$n_designed == 24L,
  y$n_success == 18L,
  y$n_fail_pcr == 0L,
  y$n_fail_no_allele == 3L,
  y$n_fail_other == 3L,
  y$n_failed == 6L,
  isTRUE(y$decomposes),
  isTRUE(y$fully_decomposed)
)

## the named classes reconstruct the denominator, both ways
stopifnot(y$n_success + y$n_failed == y$n_designed)
stopifnot(y$n_fail_pcr + y$n_fail_no_allele + y$n_fail_other == y$n_failed)

## the funnel generated from the supplement must agree
f <- ext("hrm_jsbbs2018_funnel.csv")
e3 <- f[f$experiment == "Exp. 3", ]
stopifnot(nrow(e3) == 1L, e3$n_designed == y$n_designed, e3$n_success == y$n_success)

cat("jsbbs2018 Exp. 3: 24 designed -> 18 co-dominant; 3 no polymorphism + 3 dominant\n")
