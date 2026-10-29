# =============================================================================
#  data-raw/wheat2019-conversion.R  --  Z8.1: three platforms, one denominator
# =============================================================================
#  Ayalew et al. 2019 (PLoS ONE, 10.1371/journal.pone.0217222, PMC6530864)
#  submitted 50 pre-validated wheat SNPs/INDELs to three design pipelines and
#  reports design success as 39/50 (TaqMan), 49/50 (KASP) and 49/50 (rhAmp).
#  The 11 TaqMan rejections are attributed "mainly" to INDELs longer than 6 bp,
#  but the paper gives no count of such INDELs, so the reason is qualitative
#  and the row is NOT decomposed.  KASP's and rhAmp's single rejection is not
#  explained at all.
#
#  The three rows share one denominator: same 50 loci, three chemistries.  This
#  script asserts that, so a future edit cannot make them look like three
#  independent studies and quietly turn them into a pooled rate.
#
#  Run:  Rscript data-raw/wheat2019-conversion.R
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
w <- c[grepl("^wheat2019_", c$study), ]
stopifnot(nrow(w) == 3L)

## the paper's printed design successes, by platform
stopifnot(w$n_designed[w$study == "wheat2019_taqman"] == 50L)
stopifnot(w$n_success[w$study == "wheat2019_taqman"] == 39L)
stopifnot(w$n_success[w$study == "wheat2019_kasp"] == 49L)
stopifnot(w$n_success[w$study == "wheat2019_rhamp"] == 49L)

## one denominator, three chemistries: never three independent studies
stopifnot(all(w$n_designed == 50L))
stopifnot(all(w$n_success + w$n_failed == w$n_designed))

## nothing is decomposed: the TaqMan reason is "mainly", the others unnamed
stopifnot(all(!w$decomposes), all(!w$fully_decomposed))

cat("wheat2019: 50 loci -> TaqMan 39, KASP 49, rhAmp 49 (design success, not amplification)\n")
