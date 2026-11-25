# Warm the demo cache that the test suite reads: the TAIR10 genome, the
# 4:1.5-1.56 Mb variation panel and one gene record.
#
# `tests/testthat/helper-demo-cache.R` states the contract -- a test must never
# make that fetch -- so a cold cache turns the arms, primer3, plan and batch
# tests into SKIPs instead of a network round-trip inside the suite.
#
# Best-effort by design: a getter that cannot reach Ensembl or EBI prints a
# GitHub Actions warning and the script exits 1.  CI runs it with
# `continue-on-error`, so a network failure degrades to those SKIPs rather than
# a red job.
#
# Run it on a workstation from the package root:
#   Rscript dev/warm-demo-cache.R

pkgload::load_all(".")

warm <- list(
  genome = function() ath_genome(),
  variants = function() ensembl_variants("4", 1500000L, 1560000L),
  gene = function() ensembl_gene("AT4G00650")
)

failed <- character()
for (nm in names(warm)) {
  msg <- tryCatch(
    {
      invisible(warm[[nm]]())
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (!is.null(msg)) {
    failed <- c(failed, nm)
    cat("::warning title=Demo cache warm failed::", nm, ": ",
      gsub("[\r\n]+", " ", msg), "\n",
      sep = ""
    )
  }
}

warm_now <- setdiff(names(warm), failed)
cat(sprintf(
  "demo cache: %d of %d warm%s\n", length(warm_now), length(warm),
  if (length(warm_now)) paste0(" (", paste(warm_now, collapse = ", "), ")") else ""
))
quit(status = as.integer(length(failed) > 0L))
