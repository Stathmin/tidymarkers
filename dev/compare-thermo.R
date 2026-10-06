# =============================================================================
#  dev/compare-thermo.R  --  reconcile the Tm window against the NN models
# =============================================================================
#  Runs a panel of representative oligos through several published DNA
#  nearest-neighbour parameter sets (the Primer3/oligotm family our Tm comes
#  from) and writes the Tm spread.  A spread well under the AS window
#  (KASP_LIMITS$as_tm, 10 C wide) means the threshold is model-robust; a large
#  spread means the window is model-dependent.  Also prints the provenance table
#  (kasp_priors()) so the numbers and their sources sit side by side.
#
#  Run:  Rscript dev/compare-thermo.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

# a panel: the KASP tails, a GC-balanced and an AT-rich primer, a short and a
# long one - the range of what the designer emits
seqs <- c(
  "GAAGGTGACCAAGTTCATGCT", "GAAGGTCGGAGTCAACGGATT",
  "ACGTACGTACGTACGTACGT", "TTTTGGCGCGCCGATATATAT",
  "ATATATATATATATATATATAT", "GCGCGCGCGCGCGCGCGCGC"
)
x <- thermo_crosscheck(seqs)
cat("\n== Tm across NN tables (C) ==\n")
print(as.data.frame(x), row.names = FALSE)
cat("\n== provenance ==\n")
print(as.data.frame(kasp_priors()), row.names = FALSE)

out <- fs::path(root, "results", "thermo-reconcile.csv")
fs::dir_create(fs::path_dir(out))
readr::write_csv(x, out)
cat("\nwrote", out, "\n")
