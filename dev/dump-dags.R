# =============================================================================
#  dev/dump-dags.R  --  every scheme's decision DAG as text, for diffing
# =============================================================================
#  Writes marker_decision_text() for every scheme to results/dags.txt (git-ignored)
#  so an agent can diff its own scheme before/after an edit:
#
#      Rscript dev/dump-dags.R
#      cp results/dags.txt /tmp/before.txt     # snapshot, then edit the engine
#      Rscript dev/dump-dags.R
#      git --no-pager diff --no-index /tmp/before.txt results/dags.txt
#
#  Offline and cheap: it reads the generated schema/registry only - no genome, no
#  network, no OligoArrayAux.  Run it with the renv library on R_LIBS (see
#  AGENTS.md), from the package root; `source("dev/dump-dags.R")` also works.
# =============================================================================

root <- local({
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
})
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

fs::dir_create("results")
out <- "results/dags.txt"

# one greppable block per scheme: a `#### <scheme> (<family>)` header, then the
# full text export.  Reading the family here also surfaces a scheme whose family
# gained or lost its allele-specific checks.
blocks <- lapply(marker_schemes()$scheme, function(sc) {
  c(paste0("#### ", sc, " (", marker_scheme(sc)$family, ")"), marker_decision_text(sc), "")
})
writeLines(unlist(blocks), out)
cat("wrote ", out, " for ", nrow(marker_schemes()), " schemes\n", sep = "")
