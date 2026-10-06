# =============================================================================
#  dev/render-all.R  --  (re)draw every generated figure in one pass
# =============================================================================
#  The figures are code, not art: every render script reads the engine, so the
#  diagrams follow the pipeline instead of being maintained by hand.  This is the
#  one entry point that runs them all, so `dev/check.R` can refresh them on every
#  gate run and nothing drifts.  Steps are independent and guarded - a missing
#  tool or unreachable data is a message, not a failure - so it also runs on a
#  bare machine.
#
#  Run:  Rscript dev/render-all.R [--skip-heavy] [--skip-example]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

args <- commandArgs(trailingOnly = TRUE)
heavy <- !("--skip-heavy" %in% args)
example <- !("--skip-example" %in% args)

rscript <- file.path(R.home("bin"), "Rscript")
run <- function(file, label, needs_genome = FALSE, env = character()) {
  message("== ", label, " ==")
  if (needs_genome && !nzchar(Sys.which("hybrid-min"))) {
    message("  skip (OligoArrayAux not on PATH)")
    return(invisible(FALSE))
  }
  status <- system2(rscript, file, env = env)
  ok <- identical(as.integer(status), 0L)
  if (!ok) message("  FAILED (status ", status, ")")
  invisible(ok)
}

# the decision trees are the committed, README-facing images - always refreshed
run("dev/render-decision.R", "decision trees -> README images",
  env = if (example) character() else "TIDYKASP_SKIP_EXAMPLE=1"
)
if (heavy) {
  run("dev/render-targets.R", "example-target panel", needs_genome = TRUE)
  run("dev/render-plots.R", "plot smoke", needs_genome = TRUE)
  run("dev/render-screen.R", "screen plot", needs_genome = TRUE)
  run("dev/render-real.R", "real FRI demo", needs_genome = TRUE)
}
