# =============================================================================
#  dev/check-reporter.R  --  the shared reporter for the gate scripts
# =============================================================================
#  `ok()` prints a passing check; `note()` appends a failure message to the
#  script-level `failures` vector that dev/check-docs.R and dev/check-artifacts.R
#  define.  Both run at the top level under Rscript, so the `<<-` resolves to
#  their own `failures`, and every call site stays unchanged.
#
#  This is not dev/reporter.R, which formats the test-timing report.
# =============================================================================

ok <- function(...) message("  ok   ", ...)

note <- function(...) failures <<- c(failures, paste0(...))
