# =============================================================================
#  marker_exec.R  --  external command seam
# =============================================================================
#  Everything cold that shells out goes through `run_external()`: it reports a
#  status instead of throwing on a bad exit, and a hung binary dies at
#  `timeout` seconds with the `tidymarkers_timeout` class rather than holding
#  the caller (and CI) until the runner is killed.
#
#  OligoArrayAux stays on `system2()` in `oligo_dg()`: that is the hottest loop
#  in design and a processx object costs about twice as much per oligo there.
#  `Rbowtie`/`rBLAST` do not expose their process, so their scans have no
#  timeout either; both are named limits in the docs.
# =============================================================================

#' Run an external command with a timeout
#'
#' @param command executable name or path.
#' @param args character vector of arguments.
#' @param timeout seconds before the process tree is killed.
#' @param dir working directory, or `NULL` for the current one.
#' @param env named character vector of environment overrides.
#' @param stdout,stderr `"|"` to capture, or a file path.
#' @return `processx::run()`'s list (`status`, `stdout`, `stderr`, timing).
#'   A bad exit status is returned, not thrown; only a timeout throws, with
#'   class `tidymarkers_timeout`.
#' @keywords internal
#' @noRd
run_external <- function(command, args = character(), timeout = 300,
                         dir = NULL, env = NULL, stdout = "|", stderr = "|") {
  res <- tryCatch(
    processx::run(
      command,
      args = as.character(args),
      error_on_status = FALSE,
      timeout = timeout,
      wd = dir,
      env = env,
      stdout = stdout,
      stderr = stderr
    ),
    error = function(e) stop(e)
  )
  # with error_on_status = FALSE a timeout is a result flag (status -9), not an
  # error, so it has to be turned into the typed condition here
  if (isTRUE(res$timeout)) {
    stop(structure(
      class = c("tidymarkers_timeout", "error", "condition"),
      list(
        message = paste0("'", command, "' timed out after ", timeout, " s"),
        call = NULL,
        command = command,
        timeout = timeout
      )
    ))
  }
  res
}
