# =============================================================================
#  server.R  --  stdio MCP server over the tidymarkers tool registry (YP2)
# =============================================================================
#  One JSON-RPC 2.0 message per stdin line; one response line per request (a
#  notification gets none).  Start it with:
#
#    Rscript -e 'cat(system.file("mcp/server.R", package = "tidymarkers"))'
#
#  or, from a checkout, `Rscript inst/mcp/server.R`.  Run it read-only inside a
#  container with the data mounted: see docker/README.md for the mount and the
#  network policy.
# =============================================================================

suppressPackageStartupMessages(library(tidymarkers))

# When run by `Rscript`, `--file=` names this script; fall back to the installed
# location so an `Rscript -e ...` launcher still finds the sibling tools file.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  system.file("mcp", package = "tidymarkers")
}
source(file.path(script_dir, "tools.R"))

input <- file("stdin", "r")
repeat {
  line <- trimws(readLines(input, n = 1L, warn = FALSE))
  if (!length(line)) {
    break
  }
  if (!nzchar(line)) {
    next
  }
  request <- tryCatch(
    jsonlite::fromJSON(line, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(request)) {
    next
  }
  response <- mcp_dispatch(request)
  if (!is.null(response)) {
    cat(jsonlite::toJSON(response, auto_unbox = TRUE, null = "null", force = TRUE), "\n", sep = "")
    flush(stdout())
  }
}
