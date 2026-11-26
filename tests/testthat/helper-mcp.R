# Shared helpers for the MCP tool-layer tests (test-mcp.R,
# test-mcp-validation.R).
#
# inst/mcp/tools.R is not part of the namespace, so source it into an
# environment that sees the package; that is the same lookup the server's
# `library(tidymarkers)` gives it.
.mcp_env <- function() {
  env <- new.env(parent = asNamespace("tidymarkers"))
  path <- system.file("mcp", "tools.R", package = "tidymarkers")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "mcp", "tools.R")
  }
  sys.source(path, envir = env)
  env
}

.mcp_call <- function(env, name, arguments) {
  env$mcp_dispatch(list(
    jsonrpc = "2.0", id = 1L, method = "tools/call",
    params = list(name = name, arguments = arguments)
  ))
}
