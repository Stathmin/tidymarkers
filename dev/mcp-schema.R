# =============================================================================
#  dev/mcp-schema.R  --  ZX5.3: publish the `tools/list` contract as JSON Schema
# =============================================================================
#  ADR 0004 parked inst/mcp/schema.json "until a second consumer that cannot
#  load R".  The operator decided to publish it now.  It stays a contract file:
#  nothing in the package reads it at run time, the server is still R, and
#  tests/testthat/test-mcp.R fails on drift exactly as it already does for the
#  registry itself.
#
#  The document is JSON Schema draft 2020-12 for the `tools/list` *result*:
#  `{ tools: [ { name, title, description, inputSchema }, ... ] }`, with the
#  tool count fixed and the `name` enum taken from `mcp_tool_specs()`.
#
#  Run:  Rscript dev/mcp-schema.R            # print to stdout
#        Rscript dev/mcp-schema.R --write    # write inst/mcp/schema.json
#        Rscript dev/mcp-schema.R --check    # exit 1 if the file has drifted
# =============================================================================

schema_root <- function() {
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))))
  }
  wd <- getwd()
  repeat {
    if (file.exists(file.path(wd, "DESCRIPTION")) && dir.exists(file.path(wd, "R"))) {
      return(wd)
    }
    parent <- dirname(wd)
    if (identical(parent, wd)) {
      return(getwd())
    }
    wd <- parent
  }
}

root <- schema_root()
setwd(root)

mcp_specs <- function() {
  suppressMessages(pkgload::load_all(root, quiet = TRUE))
  env <- new.env(parent = asNamespace("tidymarkers"))
  sys.source(file.path(root, "inst", "mcp", "tools.R"), envir = env)
  env$mcp_tool_specs()
}

build_schema <- function(specs) {
  n <- length(specs)
  list(
    `$schema` = "https://json-schema.org/draft/2020-12/schema",
    title = "tidymarkers MCP tools/list result",
    description = paste(
      "The result of a JSON-RPC `tools/list` call on the tidymarkers MCP",
      "server (inst/mcp/server.R). Generated from mcp_tool_specs() by",
      "dev/mcp-schema.R; tests/testthat/test-mcp.R fails on drift."
    ),
    type = "object",
    required = "tools",
    additionalProperties = FALSE,
    properties = list(tools = list(
      type = "array",
      minItems = as.integer(n),
      maxItems = as.integer(n),
      items = list(
        type = "object",
        required = c("name", "title", "description", "inputSchema"),
        additionalProperties = FALSE,
        properties = list(
          name = list(
            type = "string",
            enum = names(specs),
            description = "the registered tool name"
          ),
          title = list(type = "string", minLength = 1L),
          description = list(type = "string", minLength = 1L),
          inputSchema = list(
            type = "object",
            description = "JSON Schema draft 2020-12 for the tool's arguments"
          )
        )
      )
    ))
  )
}

render <- function(specs) {
  json <- jsonlite::toJSON(
    build_schema(specs),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  paste0(as.character(json), "\n")
}

text <- render(mcp_specs())
target <- file.path(root, "inst", "mcp", "schema.json")
args <- commandArgs(trailingOnly = TRUE)

if ("--check" %in% args) {
  if (!file.exists(target)) stop("missing ", target, "; run with --write", call. = FALSE)
  current <- readChar(target, file.info(target)$size, useBytes = TRUE)
  if (!identical(current, text)) {
    stop(
      "inst/mcp/schema.json has drifted from mcp_tool_specs(); run ",
      "Rscript dev/mcp-schema.R --write and review the diff",
      call. = FALSE
    )
  }
  cat("inst/mcp/schema.json matches mcp_tool_specs()\n")
} else if ("--write" %in% args) {
  # cat(), not writeLines(): `text` already ends in a newline and writeLines()
  # would add a second one, so --write then --check would report drift.
  cat(text, file = target, sep = "")
  cat("wrote ", target, " (", nchar(text), " bytes)\n", sep = "")
} else {
  cat(text)
}
