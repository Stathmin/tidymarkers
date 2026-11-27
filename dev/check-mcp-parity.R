# =============================================================================
#  dev/check-mcp-parity.R  --  the MCP contract gate (ADR 0011 section 3)
# =============================================================================
#  One static gate over the machine interface.  It reads files; it does not
#  start a server.  Runtime behaviour stays in test-mcp.R and
#  test-mcp-validation.R.
#
#    1. mcp_tool_specs() names == schema.json $defs.inputSchemas, each published
#       inputSchema equals the runtime one as JSON, and every spec has a handler
#    2. every schema `required` argument is read by the handler (through
#       mcp_arg/mcp_int/mcp_chr/mcp_named_chr/mcp_num, following mcp_* helpers)
#       and every argument the handler reads is declared in `properties`
#    3. every schema keyword is one the validator enforces; a top-level
#       `properties` object is closed; enum/minimum/maximum/pattern are
#       well-formed and match their declared type
#    4. annotations agree with MCP_MUTATING_TOOLS (only marker_technology is
#       non-read-only)
#    5. every package function a handler calls is an export (or an S3 method)
#       with a row in docs/user-api-surface.csv
#    6. every J1-J7 scenario has a tool path in docs/mcp-parity.csv, or is
#       listed there as human_only with a reason
#
#      Rscript dev/check-mcp-parity.R
#
#  Exit 1 on any failure.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

failures <- character()
source(file.path(root, "dev", "check-reporter.R"))

tools_path <- file.path(root, "inst", "mcp", "tools.R")
schema_path <- file.path(root, "inst", "mcp", "schema.json")
api_path <- file.path(root, "docs", "user-api-surface.csv")
parity_path <- file.path(root, "docs", "mcp-parity.csv")

if (!file.exists(tools_path)) {
  note("missing ", tools_path)
  stop("cannot read the MCP registry")
}

pkgload::load_all(root, quiet = TRUE, export_all = FALSE)
mcp_env <- new.env(parent = asNamespace("tidymarkers"))
sys.source(tools_path, envir = mcp_env)
specs <- mcp_env$mcp_tool_specs()
tool_names <- names(specs)
ok("registry: ", length(tool_names), " tool(s)")

# -- 1. registry vs published schema ----------------------------------------

schema <- jsonlite::fromJSON(schema_path, simplifyVector = FALSE)
published <- names(schema$`$defs`$inputSchemas)
if (!setequal(tool_names, published)) {
  for (nm in setdiff(tool_names, published)) note("tool not in schema.json: ", nm)
  for (nm in setdiff(published, tool_names)) note("schema.json tool not in registry: ", nm)
} else {
  ok("schema.json: ", length(published), " inputSchema(s) match the registry")
}

# The file is a contract artefact, not a runtime source (ADR 0007 section 4),
# so the two definitions have to agree byte for byte.  Compare as JSON: integer
# vs numeric and vector vs list must not blur.
published_schemas <- schema$`$defs`$inputSchemas
schema_drift <- 0L
for (nm in intersect(tool_names, published)) {
  render <- function(x) {
    as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"))
  }
  if (identical(render(specs[[nm]]$inputSchema), render(published_schemas[[nm]]))) {
    next
  }
  schema_drift <- schema_drift + 1L
  keys <- union(names(specs[[nm]]$inputSchema), names(published_schemas[[nm]]))
  bad <- keys[!vapply(keys, function(k) {
    identical(render(specs[[nm]]$inputSchema[[k]]), render(published_schemas[[nm]][[k]]))
  }, logical(1))]
  note(paste0(
    nm, ": inputSchema in schema.json differs from the runtime registry at ",
    if (length(bad)) paste(bad, collapse = ", ") else "a nested value"
  ))
}
if (!schema_drift) {
  ok("schema.json: every published inputSchema equals the runtime one")
}
missing_handler <- tool_names[!vapply(specs, function(s) is.function(s$handler), logical(1))]
for (nm in missing_handler) note("tool has no handler: ", nm)

# -- handler argument extraction --------------------------------------------
# The handlers read `args` only through the mcp_* accessors, sometimes inside
# an mcp_* helper (mcp_genome, mcp_oligo, mcp_snp).  Collect the literal names
# from the handler body and follow the helpers transitively.

ACCESSORS <- c("mcp_arg", "mcp_int", "mcp_chr", "mcp_named_chr", "mcp_num")
exprs <- parse(tools_path)
helpers <- list()
for (e in exprs) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
    is.symbol(e[[2L]]) && is.call(e[[3L]]) &&
    identical(e[[3L]][[1L]], as.name("function"))) {
    helpers[[as.character(e[[2L]])]] <- e[[3L]]
  }
}

collect_calls <- function(node, args, funcs, seen) {
  if (is.call(node)) {
    fn <- node[[1L]]
    if (is.symbol(fn)) {
      nm <- as.character(fn)
      if (nm %in% ACCESSORS && length(node) >= 3L && is.character(node[[3L]])) {
        args <- c(args, node[[3L]])
      } else if (nm %in% names(helpers) && !nm %in% seen) {
        formals <- as.list(helpers[[nm]][[2L]])
        defaults <- unlist(formals[vapply(formals, is.character, logical(1))], use.names = FALSE)
        got <- collect_calls(helpers[[nm]][[3L]], character(), character(), c(seen, nm))
        args <- c(args, defaults, got$args)
        funcs <- c(funcs, got$funcs)
      } else {
        funcs <- c(funcs, nm)
      }
    }
    for (i in seq_along(node)[-1L]) {
      got <- collect_calls(node[[i]], character(), character(), seen)
      args <- c(args, got$args)
      funcs <- c(funcs, got$funcs)
    }
  } else if (is.symbol(node)) {
    funcs <- c(funcs, as.character(node))
  } else if (is.pairlist(node) || is.list(node)) {
    for (i in seq_along(node)) {
      got <- collect_calls(node[[i]], character(), character(), seen)
      args <- c(args, got$args)
      funcs <- c(funcs, got$funcs)
    }
  }
  list(args = unique(args), funcs = unique(funcs))
}

# -- 2. declared vs read arguments -----------------------------------------

supported_keywords <- c(
  "type", "description", "title", "properties", "required",
  "additionalProperties", "items", "enum", "minimum", "maximum", "pattern"
)

check_schema <- function(node, where, tool) {
  if (!is.list(node)) {
    note(tool, ": ", where, " is not a schema object")
    return(invisible())
  }
  unknown <- setdiff(names(node), supported_keywords)
  if (length(unknown)) {
    note(tool, ": ", where, " uses keyword(s) the validator ignores: ", paste(unknown, collapse = ", "))
  }
  type <- node$type
  if (is.null(type)) {
    note(tool, ": ", where, " has no type")
  } else if (!type %in% c("string", "number", "integer", "boolean", "array", "object")) {
    note(tool, ": ", where, " has unsupported type ", type)
  }
  props <- node$properties
  req <- unlist(node$required, use.names = FALSE)
  if (length(req) && is.null(props)) {
    note(tool, ": ", where, " declares required without properties")
  } else if (length(req) && !all(req %in% names(props))) {
    for (nm in setdiff(req, names(props))) note(tool, ": ", where, " requires undeclared ", nm)
  }
  if (identical(type, "array") && is.null(node$items)) {
    note(tool, ": ", where, " is an array with no items")
  }
  if (!is.null(node$enum)) {
    vals <- unlist(node$enum, use.names = FALSE)
    if (!length(vals)) note(tool, ": ", where, " has an empty enum")
    if (identical(type, "string") && !is.character(vals)) {
      note(tool, ": ", where, " enum members do not match type string")
    }
    if (type %in% c("number", "integer") && !is.numeric(vals)) {
      note(tool, ": ", where, " enum members do not match type ", type)
    }
  }
  if (!is.null(node$minimum) && !is.numeric(node$minimum)) {
    note(tool, ": ", where, " minimum is not numeric")
  }
  if (!is.null(node$maximum) && !is.numeric(node$maximum)) {
    note(tool, ": ", where, " maximum is not numeric")
  }
  if (!is.null(node$minimum) && !is.null(node$maximum) && node$minimum > node$maximum) {
    note(tool, ": ", where, " minimum exceeds maximum")
  }
  if (!is.null(node$pattern)) {
    bad <- tryCatch(
      {
        grepl(node$pattern, "")
        FALSE
      },
      error = function(e) TRUE
    )
    if (bad) note(tool, ": ", where, " has an invalid pattern")
  }
  if (length(props)) {
    for (nm in names(props)) check_schema(props[[nm]], paste0(where, ".", nm), tool)
  }
  if (!is.null(node$items)) check_schema(node$items, paste0(where, "[]"), tool)
  invisible()
}

for (nm in tool_names) {
  spec <- specs[[nm]]
  schema_in <- spec$inputSchema
  declared <- names(schema_in$properties)
  required <- unlist(schema_in$required, use.names = FALSE)
  got <- collect_calls(body(spec$handler), character(), character(), character())
  read <- got$args
  for (a in setdiff(required, read)) {
    note(nm, ": required argument never read by the handler: ", a)
  }
  for (a in setdiff(read, declared)) {
    note(nm, ": handler reads undeclared argument: ", a)
  }
  for (a in setdiff(declared, read)) {
    note(nm, ": schema declares an argument the handler never reads: ", a)
  }
  if (!identical(schema_in$additionalProperties, FALSE)) {
    note(nm, ": inputSchema is not closed (additionalProperties must be FALSE)")
  }
  check_schema(schema_in, "inputSchema", nm)
  if (length(read)) ok(nm, ": ", length(read), " typed argument(s)")
}

# -- 4. annotations ---------------------------------------------------------

read_only <- mcp_env$MCP_ANNOTATIONS
mutating <- mcp_env$MCP_MUTATING_TOOLS
for (nm in tool_names) {
  want <- if (nm %in% mutating) mcp_env$MCP_MUTATING_ANNOTATIONS else read_only
  got <- specs[[nm]]$annotations
  if (!identical(got, want)) {
    note(nm, ": annotations disagree with MCP_MUTATING_TOOLS")
  }
}
ok("annotations: ", length(tool_names) - length(mutating), " read-only, ",
  length(mutating), " mutating")

# -- 5. every handler call is a documented export --------------------------

api <- utils::read.csv(api_path, stringsAsFactors = FALSE, check.names = FALSE)
ns <- asNamespace("tidymarkers")
exports <- getNamespaceExports(ns)
documented <- api[["function"]]
internals <- setdiff(ls(ns, all.names = TRUE), exports)

for (nm in tool_names) {
  funcs <- collect_calls(body(specs[[nm]]$handler), character(), character(), character())$funcs
  funcs <- intersect(funcs, c(documented, internals))
  if (!length(funcs)) {
    note(nm, ": handler calls no documented package function")
    next
  }
  for (f in funcs) {
    row <- api[api[["function"]] == f, , drop = FALSE]
    if (!nrow(row)) {
      note(nm, ": handler calls `", f, "` with no row in user-api-surface.csv")
    } else if (!row$kind[[1L]] %in% c("export", "s3method")) {
      note(nm, ": handler calls `", f, "` which is not an export")
    } else if (!f %in% exports && !row$kind[[1L]] == "s3method") {
      note(nm, ": handler calls `", f, "` which is absent from NAMESPACE")
    }
  }
}
ok("handlers: every called package function is a documented export")

# -- 6. scenario coverage ---------------------------------------------------

parity <- utils::read.csv(parity_path, stringsAsFactors = FALSE)
scenarios <- paste0("J", 1:7)
if (!setequal(parity$scenario, scenarios)) {
  for (s in setdiff(scenarios, parity$scenario)) note("no MCP parity row for scenario ", s)
  for (s in setdiff(parity$scenario, scenarios)) note("unknown scenario in mcp-parity.csv: ", s)
}
for (i in seq_len(nrow(parity))) {
  s <- parity$scenario[[i]]
  human <- isTRUE(as.logical(parity$human_only[[i]]))
  reason <- parity$reason[[i]]
  tools <- if (is.na(parity$tools[[i]]) || !nzchar(parity$tools[[i]])) {
    character()
  } else {
    strsplit(parity$tools[[i]], ";", fixed = TRUE)[[1L]]
  }
  if (human) {
    if (is.na(reason) || !nzchar(reason)) note(s, ": human_only without a reason")
    if (length(tools)) note(s, ": human_only but lists tools")
  } else {
    if (!length(tools)) note(s, ": no tool path and not human_only")
    for (t in setdiff(tools, tool_names)) note(s, ": unknown tool ", t)
    if (!is.na(reason) && nzchar(reason)) note(s, ": has a tool path but a human_only reason")
  }
}
ok(
  "scenarios: ", sum(!as.logical(parity$human_only)), " with a tool path, ",
  sum(as.logical(parity$human_only)), " human_only"
)

# -- summary ----------------------------------------------------------------

if (length(failures)) {
  message("\nFAIL mcp-parity: ", length(failures), " problem(s)")
  for (f in failures) message("  - ", f)
  quit(status = 1L)
}
ok("mcp-parity")
