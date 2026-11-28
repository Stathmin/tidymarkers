# Roadmap T8: `tools/call` validates its arguments against the same
# `inputSchema` that `tools/list` publishes. The shared `mcp_env()` /
# `mcp_call()` helpers live in helper-mcp.R.

.mcp_sample_value <- function(schema) {
  if (!is.null(schema$enum)) {
    return(schema$enum[[1L]])
  }
  switch(schema$type,
    string = "ACGT",
    number = 1,
    integer = 1L,
    boolean = TRUE,
    array = if (length(schema$items)) list(.mcp_sample_value(schema$items)) else list(),
    object = {
      props <- if (is.null(schema$properties)) list() else schema$properties
      out <- lapply(props, .mcp_sample_value)
      names(out) <- names(props)
      out
    },
    stop("no sample value for schema type ", schema$type)
  )
}

.mcp_valid_args <- function(schema) {
  required <- unlist(schema$required, use.names = FALSE)
  props <- schema$properties[required]
  out <- lapply(props, .mcp_sample_value)
  names(out) <- names(props)
  out
}

.mcp_wrong_value <- function(schema) {
  switch(schema$type,
    string = 42L,
    number = "not-a-number",
    integer = "not-an-integer",
    boolean = "not-a-boolean",
    array = "not-an-array",
    object = "not-an-object",
    list()
  )
}

test_that("every tool accepts a minimal valid call and rejects a wrong argument", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()
  specs <- env$mcp_tool_specs()
  expect_gt(length(specs), 30L)

  for (nm in names(specs)) {
    schema <- specs[[nm]]$inputSchema
    expect_null(env$mcp_validate_args(.mcp_valid_args(schema), schema), info = nm)

    props <- names(schema$properties)
    if (length(props)) {
      args <- .mcp_valid_args(schema)
      args[[props[[1L]]]] <- .mcp_wrong_value(schema$properties[[props[[1L]]]])
      problem <- env$mcp_validate_args(args, schema)
      expect_equal(problem$path, props[[1L]], info = nm)
    } else {
      problem <- env$mcp_validate_args(list(bogus = 1L), schema)
      expect_equal(problem$path, "bogus", info = nm)
    }
  }
})

test_that("tools/call reports invalid params as JSON-RPC -32602", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  wrong_type <- .mcp_call(env, "primer_tm", list(seq = 42L))
  expect_equal(wrong_type$error$code, -32602L)
  expect_equal(wrong_type$error$data$tool, "primer_tm")
  expect_equal(wrong_type$error$data$path, "seq")
  expect_equal(wrong_type$error$data$expected, "string")
  expect_equal(wrong_type$error$data$got, 42L)

  missing <- .mcp_call(env, "primer_tm", list())
  expect_equal(missing$error$code, -32602L)
  expect_equal(missing$error$data$path, "seq")
  expect_equal(missing$error$data$expected, "required argument")

  extra <- .mcp_call(env, "marker_schemes", list(bogus = 1L))
  expect_equal(extra$error$code, -32602L)
  expect_equal(extra$error$data$path, "bogus")

  ok <- .mcp_call(env, "primer_tm", list(seq = "ACGTACGT"))
  expect_null(ok$error)
  expect_false(isTRUE(ok$result$isError))
})

test_that("the validator reads the JSON Schema keywords it publishes", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()
  schema <- list(
    type = "object",
    additionalProperties = FALSE,
    required = list("n", "motif"),
    properties = list(
      n = list(type = "integer", minimum = 1L, maximum = 10L),
      motif = list(type = "string", pattern = "^[ACGT]+$"),
      tags = list(type = "array", items = list(type = "string", enum = c("a", "b"))),
      obj = list(type = "object", properties = list(k = list(type = "number")))
    )
  )
  valid <- function() list(n = 3L, motif = "ACGT")

  expect_null(env$mcp_validate_args(valid(), schema))

  too_small <- valid()
  too_small$n <- 0L
  expect_equal(env$mcp_validate_args(too_small, schema)$path, "n")
  expect_equal(env$mcp_validate_args(too_small, schema)$expected, ">= 1")

  too_big <- valid()
  too_big$n <- 11L
  expect_equal(env$mcp_validate_args(too_big, schema)$expected, "<= 10")

  fractional <- valid()
  fractional$n <- 1.5
  expect_equal(env$mcp_validate_args(fractional, schema)$expected, "integer")

  bad_motif <- valid()
  bad_motif$motif <- "ACGU"
  expect_equal(env$mcp_validate_args(bad_motif, schema)$expected, "pattern ^[ACGT]+$")

  bad_item <- valid()
  bad_item$tags <- list("c")
  expect_equal(env$mcp_validate_args(bad_item, schema)$path, "tags[1]")

  bad_nested <- valid()
  bad_nested$obj <- list(k = "x")
  expect_equal(env$mcp_validate_args(bad_nested, schema)$path, "obj.k")

  unknown <- valid()
  unknown$extra <- 1L
  expect_equal(env$mcp_validate_args(unknown, schema)$path, "extra")

  expect_equal(env$mcp_validate_args("not-an-object", schema)$path, "arguments")
})

# Five tools a client reaches in the documented J1-J7 paths: thermo, the scheme
# registry, batch design, screening and the one mutating tool.  These pin the
# published file to what `tools/list` actually returns, and the structured
# error to what `tools/call` actually sends.
.mcp_pinned_tools <- c(
  "primer_tm", "marker_schemes", "design_markers",
  "screen_snp_targets", "marker_technology"
)

test_that("five tools return the published schema from tools/list", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  path <- system.file("mcp", "schema.json", package = "tidymarkers")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "mcp", "schema.json")
  }
  published <- jsonlite::fromJSON(path, simplifyVector = FALSE)$`$defs`$inputSchemas
  expect_true(all(.mcp_pinned_tools %in% names(published)))

  listed <- env$mcp_dispatch(list(jsonrpc = "2.0", id = 1L, method = "tools/list"))
  runtime <- listed$result$tools
  names(runtime) <- vapply(runtime, function(tool) tool$name, character(1))

  render <- function(x) {
    as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"))
  }
  for (nm in .mcp_pinned_tools) {
    expect_identical(render(runtime[[nm]]$inputSchema), render(published[[nm]]), info = nm)
    expect_identical(runtime[[nm]]$inputSchema, env$mcp_tool_specs()[[nm]]$inputSchema, info = nm)
  }
})

test_that("five tools reject bad arguments with the four structured fields", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  for (nm in .mcp_pinned_tools) {
    error <- .mcp_call(env, nm, list(bogus = 1L))$error
    expect_equal(error$code, -32602L, info = nm)
    expect_identical(error$data$tool, nm)
    expect_true(nzchar(error$data$path), info = nm)
    expect_true(nzchar(error$data$expected), info = nm)
    expect_true("got" %in% names(error$data), info = nm)
  }
})
