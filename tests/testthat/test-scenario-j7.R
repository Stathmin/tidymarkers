# J7: the MCP surface is a scenario of its own.  It proves the protocol
# (tools/list, tools/call, a structured protocol error) and that every tool a
# scenario names in docs/mcp-parity.csv is actually served.

test_that("J7: tools/list serves every scenario tool and tools/call answers", {
  skip_if_not_installed("jsonlite")
  scenario_timer("J7")
  env <- .mcp_env()
  listed <- scenario_step(
    "j7", "mcp", "tools_list",
    env$mcp_dispatch(list(jsonrpc = "2.0", id = 1L, method = "tools/list"))
  )
  expect_identical(listed$jsonrpc, "2.0")
  expect_identical(listed$id, 1L)
  served <- vapply(listed$result$tools, `[[`, character(1), "name")
  expect_length(served, length(env$mcp_tool_specs()))

  parity <- utils::read.csv(
    test_path("..", "..", "docs", "mcp-parity.csv"),
    stringsAsFactors = FALSE
  )
  wanted <- unique(unlist(strsplit(
    parity$tools[!parity$human_only & nzchar(parity$tools)], ";",
    fixed = TRUE
  )))
  expect_true(length(wanted) > 0L)
  expect_setequal(intersect(wanted, served), wanted)

  # a tool with no engine and no write answers with typed content
  criteria <- scenario_step(
    "j7", "mcp", "criteria",
    .mcp_call(env, "marker_criteria", list(scheme = "snp"))
  )
  expect_null(criteria$error)
  expect_true(length(criteria$result$content) >= 1L)
  expect_false(isTRUE(criteria$result$isError))

  # an unknown tool is a protocol error, not a handler crash
  bad <- env$mcp_dispatch(list(
    jsonrpc = "2.0", id = 2L, method = "tools/call",
    params = list(name = "not_a_tool", arguments = list())
  ))
  expect_identical(bad$error$code, -32602L)
  expect_identical(bad$error$data$tool, "not_a_tool")
})
