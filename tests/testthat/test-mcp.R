# MCP tool layer (roadmap YP2): every registered tool is wired and callable.
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

test_that("initialize and tools/list speak JSON-RPC", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  init <- env$mcp_dispatch(list(jsonrpc = "2.0", id = 1L, method = "initialize"))
  expect_equal(init$result$protocolVersion, "2024-11-05")
  expect_equal(init$result$serverInfo$name, "tidymarkers")
  expect_type(init$result$serverInfo$version, "character")

  listed <- env$mcp_dispatch(list(jsonrpc = "2.0", id = 2L, method = "tools/list"))
  expect_length(listed$result$tools, length(env$mcp_tool_specs()))
  expect_true(all(
    c("name", "title", "description", "inputSchema") %in% names(listed$result$tools[[1L]])
  ))
})

test_that("every registered tool is wired and returns content", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()
  specs <- env$mcp_tool_specs()

  # A mounted genome and VCF: the design tools take paths, not in-memory objects.
  genome <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">Chr1", strrep("ACGT", 60L)), genome)
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "Chr1\t81\t.\tA\tG\t.\t.\t."
  ), vcf)
  primer <- "ACGTACGTACGTACGTACGT"

  # One representative call per tool.  A new tool must add its call here, so
  # "every tool" stays true as the registry grows.
  arguments <- list(
    primer_tm = list(seq = primer),
    oligo_dg = list(seq1 = primer),
    terminal_dg = list(primer = primer),
    marker_schemes = list(),
    marker_criteria = list(scheme = "snp"),
    marker_technologies = list(),
    marker_technology = list(name = "mcp_smoke", base = "snp", overwrite = TRUE),
    marker_priors = list(),
    kasp_mix_recipe = list(n_reactions = 48L),
    specificity_policy = list(scheme = "pcr"),
    genome_info = list(genome_path = genome),
    read_variants = list(variants_path = vcf),
    in_silico_pcr = list(f_primer = primer, r_primer = primer, genome_path = genome),
    design_marker = list(
      genome_path = genome, chrom = "Chr1", pos = 81L,
      ref = "A", alt = "G", schemes = list("snp")
    ),
    variant_type = list(ref = "A", alt = "G"),
    marker_site_class = list(ref = "A", alts = list("G", "T")),
    mask_positions = list(seq = "ACGTACGT", positions = list(2L, 5L)),
    asq_limits = list(),
    reaction_conditions = list(),
    scheme_conditions = list(scheme = "pcr"),
    primer_thermo = list(primers = list(primer)),
    primer_set_compatibility = list(primers = list(primer, "GCATCGATTAGCGATCGATCG")),
    dimer_thermo = list(primers = list(primer, "GCATCGATTAGCGATCGATCG")),
    hairpin_dg_thermo = list(seq = primer),
    fold_amplicon = list(seq = "ACGTACGTACGTACGTACGTACGTACGTACGT", kind = "DNA"),
    pcr_efficiency = list(
      primers = list("GAAGGTGACCAAGTTCATGCT"),
      targets = list("AGCATGAACTTGGTCACCTTC")
    ),
    wobble_mismatch = list(on_allele = "G", off_allele = "T"),
    find_palindromes = list(seq = "ACGTACGTAAAAACGTACGT"),
    has_palindrome = list(seq = "ACGTACGT"),
    find_caps_enzymes = list(
      genome_path = genome, chrom = "Chr1", pos = 81L,
      ref = "A", alt = "G", flank = 40L, min_gap = 20L
    ),
    rescue_strategy = list(),
    marker_decision_text = list(scheme = "snp"),
    screen_snp_targets = list(variants_path = vcf, genome_path = genome),
    engine_status = list(),
    engine_capabilities = list(engine = "thermo"),
    design_scheme = list(
      genome_path = genome, chrom = "Chr1", pos = 81L,
      ref = "A", alt = "G", scheme = "caps", params = list(flank = 40L)
    ),
    design_markers = list(
      genome_path = genome, variants_path = vcf, scheme = "snp", flank = 40L
    )
  )
  expect_setequal(names(specs), names(arguments))

  for (name in names(specs)) {
    response <- .mcp_call(env, name, arguments[[name]])
    expect_false(isTRUE(response$result$isError), info = name)
    expect_equal(response$result$content[[1L]]$type, "text", info = name)
    expect_true(nzchar(response$result$content[[1L]]$text), info = name)
  }
})

test_that("every inputSchema is well formed", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()
  specs <- env$mcp_tool_specs()
  expect_gt(length(specs), 30L)

  for (nm in names(specs)) {
    spec <- specs[[nm]]
    expect_true(nzchar(spec$title), info = nm)
    expect_true(nzchar(spec$description), info = nm)
    expect_type(spec$handler, "closure")

    schema <- spec$inputSchema
    expect_equal(schema$type, "object", info = nm)
    props <- names(schema$properties)
    required <- unlist(schema$required, use.names = FALSE)
    # a required name that is not declared is a schema the client cannot satisfy
    expect_true(all(required %in% props), info = nm)
    for (p in props) {
      expect_false(is.null(schema$properties[[p]]$type), info = paste(nm, p))
    }
  }
})

test_that("genome tools answer from the mounted files", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  genome <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">Chr1", strrep("ACGT", 60L)), genome)
  primer <- "ACGTACGTACGTACGTACGT"

  info <- .mcp_call(env, "genome_info", list(genome_path = genome))
  expect_false(isTRUE(info$result$isError))
  expect_equal(info$result$structuredContent$chrom, "Chr1")
  expect_equal(info$result$structuredContent$length, 240L)

  # An exact forward/reverse pair in a repeated genome still resolves to products.
  pcr <- .mcp_call(env, "in_silico_pcr", list(
    f_primer = primer, r_primer = primer, genome_path = genome
  ))
  expect_false(isTRUE(pcr$result$isError))
  expect_gt(nrow(pcr$result$structuredContent), 0L)
})

test_that("a missing path is a readable tool error", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  missing_genome <- .mcp_call(env, "genome_info", list())
  expect_true(isTRUE(missing_genome$result$isError))
  expect_match(missing_genome$result$content[[1L]]$text, "genome_path")
})

test_that("protocol edges fail readably, not silently", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  # A notification (no id) has no response.
  expect_null(env$mcp_dispatch(list(jsonrpc = "2.0", method = "notifications/initialized")))

  # An unknown method is a JSON-RPC error.
  bad <- env$mcp_dispatch(list(jsonrpc = "2.0", id = 9L, method = "no/such"))
  expect_equal(bad$error$code, -32601L)

  # An unknown tool is a tool error, not a dropped call.
  unknown <- .mcp_call(env, "no_such_tool", list())
  expect_true(isTRUE(unknown$result$isError))
})
