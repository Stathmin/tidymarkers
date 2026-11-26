# MCP tool layer (roadmap YP2/T8): every registered tool is wired and callable.
# The shared `mcp_env()` / `mcp_call()` helpers live in helper-mcp.R.

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

test_that("every tool advertises MCP annotations", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()
  specs <- env$mcp_tool_specs()
  hints <- c("readOnlyHint", "destructiveHint", "idempotentHint", "openWorldHint")

  for (nm in names(specs)) {
    annotations <- specs[[nm]]$annotations
    expect_true(all(hints %in% names(annotations)), info = nm)
    for (hint in hints) {
      expect_type(annotations[[hint]], "logical")
    }
    expect_false(annotations$destructiveHint, info = nm)
    expect_true(annotations$idempotentHint, info = nm)
    expect_false(annotations$openWorldHint, info = nm)
  }

  # Only the registry writer is not read-only; everything else has no side
  # effect beyond reading its mounted inputs.
  expect_false(specs$marker_technology$annotations$readOnlyHint)
  others <- specs[setdiff(names(specs), "marker_technology")]
  expect_true(all(vapply(others, function(s) isTRUE(s$annotations$readOnlyHint), logical(1))))

  listed <- env$mcp_dispatch(list(jsonrpc = "2.0", id = 3L, method = "tools/list"))
  expect_true(all(vapply(listed$result$tools, function(tool) {
    identical(tool$annotations, specs[[tool$name]]$annotations)
  }, logical(1))))
})

test_that("every registered tool is wired and returns content", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()
  specs <- env$mcp_tool_specs()

  # `marker_technology` registers into the package-level store; drop the smoke
  # technology so it cannot leak into another file's scheme list (test-plan).
  # `marker_schemes()` answers from the session memo, so removing the entry is
  # not enough - the memo has to be dropped too, exactly as `marker_technology()`
  # does when it changes the registry.
  on.exit(
    {
      if (exists("mcp_smoke", envir = .marker_tech_store, inherits = FALSE)) {
        rm("mcp_smoke", envir = .marker_tech_store)
      }
      session_memo_clear()
    },
    add = TRUE
  )

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
    design_chemistries = list(),
    design_scheme = list(
      genome_path = genome, chrom = "Chr1", pos = 81L,
      ref = "A", alt = "G", scheme = "caps", params = list(flank = 40L)
    ),
    design_markers = list(
      genome_path = genome, variants_path = vcf, scheme = "snp", flank = 40L
    )
  )
  expect_setequal(names(specs), names(arguments))

  # pcr_efficiency() needs DECIPHER's UNAFold call, so the tool reports an
  # error instead of a number where OligoArrayAux is absent.  The rest of the
  # registry is engine-free and still gets smoke-tested there.
  if (!has_oligoarrayaux()) {
    arguments <- arguments[setdiff(names(arguments), "pcr_efficiency")]
  }

  for (name in names(arguments)) {
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

  # Chemistries the package registers: the enum comes from marker_schemes(), not a
  # hand-written vector, so the schema cannot drift from the registry.
  marker_schemes <- env$marker_schemes

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
    # An enum (on the property or on an array's items) holds no duplicates.
    for (p in props) {
      prop <- schema$properties[[p]]
      enum <- if (is.null(prop$enum)) prop$items$enum else prop$enum
      if (!is.null(enum)) {
        expect_false(anyDuplicated(enum) > 0L, info = paste(nm, p))
      }
    }
  }

  # Every scheme-carrying property names real schemes.  The set comes from
  # marker_schemes(), so a hand-written enum cannot drift from the registry.
  scheme_props <- list(
    marker_criteria = "scheme", scheme_conditions = "scheme",
    marker_decision_text = "scheme", design_scheme = "scheme",
    design_markers = "scheme", design_marker = "schemes",
    rescue_strategy = "schemes"
  )
  for (tool in names(scheme_props)) {
    prop <- specs[[tool]]$inputSchema$properties[[scheme_props[[tool]]]]
    enum <- if (is.null(prop$enum)) prop$items$enum else prop$enum
    expect_true(all(enum %in% marker_schemes()$scheme), info = tool)
  }

  # A design enum advertises every per-locus designer: `rescue_strategy()$scheme`
  # drops `primer` (C13) and `duplex` (two loci), neither of which the
  # single-locus design tools can take.
  designable <- rescue_strategy()$scheme
  expect_setequal(specs$design_scheme$inputSchema$properties$scheme$enum, designable)
  expect_setequal(
    specs$design_marker$inputSchema$properties$schemes$items$enum,
    designable
  )
  expect_false("primer" %in% designable)
  expect_false("duplex" %in% designable)
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

test_that("a missing file is a readable tool error", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  missing_genome <- .mcp_call(env, "genome_info", list(
    genome_path = tempfile(fileext = ".fa")
  ))
  expect_true(isTRUE(missing_genome$result$isError))
  expect_true(nzchar(missing_genome$result$content[[1L]]$text))
})

test_that("protocol edges fail readably, not silently", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()

  # A notification (no id) has no response.
  expect_null(env$mcp_dispatch(list(jsonrpc = "2.0", method = "notifications/initialized")))

  # An unknown method is a JSON-RPC error.
  bad <- env$mcp_dispatch(list(jsonrpc = "2.0", id = 9L, method = "no/such"))
  expect_equal(bad$error$code, -32601L)

  # An unknown tool is invalid params, not a dropped call.
  unknown <- .mcp_call(env, "no_such_tool", list())
  expect_equal(unknown$error$code, -32602L)
  expect_equal(unknown$error$data$tool, "no_such_tool")
})

test_that("inst/mcp/schema.json matches the tool registry", {
  skip_if_not_installed("jsonlite")
  env <- .mcp_env()
  specs <- env$mcp_tool_specs()

  path <- system.file("mcp", "schema.json", package = "tidymarkers")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "mcp", "schema.json")
  }
  expect_true(file.exists(path), info = path)

  text <- readChar(path, file.info(path)$size, useBytes = TRUE)
  expect_true(jsonlite::validate(text))
  schema <- jsonlite::fromJSON(text, simplifyVector = FALSE)

  expect_equal(schema$`$schema`, "https://json-schema.org/draft/2020-12/schema")
  expect_equal(schema$type, "object")
  expect_false(isTRUE(schema$additionalProperties))

  tools <- schema$properties$tools
  # The count and the name enum are the two places a new or renamed tool shows
  # up; regenerate with `Rscript dev/mcp-schema.R --write` rather than editing.
  expect_equal(tools$minItems, length(specs))
  expect_equal(tools$maxItems, length(specs))
  expect_setequal(unlist(tools$items$properties$name$enum), names(specs))
  expect_equal(
    unlist(tools$items$required, use.names = FALSE),
    c("name", "title", "description", "inputSchema", "annotations")
  )
  expect_equal(
    unlist(tools$items$properties$annotations$required, use.names = FALSE),
    c("readOnlyHint", "destructiveHint", "idempotentHint", "openWorldHint")
  )

  # `$defs.inputSchemas` publishes the full argument schema of every tool, not
  # just `type: object`; compare it with the registry as JSON so integer vs
  # numeric and vector vs list cannot blur the difference.
  defs <- schema[["$defs"]]$inputSchemas
  expect_setequal(names(defs), names(specs))
  for (nm in names(specs)) {
    expect_equal(
      as.character(jsonlite::toJSON(defs[[nm]], auto_unbox = TRUE, null = "null")),
      as.character(jsonlite::toJSON(specs[[nm]]$inputSchema, auto_unbox = TRUE, null = "null")),
      info = nm
    )
  }

  for (nm in names(specs)) {
    input <- specs[[nm]]$inputSchema
    required <- unlist(input$required, use.names = FALSE)
    expect_true(all(required %in% names(input$properties)), info = nm)
  }
})
