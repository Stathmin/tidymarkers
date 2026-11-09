# =============================================================================
#  tools.R  --  MCP tool registry over the exported tidymarkers API (YP2)
# =============================================================================
#  One spec per call: a name, a JSON-Schema input, and a handler returning an R
#  object.  `mcp_dispatch()` answers JSON-RPC 2.0 for the `initialize`,
#  `tools/list` and `tools/call` methods, so server.R stays a thin stdin/stdout
#  loop.  Design and thermodynamics calls take the paths they need: `genome_path`
#  (FASTA), `variants_path` (VCF) and `oligoarrayaux_dir` (an OligoArrayAux prefix
#  or `bin/`).  The paths are mounted by the caller; when a path is absent the
#  call falls back to the package search (demo cache, `.tools/`, `PATH`).
#
#  This file is not part of the package namespace.  Source it after loading the
#  package; the handlers resolve the exported functions from the search path.
# =============================================================================

# `%||%` is internal to the package.  Define it here so the file also runs when
# sourced on its own (the server loads only the installed, attached package).
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# -- inputs --------------------------------------------------------------

# A NULL from JSON's `null` is "argument absent", so every handler reads its
# arguments through these two accessors rather than `$`.
mcp_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value)) NULL else value
}

mcp_int <- function(args, name) {
  value <- mcp_arg(args, name)
  if (is.null(value)) NULL else as.integer(value)
}

# JSON has no vector type: an array arrives as a list, an object as a named list.
mcp_chr <- function(args, name) {
  value <- mcp_arg(args, name)
  if (is.null(value)) NULL else as.character(unlist(value, use.names = FALSE))
}

mcp_named_chr <- function(args, name) {
  value <- mcp_arg(args, name)
  if (is.null(value)) NULL else unlist(value)
}

# -- external inputs -----------------------------------------------------

# A design call reads its genome from `genome_path`; an optional `chrom_map`
# object (defline -> name) renames sequences.  The FASTA index is built on first
# use, so a plain file is enough.
mcp_genome <- function(args) {
  path <- mcp_arg(args, "genome_path")
  if (is.null(path)) {
    stop("`genome_path` is required: a FASTA file (plain or gzipped)")
  }
  marker_genome(path = path, chrom_map = mcp_named_chr(args, "chrom_map"))
}

# Point the thermodynamics engine at a mounted OligoArrayAux prefix for one
# call.  Without `oligoarrayaux_dir` the package search decides, so the same
# tool works on the host and in the container.
mcp_oligo <- function(args, code) {
  dir <- mcp_arg(args, "oligoarrayaux_dir")
  if (is.null(dir)) {
    return(force(code))
  }
  withr::with_options(list(tidymarkers.oligoarrayaux = dir), force(code))
}

# -- serialisation -------------------------------------------------------

# One JSON string for any handler result; tibbles serialise as row objects.
mcp_json <- function(x) {
  jsonlite::toJSON(
    x,
    auto_unbox = TRUE, dataframe = "rows", null = "null",
    na = "null", digits = NA, force = TRUE
  )
}

# -- the registry --------------------------------------------------------

# Every entry: title, description, inputSchema (JSON Schema draft 2020-12), and
# a handler taking the parsed `arguments` object.
mcp_tool_specs <- function() {
  list(
    primer_tm = list(
      title = "Primer melting temperature",
      description = paste(
        "Nearest-neighbour Tm of a DNA oligo under the package's reaction",
        "conditions. `salt_method` selects the salt correction; the design",
        "default is Schildkraut2010."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          seq = list(type = "string", description = "DNA sequence, 5'->3'"),
          salt_method = list(
            type = "string",
            enum = c(
              "Schildkraut2010", "SantaLucia1998-1", "SantaLucia1998-2",
              "Owczarzy2004", "Owczarzy2008"
            )
          ),
          dnac = list(
            type = "number",
            description = "strand concentration in nM (absent = 25)"
          )
        ),
        required = list("seq")
      ),
      handler = function(args) {
        primer_tm(
          mcp_arg(args, "seq"),
          salt_method = mcp_arg(args, "salt_method"),
          dnac = mcp_arg(args, "dnac")
        )
      }
    ),
    oligo_dg = list(
      title = "Oligo hybridisation free energy",
      description = paste(
        "dG of a duplex (or a single-strand fold when `seq2` is absent) at the",
        "given sodium molarity and temperature. Returns null when the",
        "OligoArrayAux engine is not installed."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          seq1 = list(type = "string", description = "first oligo, 5'->3'"),
          seq2 = list(type = "string", description = "second oligo; absent = self-fold"),
          na = list(type = "number", description = "sodium molarity (default 0.05)"),
          temp = list(type = "number", description = "temperature in Celsius (default 37)"),
          oligoarrayaux_dir = list(
            type = "string",
            description = "OligoArrayAux prefix or bin/ (absent = package search)"
          )
        ),
        required = list("seq1")
      ),
      handler = function(args) {
        mcp_oligo(args, oligo_dg(
          mcp_arg(args, "seq1"),
          seq2 = mcp_arg(args, "seq2"),
          na = mcp_arg(args, "na") %||% 0.05,
          temp = mcp_arg(args, "temp") %||% 37
        ))
      }
    ),
    terminal_dg = list(
      title = "Terminal stack free energy",
      description = paste(
        "Free energy of the terminal `n` base pairs of a primer; the measure a",
        "3'-mismatch penalty is read against."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          primer = list(type = "string", description = "primer, 5'->3'"),
          n = list(type = "integer", description = "number of terminal bases (default 5)"),
          oligoarrayaux_dir = list(
            type = "string",
            description = "OligoArrayAux prefix or bin/ (absent = package search)"
          )
        ),
        required = list("primer")
      ),
      handler = function(args) {
        mcp_oligo(args, terminal_dg(mcp_arg(args, "primer"), n = mcp_int(args, "n") %||% 5L))
      }
    ),
    marker_schemes = list(
      title = "Marker schemes",
      description = paste(
        "The registered marker schemes: family, discrimination, readout,",
        "zygosity, amplicon class and verdict tone."
      ),
      inputSchema = list(type = "object", properties = list()),
      handler = function(args) marker_schemes()
    ),
    marker_criteria = list(
      title = "Scheme criteria",
      description = paste(
        "The gate and tier criteria a scheme evaluates, with stage, role and",
        "id. Defaults to the SNP (KASP) scheme."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          scheme = list(type = "string", description = "scheme id (default \"snp\")")
        )
      ),
      handler = function(args) {
        marker_criteria(scheme = mcp_arg(args, "scheme") %||% "snp")
      }
    ),
    marker_technologies = list(
      title = "Defined technologies",
      description = "The marker technologies defined in this session.",
      inputSchema = list(type = "object", properties = list()),
      handler = function(args) marker_technologies()
    ),
    marker_technology = list(
      title = "Define a technology",
      description = paste(
        "Compose a named marker technology from a base scheme, parameter",
        "values and filter blocks. The result is registered for the session."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          name = list(type = "string", description = "technology name"),
          base = list(type = "string", description = "base scheme id (default \"snp\")"),
          label = list(type = "string", description = "human-readable label"),
          overwrite = list(type = "boolean", description = "replace an existing name")
        ),
        required = list("name")
      ),
      handler = function(args) {
        marker_technology(
          mcp_arg(args, "name"),
          base = mcp_arg(args, "base") %||% "snp",
          label = mcp_arg(args, "label"),
          overwrite = isTRUE(mcp_arg(args, "overwrite"))
        )
      }
    ),
    marker_priors = list(
      title = "Design priors",
      description = paste(
        "The package's design priors with their source: a design coordinate,",
        "a named reference, or a package default."
      ),
      inputSchema = list(type = "object", properties = list()),
      handler = function(args) marker_priors()
    ),
    kasp_mix_recipe = list(
      title = "KASP mix recipe",
      description = paste(
        "Per-reaction and total volumes of a KASP primer mix for the given",
        "number of reactions."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          n_reactions = list(type = "integer", description = "number of reactions (default 96)"),
          reaction_uL = list(type = "number", description = "reaction volume in uL (default 10)")
        )
      ),
      handler = function(args) {
        kasp_mix_recipe(
          n_reactions = mcp_int(args, "n_reactions") %||% 96L,
          reaction_uL = mcp_arg(args, "reaction_uL") %||% 10
        )
      }
    ),
    specificity_policy = list(
      title = "Specificity policy",
      description = paste(
        "The mismatch tolerance and scoring weights a specificity scheme",
        "applies when it scans for off-target products: one of pcr, rpa,",
        "crispr, paralog, or a marker alias such as kasp, arms, caps."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          scheme = list(type = "string", description = "specificity scheme (default \"pcr\")"),
          length = list(type = "integer", description = "primer length in nt")
        )
      ),
      handler = function(args) {
        specificity_policy(
          scheme = mcp_arg(args, "scheme") %||% "pcr",
          length = mcp_int(args, "length")
        )
      }
    ),
    genome_info = list(
      title = "Genome sequences",
      description = paste(
        "Sequence names and lengths of a mounted FASTA, as the design code sees",
        "them. `chrom_map` (an object of defline -> name) renames sequences and",
        "drops the rest."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          genome_path = list(type = "string", description = "FASTA file, plain or gzipped"),
          chrom_map = list(type = "object", description = "optional defline -> name map")
        ),
        required = list("genome_path")
      ),
      handler = function(args) {
        genome <- mcp_genome(args)
        chroms <- genome_chroms(genome)
        tibble::tibble(
          chrom = chroms,
          length = as.integer(vapply(chroms, function(ch) genome_length(genome, ch), integer(1)))
        )
      }
    ),
    read_variants = list(
      title = "Read a VCF",
      description = paste(
        "Parse CHROM/POS/REF/ALT from a VCF, one row per ALT allele, optionally",
        "restricted to a region. This is the input a design call expects."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          variants_path = list(type = "string", description = "VCF file"),
          chrom = list(type = "string", description = "region chromosome"),
          start = list(type = "integer", description = "region start (1-based)"),
          end = list(type = "integer", description = "region end (inclusive)")
        ),
        required = list("variants_path")
      ),
      handler = function(args) {
        path <- mcp_arg(args, "variants_path")
        if (is.null(path)) {
          stop("`variants_path` is required: a VCF file")
        }
        read_variants_vcf(
          path,
          chrom = mcp_arg(args, "chrom"),
          start = mcp_int(args, "start"),
          end = mcp_int(args, "end")
        )
      }
    ),
    in_silico_pcr = list(
      title = "In-silico PCR",
      description = paste(
        "Every primer pair product in a mounted genome, with mismatches and",
        "penalty. `scheme` selects a specificity policy (pcr, rpa, crispr,",
        "paralog); a CRISPR call also checks `pam`."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          f_primer = list(type = "string", description = "forward primer, 5'->3'"),
          r_primer = list(type = "string", description = "reverse primer, 5'->3'"),
          genome_path = list(type = "string", description = "FASTA file, plain or gzipped"),
          chrom_map = list(type = "object", description = "optional defline -> name map"),
          max_mismatch = list(type = "integer", description = "search radius (default: policy)"),
          scheme = list(type = "string", description = "specificity scheme"),
          pam = list(type = "string", description = "PAM/PFS for a CRISPR scheme")
        ),
        required = list("f_primer", "r_primer", "genome_path")
      ),
      handler = function(args) {
        in_silico_pcr(
          mcp_arg(args, "f_primer"),
          mcp_arg(args, "r_primer"),
          mcp_genome(args),
          max_mismatch = mcp_int(args, "max_mismatch"),
          scheme = mcp_arg(args, "scheme"),
          pam = mcp_arg(args, "pam")
        )
      }
    ),
    design_marker = list(
      title = "Design one marker",
      description = paste(
        "Try the scheme ladder on one locus and report which chemistry takes it.",
        "Returns the winning scheme (null when none) and the per-scheme verdicts;",
        "the marker object itself stays on the server."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          genome_path = list(type = "string", description = "FASTA file, plain or gzipped"),
          chrom_map = list(type = "object", description = "optional defline -> name map"),
          chrom = list(type = "string", description = "locus chromosome"),
          pos = list(type = "integer", description = "locus position (1-based)"),
          ref = list(type = "string", description = "reference allele"),
          alt = list(type = "string", description = "alternate allele"),
          schemes = list(type = "array", items = list(type = "string"), description = "scheme ids to try, in order")
        ),
        required = list("genome_path", "chrom", "pos", "ref", "alt")
      ),
      handler = function(args) {
        snp <- tibble::tibble(
          chrom = as.character(mcp_arg(args, "chrom")),
          pos = as.integer(mcp_arg(args, "pos")),
          ref = as.character(mcp_arg(args, "ref")),
          alt = as.character(mcp_arg(args, "alt"))
        )
        result <- design_auto_marker(mcp_genome(args), snp, schemes = mcp_chr(args, "schemes"))
        list(
          scheme = result$scheme,
          tried = result$tried[, c("scheme", "ok"), drop = FALSE]
        )
      }
    )
  )
}

# -- dispatch ------------------------------------------------------------

mcp_pkg_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("tidymarkers")),
    error = function(e) "dev"
  )
}

# A `tools/call` result: text for a plain client, structured content for one
# that can read it. A handler error is a tool error (`isError`), not a
# protocol error, so the client sees the message instead of a dropped call.
mcp_call <- function(params, specs) {
  name <- params$name
  args <- params$arguments %||% list()
  spec <- specs[[name]]
  if (is.null(spec)) {
    return(list(
      content = list(list(type = "text", text = paste0("unknown tool: ", name))),
      isError = TRUE
    ))
  }
  value <- tryCatch(
    spec$handler(args),
    error = function(e) structure(conditionMessage(e), class = "mcp_error")
  )
  if (inherits(value, "mcp_error")) {
    return(list(
      content = list(list(type = "text", text = as.character(value))),
      isError = TRUE
    ))
  }
  list(
    content = list(list(type = "text", text = as.character(mcp_json(value)))),
    structuredContent = value
  )
}

mcp_error_response <- function(id, code, message) {
  list(jsonrpc = "2.0", id = id, error = list(code = code, message = message))
}

# `request` is one decoded JSON-RPC message. A notification (no `id`) returns
# NULL, which the server loop does not write back.
mcp_dispatch <- function(request, specs = mcp_tool_specs()) {
  method <- request$method %||% ""
  id <- request[["id"]]
  result <- switch(method,
    initialize = list(
      protocolVersion = "2024-11-05",
      serverInfo = list(name = "tidymarkers", version = mcp_pkg_version()),
      capabilities = list(tools = list())
    ),
    ping = list(),
    `notifications/initialized` = NULL,
    `tools/list` = list(tools = lapply(names(specs), function(nm) {
      spec <- specs[[nm]]
      list(
        name = nm,
        title = spec$title,
        description = spec$description,
        inputSchema = spec$inputSchema
      )
    })),
    `tools/call` = mcp_call(request$params %||% list(), specs),
    NULL
  )
  if (is.null(result)) {
    if (is.null(id)) {
      return(NULL)
    }
    return(mcp_error_response(id, -32601L, paste0("method not found: ", method)))
  }
  list(jsonrpc = "2.0", id = id, result = result)
}
