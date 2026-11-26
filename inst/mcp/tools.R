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
#  The registry covers the workflow the vignettes show: read a genome/VCF,
#  classify a variant, screen targets, design one marker or a batch by scheme,
#  and read the verdict.  Thermo/QC helpers and the engine probes are exposed
#  as such.  Calls that cannot cross the wire are deliberately absent: the
#  plotting calls and the pipeline objects from marker_plan()/marker_solve()
#  (both are R structures, not JSON), the writers (marker_export(),
#  write_marker_results()) because an MCP call should not write outside its
#  mount, and the Ensembl lookups (ensembl_gene()/ensembl_variants()) because
#  they need the network and the server is offline by default.
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

mcp_num <- function(args, name) {
  value <- mcp_arg(args, name)
  if (is.null(value)) NULL else as.numeric(value)
}

# A primer set arrives as an array (labelled p1, p2, ...) or as an object whose
# names are the labels the caller wants back.
mcp_primers <- function(args, name = "primers") {
  value <- mcp_named_chr(args, name)
  if (is.null(value)) {
    value <- mcp_chr(args, name)
  }
  if (is.null(value)) {
    return(NULL)
  }
  if (is.null(names(value)) || any(!nzchar(names(value)))) {
    names(value) <- paste0("p", seq_along(value))
  }
  value
}

# One locus as the one-row tibble the design calls expect.
mcp_snp <- function(args) {
  tibble::tibble(
    chrom = as.character(mcp_arg(args, "chrom")),
    pos = as.integer(mcp_arg(args, "pos")),
    ref = as.character(mcp_arg(args, "ref")),
    alt = as.character(mcp_arg(args, "alt"))
  )
}

# Everything a caller needs from a designed marker that survives JSON: the
# class, the pass verdict, the primers, risk, the check table and what the
# scheme does or does not simulate. An extractor that does not apply to this
# marker returns NULL, which means "not applicable", not a failed call.
mcp_marker_summary <- function(marker) {
  if (is.null(marker)) {
    return(NULL)
  }
  take <- function(f) tryCatch(f(marker), error = function(e) NULL)
  primers <- take(marker_primers)
  list(
    class = class(marker)[1L],
    passes = tryCatch(isTRUE(passes_all(marker)), error = function(e) NA),
    primers = if (is.data.frame(primers) || is.character(primers)) primers else NULL,
    risk = take(marker_risk),
    checks = take(marker_explain),
    readout = take(marker_readout)
  )
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
  marker_genome(
    path = path,
    chrom_map = mcp_named_chr(args, "chrom_map"),
    blast_db = mcp_arg(args, "blast_db"),
    fm_index = mcp_arg(args, "fm_index")
  )
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

# The chemistries the package registers, computed from the registry so the
# advertised enum cannot drift.  A design enum is the per-locus ladder
# (`rescue_strategy()$scheme`): it drops `primer` (no locus designer, C13) and
# `duplex` (two loci), which the single-locus design tools cannot take.  A
# criteria/decision-text enum keeps every scheme, because those calls evaluate
# any of them.
mcp_scheme_ids <- function(design = FALSE) {
  if (design) rescue_strategy()$scheme else marker_schemes()$scheme
}

# Every entry: title, description, inputSchema (JSON Schema draft 2020-12), and
# a handler taking the parsed `arguments` object.
mcp_tool_specs <- function() {
  specs <- list(
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
          scheme = list(
            type = "string",
            enum = mcp_scheme_ids(),
            description = "scheme id (default \"snp\")"
          )
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
          schemes = list(
            type = "array",
            items = list(type = "string", enum = mcp_scheme_ids(design = TRUE)),
            description = "scheme ids to try, in order (one per registered designer; \"primer\" is excluded)"
          )
        ),
        required = list("genome_path", "chrom", "pos", "ref", "alt")
      ),
      handler = function(args) {
        result <- design_auto_marker(
          mcp_genome(args), mcp_snp(args),
          schemes = mcp_chr(args, "schemes")
        )
        list(
          scheme = result$scheme,
          tried = result$tried[, c("scheme", "ok"), drop = FALSE]
        )
      }
    ),
    variant_type = list(
      title = "Variant type",
      description = "Classify a REF/ALT pair: SNP, insertion, deletion, MNP or complex.",
      inputSchema = list(
        type = "object",
        properties = list(
          ref = list(type = "string", description = "reference allele"),
          alt = list(type = "string", description = "alternate allele")
        ),
        required = list("ref", "alt")
      ),
      handler = function(args) {
        variant_type(mcp_arg(args, "ref"), mcp_arg(args, "alt"))
      }
    ),
    marker_site_class = list(
      title = "Site class of a multi-allelic locus",
      description = paste(
        "Group a REF and its ALT alleles by the IUPAC code a primer must carry",
        "at the 3' end to cover all of them."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          ref = list(type = "string", description = "reference allele"),
          alts = list(type = "array", items = list(type = "string"), description = "alternate alleles")
        ),
        required = list("ref", "alts")
      ),
      handler = function(args) {
        marker_site_class(mcp_arg(args, "ref"), mcp_chr(args, "alts"))
      }
    ),
    mask_positions = list(
      title = "Mask variant positions",
      description = paste(
        "Replace the given 1-based positions of a DNA sequence with N, so a",
        "neighbouring variant cannot leak into a designed primer."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          seq = list(type = "string", description = "DNA sequence"),
          positions = list(
            type = "array", items = list(type = "integer"),
            description = "1-based positions to mask"
          )
        ),
        required = list("seq", "positions")
      ),
      handler = function(args) {
        as.character(mask_positions(
          mcp_arg(args, "seq"),
          as.integer(unlist(mcp_arg(args, "positions")))
        ))
      }
    ),
    asq_limits = list(
      title = "Allele-specific qPCR limits",
      description = paste(
        "The amplicon and Tm windows an ASQ/allele-specific qPCR design must sit",
        "in, and the flags each window carries."
      ),
      inputSchema = list(type = "object", properties = list()),
      handler = function(args) asq_limits()
    ),
    reaction_conditions = list(
      title = "Reaction conditions",
      description = paste(
        "The salt, dNTP, DMSO and temperature values the thermodynamics engines",
        "read. Absent fields keep the package defaults."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          Na = list(type = "number", description = "sodium, mM"),
          Mg = list(type = "number", description = "magnesium, mM"),
          dNTPs = list(type = "number", description = "dNTPs, mM"),
          DMSO = list(type = "number", description = "DMSO, percent"),
          temperature = list(type = "number", description = "temperature, Celsius")
        )
      ),
      handler = function(args) {
        over <- list(
          Na = mcp_num(args, "Na"), Mg = mcp_num(args, "Mg"),
          dNTPs = mcp_num(args, "dNTPs"), DMSO = mcp_num(args, "DMSO"),
          temperature = mcp_num(args, "temperature")
        )
        do.call(reaction_conditions, over[!vapply(over, is.null, logical(1))])
      }
    ),
    scheme_conditions = list(
      title = "Scheme reaction conditions",
      description = "The reaction conditions a marker scheme defaults to.",
      inputSchema = list(
        type = "object",
        properties = list(
          # No enum: the handler documents a scheme id *or* a temperature
          # family, and treats an unregistered id as PCR-style.
          scheme = list(
            type = "string",
            description = "scheme id or pcr/lamp/rpa/ambient (default \"pcr\")"
          )
        )
      ),
      handler = function(args) {
        scheme_conditions(mcp_arg(args, "scheme") %||% "pcr")
      }
    ),
    primer_thermo = list(
      title = "Primer thermodynamics",
      description = paste(
        "Length, Tm, GC, hairpin dG and terminal dG for each primer, in one",
        "table. The dG columns are null without the OligoArrayAux engine."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          primers = list(type = "array", items = list(type = "string"), description = "primers, 5'->3'")
        ),
        required = list("primers")
      ),
      handler = function(args) primer_thermo(mcp_chr(args, "primers"))
    ),
    primer_set_compatibility = list(
      title = "Primer set compatibility",
      description = paste(
        "Pairwise complementarity of a primer set: self versus cross runs, the",
        "cross-dimer dG, and whether each pair passes the tube thresholds."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          primers = list(type = "array", items = list(type = "string"), description = "primers, 5'->3'"),
          min_cross_run = list(type = "integer", description = "minimum cross complementarity run (default 6)"),
          max_selfcomp = list(type = "integer", description = "maximum self complementarity run (default 5)")
        ),
        required = list("primers")
      ),
      handler = function(args) {
        primer_set_compatibility(
          mcp_primers(args),
          min_cross_run = mcp_int(args, "min_cross_run") %||% 6L,
          max_selfcomp = mcp_int(args, "max_selfcomp") %||% 5L
        )
      }
    ),
    dimer_thermo = list(
      title = "Dimer free energy",
      description = paste(
        "dG of every distinct pair in a primer set (self pairs included). Primers",
        "may arrive as an array or as an object with labels."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          primers = list(type = "array", items = list(type = "string"), description = "primers, 5'->3'")
        ),
        required = list("primers")
      ),
      handler = function(args) dimer_thermo(mcp_primers(args))
    ),
    hairpin_dg_thermo = list(
      title = "Hairpin free energy",
      description = "The hairpin dG of one oligo at the reaction conditions.",
      inputSchema = list(
        type = "object",
        properties = list(
          seq = list(type = "string", description = "DNA sequence, 5'->3'")
        ),
        required = list("seq")
      ),
      handler = function(args) hairpin_dg_thermo(mcp_arg(args, "seq"))
    ),
    fold_amplicon = list(
      title = "Fold an amplicon",
      description = paste(
        "Minimum free energy and dot-bracket structure of a single strand, or",
        "status `unknown` with the reason when no folding engine is present."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          seq = list(type = "string", description = "sequence to fold"),
          kind = list(type = "string", enum = c("DNA", "RNA"), description = "nucleic acid (default DNA)")
        ),
        required = list("seq")
      ),
      handler = function(args) {
        fold_amplicon(mcp_arg(args, "seq"), kind = mcp_arg(args, "kind") %||% "DNA")
      }
    ),
    pcr_efficiency = list(
      title = "PCR amplification efficiency",
      description = paste(
        "Primer/target amplification efficiency (0-1) from the thermodynamics",
        "model. Requires OligoArrayAux and DECIPHER."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          primers = list(type = "array", items = list(type = "string"), description = "primers, 5'->3'"),
          targets = list(
            type = "array", items = list(type = "string"),
            description = "reverse-complement of each binding site, same length as primers"
          ),
          temp = list(type = "number", description = "annealing temperature, Celsius"),
          cores = list(type = "integer", description = "processors for the vectorised call")
        ),
        required = list("primers", "targets")
      ),
      handler = function(args) {
        extra <- list(temp = mcp_num(args, "temp"), cores = mcp_int(args, "cores"))
        do.call(pcr_efficiency, c(
          list(mcp_chr(args, "primers"), mcp_chr(args, "targets")),
          extra[!vapply(extra, is.null, logical(1))]
        ))
      }
    ),
    wobble_mismatch = list(
      title = "Wobble mismatch",
      description = paste(
        "Does an on/off-allele pair form the worst-case G:T or A:C wobble, which",
        "does not stop extension?"
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          on_allele = list(type = "string", description = "on-target allele base"),
          off_allele = list(type = "string", description = "off-target allele base")
        ),
        required = list("on_allele", "off_allele")
      ),
      handler = function(args) {
        wobble_mismatch(mcp_arg(args, "on_allele"), mcp_arg(args, "off_allele"))
      }
    ),
    find_palindromes = list(
      title = "Inverted repeats",
      description = paste(
        "Inverted repeats in a sequence that can fold into a hairpin, with the",
        "arm length and loop constraints."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          seq = list(type = "string", description = "DNA sequence"),
          min_armlength = list(type = "integer", description = "minimum arm length (default 4)"),
          max_looplength = list(type = "integer", description = "maximum loop length (default 6)"),
          max_mismatch = list(type = "integer", description = "allowed arm mismatches (default 0)")
        ),
        required = list("seq")
      ),
      handler = function(args) {
        find_palindromes(
          mcp_arg(args, "seq"),
          min_armlength = mcp_int(args, "min_armlength") %||% 4L,
          max_looplength = mcp_int(args, "max_looplength") %||% 6L,
          max_mismatch = mcp_int(args, "max_mismatch") %||% 0L
        )[[1L]]
      }
    ),
    has_palindrome = list(
      title = "Has an inverted repeat",
      description = "Whether a sequence folds into a hairpin under the given constraints.",
      inputSchema = list(
        type = "object",
        properties = list(
          seq = list(type = "string", description = "DNA sequence"),
          min_armlength = list(type = "integer", description = "minimum arm length (default 4)"),
          max_looplength = list(type = "integer", description = "maximum loop length (default 6)"),
          max_mismatch = list(type = "integer", description = "allowed arm mismatches (default 0)")
        ),
        required = list("seq")
      ),
      handler = function(args) {
        has_palindrome(
          mcp_arg(args, "seq"),
          min_armlength = mcp_int(args, "min_armlength") %||% 4L,
          max_looplength = mcp_int(args, "max_looplength") %||% 6L,
          max_mismatch = mcp_int(args, "max_mismatch") %||% 0L
        )
      }
    ),
    find_caps_enzymes = list(
      title = "CAPS/dCAPS enzymes",
      description = paste(
        "Restriction enzymes whose site is gained or lost by the ALT allele and",
        "whose digest fragments are resolvable on a gel. Needs a mounted genome."
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
          flank = list(type = "integer", description = "window around the locus (default 250)"),
          min_gap = list(type = "integer", description = "minimum resolvable band gap (default 50)")
        ),
        required = list("genome_path", "chrom", "pos", "ref", "alt")
      ),
      handler = function(args) {
        find_caps_enzymes(
          mcp_genome(args), mcp_snp(args),
          flank = mcp_int(args, "flank") %||% 250L,
          min_gap = mcp_int(args, "min_gap") %||% 50L
        )
      }
    ),
    rescue_strategy = list(
      title = "Rescue strategy",
      description = paste(
        "The ordered ladder of marker schemes to try for a locus, with each",
        "scheme's discrimination and readout."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          schemes = list(
            type = "array",
            items = list(type = "string", enum = mcp_scheme_ids(design = TRUE)),
            description = "scheme ids to order (one per registered designer; \"primer\" is excluded)"
          )
        )
      ),
      handler = function(args) rescue_strategy(mcp_chr(args, "schemes"))
    ),
    marker_decision_text = list(
      title = "Decision graph as text",
      description = paste(
        "A scheme's decision graph as text, one node per line. `evidence`",
        "attaches each criterion's evidence rung."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          scheme = list(
            type = "string",
            enum = mcp_scheme_ids(),
            description = "scheme id (default \"snp\")"
          ),
          evidence = list(type = "boolean", description = "include the evidence rung")
        )
      ),
      handler = function(args) {
        marker_decision_text(
          mcp_arg(args, "scheme") %||% "snp",
          evidence = isTRUE(mcp_arg(args, "evidence"))
        )
      }
    ),
    screen_snp_targets = list(
      title = "Screen SNP targets",
      description = paste(
        "For every variant in a VCF, count the flanking primer candidates that",
        "clear the 3' cleanliness, length and amplicon windows. Needs a genome."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          variants_path = list(type = "string", description = "VCF file"),
          genome_path = list(type = "string", description = "FASTA file, plain or gzipped"),
          chrom_map = list(type = "object", description = "optional defline -> name map")
        ),
        required = list("variants_path", "genome_path")
      ),
      handler = function(args) {
        path <- mcp_arg(args, "variants_path")
        if (is.null(path)) {
          stop("`variants_path` is required: a VCF file")
        }
        screen_snp_targets(read_variants_vcf(path), mcp_genome(args))
      }
    ),
    engine_status = list(
      title = "Engine availability",
      description = "Whether the external engines (OligoArrayAux, BLAST+) are installed.",
      inputSchema = list(type = "object", properties = list()),
      handler = function(args) {
        list(oligoarrayaux = has_oligoarrayaux(), blast = has_blast())
      }
    ),
    engine_capabilities = list(
      title = "Engine capabilities",
      description = paste(
        "The capability flags of an engine: Tm, salt, duplex/hairpin dG, fold for",
        "thermo; in-silico PCR, seed counts, IUPAC, mismatch and BLAST for",
        "specificity."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          engine = list(type = "string", enum = c("thermo", "specificity"), description = "engine (default thermo)")
        )
      ),
      handler = function(args) {
        capabilities(marker_engines()[[mcp_arg(args, "engine") %||% "thermo"]])
      }
    ),
    design_chemistries = list(
      title = "Chemistries available for design",
      description = paste(
        "The selection object a design call accepts: the per-locus ladder",
        "(`per_locus`, the valid `scheme`/`schemes` ids), every registered scheme,",
        "and the relative costs used to order the ladder."
      ),
      inputSchema = list(type = "object", properties = list()),
      handler = function(args) {
        list(
          per_locus = rescue_strategy()$scheme,
          all = marker_schemes()$scheme,
          costs = marker_costs()
        )
      }
    ),
    design_scheme = list(
      title = "Design a marker with one scheme",
      description = paste(
        "Design one locus with a named scheme and return the verdict, primers,",
        "risk, check table and readout. `params` is passed to the scheme's",
        "designer unchanged (e.g. `{ \"flank\": 200, \"cas\": \"cas12a\" }`)."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          genome_path = list(type = "string", description = "FASTA file, plain or gzipped"),
          chrom_map = list(type = "object", description = "optional defline -> name map"),
          blast_db = list(type = "string", description = "optional BLAST database prefix"),
          chrom = list(type = "string", description = "locus chromosome"),
          pos = list(type = "integer", description = "locus position (1-based)"),
          ref = list(type = "string", description = "reference allele"),
          alt = list(type = "string", description = "alternate allele"),
          scheme = list(
            type = "string",
            enum = mcp_scheme_ids(design = TRUE),
            description = "scheme id (default \"snp\"; \"primer\" has no locus-level designer)"
          ),
          params = list(type = "object", description = "extra arguments for the scheme's designer")
        ),
        required = list("genome_path", "chrom", "pos", "ref", "alt")
      ),
      handler = function(args) {
        params <- mcp_arg(args, "params")
        if (is.null(params)) {
          params <- list()
        }
        result <- do.call(design_auto_marker, c(
          list(mcp_genome(args), mcp_snp(args), schemes = mcp_arg(args, "scheme") %||% "snp"),
          params
        ))
        list(
          scheme = result$scheme,
          tried = result$tried[, c("scheme", "ok"), drop = FALSE],
          marker = mcp_marker_summary(result$marker)
        )
      }
    ),
    design_markers = list(
      title = "Design a batch of markers",
      description = paste(
        "Design every variant in a VCF, or a region of it, with one scheme, and",
        "return one verdict row per locus."
      ),
      inputSchema = list(
        type = "object",
        properties = list(
          genome_path = list(type = "string", description = "FASTA file, plain or gzipped"),
          chrom_map = list(type = "object", description = "optional defline -> name map"),
          variants_path = list(type = "string", description = "VCF file"),
          chrom = list(type = "string", description = "region chromosome"),
          start = list(type = "integer", description = "region start (1-based)"),
          end = list(type = "integer", description = "region end (inclusive)"),
          scheme = list(
            type = "string",
            enum = mcp_scheme_ids(design = TRUE),
            description = "scheme id (default \"snp\"; \"primer\" has no locus-level designer)"
          ),
          flank = list(type = "integer", description = "window around each locus")
        ),
        required = list("genome_path", "variants_path")
      ),
      handler = function(args) {
        path <- mcp_arg(args, "variants_path")
        if (is.null(path)) {
          stop("`variants_path` is required: a VCF file")
        }
        loci <- read_variants_vcf(
          path,
          chrom = mcp_arg(args, "chrom"),
          start = mcp_int(args, "start"),
          end = mcp_int(args, "end")
        )
        extra <- list()
        flank <- mcp_int(args, "flank")
        if (!is.null(flank)) {
          extra$flank <- flank
        }
        do.call(design_markers, c(
          list(mcp_genome(args), loci, scheme = mcp_arg(args, "scheme") %||% "snp"),
          extra
        ))
      }
    )
  )
  mcp_normalise_specs(mcp_annotate_specs(specs))
}

# MCP annotation hints. The registry reads files and computes values; only
# `marker_technology` changes session state. No tool opens the network, and a
# repeated call with the same arguments returns the same value.
MCP_ANNOTATIONS <- list(
  readOnlyHint = TRUE,
  destructiveHint = FALSE,
  idempotentHint = TRUE,
  openWorldHint = FALSE
)
MCP_MUTATING_ANNOTATIONS <- list(
  readOnlyHint = FALSE,
  destructiveHint = FALSE,
  idempotentHint = TRUE,
  openWorldHint = FALSE
)
MCP_MUTATING_TOOLS <- c("marker_technology")

mcp_annotate_specs <- function(specs) {
  for (nm in names(specs)) {
    specs[[nm]]$annotations <- if (nm %in% MCP_MUTATING_TOOLS) {
      MCP_MUTATING_ANNOTATIONS
    } else {
      MCP_ANNOTATIONS
    }
  }
  specs
}

# Every published argument object is closed (`additionalProperties = FALSE`)
# unless a schema opts out; the tool then never receives a name a client could
# not see in `properties`.
mcp_normalise_specs <- function(specs) {
  for (nm in names(specs)) {
    schema <- specs[[nm]]$inputSchema
    if (is.null(schema$additionalProperties)) {
      schema$additionalProperties <- FALSE
      specs[[nm]]$inputSchema <- schema
    }
  }
  specs
}

# -- validation ----------------------------------------------------------

# The runtime validator reads the same `inputSchema` that `tools/list`
# publishes, so the server never hands a handler an argument the client was
# told to send.  `path` names the failing argument, `got` is what arrived.
mcp_type_of <- function(value) {
  if (is.null(value)) {
    return("null")
  }
  if (is.list(value)) {
    return(if (is.null(names(value))) "array" else "object")
  }
  if (is.logical(value)) {
    return("boolean")
  }
  if (is.numeric(value)) {
    return("number")
  }
  if (is.character(value)) {
    return("string")
  }
  typeof(value)
}

mcp_violation <- function(path, expected, got) {
  list(path = path, expected = expected, got = got)
}

mcp_child_path <- function(path, name) {
  if (nzchar(path)) paste0(path, ".", name) else name
}

mcp_validate_value <- function(value, schema, path) {
  type <- schema$type
  if (!is.null(type)) {
    actual <- mcp_type_of(value)
    ok <- switch(type,
      string = is.character(value) && length(value) == 1L,
      number = is.numeric(value) && length(value) == 1L,
      integer = is.numeric(value) && length(value) == 1L && value == trunc(value),
      boolean = is.logical(value) && length(value) == 1L,
      # An empty JSON object and an empty JSON array are both `list()`; accept
      # either as an empty object rather than guess.
      object = is.list(value) &&
        (is.null(names(value)) && length(value) == 0L || !is.null(names(value))),
      array = is.list(value) && is.null(names(value)),
      actual == type
    )
    if (!ok) {
      return(mcp_violation(path, type, value))
    }
  }

  if (!is.null(schema$enum)) {
    choices <- as.list(schema$enum)
    if (!any(vapply(choices, identical, logical(1), value))) {
      return(mcp_violation(path, paste0("one of: ", paste(unlist(choices), collapse = ", ")), value))
    }
  }

  if (is.numeric(value) && length(value) == 1L) {
    if (!is.null(schema$minimum) && value < schema$minimum) {
      return(mcp_violation(path, paste0(">= ", schema$minimum), value))
    }
    if (!is.null(schema$maximum) && value > schema$maximum) {
      return(mcp_violation(path, paste0("<= ", schema$maximum), value))
    }
  }

  if (is.character(value) && length(value) == 1L && !is.null(schema$pattern)) {
    if (!grepl(schema$pattern, value)) {
      return(mcp_violation(path, paste0("pattern ", schema$pattern), value))
    }
  }

  object_like <- is.list(value) && !is.null(names(value))
  empty_object <- is.list(value) && length(value) == 0L && identical(type, "object")
  if (object_like || empty_object) {
    required <- unlist(schema$required %||% list(), use.names = FALSE)
    for (nm in required) {
      if (is.null(value[[nm]])) {
        return(mcp_violation(mcp_child_path(path, nm), "required argument", NULL))
      }
    }
    props <- schema$properties %||% list()
    if (isFALSE(schema$additionalProperties)) {
      extra <- setdiff(names(value), names(props))
      if (length(extra)) {
        return(mcp_violation(
          mcp_child_path(path, extra[[1L]]),
          paste0("one of: ", paste(names(props), collapse = ", ")),
          value[[extra[[1L]]]]
        ))
      }
    }
    for (nm in intersect(names(value), names(props))) {
      problem <- mcp_validate_value(value[[nm]], props[[nm]], mcp_child_path(path, nm))
      if (!is.null(problem)) {
        return(problem)
      }
    }
    return(NULL)
  }

  if (is.list(value) && is.null(names(value)) && !is.null(schema$items)) {
    for (i in seq_along(value)) {
      problem <- mcp_validate_value(value[[i]], schema$items, paste0(path, "[", i, "]"))
      if (!is.null(problem)) {
        return(problem)
      }
    }
  }

  NULL
}

# `args` is the decoded JSON `arguments` object: a named list (or `list()` for
# an absent/empty object).
mcp_validate_args <- function(args, schema) {
  if (!is.list(args)) {
    return(mcp_violation("arguments", "object", args))
  }
  mcp_validate_value(args, schema, "")
}

mcp_protocol_error <- function(message, tool, path, expected, got) {
  structure(
    list(
      message = message,
      data = list(tool = tool, path = path, expected = expected, got = got)
    ),
    class = "mcp_protocol_error"
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
# An unknown tool or an argument the schema rejects is JSON-RPC `-32602`
# ("invalid params"), with `data` naming the tool, path, expectation and value.
mcp_call <- function(params, specs) {
  if (!is.list(params)) {
    return(mcp_protocol_error("invalid params: `params` must be an object", "", "params", "object", params))
  }
  name <- params$name
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    return(mcp_protocol_error(
      "invalid params: `name` is required", "", "name", "string", params$name
    ))
  }
  spec <- specs[[name]]
  if (is.null(spec)) {
    return(mcp_protocol_error(
      paste0("unknown tool: ", name), name, "name", "a registered tool name", name
    ))
  }
  args <- params$arguments %||% list()
  problem <- mcp_validate_args(args, spec$inputSchema)
  if (!is.null(problem)) {
    label <- if (nzchar(problem$path)) problem$path else "arguments"
    return(mcp_protocol_error(
      paste0("invalid arguments for ", name, ": ", label),
      name, label, problem$expected, problem$got
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

mcp_error_response <- function(id, code, message, data = NULL) {
  error <- list(code = code, message = message)
  if (!is.null(data)) {
    error$data <- data
  }
  list(jsonrpc = "2.0", id = id, error = error)
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
        inputSchema = spec$inputSchema,
        annotations = spec$annotations
      )
    })),
    `tools/call` = mcp_call(request$params %||% list(), specs),
    NULL
  )
  if (inherits(result, "mcp_protocol_error")) {
    if (is.null(id)) {
      return(NULL)
    }
    return(mcp_error_response(id, -32602L, result$message, result$data))
  }
  if (is.null(result)) {
    if (is.null(id)) {
      return(NULL)
    }
    return(mcp_error_response(id, -32601L, paste0("method not found: ", method)))
  }
  list(jsonrpc = "2.0", id = id, result = result)
}
