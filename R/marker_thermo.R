# =============================================================================
#  marker_thermo.R  --  Tm, dG and PCR efficiency (packages)
# =============================================================================
# =============================================================================
#  2. Thermodynamics via packages
# =============================================================================

#' Give a function a new enclosing environment, keeping its formals and body
#'
#' The result is a fresh closure, so rebinding it cannot mutate the caller's
#' copy of `FUN`.  `base::` is explicit on purpose: here the unqualified
#' `body<-` is not found, and the assignment then silently leaves `FUN`
#' unchanged -- which looks exactly like a working rebinding.
#' @param FUN a closure.
#' @param env its new enclosing environment.
#' @return `FUN` with `environment()` set to `env`.
#' @keywords internal
#' @noRd
rebind <- function(FUN, env) {
  args <- formals(FUN)
  out <- function(...) NULL
  formals(out) <- args
  base::body(out) <- base::body(FUN)
  environment(out) <- env
  out
}

#' Re-home the package's functions into one plain environment
#'
#' Serializing a function whose enclosing environment is this namespace ships
#' only its *name*, and the worker rebinds the name to the installed namespace.
#' A function re-homed into a plain environment that is not a namespace or global
#' environment is serialized *with* that environment's contents instead, and the
#' contents are this package's own definitions.
#'
#' The names the task closure uses but does not bind itself (`genome`, `radius`
#' in a scan over a genome) are read from its own enclosing environment, and
#' travel as themselves; a name the task does bind -- a local `scan_one` -- is
#' rebuilt here, so it sees the carrier too.
#' @param FUN the task closure.
#' @return a copy of `FUN` closed over the carrier.
#' @keywords internal
#' @noRd
source_carrier <- function(FUN) {
  ns <- parent.env(environment())
  here <- environment(FUN)
  carrier <- list2env(as.list(here, all.names = TRUE), parent = globalenv())
  for (n in sort(ls(ns, all.names = TRUE))) {
    if (exists(n, envir = carrier, inherits = FALSE)) {
      next
    }
    obj <- get(n, envir = ns)
    # functions and plain data only: everything else here is an environment
    # (this very namespace, the S3 table), which serializes by reference -- a
    # worker would then reload the installed one and be back to stale code
    if (is.function(obj)) {
      assign(n, rebind(obj, carrier), envir = carrier)
    } else if (is.atomic(obj) || is.list(obj)) {
      assign(n, obj, envir = carrier)
    }
  }
  for (n in ls(here, all.names = TRUE)) {
    if (is.function(get(n, envir = here))) {
      assign(n, rebind(get(n, envir = here), carrier), envir = carrier)
    }
  }
  carrier
}

#' GC fraction and homopolymer-run test
#' @keywords internal
#' @noRd
gc_frac <- function(seq) unname(Biostrings::letterFrequency(Biostrings::DNAStringSet(seq), "GC", as.prob = TRUE)[, 1])

#' Regex for a tandem repeat of `unit`-length units, at least `min_copies` times
#'
#' The one place the homopolymer / SSR pattern is built, so the run screen, the
#' SSR screen and the screening scan all agree.
#' @keywords internal
#' @noRd
repeat_regex <- function(unit, min_copies = 3L) {
  # sprintf, not glue: this is called once per sequence per unit on the candidate
  # grid, and glue's data masking dominated the design profile
  sprintf("([ACGT]{%d})\\1{%d,}", as.integer(unit), as.integer(min_copies) - 1L)
}

has_run <- function(seq, n = 4L) stringr::str_detect(seq, repeat_regex(1L, n))

#' Make a parallel worker able to dispatch on Biostrings objects
#'
#' A `SnowParam`/SOCK worker is a fresh R process.  It resolves the package's
#' closures, but unless `Biostrings` is loaded there, S4 dispatch on a serialized
#' `DNAStringSet` (e.g. `genome[[i]]`, `subseq()`) falls back to base and fails
#' with "this S4 class is not subsettable".  Loading the namespace in the worker
#' fixes the dispatch; fork-based backends ([BiocParallel::MulticoreParam()])
#' inherit loaded namespaces and are unaffected.  Called at the top of every
#' worker closure that touches the genome.
#' @keywords internal
#' @noRd
worker_biostrings <- function() requireNamespace("Biostrings", quietly = TRUE)

#' Parallel map with an interactive progress bar
#'
#' A thin wrapper so every parallel map shows a progress bar only when
#' interactive.  Reproducible RNG is the backend's job: set `RNGseed` on the
#' `BPPARAM` (`MulticoreParam(RNGseed = 1)`); the design kernels here are
#' deterministic, so nothing needs seeding (`bpoptions(RNGseed = )` is not
#' supported by BiocParallel).
#'
#' On a socket backend the worker resolves a package closure against the
#' *installed* copy when one is present: BiocParallel exports the enclosing
#' environment of the task, and that environment's parent is the namespace the
#' function was found in.  A `load_all()` session therefore runs stale installed
#' code in the worker -- the same call that works serially fails there with
#' "unused argument".  For a **serial** backend the task runs in this session,
#' where the namespace is already current, so nothing is rewritten.
#' @keywords internal
#' @noRd
bp_lapply <- function(X, FUN, BPPARAM) {
  if (!serial_backend(BPPARAM)) {
    # a worker process resolves every package name against the *installed* copy;
    # hand it this session's definitions instead (see source_carrier())
    FUN <- rebind(FUN, source_carrier(FUN))
  }
  BiocParallel::bplapply(
    X, FUN,
    BPPARAM = BPPARAM,
    BPOPTIONS = BiocParallel::bpoptions(progressbar = interactive())
  )
}

#' Does this backend run the task in this session?
#' @param BPPARAM a [BiocParallel::BiocParallelParam()] object.
#' @return `TRUE` for the serial and fork backends, which need no carrier;
#'   `FALSE` for a socket backend, which starts fresh R processes.
#' @keywords internal
#' @noRd
serial_backend <- function(BPPARAM) {
  methods::is(BPPARAM, "SerialParam") || methods::is(BPPARAM, "MulticoreParam")
}


#' Tandem repeats (simple sequence repeats) in a sequence
#'
#' EasyKASP screens candidate primers for SSR repeat units, not just
#' homopolymers.  Returns, for each input, the repeated units of length `units`
#' occurring at least `min_copies` times in tandem.
#' @param seq character vector of DNA sequences.
#' @param units repeat-unit lengths to test.
#' @param min_copies tandem copies required for a hit.
#' @return a list, one character vector per input, of the repeated units (unique,
#'   possibly empty).
#' @examples
#' find_repeats(c("ACACACAC", "ACGTACGT"))
#' @family KASP QC
#' @export
find_repeats <- function(seq, units = 2:4, min_copies = 3L) {
  pats <- repeat_regex(units, min_copies)
  purrr::map(seq, function(s) {
    purrr::map(pats, function(p) stringr::str_match_all(s, p)[[1L]][, 2L]) |>
      purrr::list_c() |>
      unique()
  })
}

#' Does a sequence contain an SSR? (see [find_repeats()])
#' @param seq character vector of DNA sequences.
#' @param units repeat-unit lengths to test.
#' @param min_copies tandem copies required for a hit.
#' @return a logical vector, one per input sequence.
#' @examples
#' has_repeat(c("ACACACAC", "ACGT"))
#' @family KASP QC
#' @export
has_repeat <- function(seq, units = 2:4, min_copies = 3L) {
  if (!length(units)) {
    return(rep(FALSE, length(seq)))
  }
  # a boolean question must not extract every match: one vectorised detect per
  # unit over the whole candidate vector
  pats <- repeat_regex(units, min_copies)
  Reduce(`|`, lapply(pats, function(p) stringr::str_detect(seq, p)))
}

#' Inverted repeats / palindromes in a primer (library screen)
#'
#' `Biostrings::findPalindromes()` finds reverse-complementary arms - the
#' structural signature of a hairpin.  This complements the thermodynamic
#' hairpin screen in [primer_thermo()] with a sequence-level flag that does not
#' need OligoArrayAux.  The defaults are deliberately strict (arms >= 4 bp, loop
#' <= 6 bp, no mismatches): a well-designed primer returns no hits, whereas a
#' repetitive or self-symmetric one returns several.
#'
#' @param seq character vector of DNA sequences.
#' @param min_armlength minimum arm length of a palindromic hit.
#' @param max_looplength maximum loop (spacer) length.
#' @param max_mismatch mismatches allowed in the arms.
#' @return a list, one tibble per input, with `start`, `end`, `width` of each
#'   palindromic view (empty tibble when none).
#' @examples
#' find_palindromes(c(hit = "ACGTACGT", none = "AACCGG"))
#' @family KASP QC
#' @export
find_palindromes <- function(seq, min_armlength = 4L, max_looplength = 6L, max_mismatch = 0L) {
  nm <- names(seq)
  hits <- purrr::map(as.character(seq), function(s) {
    h <- Biostrings::findPalindromes(Biostrings::DNAString(s),
      min.armlength = min_armlength,
      max.looplength = max_looplength,
      max.mismatch = max_mismatch
    )
    tibble::tibble(
      start = as.integer(IRanges::start(h)),
      end = as.integer(IRanges::end(h)),
      width = as.integer(IRanges::width(h))
    )
  })
  rlang::set_names(hits, nm)
}

#' Does a primer contain an inverted repeat? (see [find_palindromes()])
#' @param seq character vector of DNA sequences.
#' @param min_armlength minimum arm length of a palindromic hit.
#' @param max_looplength maximum loop (spacer) length.
#' @param max_mismatch mismatches allowed in the arms.
#' @return a logical vector, one per input sequence.
#' @examples
#' has_palindrome(c("ACGTACGT", "AACCGG"))
#' @family KASP QC
#' @export
has_palindrome <- function(seq, min_armlength = 4L, max_looplength = 6L, max_mismatch = 0L) {
  purrr::map_int(find_palindromes(seq, min_armlength, max_looplength, max_mismatch), nrow) > 0L
}

#' Nearest-neighbour Tm of oligos (perfect duplex), salt-corrected (TmCalculator)
#'
#' `res$gr` uses TmCalculator's `$` accessor for the result GRanges; `Tm` is a
#' metadata column on it.  `nn_table` selects the nearest-neighbour parameter set
#' (TmCalculator ships several: SantaLucia 2004 is the package default, and
#' Breslauer / Allawi / Sugimoto / Ghosh are available for an independent
#' thermodynamic cross-check, see [thermo_crosscheck()]).
#' @param seq oligo(s) 5'->3'.
#' @param DMSO DMSO concentration (%), as in a KASP master mix; it lowers Tm.
#'   `NULL` (the default) uses `conditions$DMSO`.
#' @param nn_table nearest-neighbour parameter table (see details).
#' @param conditions reaction conditions from [reaction_conditions()] (or
#'   [scheme_conditions()], or a partial named list); its `Na`, `Mg`, `dNTPs` and
#'   `DMSO` fields are used.  `temperature` does not enter the nearest-neighbour
#'   Tm calculation (the NN model gives Tm directly), so passing
#'   `scheme_conditions("lamp")` shifts the Tm only through the buffer fields.
#'   `NULL` uses [reaction_conditions()], i.e. the package's historical numbers.
#'   For an isothermal method whose buffer this model does not represent (RPA's
#'   crowded recombinase buffer) the Tm is informational, not a gate.
#' @param salt_method monovalent-salt correction passed to
#'   `TmCalculator::tm_calculate()`.  `NULL` keeps the engine default
#'   (`"Schildkraut2010"`, Schildkraut & Lifson 1965, `16.6*log10[Na+]`), so
#'   existing designs do not move.  A caller who has a measured buffer can ask
#'   for a modern form: `"Owczarzy2004"` / `"Owczarzy2008"` (Owczarzy et al.
#'   2004, 2008, which correct `1/T` and include a duplex-length term) or
#'   `"SantaLucia1998-1"` / `"SantaLucia1998-2"` (SantaLucia 1998,
#'   `11.7*log10[Na+]`).  This is a level shift, not a precision gain: on the
#'   mm2020 oracle (500 mM Na equivalent, no Mg) the default reads 16.9 C low
#'   and Owczarzy 2008 reads 9.3 C low, with the same correlation; see
#'   `.tmp/sdd/report-tmfix.md`.
#' @return a numeric vector of Tm (C), one per oligo.
#' @examples
#' primer_tm("GAAGGTGACCAAGTTCATGCT")
#' primer_tm(c("ACGTACGTACGTACGTACGT", "TTTTGGCGCGCCGATATATAT"))
#' primer_tm("GAAGGTGACCAAGTTCATGCT", conditions = reaction_conditions(Na = 100))
#' @family KASP thermodynamics
#' @export
primer_tm <- function(seq, DMSO = NULL, nn_table = NN_TABLE,
                      conditions = NULL, salt_method = NULL) {
  conditions <- as_reaction_conditions(conditions)
  DMSO <- DMSO %||% conditions$DMSO
  seq <- as.character(seq)
  gr <- TmCalculator::to_genomic_ranges(seq)
  res <- TmCalculator::tm_calculate(gr,
    method = "tm_nn", nn_table = nn_table,
    Na = conditions$Na, Mg = conditions$Mg, dNTPs = conditions$dNTPs, DMSO = DMSO,
    salt_method = salt_method %||% "Schildkraut2010"
  )
  as.numeric(S4Vectors::mcols(res$gr)$Tm)
}

#' Is OligoArrayAux (`hybrid-min`) available?
#'
#' DECIPHER's thermodynamics runs through OligoArrayAux.  It is looked up as
#' `options(tidymarkers.oligoarrayaux)`, a project `.tools/`, the user cache (see
#' [marker_install_tools()]) and finally the `PATH`.
#'
#' @return a single logical.
#' @examples
#' has_oligoarrayaux()
#' @family KASP utilities
#' @export
has_oligoarrayaux <- function() !is.null(oligoarrayaux_paths())

#' Point OligoArrayAux at its data tables for the duration of `code`
#'
#' The thermodynamics engine (DECIPHER via OligoArrayAux) needs `UNAFOLDDAT` to
#' find its parameter tables.  This temporarily sets it, so you need not fiddle
#' with the environment by hand.
#'
#' @param path directory with the OligoArrayAux data tables.
#' @param code expression to evaluate (e.g. `assemble_kasp_assay(...)`).
#' @return the value of `code`, evaluated with `UNAFOLDDAT` set to `path`.
#' @examples
#' \dontrun{
#' with_unafloddat("/path/to/oligoarrayaux", assemble_kasp_assay(genome, snp))
#' }
#' @family KASP utilities
#' @export
with_unafloddat <- function(path, code) {
  if (!fs::dir_exists(path)) cli::cli_abort(gettext("Directory not found: {path}"))
  withr::with_envvar(c(UNAFOLDDAT = path), code)
}

#' Turn a core count into a `BiocParallel` backend (or `NULL` for one core)
#'
#' The single place `cores` becomes a `BPPARAM`: fork on Unix, a worker pool
#' elsewhere, and `NULL` for a single core - the CRAN-friendly default.  Pass a
#' `BPPARAM` explicitly to override, which is also how a caller avoids nesting
#' (an outer locus-level backend should force the inner one to `NULL`).
#' @keywords internal
#' @noRd
marker_bpparam <- function(cores = getOption("tidymarkers.cores", 1L)) {
  n <- max(1L, as.integer(cores))
  if (n <= 1L) {
    return(NULL)
  }
  if (.Platform$OS.type == "unix") BiocParallel::MulticoreParam(n) else BiocParallel::SnowParam(n)
}

#' Worker count of a `BPPARAM` (1 when `NULL`)
#' @keywords internal
#' @noRd
marker_ncores <- function(BPPARAM = NULL) {
  if (is.null(BPPARAM)) 1L else max(1L, BiocParallel::bpnworkers(BPPARAM))
}

#' OligoArrayAux dG for many oligos, optionally in parallel
#'
#' A vectorised `oligo_dg()` (duplex when `seq2` is given, else hairpin folding).
#' Each pair is its own `hybrid-min` subprocess, so with a `BPPARAM` the batch
#' runs across cores - the COM candidates and the QC maps are the per-SNP hot
#' spot.
#' @keywords internal
#' @noRd
oligo_dg_vec <- function(seq1, seq2 = NULL, BPPARAM = NULL,
                         na = 0.05, temp = 37) {
  n <- length(seq1)
  if (n == 0L) {
    return(numeric(0))
  }
  one <- function(i) {
    if (is.null(seq2)) {
      oligo_dg(seq1[[i]], na = na, temp = temp)
    } else {
      oligo_dg(seq1[[i]], seq2[[i]], na = na, temp = temp)
    }
  }
  if (is.null(BPPARAM)) {
    vapply(seq_len(n), one, numeric(1))
  } else {
    purrr::list_c(bp_lapply(seq_len(n), one, BPPARAM))
  }
}

#' Hairpin (self-folding) dG at a scheme's reaction temperature
#'
#' The hairpin counterpart of a temperature-aware [dimer_thermo()]: the same
#' OligoArrayAux folding engine, evaluated at the temperature the chemistry runs
#' at ([scheme_conditions()]) instead of a fixed 37 C, and with the salt taken
#' from the conditions object.  A higher temperature destabilises the fold, so
#' `dg` is never more negative than at 37 C for the same sequence.
#'
#' The duplex counterpart's signature (`dimer_thermo(primers, BPPARAM,
#' conditions)`) lives in `R/marker_qc.R`; this helper and [duplex_dg_thermo()]
#' exist so callers that need the reaction temperature do not have to go through a
#' global.  Passing a conditions object with its default `temperature` (37 C)
#' reproduces the historical numbers exactly.
#'
#' @param seq oligo(s) 5'->3'.
#' @param conditions reaction conditions from [reaction_conditions()] or
#'   [scheme_conditions()]; `Na` (mM) and `temperature` (C) are used.  `NULL` uses
#'   [reaction_conditions()], i.e. 50 mM / 37 C, the package's historical numbers.
#' @param BPPARAM optional `BiocParallel` backend.
#' @return a numeric vector of dG (kcal/mol), one per oligo; more negative = a
#'   more stable fold.  `NA_real_` without OligoArrayAux.
#' @examples
#' \dontrun{
#' hairpin_dg_thermo("GAAGGTGACCAAGTTCATGCT", conditions = scheme_conditions("lamp"))
#' }
#' @seealso [primer_tm()], [duplex_dg_thermo()], [scheme_conditions()]
#' @family KASP thermodynamics
#' @export
hairpin_dg_thermo <- function(seq, conditions = NULL, BPPARAM = NULL) {
  conditions <- as_reaction_conditions(conditions)
  oligo_dg_vec(
    as.character(seq),
    BPPARAM = BPPARAM,
    na = conditions$Na / 1000, temp = conditions$temperature
  )
}

#' Duplex (oligo pair) dG at a scheme's reaction temperature
#'
#' The temperature-aware counterpart of [dimer_thermo()] for named oligo pairs:
#' the same `hybrid-min` engine, with salt and temperature taken from the
#' conditions object.  [dimer_thermo()] itself already honours a
#' `conditions$temperature` (its default is 37 C, the historical behaviour); this
#' vector helper is the one-oligo-per-element form the per-candidate loops and
#' their tests use, so a scheme can evaluate its dimers at its own temperature
#' without changing any existing caller.
#'
#' @param seq1,seq2 oligo vectors 5'->3', paired element-wise.
#' @param conditions reaction conditions from [reaction_conditions()] or
#'   [scheme_conditions()]; `Na` (mM) and `temperature` (C) are used.  `NULL` uses
#'   [reaction_conditions()], i.e. 50 mM / 37 C.
#' @param BPPARAM optional `BiocParallel` backend.
#' @return a numeric vector of dG (kcal/mol); more negative = a stronger duplex.
#'   `NA_real_` without OligoArrayAux.
#' @examples
#' \dontrun{
#' duplex_dg_thermo(
#'   "GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC",
#'   conditions = scheme_conditions("pcr")
#' )
#' }
#' @seealso [dimer_thermo()], [hairpin_dg_thermo()], [scheme_conditions()]
#' @family KASP thermodynamics
#' @export
duplex_dg_thermo <- function(seq1, seq2, conditions = NULL, BPPARAM = NULL) {
  conditions <- as_reaction_conditions(conditions)
  oligo_dg_vec(
    as.character(seq1), as.character(seq2),
    BPPARAM = BPPARAM,
    na = conditions$Na / 1000, temp = conditions$temperature
  )
}

#' PCR amplification efficiency of primer/target pairs (DECIPHER)
#'
#' Models hybridisation (incl. side reactions) times Taq elongation efficiency,
#' so a 3'-terminal mismatch is strongly penalised - the physical basis of KASP.
#'
#' @param primers character vector, 5'->3'.
#' @param targets character vector (same length) with the reverse-complement of
#'   the binding site (as in `?DECIPHER::CalculateEfficiencyPCR`).
#' @param temp annealing temperature (C).
#' @param P primer concentration (M).
#' @param ions monovalent + divalent salt (M).
#' @param cores processors for the vectorised DECIPHER call (default
#'   `getOption("tidymarkers.cores", 1L)`; the many OligoArrayAux subprocesses are
#'   the hot spot, so >1 helps a lot, and 1 keeps CRAN happy).
#' @return a numeric vector of efficiencies (0-1), one per primer/target pair.
#' @examples
#' \dontrun{
#' pcr_efficiency("GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC")
#' }
#' @family KASP thermodynamics
#' @export
pcr_efficiency <- function(primers, targets, temp = PCR_EFF$temp,
                           P = PCR_EFF$P, ions = PCR_EFF$ions,
                           cores = getOption("tidymarkers.cores", 1L)) {
  rlang::check_installed("DECIPHER")
  if (!has_oligoarrayaux()) {
    cli::cli_abort(c(gettext("OligoArrayAux ('hybrid-min') is required for PCR efficiencies."),
      i = gettext("run marker_install_tools(), or build it with data-raw/install-oligoarrayaux.R")
    ))
  }
  primers <- as.character(primers)
  targets <- as.character(targets)
  with_oligoarrayaux(as.numeric(DECIPHER::CalculateEfficiencyPCR(
    primer = Biostrings::DNAStringSet(primers),
    target = Biostrings::DNAStringSet(targets),
    temp = temp, P = P, ions = ions, taqEfficiency = TRUE,
    processors = max(1L, as.integer(cores))
  )))
}
