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

#' Environments whose contents a socket worker must see
#'
#' The namespace, its `.__NAMESPACE__.` and `.__S3MethodsTable__.` must stay
#' by-reference (a worker would otherwise reload the installed copy, the bug
#' `source_carrier()` exists to avoid), and the option/cache state in
#' `.marker_state` is session-local.  A *registry* is neither: `.marker_tech_store`
#' holds the technologies composed with [marker_technology()], which are part of
#' what the task means, so its contents travel with the task.
#' @return a character vector of environment names to carry.
#' @keywords internal
#' @noRd
carried_registries <- function() c(".marker_tech_store")

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
    # functions and plain data; a registry environment is added below, and the
    # rest (this very namespace, the S3 table, the option/cache state) is left
    # out -- it must serialize by reference or not travel at all
    if (is.function(obj)) {
      assign(n, rebind(obj, carrier), envir = carrier)
    } else if (is.atomic(obj) || is.list(obj)) {
      assign(n, obj, envir = carrier)
    } else if (is.environment(obj) && n %in% carried_registries()) {
      # a runtime registry holds definitions the task needs, not session state.
      # Its parent is emptyenv(), so serializing it by value ships exactly its
      # own contents; the functions inside still resolve through the namespace.
      assign(n, obj, envir = carrier)
    }
  }
  for (n in ls(here, all.names = TRUE)) {
    if (is.function(get(n, envir = here))) {
      assign(n, rebind(get(n, envir = here), carrier), envir = carrier)
    }
  }
  # `.marker_state` is deliberately not carried (session state, above), but the
  # re-homed `session_memo()` still looks it up: without this the worker dies on
  # "object '.marker_state' not found".  A fresh, empty state per task is what a
  # worker should get anyway - the memo is cold there either way.
  assign(".marker_state", marker_state_new(), envir = carrier)
  carrier
}

#' GC fraction and homopolymer-run test
#' @keywords internal
#' @noRd
gc_frac <- function(seq) {
  # `letterFrequency(DNAStringSet(...))` spent 1.45 ms per call on a 27-mer, and
  # the count itself is one `gsub`.  The alphabet check is what the S4 call gave
  # for free and is kept: a silent `gc` for a bad character would reach a gate.
  # The denominator is `nchar()`, which is what `as.prob = TRUE` divides by (a
  # `+`/`-`/`.` in the DNA alphabet counts in it but is never GC).
  seq <- as.character(seq)
  n <- nchar(seq, keepNA = TRUE)
  if (anyNA(n) || any(grepl("[^ACGTMRWSYKVHDBNacgtmrwsykvhdbn+.-]", seq))) {
    abort_msg("`seq` must be DNA (ACGT plus IUPAC ambiguity codes).")
  }
  unname(nchar(gsub("[^GC]", "", toupper(seq)), keepNA = TRUE) / n)
}

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
#'
#' `chunksize` groups the input into tasks of that size, which is the unit a
#' worker receives; the result is flattened back into input order.  With
#' `on_error = "collect"` a failing element becomes a `tidymarkers_task_error`
#' (its `index` and `message`) and the rest still run; the default `"stop"`
#' keeps the historical fail-fast behaviour.  `RNGseed` is set on the backend
#' when given, because `bpoptions(RNGseed = )` is not supported by BiocParallel.
#' @param X input vector/list.
#' @param FUN function of one element.
#' @param BPPARAM a [BiocParallel::BiocParallelParam()] object.
#' @param chunksize elements per task; `NULL` keeps one task per element.
#' @param on_error `"stop"` (default) or `"collect"`.
#' @param RNGseed optional backend RNG seed.
#' @return a list in input order.
#' @keywords internal
#' @noRd
bp_lapply <- function(X, FUN, BPPARAM, chunksize = NULL,
                      on_error = c("stop", "collect"),
                      RNGseed = NULL) { # nolint: object_name_linter.
  on_error <- match.arg(on_error)
  if (length(X) == 0L) {
    return(list())
  }
  if (!is.null(RNGseed)) {
    BiocParallel::bpRNGseed(BPPARAM) <- as.integer(RNGseed)
  }
  run_one <- if (on_error == "stop") {
    FUN
  } else {
    function(x) {
      tryCatch(
        FUN(x),
        error = function(e) {
          structure(
            list(message = conditionMessage(e)),
            class = "tidymarkers_task_error"
          )
        }
      )
    }
  }
  if (!serial_backend(BPPARAM)) {
    # a worker process resolves every package name against the *installed* copy;
    # hand it this session's definitions instead (see source_carrier())
    run_one <- rebind(run_one, source_carrier(run_one))
  }
  n <- length(X)
  chunks <- if (is.null(chunksize)) {
    as.list(seq_len(n))
  } else {
    chunksize <- as.integer(chunksize)
    if (is.na(chunksize) || chunksize < 1L) {
      abort_msg("`chunksize` must be a positive integer.")
    }
    unname(split(seq_len(n), ceiling(seq_len(n) / chunksize)))
  }
  parts <- BiocParallel::bplapply(
    chunks, function(ii) lapply(ii, function(i) run_one(X[[i]])),
    BPPARAM = BPPARAM,
    BPOPTIONS = BiocParallel::bpoptions(progressbar = interactive())
  )
  unlist(parts, recursive = FALSE, use.names = FALSE)
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
#' @family QC
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
#' @family QC
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
#' @family QC
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
#'
#' IUPAC ambiguity codes are wildcards: two arm positions pair when the base
#' sets of their codes intersect (`N` pairs with anything, `R` with `Y`, ...).
#' A degenerate oligo is a mixture of concrete oligos, so if any member of the
#' mixture can form the inverted repeat the primer carries a real hairpin risk
#' and the QC gate must see it.  Plain-ACGT input is unchanged: it keeps the
#' exact substring fast path, and a wildcard comparison on ACGT is the same
#' test.  Wildcard matching needs `max_mismatch = 0`; a degenerate oligo asked
#' for mismatches is still screened by `findPalindromes()`, symbol-exactly.
#' @param seq character vector of DNA sequences.
#' @param min_armlength minimum arm length of a palindromic hit, in nt (one
#'   positive whole number).
#' @param max_looplength maximum loop (spacer) length, in nt (one non-negative
#'   whole number).
#' @param max_mismatch mismatches allowed in the arms (one non-negative whole
#'   number).
#' @return a logical vector, one per input sequence.
#' @examples
#' has_palindrome(c("ACGTACGT", "AACCGG"))
#' # wildcard: the `N` arms can be completed to ACGT...ACGT
#' has_palindrome("ACGTNNACGT")
#' @family QC
#' @export
has_palindrome <- function(seq, min_armlength = 4L, max_looplength = 6L, max_mismatch = 0L) {
  # `findPalindromes()` plus a tibble per sequence cost ~1.2 ms (measured); the
  # question is only "does any arm pair exist".  For a plain-ACGT oligo that is
  # a few vectorised substring comparisons (~60 us); an IUPAC one is the same
  # scan on 4-bit base masks.  A mismatch tolerance still takes the
  # Biostrings path, so the fast path cannot change an answer it does not
  # cover.
  seq <- as.character(seq)
  max_mismatch <- check_scan_radius(max_mismatch, "max_mismatch")
  min_armlength <- check_scan_radius(min_armlength, "min_armlength", min = 1L)
  max_looplength <- check_scan_radius(max_looplength, "max_looplength", min = 0L)
  fast <- identical(max_mismatch, 0L) && min_armlength >= 2L
  # unname(): the old `map_int(...) > 0L` dropped names because a comparison
  # strips attributes; callers must not start seeing names on this vector
  unname(vapply(seq, function(s) {
    if (is.na(s)) {
      # a missing sequence cannot form an inverted repeat; Biostrings aborts on it
      return(FALSE)
    }
    if (!fast) {
      found <- Biostrings::findPalindromes(
        Biostrings::DNAString(s),
        min.armlength = min_armlength,
        max.looplength = max_looplength,
        max.mismatch = max_mismatch
      )
      return(length(found) > 0L)
    }
    if (grepl("[^ACGT]", s)) {
      wildcard_arms(s, as.integer(min_armlength), as.integer(max_looplength))
    } else {
      palindromic_arms(s, as.integer(min_armlength), as.integer(max_looplength))
    }
  }, logical(1)))
}

#' Do two reverse-complementary arms of `k` bases sit within `L` of each other?
#'
#' `findPalindromes()` is an inverted-repeat finder: the arms match reverse-
#' complementarily and the loop between them is arbitrary (checked: it reports
#' `ACGTAAAAACGT` at `k = 4, L = 4`).  A view with arms `A >= k` and loop `d`
#' contains a `k`-arm pair with the same loop `d`, and any `k`-arm pair with
#' loop `d <= L` is such a view, so "any hit" is exactly "any `k`-arm pair".
#'
#' The reverse complement is taken once for the whole sequence and read back
#' through mirrored coordinates, so no window is complemented on its own.
#' @param s a plain-ACGT sequence.
#' @param k arm length.
#' @param L maximum loop length (minimum loop length is 0).
#' @return a single logical.
#' @keywords internal
#' @noRd
palindromic_arms <- function(s, k, L) {
  n <- nchar(s)
  if (n < 2L * k) {
    return(FALSE)
  }
  rc <- intToUtf8(rev(utf8ToInt(chartr("ACGT", "TGCA", s))))
  for (d in 0L:L) {
    npos <- n - 2L * k - d + 1L
    if (npos < 1L) break
    i <- seq_len(npos)
    left <- substring(s, i, i + k - 1L)
    # the right arm starts at i + k + d, so its reverse complement starts at
    # n - (i + k + d + k - 1) + 1 in the genome's reverse complement
    rs <- n - i - 2L * k - d + 2L
    if (any(left == substring(rc, rs, rs + k - 1L))) {
      return(TRUE)
    }
  }
  FALSE
}

# 4-bit base mask, A=1 C=2 G=4 T=8 (the bit order of `IUPAC_BASES`).  A code is
# the union of its bases, so "the base sets intersect" is a bitwise AND.  Any
# other character is -1, i.e. it pairs with nothing.  No names, so `[` lookup
# does not warn on a character outside the table.
IUPAC_MASKS <- c(
  A = 1L, C = 2L, G = 4L, T = 8L,
  R = 5L, Y = 10L, S = 6L, W = 9L, K = 12L, M = 3L,
  B = 14L, D = 13L, H = 11L, V = 7L, N = 15L
)

#' Base masks of a sequence, one per position
#' @param s a sequence of IUPAC codes, any case.
#' @return an integer vector; any character outside the table pairs with nothing.
#' @keywords internal
#' @noRd
iupac_masks <- function(s) {
  IUPAC_MASKS[toupper(stringr::str_split(s, "")[[1L]])]
}

#' Complement a sequence of IUPAC codes (A<->T, C<->G, R<->Y, ...)
#' @param s a sequence of IUPAC codes, any case.
#' @return `s` with each code replaced by its complement; unknown characters
#'   are left as they are.
#' @keywords internal
#' @noRd
iupac_complement <- function(s) {
  # `chartr` works on characters, so translate the upper form and keep case
  ch <- stringr::str_to_upper(s)
  comp <- chartr("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN", ch)
  out <- ifelse(ch == s, comp, stringr::str_to_lower(comp))
  paste0(out, collapse = "")
}

#' Does any wildcard-matched arm pair exist?
#'
#' Same offsets as `palindromic_arms()`: loop `d`, left arm `i..i+k-1`, right
#' arm `i+k+d..i+2k+d-1`.  Position `p` of each pairs when the two masks
#' intersect, i.e. some concrete base satisfies both codes.  The right arm is
#' read through the reverse complement, as `palindromic_arms()` does.
#' @param s a sequence containing at least one IUPAC code.
#' @param k arm length.
#' @param L maximum loop length.
#' @return a single logical.
#' @keywords internal
#' @noRd
wildcard_arms <- function(s, k, L) {
  n <- nchar(s)
  if (n < 2L * k) {
    return(FALSE)
  }
  # integer arithmetic: with `2L * k` double-typed the offsets become double and
  # `m[0 + p]` truncates to index 1, silently pairing the wrong positions
  k <- as.integer(k)
  m <- iupac_masks(s)
  # mask of the reverse complement, read with the same mirrored coordinates as
  # `palindromic_arms()`; `iupac_complement()` complements codes, not bases
  mr <- iupac_masks(intToUtf8(rev(utf8ToInt(iupac_complement(s)))))
  for (d in seq.int(0L, as.integer(L))) {
    npos <- n - 2L * k - d + 1L
    if (npos < 1L) break
    i <- seq_len(npos)
    ok <- rep(TRUE, npos)
    for (p in seq.int(0L, k - 1L)) {
      # left offset p against the reverse-complement arm start, rs above
      ok <- ok & bitwAnd(m[i + p], mr[n - i - 2L * k - d + 2L + p]) > 0L
      if (!any(ok)) break
    }
    if (any(ok)) {
      return(TRUE)
    }
  }
  FALSE
}

#' Nearest-neighbour Tm of oligos (perfect duplex), salt-corrected (TmCalculator)
#'
#' `res$gr` uses TmCalculator's `$` accessor for the result GRanges; `Tm` is a
#' metadata column on it.  `nn_table` selects the nearest-neighbour parameter set
#' (TmCalculator ships several: SantaLucia 2004 is the package default, and
#' Breslauer / Allawi / Sugimoto / Ghosh are available for an independent
#' thermodynamic cross-check, see [thermo_crosscheck()]).
#' @param seq oligo(s) 5'->3'.
#' @param DMSO DMSO concentration (%); it lowers Tm.  `NULL` (the default) uses
#'   `conditions$DMSO`, which is 0 for every scheme here: no design path in this
#'   package adds DMSO, including KASP.
#' @param nn_table nearest-neighbour parameter table (see details).
#' @param conditions reaction conditions from [reaction_conditions()] (or
#'   [scheme_conditions()], or a partial named list); its `Na`, `Mg`, `dNTPs` and
#'   `DMSO` fields are used.  `temperature` does not enter the nearest-neighbour
#'   Tm calculation (the NN model gives Tm directly), so passing
#'   `scheme_conditions("lamp")` changes the Tm only if that scheme's buffer
#'   differs from `PCR_COND`.  No scheme differs today, so it changes nothing;
#'   its `unmodelled` attribute names the additives neither model contains.
#'   `NULL` uses [reaction_conditions()].
#'   For an isothermal method whose buffer this model does not represent (RPA's
#'   crowded recombinase buffer) the Tm is informational, not a gate.
#' @param salt_method monovalent-salt correction passed to
#'   `TmCalculator::tm_calculate()`.  `NULL` uses the design coordinate,
#'   `TM_SALT_METHOD` (`"Schildkraut2010"`, Schildkraut & Lifson 1965,
#'   `16.6*log10[Na+]`).  The design windows (KASP 55-65/58-63 C, dPCR 44-50 C)
#'   are stated on that coordinate, so it does not move on its own; see ZW3.4.
#'   For a value comparable with Primer3/`oligotm` and Biopython, pass
#'   `"SantaLucia1998-2"`, the NN-consistent entropy
#'   correction of the same SantaLucia 1998 paper as the NN table; it reproduces
#'   a printed Primer3 panel to 0.08 C (ZW3.5).  `"Owczarzy2004"` /
#'   `"Owczarzy2008"` are available for a cross-check.  On the mm2020 measured
#'   set (50 mM Na, no Mg) the default reads 16.9 C low and the reference model
#'   9.6 C low at the design strand concentration; matching the set's own 1 uM
#'   total strand with `dnac` closes a further ~4.2 C, and the remainder is the
#'   fluorescent-label / first-derivative reference state, recorded and not
#'   removed (ZW3.4, `test-tm-mm2020.R`).
#' @param dnac strand concentration (nM, each strand) passed to
#'   `TmCalculator::tm_calculate()`.  `NULL` keeps the engine/design default of
#'   25 nM each (50 nM total), the Biopython `MeltingTemp` and Primer3 design
#'   convention.  Match a measured set's own strand concentration here when
#'   comparing against it: the dimerisation term makes Tm rise ~1.4 C per
#'   e-fold of concentration (ZW3.4).
#' @return a numeric vector of Tm (C), one per oligo.
#' @examples
#' primer_tm("GAAGGTGACCAAGTTCATGCT")
#' primer_tm(c("ACGTACGTACGTACGTACGT", "TTTTGGCGCGCCGATATATAT"))
#' primer_tm("GAAGGTGACCAAGTTCATGCT", conditions = reaction_conditions(Na = 100))
#' primer_tm("GAAGGTGACCAAGTTCATGCT", dnac = 500)
#' @family thermodynamics
#' @export
primer_tm <- function(seq, DMSO = NULL, nn_table = NN_TABLE,
                      conditions = NULL, salt_method = NULL, dnac = NULL) {
  conditions <- as_reaction_conditions(conditions)
  DMSO <- DMSO %||% conditions$DMSO
  seq <- as.character(seq)
  # Zero candidates is a normal state for a `mutate()` over an empty filter
  # (a very short template leaves design_common_primers() no window).  Passing
  # it on raised TmCalculator's "Input sequence cannot be NULL or empty" before
  # the caller could reach its own, clearer abort.
  if (!length(seq)) {
    return(numeric(0))
  }
  gr <- TmCalculator::to_genomic_ranges(seq)
  dnac <- dnac %||% formals(TmCalculator::tm_nn)$dnac_high
  res <- TmCalculator::tm_calculate(gr,
    method = "tm_nn", nn_table = nn_table,
    Na = conditions$Na, Mg = conditions$Mg, dNTPs = conditions$dNTPs, DMSO = DMSO,
    dnac_high = dnac, dnac_low = dnac,
    salt_method = salt_method %||% TM_SALT_METHOD
  )
  as.numeric(S4Vectors::mcols(res$gr)$Tm)
}

#' Is OligoArrayAux (`hybrid-min`) available?
#'
#' DECIPHER's thermodynamics runs through OligoArrayAux.  It is looked up as
#' `options(tidymarkers.oligoarrayaux)`, a project `.tools/`, the user cache (see
#' [marker_install_tools()]) and finally the `PATH`.
#'
#' This is the probe, not the contract: design code asks [capabilities()] /
#' `engine_has()` for the capability it needs, and only the install hint and the
#' absent-tool abort name this function.
#'
#' @return a single logical.
#' @examples
#' has_oligoarrayaux()
#' @family utilities
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
#' # the env var is restored when the expression finishes
#' with_unafloddat(tempdir(), Sys.getenv("UNAFOLDDAT"))
#' @family utilities
#' @export
with_unafloddat <- function(path, code) {
  if (!fs::dir_exists(path)) abort_msg("Directory not found: {path}")
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
  # duplex and hairpin through the configured engine: the internal adapters are
  # vectorised and need no subprocess (the hairpin one is a stem-only screen)
  engine <- getOption("tidymarkers.thermo_engine", "auto")
  if (identical(engine, "internal")) {
    if (is.null(seq2)) {
      return(internal_hairpin_dg(seq1, na = na, temp = temp))
    }
    return(internal_duplex_dg(seq1, seq2, na = na, temp = temp))
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
#' @family thermodynamics
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
#' @family thermodynamics
#' @export
duplex_dg_thermo <- function(seq1, seq2, conditions = NULL, BPPARAM = NULL) {
  conditions <- as_reaction_conditions(conditions)
  oligo_dg_vec(
    as.character(seq1), as.character(seq2),
    BPPARAM = BPPARAM,
    na = conditions$Na / 1000, temp = conditions$temperature
  )
}

#' Fold a whole amplicon/target, not just a primer
#'
#' Nothing is scored here: it returns the engine's raw answer plus the caveats a
#' scheme has to live with.  (1) The whole `seq` folds.  UNAFold's
#' `hybrid-ss-min` folds the pairwise part only, so a long two-helix hairpin
#' comes out higher than with the full base-pair-probability model, and the run
#' is reported negative: that is not a structural artefact.  (2) `kind` says
#' whether the fold was DNA or RNA.  The engine can fold both (`-n RNA`), but
#' this OligoArrayAux build has no RNA energy tables, so an RNA target is tried
#' as RNA, reported as `kind = "RNA"` when it works, and otherwise `unknown`;
#' it is never folded with the DNA tables - RNA:RNA G-U pairs would be dropped
#' and the structure would be wrong, not merely imprecise.  (3) `salt_modelled`
#' is `FALSE` for every RNA fold: UNAFold ignores `[Na+]`/`[Mg2+]` for RNA.
#'
#' an absent engine returns `status = "unknown"` (via [capabilities()]), never a
#' number.  Oligonucleotide concentration is not modelled, so this is the
#' intrinsic (temperature-independent) fold - use the DNA hairpin of
#' [hairpin_dg_thermo()] when the question is "does the primer itself fold at the
#' reaction temperature".  An amplicon that folds closed does not anneal a probe:
#' an MFE that keeps most of the target single-stranded is the signal for HRM.
#'
#' @param seq one amplicon/target, 5'->3'.
#' @param kind `"DNA"` or `"RNA"` sequence alphabet; default `"DNA"`.
#' @param conditions reaction conditions from [reaction_conditions()] or
#'   [scheme_conditions()]; `Na` (M = mM / 1000) and the temperature are what the
#'   engine is given.  `NULL` uses [reaction_conditions()].
#' @return a list: `status` (`"ok"`/`"unknown"`), `kind` (`"DNA"`/`"RNA"`),
#'   `length`, `mfe` (kcal/mol, `NA` unless `"ok"`), `structure` (dot-bracket,
#'   `NA` unless `"ok"`), `paired`, `single_stranded` (fractions), `salt_modelled`,
#'   `engine` and `note`.
#' @examples
#' \dontrun{
#' fold_amplicon(paste0("GCGCGCGCGCGGGGCCCCCAATT", "GCGCGCGCGCGGGGCCCCCAATT"))
#' fold_amplicon("ACGUACGUACGUGGCAUGCUAGC", kind = "RNA")
#' }
#' @seealso [hairpin_dg_thermo()], [oligo_dg()], [capabilities()]
#' @family thermodynamics
#' @export
fold_amplicon <- function(seq, kind = c("DNA", "RNA"), conditions = NULL) {
  kind <- match.arg(kind)
  conditions <- as_reaction_conditions(conditions)
  seq <- gsub("[[:space:]]", "", paste(as.character(seq), collapse = ""))
  if (!nzchar(seq)) abort_msg("`seq` must hold a sequence.")
  # `fold` is hybrid-ss-min here: hybrid-min folds pairs of strands, not one
  avail <- capabilities(marker_engines()$thermo, "fold")
  if (!isTRUE(avail)) {
    return(amplicon_fold_result(seq, kind,
      status = "unknown", engine = "OligoArrayAux",
      note = gettext("no folding engine (OligoArrayAux `hybrid-ss-min`): fold unknown")
    ))
  }
  dir <- withr::local_tempdir()
  # hybrid-ss-min writes <prefix>.ct (structure + dG) into the working directory,
  # and -q would skip the file and print 999.999 ("no structure") instead
  writeLines(c(">amplicon", seq), fs::path(dir, "fold.fa"))
  args <- c("-n", kind, "-t", conditions$temperature, "-T", conditions$temperature)
  # shortcut: UNAFold prints "salt concentrations ignored for RNA" and its `-q`
  # dG stays the only machine-readable output, so RNA goes without -N / -M.
  if (kind == "DNA") args <- c(args, "-N", conditions$Na / 1000, "-M", conditions$Mg / 1000)
  res <- withr::with_dir(dir, with_oligoarrayaux(system2(
    as.character(fs::path(oligoarrayaux_paths()$bin, "hybrid-ss-min")),
    c(args, "-o", "fold", "fold.fa"),
    stdout = TRUE, stderr = TRUE
  )))
  status <- attr(res, "status")
  ct <- fs::path(dir, "fold.ct")
  lines <- if (is.null(status) && fs::file_exists(ct)) readLines(ct, warn = FALSE) else character(0)
  mfe <- if (length(lines)) amplicon_mfe(lines) else NA_real_
  if (!is.finite(mfe)) {
    note <- if (kind == "RNA") {
      gettext("this OligoArrayAux build cannot fold RNA with its tables: RNA fold unknown (not folded as DNA)")
    } else {
      gettext("the folding engine returned no structure: fold unknown")
    }
    return(amplicon_fold_result(seq, kind, status = "unknown", engine = "OligoArrayAux", note = note))
  }
  amplicon_fold_result(seq, kind, mfe = mfe, structure = amplicon_structure(lines))
}

#' Assemble [fold_amplicon()]'s answer (single place for the shape)
#' @keywords internal
#' @noRd
amplicon_fold_result <- function(seq, kind, status = "ok", mfe = NA_real_,
                                 structure = NA_character_, engine = "OligoArrayAux",
                                 note = NA_character_) {
  ok <- identical(status, "ok")
  if (!ok) note <- paste0(note, "; the intrinsic dG is not known")
  list(
    status = status, kind = kind, length = nchar(seq),
    mfe = mfe, structure = structure,
    paired = if (ok) sum(strsplit(structure, "")[[1L]] != ".") / nchar(seq) else NA_real_,
    single_stranded = if (ok) sum(strsplit(structure, "")[[1L]] == ".") / nchar(seq) else NA_real_,
    salt_modelled = ok && kind == "DNA",
    engine = engine, note = note
  )
}

#' MFE dG (kcal/mol) from a UNAFold `.ct` file
#'
#' Line 1 is `N dG = <value> <name>`; `NA` when the value is the engine's
#' "no structure" sentinel or missing.
#' @keywords internal
#' @noRd
amplicon_mfe <- function(lines) {
  line <- lines[[1L]]
  value <- suppressWarnings(as.numeric(sub(".*dG[[:space:]]*=[[:space:]]*([-0-9.]+).*", "\\1", line)))
  if (!is.finite(value) || value > 100) {
    return(NA_real_)
  }
  value
}

#' Dot-bracket from a UNAFold `.ct` file
#'
#' The partner of base `i` is column 5 of its row, and bases are numbered from 1.
#' @keywords internal
#' @noRd
amplicon_structure <- function(lines) {
  # line 1 is "N dG = x name"; a base row starts with its index, not a length
  rows <- strsplit(trimws(lines[-1L]), "[[:space:]]+")
  n <- length(rows)
  pair <- if (n) vapply(rows, function(r) as.integer(r[[5L]]), integer(1)) else integer(0)
  if (n && (length(pair) != n || anyNA(pair) || any(pair > n))) {
    return(NA_character_)
  }
  s <- rep(".", n)
  for (i in seq_along(pair)) {
    if (pair[[i]] > i) {
      s[[i]] <- "("
      s[[pair[[i]]]] <- ")"
    }
  }
  paste(s, collapse = "")
}

# Z10.1: `pass` above this share of the amplicon left single-stranded.  It is a
# readout heuristic (structure, not energy), not a measured threshold: an MFE
# that keeps most of the target open is the HRM/reverse-transcription signal.
AMPLICON_FOLD_MIN_SS <- 0.5

#' Amplicon/target folding check (whole amplicon, not the primer)
#'
#' An opt-in structure check that folds the WHOLE amplicon with [fold_amplicon()]
#' and grades the intrinsic fold as `pass`/`warn`/`unknown`.  It is additive: no
#' existing gate reads it, and no gate result moves.
#'
#' **A primer's folding dG is not an amplicon's.**  [hairpin_dg_thermo()] (and
#' [oligo_dg()]) fold a short oligo at the reaction temperature: the question is
#' whether that oligo closes on itself before it anneals.  This check folds the
#' product, which is longer, has extra pairing partners, and is scored at the
#' engine's default temperature as an intrinsic (concentration-free) MFE.  The
#' two numbers are not interchangeable, and the oligo dG is not a safe proxy for
#' the amplicon.
#'
#' **What stays `unknown`.**  (1) An RNA target (`kind = "RNA"`, or a chemistry
#' whose target is RNA - Cas13, or a LAMP/RPA assay run on an RNA template): the
#' package's DNA hairpin engine is not an RNA folding model, and a DNA number
#' would drop G-U pairs, so the answer is `unknown` unless [fold_amplicon()]
#' itself reports an RNA fold.  (2) No engine: `capabilities("fold")` is `FALSE`
#' or the engine call fails.  (3) A long amplicon (sanger, 400-800 bp): UNAFold's
#' folding is quadratic-ish in length, so it is not run past `max_length` and is
#' reported `unknown` rather than scored - LAMP (120-300 bp) is the longest
#' chemistry that is scored here.
#'
#' @param seq one amplicon/target, 5'->3' (a vector is folded element-wise).
#' @param kind `"DNA"` or `"RNA"` alphabet of the target.  A chemistry whose
#'   target is RNA is the caller's to declare; it is never inferred.
#' @param conditions reaction conditions from [reaction_conditions()] or
#'   [scheme_conditions()]; `NULL` uses [reaction_conditions()].
#' @param min_single_stranded the `single_stranded` fraction at which the fold
#'   passes; below it the status is `"warn"`.
#' @param max_length longest amplicon folded; longer is `unknown`.
#' @return a tibble, one row per `seq`: `length`, `kind`, `status`
#'   (`"pass"`/`"warn"`/`"unknown"`), `mfe` (kcal/mol, `NA` unless scored),
#'   `single_stranded`, `engine` and `note` (why it is `unknown`/`warn`).
#' @examples
#' amplicon_fold_check(paste0("GCGCGCGCGCGGGGCCCCCAATT", "GCGCGCGCGCGGGGCCCCCAATT"))
#' amplicon_fold_check("ACGUACGUACGUGGCAUGCUAGC", kind = "RNA")
#' @seealso [fold_amplicon()], [hairpin_dg_thermo()], [capabilities()]
#' @family thermodynamics
#' @export
amplicon_fold_check <- function(seq, kind = c("DNA", "RNA"), conditions = NULL,
                                min_single_stranded = AMPLICON_FOLD_MIN_SS,
                                max_length = 300L) {
  kind <- match.arg(kind)
  seq <- gsub("[[:space:]]", "", as.character(seq))
  if (!length(seq) || !all(nzchar(seq))) abort_msg("`seq` must hold a sequence.")
  purrr::map(seq, \(s) amplicon_fold_row(s, kind, conditions, min_single_stranded, max_length)) |>
    purrr::list_rbind()
}

#' One [amplicon_fold_check()] row, or the reason it is `unknown`
#' @keywords internal
#' @noRd
amplicon_fold_row <- function(seq, kind, conditions, min_single_stranded, max_length) {
  unknown <- function(note, engine = "OligoArrayAux") {
    tibble::tibble(
      length = nchar(seq), kind = kind, status = "unknown", mfe = NA_real_,
      single_stranded = NA_real_, engine = engine, note = note
    )
  }
  if (!isTRUE(capabilities(marker_engines()$thermo, "fold"))) {
    return(unknown(gettext("no folding engine (OligoArrayAux `hybrid-ss-min`): fold unknown")))
  }
  if (nchar(seq) > max_length) {
    return(unknown(gettextf(
      "amplicon longer than %d bp is not folded here: fold unknown", max_length
    )))
  }
  # the engine call itself may fail (a build without the RNA tables reports that
  # through `status = "unknown"`; an I/O failure is caught the same way)
  res <- tryCatch(
    fold_amplicon(seq, kind = kind, conditions = conditions),
    error = function(e) list(status = "unknown", note = conditionMessage(e))
  )
  if (!identical(res$status, "ok")) {
    note <- res$note
    if (!is.character(note) || length(note) != 1L || is.na(note)) {
      note <- gettext("the folding engine returned no structure: fold unknown")
    }
    return(unknown(note))
  }
  # an RNA fold is only reported as RNA; `salt_modelled` is FALSE for RNA by
  # construction, so a DNA fold is the only one that carries a salted number
  ss <- res$single_stranded
  pass <- is.finite(ss) && ss >= min_single_stranded
  tibble::tibble(
    length = res$length, kind = kind, status = if (pass) "pass" else "warn",
    mfe = res$mfe, single_stranded = ss, engine = res$engine,
    note = if (pass) {
      NA_character_
    } else {
      gettextf("most of the amplicon is paired (%.0f%% single-stranded)", 100 * ss)
    }
  )
}

#' PCR amplification efficiency of primer/target pairs (DECIPHER)
#'
#' Models hybridisation (incl. side reactions) times Taq elongation efficiency,
#' so a 3'-terminal mismatch is strongly penalised - the physical basis of KASP.
#' Needs OligoArrayAux through DECIPHER, i.e. the `pcr_efficiency` capability of
#' the thermo adapter ([marker_engines()]); callers gate on that, not on the tool.
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
#' @family thermodynamics
#' @export
pcr_efficiency <- function(primers, targets, temp = PCR_EFF$temp,
                           P = PCR_EFF$P, ions = PCR_EFF$ions,
                           cores = getOption("tidymarkers.cores", 1L)) {
  rlang::check_installed("DECIPHER")
  if (!has_oligoarrayaux()) {
    abort_msg("OligoArrayAux ('hybrid-min') is required for PCR efficiencies.",
      i = "run marker_install_tools(), or build it with data-raw/install-oligoarrayaux.R"
    )
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
