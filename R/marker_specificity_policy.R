# =============================================================================
#  marker_specificity_policy.R  --  chemistry tolerance and scan grading
# =============================================================================
#  Split out of marker_specificity.R (T7, cycle2/tidy-20261126).  This is the
#  contract side of a specificity scan: what each chemistry tolerates
#  (`specificity_policy()`), which tier answered (`specificity_backend()`), the
#  PAM rule, and how a hit becomes a verdict (`specificity_verdict()`).
#  Function names are unchanged by the split.
# =============================================================================

# =============================================================================
#  1. Scheme-specific mismatch tolerance (heuristic, position-aware)
# =============================================================================

# The tolerance table.  One row per chemistry; `specificity_policy()` returns
# the row as an object.  Every number is a heuristic (calibrated = FALSE) taken
# from the thermodynamics the markers rely on, not fitted to data:
#
#  * pcr/KASP: polymerase extension starts at the 3' terminus, so a 3'-terminal
#    mismatch is not a "mismatch" but loss of extension.  A mismatch in the last
#    ~5 nt is destabilising enough to swallow ordinary allele discrimination, so
#    it is counted separately and capped at 0.  One mismatch elsewhere (the
#    deliberate-mismatch and 5'/middle positions) is tolerated.
#    KASP allele-specific PCR relies on the AS primer's 3' base, so an
#    exact-hit-only screen and a 15-nt 3' seed both miss productive off-targets
#    when the mismatch sits in the middle - hence `max_mismatch = 1` here rather
#    than 0.
#  * rpa: recombinase filaments tolerate mismatches and strand exchange is not
#    proofreading, so several mismatches (up to ~3-5) can still amplify - an
#    exact-hit-only in-silico RPA is wrong.  The 3' end still matters (the
#    strand-displacement polymerase needs a paired terminus) but less strictly
#    than for PCR: one mismatch in the last 5 nt, no terminal mismatch.
#  * crispr: a Cas12a guide's PAM-proximal ~8-12 nt seed is the mismatch-critical
#    region (0-1 seed mismatch is tolerated); Cas13's critical region is central,
#    reported in `critical_from`/`critical_to` for clarity.  A changed PAM is a
#    HARD loss of target (a mismatch in the PAM is not a tolerated primer
#    mismatch and there is no partial credit).  Cas13 has no DNA PAM but a weak
#    PFS preference, carried in `pfs = "weak"` - Cas13 is not PAM-free.
#  * paralog/homeolog: no exact-hit assumption at all - a homeolog may differ by
#    several bases yet still be a source of cross-hybridisation - so it is the
#    most permissive DNA mode, with the 3' end and seed still capped.
SPECIFICITY_POLICIES <- list(
  pcr = list(
    mode = "pcr", max_mismatch = 1L, seed_len = 5L, seed_max = 0L,
    three_prime_max = 0L, requires_pam = FALSE, pam = NA_character_, pfs = NA_character_,
    weight_middle = 1, weight_seed = 3, weight_three_prime = 8, weight_terminal = 12
  ),
  rpa = list(
    mode = "rpa", max_mismatch = 4L, seed_len = 5L, seed_max = 2L,
    three_prime_max = 1L, requires_pam = FALSE, pam = NA_character_, pfs = NA_character_,
    weight_middle = 1, weight_seed = 2, weight_three_prime = 5, weight_terminal = 12
  ),
  crispr = list(
    mode = "crispr", max_mismatch = 5L, seed_len = 8L, seed_max = 1L,
    three_prime_max = 2L, requires_pam = TRUE, pam = "TTTV", pfs = "weak",
    critical_from = 1L, critical_to = 8L,
    weight_middle = 1, weight_seed = 4, weight_three_prime = 8, weight_terminal = 12,
    spacer_len_max_mismatch = c(`16` = 0L, `20` = 1L, `23` = 2L, `28` = 0L),
    spacer_len_central = list(`20` = c(7L, 9L))
  ),
  paralog = list(
    mode = "paralog", max_mismatch = 4L, seed_len = 8L, seed_max = 1L,
    three_prime_max = 1L, requires_pam = FALSE, pam = NA_character_, pfs = NA_character_,
    weight_middle = 1, weight_seed = 2, weight_three_prime = 5, weight_terminal = 12
  )
)

# Chemistry aliases -> a policy name (the scheme a marker object would carry)
SPECIFICITY_SCHEME_ALIASES <- c(
  kap = "pcr", kas = "pcr", kasp = "pcr", tas = "pcr", arms = "pcr",
  caps = "pcr", scar = "pcr", ssr = "pcr", dcap = "pcr", dpcr = "pcr",
  as = "pcr", hrm = "pcr", pcrr = "pcr",
  rpa = "rpa", lamp = "rpa", aslamp = "rpa",
  cas9 = "crispr", cas12 = "crispr", cas12a = "crispr", cas13 = "crispr",
  cas13a = "crispr", cas13b = "crispr", sherlock = "crispr", detectr = "crispr",
  paralog = "paralog", homeolog = "paralog", homeologue = "paralog",
  # a homoeolog is the same concept: a duplicated subgenome copy, not a gene
  # family member.  Spellings follow the literature (wheat papers use
  # "homoeolog"); the tables keep the shorter `homeolog` column name.
  homoeolog = "paralog", homoeologue = "paralog",
  homolog = "paralog", homologe = "paralog"
)

#' Per-position mismatch weight for a primer hit (heuristic, not calibrated)
#'
#' A mismatch near the primer's 3' terminus destabilises the annealed end far
#' more than one in the middle, because extension starts there; the seed region
#' (the 3'-terminal `seed_len` nt, the mismatch-critical region of an
#' allele-specific primer) sits just behind it.  The ladder is the scheme's own
#' (`weight_seed`, `weight_three_prime`, `weight_terminal` of
#' [specificity_policy()]), so RPA's tolerant filaments are not scored like a
#' proofreading polymerase.  It is a monotone heuristic, not a free-energy
#' calculation.
#' @param n_nt primer length.
#' @param scheme chemistry (see [specificity_policy()]).
#' @param seed_len seed length (default: the policy's).
#' @return numeric vector of length `n_nt`, weight per position (5'->3').
#' @keywords internal
#' @noRd
mismatch_weights <- function(n_nt, scheme = "pcr", seed_len = NULL) {
  pol <- specificity_policy(scheme)
  seed_len <- seed_len %||% pol$seed_len
  n_nt <- as.integer(n_nt)
  if (!n_nt) {
    return(numeric())
  }
  w <- rep(pol$weight_middle, n_nt)
  seed <- seq.int(max(1L, n_nt - seed_len + 1L), n_nt)
  w[seed] <- pol$weight_seed
  tail5 <- seq.int(max(1L, n_nt - pol$three_prime_window + 1L), n_nt)
  w[tail5] <- pol$weight_three_prime
  w[tail_terminal(n_nt)] <- pol$weight_terminal
  w
}

#' Last `THREE_PRIME_TERMINAL` positions of `n_nt` (5'->3', 1-based)
#' @keywords internal
#' @noRd
tail_terminal <- function(n_nt) {
  seq.int(max(1L, n_nt - THREE_PRIME_TERMINAL + 1L), n_nt)
}

#' Mismatch tolerance for a specificity scheme
#'
#' Encodes the *known biology* each chemistry imposes on a target site: how many
#' mismatches can still amplify, where they may sit, and whether a licensing
#' element (a CRISPR PAM) must be intact.  The tolerance numbers are
#' **heuristics** (`calibrated = FALSE`): they follow the qualitative rule the
#' literature states, they are not fitted to an assay dataset.
#'
#' | scheme | max_mismatch | seed_len | seed_max | three_prime_max | requires_pam |
#' |---|---|---|---|---|---|
#' | `pcr` (KASP/ARMS/CAPS/SSR/dPCR) | 1 | 5 | 0 | 0 | no |
#' | `rpa` (RPA/LAMP) | 4 | 5 | 2 | 1 | no |
#' | `crispr` (Cas9/Cas12/Cas13) | 5 | 8 | 1 | 2 | yes (Cas12a `TTTV`) |
#' | `paralog` (homeolog/homolog) | 4 | 8 | 1 | 1 | no |
#'
#' ## CRISPR spacer length
#'
#' The `crispr` count above is spacer-independent and stays the fallback: an
#' effector's mismatch budget also depends on the guide's length, and the
#' published series gives an endpoint per length.
#'
#' | spacer_len (nt) | max_mismatch | evidence |
#' |---|---|---|
#' | 16 | 0 | a 16-nt LbuCas13a spacer still required a perfect match (PMC10810210) |
#' | 20 | 1 | a single middle mismatch (positions +7 to +9) costs up to 70 % activity (PMC10810210) |
#' | 23 | 2 | 2 central mismatches still silence the K-RAS G12D transcript 54.8 %; 3 leave 11 % (PMC11654686) |
#' | 28 | 0 | a 28-nt LbuCas13a spacer has no single-mismatch discriminatory capacity (PMC10810210) |
#'
#' The 20-nt row is stricter than "anywhere": its one tolerated mismatch must
#' sit outside spacer positions 7-9 (1-based from the 5' end), the middle window
#' the paper measured.  These numbers stay heuristics (`calibrated = FALSE`):
#' they are published endpoints, not a fit.
#'
#' Position ladder (per-position weight of a mismatch, 5'->3' of the primer):
#'
#' | scheme | `weight_middle` | `weight_seed` | `weight_three_prime` | `weight_terminal` |
#' |---|---|---|---|---|
#' | `pcr` | 1 | 3 | 8 | 12 |
#' | `rpa` | 1 | 2 | 5 | 12 |
#' | `crispr` | 1 | 4 | 8 | 12 |
#' | `paralog` | 1 | 2 | 5 | 12 |
#'
#' `weight_seed` applies to the 3'-terminal `seed_len` nt, `weight_three_prime`
#' to the last `three_prime_window` (5) nt and `weight_terminal` to the last 3 nt
#' - the 3'-terminal mismatch that stops extension.  RPA/LAMP tolerate several
#' mismatches, so a middle mismatch costs less there than for a proofreading
#' polymerase; the terminal weight is the same everywhere, because no chemistry
#' extends from a mismatched 3' end.  The pcr row is the historical scoring, so
#' `pcr`/`kasp` callers' scores are unchanged.
#'
#' Biological basis:
#' * **PCR/KASP** - extension begins at the 3' terminus: a 3'-terminal mismatch
#'   is a hard loss (never a tolerated mismatch), 0 mismatches are allowed in the
#'   last ~5 nt, and 0-1 outside them.  An exact-hit-only screen is not
#'   sufficient, because a middle mismatch can still amplify.
#' * **RPA** - the recombinase filament and strand-displacement polymerase are
#'   mismatch-tolerant: several mismatches (up to ~3-5) can still amplify, so
#'   exact-hit-only is wrong.  The 3' end still matters, so at most one mismatch
#'   in the last 5 nt and none terminal.
#' * **CRISPR** - the Cas12a seed is PAM-proximal (~8-12 nt from the PAM): 0-1
#'   seed mismatch, and a changed PAM is a hard loss of target, not a tolerated
#'   mismatch.  Cas13 has no DNA PAM, but its PFS is a weak preference, carried
#'   in `pfs = "weak"` - Cas13 is never PAM-free.
#' * **paralog/homeolog** - mismatch-aware by construction: a homeolog differs
#'   from the assayed copy by several bases and still cross-hybridises, so an
#'   exact-hit assumption understates off-targets.
#'
#' @param scheme a scheme name (`"pcr"`, `"rpa"`, `"crispr"`, `"paralog"`) or a
#'   chemistry alias (`"kasp"`, `"lamp"`, `"cas12a"`, `"homeolog"`, ...).
#' @param length optional target-site length; when given, the returned object
#'   carries a named `weight` vector with the per-position ladder (5'->3').
#' @param spacer_len spacer (guide) length in nt for `"crispr"`; one of the
#'   published series (`16`, `20`, `23`, `28`) sets the mismatch budget and, for
#'   `20`, an intolerant middle window.  An unlisted length keeps the scheme's
#'   `max_mismatch` and reports `spacer_len_known = FALSE`; `NULL` keeps the
#'   length-independent policy.
#' @return a `specificity_policy` list: `mode`, `max_mismatch`, `seed_len`,
#'   `seed_max`, `three_prime_max`, `three_prime_window`, `requires_pam`, `pam`,
#'   `pfs`, `critical_from`/`critical_to` (CRISPR), the position ladder
#'   (`weight_middle`, `weight_seed`, `weight_three_prime`, `weight_terminal`),
#'   `priors` (the tolerance row as a one-row tibble) and `calibrated` (`FALSE`).
#'   A CRISPR policy given a `spacer_len` also carries `spacer_len`,
#'   `spacer_len_known` and `spacer_len_central` (the intolerant middle window,
#'   or `NULL` when that length has none).
#' @examples
#' specificity_policy("kasp")$max_mismatch
#' specificity_policy("rpa")$max_mismatch
#' specificity_policy("cas12a")$requires_pam
#' specificity_policy("homeolog")$mode
#' @family specificity
#' @export
specificity_policy <- function(scheme = "pcr", length = NULL, spacer_len = NULL) {
  if (is.null(scheme)) {
    scheme <- "pcr"
  }
  if (!is.character(scheme) || length(scheme) != 1L) {
    abort_msg("`scheme` must be a single scheme name.")
  }
  if (!is.null(spacer_len)) {
    n <- if (is.numeric(spacer_len)) spacer_len else NA_real_
    whole <- length(n) == 1L && is.finite(n) && n >= 1 &&
      n <= .Machine$integer.max && n == round(n)
    if (!whole) {
      abort_msg(
        "`spacer_len` must be a single positive whole number of nucleotides."
      )
    }
  }
  key <- gsub("[-_ ]", "", tolower(scheme))
  alias <- match(key, names(SPECIFICITY_SCHEME_ALIASES))
  key <- if (!is.na(alias)) {
    unname(SPECIFICITY_SCHEME_ALIASES[[alias]])
  } else if (key %in% names(SPECIFICITY_POLICIES)) {
    key
  } else {
    cli::cli_abort(c(
      gettext("Unknown specificity scheme {.val {scheme}}."),
      i = gettext(paste(
        "Use one of {.val {names(SPECIFICITY_POLICIES)}} or an alias",
        "({.val {names(SPECIFICITY_SCHEME_ALIASES)}})."
      ))
    ))
  }
  pol <- SPECIFICITY_POLICIES[[key]]
  # a published spacer length overrides the scheme budget; an unlisted length
  # keeps it and says so, so a length the evidence does not cover is never
  # guessed at, and `priors` below then carries the overriding number
  if (identical(key, "crispr") && !is.null(spacer_len)) {
    pol$spacer_len <- as.integer(spacer_len)
    name <- as.character(pol$spacer_len)
    if (name %in% names(pol[["spacer_len_max_mismatch"]])) {
      pol$max_mismatch <- pol[["spacer_len_max_mismatch"]][[name]]
      pol$spacer_len_known <- TRUE
    } else {
      pol$spacer_len_known <- FALSE
    }
    pol$spacer_len_central <- if (!is.null(pol[["spacer_len_central"]])) {
      pol[["spacer_len_central"]][[name]]
    } else {
      NULL
    }
  }
  pol$scheme <- scheme
  pol$three_prime_window <- THREE_PRIME_WINDOW
  pol$calibrated <- FALSE
  pol$priors <- tibble::tibble(
    mode = pol$mode, max_mismatch = pol$max_mismatch, seed_len = pol$seed_len,
    seed_max = pol$seed_max, three_prime_max = pol$three_prime_max,
    requires_pam = pol$requires_pam, calibrated = FALSE
  )
  if (!is.null(length)) {
    pol$weight <- mismatch_weights(as.integer(length), pol$mode, pol$seed_len)
  }
  structure(pol, class = c("specificity_policy", "list"))
}

#' @export
print.specificity_policy <- function(x, ...) {
  cat("<specificity_policy:", x$mode, "(heuristic, calibrated =", x$calibrated, ")>\n")
  cat("  max_mismatch", x$max_mismatch, "| seed", x$seed_max, "of", x$seed_len,
    "| last", x$three_prime_window, "nt", x$three_prime_max, "\n",
    sep = " "
  )
  cat("  weights: middle", x$weight_middle, "| seed", x$weight_seed,
    "| last", x$three_prime_window, "nt", x$weight_three_prime,
    "| terminal", x$weight_terminal, "\n",
    sep = " "
  )
  if (!is.null(x[["spacer_len"]])) {
    cat("  spacer", x[["spacer_len"]], "nt | max mismatch", x$max_mismatch,
      if (!is.null(x[["spacer_len_central"]])) {
        c(
          "| no mismatch in positions",
          paste0(x[["spacer_len_central"]], collapse = "-")
        )
      },
      "\n",
      sep = " "
    )
  }
  if (x$requires_pam) cat("  PAM", x$pam, "required (a changed PAM = loss of target)\n")
  if (!is.na(x$pfs)) cat("  PFS preference:", x$pfs, "\n")
  invisible(x)
}

#' Is there a genome/specificity engine able to answer for this genome?
#'
#' The in-memory scan needs a `DNAStringSet` (or an indexed `marker_genome` read
#' through Rsamtools).  An indexed genome without a BLAST database on this
#' machine has **no** genome-wide engine, and that is reported as missing - the
#' caller must not read it as a pass.  A Bowtie/BWA/crisprVerse backend is not
#' implemented here; when one is needed (CRISPR genome-wide), it is reported
#' missing too, with its capability flag.
#' @param genome a `DNAStringSet`, a `marker_genome`, or `NULL`.
#' @param scheme specificity scheme ([specificity_policy()]).
#' @return `TRUE`/`FALSE` for `engine` (a usable in-memory/indexed scan) and
#'   named capability flags (`exact_scan`, `mismatch_scan`, `blast`,
#'   `crispr_offtarget`, `bowtie`, `bwa`, `crisprverse`).
#' @keywords internal
#' @noRd
specificity_backend <- function(genome = NULL, scheme = "pcr") {
  pol <- specificity_policy(scheme)
  in_memory <- inherits(genome, "DNAStringSet")
  blast_ok <- has_blast() && isTRUE(requireNamespace("rBLAST", quietly = TRUE))
  if (inherits(genome, "genome_index")) {
    # the declared strategy is the answer: a capability the strategy does not
    # have is a missing engine, never an inferred one
    caps <- genome$capabilities
    return(list(
      engine = switch(genome$strategy,
        memory = TRUE,
        bowtie = "index",
        blast = "blast",
        fm = "fm"
      ),
      capabilities = c(
        exact_scan = isTRUE(caps$seed_scan),
        mismatch_scan = isTRUE(caps$mismatch_tolerant),
        blast = blast_ok,
        crispr_offtarget = FALSE,
        bowtie = isTRUE(caps$whole_primer_scan),
        bwa = nzchar(Sys.which("bwa")),
        crisprverse = namespace_available("crisprVerse")
      ),
      scheme = pol$mode,
      reason = NA_character_
    ))
  }
  if (inherits(genome, "marker_genome")) {
    in_memory <- !is.null(genome$blast_db) && blast_ok
  }
  caps <- c(
    exact_scan = in_memory,
    mismatch_scan = in_memory,
    blast = blast_ok,
    crispr_offtarget = FALSE,
    bowtie = nzchar(Sys.which("bowtie")),
    bwa = nzchar(Sys.which("bwa")),
    crisprverse = namespace_available("crisprVerse")
  )
  list(
    engine = in_memory,
    capabilities = caps,
    scheme = pol$mode,
    reason = if (in_memory) {
      NA_character_
    } else if (is.null(genome)) {
      gettext("No genome was supplied.")
    } else {
      gettext(paste(
        "No genome-wide engine for this genome: it is not an in-memory",
        "DNAStringSet and it has no usable BLAST database."
      ))
    }
  )
}

# Six-frame translation need not be rigorous here: a protein PAM is a 3-codon
# string.  (The guide-level PAM check uses `pam_state()` below.)

#' Is a CRISPR PAM/PFS intact?
#'
#' `"match"` when a PAM matching the scheme's pattern (case-insensitive, IUPAC
#' `V`/`D` expanded) was observed; `"pam_loss"` when it does not - a hard loss of
#' target, never a tolerated primer mismatch; `"absent"` when no PAM was
#' supplied.  `"PFS"` (a Cas13 guide) is checked only for its weak documented
#' rule: a `G` immediately 3' of the protospacer is disfavoured, any other base
#' is `"match"`.  That weak PFS is a preference, not a PAM, and is never
#' reported as "PAM-free".
#' @param observed the PAM/PFS sequence read from the target, or `"PFS"`.
#' @param pattern the licensed PAM for the effector (default Cas12a's `TTTV`).
#' @keywords internal
#' @noRd
pam_state <- function(observed, pattern = "TTTV") {
  if (is.null(observed) || is.na(observed) || !nzchar(observed)) {
    return("absent")
  }
  obs <- toupper(observed)
  if (identical(obs, "PFS")) {
    return("match")
  }
  if (nchar(obs) != nchar(pattern)) {
    return("pam_loss")
  }
  ok <- mapply(function(code, base) iupac_covers(code, base),
    strsplit(pattern, "", fixed = TRUE)[[1L]],
    strsplit(obs, "", fixed = TRUE)[[1L]],
    USE.NAMES = FALSE
  )
  if (all(ok)) "match" else "pam_loss"
}

#' Grade in-silico PCR products against a scheme's mismatch policy
#'
#' The structured answer to "is this specific enough?": `pass` when at least one
#' product is tolerated *and* it is the only tolerated product;
#' `fail` when every product breaks the policy (including a changed PAM),
#' more than one tolerated product exists, or the backend is primer-level; `unknown`
#' when there is no genome-wide engine - a missing backend is never a silent pass.
#' The result carries the products, the tolerated flag, the best product and the
#' policy it was graded against, so the evidence travels with the verdict, plus
#' `tier`, the evidence level behind the verdict:
#'
#' 1. `"pair"` - pair-level in-silico PCR with product-size bounds, graded here;
#' 2. `"blast"` (`"fm"`) - a primer-level BLAST database / FM-index scan.  The
#'    per-product `n_mismatch` is this pair's primer evidence, so at primer level
#'    it cannot stand in for a genome-wide off-target count: the verdict grades
#'    product formation, and a changed PAM still fails;
#' 3. `"none"` - nothing genome-wide was consulted, so the verdict is `unknown`.
#'
#' [specificity_seed()] never produces a verdict: it is a prefilter, and its
#' fixed false-negative rate is not genome-wide evidence.
#' @param hits in-silico PCR products ([in_silico_pcr()], as a tibble: one row
#'   per product, `chrom`/`start`/`end`/`width` and the mismatch evidence;
#'   extra columns - e.g. a `guide` id from a multi-guide table - are ignored.
#'   An optional `scheme`/`tolerated` pair is recomputed here, so a table from
#'   `in_silico_pcr(scheme = )` is graded exactly as a bare one.
#' @param scheme specificity scheme ([specificity_policy()]).
#' @param genome optional genome; when supplied, a missing engine turns the
#'   answer into `unknown` (see `specificity_backend()`).
#' @param pam optional PAM/PFS observation for `"crispr"` (see `pam_state()`).
#' @param spacer_len optional spacer (guide) length in nt for `"crispr"`; a
#'   length in the published series sets the mismatch budget and, for `20`, the
#'   intolerant middle window (see [specificity_policy()]).
#' @return a `specificity_verdict` list: `verdict` (`"pass"`/`"fail"`/
#'   `"unknown"`), `reason` (an English sentence), `scheme`, `policy` (the policy
#'   graded against, carrying `spacer_len`/`spacer_len_known` for a CRISPR
#'   guide), `n_products`, `n_tolerated`, `n_mismatch` (best product), `best`
#'   (one-row product tibble), `products` (the graded tibble), `pam` (state),
#'   `backend` (the capability list, or `NA` when no genome was supplied) and
#'   `tier` (`"pair"`/`"index"`/`"blast"`/`"fm"`/`"none"`).
#' @examples
#' genome <- Biostrings::DNAStringSet(c(chr1 = paste0(
#'   strrep("A", 20L), "ACGTACGTACGTACGTACGT", strrep("T", 20L),
#'   "ACGTACGTACGTACGTACGT", strrep("A", 20L)
#' )))
#' hits <- in_silico_pcr("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT", genome)
#' # `unknown` when the genome is missing: a missing engine is never a pass
#' specificity_verdict(hits, "kasp", genome = NULL)$verdict
#' specificity_verdict(hits, "kasp", genome = NULL)$tier
#' @family specificity
#' @export
specificity_verdict <- function(hits, scheme = "pcr", genome = NULL, pam = NULL,
                                spacer_len = NULL) {
  pol <- specificity_policy(scheme, spacer_len = spacer_len)
  backend <- if (is.null(genome)) NA else specificity_backend(genome, scheme)
  # a table may answer one question per genome, never two for the same genome:
  # a count cannot be both a seed count and a whole-primer count
  if (!is.null(attr(hits, "scan"))) {
    q <- scan_question(hits)
    by_genome <- split(unname(q), names(q) %||% rep("", length(q)))
    mixed <- any(vapply(by_genome, function(v) length(unique(v)) > 1L, logical(1)))
    if (mixed) {
      abort_msg(
        "This table mixes \"3'-seed\" and \"whole-primer\" counts for one genome.",
        i = "Keep one question per genome: the two counts are not comparable.",
        class = "tidymarkers_scan_mixed_error"
      )
    }
  }
  # one verdict answers for one build: a table that spans builds has no single
  # meaning, so refuse it instead of pooling A/B/D
  seen <- if (!is.null(hits) && "build" %in% names(hits)) {
    unique(as.character(hits$build))
  } else {
    as.character(attr(hits, "genome_build") %||% character(0))
  }
  seen <- seen[!is.na(seen)]
  if (length(seen) > 1L) {
    abort_msg(
      "This table mixes builds {.val {seen}}; a specificity verdict answers for one build.",
      i = "Grade each build separately: one specificity_of() call and one verdict per build.",
      class = "tidymarkers_mixed_build_error"
    )
  }
  products <- if (is.null(hits)) {
    empty_pcr_products()
  } else {
    tibble::as_tibble(hits)[intersect(names(empty_pcr_products()), names(hits))]
  }
  pam_given <- !is.null(pam)
  if (!all(c("n_mismatch", "seed_mismatch", "three_prime_mismatch") %in% names(products))) {
    # a legacy exact table: it carries no evidence, so it cannot be graded
    products$n_mismatch <- 0L
    products$seed_mismatch <- 0L
    products$three_prime_mismatch <- 0L
    products$three_prime_terminal <- FALSE
  }
  if (!"penalty" %in% names(products)) products$penalty <- 0
  # A CRISPR PAM is licensing, not a mismatch: when it was observed, a changed
  # PAM removes the tolerated status of every product, even a perfect one.
  pam_val <- if (pam_given && identical(pol$mode, "crispr")) {
    pam_state(pam, pol$pam %||% "TTTV")
  } else {
    NA_character_
  }
  central_bad <- central_violated(products, pol)
  products$tolerated <- pcr_products_tolerated(products, pol) & !central_bad
  products$violations <- mapply(pcr_violations,
    n_mismatch = products$n_mismatch, seed_mismatch = products$seed_mismatch,
    three_prime_mismatch = products$three_prime_mismatch,
    three_prime_terminal = products$three_prime_terminal,
    MoreArgs = list(policy = pol), USE.NAMES = FALSE
  )
  if (any(central_bad)) {
    products$violations[central_bad] <- lapply(
      products$violations[central_bad], function(v) union(v, "central")
    )
  }
  if (identical(pam_val, "pam_loss")) {
    products$tolerated <- rep(FALSE, nrow(products))
    products$violations <- lapply(products$violations, function(v) union(v, "pam"))
  }
  # The best evidence: a tolerated product first, then the smallest penalty.  At
  # equal penalty the complete (widest) product wins over an end-padded one, so
  # the scan radius cannot manufacture the headline number.
  best <- if (nrow(products)) {
    products[order(
      !products$tolerated, products$penalty,
      -products$width, products$start
    ), ][1L, ]
  } else {
    NULL
  }
  n_tol <- sum(products$tolerated)
  # the evidence level: a genome-wide index (BLAST/FM/bowtie) is primer-level
  # (it can fail a primer genome-wide), an in-memory genome is pair-level here,
  # and no engine at all is no tier
  engine <- if (is.list(backend)) backend$engine else FALSE
  tier <- if (identical(engine, "blast")) {
    "blast"
  } else if (identical(engine, "fm")) {
    "fm"
  } else if (identical(engine, "index")) {
    "index"
  } else if (isTRUE(engine)) {
    "pair"
  } else {
    "none"
  }
  primer_level <- tier %in% c("blast", "fm", "index")
  # a backend that is not a list at all (no genome given) carries no engine
  if (identical(tier, "none")) {
    out <- list(
      verdict = "unknown",
      reason = if (is.null(genome)) {
        gettext("No genome-wide engine was consulted; a missing engine is never a pass.")
      } else if (nzchar(backend$reason)) {
        backend$reason
      } else {
        gettext("No genome-wide engine for this genome; a missing engine is never a pass.")
      },
      tier = tier,
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = n_tol, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else if (n_tol == 0L) {
    # a CRISPR guide whose PAM changed has lost its target: there is no
    # tolerated product to be specific about, so this is not a specificity pass
    reason <- if (!nrow(products)) {
      gettext("No in-silico PCR product at all: the primers do not amplify.")
    } else if (identical(pam_val, "pam_loss")) {
      gettext("The guide's PAM changed: a hard loss of target, not a tolerated mismatch.")
    } else {
      gettext("Every product breaks the scheme's mismatch policy.")
    }
    out <- list(
      verdict = "fail",
      reason = reason,
      tier = tier,
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = 0L, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else if (n_tol > 1L) {
    out <- list(
      verdict = "fail",
      reason = gettextf(
        "%d tolerated in-silico PCR products: not specific.", n_tol
      ),
      tier = tier,
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = n_tol, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else if (primer_level) {
    # the product data is this pair's primer evidence, not an off-target count:
    # promote it to a pass only if a genome-wide primer scan was consulted
    out <- list(
      verdict = "pass",
      reason = gettext(
        "One tolerated product and a genome-wide primer scan: no other binding site."
      ),
      tier = tier,
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = 1L, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  } else {
    out <- list(
      verdict = "pass",
      reason = gettext(paste(
        "Exactly one tolerated in-silico PCR product; pair-level evidence only,",
        "not a genome-wide primer scan."
      )),
      tier = tier,
      scheme = pol$mode, policy = pol, n_products = nrow(products),
      n_tolerated = 1L, n_mismatch = best$n_mismatch %||% NA_integer_,
      best = best, products = products, pam = pam_val, backend = backend
    )
  }
  structure(out, class = c("specificity_verdict", "list"))
}

#' Which caps of a policy a product breaks
#' @keywords internal
#' @noRd
pcr_violations <- function(n_mismatch, seed_mismatch, three_prime_mismatch,
                           three_prime_terminal, policy) {
  v <- character()
  if (n_mismatch > policy$max_mismatch) v <- c(v, "total")
  if (seed_mismatch > policy$seed_max) v <- c(v, "seed")
  if (three_prime_mismatch > policy$three_prime_max) v <- c(v, "three_prime")
  if (isTRUE(three_prime_terminal)) v <- c(v, "terminal")
  v
}

#' @export
print.specificity_verdict <- function(x, ...) {
  cat("<specificity_verdict:", x$verdict, ">\n")
  cat("  scheme:", x$scheme, "| tier:", x$tier, "| products:", x$n_products,
    "| tolerated:", x$n_tolerated, "\n",
    sep = " "
  )
  cat("  reason:", x$reason, "\n")
  invisible(x)
}
