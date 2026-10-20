# =============================================================================
#  marker_duplex.R  --  two loci in one tube (the duplex chemistry engine)
# =============================================================================
#  The scheme registry described the duplex chemistries (dKASP, mxKASP,
#  Multi-KASP) declaratively; this is their *engine*.  A duplex assay is the two
#  per-locus KASPs plus a compatibility stage across every primer of both loci:
#  harmonised Tm, no inter-pair complementarity, a bounded self-complementarity,
#  and the amplicon budget.  The readout is either co-dominant (dKASP: two
#  cassettes, one per locus) or presence/absence (mxKASP / Multi-KASP: a dominant
#  call from the combination), which is why `homologous` matters.
#
#  Design rules (researched, not invented).  The compatibility screen is the one
#  published for multiplex KASP - Wang Y, Yang J, Gao Y, Ma H, Dai Y, Ma H (2025)
#  "Development and application of a cost-effective multiplex Kompetitive
#  Allele-Specific polymerase chain reaction assay for pyramiding resistant genes
#  of fusarium head blight and powdery mildew in wheat", BMC Plant Biol 25:963,
#  doi:10.1186/s12870-025-07005-8, whose primer screen is stated verbatim:
#    (1) "no individual primer should contain more than five consecutive
#        self-complementary bases, and the 3' end of each primer should not have
#        four or more consecutive bases complementary to any internal region of
#        the same primer";
#    (2) "any two primers within the multiplex panel should avoid having six or
#        more consecutive complementary bases" (inter-primer dimer);
#    (3) "the Tm of all primers should be between 59 C and 65 C";
#    (4) "the expected amplicon length should not exceed 150 bp".
#  The two topologies are the published ones:
#   * dKASP - Jiang P et al. (2023) "Cost-effective duplex Kompetitive Allele
#     Specific PCR markers for homologous genes facilitating wheat breeding",
#     BMC Plant Biol 23:119, doi:10.1186/s12870-023-04116-y: two *homologous*
#     genes share one common primer, one allele-specific primer per gene, each
#     carrying its own cassette; the call is dominant and the 2x2 combination
#     encodes both loci.
#   * mxKASP - Meilink WRM et al. (2025) "Determining Zygosity With Multiplex
#     Kompetitive Allele-Specific PCR (mxKASP) Genotyping", Ecol Evol 15:e71642,
#     doi:10.1002/ece3.71642: two *non-homologous* markers, each with its own
#     primer pair, one dye per locus; presence/absence per locus gives zygosity.
#  Only (1)'s first half and (2)-(4) are sequence criteria that can be *computed*
#  from a candidate; (3) and (4) are enforced here, (4) more strictly through the
#  per-locus KASP amplicon window (see `DUPLEX_MAX_AMPLICON` below).
#
#  The two orders the stage declares (C6): gates in application order (grave /
#  cheap first), then the relaxation ladder (least grave first).
# =============================================================================

# ---- the duplex selection stage ---------------------------------------------
# Order = grave / cheap first: the declared topology (a wrong readout cannot be
# rescued), then inter-primer cross-talk (fatal for a multiplex), then the shared
# Tm window, then the mild self-complementarity cap.  Only the last is relaxable.
KASP_GATES_DUPLEX <- c(
  "gate_duplex_targets", "gate_duplex_complement", "gate_duplex_tm",
  "gate_duplex_selfcomp"
)
# The ladder relaxes the *least grave* criterion: one extra self-complementary
# base (never the cross-talk rules - Wang 2025 rule 2 is fatal for a multiplex).
KASP_RELAX_DUPLEX <- "gate_duplex_selfcomp"
# Tiers rank the survivors; both are margins *of this stage's own soft gate* and
# of an ungated budget, so neither merely restates a hard gate (C14).
KASP_TIERS_DUPLEX <- c(
  tier1 = "self-complementarity headroom", tier2 = "amplicon budget headroom"
)

# The published limits, named once so the labels and the engine cannot drift.
DUPLEX_TM_WINDOW <- c(59, 65) # C     (Wang 2025 rule 3)
DUPLEX_MAX_SELFCOMP <- 5L # bases (Wang 2025 rule 1, first half)
DUPLEX_MIN_CROSS_RUN <- 6L # bp    (Wang 2025 rule 2)
DUPLEX_MAX_AMPLICON <- 150L # bp   (Wang 2025 rule 4; the per-locus KASP window
#                                   [50, 100] is stricter, so this is a budget,
#                                   not a second gate - see KASP_LIMITS$amplicon)
DUPLEX_SELFCOMP_FLOOR <- -2L # one rung: up to 7 self-complementary bases

# Numeric floors per soft gate, in margin units (0 = hard).  Mirrors the
# KASP_GATE_FLOOR_AS / _PAIR convention in R/marker_select.R; the plan flagged that
# the duplex stage declared a ladder with no floor behind it.
KASP_GATE_FLOOR_DUPLEX <- c(gate_duplex_selfcomp = DUPLEX_SELFCOMP_FLOOR)

# This scheme's entries in the shared gate-label lookups (assembled in
# R/marker_select.R), kept beside the gates they name.
KASP_GATE_LABEL_DUPLEX <- c(
  gate_duplex_targets = "the readout does not match the loci's relationship",
  gate_duplex_complement = "inter-pair primer complementarity",
  gate_duplex_tm = "primers outside the shared Tm window",
  gate_duplex_selfcomp = "primer self-complementarity above the duplex cap"
)
KASP_GATE_SHORT_DUPLEX <- c(
  gate_duplex_targets = "duplex targets",
  gate_duplex_complement = "inter-pair complementarity",
  gate_duplex_tm = "Tm window (all primers)",
  gate_duplex_selfcomp = "self-complementarity"
)
KASP_GATE_LIMIT_DUPLEX <- c(
  gate_duplex_targets = "homologous loci (dKASP) or non-homologous (mxKASP)",
  gate_duplex_complement = glue::glue("< {DUPLEX_MIN_CROSS_RUN} bp inter-pair complementarity"),
  gate_duplex_tm = glue::glue("all primers Tm in [{DUPLEX_TM_WINDOW[1L]}, {DUPLEX_TM_WINDOW[2L]}] C"),
  gate_duplex_selfcomp = glue::glue(
    "<= {DUPLEX_MAX_SELFCOMP} self-complementary bases, no >= 4 bp 3' fold-back"
  )
)
KASP_GATE_FLOOR_LABEL_DUPLEX <- c(
  gate_duplex_selfcomp = glue::glue("up to {-DUPLEX_SELFCOMP_FLOOR} extra self-complementary bases")
)
KASP_GATE_HARD_LABEL_DUPLEX <- c(
  gate_duplex_targets = "the readout must match the loci's relationship",
  gate_duplex_complement = "primer cross-talk defeats multiplexing",
  gate_duplex_tm = "all primers must share one annealing window"
)

#' Longest complementary run between two oligos
#'
#' The longest stretch over which `a` matches the reverse complement of `b` (i.e.
#' the longest run of Watson-Crick complementarity), via the classic longest
#' common substring dynamic program.  This is the base-pair-count rule multiplex
#' design uses ("no >= 6 bp inter-primer complementarity"; Wang 2025 rule 2); the
#' thermodynamic dimer dG ([dimer_thermo()]) is the orthogonal, stricter measure.
#' @param a,b oligos, 5'->3'.
#' @return a single integer (base pairs).
#' @keywords internal
#' @noRd
complementary_run <- function(a, b) {
  x <- strsplit(as.character(a)[1L], "")[[1L]]
  y <- strsplit(
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(as.character(b)[1L]))), ""
  )[[1L]]
  prev <- integer(length(y) + 1L)
  best <- 0L
  for (i in seq_along(x)) {
    cur <- integer(length(y) + 1L)
    for (j in seq_along(y)) {
      if (identical(x[i], y[j])) {
        cur[j + 1L] <- prev[j] + 1L
        if (cur[j + 1L] > best) best <- cur[j + 1L]
      }
    }
    prev <- cur
  }
  best
}

#' A primer's self-complementarity, both halves of the multiplex rule
#'
#' Wang 2025's first screen has two parts that are not the same test: the
#' *longest* self-complementary run anywhere in the oligo must be at most five
#' bases, **and** the 3' end must not have four or more consecutive bases
#' complementary to an internal region (a 3'-anchored fold-back).  A 4-base 3'
#' fold passes the first test (4 <= 5) but fails the second, so both are reported.
#' @param primer an oligo, 5'->3'.
#' @return list with `run` (longest self-complementary run), `fold` (the 3' end
#'   folds back >= 4 bp) and `observed` (the run the cap is compared against).
#' @keywords internal
#' @noRd
duplex_selfcomp <- function(primer) {
  p <- as.character(primer)[1L]
  n <- nchar(p)
  run <- complementary_run(p, p)
  fold <- FALSE
  if (n > 4L) {
    tail4 <- substr(p, n - 3L, n)
    tail_rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(tail4)))
    # an *internal* region excludes the 3'-terminal 4 bases themselves
    fold <- grepl(tail_rc, substr(p, 1L, n - 4L), fixed = TRUE)
  }
  list(run = run, fold = fold, observed = run)
}

#' Raw compatibility metrics of a duplex pairing
#'
#' One evaluation of the two-locus stage, shared by the gate table
#' (`duplex_compatibility()`) and the relaxation ladder (`duplex_selection()`) so
#' the two can never disagree.
#' @param assays a list of two assembled assays (one per locus).
#' @param readout `"co-dominant"` (one cassette per locus) or `"presence-absence"`.
#' @param homologous are the two loci homologous? (see `duplex_compatibility()`).
#' @param tm_window,max_amplicon,max_selfcomp,min_cross_run the duplex limits.
#' @return a list of scalars/vectors (see the body).
#' @keywords internal
#' @noRd
duplex_metrics <- function(assays, readout, homologous, tm_window, max_amplicon,
                           max_selfcomp, min_cross_run) {
  primers <- lapply(assays, function(a) {
    c(a$assay$as_ref[1L], a$assay$as_alt[1L], a$assay$com[1L])
  })
  all_pr <- unlist(primers, use.names = FALSE)
  tm <- primer_tm(all_pr)
  # the declared topology and the loci's relationship must agree: a shared-cassette
  # duplex (co-dominant, one cassette per locus) is the dKASP build and pairs
  # *homologous* genes, so one common primer serves both (Jiang 2023); the
  # independent-pair builds (mxKASP / Multi-KASP) target *non-homologous* loci
  # (Meilink 2025, Wang 2025).  An undeclared relationship is noted, not assumed.
  tgt_pass <- !(readout == "co-dominant" && isFALSE(homologous))
  # cross-complementarity only *between* the loci's primers (within a locus the
  # per-locus dimer screen already covers it)
  cross <- unlist(lapply(seq_along(primers), function(i) {
    lapply(seq_along(primers), function(j) {
      if (i >= j) {
        return(0L)
      }
      outer(primers[[i]], primers[[j]], Vectorize(complementary_run))
    })
  }))
  self <- lapply(all_pr, duplex_selfcomp)
  amp <- vapply(assays, function(a) as.numeric(a$assay$amplicon[1L]), numeric(1))
  list(
    primers = all_pr, tm = tm, tgt_pass = tgt_pass, cross = max(cross),
    self_run = max(vapply(self, `[[`, integer(1), "run")),
    self_fold = any(vapply(self, `[[`, logical(1), "fold")),
    amplicon = max(amp),
    tm_window = tm_window, max_amplicon = max_amplicon,
    max_selfcomp = max_selfcomp, min_cross_run = min_cross_run
  )
}

#' The compatibility gates of a duplex pairing
#'
#' Evaluates the two-locus stage of [marker_schemes()] `"duplex"` on the primers of
#' two assembled per-locus assays, using the same gate ids the decision registry
#' declares (`KASP_GATES_DUPLEX`), so the selector, the tree and this engine share
#' one vocabulary.  Every gate here is computed from the candidate primers: the
#' declared topology, the inter-pair complementarity, the all-primer Tm window and
#' the self-complementarity cap.  The multiplex amplicon budget (Wang 2025 rule 4)
#' is *not* a gate - the per-locus KASP amplicon window is already stricter - it
#' is reported as the `amplicon budget headroom` tier.
#' @param assays a list of two assembled assays (one per locus).
#' @param readout `"co-dominant"` (dKASP: a cassette per locus) or
#'   `"presence-absence"` (mxKASP / Multi-KASP: a dominant call).
#' @param homologous are the two loci homologous?  `TRUE` for dKASP, `FALSE` for
#'   mxKASP/Multi-KASP, `NA` when undeclared (the gate then notes it).
#' @param tm_window,max_amplicon,max_selfcomp,min_cross_run the duplex limits (see
#'   the published rules in the file header).
#' @return a tidy tibble `gate`, `observed`, `limit`, `pass`, `detail`.
#' @keywords internal
#' @noRd
duplex_compatibility <- function(assays, readout = c("co-dominant", "presence-absence"),
                                 homologous = NA, tm_window = DUPLEX_TM_WINDOW,
                                 max_amplicon = DUPLEX_MAX_AMPLICON,
                                 max_selfcomp = DUPLEX_MAX_SELFCOMP,
                                 min_cross_run = DUPLEX_MIN_CROSS_RUN) {
  readout <- match.arg(readout)
  m <- duplex_metrics(
    assays, readout, homologous, tm_window, max_amplicon, max_selfcomp, min_cross_run
  )
  gate <- c(
    "gate_duplex_targets", "gate_duplex_complement", "gate_duplex_tm",
    "gate_duplex_selfcomp"
  )
  obs <- c(
    dplyr::if_else(is.na(homologous), NA_real_, as.numeric(homologous)),
    m$cross, max(m$tm) - min(m$tm), m$self_run
  )
  lim <- c(
    if (readout == "co-dominant") "homologous" else "any",
    glue::glue("< {min_cross_run} bp"),
    glue::glue("<= {tm_window[2L] - tm_window[1L]} C spread, all in [{tm_window[1L]}, {tm_window[2L]}] C"),
    glue::glue("<= {max_selfcomp} bp, no >= 4 bp 3' fold")
  )
  pass <- c(
    m$tgt_pass,
    m$cross < min_cross_run,
    all(dplyr::between(m$tm, tm_window[1L], tm_window[2L])),
    m$self_run <= max_selfcomp && !m$self_fold
  )
  detail <- c(
    if (is.na(homologous)) "homology undeclared" else if (homologous) "homologous" else "non-homologous",
    NA, NA,
    if (m$self_fold) "3' end folds back (>= 4 bp)" else NA
  )
  tibble::tibble(gate = gate, observed = obs, limit = lim, pass = pass, detail = detail)
}

#' Pick the duplex pairing through the declared relaxation ladder
#'
#' Turns the compatibility metrics into the *margins* the shared selector ranks
#' (`select_with_ladder()`): a hard gate (floor 0) dooms, a soft one
#' (`gate_duplex_selfcomp`) may be lowered as far as its published floor, and the
#' winner is the (single) pairing the caller assembled.  This is what makes the
#' stage's declared ladder real (C7) and a hard failure DOOMED rather than a mere
#' warning (C4).
#' @param assays a list of two assembled assays (one per locus).
#' @param readout,homologous,tm_window,max_amplicon,max_selfcomp,min_cross_run see
#'   `duplex_compatibility()`.
#' @return list with `pick` (the `select_with_ladder()` result) and `metrics`.
#' @keywords internal
#' @noRd
duplex_selection <- function(assays, readout = c("co-dominant", "presence-absence"),
                             homologous = NA, tm_window = DUPLEX_TM_WINDOW,
                             max_amplicon = DUPLEX_MAX_AMPLICON,
                             max_selfcomp = DUPLEX_MAX_SELFCOMP,
                             min_cross_run = DUPLEX_MIN_CROSS_RUN) {
  readout <- match.arg(readout)
  m <- duplex_metrics(
    assays, readout, homologous, tm_window, max_amplicon, max_selfcomp, min_cross_run
  )
  # one row per candidate pairing: the two best per-locus designs are fixed, so
  # the duplex stage has exactly one candidate and the ladder bounds how far its
  # soft cap may be lowered before the pairing is doomed
  data <- tibble::tibble(
    m_targets = if (m$tgt_pass) 0 else -1,
    m_complement = min_cross_run - m$cross,
    m_tm = min(min(m$tm) - tm_window[1L], tm_window[2L] - max(m$tm)),
    m_selfcomp = max_selfcomp - m$self_run,
    tier1 = max_selfcomp - m$self_run,
    tier2 = max_amplicon - m$amplicon
  )
  pick <- select_with_ladder(
    data,
    gates = KASP_GATES_DUPLEX,
    relax_order = KASP_RELAX_DUPLEX,
    tiers = names(KASP_TIERS_DUPLEX),
    ties = character(0),
    gate_margin = c(
      gate_duplex_targets = "m_targets", gate_duplex_complement = "m_complement",
      gate_duplex_tm = "m_tm", gate_duplex_selfcomp = "m_selfcomp"
    ),
    gate_floor = KASP_GATE_FLOOR_DUPLEX
  )
  list(pick = pick, metrics = m)
}

#' Design a duplex assay: two loci in one tube
#'
#' Assembles each locus on its own ([assemble_kasp_assay()]) and then evaluates
#' the two-locus **compatibility stage** (`duplex_compatibility()`) across every
#' primer of both loci - the declared topology, harmonised Tm, no inter-pair
#' complementarity and a bounded self-complementarity - through the declared
#' relaxation ladder (`duplex_selection()`).  A failed *hard* gate dooms the
#' pairing (`verdict = "fail"`); a relaxed soft gate warns.  The readout is either
#' co-dominant (`"dKASP"`: a cassette per locus, `homologous = TRUE`) or
#' presence/absence (`"mxKASP"`/`"Multi-KASP"`: a dominant call from the
#' combination, `homologous = FALSE`).
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param loci two-row tibble (`chrom`, `pos`, `ref`, `alt`).
#' @param readout `"co-dominant"` or `"presence-absence"`.
#' @param homologous are the two loci homologous? (see `duplex_compatibility()`).
#' @param ... passed to [assemble_kasp_assay()].
#' @return a `marker_duplex` list (`loci`, `assays`, `readout`, `homologous`,
#'   `compatibility`, `selection`, `verdict`).
#' @examples
#' \dontrun{
#' loci <- marker_example()[1:2, c("chrom", "pos", "ref", "alt")]
#' design_duplex(ath_genome(), loci, readout = "presence-absence", homologous = FALSE)
#' }
#' @family KASP design
#' @export
design_duplex <- function(genome, loci, readout = c("co-dominant", "presence-absence"),
                          homologous = NA, ...) {
  readout <- match.arg(readout)
  if (nrow(loci) != 2L) {
    cli::cli_abort(gettext("`loci` must have exactly two rows for a duplex assay."))
  }
  assays <- purrr::map(seq_len(2L), function(i) {
    suppressWarnings(assemble_kasp_assay(
      genome, loci[i, , drop = FALSE], ...
    ))
  })
  # a duplex needs *both* loci: one doomed locus dooms the pairing (there is no
  # half-duplex), so this is a fail, never a warning
  ok <- !vapply(assays, inherits, logical(1), "marker_doomed")
  compat <- NULL
  sel <- NULL
  verdict <- "fail"
  if (all(ok)) {
    compat <- duplex_compatibility(assays, readout, homologous)
    sel <- duplex_selection(assays, readout, homologous)
    # a hard gate has floor 0 and is never lowered, so `best = NULL` means a hard
    # failure (DOOMED); a non-empty `violated` means the soft cap was relaxed
    verdict <- if (is.null(sel$pick$best)) {
      "fail"
    } else if (length(sel$pick$violated)) {
      "warn"
    } else {
      "pass"
    }
    if (identical(verdict, "warn")) {
      cli::cli_warn(c(
        gettext("Duplex pairing accepted only after relaxing a soft compatibility cap."),
        rlang::set_names(
          gettextf("relaxed: %s", paste(sel$pick$violated, collapse = ", ")), "*"
        )
      ))
    }
  }
  out <- list(
    loci = loci, assays = assays, readout = readout, homologous = homologous,
    compatibility = compat, selection = sel, verdict = verdict
  )
  structure(out, class = "marker_duplex")
}

#' @export
print.marker_duplex <- function(x, ...) {
  cat("<marker_duplex>", x$readout, "| verdict", x$verdict, "\n")
  for (i in seq_along(x$assays)) {
    a <- x$assays[[i]]
    cat("  locus", i, ":", x$loci$chrom[i], ":", x$loci$pos[i], " -> ",
      if (inherits(a, "marker_doomed")) paste("doomed:", a$reason) else "ok",
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$compatibility)) {
    bad <- x$compatibility$gate[!x$compatibility$pass]
    cat("  compatibility:", if (length(bad)) paste("failed", paste(bad, collapse = ", ")) else "ok", "\n")
  }
  if (identical(x$verdict, "warn") && length(x$selection$pick$violated)) {
    cat("  relaxed:", paste(x$selection$pick$violated, collapse = ", "), "\n")
  }
  invisible(x)
}
