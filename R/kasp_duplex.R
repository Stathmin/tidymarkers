# =============================================================================
#  kasp_duplex.R  --  two loci in one tube (the duplex chemistry engine)
# =============================================================================
#  The scheme registry described the duplex chemistries (dKASP, mxKASP,
#  Multi-KASP) declaratively; this is their *engine*.  A duplex assay is the two
#  per-locus KASPs plus a compatibility stage across every primer of both loci:
#  harmonised Tm, no inter-pair complementarity (cross-dimers grow quadratically
#  with primer number and are the usual multiplex failure - multiplex-PCR design
#  guidance), a bounded self-complementarity, and an amplicon budget.  The readout
#  is either co-dominant (dKASP: two cassettes, one per locus) or presence/absence
#  (Multi-KASP: a dominant call from the combination), which is why `homologous`
#  is only meaningful for the co-dominant reading.
# =============================================================================

#' Longest complementary run between two oligos
#'
#' The longest stretch over which `a` matches the reverse complement of `b` (i.e.
#' the longest run of Watson-Crick complementarity), via the classic longest
#' common substring dynamic program.  This is the base-pair-count rule multiplex
#' design uses ("no >= 6 bp inter-primer complementarity"); the thermodynamic
#' dimer dG ([dimer_thermo()]) is the orthogonal, stricter measure.
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

#' The compatibility gates of a duplex pairing
#'
#' Evaluates the two-locus stage of [kasp_schemes()] `"duplex"` on the primers of
#' two assembled per-locus assays, using the same gate ids the decision registry
#' declares (`KASP_GATES_DUPLEX`), so the selector, the tree and this engine share
#' one vocabulary.
#' @param assays a list of two assembled assays (one per locus).
#' @param readout `"co-dominant"` (dKASP: a cassette per locus) or
#'   `"presence-absence"` (Multi-KASP: a dominant call from the combination).
#' @param homologous are the two loci homologous?  `TRUE` for dKASP, `FALSE` for
#'   mxKASP/Multi-KASP, `NA` when undeclared (the gate then notes it).
#' @param tm_window,max_amplicon,max_selfcomp,min_cross_run the duplex limits (see
#'   the published rules in `kasp_schemes()`).
#' @return a tidy tibble `gate`, `observed`, `limit`, `pass`, `detail`.
#' @keywords internal
#' @noRd
duplex_compatibility <- function(assays, readout = c("co-dominant", "presence-absence"),
                                 homologous = NA, tm_window = c(59, 65),
                                 max_amplicon = 150L, max_selfcomp = 5L, min_cross_run = 6L) {
  readout <- match.arg(readout)
  primers <- lapply(assays, function(a) {
    c(a$assay$as_ref[1L], a$assay$as_alt[1L], a$assay$com[1L])
  })
  all_pr <- unlist(primers, use.names = FALSE)
  tm <- primer_tm(all_pr)
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
  self <- vapply(all_pr, function(p) complementary_run(p, p), integer(1))
  amp <- vapply(assays, function(a) as.numeric(a$assay$amplicon[1L]), numeric(1))
  gate <- c(
    "gate_duplex_targets", "gate_duplex_tm", "gate_duplex_complement",
    "gate_duplex_selfcomp", "gate_duplex_amplicon"
  )
  obs <- c(
    ifelse(is.na(homologous), NA, as.numeric(homologous)), max(tm) - min(tm),
    max(cross), max(self), max(amp)
  )
  lim <- c(
    if (readout == "co-dominant") "homologous" else "any",
    glue::glue("[{tm_window[1]}, {tm_window[2]}] C"),
    glue::glue("< {min_cross_run} bp"), glue::glue("<= {max_selfcomp} bp"),
    glue::glue("<= {max_amplicon} bp")
  )
  pass <- c(
    tgt_pass,
    all(dplyr::between(tm, tm_window[1L], tm_window[2L])),
    max(cross) < min_cross_run, max(self) <= max_selfcomp, max(amp) <= max_amplicon
  )
  detail <- c(
    if (is.na(homologous)) "homology undeclared" else if (homologous) "homologous" else "non-homologous",
    NA, NA, NA, NA
  )
  tibble::tibble(gate = gate, observed = obs, limit = lim, pass = pass, detail = detail)
}

#' Design a duplex assay: two loci in one tube
#'
#' Assembles each locus on its own ([assemble_kasp_assay()]) and then evaluates
#' the two-locus **compatibility stage** ([duplex_compatibility()]) across every
#' primer of both loci - harmonised Tm, no inter-pair complementarity, a bounded
#' self-complementarity and the amplicon budget.  The readout is either
#' co-dominant (`"dKASP"`: a cassette per locus, `homologous = TRUE`) or
#' presence/absence (`"Multi-KASP"`: a dominant call from the combination,
#' `homologous = FALSE`); a failed locus or gate makes the verdict `"fail"`.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param loci two-row tibble (`chrom`, `pos`, `ref`, `alt`).
#' @param readout `"co-dominant"` or `"presence-absence"`.
#' @param homologous are the two loci homologous? (see [duplex_compatibility()]).
#' @param ... passed to [assemble_kasp_assay()].
#' @return a `kasp_duplex` list (`loci`, `assays`, `readout`, `homologous`,
#'   `compatibility`, `verdict`).
#' @examples
#' \dontrun{
#' loci <- kasp_example()[1:2, c("chrom", "pos", "ref", "alt")]
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
  ok <- !vapply(assays, inherits, logical(1), "kasp_doomed")
  compat <- if (all(ok)) duplex_compatibility(assays, readout, homologous) else NULL
  verdict <- if (all(ok) && all(compat$pass)) "pass" else if (any(ok)) "warn" else "fail"
  out <- list(
    loci = loci, assays = assays, readout = readout, homologous = homologous,
    compatibility = compat, verdict = verdict
  )
  structure(out, class = "kasp_duplex")
}

#' @export
print.kasp_duplex <- function(x, ...) {
  cat("<kasp_duplex>", x$readout, "| verdict", x$verdict, "\n")
  for (i in seq_along(x$assays)) {
    a <- x$assays[[i]]
    cat("  locus", i, ":", x$loci$chrom[i], ":", x$loci$pos[i], " -> ",
      if (inherits(a, "kasp_doomed")) paste("doomed:", a$reason) else "ok",
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$compatibility)) {
    bad <- x$compatibility$gate[!x$compatibility$pass]
    cat("  compatibility:", if (length(bad)) paste("failed", paste(bad, collapse = ", ")) else "ok", "\n")
  }
  invisible(x)
}
