# =============================================================================
#  marker_panel.R  --  GWAS-to-panel end to end (roadmap ZY1)
# =============================================================================
#  One wrapper over the bricks that already exist.  Input: a genome, the variants
#  and phased haplotypes of the population, the chemistries a lab can run, and
#  the targets (or a region).  Output: a marker_assay per target - or the cause
#  chain that rejected every chemistry - plus the one-pot multiplex estimate and
#  the decision artefact.
#
#  Nothing is designed or grouped here: marker_plan() + marker_solve() do the
#  design, marker_one_pot() does the grouping, marker_artefact() records the
#  decisions.  A target that no chemistry designs is not dropped: its row carries
#  every attempt's cause, and a panel with no design at all aborts with them.
#
#  The tolerance column is `footprint_dropout`, the frequency-weighted fraction of
#  templates the ALT-primer footprint misses (KASP_LIMITS$footprint_dropout is the
#  acceptance edge).  It is NA when the designer reports no haplotype binding
#  ("not measured"), never 0 as a silent pass.
# =============================================================================

#' Design a marker panel from variants, haplotypes and chemistries
#'
#' The whole pipeline as one call: plan the targets, solve them cheapest-first
#' across the allowed chemistries with population context, group the designed
#' assays into one-pot mixes, and write the decision artefact.  Every target ends
#' up in exactly one row of `targets`: a design (`assay`, `scheme`,
#' `footprint_dropout`) or the cause chain of the chemistry attempts that failed.
#'
#' @param genome the reference sequence: a `Biostrings::DNAStringSet`, or a
#'   `marker_genome` (see [marker_genome()]).
#' @param variants optional variant context, as a tidy table (`chrom`, `pos`, and
#'   either `alt` or `ref_genome`; optional `pop_freq`) or a VCF path read with
#'   [read_variants_vcf()].  It is how many other variants sit under a primer.
#' @param haplotypes optional phased haplotypes, as a tidy table (`haplotype`,
#'   `freq`, `chrom`, `pos`, `ref`, `alt`) or a CSV path.  They make selection
#'   population-aware: a primer is scored on the haplotypes it actually binds, not
#'   on an allele frequency alone (see [marker_context()]).
#' @param chemistries scheme ids to allow, e.g. `c("snp", "hrm")` (see
#'   [marker_schemes()]).  The solver walks them cheapest-first and stops at the
#'   first that designs each target.
#' @param targets optional target table (`chrom`, `pos`, `ref`, `alt`; optional
#'   `label`, `type`, `alt2`, `strand`, `id`).
#' @param region optional single region (`GRanges`, or a data frame with
#'   `chrom`, `start`, `end`); it expands to the variants the context places
#'   inside it ([marker_plan()]).
#' @param dir directory the artefact is written to (see [marker_artefact()]).
#' @param ld optional LD step for the target table: `NULL` (not applied), a
#'   numeric `r2` threshold computed from `haplotypes`, or a linkage table with
#'   `target1`/`target2` (labels).  A linked pair keeps the first target; the
#'   record is written to the artefact (`"not_applied"` when `NULL`).
#' @param costs optional chemistry costs passed to [marker_tech()], so the
#'   cheapest-first order is a seam (see [marker_costs()]).
#' @param created optional artefact timestamp; `NULL` reads
#'   `SOURCE_DATE_EPOCH`, then the current time (see [marker_artefact()]).
#' @param ... passed to the per-scheme designer through [marker_solve()] (e.g.
#'   `flank`, `cores`).
#' @return A `marker_panel` list:
#' * `targets` - one row per input target: `target` (the label), `assay` (the
#'   designed `marker_assay`, or `NULL`), `scheme` (the chemistry that designed
#'   it, or `NA`), `footprint_dropout` (frequency-weighted fraction of haplotypes
#'   the ALT primer misses; `NA` when the designer did not report binding),
#'   `cause` (a list column: `"<scheme>: <cause>"` per failed attempt, empty when
#'   designed) and `ok` (the design's logical verdict);
#' * `one_pot` - the [marker_one_pot()] grouping of the designed assays (both
#'   alleles are read on the KASP channels, which is not modelled, so `channel`
#'   is `NA` and the amplicon gap is the conservative rule);
#' * `artefact` - the [marker_artefact()] decision record;
#' * `solved` - the underlying `marker_solved`;
#' * `experimental` - `TRUE`: the one-pot grouping is a bench proposal, not a
#'   validated mix (reagent compatibility and channel overlap are not modelled).
#'
#' If no chemistry designs any target the call aborts with every target and every
#' attempt's cause.
#' @examples
#' # the argument contract is checked before any designer runs
#' genome <- Biostrings::DNAStringSet(c(Chr1 = "ACGTACGTACGTACGT"))
#' tryCatch(
#'   marker_design_panel(genome, chemistries = "snp"),
#'   error = conditionMessage
#' )
#'
#' \dontrun{
#' panel <- marker_design_panel(
#'   ath_genome(),
#'   haplotypes = "population_haplotypes.csv",
#'   chemistries = c("snp", "hrm"),
#'   targets = marker_example("SNP")[, c("chrom", "pos", "ref", "alt")]
#' )
#' panel
#' }
#' @family plan
#' @export
marker_design_panel <- function(genome, variants = NULL, haplotypes = NULL,
                                chemistries, targets = NULL, region = NULL,
                                ld = NULL, costs = NULL, created = NULL,
                                dir = tempdir(), ...) {
  if (!inherits(genome, "DNAStringSet") && !inherits(genome, "marker_genome")) {
    abort_msg(
      "`genome` must be a DNAStringSet or a marker_genome (see marker_genome())."
    )
  }
  if (
    !is.character(chemistries) || !length(chemistries) ||
      anyNA(chemistries) || !all(nzchar(chemistries))
  ) {
    abort_msg(
      "`chemistries` must be one or more scheme ids, e.g. c(\"snp\", \"hrm\")."
    )
  }
  if (is.null(targets) && is.null(region)) {
    abort_msg(
      "Give `targets` (a chrom/pos/ref/alt table) or a `region` to design for."
    )
  }
  variants <- panel_variants(variants)
  haplotypes <- panel_haplotypes(haplotypes)

  # resolve the target rows before the LD and haplotype seams, so each seam
  # checks the same table the solver will design
  target_tbl <- if (!is.null(targets)) normalise_targets(targets) else NULL
  if (is.null(target_tbl) && !is.null(region)) {
    seed <- marker_plan(genome, region = region)
    if (!is.null(variants) || !is.null(haplotypes)) {
      seed <- seed + marker_context(variants = variants, haplotypes = haplotypes)
    }
    target_tbl <- region_targets(seed)
  }
  if (is.null(target_tbl) || !nrow(target_tbl)) {
    abort_msg("`targets` or `region` must resolve to at least one target.")
  }
  panel_duplicate_guard(target_tbl)
  if (!is.null(haplotypes)) {
    haplotypes <- panel_haplotype_contract(haplotypes, target_tbl, variants)
  }
  ldres <- panel_ld(target_tbl, haplotypes, ld)
  target_tbl <- ldres$targets
  # a region stays a region only when nothing pruned its targets; after an LD
  # prune the survivors are solved explicitly
  use_region <- if (is.null(targets) && is.null(ld)) region else NULL
  use_targets <- if (is.null(use_region)) target_tbl else NULL

  plan <- marker_plan(genome, targets = use_targets, region = use_region) +
    marker_tech(schemes = chemistries, costs = costs)
  if (!is.null(variants) || !is.null(haplotypes)) {
    plan <- plan + marker_context(variants = variants, haplotypes = haplotypes)
  }
  plan$ld <- ldres$record

  solved <- marker_solve(plan, progress = FALSE, ...)
  targets_tbl <- panel_targets(solved)
  if (!any(targets_tbl$ok)) {
    cli::cli_abort(c(
      gettext("No chemistry designed any target; every attempt and its cause:"),
      panel_failure_lines(targets_tbl)
    ))
  }

  structure(
    list(
      targets = targets_tbl,
      one_pot = panel_one_pot(targets_tbl),
      artefact = marker_artefact(solved, dir = dir, created = created),
      solved = solved,
      ld = ldres$record,
      experimental = TRUE
    ),
    class = "marker_panel"
  )
}

#' @export
#' @method print marker_panel
#' @family plan
print.marker_panel <- function(x, ...) {
  t <- x$targets
  cat(glue::glue(
    "<marker_panel> {nrow(t)} target(s), {sum(t$ok)} designed, {sum(!t$ok)} failed; ",
    "{x$one_pot$n_groups} pot(s)\n"
  ))
  print(utils::head(
    tibble::tibble(
      target = t$target, scheme = t$scheme,
      footprint_dropout = t$footprint_dropout, ok = t$ok
    ),
    12L
  ))
  if (any(!t$ok)) {
    cat("  no design:\n")
    for (i in which(!t$ok)) {
      cat(glue::glue("    {t$target[[i]]}: {paste(t$cause[[i]], collapse = ' | ')}\n"))
    }
  }
  if (isTRUE(x$experimental)) {
    cat("  experimental: the one-pot grouping is a proposal, not a validated wet-lab mix\n")
  }
  invisible(x)
}

# ---- internals ---------------------------------------------------------------

#' Read a variants input: VCF path or table, `NULL` for none
#' @keywords internal
#' @noRd
panel_variants <- function(variants) {
  if (is.null(variants)) {
    return(NULL)
  }
  if (is.character(variants) && length(variants) == 1L && !is.na(variants)) {
    return(read_variants_vcf(variants))
  }
  variants
}

#' Read and normalise a haplotypes input: CSV path or table, `NULL` for none
#' @keywords internal
#' @noRd
panel_haplotypes <- function(haplotypes) {
  if (is.null(haplotypes)) {
    return(NULL)
  }
  if (is.character(haplotypes) && length(haplotypes) == 1L && !is.na(haplotypes)) {
    haplotypes <- readr::read_csv(haplotypes, show_col_types = FALSE)
  }
  normalise_haplotypes(haplotypes)
}

#' The panel's one-row-per-target table, with the cause chain per failure
#' @keywords internal
#' @noRd
panel_targets <- function(solved) {
  report <- solved$report
  attempts <- solved$attempts
  tibble::tibble(
    target = as.character(unname(report$label)),
    assay = unname(report$marker),
    scheme = unname(report$scheme),
    footprint_dropout = unname(vapply(report$marker, panel_footprint_dropout, numeric(1))),
    # a designed row succeeded, so it has nothing to explain: the chain exists
    # only where the target has no design
    cause = lapply(seq_len(nrow(report)), function(i) {
      panel_causes(attempts, report$label[[i]], report$ok[[i]])
    }),
    ok = unname(report$ok)
  )
}

#' The `"<scheme>: <cause>"` strings of a target that no chemistry designed
#' @keywords internal
#' @noRd
panel_causes <- function(attempts, label, ok) {
  if (isTRUE(ok)) {
    return(character(0))
  }
  at <- attempts[attempts$label == label, , drop = FALSE]
  paste0(at$scheme, ": ", at$cause)
}

#' Frequency-weighted fraction of haplotypes the ALT primer does not bind
#' @keywords internal
#' @noRd
panel_footprint_dropout <- function(marker) {
  hits <- marker$haplotype_hits
  if (is.null(hits) || !nrow(hits) || !"binds" %in% names(hits)) {
    return(NA_real_)
  }
  total <- sum(hits$freq)
  if (!is.finite(total) || total <= 0) {
    return(NA_real_)
  }
  bound <- !is.na(hits$binds) & hits$binds
  1 - sum(hits$freq[bound]) / total
}

#' One field of a designed assay's first primer row, `NA` when it has none
#' @keywords internal
#' @noRd
panel_assay_field <- function(marker, field) {
  a <- marker$assay
  if (is.null(a) || !nrow(a) || !field %in% names(a)) {
    return(NA_real_)
  }
  as.numeric(a[[field]][1L])
}

#' The one-pot estimate for the assays that designed
#' @keywords internal
#' @noRd
panel_one_pot <- function(targets_tbl) {
  designed <- targets_tbl[targets_tbl$ok, , drop = FALSE]
  marker_one_pot(tibble::tibble(
    label = designed$target,
    scheme = designed$scheme,
    tm = vapply(designed$assay, panel_assay_field, numeric(1), field = "com_tm"),
    amplicon = vapply(designed$assay, panel_assay_field, numeric(1), field = "amplicon"),
    # KASP reads both alleles on FAM/HEX; which channel an assay uses is optics,
    # not design, so it stays unassigned and the conservative gap rule applies
    channel = NA_character_
  ))
}

#' One abort bullet per (target, failed attempt)
#' @keywords internal
#' @noRd
panel_failure_lines <- function(targets_tbl) {
  failed <- targets_tbl[!targets_tbl$ok, , drop = FALSE]
  unlist(Map(
    function(target, causes) {
      if (!length(causes)) {
        return(glue::glue("{target}: no cause was recorded"))
      }
      glue::glue("{target}: {causes}")
    },
    failed$target, failed$cause
  ))
}

#' Reject a target table with two rows for the same locus
#' @keywords internal
#' @noRd
panel_duplicate_guard <- function(targets) {
  key <- paste(targets$chrom, targets$pos, targets$ref, targets$alt, sep = "\r")
  dup <- duplicated(key)
  if (any(dup)) {
    i <- which(dup)[[1L]]
    abort_msg(
      "Duplicate target {targets$label[[i]]}: every target must appear once.",
      class = "tidymarkers_duplicate_target"
    )
  }
  invisible(targets)
}

#' The phased-haplotype contract of the panel
#'
#' Checks the input before it reaches a designer: one allele per
#' `(haplotype, chrom, pos)`, `sum(freq) <= 1` per chromosome, a target-position
#' row whose `ref`/`alt` match the target, and (when a variant table is given)
#' every other row present in that table.  Absence at a target is reference; the
#' number of such calls is attached as `reference_calls`.
#' @keywords internal
#' @noRd
panel_haplotype_contract <- function(haplotypes, targets, variants = NULL) {
  h <- normalise_haplotypes(haplotypes)
  key <- paste(h$haplotype, h$chrom, h$pos, sep = "\r")
  dup <- duplicated(key)
  if (any(dup)) {
    i <- which(dup)[[1L]]
    abort_msg(
      "Haplotype {h$haplotype[[i]]} has two alleles at {h$chrom[[i]]}:{h$pos[[i]]}; the input is not phased.",
      class = "tidymarkers_haplotype_error"
    )
  }
  freq <- unique(h[, c("haplotype", "chrom", "freq")])
  total <- tapply(freq$freq, freq$chrom, sum)
  over <- names(total)[total > 1 + 1e-8]
  if (length(over)) {
    chrom <- over[[1L]]
    first <- freq[freq$chrom == chrom, , drop = FALSE][1L, ]
    pos <- h$pos[h$haplotype == first$haplotype & h$chrom == chrom][[1L]]
    abort_msg(
      "Haplotype frequencies sum to {round(total[[chrom]], 3)} on {chrom} (> 1); first offender {first$haplotype} at {chrom}:{pos}.",
      class = "tidymarkers_haplotype_error"
    )
  }
  tg <- tibble::as_tibble(targets)
  for (i in seq_len(nrow(tg))) {
    t <- tg[i, , drop = FALSE]
    at <- h[h$chrom == t$chrom & h$pos == t$pos, , drop = FALSE]
    if (!nrow(at)) {
      next
    }
    allowed <- c(t$alt, if ("alt2" %in% names(t)) t$alt2[!is.na(t$alt2)] else character(0))
    bad <- at$ref != t$ref | !at$alt %in% allowed
    if (any(bad)) {
      j <- which(bad)[[1L]]
      abort_msg(
        "Haplotype {at$haplotype[[j]]} has {at$ref[[j]]}>{at$alt[[j]]} at {at$chrom[[j]]}:{at$pos[[j]]}, but the target is {t$ref}>{t$alt}.",
        class = "tidymarkers_haplotype_error"
      )
    }
  }
  if (!is.null(variants)) {
    v <- tibble::as_tibble(variants)
    need_cols(v, c("chrom", "pos", "ref", "alt"), arg = "variants")
    target_key <- paste(tg$chrom, tg$pos, sep = "\r")
    elsewhere <- !paste(h$chrom, h$pos, sep = "\r") %in% target_key
    if (any(elsewhere)) {
      vk <- paste(v$chrom, as.integer(v$pos), toupper(v$ref), toupper(v$alt), sep = "\r")
      hk <- paste(
        h$chrom[elsewhere], h$pos[elsewhere],
        h$ref[elsewhere], h$alt[elsewhere],
        sep = "\r"
      )
      bad <- !hk %in% vk
      if (any(bad)) {
        i <- which(elsewhere)[which(bad)[[1L]]]
        abort_msg(
          "Haplotype {h$haplotype[[i]]} carries {h$ref[[i]]}>{h$alt[[i]]} at {h$chrom[[i]]}:{h$pos[[i]]}, which is not in the variant table.",
          class = "tidymarkers_haplotype_error"
        )
      }
    }
  }
  ref_calls <- 0L
  for (i in seq_len(nrow(tg))) {
    t <- tg[i, , drop = FALSE]
    on_chrom <- unique(h$haplotype[h$chrom == t$chrom])
    allowed <- c(t$alt, if ("alt2" %in% names(t)) t$alt2[!is.na(t$alt2)] else character(0))
    carriers <- unique(h$haplotype[
      h$chrom == t$chrom & h$pos == t$pos & h$alt %in% allowed
    ])
    ref_calls <- ref_calls + length(setdiff(on_chrom, carriers))
  }
  attr(h, "reference_calls") <- ref_calls
  h
}

#' The opt-in LD step: collapse targets linked at or above a threshold
#'
#' `ld = NULL` records `"not_applied"`.  A numeric threshold computes pairwise
#' `r2` from the phased haplotypes (carrier indicator per haplotype); a linkage
#' table names the pairs (`target1`, `target2` labels) to collapse.  The first
#' target of a pair is kept, so the order is the input order.
#' @keywords internal
#' @noRd
panel_ld <- function(targets, haplotypes = NULL, ld = NULL) {
  if (is.null(ld)) {
    return(list(targets = targets, record = "not_applied"))
  }
  keep <- rep(TRUE, nrow(targets))
  if (is.data.frame(ld)) {
    need_cols(ld, c("target1", "target2"), arg = "ld")
    for (i in seq_len(nrow(ld))) {
      first <- match(as.character(ld$target1[[i]]), targets$label)
      second <- match(as.character(ld$target2[[i]]), targets$label)
      if (is.na(first) || is.na(second)) {
        abort_msg("The LD table names a target that is not in `targets`.")
      }
      keep[second] <- FALSE
    }
    record <- list(table = ld)
  } else if (is.numeric(ld) && length(ld) == 1L && is.finite(ld) && ld >= 0 && ld <= 1) {
    if (is.null(haplotypes)) {
      abort_msg("A numeric `ld` threshold needs `haplotypes` to compute r2.")
    }
    if (nrow(targets) > 1L) {
      hap <- unique(haplotypes$haplotype)
      dosage <- vapply(seq_len(nrow(targets)), function(i) {
        t <- targets[i, , drop = FALSE]
        allowed <- c(t$alt, if ("alt2" %in% names(t)) t$alt2[!is.na(t$alt2)] else character(0))
        carriers <- unique(haplotypes$haplotype[
          haplotypes$chrom == t$chrom & haplotypes$pos == t$pos &
            haplotypes$alt %in% allowed
        ])
        as.numeric(hap %in% carriers)
      }, numeric(length(hap)))
      dosage <- matrix(dosage, nrow = length(hap))
      r2 <- suppressWarnings(stats::cor(dosage)^2)
      for (j in seq_len(nrow(targets))) {
        linked <- which(keep & seq_len(nrow(targets)) < j)
        if (length(linked) && any(r2[linked, j] >= ld, na.rm = TRUE)) {
          keep[[j]] <- FALSE
        }
      }
    }
    record <- list(threshold = ld, dropped = targets$label[!keep])
  } else {
    abort_msg("`ld` must be NULL, one `r2` threshold in [0, 1], or a linkage table.")
  }
  list(targets = targets[keep, , drop = FALSE], record = record)
}
