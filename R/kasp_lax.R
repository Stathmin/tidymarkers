# =============================================================================
#  kasp_lax.R  --  lax KASP: one ALT primer per haplotype (linked variants OK)
# =============================================================================
#  Pure KASP separates the alleles at a *single* 3'-terminal base, so a site whose
#  ALT is really a whole haplotype - several linked variants riding together - is
#  taken only one base at a time.  A "lax" KASP relaxes exactly that: the ALT
#  primer is allowed to match a haplotype (multiple SNPs, and linked indels), and
#  the site gets *one ALT primer per ALT haplotype*.  Discrimination is then the
#  combination of mismatches to REF across the primer, not one base.
#
#  Nothing else changes: each per-haplotype primer is assembled by
#  assemble_kasp_assay(), so the full KASP QC still applies - thermodynamics,
#  complementarity, internal repeats, off-targets, the decision tree.  Only the
#  binding is relaxed (the variants are *linked*, so baking them into the ALT
#  primer is a feature, not a contested 3' window).
# =============================================================================

#' Design a lax KASP: one ALT primer per haplotype
#'
#' The haplotype-aware member of the KASP family.  Given a tidy haplotype table
#' ([assemble_kasp_assay()]'s `haplotypes=` shape), it designs **one ALT primer
#' per ALT haplotype**, each matching that haplotype's own linked variants (baked
#' in via [design_as_primers()]'s `linked=`), so a target whose ALT is a multi-SNP
#' (or linked-indel) haplotype still makes a KASP.  The REF primer is the ordinary
#' one.  Every per-haplotype primer runs through the usual [assemble_kasp_assay()]
#' pipeline, so thermodynamics, complementarity, internal repeats and off-targets
#' are checked as for pure KASP.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param haplotypes tidy haplotype table (`haplotype`, `freq`, `chrom`, `pos`,
#'   `ref`, `alt`); the ALT haplotypes define the primers.  `NULL` degenerates to
#'   a single pure-KASP primer (no linked variants).
#' @param ... passed to [assemble_kasp_assay()].
#' @return a `kasp_lax` list: `snp`, `haplotypes`, `alt` (the ALT haplotype
#'   names), `assays` (one assembled assay per ALT haplotype), `designed` (those
#'   that designed) and `verdict`.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_lax_marker(genome, snp, haplotypes)
#' }
#' @family KASP design
#' @export
design_lax_marker <- function(genome, snp, haplotypes = NULL, ...) {
  if (is.null(haplotypes)) {
    # no haplotype table: a single, ordinary ALT primer (pure-KASP behaviour)
    a <- suppressWarnings(assemble_kasp_assay(genome, snp, ...))
    out <- list(
      snp = snp, haplotypes = NULL, alt = "ALT", assays = list(ALT = a),
      designed = if (is_target_designable(a)) "ALT" else character(0),
      verdict = if (is_target_designable(a)) "pass" else "fail"
    )
    return(structure(out, class = "kasp_lax"))
  }
  hv <- haplotype_variants(haplotypes, snp)
  alt_hap <- hv$alt
  if (!length(alt_hap)) {
    cli::cli_inform(gettext("No ALT haplotype carries the target allele; nothing to design."))
    return(NULL)
  }
  assays <- purrr::map(alt_hap, function(hp) {
    lk <- NULL
    if (!is.null(hv$linked)) {
      d <- hv$linked[hv$linked$haplotype == hp, c("target", "linked"), drop = FALSE]
      if (nrow(d)) lk <- d
    }
    suppressWarnings(assemble_kasp_assay(genome, snp, variants = hv$variants, linked = lk, ...))
  })
  names(assays) <- alt_hap
  ok <- vapply(assays, is_target_designable, logical(1))
  out <- list(
    snp = snp, haplotypes = hv$haplotypes, alt = alt_hap, assays = assays,
    designed = alt_hap[ok],
    verdict = if (all(ok)) "pass" else if (any(ok)) "warn" else "fail"
  )
  structure(out, class = "kasp_lax")
}

#' @export
print.kasp_lax <- function(x, ...) {
  cat(
    "<kasp_lax>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt,
    " | ", length(x$alt), " ALT haplotype(s) | verdict ", x$verdict, "\n",
    sep = ""
  )
  for (hp in x$alt) {
    a <- x$assays[[hp]]
    cat("  haplotype ", hp, " -> ",
      if (inherits(a, "kasp_doomed")) paste("doomed:", a$reason) else "ok",
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' The lax-KASP scheme spec (see [kasp_scheme_registry()])
#' @keywords internal
#' @noRd
scheme_def_lax <- function() {
  gates <- c("gate_lax_coverage", "gate_lax_discrimination")
  list(
    scheme = "lax", family = "kasp",
    label = "lax KASP (one primer per haplotype)",
    discrimination = "haplotype block (multiple linked variants), not one 3' base",
    readout = "co-dominant two-colour", zygosity = "co-dominant",
    amplicon = "critical", loci = 1L, verdict = "usable lax KASP assay",
    gates = gates, relax = gates, tiers = c(tier1 = "haplotype coverage", tier2 = "discrimination", tier3 = "Tm"),
    gate_short = c(gate_lax_coverage = "haplo coverage", gate_lax_discrimination = "haplo discrim."),
    gate_limit = c(
      gate_lax_coverage = "one binding ALT primer per ALT haplotype",
      gate_lax_discrimination = "each ALT primer mismatches REF at the haplotype"
    ),
    gate_floor_label = c(
      gate_lax_coverage = "leave a haplotype uncovered (warn)",
      gate_lax_discrimination = "accept a weaker per-haplotype mismatch"
    ),
    gate_hard_label = character(0),
    stage_label = "lax haplotype stage"
  )
}
