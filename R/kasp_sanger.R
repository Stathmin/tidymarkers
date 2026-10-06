# =============================================================================
#  kasp_sanger.R  --  Sanger sequencing (dideoxy read-out) as a marker scheme
# =============================================================================
#  Sanger reads the variant as a *base in a chromatogram*, not as a ratio or a
#  band, so there is no allele-specific primer and no restriction site: one
#  flanking primer pair delimits an amplicon that is read with a sequencing
#  primer.  The design rules below are the consensus of the sequencing-core
#  guidance surveyed for this file:
#
#   * amplicon ~400-800 bp around the variant, with the variant >= ~50 bp from
#     the primer binding sites so it sits inside the high-quality read and never
#     under a primer;
#     https://duckduckgo.com/html/?q=Sanger+sequencing+amplicon+design+primer+400-800+bp+variant+distance+from+primer
#     (Sanger validation guide: amplicon 400-800 bp, variant >= 50 bp from the
#     primer binding sites, specificity checks against the reference genome).
#   * primer length 17-25 nt, Tm ~55-65 C (50-70 acceptable), GC ~50% (we keep
#     the package GC window 40-60%), no hairpins / dimers, avoid long runs and
#     >2 identical 3' bases (Cornell BRC Sanger primer-design guide);
#   * the first ~50 bases of a read are low quality (a 650 bp read with 50 bp
#     low-quality edges leaves a 550 bp usable window), which is the ">= 50 bp
#     from the read start" rule above;
#     https://molbiologytools.com/sanger_primer_designer_help
#
#  The read-out is not simulated (no chromatogram model), so this returns a
#  well-formed, unique amplicon with the variant correctly placed, not a trace.
# =============================================================================

# The Sanger acceptance windows (the constraints the gates below encode).  Tm is
# the sequencing-core 55-65 C rather than the tighter KASP COM window; min_dist
# is the low-quality lead of a read.
KASP_SANGER_LIMITS <- list(
  amplicon = c(400L, 800L),
  tm = c(55, 65),
  gc = c(0.40, 0.60),
  min_dist = 50L
)

# The Sanger selection gates, in pipeline order.  `dist` and `unique` are hard:
# a variant under a primer or a non-unique amplicon cannot be read at all.
KASP_GATES_SANGER <- c("gate_sanger_amplicon", "gate_sanger_tm", "gate_sanger_dist", "gate_sanger_unique")
KASP_RELAX_SANGER <- c("gate_sanger_amplicon", "gate_sanger_tm")
KASP_TIERS_SANGER <- c(tier1 = "unique amplicon", tier2 = "variant centrality", tier3 = "Tm balance")

# Short names for the decision graph.
KASP_GATE_SHORT_SANGER <- c(
  gate_sanger_amplicon = "amplicon length (Sanger)",
  gate_sanger_tm = "primer Tm (Sanger)",
  gate_sanger_dist = "variant from read start",
  gate_sanger_unique = "unique amplicon"
)

# Human-readable relaxation floor and hard-stop reason per gate.
KASP_GATE_FLOOR_LABEL_SANGER <- c(
  gate_sanger_amplicon = "widen the amplicon window",
  gate_sanger_tm = "widen the Tm window"
)
KASP_GATE_HARD_LABEL_SANGER <- c(
  gate_sanger_dist = "a variant under a primer or in the low-quality lead is unreadable",
  gate_sanger_unique = "a non-unique amplicon cannot be called"
)

#' Design a Sanger sequencing marker at a variant
#'
#' The amplicon analogue of [assemble_kasp_assay()]: a single flanking pair
#' delimits a ~400-800 bp product, the variant sits at least `min_dist` bp
#' inside both ends (clear of the low-quality read lead and of either primer),
#' and the pair is checked for genome-wide uniqueness with [in_silico_pcr()].
#' There is no allele-specific primer - Sanger reads the variant directly - so
#' the design optimises the amplicon window, the variant placement and the Tm
#' balance, then gates on off-target uniqueness.
#'
#' It does **not** simulate the chromatogram, so it returns a well-formed,
#' unique amplicon, not a predicted trace.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width (must leave room for the amplicon).
#' @param min_dist minimal distance from the variant to each amplicon end (the
#'   read start / primer binding site).
#' @param min_product,max_product the accepted amplicon window.
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows.
#' @param check_unique when `TRUE`, require exactly one [in_silico_pcr()] product.
#' @param ... passed to [in_silico_pcr()] (e.g. `BPPARAM`).
#' @return a `kasp_sanger` list (`snp`, `primers`, `amplicon`), or `NULL` (with a
#'   message) when no unique pair passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_sanger_marker(genome, kasp_example("SNP")[1, ], flank = 500)
#' }
#' @family KASP Sanger
#' @export
design_sanger_marker <- function(genome, snp, flank = 500L,
                                 min_dist = KASP_SANGER_LIMITS$min_dist,
                                 min_product = KASP_SANGER_LIMITS$amplicon[1L],
                                 max_product = KASP_SANGER_LIMITS$amplicon[2L],
                                 lengths = 18:25, tm_range = KASP_SANGER_LIMITS$tm,
                                 gc_range = KASP_SANGER_LIMITS$gc,
                                 check_unique = TRUE, ...) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  si <- tmpl$snp_i
  # A guard band around the variant: design_flanking_pair() then keeps every
  # forward primer entirely left of `lo` and every reverse primer right of `hi`,
  # which is exactly the "variant >= min_dist from each read start" rule.
  lo <- si - as.integer(min_dist) + 1L
  hi <- si + as.integer(min_dist) - 1L
  primers <- design_flanking_pair(tmpl$seq, lo, hi, lengths, min_product, max_product, tm_range, gc_range)
  if (is.null(primers)) {
    cli::cli_inform(gettext("No flanking pair kept the variant clear of both read starts for Sanger."))
    return(NULL)
  }
  if (check_unique) {
    hits <- in_silico_pcr(primers$forward, primers$reverse, genome, ...)
    if (nrow(hits) != 1L) {
      cli::cli_inform(gettext("The Sanger pair is not a unique amplicon ({nrow(hits)} in-silico products)."))
      return(NULL)
    }
  }
  out <- list(
    snp = snp, primers = primers,
    amplicon = tibble::tibble(
      start = snp$pos + primers$start - si,
      end = snp$pos + primers$end - si
    )
  )
  structure(out, class = "kasp_sanger")
}

#' @export
print.kasp_sanger <- function(x, ...) {
  cat(
    "<kasp_sanger>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  amplicon:", x$primers$product, "bp (Sanger read)\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  invisible(x)
}

#' The Sanger scheme definition (internal, wired into the shared registry)
#'
#' One scheme's declarative row plus its decision-tree vocabulary, in the shape
#' the registry expects: the chemistry's attributes, its selection gates (hard /
#' soft), its ranking tiers, and the short / limit / floor / hard labels the
#' decision graph and text export read.
#' @return a list (`scheme`, `family`, `label`, `discrimination`, `readout`,
#'   `zygosity`, `amplicon`, `loci`, `verdict`, `gates`, `relax`, `tiers`,
#'   `gate_short`, `gate_limit`, `gate_floor_label`, `gate_hard_label`,
#'   `stage_label`).
#' @keywords internal
#' @noRd
scheme_def_sanger <- function() {
  amp <- KASP_SANGER_LIMITS$amplicon
  tm <- KASP_SANGER_LIMITS$tm
  dist <- KASP_SANGER_LIMITS$min_dist
  list(
    scheme = "sanger",
    family = "sanger",
    label = "Sanger sequencing (dideoxy read-out)",
    discrimination = "direct base call of the variant",
    readout = "dideoxy chromatogram (capillary)",
    zygosity = "co-dominant (direct read)",
    amplicon = "400-800 bp",
    loci = 1L,
    verdict = "usable Sanger amplicon",
    gates = KASP_GATES_SANGER,
    relax = KASP_RELAX_SANGER,
    tiers = KASP_TIERS_SANGER,
    gate_short = KASP_GATE_SHORT_SANGER,
    gate_limit = c(
      gate_sanger_amplicon = glue::glue("amplicon within [{amp[1]}, {amp[2]}] bp"),
      gate_sanger_tm = glue::glue("both primers Tm within [{tm[1]}, {tm[2]}] C"),
      gate_sanger_dist = glue::glue("variant >= {dist} bp from each read start"),
      gate_sanger_unique = "exactly one in-silico amplicon"
    ),
    gate_floor_label = KASP_GATE_FLOOR_LABEL_SANGER,
    gate_hard_label = KASP_GATE_HARD_LABEL_SANGER,
    stage_label = "Sanger amplicon"
  )
}
