# =============================================================================
#  kasp_hrm.R  --  HRM (high-resolution melting) as a marker scheme
# =============================================================================
#  HRM reads the genotype from the *melting curve of a short amplicon*: the
#  alleles differ in how the duplex melts, and heterozygotes show a heteroduplex
#  shape.  It is the scheme for the loci KASP cannot take - a variant in the AS
#  3' window, a poly-A junction - because it anchors the variant in the
#  *amplicon*, not in a primer, and *tolerates neighbours* (the curve separates
#  them).  So there is no allele-specific primer: one flanking pair straddles the
#  variant, the amplicon is kept SHORT (a short melt domain sharpens the curve -
#  Premier Biosoft / Thermo Fisher HRM guidance), and the pair is Tm-balanced.
#  The melt curve itself is not simulated (that needs uMELT/PolyMelt-class
#  models), so this returns a well-formed amplicon, not a predicted profile.
# =============================================================================

#' Design an HRM marker at a variant
#'
#' The melting-curve analogue of [assemble_kasp_assay()]: a single flanking pair
#' straddles the variant, which sits in the melt domain and under neither primer.
#' The design optimises the two things that govern HRM - a **short** amplicon and
#' a **Tm-balanced** pair - and, because HRM tolerates neighbouring variants,
#' imposes no clean 3' window.  This is why HRM is the fallback for the loci KASP
#' dooms (a variant in the AS 3' window, an AT-rich / poly-A junction).
#'
#' It does **not** simulate the melting curve (that needs uMELT/PolyMelt-class
#' models), so it returns a well-formed amplicon, not a predicted melt profile.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param flank template half-width.
#' @param min_product,max_product the accepted amplicon window (HRM wants short).
#' @param lengths primer lengths to try.
#' @param tm_range,gc_range accepted primer windows.
#' @return a `kasp_hrm` list (`snp`, `primers`, `amplicon`), or `NULL` (with a
#'   message) when no flanking pair passes the filters.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_hrm_marker(genome, kasp_example("SNP")[1, ])
#' }
#' @family KASP HRM
#' @export
design_hrm_marker <- function(genome, snp, flank = 150L,
                              min_product = 40L, max_product = 150L, lengths = 18:25,
                              tm_range = KASP_LIMITS$com_tm, gc_range = KASP_LIMITS$gc) {
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  # the variant sits *between* the primers, so the whole of it is in the melt
  # domain and never under a primer (lo = hi = the variant position)
  si <- tmpl$snp_i
  primers <- design_flanking_pair(
    tmpl$seq, si, si, lengths, min_product, max_product, tm_range, gc_range
  )
  if (is.null(primers)) {
    cli::cli_inform(gettext("No short flanking primer pair passed the Tm/GC filters for HRM."))
    return(NULL)
  }
  out <- list(
    snp = snp, primers = primers,
    amplicon = tibble::tibble(
      start = snp$pos + primers$start - si,
      end = snp$pos + primers$end - si
    )
  )
  structure(out, class = "kasp_hrm")
}

#' @export
print.kasp_hrm <- function(x, ...) {
  cat(
    "<kasp_hrm>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  amplicon:", x$primers$product, "bp (melt domain)\n")
  cat("  F:", x$primers$forward, "\n  R:", x$primers$reverse, "\n")
  invisible(x)
}
