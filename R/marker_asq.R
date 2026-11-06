# =============================================================================
#  marker_asq.R  --  ASQ (allele-specific qPCR): a universal-probe readout of KASP
# =============================================================================
#  ASQ reads the same competitive allele-specific extension as KASP, but its
#  distinguishing feature is the geometry of the allele-specific primer:
#   * the variant sits at the penultimate (or antepenultimate) 3' base of the
#     allele-specific primer, not the terminal base KASP uses; and
#   * the reporter is a set of universal probes (UP1-UP4, FAM/HEX/Cy3/Cy5) read
#     against a *single* common quencher (Uni-Q), i.e. a master-mix / instrument
#     (panel) property, up to 4-plex.
#
#  Design rules (researched, not invented).  Kalendar R, Shustov AV, Akhmetollayev
#  I, Kairov U (2022) "Designing Allele-Specific Competitive-Extension PCR-Based
#  Assays for High-Throughput Genotyping and Gene Characterization", Front Mol
#  Biosci 9:773956, doi:10.3389/fmolb.2022.773956 (open access, PMC8921500):
#   * "each allele-specific primer is designed such that the penultimate base
#     before the primer's 3' end base is positioned at the SNP site" - the tool
#     "automatically places a variable base at the second (penultimate) position
#     of the primer's 3' terminus", or the third, chosen thermodynamically;
#   * the FRET cassette becomes "two or more universal probes (UP-1-4) with
#     corresponding tags and different fluorescent dyes at the 5' end and a single
#     common UP with a quencher at the 3' end (Uni-Q)" - "4 in ASQ, 2 in KASP";
#   * the allele-specific primers "were selected to obtain a Tm close to 63 C
#     (range 62-64 C)"; the universal probes are designed higher (66-68 C) and the
#     quencher lower (53-55 C), the signal being read at a second 55 C step;
#   * "amplicons ranged in size from 100 to 200 bases" and encompassed the SNP.
#
#  The engine is the KASP one, re-parameterised: [design_asq_marker()] calls
#  [assemble_kasp_assay()] with `variant_offset = 1` (the geometric difference),
#  ASQ's own acceptance windows ([asq_limits()]) and the 63 C AS-Tm target.  The
#  shared AS/pair selection stages therefore grade ASQ with ASQ's thresholds, and
#  the penultimate placement is enforced by the geometry itself, not by a
#  decorative gate.
#
#  The reporter/plexity rules are panel and format properties - they are not
#  computable from a candidate, so they are not selection gates (C7/C8).  The
#  `asq` stage in the decision schema is documentation of the geometry rule and
#  its 63 C ranking dimension; it declares no gate of its own because the shared
#  AS/pair stages already compute ASQ's real gates under ASQ's limits.
# =============================================================================

# ---- ASQ acceptance windows (Kalendar 2022) ---------------------------------
ASQ_AMPLICON <- c(100, 200) # bp   ("amplicons ranged in size from 100 to 200")
ASQ_AS_TM_LIMITS <- c(62, 64) # C  ("a Tm close to 63 C (range 62-64 C)")
ASQ_AS_TM_TARGET <- 63 # C       (the 62-64 C window's optimum)
# The COM primer is not printed in the paper; it only has to clear the shared
# "COM Tm >= AS Tm" gate.  The window is the ASQ AS window extended upward by the
# universal-probe margin, so a COM primer reached by the engine is always inside
# the ASQ AS-Tm range or above it.
ASQ_COM_TM_LIMITS <- c(62, 68) # C

#' ASQ acceptance thresholds (the ASQ twin of `KASP_LIMITS`)
#'
#' The KASP thresholds with ASQ's own amplicon and AS-primer-Tm windows: 100-200
#' bp and 62-64 C (Kalendar et al. 2022).  An assay assembled by
#' [design_asq_marker()] carries this object as `$limits`, and every acceptance
#' check and verdict reads it, so ASQ is never graded against KASP's 50-100 bp /
#' 55-65 C optimum.
#' @param amplicon the amplicon window (bp).
#' @param as_tm the AS-primer Tm window (C).
#' @param com_tm the COM-primer Tm window (C).
#' @return a named list, the same shape as `KASP_LIMITS`.
#' @examples
#' asq_limits()
#' @family design
#' @export
asq_limits <- function(amplicon = ASQ_AMPLICON, as_tm = ASQ_AS_TM_LIMITS,
                       com_tm = ASQ_COM_TM_LIMITS) {
  out <- KASP_LIMITS
  out$amplicon <- as.numeric(amplicon)
  out$as_tm <- as.numeric(as_tm)
  out$com_tm <- as.numeric(com_tm)
  out
}

#' Design an ASQ / PACE marker (variant at the penultimate 3' base)
#'
#' ASQ (allele-specific qPCR, also marketed as PACE) uses the same AS + common
#' primer pair as KASP, but places the variant at the **penultimate** 3' base of
#' the allele-specific primer (Kalendar et al. 2022).  This is the KASP engine
#' with that geometry, ASQ's windows ([asq_limits()]) and the 63 C AS-Tm target,
#' so it returns the same `marker_assay` object every other chemistry returns.
#'
#' The returned assay carries `as_variant_offset = 1L`: the variant block sits one
#' base inside the AS primer's 3' terminus, and both allele primers share their
#' terminal base.  Only substitutions are designed (the penultimate placement of
#' a length-changing indel junction is a different geometry); an indel target is
#' doomed rather than silently designed as a KASP.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp a one-row tibble with `chrom`, `pos`, `ref`, `alt` (and optionally
#'   `alt2`).
#' @param variants optional variant table (see [assemble_kasp_assay()]).
#' @param haplotypes optional tidy haplotype table (see [assemble_kasp_assay()]).
#' @param flank template half-width (bp); ASQ products reach 200 bp, so the
#'   default is larger than KASP's 150.
#' @param BPPARAM optional `BiocParallel` backend for the genome-wide steps
#'   (see [assemble_kasp_assay()]); `NULL` runs them serially.
#' @param ... further arguments passed to [assemble_kasp_assay()].
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' snp <- marker_example("SNP")[1, c("chrom", "pos", "ref", "alt")]
#' assay <- design_asq_marker(genome, snp)
#' assay$as_variant_offset
#' }
#' @family design
#' @export
design_asq_marker <- function(genome, snp, variants = NULL, haplotypes = NULL,
                              flank = 200L, BPPARAM = NULL, ...) {
  # the penultimate placement is defined for a substitution block; a
  # length-changing indel junction is a different geometry, so doom it rather
  # than silently designing a terminal (KASP) indel under the ASQ name
  if (nchar(as.character(snp$ref)[1L]) != nchar(as.character(snp$alt)[1L])) {
    tmpl <- snp_template(genome, snp$chrom[1L], snp$pos[1L], flank)
    return(marker_doomed(
      gettext("ASQ designs substitutions; a length-changing indel is not placed at the penultimate base."),
      "as", snp, tmpl
    ))
  }
  lim <- asq_limits()
  assemble_kasp_assay(
    genome, snp,
    variants = variants, haplotypes = haplotypes, flank = flank,
    scheme = "asq", limits = lim, variant_offset = 1L,
    tm_range = lim$as_tm, com_tm_range = lim$com_tm, com_gaps = 18:200L,
    as_tm_target = ASQ_AS_TM_TARGET, com_tm_target = mean(lim$com_tm),
    # the 62-64 C AS window is reached by longer primers than KASP's 20-27 nt
    as_lengths = 22:35L,
    BPPARAM = BPPARAM, ...
  )
}

# ---- the ASQ selection stage (documentation-only) ---------------------------
# The real gates are the shared AS/pair stages, now graded with asq_limits().
# This stage adds no gate of its own: the geometry rule (penultimate placement)
# is enforced by the engine, and the only ASQ-specific ranking dimension is the
# AS primer Tm toward the 63 C optimum.  No relaxation ladder either (C6's
# relaxation order is empty): there is nothing to relax.
KASP_GATES_ASQ <- character(0)
KASP_RELAX_ASQ <- character(0)
KASP_TIERS_ASQ <- c(tier1 = "AS primer Tm (toward the 63 C ASQ optimum)")

# This stage contributes no gate, so it has no entry in the shared gate-label
# lookups (assembled in R/marker_select.R); the empty vectors keep those
# concatenations valid without a dangling label for a gate that does not exist.
KASP_GATE_LABEL_ASQ <- character(0)
KASP_GATE_SHORT_ASQ <- character(0)
KASP_GATE_LIMIT_ASQ <- character(0)
KASP_GATE_HARD_LABEL_ASQ <- character(0)
