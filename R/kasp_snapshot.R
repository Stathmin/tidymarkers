# =============================================================================
#  kasp_snapshot.R  --  SNaPshot / single-base extension as a marker scheme
# =============================================================================
#  SNaPshot reads the genotype from a *single-base extension*: an extension (or
#  "minisequencing") primer anneals so that its 3' terminus sits immediately 5'
#  of the variant, and the polymerase adds exactly one fluorescent ddNTP - the
#  allele's base - before the products are sized by capillary electrophoresis.
#  There is no allele-specific primer and no 3'-terminal mismatch: the variant is
#  read by *base identity*, so the method tolerates nearly any locus and is
#  co-dominant (homozygotes give one peak, heterozygotes two).  The design
#  therefore optimises two things - the exact 3' placement of the extension
#  primer on the variant, and a clean flanking pair that makes the template.
#
#  Design rules (web-researched; do not invent thresholds):
#    * Extension primer 3' end immediately 5' of the variant, so the first ddNTP
#      added reports the allele.
#      - Namdar Aligoodarzi P, et al. Iran Biomed J 2023;27(1):46-57,
#        doi:10.52547/ibj.3856 (PMC9971712): "A primer melting temperature was
#        calculated ... by annealing one base before the SNP to be analyzed.
#        Therefore, the polymorphisms would be identified by the base after the
#        minisequencing primer for each SNP."
#      - Fanis P, et al. BMC Genomics 2014;15:108, doi:10.1186/1471-2164-15-108
#        (PMC3922441): "Minisequencing primers were designed with their 3' base
#        corresponding to the base immediately before the SNP to be investigated."
#    * Extension-primer length/Tm: 15-25 nt at Tm 49-53 C (Namdar Aligoodarzi);
#      "there is little flexibility in the positioning of the primer, although
#      parameters like the length of the primer and its orientation are more
#      amenable to manipulation" (Fanis).  ASSUMPTION: we take the task's window
#      (Tm 50-55 C, ~18-25 nt) as the gate and widen it by 5 C on each side (the
#      ladder's tm floor: [45, 60] C) when nothing clears, because the package
#      primer_tm is salt-corrected (Mg/dNTPs) and runs hotter than the Promega
#      formula those papers quote.
#    * The extension primer must be free of stable hairpins/self-dimers: "stable
#      oligonucleotide secondary structures can affect signals" (Fanis); hairpins
#      and dimers "were checked" (Namdar Aligoodarzi).
#    * The extension primer must be *unique*: a second perfect binding site
#      extends elsewhere and adds a spurious peak to the electropherogram, so the
#      design requires exactly one perfect genomic site for it.
#    * Amplification: a flanking pair over a short amplicon.  Published amplicons
#      span ~60-882 bp (Fanis 60-649; Namdar Aligoodarzi 143-882).  ASSUMPTION:
#      the task's 100-300 bp target is a practical default, not a hard rule - hence
#      a *soft* gate whose floor widens it to [50, 350] bp.
#    * Readout: capillary electrophoresis of the extended products; the kit
#      multiplexes up to 10 SNPs (Thermo Fisher, SNaPshot Multiplex Kit,
#      cat. 4323151, https://www.thermofisher.com/order/catalog/product/4323151).
#  The literature hit-list came from the Europe PMC REST search
#  (https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=SNaPshot
#  %20minisequencing%20primer%20design&format=json); DuckDuckGo's HTML endpoint
#  (https://duckduckgo.com/html/?q=SNaPshot+minisequencing+primer+design) was
#  bot-blocked at the time.
#
#  Every declared gate is evaluated by the designer (C7): the candidate scans
#  enforce the hard gates and `grade_one_marker()` records them through the shared
#  relaxation ladder (least grave first: amplicon, tm, structure), so the
#  scheme's decision graph and its engine read the same thresholds.
# =============================================================================

# Default design targets for a SNaPshot marker (sources in the file header).
SNAPSHOT_AMPLICON <- c(100L, 300L)
SNAPSHOT_AMP_LENGTHS <- 18:25
SNAPSHOT_EXT_LENGTHS <- 18:25
SNAPSHOT_EXT_TM <- c(50, 55)
SNAPSHOT_EXT_GC <- c(0.30, 0.70)

# How far the ladder may lower the three soft SNaPshot gates (margin units): the
# amplicon window by this many bp on each side, the extension Tm window by this
# many C on each side, and this much extra probe/extension self-structure.
SNAPSHOT_AMPLICON_RELAX <- 50L
SNAPSHOT_TM_RELAX <- 5
SNAPSHOT_STRUCT_RELAX <- 3

#' Design a SNaPshot / single-base-extension marker at a variant
#'
#' The minisequencing analogue of [assemble_kasp_assay()]: a flanking pair makes
#' an amplicon over the variant, and an **extension primer** is placed so that
#' its 3' terminus sits exactly one base 5' of the variant.  The first ddNTP
#' added therefore reports the allele by base identity (co-dominant), without the
#' 3'-terminal mismatch KASP relies on.  The extension primer is kept short with
#' a low Tm (so it anneals but extends specifically), free of stable hairpins and
#' self-dimers, and the amplicon is kept short.  Single-base substitutions only.
#' @param genome a `kasp_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`); single-base only.
#' @param flank template half-width.
#' @param amp_lengths primer lengths to try for the flanking pair.
#' @param min_amplicon,max_amplicon the accepted amplicon window (bp).
#' @param amp_tm,amp_gc accepted windows for the flanking pair.
#' @param ext_lengths extension-primer lengths to try.
#' @param ext_tm accepted extension-primer Tm window (C); widened once if no
#'   candidate clears it (see the file header).
#' @param ext_gc accepted extension-primer GC window (fraction).
#' @param limits thresholds list (see `KASP_LIMITS`); its `hairpin_dg` /
#'   `homodimer_dg` floor the structure gate.
#' @return a `kasp_snapshot` list (`snp`, `amplification`, `extension`,
#'   `amplicon`, `compatibility` - the joint dimer/hairpin screen over the pair
#'   and the extension primer), or `NULL` (with a message) when no pair/primer
#'   passes.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_snapshot_marker(genome, kasp_example("SNP")[1, ])
#' }
#' @family KASP SNaPshot
#' @export
design_snapshot_marker <- function(genome, snp, flank = 250L,
                                   amp_lengths = SNAPSHOT_AMP_LENGTHS,
                                   min_amplicon = SNAPSHOT_AMPLICON[1L],
                                   max_amplicon = SNAPSHOT_AMPLICON[2L],
                                   amp_tm = KASP_LIMITS$com_tm,
                                   amp_gc = KASP_LIMITS$gc,
                                   ext_lengths = SNAPSHOT_EXT_LENGTHS,
                                   ext_tm = SNAPSHOT_EXT_TM,
                                   ext_gc = SNAPSHOT_EXT_GC,
                                   limits = KASP_LIMITS) {
  ref <- as.character(snp$ref)[1L]
  alt <- as.character(snp$alt)[1L]
  if (nchar(ref) != 1L || nchar(alt) != 1L) {
    cli::cli_inform(gettext("SNaPshot single-base extension is implemented for single-base substitutions only."))
    return(NULL)
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  si <- tmpl$snp_i
  # the amplicon rung of the ladder: the declared window first, widened by the
  # floor only when nothing clears it (so the soft gate is genuinely reachable)
  amp_in <- function(lo, hi) {
    design_flanking_pair(tmpl$seq, si, si, amp_lengths, lo, hi, amp_tm, amp_gc)
  }
  amp <- amp_in(min_amplicon, max_amplicon) %||%
    amp_in(min_amplicon - SNAPSHOT_AMPLICON_RELAX, max_amplicon + SNAPSHOT_AMPLICON_RELAX)
  if (is.null(amp)) {
    cli::cli_inform(gettext("No amplification pair passed the Tm/GC and amplicon-length filters."))
    return(NULL)
  }
  ext <- design_extension_primer(tmpl$seq, si, ext_lengths, ext_tm, ext_gc, limits)
  if (is.null(ext)) {
    # the tm rung: the same documented relaxation, declared as the gate's floor
    ext <- design_extension_primer(
      tmpl$seq, si, ext_lengths,
      ext_tm + c(-SNAPSHOT_TM_RELAX, SNAPSHOT_TM_RELAX),
      ext_gc, limits
    )
  }
  if (is.null(ext)) {
    cli::cli_inform(gettext("No extension primer ended immediately 5' of the variant within the Tm/GC filters."))
    return(NULL)
  }
  off <- snp$pos - si
  out <- list(
    snp = snp,
    amplification = list(
      forward = amp$forward, reverse = amp$reverse,
      start = amp$start + off, end = amp$end + off, product = amp$product,
      tm_f = amp$tm_f, tm_r = amp$tm_r
    ),
    extension = list(
      primer = ext$primer, strand = ext$strand,
      start = ext$start + off, end = ext$end + off,
      tm = ext$tm, gc = ext$gc, placement = ext$placement,
      hairpin_dg = ext$hairpin_dg, homodimer_dg = ext$homodimer_dg,
      structure_ok = ext$structure_ok
    ),
    amplicon = tibble::tibble(start = amp$start + off, end = amp$end + off)
  )
  # pair plus extension primer share the tube: screen the joint set
  out$compatibility <- tube_compatibility(c(
    amp_forward = out$amplification$forward, amp_reverse = out$amplification$reverse,
    extension = out$extension$primer
  ))

  # ---- C7: the declared gates are all evaluated on the winner ----------------
  d <- scheme_def_snapshot()
  # placement is structural, but verify it on the winner rather than assume it
  placement_ok <- if (identical(ext$strand, "+")) ext$end == si - 1L else ext$start == si + 1L
  ext_gc_obs <- gc_frac(ext$primer)
  comp_ok <- ext_gc_obs >= SNAPSHOT_EXT_GC[1L] && ext_gc_obs <= SNAPSHOT_EXT_GC[2L] &&
    !has_run(ext$primer) && !has_repeat(ext$primer)
  tm_win <- d$ext_tm_window
  amp_win <- d$amplicon_window
  graded <- tibble::tibble(
    m_placement = if (placement_ok) 1 else -1,
    m_comp = if (comp_ok) 1 else -1,
    m_unique = if (exact_oligo_hits(ext$primer, genome) == 1L) 1 else -1,
    m_structure = pmin(ext$hairpin_dg - limits$hairpin_dg, ext$homodimer_dg - limits$homodimer_dg),
    m_tm = min(ext$tm - tm_win[1L], tm_win[2L] - ext$tm),
    m_amplicon = pmin(amp$product - amp_win[1L], amp_win[2L] - amp$product)
  )
  gate_margin <- c(
    gate_snapshot_placement = "m_placement", gate_snapshot_comp = "m_comp",
    gate_snapshot_unique = "m_unique", gate_snapshot_structure = "m_structure",
    gate_snapshot_tm = "m_tm", gate_snapshot_amplicon = "m_amplicon"
  )
  sel <- grade_one_marker(
    graded,
    gates = d$gates, relax_order = d$relax,
    tiers = c("m_tm", "m_structure", "m_amplicon"), # tier1..3 in scheme_def_snapshot()
    gate_margin = gate_margin, gate_floor = d$gate_floor
  )
  if (is.null(sel)) {
    cli::cli_inform(gettext(
      "The SNaPshot candidate failed a hard gate (extension placement, composition or uniqueness)."
    ))
    return(NULL)
  }
  out$selection <- sel
  structure(out, class = "kasp_snapshot")
}

#' @export
print.kasp_snapshot <- function(x, ...) {
  end3 <- if (x$extension$strand == "+") x$extension$end else x$extension$start
  cat(
    "<kasp_snapshot>", x$snp$chrom, ":", x$snp$pos, " ",
    x$snp$ref, ">", x$snp$alt, "\n",
    sep = ""
  )
  cat("  amplicon:", x$amplification$product, "bp\n")
  cat("  F:", x$amplification$forward, "\n  R:", x$amplification$reverse, "\n")
  cat(
    "  extension (", x$extension$strand, "): ", x$extension$primer,
    " (Tm ", round(x$extension$tm, 1), " C, 3' end at ", end3, ")\n",
    sep = ""
  )
  invisible(x)
}

#' The best extension primer whose 3' end abuts the variant
#'
#' Placement fixes the candidate set: on the plus strand the primer ends at
#' `si - 1`, on the minus strand (reverse complement) it starts at `si + 1`.
#' Candidates are GC/Tm-screened by [primer_candidates()], structure-screened
#' with OligoArrayAux when available, then ranked by Tm centrality.
#' @keywords internal
#' @noRd
design_extension_primer <- function(seq, si, lengths, tm_range, gc_range, limits) {
  s <- as.character(seq)
  n <- nchar(s)
  short <- min(lengths)
  seg <- function(a, b) if (a >= 1L && b <= n) stringr::str_sub(s, a, b) else ""
  # The placement pins each orientation to a nested candidate set (all lengths
  # share the 3' end), so if the *shortest* primer is a run/repeat every longer
  # one is too and primer_candidates() would be handed an empty vector - skip it.
  plus_seg <- seg(si - short, si - 1L)
  minus_seg <- seg(si + 1L, si + short)
  plus_ok <- nzchar(plus_seg) && !has_run(plus_seg) && !has_repeat(plus_seg)
  minus_ok <- nzchar(minus_seg) && !has_run(minus_seg) && !has_repeat(minus_seg)
  empty <- tibble::tibble(
    start = integer(), len = integer(), end = integer(),
    primer = character(), tm = numeric(), gc = numeric(), strand = character()
  )
  plus <- if (plus_ok) {
    primer_candidates(
      seq, seq.int(max(1L, si - max(lengths)), max(1L, si - min(lengths))),
      lengths, "+", tm_range, gc_range
    ) |>
      dplyr::filter(.data$end == si - 1L) |>
      dplyr::mutate(strand = "+")
  } else {
    empty
  }
  minus <- if (minus_ok) {
    primer_candidates(seq, si + 1L, lengths, "-", tm_range, gc_range) |>
      dplyr::mutate(strand = "-")
  } else {
    empty
  }
  cand <- dplyr::bind_rows(plus, minus)
  if (!nrow(cand)) {
    return(NULL)
  }
  if (has_oligoarrayaux()) {
    cand$hairpin_dg <- oligo_dg_vec(cand$primer)
    cand$homodimer_dg <- oligo_dg_vec(cand$primer, cand$primer)
  } else {
    cand$hairpin_dg <- NA_real_
    cand$homodimer_dg <- NA_real_
  }
  cand$structure_ok <- (is.na(cand$hairpin_dg) | cand$hairpin_dg >= limits$hairpin_dg) &
    (is.na(cand$homodimer_dg) | cand$homodimer_dg >= limits$homodimer_dg)
  if (any(cand$structure_ok)) {
    cand <- cand[cand$structure_ok, , drop = FALSE]
  }
  cand$placement <- 1L # the 3' end sits one base from the variant
  cand$dtm <- abs(cand$tm - mean(tm_range))
  cand[order(cand$dtm, cand$len, cand$strand), , drop = FALSE][1L, , drop = FALSE]
}

#' The SNaPshot scheme definition (gates, tiers, labels)
#' @keywords internal
#' @noRd
scheme_def_snapshot <- function() {
  gates <- c(
    "gate_snapshot_placement", "gate_snapshot_comp", "gate_snapshot_unique",
    "gate_snapshot_structure", "gate_snapshot_tm", "gate_snapshot_amplicon"
  )
  list(
    scheme = "snapshot",
    family = "snapshot",
    label = "SNaPshot / single-base extension (minisequencing)",
    discrimination = "one fluorescent ddNTP added at the variant (base identity)",
    readout = "capillary electrophoresis (dye colour + product size)",
    zygosity = "co-dominant",
    amplicon = "bounded (<= 300 bp)",
    loci = 1L,
    verdict = "usable SNaPshot marker",
    gates = gates,
    # least grave first: the amplicon window is a practical default, the Tm
    # window is salt-corrected, and a hairpin/dimer is the costliest deviation
    relax = c("gate_snapshot_amplicon", "gate_snapshot_tm", "gate_snapshot_structure"),
    tiers = c(
      tier1 = "extension-primer Tm centrality",
      tier2 = "extension self-structure (hairpin / homodimer)",
      tier3 = "amplicon length (near the ideal)"
    ),
    gate_short = c(
      gate_snapshot_placement = "extension 3' placement",
      gate_snapshot_comp = "extension composition",
      gate_snapshot_unique = "extension primer uniqueness",
      gate_snapshot_structure = "extension structure",
      gate_snapshot_tm = "extension primer Tm",
      gate_snapshot_amplicon = "amplicon length"
    ),
    gate_limit = c(
      gate_snapshot_placement = "extension 3' end immediately 5' of the variant",
      gate_snapshot_comp = glue::glue(
        "GC {SNAPSHOT_EXT_GC[1L] * 100}-{SNAPSHOT_EXT_GC[2L] * 100} %, no run >= 4, no SSR"
      ),
      gate_snapshot_unique = "exactly one perfect extension-primer site in the genome",
      gate_snapshot_structure = glue::glue("hairpin/homodimer dG37 >= {KASP_LIMITS$hairpin_dg}"),
      gate_snapshot_tm = glue::glue("extension primer Tm within [{SNAPSHOT_EXT_TM[1]}, {SNAPSHOT_EXT_TM[2]}] C"),
      gate_snapshot_amplicon = glue::glue("amplicon within [{SNAPSHOT_AMPLICON[1]}, {SNAPSHOT_AMPLICON[2]}] bp")
    ),
    # the numeric windows the engine scores against (so gate margins and the
    # declared limits cannot drift) and the ladder's floors in margin units
    amplicon_window = SNAPSHOT_AMPLICON,
    ext_tm_window = SNAPSHOT_EXT_TM,
    ext_gc_window = SNAPSHOT_EXT_GC,
    gate_floor = c(
      gate_snapshot_amplicon = -SNAPSHOT_AMPLICON_RELAX,
      gate_snapshot_tm = -SNAPSHOT_TM_RELAX,
      gate_snapshot_structure = -SNAPSHOT_STRUCT_RELAX
    ),
    gate_floor_label = c(
      gate_snapshot_amplicon = glue::glue("widen the amplicon window by {SNAPSHOT_AMPLICON_RELAX} bp"),
      gate_snapshot_tm = glue::glue("widen the Tm window by {SNAPSHOT_TM_RELAX} C each side"),
      gate_snapshot_structure = glue::glue("allow {SNAPSHOT_STRUCT_RELAX} kcal/mol more stable")
    ),
    gate_hard_label = c(
      gate_snapshot_placement = "the first ddNTP reads the variant only if the 3' end abuts it",
      gate_snapshot_comp = "a run or extreme GC content stops the extension primer annealing",
      gate_snapshot_unique = "a second binding site extends elsewhere and adds a spurious peak"
    ),
    stage_label = "SNaPshot extension"
  )
}
