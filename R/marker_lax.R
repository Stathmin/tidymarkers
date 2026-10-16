# =============================================================================
#  marker_lax.R  --  lax KASP: one ALT primer per haplotype (linked variants OK)
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
#
#  Design rules (researched, not invented).  The method basis is the KASP manual
#  (README section 3, "Краткий мануал по подбору KASP": LGC/KBioscience practice
#  plus EasyKASP, Zhang, Yang & Wen, BMC Bioinformatics 2025, 26:292, where the
#  optimal ranges come from 500 validated KASP tests), this file's own header, and
#  the linkage section of ROADMAP section H.
#
#    * What discriminates REF from an ALT haplotype: the AS primer's 3'-terminal
#      base, which Taq cannot extend past a mismatch, plus a deliberate second
#      mismatch 1-3 nt in from the 3' end that destabilises the wrong allele
#      further (README section 3, "Аллель-специфичные праймеры (AS)").  A *lax*
#      ALT primer keeps that 3' anchor and additionally carries the ALT
#      haplotype's linked alleles in its body, so against REF it mismatches at the
#      anchor *and* at every baked linked position.
#    * Keeping the promise is therefore a *coverage* question: "the site gets one
#      ALT primer per ALT haplotype" (this file's header, lines 6-9).  A haplotype
#      that no designed primer binds cannot be called, so it is a HARD gate (C8).
#    * Linkage is only a licence to bake a variant into the ALT primer while it
#      really travels with the target.  README section 3 disqualifies a primer
#      whose last ~5 nt cover "другой частый SNP/индел" (another *common*
#      SNP/indel), and ROADMAP section H states the premise: a variant in the AS
#      zone *together with* the target is an ally, not a sin.  The assembler acts
#      on that licence by excusing the baked positions from the AS footprint
#      (`as_exclude`, R/marker_assemble.R), so the designer must verify it: an
#      allele the *wild type* also carries makes the REF primer mismatch those
#      templates (wild-type dropout) - the second gate.
#    * The manual screens EVERY pair in the tube (README section 3, "Вторичная
#      структура и комплементарность": "Отбраковываем сильные димеры любой пары
#      (AS1-COM, AS2-COM, AS1-AS2, само-CM)") but tolerates a strong AS1-AS2 pair
#      because those primers compete for one site.  For lax the haplotype ALT
#      primers all share the target 3' base, so the risk the dimer screen cannot
#      see is *cross-priming*: a primer is only usable if its 3'-proximal window
#      separates it from the other groups, or every group fires on every ALT
#      template and the haplotype resolution collapses - the third gate.
# =============================================================================

# The lax decision stage (declared here, read by the shared registry through
# scheme_def_lax()).  It is a *set-level* stage: every per-haplotype ALT primer
# has already been through the shared `as`/`pair` stages (one
# assemble_kasp_assay() per group), so what is left to decide is whether the
# *tube* keeps the scheme's promise.  Gate order = grave / cheap first:
#
#   1. gate_lax_coverage      HARD  every ALT haplotype is bound by a designed ALT
#                                   primer (no primer = no call; C8)
#   2. gate_lax_alt_specific  soft  the linked alleles baked into those primers
#                                   are ALT-specific: no wild-type (REF) haplotype
#                                   carries them, or the REF primer drops those
#                                   templates out
#   3. gate_lax_cross_priming soft  the ALT primers are separable from one another
#                                   inside their 3'-proximal window, so they do not
#                                   cross-prime each other's haplotypes
#
# Relaxation order is the opposite (least grave first): losing a haplotype group
# to cross-priming still calls REF/ALT, while a wild-type dropout corrupts the
# call itself, so the milder risk is tried first (C6 - the two orders differ).
KASP_GATES_LAX <- c("gate_lax_coverage", "gate_lax_alt_specific", "gate_lax_cross_priming")
KASP_RELAX_LAX <- c("gate_lax_cross_priming", "gate_lax_alt_specific")

# The one ranking dimension the set level adds: the population fraction of ALT
# haplotypes the tube can actually *resolve*.  It is a margin, not a gate (the
# cross-priming gate above is the count that fails), and it does not restate the
# hard coverage gate - a covered haplotype can still be unresolvable when its
# group shares its 3'-proximal window with another (C14).  One tube is one
# candidate, so the tier reports quality rather than breaking a tie.
KASP_TIERS_LAX <- c(tier1 = "resolvable ALT haplotype frequency")

# The tolerance a baked linked allele may also have in the wild type: the same bar
# as a footprint variant's expected dropout (KASP_LIMITS), because that is exactly
# what the REF primer then loses.
KASP_LAX_ALT_SPECIFIC_TOL <- KASP_LIMITS$footprint_dropout %||% 0.05

# Floor per soft gate (margin units; 0 = hard).  The wild-type-leakage floor is one
# more tolerance, the cross-priming floor tolerates exactly one collapsible group
# pair - both bounded, never -Inf (C8): a relaxed winner is still a warned, usable
# tube, not a silent one.
KASP_GATE_FLOOR_LAX <- c(
  gate_lax_coverage = 0,
  gate_lax_alt_specific = -KASP_LAX_ALT_SPECIFIC_TOL,
  gate_lax_cross_priming = -1
)

KASP_GATE_SHORT_LAX <- c(
  gate_lax_coverage = "haplotype coverage",
  gate_lax_alt_specific = "linked alleles ALT-specific",
  gate_lax_cross_priming = "haplotype cross-priming"
)
KASP_GATE_LIMIT_LAX <- c(
  gate_lax_coverage = "every ALT haplotype is bound by a designed ALT primer",
  gate_lax_alt_specific = glue::glue(
    "no wild-type haplotype carries a baked linked allele (leakage <= {KASP_LAX_ALT_SPECIFIC_TOL})"
  ),
  gate_lax_cross_priming = "each ALT primer is separable from the others in its 3'-proximal window"
)
KASP_GATE_FLOOR_LABEL_LAX <- c(
  gate_lax_alt_specific = "tolerate wild-type leakage up to 2x tolerance",
  gate_lax_cross_priming = "tolerate one collapsible haplotype-group pair"
)
KASP_GATE_HARD_LABEL_LAX <- c(
  gate_lax_coverage = "a haplotype no ALT primer binds cannot be called"
)
KASP_GATE_LABEL_LAX <- c(
  gate_lax_coverage = "an ALT haplotype is left without an ALT primer",
  gate_lax_alt_specific = "a baked linked allele is also carried by the wild type",
  gate_lax_cross_priming = "two haplotype ALT primers cross-prime each other"
)

# The numbers behind the stage, declared next to it (C9: real, parameterised
# values, never a `[lo, hi]` placeholder).
KASP_LAX_CROSS_PRIMING_FLOOR <- 1L # collapsible group pairs the soft floor tolerates

#' The fraction of wild-type (REF) templates a baked linked set drops out
#'
#' The lax engine bakes an ALT haplotype's co-segregating alleles into its ALT
#' primer; `assemble_kasp_assay()` then excuses those positions from the AS
#' footprint and 3'-window counts (`as_exclude`), because they are assumed to ride
#' with the ALT haplotype.  That licence only holds while the *wild type* does not
#' carry them too: a REF haplotype with the baked allele is a template its REF
#' primer now mismatches - wild-type dropout, and a false REF/ALT reading.  This is
#' that check for one tube.
#' @param haplotypes a tidy haplotype table (`haplotype`, `freq`, `chrom`, `pos`,
#'   `ref`, `alt`).
#' @param snp one-row target (`pos`, `ref`, `alt`).
#' @param linked_positions the variant positions the kept ALT primers baked in.
#' @return a fraction in `[0, 1]` (`0` when no REF haplotype carries one).
#' @keywords internal
#' @noRd
lax_ref_leakage <- function(haplotypes, snp, linked_positions) {
  h <- normalise_haplotypes(haplotypes)
  linked_positions <- unique(as.integer(linked_positions))
  linked_positions <- linked_positions[!is.na(linked_positions)]
  if (!length(linked_positions)) {
    return(0)
  }
  alt_hap <- unique(h$haplotype[h$pos == snp$pos & h$alt == snp$alt])
  ref_hap <- setdiff(unique(h$haplotype), alt_hap)
  if (!length(ref_hap)) {
    return(0)
  }
  freqs <- h |>
    dplyr::distinct(.data$haplotype, .data$freq) |>
    tibble::deframe()
  # the alleles the ALT primers bake: each ALT carrier's own allele at a linked
  # position - exactly what design_as_primers(linked=) writes into the primer
  baked <- h[h$haplotype %in% alt_hap & h$pos %in% linked_positions, c("pos", "alt")]
  baked <- baked[!duplicated(baked), , drop = FALSE]
  if (!nrow(baked)) {
    return(0)
  }
  ref <- h[h$haplotype %in% ref_hap & h$pos %in% linked_positions, c("haplotype", "pos", "alt")]
  leaky <- unique(ref$haplotype[paste(ref$pos, ref$alt) %in% paste(baked$pos, baked$alt)])
  min(1, sum(freqs[leaky]))
}

#' Haplotype-group pairs whose ALT primers cross-prime each other
#'
#' Every haplotype ALT primer shares the target 3' base (all ALT haplotypes carry
#' the ALT there), so what tells the groups apart is the primer *body*.  Two groups
#' whose designs differ only outside the 3'-proximal window have primers that
#' anneal identically at their 3' ends (Taq extends from a matched 3' end): each
#' fires on the other's template, so the pooled tube cannot resolve the groups and
#' they compete for one site.  This counts the unordered pairs of groups that share
#' that window.  [primer_set_compatibility()] cannot see it - the primers are
#' near-identical, not complementary.
#' @param haplotypes a tidy haplotype table.
#' @param snp one-row target.
#' @param designs the kept designs (see [haplotype_cover()]), each with `rep` (a
#'   haplotype name) and `design` (the assembler's `marker_assay`).
#' @return an integer count of collapsible group pairs.
#' @keywords internal
#' @noRd
lax_cross_priming_pairs <- function(haplotypes, snp, designs) {
  if (length(designs) < 2L) {
    return(0L)
  }
  h <- normalise_haplotypes(haplotypes)
  # every design sits on the base strand (design_lax_marker forces `strands`), so
  # the AS primers share one 3'-proximal interval: one window signature per group
  a <- designs[[1L]]$design$assay[1L, ]
  win <- vapply(designs, function(d) {
    hap_signature(h, d$rep, snp, a$win3p_start, a$win3p_end)
  }, character(1))
  counts <- as.numeric(table(win))
  as.integer(sum(counts * (counts - 1) / 2))
}

#' The linked positions the kept designs bake into their ALT primers
#'
#' The designer hands each group's representative to `assemble_kasp_assay()`
#' (`linked =`), so the positions that actually rewrite that group's ALT primer are
#' its linked set intersected with the primer's own footprint.
#' @param linked the haplotype engine's linked table (see [haplotype_variants()]);
#'   `NULL` when there is none.
#' @param designs the kept designs (see [haplotype_cover()]), each with `rep` and
#'   `fp` (its AS footprint).
#' @return a unique integer vector of positions (`integer(0)` when none).
#' @keywords internal
#' @noRd
lax_baked_positions <- function(linked, designs) {
  if (is.null(linked) || !nrow(linked) || !length(designs)) {
    return(integer(0))
  }
  unique(unlist(lapply(designs, function(d) {
    lp <- linked$linked[linked$haplotype == d$rep]
    lp[lp >= d$fp[1L] & lp <= d$fp[2L]]
  })))
}

#' The population frequency of ALT haplotypes the tube can resolve
#'
#' The ranking counterpart of [lax_cross_priming_pairs()]: a group whose
#' 3'-proximal window is unique among the designs is separable (its ALT primer
#' fires only on its own templates), while groups that share a window cannot be
#' told apart.  This weights that by the population each covers - the scheme's
#' actual quality, not a pass/fail.
#' @param haplotypes a tidy haplotype table.
#' @param snp one-row target.
#' @param cover the [haplotype_cover()] result (`groups` + `designs`).
#' @return a fraction in `[0, 1]`.
#' @keywords internal
#' @noRd
lax_resolvable_freq <- function(haplotypes, snp, cover) {
  designs <- cover$designs
  if (!length(designs)) {
    return(0)
  }
  h <- normalise_haplotypes(haplotypes)
  a <- designs[[1L]]$design$assay[1L, ]
  win <- vapply(designs, function(d) {
    hap_signature(h, d$rep, snp, a$win3p_start, a$win3p_end)
  }, character(1))
  uniq <- !duplicated(win) & !duplicated(win, fromLast = TRUE)
  # the cover table labels each design's group by its representative's footprint
  # signature, so a design's templates are the rows sharing its label
  grp <- vapply(designs, function(d) {
    g <- cover$groups$group[cover$groups$haplotype == d$rep]
    if (length(g)) g[1L] else NA_character_
  }, character(1))
  haps <- cover$groups$haplotype[cover$groups$group %in% grp[uniq]]
  freqs <- h |>
    dplyr::distinct(.data$haplotype, .data$freq) |>
    tibble::deframe()
  min(1, sum(freqs[haps]))
}

#' Run one tube through the lax decision stage
#'
#' The set-level twin of [select_with_ladder()]: the per-haplotype primers are
#' already chosen, so the table has exactly one row (the tube) and the shared
#' ladder only has to decide whether its three gate margins clear.  Running this
#' here is what makes the declared lax gates real rather than decorative (C7): the
#' designer records the outcome and the observed metrics on the marker.
#' @param n_uncovered ALT haplotypes no kept design covers.
#' @param leakage wild-type leakage fraction (see [lax_ref_leakage()]).
#' @param n_cross collapsible group pairs (see [lax_cross_priming_pairs()]).
#' @param resolvable the tube's ranking margin (see [lax_resolvable_freq()]).
#' @return `list(pick, metrics)`: the ladder's outcome (`best`, `relaxed`,
#'   `violated`) and the observed numbers.
#' @keywords internal
#' @noRd
lax_stage_pick <- function(n_uncovered = 0L, leakage = 0, n_cross = 0L, resolvable = 0) {
  data <- tibble::tibble(
    m_coverage = -as.numeric(n_uncovered),
    m_alt_specific = KASP_LAX_ALT_SPECIFIC_TOL - leakage,
    m_cross_priming = -as.numeric(n_cross),
    tier1 = resolvable
  )
  pick <- select_with_ladder(
    data,
    gates = KASP_GATES_LAX, relax_order = KASP_RELAX_LAX,
    tiers = names(KASP_TIERS_LAX), ties = character(0),
    gate_margin = c(
      gate_lax_coverage = "m_coverage",
      gate_lax_alt_specific = "m_alt_specific",
      gate_lax_cross_priming = "m_cross_priming"
    ),
    gate_floor = KASP_GATE_FLOOR_LAX
  )
  warn_relaxed(pick$violated)
  list(
    pick = pick,
    metrics = list(
      coverage = n_uncovered, alt_specific = leakage,
      cross_priming = n_cross, resolvable = resolvable
    )
  )
}

#' The lax stage's record on the marker (`$selection`)
#' @param lax a [lax_stage_pick()] result.
#' @param designs the kept designs, whose own `as`/`pair` outcomes are rolled up.
#' @return a list (`relaxed`/`violated`/`doomed`/the metrics + the rolled-up
#'   `as`/`pair` gate ids) suitable as a marker's `$selection`.
#' @keywords internal
#' @noRd
lax_selection <- function(lax, designs = list()) {
  roll <- function(field) {
    unique(unlist(lapply(designs, function(d) d$design$selection[[field]] %||% character(0))))
  }
  list(
    # the lax (set-level) stage's ladder outcome, and the metrics behind it
    relaxed = lax$pick$relaxed,
    violated = lax$pick$violated,
    doomed = is.null(lax$pick$best),
    coverage = lax$metrics$coverage,
    alt_specific = lax$metrics$alt_specific,
    cross_priming = lax$metrics$cross_priming,
    resolvable = lax$metrics$resolvable,
    # the per-haplotype `as`/`pair` stages, rolled up over the kept designs (the
    # marker is a tube, not one assay, so these are the union of its groups')
    as = roll("as"),
    pair = roll("violated")
  )
}

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
#'
#' The tube is then put through the scheme's own *set-level* stage, so the
#' declared lax gates are evaluated and recorded in `$selection`: every ALT
#' haplotype must be covered (hard), the baked linked alleles must be ALT-specific
#' ([lax_ref_leakage()]) and the haplotype ALT primers must not cross-prime each
#' other ([lax_cross_priming_pairs()]).  The last two are soft, with bounded floors
#' (see `KASP_GATE_FLOOR_LAX`), so a relaxed tube comes back with a warning and a
#' `"warn"` verdict.
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param snp one-row target (`chrom`, `pos`, `ref`, `alt`).
#' @param haplotypes tidy haplotype table (`haplotype`, `freq`, `chrom`, `pos`,
#'   `ref`, `alt`); the ALT haplotypes define the primers.  `NULL` degenerates to
#'   a single pure-KASP primer (no linked variants).
#' @param ... passed to [assemble_kasp_assay()].
#' @return a `marker_lax` list: `snp`, `haplotypes`, `alt` (the ALT haplotype
#'   names), `base` (the shared backbone - one ordinary KASP assembly whose REF
#'   primer and COM are reused), `alt_primers` (one ALT primer per haplotype, on
#'   the base's strand; `NA` when a haplotype did not design), `pool` +
#'   `compatibility` (the *joint* check of the shared REF + shared COM + every ALT
#'   primer, see [primer_set_compatibility()]), `selection` (the lax stage's gate
#'   outcome and the metrics behind it, see [lax_stage_pick()]) and `verdict`
#'   (`"fail"` when no haplotype designs or the hard coverage gate blocks,
#'   `"warn"` when a soft gate was relaxed or the pool is not clean, else
#'   `"pass"`).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' design_lax_marker(genome, snp, haplotypes)
#' }
#' @family KASP design
#' @export
design_lax_marker <- function(genome, snp, haplotypes = NULL, ...) {
  if (is.null(haplotypes)) {
    # no haplotype table: a single, ordinary ALT primer (pure-KASP behaviour).
    # There is no haplotype set to cover, no linked allele to bake and no second
    # ALT primer to compete with, so the lax stage's risks are vacuous - run it
    # anyway, so every declared gate is evaluated on every path (C7)
    a <- suppressWarnings(assemble_kasp_assay(genome, snp, ...))
    lax <- lax_stage_pick()
    out <- list(
      snp = snp, haplotypes = NULL, alt = "ALT", assays = list(ALT = a),
      designed = if (is_target_designable(a)) "ALT" else character(0),
      selection = lax_selection(lax),
      verdict = if (is_target_designable(a) && !is.null(lax$pick$best)) "pass" else "fail"
    )
    return(structure(out, class = "marker_lax"))
  }
  hv <- haplotype_variants(haplotypes, snp)
  alt_hap <- hv$alt
  if (!length(alt_hap)) {
    cli::cli_inform(gettext("No ALT haplotype carries the target allele; nothing to design."))
    return(NULL)
  }
  # one shared backbone for the whole tube: the REF primer and the COM come from a
  # single ordinary KASP assembly, and only the ALT primer varies per haplotype.
  # Forcing every ALT primer onto the base's strand keeps them compatible with the
  # shared COM (they are the forward/reverse partner of the same reverse primer).
  base <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = hv$variants, ...))
  if (!is_target_designable(base)) {
    # the backbone is doomed, so no group can be covered: say so through the same
    # stage the rest of the marker uses, rather than a bare "fail"
    lax <- lax_stage_pick(n_uncovered = length(alt_hap))
    return(structure(list(
      snp = snp, haplotypes = hv$haplotypes, alt = alt_hap, base = base,
      alt_primers = NULL, pool = NULL, compatibility = NULL,
      selection = lax_selection(lax), verdict = "fail"
    ), class = "marker_lax"))
  }
  strand <- base$strand
  h <- hv$haplotypes
  # the shared two-pass engine: one ALT primer per primer-window group, merged by
  # the footprint each primer actually covers (see haplotype_cover())
  cover <- haplotype_cover(h, snp, function(hpk) {
    lk <- if (!is.null(hv$linked)) {
      d <- hv$linked[hv$linked$haplotype == hpk, c("target", "linked"), drop = FALSE]
      if (nrow(d)) d else NULL
    } else {
      NULL
    }
    suppressWarnings(assemble_kasp_assay(
      genome, snp,
      variants = hv$variants, linked = lk, strands = strand, ...
    ))
  })
  designs <- cover$designs
  alt_primers <- stats::setNames(
    vapply(cover$groups$design, function(d) {
      if (is.null(d)) NA_character_ else d$assay$as_alt[1L]
    }, character(1)),
    cover$groups$haplotype
  )
  ok <- !is.na(alt_primers)
  # ---- the lax stage: the scheme's own, set-level decision (C7) --------------
  n_uncovered <- length(cover$uncovered)
  leakage <- lax_ref_leakage(h, snp, lax_baked_positions(hv$linked, designs))
  n_cross <- lax_cross_priming_pairs(h, snp, designs)
  lax <- lax_stage_pick(n_uncovered, leakage, n_cross, lax_resolvable_freq(h, snp, cover))
  doomed <- is.null(lax$pick$best)
  if (doomed) {
    reasons <- c(
      if (n_uncovered > 0L) {
        gettext("no ALT primer covers {n_uncovered} ALT haplotype(s)")
      },
      if (leakage > 2 * KASP_LAX_ALT_SPECIFIC_TOL) {
        gettext("the wild type also carries a baked linked allele (leakage {round(leakage, 2)})")
      },
      if (n_cross > KASP_LAX_CROSS_PRIMING_FLOOR) {
        gettext("{n_cross} haplotype-group pairs cross-prime")
      }
    )
    cli::cli_warn(c(
      gettext("The lax tube fails the haplotype stage."),
      rlang::set_names(reasons, rep("*", length(reasons)))
    ))
  }
  # the *joint* pool: the shared REF primer + COM + one ALT primer per group
  prim <- c(
    AS_ref = base$assay$as_ref[1L], COM = base$assay$com[1L],
    stats::setNames(alt_primers[ok], paste0("AS_alt_", alt_hap[ok]))
  )
  prim <- prim[!duplicated(prim)]
  compat <- if (length(prim) > 1L) primer_set_compatibility(prim) else NULL
  pool_ok <- is.null(compat) || all(compat$pass)
  verdict <- if (doomed || !any(ok)) {
    "fail"
  } else if (pool_ok && !length(lax$pick$violated)) {
    "pass"
  } else {
    "warn"
  }
  out <- list(
    snp = snp, haplotypes = h, alt = alt_hap, base = base,
    haplotype_groups = tibble::tibble(group = cover$groups$group, haplotype = cover$groups$haplotype),
    alt_primers = alt_primers, pool = prim, compatibility = compat,
    selection = lax_selection(lax, designs),
    verdict = verdict
  )
  structure(out, class = "marker_lax")
}

#' @export
print.marker_lax <- function(x, ...) {
  cat(
    "<marker_lax>", x$snp$chrom, ":", x$snp$pos, " ", x$snp$ref, ">", x$snp$alt,
    " | ", length(x$alt), " ALT haplotype(s) | verdict ", x$verdict, "\n",
    sep = ""
  )
  if (is.null(x$alt_primers)) {
    cat("  plain ALT primer (no haplotype table given)\n")
  } else {
    for (hp in x$alt) {
      p <- x$alt_primers[[hp]]
      cat("  haplotype ", hp, " ALT primer ",
        if (is.null(p) || is.na(p)) "not designed" else p,
        "\n",
        sep = ""
      )
    }
  }
  if (!is.null(x$pool)) {
    cat("  shared: REF ", x$pool[["AS_ref"]], " COM ", x$pool[["COM"]], "\n", sep = "")
  }
  if (!is.null(x$compatibility)) {
    bad <- x$compatibility[!x$compatibility$pass, , drop = FALSE]
    cat("  pool (", nrow(x$compatibility), " pairs): ",
      if (nrow(bad)) paste("failed", paste(paste(bad$a, bad$b, sep = "-"), collapse = ", ")) else "clean",
      "\n",
      sep = ""
    )
  }
  if (!is.null(x$selection)) {
    cat(
      "  lax stage: ", x$selection$coverage, " uncovered, wild-type leakage ",
      round(x$selection$alt_specific, 3L), ", cross-priming pairs ",
      x$selection$cross_priming, "\n",
      sep = ""
    )
    if (length(x$selection$violated)) {
      cat("    relaxed: ", paste(x$selection$violated, collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}

#' The lax-KASP scheme spec (discovered by the shared registry, see [marker_schemes()])
#' @keywords internal
#' @noRd
scheme_def_lax <- function() {
  list(
    scheme = "lax", family = "kasp",
    label = "lax KASP (one primer per haplotype)",
    discrimination = "haplotype block (multiple linked variants), not one 3' base",
    readout = "co-dominant two-colour", zygosity = "co-dominant",
    amplicon = "critical", loci = 1L, verdict = "usable lax KASP assay",
    gates = KASP_GATES_LAX, relax = KASP_RELAX_LAX,
    # one tube is one candidate, so the single tier reports the tube's resolvable
    # haplotype frequency rather than breaking a tie between candidates (C14)
    tiers = KASP_TIERS_LAX,
    gate_short = KASP_GATE_SHORT_LAX,
    gate_limit = KASP_GATE_LIMIT_LAX,
    gate_floor_label = KASP_GATE_FLOOR_LABEL_LAX,
    gate_hard_label = KASP_GATE_HARD_LABEL_LAX,
    gate_label = KASP_GATE_LABEL_LAX,
    stage_label = "lax haplotype set"
  )
}
