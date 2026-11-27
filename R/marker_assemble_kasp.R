# =============================================================================
#  marker_assemble_kasp.R  --  assemble_kasp_assay(): AS x COM x contract
# =============================================================================
#  Split from marker_assemble.R (roadmap ZT T7); the AS candidate stages stay in
#  marker_assemble.R, the fixed-primer constructor in marker_assemble_fixed.R.
# =============================================================================
#' Assemble and score the best KASP assay for a target SNP
#'
#' When `variants` is supplied the selection is *variant-aware*: among primers
#' that pass the physico-chemical filters we minimise the number of other
#' variants falling under the primer footprint (the 3'-terminal 5 nt are handled
#' earlier by [pick_target_snp()]), then maximise discrimination.
#'
#' @section Selection policy:
#' The AS placement is chosen lexicographically by:
#' 1. cleanest 3' window (`n_3p`, best of the two orientations);
#' 2. fewest other variants under the primer footprint (`n_var`);
#' 3. highest discrimination (DECIPHER);
#' 4. not a weak G:T / A:C wobble mismatch;
#' 5. quieter 3'-seed genome-wide (`seed_hits`, best of the two strands) -- a
#'    specificity tie-break computed only for the top `specificity_candidates`
#'    placements, mirroring [pick_target_snp()].  In the returned
#'    `as_candidates` table, `seed_hits` is `NA` for every placement outside that
#'    shortlist: `NA` means *not screened*, not *not specific*, so it is a
#'    partial tie-break and not a complete ranking key;
#' 6. then `|Tm - target|`, 3'-stability and length.
#' The common primer is then placed on the winning strand and ranked by
#' amplicon length, `COM Tm >= AS Tm`, dimers/hairpin, footprint load and Tm
#' proximity (see [design_common_primers()]).
#' @section KASP design pipeline:
#' Typical use is `ath_genome()` -> `ensembl_gene()`/`ensembl_variants()` ->
#' `pick_target_snp()` -> `assemble_kasp_assay()` -> `marker_design_checks()`.
#' See `vignette("kasp-intro")`.
#'
#' @param genome a `marker_genome` or `DNAStringSet` for any genome; the examples
#'   use the bundled Arabidopsis TAIR10 demo ([ath_genome()]).
#' @param snp    one-row tibble (chrom, pos, ref, alt) from [pick_target_snp()].
#' @param variants variant table (for the same chromosome) to count overlaps.
#' @param strands which AS-primer orientations to search (see [design_as_primers()]).
#' @param tm_range accept AS primers whose Tm falls in this window (C).
#' @param clean_3p size of the "clean" AS 3' window used to count other variants.
#' @param max_footprint_variants hard cap on *other* variants under the AS primer
#'   footprint ("N-masking" as a constraint): `0` allows none, the default `Inf`
#'   only penalises them in the ranking.  Use [mask_positions()] to blank such
#'   positions in your own template if you need the literal `N`s.
#' @param min_disc minimum discrimination (on- minus off-target efficiency).
#' @param specificity_candidates how many leading AS placements to screen for
#'   3'-seed genome-wide uniqueness when breaking ties (mirrors
#'   [pick_target_snp()]; `0` disables the screen).
#' @param pair_screen how many leading pairings (by the genome-free ranking) to
#'   genome-screen for off-target exact primer hits, a second in-silico PCR
#'   product and an inverted repeat; the first clean one wins.  These need a
#'   genome/PCR pass, so only the leaders are screened (`0` disables; the ranked
#'   winner is then kept).
#' @param on_offtarget the specificity policy (roadmap Y3): `"doom"` - a pairing
#'   that fails the genome-level screen is not returned, the target is doomed with
#'   the reason; `"warn"` - keep the warned ranked leader.  `NULL` (default)
#'   follows the scheme's own policy (`marker_specificity()`): `unique` /
#'   `homeolog-aware` doom, the tolerant default `multicopy` / `off` warn.
#' @param run_blast also run a whole-primer genome scan (needs NCBI BLAST+, or
#'   a `bowtie_index` on the genome).  On the indexed path the result is an exact
#'   whole-primer scan (`scan` attribute `"whole-primer"`) of the AS primers
#'   *with* their deliberate mismatch; `specificity` is the seed tier and is not
#'   the same table.
#' @param cores processors for the per-SNP subprocess work (the DECIPHER PCR
#'   efficiencies and the OligoArrayAux dG maps); `1` (the default, from
#'   `getOption("tidymarkers.cores")`) keeps CRAN happy, and larger values give a
#'   near-linear speed-up of a single assay.  When `BPPARAM` is set it wins.
#' @param BPPARAM optional `BiocParallel` backend for the genome-wide steps
#'   (seed specificity and in-silico PCR; split by chromosome) and for the
#'   subprocess work.  Use it to avoid nesting: inside a locus-level backend pass
#'   `cores = 1` / a serial `BPPARAM`.
#' @param diagnostic_positions the subgenome-diagnostic mask (roadmap Y3.1):
#'   `NULL`, a data frame / list with `chrom` and `pos`, a `GRanges`, or a named
#'   list `chrom -> integer positions`.  It is the PolyMarker mask derived by
#'   comparing subgenomes, so the candidate never infers it from one genome; with
#'   `NULL` the subgenome verdict is `NA`, with a supplied mask the common
#'   primer's 3' base is checked against it and a miss warns and reports
#'   `multicopy`.
#' @return a `marker_assay` object: a list holding the `snp`, the template, the
#'   chosen `strand`, ranked `assay`/`as_candidates` tibbles, `efficiency`,
#'   `thermo`, `complementarity`, `specificity`, `amplicons`, `blast` results and
#'   the tailed `tails` primers (plus `haplotypes`/`haplotype_hits` when a
#'   haplotype table was given), and `evaluation` (the [marker_design_checks()]
#'   table for this assay, stored once at assembly time).  When no candidate clears the gates it returns a
#'   `marker_doomed` object instead (the target, template, `reason`, `stage` and the
#'   candidate table), so a batch or a plot can report the target rather than
#'   abort.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' v <- ensembl_variants("Chr4", 268901, 271503)
#' snp <- pick_target_snp(v, genome, ensembl_gene("AT4G00650"))
#' assemble_kasp_assay(genome, snp)
#' }
#' @section External tool:
#' The thermodynamics engine (DECIPHER via OligoArrayAux) needs `hybrid-min` on
#' the `PATH`; if the data tables are not found automatically, point them out:
#' `with_unafloddat("/path/to/share/oligoarrayaux", assemble_kasp_assay(...))`.
#' @param flank template half-width (bp).
#' @param linked optional two-column table (`target`, `linked` positions) of
#'   co-segregating substitutions; those riding with `snp$pos` are baked into the
#'   ALT primer so it mismatches REF twice (see [design_as_primers()]).
#' @param haplotypes optional tidy haplotype table
#'   (`haplotype`, `freq`, `chrom`, `pos`, `ref`, `alt`).  Its per-position
#'   frequencies become the `pop_freq` weights and its co-segregating positions
#'   the `linked` set (both derived, see `haplotype_variants()`); the object also
#'   carries `haplotypes` and a per-haplotype `haplotype_hits` table (which AS
#'   primer each is expected to bind).
#' @param min_gc,max_gc accepted AS GC window.
#' @param limits acceptance thresholds, a named list like `KASP_LIMITS`.  The
#'   same object grades the candidates and is stored on the object (see
#'   [marker_design_checks()]); [design_asq_marker()] passes ASQ's own windows.
#' @param variant_offset bases between the variant block and the primer's 3'
#'   terminus: `0` puts the variant at the terminal base (KASP), `1` one base
#'   inside it (ASQ's penultimate placement, see [design_asq_marker()]).
#' @param as_tm_target preferred AS-primer Tm (C) the lexicographic grader
#'   ranks towards; the accepted window is still `limits$as_tm` / `tm_range`.
#' @param com_tm_target preferred common-primer Tm (C) the pair grader ranks
#'   towards (default `COM_TM_TARGET`).
#' @param com_tm_range optional common-primer Tm window; `NULL` reuses
#'   `tm_range`, so both primers are held to the same window.
#' @param com_gaps distances (bp) from the variant anchor to the common primer
#'   that the search tries; passed to [design_common_primers()]'s `gaps`.
#' @param as_lengths AS-primer lengths to try (passed to
#'   [design_as_primers()]'s `lengths`; ASQ needs longer primers than KASP's
#'   20-27 nt to reach its Tm window).
#' @param scheme the marker chemistry this assay is assembled *as* (see
#'   [marker_schemes()]): stored on the object, and used to build the matching
#'   decision graph and to grade the acceptance checks, so a scheme never shows
#'   (or is graded by) another chemistry's tree (C12).
#' @param engines optional `marker_engines` adapter list (see `engine_has()`);
#'   `NULL` reads the session default.  It gates the degenerate-primer efficiency
#'   screen (and the fixed-primer `assemble_primer_assay()` it delegates to).
#' @family design
#' @export
assemble_kasp_assay <- function(genome, snp, variants = NULL, flank = 150L,
                                strands = c("+", "-"),
                                min_gc = 0.35, max_gc = 0.70,
                                tm_range = KASP_LIMITS$as_tm, clean_3p = 5L,
                                max_footprint_variants = Inf,
                                min_disc = 0.5, specificity_candidates = 10L,
                                pair_screen = 3L,
                                linked = NULL, haplotypes = NULL, run_blast = FALSE,
                                scheme = "snp", on_offtarget = NULL,
                                cores = getOption("tidymarkers.cores", 1L), BPPARAM = NULL,
                                diagnostic_positions = NULL,
                                limits = KASP_LIMITS, variant_offset = 0L,
                                as_tm_target = AS_TM_TARGET,
                                com_tm_target = COM_TM_TARGET,
                                com_tm_range = NULL, com_gaps = 18:80,
                                as_lengths = KASP_AS_LENGTHS, engines = NULL) {
  call <- sys.call()
  if (contracts_enabled()) contract_assay_input(genome, snp, call)
  if (is.null(BPPARAM)) BPPARAM <- marker_bpparam(cores)
  mask <- normalise_diagnostic_positions(diagnostic_positions)
  # the scheme's own reaction conditions, so every dG this assay reports is
  # evaluated at the temperature the chemistry actually runs at (Z3)
  conditions <- scheme_conditions(scheme)
  # Y3: the specificity policy.  A `unique` chemistry never returns a marker whose
  # leading pairing fails the genome-level screen - it is doomed with the reason
  # (no silent warn-and-return); `multicopy`/`off`/`homeolog-aware` keep the
  # warned leader.  `on_offtarget` overrides the scheme's own policy.
  on_offtarget <- on_offtarget %||%
    if (marker_specificity(scheme) %in% KASP_SPECIFICITY_STRICT) "doom" else "warn"
  on_offtarget <- match.arg(on_offtarget, c("warn", "doom"))
  # one-line install UX: a missing tool is a hint (once per session), not a silent NA
  if (!has_oligoarrayaux() && !isTRUE(.marker_state$tool_hint_shown)) {
    .marker_state$tool_hint_shown <- TRUE
    cli::cli_inform(c(
      gettext("OligoArrayAux was not found: thermodynamic and structural checks will be NA."),
      i = gettext("run tidymarkers::marker_install_tools() once to build it.")
    ))
  }
  # the target's allele(s): a degenerate all-SNP site designs one consensus ALT primer
  alts <- marker_snp_alts(snp)
  alts_in <- alts$alts_in
  degenerate <- alts$degenerate
  alt_pass <- alts$alt_pass
  snp <- alts$snp
  tmpl <- snp_template(genome, snp$chrom, snp$pos, flank)
  hap <- NULL
  if (!is.null(haplotypes)) {
    hap <- haplotype_variants(haplotypes, snp)
    if (is.null(variants)) variants <- hap$variants
    if (is.null(linked)) linked <- hap$linked
  }
  off <- snp$pos - flank # template coordinate -> genomic coordinate
  # the variant table's alleles/weights, and the mapping to template coordinates
  vc <- marker_variant_alleles(variants, snp, tmpl$seq, off)
  vpos <- vc$vpos
  vweight <- vc$vweight
  vsub <- vc$vsub
  valt <- vc$valt
  ref_len <- nchar(snp$ref)
  alt_len <- nchar(snp$alt)
  is_indel <- ref_len != alt_len
  vblock <- seq.int(snp$pos, snp$pos + max(ref_len, alt_len) - 1L) # the target's own bases

  # one service for both critical zones (AS 3' window, COM footprint): DECIPHER
  # turns a footprint variant into an expected dropout (NA without it / on error)
  use_eff <- engine_has("pcr_efficiency", engines)
  batch_or_safe <- marker_dropout_batch(tmpl, off, vpos, valt, vweight, use_eff, BPPARAM)

  # co-segregating substitutions to bake into the ALT primer (see design_as_primers)
  lo <- marker_linked_offsets(linked, snp, vpos, valt, vsub, tmpl$seq, off, tmpl$snp_i, vblock)
  linked_off <- lo$linked_off
  linked_skipped <- lo$linked_skipped
  as_exclude <- lo$as_exclude

  # all primer geometry comes from the shared helpers (plus-strand coordinates)
  com_iv <- function(strand, L, gap) {
    primer_interval(as_anchor(snp$pos, ref_len, strand), strand, L, gap = gap, role = "COM")
  }

  as_cand <- marker_as_candidates(
    tmpl, snp, alt_pass, is_indel, off, strands, tm_range, min_gc, max_gc,
    linked_off, ref_len, alt_len, BPPARAM,
    conditions = conditions, variant_offset = variant_offset, as_lengths = as_lengths
  )
  if (is_kasp_doomed(as_cand)) {
    return(as_cand)
  }

  as_cand <- marker_as_zone_metrics(
    as_cand, clean_3p, vpos, vweight, as_exclude, batch_or_safe,
    max_footprint_variants, snp, tmpl,
    variant_offset = variant_offset
  )
  if (is_kasp_doomed(as_cand)) {
    return(as_cand)
  }

  # prioritise with the shared grader: gate the grave sins, then rank
  # lexicographically on the margins (the genome 3'-seed screen is a tie-break)
  as_pick <- marker_pick_as(
    as_cand, genome, min_disc, specificity_candidates, BPPARAM, snp, tmpl,
    limits = limits, tm_target = as_tm_target
  )
  if (is_kasp_doomed(as_pick)) {
    return(as_pick)
  }
  best_as <- as_pick$best_as
  pick_as <- as_pick$pick_as
  # the ranked candidate table for inspection, minus the grader's helper columns
  as_cand <- as_pick$as_cand

  # name the soft failures and the linkage dependence (see marker_as_warnings())
  linked_dependent <- marker_as_warnings(best_as, snp, clean_3p, min_disc, lo$link_pos)

  # COM: design on the SAME strand as the winning AS primer so the amplicon
  # closes.  design_common_primers() works in TEMPLATE coordinates, so its
  # anchor is the template-local variant anchor (not the genomic one).
  common <- design_common_primers(tmpl$seq, as_anchor(tmpl$snp_i, ref_len, best_as$strand),
    strand = best_as$strand, gaps = com_gaps,
    tm_range = com_tm_range %||% tm_range, n = 60L,
    conditions = conditions
  ) |>
    dplyr::mutate(
      com_homodimer = if (engine_has("duplex_dg", engines)) {
        oligo_dg_vec(.data$com, .data$com, BPPARAM = BPPARAM)
      } else {
        NA_real_
      },
      com_hairpin = if (engine_has("hairpin_dg", engines)) {
        oligo_dg_vec(.data$com, BPPARAM = BPPARAM)
      } else {
        NA_real_
      }
    )

  as_tm_mean <- mean(best_as$tm)

  assays <- common |>
    dplyr::mutate(
      strand = best_as$strand,
      com_start = purrr::pmap_int(
        list(.data$strand, .data$len, .data$gap),
        function(st, L, g) IRanges::start(com_iv(st, L, g))
      ),
      com_end = purrr::pmap_int(
        list(.data$strand, .data$len, .data$gap),
        function(st, L, g) IRanges::end(com_iv(st, L, g))
      ),
      n_var = purrr::map2_int(.data$gap, .data$len, function(g, L) {
        count_variants_in(com_iv(best_as$strand, L, g), vpos, vblock)
      }),
      w_var = purrr::map2_dbl(.data$gap, .data$len, function(g, L) {
        variant_load_in(com_iv(best_as$strand, L, g), vpos, vweight, vblock)
      }),
      # expected dropout of the COM footprint (position- and identity-aware)
      com_dropout = batch_or_safe(
        .data$com, .data$com_start, .data$com_end,
        .data$strand != "+", .data$n_var, vblock
      ),
      # the COM primer's 3'-terminal dG (README section 3: a low 3'-end
      # stability), ranked by the pair stage's "COM 3' end weakness" tier
      com_dg_3p = terminal_dg(.data$com, BPPARAM = BPPARAM),
      as_ref = best_as$ref_primer, as_alt = best_as$alt_primer,
      as_alt_degen = best_as$alt_degen, as_dropout = best_as$as_dropout,
      as_foot_dropout = best_as$as_foot_dropout,
      discrimination = best_as$discrimination, as_gc = best_as$gc,
      as_len = best_as$as_length, as_n_var = best_as$n_var, as_n_3p = best_as$n_3p,
      # worst-case product over both alleles (indels differ by delta)
      amplicon = pmax(
        pmax(best_as$as_end, .data$com_end) - pmin(best_as$as_start, .data$com_start) + 1L,
        pmax(best_as$as_alt_end, .data$com_end) - pmin(best_as$as_alt_start, .data$com_start) + 1L
      ),
      # store the geometry once, for the plots to read instead of recomputing
      as_start = best_as$as_start, as_end = best_as$as_end,
      as_alt_start = best_as$as_alt_start, as_alt_end = best_as$as_alt_end,
      win3p_start = best_as$win3p_start, win3p_end = best_as$win3p_end
    )

  # choose among the pairings: gate the grave sins, then rank lexicographically
  # on tier margins, relaxing gates one at a time only when nothing clears them
  graded <- grade_assay_pairs(assays, limits, as_tm_mean, com_tm_target = com_tm_target)
  # name the gate(s) that could not be evaluated at all (all margins NA), so a
  # missing engine is reported as unknown rather than as a plain gate failure
  unknown <- names(graded$gate_margin)[vapply(graded$gate_margin, function(col) {
    m <- graded$data[[col]]
    length(m) > 0L && all(is.na(m))
  }, logical(1))]
  if (length(unknown)) {
    cli::cli_inform(c(
      gettext("The assay pairing could not be graded: {paste(unknown, collapse = ', ')} has no value."),
      i = gettext("A gate whose margin is NA is unknown, not a pass; OligoArrayAux supplies the structure dG.")
    ))
  }
  pick <- select_with_ladder(
    graded$data, graded$gates, KASP_RELAX_PAIR, graded$tiers, graded$ties,
    graded$gate_margin, graded$gate_floor
  )
  if (is.null(pick$best)) {
    # the doom must carry the reason: an unevaluable gate is not a plain failure
    reason <- if (length(unknown)) {
      glue::glue(gettext("No candidate assay could be graded: {paste(unknown, collapse = ', ')} has no value."))
    } else {
      gettext("No candidate assay survived even the relaxed constraints.")
    }
    return(marker_doomed(reason, "pair", snp, tmpl, as_cand))
  }
  best <- pick$best
  # keep the chosen pairing first: the rest of the package (checks, plots) reads
  # assay[1, ] as "the" assay, so the ranking must be reflected in the row order
  ord <- c(pick$best$.row, setdiff(seq_len(nrow(assays)), pick$best$.row))
  assays <- assays[ord, , drop = FALSE]
  assays$rank <- seq_len(nrow(assays))
  warn_relaxed(pick$violated)

  # The ladder ranks on genome-free margins, so a second in-silico PCR product, an
  # off-target exact primer hit or an inverted repeat can still sit on the winner.
  # Those cost a genome/PCR pass, so only the leading `pair_screen` pairings are
  # screened; the first clean one wins (the ranked winner is kept when none is).
  if (!is.null(genome) && pair_screen > 1L && nrow(assays) > 1L) {
    n_screen <- min(as.integer(pair_screen), nrow(assays))
    # The screen costs one genome pass per pairing (its COM primer, and the two
    # alleles of its virtual PCR); `pair_screen` pairings therefore cost
    # `pair_screen` passes of the same radius.  Each pass is per primer, not per
    # pattern, so the whole screen is scanned once up front instead: every COM
    # primer, every AVS pattern of the screened pairings (`as_hits` included,
    # which is scanned below anyway) and the reverse complement of each (a primer
    # binds either strand) is searched in one grouped pass, and `hit_count()`
    # reads that scan instead of rescanning; the `in_silico_pcr()` calls it
    # replaces are not re-issued.
    as_primers <- c(AS_ref = assays$as_ref[1L], AS_alt = assays$as_alt[1L])
    screened <- assays[seq_len(n_screen), , drop = FALSE]
    rc_pattern <- function(p) {
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(p)))
    }
    screen_template <- function(i, side) {
      row <- screened[i, ]
      # the same pairs `as_pcr_pairs()` builds, so the scan covers what the PCR asks
      fr <- as_pcr_pairs(list(assay = row, template = tmpl, snp = snp))
      vapply(fr, function(x) x[[side]], character(1))
    }
    screen_patterns <- c(
      as_primers,
      vapply(as_primers, rc_pattern, character(1)),
      screened$com,
      vapply(screened$com, rc_pattern, character(1)),
      unlist(lapply(seq_len(n_screen), function(i) screen_template(i, "F")), use.names = FALSE),
      # the reverse primers are searched as their reverse complements
      unlist(lapply(seq_len(n_screen), function(i) screen_template(i, "R")), use.names = FALSE) |>
        lapply(rc_pattern) |>
        unlist(use.names = FALSE)
    )
    screen_hits <- tryCatch(
      pattern_hits_batch(screen_patterns, genome, 0L),
      error = function(e) NULL
    )
    hit_count <- function(primers) {
      # One row per (primer, chromosome) that holds an exact hit on either strand:
      # the shape `primer_hits()` returns and the only thing the guards below read
      # (`table(hits$primer)`).  A table without the `primer` column made that
      # count NULL at every guard, so the screen kept an off-target pairing.
      if (is.null(screen_hits)) {
        return(tryCatch(primer_hits(primers, genome, BPPARAM = BPPARAM), error = function(e) NULL))
      }
      purrr::imap(primers, function(p, nm) {
        back <- rc_pattern(p)
        purrr::map(names(screen_hits[[p]]), function(ch) {
          n <- length(screen_hits[[p]][[ch]]) + length(screen_hits[[back]][[ch]])
          if (!n) NULL else tibble::tibble(primer = nm, chrom = ch, n_hits = n)
        }) |>
          purrr::list_rbind()
      }) |>
        purrr::list_rbind()
    }
    as_hits <- hit_count(as_primers)
    as_bad <- !is.null(as_hits) && nrow(as_hits) && max(table(as_hits$primer)) > 1L
    # the virtual PCR of every screened pairing is one scan of their pooled
    # patterns (`in_silico_pcr_pool()`), not one per pairing.  the pool is scanned
    # exactly: `near = TRUE` adds a radius-1 pass (pooled near call: 130.05 s
    # on TAIR10 before the `PDict` fix, 1.99 s after; 193 s on a differently
    # loaded machine) for a warning this screen does not need
    screen_pcr <- tryCatch(
      lapply(seq_len(n_screen), function(i) {
        list(assay = screened[i, ], template = tmpl, snp = snp)
      }) |>
        lapply(as_pcr_pairs) |>
        unlist(recursive = FALSE) |>
        stats::setNames(paste0(rep(c("AS_ref", "AS_alt"), n_screen), "_", rep(seq_len(n_screen), each = 2L))) |>
        in_silico_pcr_pool(genome),
      error = function(e) NULL
    )
    screen_pair <- function(i) {
      row <- screened[i, ]
      if (as_bad || any(has_palindrome(c(as_primers, row$com)))) {
        return(FALSE)
      }
      com_hits <- hit_count(c(COM = row$com))
      if (!is.null(com_hits) && nrow(com_hits) && max(table(com_hits$primer)) > 1L) {
        return(FALSE)
      }
      # the pool names each pair's products by its pool position; rename the slice
      # back to the AS primers `assay_amplicon_counts()` reports
      prod <- if (is.null(screen_pcr)) {
        NULL
      } else {
        stats::setNames(
          screen_pcr[paste0(c("AS_ref", "AS_alt"), "_", i)],
          c("AS_ref", "AS_alt")
        )
      }
      if (is.null(prod) || any(vapply(prod, is.null, logical(1)))) {
        return(TRUE) # unknown -> do not block
      }
      cnt <- assay_amplicon_counts(prod)
      ref <- cnt$n_products[cnt$primer == "AS_ref"]
      alt <- cnt$n_products[cnt$primer == "AS_alt"]
      length(ref) == 1L && length(alt) == 1L && ref == 1L && alt == 0L
    }
    clean <- purrr::map_lgl(seq_len(n_screen), screen_pair)
    win <- which(clean)[1L]
    if (is.na(win) && identical(on_offtarget, "doom")) {
      return(marker_doomed(
        gettext("No screened pairing is genome-specific (an off-target hit, a second product or an inverted repeat)."),
        "pair", snp, tmpl, as_cand
      ))
    }
    if (!is.na(win) && win > 1L) {
      cli::cli_warn(c(
        gettext(
          "Top-ranked pairing failed genome-level checks; kept the next clean one ({win} of {n_screen} screened)."
        ),
        i = gettext("Ranking is genome-free; off-target hits, in-silico PCR and inverted repeats are screened here.")
      ))
      ord <- c(win, setdiff(seq_len(nrow(assays)), win))
      assays <- assays[ord, , drop = FALSE]
      assays$rank <- seq_len(nrow(assays))
    }
  }
  best <- assays[1L, ]
  primers <- c(AS_ref = best$as_ref, AS_alt = best$as_alt, COM = best$com)

  # specificity uses the template-matched allele for the AS primers (no
  # deliberate mismatch): that answers "is the target site unique?".
  # make_as_primer() works in TEMPLATE coordinates, hence the template anchor.
  tmpl_anchor <- as_anchor(tmpl$snp_i, ref_len, best$strand, variant_offset)
  seed_primers <- c(
    AS_ref = make_as_primer(tmpl$seq, tmpl_anchor, best$as_len, snp$ref,
      strand = best$strand, variant_offset = variant_offset
    ),
    AS_alt = make_as_primer(tmpl$seq, tmpl_anchor, best$as_len, snp$alt,
      strand = best$strand, variant_offset = variant_offset
    ),
    COM = best$com
  )

  # in-silico PCR: the REF primer must give one product, the ALT primer none.
  # The products are computed once and reused for the amplicon counts, the
  # ideogram below and (via `out$ideogram`) the plots, so the genome is scanned
  # a single time per assay instead of three.
  stub <- list(assay = best, template = tmpl, snp = snp)
  pcr_products <- tryCatch(
    assay_pcr_products(stub, genome, BPPARAM = BPPARAM),
    error = function(e) {
      warn_msg("in-silico PCR failed: {conditionMessage(e)}")
      NULL
    }
  )
  amplicons <- if (is.null(pcr_products)) NULL else assay_amplicon_counts(pcr_products)

  # whole-primer genome scan (only when requested): an attached bowtie index is
  # queried first (index-based, flat cost in genome and primer count), otherwise
  # a BLAST database
  blast <- NULL
  if (isTRUE(run_blast)) {
    blast <- if (inherits(genome, "marker_genome") && !is.null(genome$bowtie_index)) {
      tryCatch(
        specificity_bowtie(primers, genome$bowtie_index,
          nthreads = if (is.null(BPPARAM)) 4L else 1L
        ),
        error = function(e) {
          warn_msg("bowtie scan failed: {conditionMessage(e)}")
          NULL
        }
      )
    } else {
      db <- tryCatch(build_blast_db(genome), error = function(e) {
        warn_msg("BLAST db unavailable: {conditionMessage(e)}")
        NULL
      })
      blast_primers(primers, db)
    }
    # whole-primer scan either way: the AS primers here carry their deliberate
    # mismatch, which the 3'-seed `specificity` table above does not (D3)
    if (!is.null(blast)) attr(blast, "scan") <- "whole-primer"
  }

  # a per-call `engines` override must reach the numbers the evaluation grades:
  # without the structure engine the dG columns are NA, not a silent pass
  thermo <- purrr::imap(primers, ~ dplyr::mutate(primer_thermo(.x, BPPARAM = BPPARAM), name = .y)) |>
    purrr::list_rbind()
  if (!engine_has("hairpin_dg", engines)) thermo$hairpin_dg <- NA_real_
  if (!engine_has("duplex_dg", engines)) thermo$terminal_dg <- NA_real_
  complementarity <- dimer_thermo(primers, BPPARAM = BPPARAM, conditions = conditions)
  if (!engine_has("duplex_dg", engines)) complementarity$dg <- NA_real_

  out <- list(
    scheme = scheme,
    snp = snp, template = tmpl, strand = best$strand, assay = assays,
    limits = limits, as_variant_offset = variant_offset,
    as_candidates = as_cand,
    efficiency = best_as |> dplyr::select("ref_on", "ref_off", "alt_on", "alt_off", "discrimination"),
    thermo = thermo,
    complementarity = complementarity,
    specificity = specificity_of(seed_primers, genome, BPPARAM = BPPARAM),
    amplicons = amplicons,
    selection = list(
      relaxed = pick$relaxed, violated = pick$violated, as = pick_as$violated,
      linked = linked_off, linked_skipped = linked_skipped,
      linked_dependent = linked_dependent
    ),
    blast = blast,
    haplotypes = if (is.null(hap)) NULL else hap$haplotypes,
    haplotype_hits = if (is.null(hap)) {
      NULL
    } else {
      haplotype_hits(
        haplotypes, snp, best$as_start, best$as_end,
        alt_hap = hap$alt,
        alt_start = best$as_alt_start, alt_end = best$as_alt_end,
        alt_baked = lo$link_pos
      )
    },
    tails = c(
      AS_ref = paste0(KASP_TAIL_FAM, best$as_ref),
      AS_alt = paste0(KASP_TAIL_HEX, best$as_alt),
      COM = best$com
    )
  )
  # precompute the whole-genome ideogram inputs once, so plotting only reads the
  # object instead of rescanning the genome.  `primer_hits()` would scan each
  # primer and its reverse complement on its own genome pass; batching them into
  # one pass cuts the ideogram's six passes to one.
  ideogram_hits <- tryCatch(
    batch_primer_hits(primers, genome),
    error = function(e) NULL
  )
  out$ideogram <- tryCatch(
    assay_ideogram_data(out, genome, BPPARAM = BPPARAM, hits = ideogram_hits, products = pcr_products),
    error = function(e) {
      warn_msg("genome ideogram failed: {conditionMessage(e)}")
      NULL
    }
  )
  # the target site's KASP class (biallelic / degenerate / complementary / indel),
  # from every ALT seen at the position in the variant table
  site_alts <- alts_in
  if (!is.null(variants) && all(c("chrom", "pos") %in% names(variants))) {
    rows <- variants[variants$chrom == snp$chrom & variants$pos == snp$pos, , drop = FALSE]
    if (nrow(rows)) {
      a <- if ("ref_genome" %in% names(rows)) {
        # `ref_genome` is a LIST column (one allele vector per variant): unlist
        # first - `as.character()` on a list would stringify the whole vector
        as.character(unlist(rows$ref_genome, use.names = FALSE))
      } else if ("alt" %in% names(rows)) {
        as.character(rows$alt)
      } else {
        character(0)
      }
      a <- a[nzchar(a) & grepl("^[ACGT]+$", a)]
      if (length(a)) site_alts <- unique(c(site_alts, a))
    }
  }
  out$site <- marker_site_class(snp$ref, site_alts)
  # how many of the site's ALTs the single ALT primer actually covers: all of them
  # for a degenerate all-SNP site (one IUPAC primer), otherwise just the first
  out$site$alt_covered <- if (degenerate) out$site$n_alt[[1L]] else 1L
  # a multi-allelic site is reported as such.  When it is all-SNP it is designed
  # as a degenerate primer covering every ALT; otherwise the ALT primer still
  # covers the first ALT only - say so, do not pretend.
  if (out$site$n_alt[1L] > 1L) {
    if (degenerate) {
      cli::cli_inform(c(
        gettext(
          "Target {snp$chrom}:{snp$pos} is multi-allelic ({out$site$n_alt[1L]} ALTs, class {out$site$class[1L]})."
        ),
        i = gettext(
          "The ALT primer is degenerate (3' base {out$site$code[1L]}): it covers every ALT but cannot tell them apart."
        )
      ))
    } else {
      cli::cli_warn(c(
        gettext(
          "Target {snp$chrom}:{snp$pos} is multi-allelic ({out$site$n_alt[1L]} ALTs, class {out$site$class[1L]})."
        ),
        i = gettext(
          "The ALT primer is built for {snp$alt} only; the other ALTs are not covered."
        )
      ))
    }
  }
  # evaluate every registered acceptance check once, so plots/reports read the
  # results instead of recomputing them
  out$evaluation <- tryCatch(marker_design_checks(out, limits = limits, engines = engines), error = function(e) {
    # A failed evaluator is not "no opinion": say so, or the gap looks like a
    # design that has no checks.
    cli::cli_warn(c(
      gettext("The acceptance checks could not be evaluated: {conditionMessage(e)}"),
      i = gettext("assay$evaluation is NULL; the verdict is unknown, not usable.")
    ))
    NULL
  })
  # structural / specificity checks that the selection gates do not cover (they
  # need the genome or a PCR pass, computed only after the pair is chosen):
  # rather than return an ugly assay in silence, name the failures.  Which checks
  # those are is per-scheme (C10) - a scheme that already gates one is not told
  # it is ungated.
  if (!is.null(out$evaluation)) {
    ungated <- KASP_UNGATED_CHECKS(scheme)
    hit <- !out$evaluation$pass & out$evaluation$criterion %in% ungated
    bad <- out$evaluation$criterion[hit]
    if (length(bad)) {
      cli::cli_warn(c(
        gettext(
          "The chosen assay fails structural/specificity checks: {paste(bad, collapse = '; ')}."
        ),
        i = gettext("Selection does not gate these yet - inspect assay$evaluation and consider a trade-off.")
      ))
    }
  }
  # the selection policy of *this* scheme as a graph, with the assay's runtime
  # status stamped on (C12); the object carries the scheme, so a plot or a report
  # can rebuild/match it without guessing
  out$decision <- tryCatch(marker_decision_graph(out, scheme = scheme),
    error = function(e) NULL
  )
  # Y3.1: the two 3' axes, explicit.  The allele axis is structural (always
  # filled); the subgenome axis is read against the caller-supplied mask and is
  # never guessed from the one genome.
  out <- set_homeolog_axes(out, mask)
  contract_assay_output(structure(out, class = c("marker_assay", "marker")), call)
}
