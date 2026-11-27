# =============================================================================
#  marker_assemble_fixed.R  --  assemble_primer_assay() from given primers
# =============================================================================
#  Split from marker_assemble.R (roadmap ZT T7).
# =============================================================================
# =============================================================================
#  7c. Assemble an assay from FIXED (given) primers
# =============================================================================

# The population / footprint checks need a variant table; the fixed-primer
# constructor never invents one, so they are reported `unknown`, never passed.
FIXED_PRIMER_NO_VARIANTS <- c(
  "linked variants skipped (n)", "design relies on linkage",
  "AS 3' window variants (n)", "AS footprint variants (n)",
  "COM footprint variants (n)", "AS 3' window dropout",
  "AS footprint dropout", "COM footprint dropout"
)
# the subset whose pass condition treats an NA dropout as a pass
FIXED_PRIMER_DROPOUTS <- c(
  "AS 3' window dropout", "AS footprint dropout", "COM footprint dropout"
)

#' @keywords internal
#' @noRd
normalise_fixed_snp <- function(snp) {
  if (is.list(snp) && !is.data.frame(snp)) snp <- tibble::as_tibble(snp)
  if (!is.data.frame(snp) || nrow(snp) != 1L) {
    abort_msg("`snp` must be a one-row table (or list) with `chrom`, `pos`, `ref`, `alt`.")
  }
  need_cols(snp, c("chrom", "pos", "ref", "alt"), arg = "snp")
  out <- tibble::tibble(
    chrom = as.character(snp$chrom[[1L]]),
    pos = as.integer(snp$pos[[1L]]),
    ref = toupper(as.character(snp$ref[[1L]])),
    alt = toupper(as.character(snp$alt[[1L]]))
  )
  if (is.na(out$pos) || !nzchar(out$chrom)) {
    abort_msg("`snp` needs a non-missing chromosome and a positive integer position.")
  }
  if (!grepl("^[ACGT]$", out$ref) || !grepl("^[ACGT]$", out$alt)) {
    abort_msg("Fixed-primer assembly currently handles single-base substitutions only.")
  }
  if (identical(out$ref, out$alt)) {
    abort_msg("`snp$ref` and `snp$alt` must differ.")
  }
  out
}

#' @keywords internal
#' @noRd
normalise_fixed_primers <- function(primers) {
  aliases <- c(
    AS_ref = "AS_ref", as_ref = "AS_ref", allele1 = "AS_ref", ref = "AS_ref",
    AS_alt = "AS_alt", as_alt = "AS_alt", allele2 = "AS_alt", alt = "AS_alt",
    COM = "COM", com = "COM", common = "COM"
  )
  if (is.data.frame(primers)) {
    role_col <- intersect(c("role", "name"), names(primers))[1L]
    # `body` wins when both are present: it is the genome-matching part, while a
    # `sequence` column may still carry the 5' KASP tail
    seq_col <- intersect(c("body", "sequence", "seq"), names(primers))[1L]
    if (is.na(role_col) || is.na(seq_col)) {
      abort_msg("`primers` needs a `role` (or `name`) and a `sequence` (or `body`) column.")
    }
    role <- as.character(primers[[role_col]])
    seqs <- as.character(primers[[seq_col]])
  } else if (is.character(primers) && !is.null(names(primers))) {
    role <- names(primers)
    seqs <- unname(primers)
  } else {
    abort_msg("`primers` must be a table with `role`/`sequence`, or a named character vector.")
  }
  idx <- match(role, names(aliases))
  role[!is.na(idx)] <- aliases[idx[!is.na(idx)]]
  seqs <- toupper(gsub("\\s+", "", seqs))
  if (length(seqs) != 3L || anyDuplicated(role)) {
    abort_msg("`primers` must be exactly three distinct roles: AS_ref, AS_alt and COM.")
  }
  if (!setequal(role, c("AS_ref", "AS_alt", "COM"))) {
    abort_msg("`primers` roles must be AS_ref / AS_alt / COM (aliases: allele1 / allele2 / common).")
  }
  p <- rlang::set_names(seqs, role)[c("AS_ref", "AS_alt", "COM")]
  if (any(!nzchar(p)) || any(grepl("[^ACGT]", p))) {
    abort_msg("Every primer must be a non-empty ACGT sequence.")
  }
  p
}

#' Assemble an assay from fixed (already given) primers
#'
#' [assemble_kasp_assay()] *designs* the oligos.  This constructor starts from a
#' primer set that already exists and turns it into the same `marker_assay`
#' object, so [marker_design_checks()], [graphics::plot()] / [plot_primers()] and the
#' decision graph read it like any designed assay.  Nothing is taken on trust:
#' every oligo is located by an exact full-length genome hit, its binding strand
#' and the REF/ALT orientation are read from the genome, and the 3'-terminal base
#' is checked against the target allele.  The amplicon, the in-silico PCR and the
#' off-target counts are then computed, not copied from a document.
#'
#' The primers are the *bodies*: they must match the genome exactly (a 5' tail by
#' definition does not, so strip it first).  `AS_ref` and `AS_alt` must be the
#' same length and differ at exactly their 3' base, which is the allele's base on
#' the binding strand.
#'
#' The panel has no population variant table, so the footprint / linkage checks
#' are reported `unknown` unless `variants` is supplied; they are never silently
#' passed.
#'
#' @param genome a `marker_genome` or a `DNAStringSet`.
#' @param snp a one-row table with `chrom`, `pos`, `ref` and `alt` (the
#'   plus-strand alleles).
#' @param primers a table with a `role` column (`AS_ref` / `AS_alt` / `COM`, or the
#'   aliases `allele1` / `allele2` / `common`) and a `sequence` (or `body`)
#'   column, or a named character vector of the three bodies.
#' @param variants optional variant table (`chrom`, `pos`, `alt`, and optionally
#'   `pop_freq`); it enables the footprint / linkage checks.  `NULL` leaves them
#'   `unknown`.
#' @param flank minimum template half-width (bp); widened automatically to cover
#'   the amplicon.
#' @param scheme the marker chemistry to assemble *as* (default `"snp"`), which
#'   selects the reaction conditions and the graded criterion set.
#' @param salt_method salt correction for `Tm`.  The default is the design
#'   coordinate `TM_SALT_METHOD` (`"Schildkraut2010"`), the same scale the
#'   KASP/AS/COM gate windows are stated on, so the assay is graded on the
#'   coordinate its limits live in.  Pass `"SantaLucia1998-2"` to reproduce a
#'   printed Primer3 panel (ROADMAP `ZW3.5`, `ZW3.10`).
#' @param cores worker count for the genome/thermodynamics batches.
#' @param BPPARAM optional `BiocParallel` backend (overrides `cores`).
#' @param diagnostic_positions the subgenome-diagnostic mask (roadmap Y3.1):
#'   `NULL`, a data frame / list with `chrom` and `pos`, a `GRanges`, or a named
#'   list `chrom -> integer positions`.  See [assemble_kasp_assay()].
#' @param engines optional `marker_engines` adapter list (see `engine_has()`);
#'   `NULL` reads the session default.  It gates the footprint-dropout screen
#'   and the discrimination score, both of which report `"unknown"` without the
#'   PCR-efficiency engine.
#' @return a `marker_assay`.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' tmpl <- snp_template(genome, "Chr4", 2844L, 150L)
#' as_ref <- as.character(Biostrings::subseq(tmpl$seq, tmpl$snp_i - 24L, tmpl$snp_i))
#' com <- as.character(Biostrings::reverseComplement(
#'   Biostrings::subseq(tmpl$seq, tmpl$snp_i + 60L, tmpl$snp_i + 84L)
#' ))
#' p <- tibble::tibble(
#'   role = c("AS_ref", "AS_alt", "COM"),
#'   sequence = c(as_ref, paste0(substr(as_ref, 1L, 24L), "C"), com)
#' )
#' snp <- tibble::tibble(chrom = "Chr4", pos = 2844L, ref = "T", alt = "C")
#' assay <- assemble_primer_assay(genome, snp, p)
#' assay$assay[, c("as_ref", "as_alt", "com", "amplicon")]
#' @family design
#' @export
assemble_primer_assay <- function(genome, snp, primers, variants = NULL,
                                  flank = 150L, scheme = "snp",
                                  salt_method = TM_SALT_METHOD,
                                  cores = getOption("tidymarkers.cores", 1L),
                                  BPPARAM = NULL,
                                  diagnostic_positions = NULL, engines = NULL) {
  if (is.null(BPPARAM)) BPPARAM <- marker_bpparam(cores)
  mask <- normalise_diagnostic_positions(diagnostic_positions)
  conditions <- scheme_conditions(scheme)
  snp <- normalise_fixed_snp(snp)
  p <- normalise_fixed_primers(primers)
  as_ref <- p[["AS_ref"]]
  as_alt <- p[["AS_alt"]]
  com <- p[["COM"]]
  n_as <- nchar(as_ref)
  if (nchar(as_alt) != n_as) {
    abort_msg("AS_ref and AS_alt must be the same length (they are one allele-specific primer pair).")
  }
  n_pre <- n_as - 1L
  same_prefix <- substr(as_ref, 1L, n_pre) == substr(as_alt, 1L, n_pre)
  same_last <- substr(as_ref, n_as, n_as) == substr(as_alt, n_as, n_as)
  if (!same_prefix || same_last) {
    abort_msg("AS_ref and AS_alt must differ at exactly their 3'-terminal base.")
  }

  # --- locate the oligos on the genome (never accept coordinates from prose) ---
  hits <- exact_hit_table(p, genome, BPPARAM = BPPARAM)
  h <- hits[hits$primer == "AS_ref", , drop = FALSE]
  anchor <- ifelse(h$strand == "+", h$end, h$start)
  cand <- h[h$chrom == snp$chrom & anchor == snp$pos, , drop = FALSE]
  if (nrow(cand) != 1L) {
    abort_msg(
      "The AS_ref primer must have exactly one exact genome hit whose 3' end sits on the target.",
      x = "{snp$chrom}:{snp$pos}: found {nrow(cand)} such hit(s).",
      i = "AS_ref (5'->3'): {as_ref}"
    )
  }
  strand <- cand$strand[[1L]]
  as_start <- cand$start[[1L]]
  as_end <- cand$end[[1L]]

  # the 3' base pairs with the template: on "+" it is the plus-strand allele, on
  # "-" it is its complement
  last_base <- function(x) substr(x, nchar(x), nchar(x))
  expected <- function(allele) {
    if (identical(strand, "+")) allele else as.character(Biostrings::complement(Biostrings::DNAString(allele)))
  }
  if (!identical(last_base(as_ref), expected(snp$ref))) {
    abort_msg(
      "AS_ref's 3' base does not match the reference allele at the target.",
      x = "target {snp$chrom}:{snp$pos} has {snp$ref}; strand {strand} expects {expected(snp$ref)}."
    )
  }
  if (!identical(last_base(as_alt), expected(snp$alt))) {
    abort_msg("AS_alt's 3' base does not match the ALT allele ({snp$alt}) on strand {strand}.")
  }

  # the common primer: an opposite-strand hit on the target chromosome that
  # closes the amplicon (the shortest span that contains the SNP)
  hc <- hits[hits$primer == "COM", , drop = FALSE]
  hc <- hc[hc$chrom == snp$chrom, , drop = FALSE]
  if (nrow(hc)) {
    spans <- pmin(as_start, hc$start) <= snp$pos & pmax(as_end, hc$end) >= snp$pos
    closes <- if (identical(strand, "+")) hc$start > as_end else hc$end < as_start
    hc <- hc[hc$strand != strand & spans & closes, , drop = FALSE]
  }
  if (!nrow(hc)) {
    abort_msg(
      "The COM primer has no exact opposite-strand hit that closes an amplicon on {snp$chrom}.",
      i = "COM (5'->3'): {com}"
    )
  }
  span_len <- pmax(as_end, hc$end) - pmin(as_start, hc$start) + 1L
  best <- which(span_len == min(span_len))
  if (length(best) > 1L) {
    abort_msg("The COM primer has more than one equally short amplicon on {snp$chrom}; ambiguous.")
  }
  com_hit <- hc[best, , drop = FALSE]
  com_start <- com_hit$start[[1L]]
  com_end <- com_hit$end[[1L]]
  amplicon <- max(as_end, com_end) - min(as_start, com_start) + 1L

  # --- template (symmetric, so `design_window()` reads the right flank) -------
  needed <- max(abs(c(as_start, as_end, com_start, com_end) - snp$pos)) + 10L
  half <- max(as.integer(flank), as.integer(needed))
  if (snp$pos - half < 1L) {
    abort_msg("The amplicon runs off the start of {snp$chrom}; cannot build a template around {snp$pos}.")
  }
  tmpl <- snp_template(genome, snp$chrom, snp$pos, half)
  if (!identical(as.character(Biostrings::subseq(tmpl$seq, tmpl$snp_i, tmpl$snp_i)), snp$ref)) {
    abort_msg("The genome at {snp$chrom}:{snp$pos} is not the supplied REF allele {snp$ref}.")
  }

  win3p <- as_3p_interval(snp$pos, strand)
  win3p_start <- IRanges::start(win3p)
  win3p_end <- IRanges::end(win3p)

  # --- thermodynamics (the same engine the designer uses) ---------------------
  thermo <- purrr::imap(p, function(s, nm) {
    dplyr::mutate(
      primer_thermo(s, BPPARAM = BPPARAM),
      name = nm,
      tm = primer_tm(s, salt_method = salt_method)
    )
  }) |>
    purrr::list_rbind()
  com_tm <- thermo$tm[thermo$name == "COM"]

  # --- discrimination over the AS binding window carrying REF vs ALT ----------
  seg <- as.character(genome_seq(genome, snp$chrom, as_start, as_end))
  anchor_seg <- if (identical(strand, "+")) nchar(seg) else 1L
  if (substr(seg, anchor_seg, anchor_seg) != snp$ref) {
    abort_msg("The AS_ref hit does not place {snp$ref} at its 3' base.")
  }
  seg_alt <- seg
  substr(seg_alt, anchor_seg, anchor_seg) <- snp$alt
  t_ref <- as.character(as_target(Biostrings::DNAString(seg), strand))
  t_alt <- as.character(as_target(Biostrings::DNAString(seg_alt), strand))
  # unknown, not a pass, when the engine is absent: pcr_efficiency() aborts
  use_eff <- engine_has("pcr_efficiency", engines)
  eff <- if (use_eff) {
    pcr_efficiency(
      c(as_ref, as_ref, as_alt, as_alt),
      c(t_ref, t_alt, t_alt, t_ref)
    )
  }
  discrimination <- if (use_eff) {
    ((eff[[1L]] - eff[[2L]]) + (eff[[3L]] - eff[[4L]])) / 2
  } else {
    NA_real_
  }

  # --- footprint / linkage metrics (unknown without a variant table) ----------
  if (is.null(variants)) {
    as_n_3p <- NA_integer_
    as_n_var <- NA_integer_
    n_var <- NA_integer_
    as_dropout <- NA_real_
    as_foot_dropout <- NA_real_
    com_dropout <- NA_real_
  } else {
    off <- snp$pos - half
    vc <- marker_variant_alleles(variants, snp, tmpl$seq, off)
    vpos <- vc$vpos
    vblock <- snp$pos
    batch_or_safe <- marker_dropout_batch(tmpl, off, vpos, vc$valt, vc$vweight, use_eff, BPPARAM)
    as_n_3p <- count_variants_in(win3p, vpos, vblock)
    as_n_var <- count_variants_in(IRanges::IRanges(as_start, as_end), vpos, vblock)
    n_var <- count_variants_in(IRanges::IRanges(com_start, com_end), vpos, vblock)
    as_dropout <- batch_or_safe(as_ref, win3p_start, win3p_end, strand == "+", as_n_3p, vblock)
    as_foot_dropout <- batch_or_safe(as_ref, as_start, as_end, strand == "+", as_n_var, vblock)
    com_dropout <- batch_or_safe(com, com_start, com_end, strand != "+", n_var, vblock)
  }

  assay_row <- tibble::tibble(
    amplicon = amplicon,
    strand = strand,
    com_tm = com_tm,
    as_ref = as_ref, as_alt = as_alt, com = com,
    as_gc = gc_frac(as_ref), com_gc = gc_frac(com),
    as_start = as_start, as_end = as_end,
    as_alt_start = as_start, as_alt_end = as_end,
    com_start = com_start, com_end = com_end,
    win3p_start = win3p_start, win3p_end = win3p_end,
    as_len = n_as,
    as_n_3p = as_n_3p, as_dropout = as_dropout,
    as_n_var = as_n_var, as_foot_dropout = as_foot_dropout,
    n_var = n_var, com_dropout = com_dropout
  )

  site <- marker_site_class(snp$ref, snp$alt)
  site$alt_covered <- 1L

  out <- list(
    scheme = scheme,
    snp = snp,
    template = tmpl,
    strand = strand,
    assay = assay_row,
    as_candidates = NULL,
    efficiency = tibble::tibble(
      ref_on = eff[[1L]], ref_off = eff[[2L]],
      alt_on = eff[[3L]], alt_off = eff[[4L]],
      discrimination = discrimination
    ),
    thermo = thermo,
    complementarity = dimer_thermo(p, BPPARAM = BPPARAM, conditions = conditions),
    amplicons = NULL,
    selection = NULL,
    blast = NULL,
    haplotypes = NULL,
    haplotype_hits = NULL,
    site = site
  )

  # in-silico PCR and the whole-genome ideogram (computed once, reused)
  stub <- list(assay = assay_row, template = tmpl, snp = snp)
  pcr_products <- tryCatch(
    assay_pcr_products(stub, genome, BPPARAM = BPPARAM),
    error = function(e) {
      warn_msg("in-silico PCR failed: {conditionMessage(e)}")
      NULL
    }
  )
  out$amplicons <- if (is.null(pcr_products)) NULL else assay_amplicon_counts(pcr_products)
  out$ideogram <- tryCatch(
    assay_ideogram_data(out, genome, BPPARAM = BPPARAM, products = pcr_products),
    error = function(e) {
      warn_msg("genome ideogram failed: {conditionMessage(e)}")
      NULL
    }
  )

  # grade the assay with the same registry the designer uses
  out$evaluation <- tryCatch(
    marker_design_checks(out, scheme = scheme, engines = engines),
    error = function(e) {
      cli::cli_warn(c(
        gettext("The acceptance checks could not be evaluated: {conditionMessage(e)}"),
        i = gettext("assay$evaluation is NULL; the verdict is unknown, not usable.")
      ))
      NULL
    }
  )
  if (!is.null(out$evaluation)) {
    # a check that could not be computed must never read as a pass (an NA
    # dropout would otherwise satisfy `is.na(dropout) | ...`)
    unknown <- if (is.null(variants)) {
      FIXED_PRIMER_NO_VARIANTS
    } else if (!use_eff) {
      FIXED_PRIMER_DROPOUTS
    } else {
      character(0)
    }
    # an absent engine also leaves the discrimination unknown
    if (!use_eff) unknown <- c(unknown, "discrimination")
    if (length(unknown)) {
      i <- out$evaluation$criterion %in% unknown
      out$evaluation$observed[i] <- NA_real_
      # `pass` is a view of `status` (marker_evaluation_table), so it is rebuilt,
      # not set to NA here
      out$evaluation$pass[i] <- FALSE
      out$evaluation$status[i] <- "unknown"
      no_var <- out$evaluation$criterion %in% FIXED_PRIMER_NO_VARIANTS
      out$evaluation$detail[i] <- ifelse(no_var[i],
        "not evaluated: no population variant table was supplied",
        "not evaluated: OligoArrayAux is unavailable"
      )
    }
  }
  # no selection ladder ran on a fixed primer set, so the graph is left neutral
  # rather than pretending every gate was cleared
  out$decision <- tryCatch(marker_decision_graph(out, scheme = scheme), error = function(e) NULL)
  # Y3.1: the two 3' axes, explicit, read from the fixed primers' own geometry
  # (com_start/com_end) rather than from the designer's pair search.
  out <- set_homeolog_axes(out, mask)
  structure(out, class = c("marker_assay", "marker"))
}
