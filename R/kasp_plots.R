# =============================================================================
#  kasp_plots.R  --  library-driven visualisations & oligo ordering helpers
# =============================================================================
#  No hand-rolled geometry: every panel is drawn by a dedicated library.
#
#    - ggmsa (Bioconductor)  : DNA multiple-sequence alignments (+ seq logo)
#    - gggenes               : strand-aware arrows (primers on the template,
#                              gene model + variant track)
#    - ggrepel               : non-overlapping labels
#    - OligoArrayAux         : DNA duplex / hairpin free energies (via oligo_dg())
#
#  Hairpins / self-complementarity are shown the DNA way: an alignment of each
#  primer against its own reverse complement (ggmsa). RNA-only tools (R4RNA /
#  ggmsa::geom_helix) are deliberately NOT used - we work with DNA.
#
#  export_primers() writes tidy order sheets for IDT (bulk oligo input) and
#  LGC/KASP (assay order).
# =============================================================================

#' The three assay primers as a named vector
#' @keywords internal
#' @noRd
assay_primers <- function(assay) {
  c(AS_ref = assay$assay$as_ref[1], AS_alt = assay$assay$as_alt[1], COM = assay$assay$com[1])
}

#' Design window (plus-strand genomic coords) as c(start, end)
#' @keywords internal
#' @noRd
design_window <- function(assay, snp) {
  flank <- (length(assay$template$seq) - 1L) / 2L
  c(snp$pos - flank, snp$pos + flank)
}

#' Variants inside a region, as a tidy tibble
#' @keywords internal
#' @noRd
variants_in_window <- function(variants, chrom, win) {
  variants |>
    dplyr::filter(.data$chrom == !!chrom,
                  dplyr::between(.data$pos, win[1], win[2])) |>
    dplyr::mutate(consequence = forcats::fct_lump_n(.data$consequence, 6L))
}


# =============================================================================
#  Alignments (ggmsa)
# =============================================================================

#' AS primers aligned against the template; SNP + mismatches visible
#'
#' A DNA MSA rendered by ggmsa, with a sequence logo on top.  The two template
#' rows are the primer-binding sequences carrying the REF/ALT allele (without the
#' deliberate stabilising mismatch); the two AS rows are the real primers, so a
#' column that differs between a template row and its AS row is the injected
#' secondary mismatch.  Works for either assay orientation.
#' @export
plot_alignment <- function(assay) {
  tmpl <- assay$template$seq; snp_i <- assay$template$snp_i
  p <- assay_primers(assay); L <- nchar(p[["AS_ref"]])
  strand <- assay$strand %||% "+"
  anchor <- as_anchor(snp_i, nchar(assay$snp$ref), strand)
  # template rows in the SAME 5'->3' orientation as the primers, one per allele
  t_ref <- make_as_primer(tmpl, anchor, L, assay$snp$ref, strand = strand)
  t_alt <- make_as_primer(tmpl, anchor, L, assay$snp$alt, strand = strand)
  msa <- Biostrings::DNAStringSet(c(template_ref = t_ref, template_alt = t_alt,
                                    AS_ref = p[["AS_ref"]], AS_alt = p[["AS_alt"]]))
  suppressWarnings(
    ggmsa::ggmsa(msa, 1, L, color = "Chemistry_NT", seq_name = TRUE, char_width = 0.6) +
      ggmsa::geom_seqlogo(color = "Chemistry_NT"))
}

#' Self-complementarity of each primer (DNA): primer vs its reverse complement
#'
#' Columns that match are self-complementary (hairpin-prone stems).
#' @export
plot_self_complementarity <- function(assay, title = "Self-complementarity (DNA): matching columns = hairpin-prone stems") {
  plots <- purrr::imap(assay_primers(assay), function(s, nm) {
    rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
    suppressWarnings(
      ggmsa::ggmsa(Biostrings::DNAStringSet(c(primer = s, revcomp = rc)),
                   1, nchar(s), color = "Chemistry_NT", seq_name = TRUE, char_width = 0.5)) +
      ggplot2::theme(plot.margin = ggplot2::margin(1, 4, 1, 4))
  })
  p <- patchwork::wrap_plots(plots, ncol = 1)
  if (!is.null(title)) p <- p + patchwork::plot_annotation(title = title)
  p
}


# =============================================================================
#  Candidate ranking
# =============================================================================

#' Candidate AS primers ranked by discrimination
#'
#' One bar chart, shared by the driver and the report so the two cannot drift.
#' @export
plot_discrimination <- function(assay, n = 15L) {
  assay$as_candidates |>
    dplyr::slice_max(.data$discrimination, n = n) |>
    dplyr::mutate(label = paste0("strand ", .data$strand, " / len ", .data$as_length,
                                 " / mm ", .data$stab_at, " / ", .data$stab_base)) |>
    ggplot2::ggplot(ggplot2::aes(.data$discrimination,
                                 stats::reorder(.data$label, .data$discrimination))) +
    ggplot2::geom_col(fill = "darkorange") +
    ggplot2::labs(x = "discrimination (on- minus off-target efficiency)", y = NULL,
                  title = "Allele-specific discrimination of candidate AS primers") +
    ggplot2::theme_minimal()
}


# =============================================================================
#  Gene / primer maps (gggenes)
# =============================================================================

#' Primer annealing on the amplicon + local 1001-Genomes variant track
#'
#' Top: gggenes arrows on the design window (genomic coordinates), pointing in
#' the direction of extension.  Bottom: the variants that fall inside that
#' window, on the same x axis.  The assay orientation (`+`/`-`) decides which
#' side of the SNP the primers sit on.
#' @export
plot_primers <- function(assay, variants, snp) {
  win <- design_window(assay, snp)
  aa  <- assay$assay[1, ]
  strand <- aa$strand %||% "+"
  plus   <- identical(strand, "+")

  # coordinates come straight from the assay (single source of truth)
  d <- tibble::tibble(
    primer  = factor(c("AS-ref", "AS-alt", "COM"), levels = c("COM", "AS-alt", "AS-ref")),
    start   = c(aa$as_start, aa$as_start, aa$com_start),
    end     = c(aa$as_end,   aa$as_end,   aa$com_end),
    forward = c(plus, plus, !plus))

  p_arrows <- ggplot2::ggplot(d, ggplot2::aes(xmin = .data$start, xmax = .data$end,
                                              y = .data$primer, forward = .data$forward,
                                              fill = .data$primer)) +
    gggenes::geom_gene_arrow() +
    ggplot2::geom_vline(xintercept = snp$pos, linetype = "dashed", colour = "firebrick") +
    ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    gggenes::theme_genes()

  vw <- variants_in_window(variants, snp$chrom, win)
  reg3p    <- c(aa$win3p_start, aa$win3p_end)
  as_span  <- c(aa$as_start,   aa$as_end)
  com_span <- c(aa$com_start,  aa$com_end)
  foot <- dplyr::filter(vw, (.data$pos >= as_span[1] & .data$pos <= as_span[2]) |
                              (.data$pos >= com_span[1] & .data$pos <= com_span[2]))

  p_var <- ggplot2::ggplot(vw, ggplot2::aes(.data$pos, 1)) +
    ggplot2::annotate("rect", xmin = reg3p[1] - 0.5, xmax = reg3p[2] + 0.5,
                      ymin = 0.4, ymax = 1.6, fill = "firebrick", alpha = 0.12) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$consequence), size = 1.7) +
    ggplot2::geom_point(data = foot, ggplot2::aes(.data$pos, 1),
                        shape = 1, colour = "black", size = 2.6) +
    ggplot2::geom_point(data = snp, ggplot2::aes(.data$pos, 1),
                        shape = 21, fill = "firebrick", colour = "white", size = 3) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "consequence") +
    ggplot2::scale_y_continuous(breaks = 1, labels = "variants") +
    ggplot2::labs(x = glue::glue("design window, {snp$chrom} position (bp)"), y = NULL,
                  caption = "shaded band = AS 3' 5 nt; ringed = other variants under AS/COM") +
    ggplot2::theme_minimal()

  patchwork::wrap_plots(p_arrows, p_var, ncol = 1, heights = c(2, 1)) +
    patchwork::plot_annotation(
      title = "KASP primers + local variants (1001 Genomes)",
      subtitle = glue::glue("target {snp$chrom}:{snp$pos} ({snp$ref}>{snp$alt}) on strand ",
                            "{strand}; other variants: AS 3'={aa$as_n_3p %||% NA}, ",
                            "AS footprint={aa$as_n_var}, COM={aa$n_var}"))
}

#' Gene model with the 1001-Genomes variant track + inset of the target region
#' @export
plot_variant_landscape <- function(gene, variants, snp, assay) {
  gene_df <- tibble::tibble(molecule = gene$symbol, start = gene$start, end = gene$end,
                            strand = gene$strand)
  v <- variants |>
    dplyr::mutate(consequence = forcats::fct_lump_n(.data$consequence, 6L))

  main <- ggplot2::ggplot() +
    gggenes::geom_gene_arrow(data = gene_df,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = 1, forward = .data$strand > 0),
      fill = "grey80") +
    ggplot2::geom_point(data = v, ggplot2::aes(.data$pos, 2, colour = .data$consequence),
                        size = 1.3, shape = 108) +
    ggplot2::geom_point(data = snp, ggplot2::aes(.data$pos, 2),
                        shape = 21, fill = "firebrick", colour = "white", size = 3) +
    ggrepel::geom_text_repel(data = snp, ggplot2::aes(.data$pos, 2, label = .data$id),
                             size = 2.6, nudge_y = 0.6, min.segment.length = 0) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "consequence") +
    ggplot2::scale_y_continuous(limits = c(0.5, 3), breaks = c(1, 2),
                                labels = c(gene$symbol, "variants")) +
    ggplot2::labs(x = glue::glue("{gene$chrom} position (bp)"), y = NULL,
                  title = glue::glue("Variation across {gene$symbol} (1001 Genomes)")) +
    ggplot2::theme_minimal()

  win <- design_window(assay, snp)
  vw  <- variants_in_window(variants, snp$chrom, win)
  inset <- ggplot2::ggplot(vw, ggplot2::aes(.data$pos, 1)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$consequence), size = 1.1) +
    ggplot2::geom_point(data = snp, ggplot2::aes(.data$pos, 1),
                        shape = 21, fill = "firebrick", colour = "white", size = 2) +
    ggplot2::scale_colour_brewer(palette = "Dark2", guide = "none") +
    ggplot2::scale_y_continuous(breaks = 1, labels = "variants") +
    ggplot2::labs(x = glue::glue("design window ({win[1]}-{win[2]} bp)"), y = NULL,
                  title = "target region") +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = "grey70"))

  main + patchwork::inset_element(inset, left = 0.53, bottom = 0.42, right = 1, top = 1)
}


# =============================================================================
#  Cross-complementarity network (ggraph / tidygraph)
# =============================================================================

#' Dimer QC: pairwise duplex dG heatmap + hairpin dG (OligoArrayAux, DNA)
#'
#' Readable replacement for a 3-node network: a symmetric matrix of duplex
#' free energies (kcal/mol, 37 C; more negative = stronger dimer) plus a bar
#' panel of intramolecular folding (hairpin) free energies.
#' @export
plot_dimer_qc <- function(assay) {
  d <- assay$complementarity          # already a full symmetric matrix

  heat <- ggplot2::ggplot(d, ggplot2::aes(.data$a, .data$b, fill = .data$dg)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = round(.data$dg, 1)), size = 4) +
    ggplot2::scale_fill_gradient(low = "firebrick", high = "white",
                                 limits = c(min(d$dg, na.rm = TRUE), 0),
                                 name = "dG37\n(kcal/mol)") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, title = "Duplex dG37 (kcal/mol)") +
    ggplot2::theme_minimal()

  hair <- ggplot2::ggplot(assay$thermo,
                          ggplot2::aes(.data$hairpin_dg, stats::reorder(.data$name, .data$hairpin_dg))) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_vline(xintercept = HAIRPIN_DG_WARN, linetype = "dashed", colour = "firebrick") +
    ggplot2::labs(x = "hairpin dG37 (kcal/mol)", y = NULL,
                  title = "Hairpin dG37 (kcal/mol)") +
    ggplot2::theme_minimal()

  sc <- plot_self_complementarity(assay, title = NULL)
  (heat | hair) / sc +
    patchwork::plot_layout(heights = c(1, 1.1)) +
    patchwork::plot_annotation(title = "Primer QC: duplex/hairpin dG37 + self-complementarity")
}


# =============================================================================
#  Oligo ordering (IDT / LGC)
# =============================================================================

#' Tidy table of the assay oligos (with KASP tails)
#' @export
kasp_primer_table <- function(assay, gene) {
  nm <- glue::glue("{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}_{assay$snp$ref}to{assay$snp$alt}")
  roles <- tibble::tibble(
    name        = c("AS_ref", "AS_alt", "COM"),
    role        = c("AS-ref (FAM)", "AS-alt (HEX/VIC)", "common (reverse)"),
    oligo_id    = c("AS1_FAM", "AS2_HEX", "COM"),
    tail        = c(KASP_TAIL_FAM, KASP_TAIL_HEX, ""),
    seq_no_tail = c(assay$assay$as_ref[1], assay$assay$as_alt[1], assay$assay$com[1]))

  assay$thermo |>
    dplyr::select("name", "tm", "gc") |>
    dplyr::inner_join(roles, by = "name") |>               # join, not positional match()
    dplyr::mutate(assay    = nm,
                  sequence = paste0(.data$tail, .data$seq_no_tail),  # with 5' tail
                  gc       = round(.data$gc, 2),
                  length   = nchar(.data$sequence)) |>
    dplyr::select("assay", "role", "oligo_id", "sequence", "seq_no_tail",
                  "tail", "tm", "gc", "length")
}

#' Write IDT, LGC and tidy order sheets for the assay primers
#'
#' - `*_idt.csv`  : IDT bulk oligo input (Name, Sequence, Scale, Purification)
#' - `*_lgc.csv`  : LGC/KASP assay order (one row per assay, wide)
#' - `*_tidy.csv` : the long tidy table
#' @export
export_primers <- function(assay, gene, dir = here::here("results")) {
  fs::dir_create(dir)
  tbl <- kasp_primer_table(assay, gene)
  base <- glue::glue("kasp_{gene$symbol}_{assay$snp$chrom}_{assay$snp$pos}")

  tbl |>
    dplyr::transmute(Name = glue::glue("{assay}_{oligo_id}"),
                     Sequence = .data$sequence, Scale = "25nm", Purification = "STD") |>
    readr::write_csv(fs::path(dir, glue::glue("{base}_idt.csv")))

  tibble::tibble(
    Assay      = unique(tbl$assay),
    AS1_FAM    = tbl$sequence[tbl$oligo_id == "AS1_FAM"],
    AS2_HEX    = tbl$sequence[tbl$oligo_id == "AS2_HEX"],
    Common     = tbl$sequence[tbl$oligo_id == "COM"]) |>
    readr::write_csv(fs::path(dir, glue::glue("{base}_lgc.csv")))

  readr::write_csv(tbl, fs::path(dir, glue::glue("{base}_tidy.csv")))
  invisible(tbl)
}
