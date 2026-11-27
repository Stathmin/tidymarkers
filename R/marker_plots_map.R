# =============================================================================
#  marker_plots_map.R  --  gene/primer maps, haplotype layers, plot_primers()
# =============================================================================
#  Split from marker_plots.R (roadmap ZT T7); the plotting vocabulary and the
#  ideogram live in marker_plots.R, the order/output helpers in
#  marker_plots_output.R.
# =============================================================================
# =============================================================================
#  Gene / primer maps (gggenes)
# =============================================================================

#' Normalise an annotation table to the feature row's shape
#'
#' Accepts a data frame with `start`/`end` (and optional `strand`, `name`/`symbol`),
#' a path to an annotation file ([marker_features()]), or a gene record, keeps only
#' features overlapping `span` (and, when a `chrom` column is present, `chrom`),
#' and returns `start`/`end`/`forward`/`name` (or `NULL`).  One place so the doomed
#' and the designed views cannot drift.
#' @keywords internal
#' @noRd
normalise_features <- function(features, span, chrom = NULL) {
  if (is.character(features) && length(features) == 1L && fs::file_exists(features)) {
    features <- marker_features(features, chrom = chrom)
  }
  if (is.null(features)) {
    return(NULL)
  }
  features <- tibble::as_tibble(features)
  if (!nrow(features) || !all(c("start", "end") %in% names(features))) {
    return(NULL)
  }
  # the feature file is genome-wide; keep the chromosome being drawn
  if (!is.null(chrom) && "chrom" %in% names(features)) {
    features <- features[
      !is.na(features$chrom) & as.character(features$chrom) == as.character(chrom), ,
      drop = FALSE
    ]
  }
  nm <- if ("name" %in% names(features)) {
    features$name
  } else if ("symbol" %in% names(features)) {
    features$symbol
  } else {
    NA_character_
  }
  feat <- tibble::tibble(
    start = as.integer(features$start), end = as.integer(features$end),
    forward = if ("strand" %in% names(features)) features$strand != "-" else TRUE,
    name = nm
  )
  feat <- feat[feat$end >= span[1L] & feat$start <= span[2L], , drop = FALSE]
  if (nrow(feat) == 0L) NULL else feat
}

#' Context plot for a doomed target
#'
#' A target no assay clears the gates for still deserves a picture: the local
#' window, its variants and the feature annotation, on the same x axis as
#' [plot_primers()], with the reason it was doomed.  Drawn by [plot_primers()]
#' when it is handed a `marker_doomed` object.
#' @keywords internal
#' @noRd
plot_doomed_target <- function(assay, variants, snp, features = NULL) {
  win <- design_window(assay, snp)
  vw <- variants[
    variants$chrom == snp$chrom & variants$pos >= win[1L] & variants$pos <= win[2L], ,
    drop = FALSE
  ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = tibble::tibble(x = win[1L], xend = win[2L]),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = 1, yend = 1),
      colour = "grey72", linewidth = 0.3
    )
  feat <- normalise_features(features, win, chrom = snp$chrom)
  if (is.null(feat)) {
    p <- p + ggplot2::annotate("text",
      x = mean(win), y = 1, label = "no feature annotation",
      size = 2.1, colour = "grey55", fontface = "italic"
    )
  } else {
    p <- p + ggplot2::geom_rect(
      data = feat,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = 0.9, ymax = 1.1),
      fill = "grey72", colour = "grey40"
    )
    lab <- feat[!is.na(feat$name), , drop = FALSE]
    if (nrow(lab)) {
      lab$x <- (pmax(lab$start, win[1L]) + pmin(lab$end, win[2L])) / 2
      p <- p + ggplot2::geom_text(
        data = lab, ggplot2::aes(x = .data$x, y = 1, label = .data$name),
        size = 2.1, colour = "black", fontface = "bold"
      )
    }
  }

  p +
    ggplot2::geom_point(data = vw, ggplot2::aes(.data$pos, 0.45), colour = "grey65", size = 0.8) +
    ggplot2::geom_point(
      data = snp, ggplot2::aes(.data$pos, 0.45),
      shape = 21, fill = "firebrick", colour = "white", size = 2.6
    ) +
    ggplot2::geom_vline(xintercept = snp$pos, linetype = "dashed", colour = "grey45") +
    ggplot2::coord_cartesian(xlim = win, ylim = c(0.2, 1.25)) +
    ggplot2::labs(
      x = glue::glue("{snp$chrom} position (bp)"), y = NULL,
      title = glue::glue("Doomed target {snp$chrom}:{snp$pos} {snp_label(snp)}"),
      subtitle = glue::glue("no assay passed the gates - {assay$reason}")
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(colour = "#d7301f", face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
}

#' Layered haplotype track: one row per haplotype, its variants and frequency
#'
#' Drawn only when the assay carries a haplotype table.  Each row is a haplotype,
#' labelled with its population frequency and the AS primer it is expected to
#' bind; its own variants are marked (the target in red), so it is visible which
#' haplotype a primer catches and which variant rides with it.
#' @keywords internal
#' @noRd
haplotype_layers <- function(assay, snp, span) {
  h <- assay$haplotypes
  if (is.null(h) || !nrow(h)) {
    return(NULL)
  }
  haps <- unique(h$haplotype)
  h$y <- match(h$haplotype, haps)
  vari <- h[h$pos != snp$pos, , drop = FALSE]
  tgt <- tibble::tibble(x = snp$pos, y = seq_along(haps))
  hits <- assay$haplotype_hits
  lab <- if (is.null(hits)) {
    haps
  } else {
    hits$y <- match(hits$haplotype, haps)
    glue::glue("{hits$haplotype} {round(hits$freq, 2)} -> {hits$primer}")
  }
  base <- tibble::tibble(x = span[1L], xend = span[2L], y = seq_along(haps))
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = base,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      colour = "grey92"
    ) +
    ggplot2::geom_point(
      data = vari, ggplot2::aes(.data$pos, .data$y),
      shape = 21, fill = "#e6ab02", colour = "white", size = 1.8
    ) +
    ggplot2::geom_text(
      data = vari, ggplot2::aes(.data$pos, .data$y, label = .data$alt),
      vjust = -0.9, size = 1.9, family = "mono", colour = "#7f6000"
    ) +
    ggplot2::geom_point(
      data = tgt, ggplot2::aes(.data$x, .data$y),
      shape = 21, fill = "firebrick", colour = "white", size = 2.0
    ) +
    ggplot2::geom_vline(xintercept = snp$pos, linetype = "dashed", colour = "grey65") +
    ggplot2::scale_y_reverse(breaks = seq_along(haps), labels = lab) +
    ggplot2::coord_cartesian(xlim = span) +
    ggplot2::labs(x = NULL, y = NULL, title = "haplotypes (freq -> AS primer)") +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black")
    )
}

#' Local assay map: styled primers, Tm, and the acceptance checks
#'
#' The local view of the assay on the design window (genomic coordinates).  Each
#' primer is drawn as its sequence of letters: neutral black by default, **bold**
#' over the 3' GC-clamp window, and **bright red** on the bases of an SSR or a
#' hairpin stem; a ring marks an injected mismatch / the 3' discriminating base.
#' The Tm sits just past the 3' end, with its deviation from the target as a
#' green/red tumbler.  Below: the amplicon's own sequence, the local variant
#' track and a feature-annotation row (same x axis), then - when `genome` is
#' given - a thin genome-wide ideogram (`assay_ideogram()`) of exact primer hits
#' and in-silico PCR products.  Right: a thin column of per-criterion tumblers.
#' @param assay assembled assay, or a `marker_doomed` object (then only the target
#'   context and the reason are drawn).
#' @param variants variant table.
#' @param snp target variant (one row).
#' @param genome optional `marker_genome` or `DNAStringSet`; when given, the
#'   genome-wide ideogram panel is added.
#' @param features optional annotation to draw on the feature row: a data frame
#'   with `start`/`end` (and optional `strand`, `name`, `chrom`), a path to an
#'   NCBI feature table or GFF3 (read by [marker_features()]), or a gene record
#'   from [ensembl_gene()].  A KASP marker should sit on an annotated feature, so
#'   the row says `"no feature annotation"` when nothing (or no overlap) is given.
#' @return a `patchwork` object: the local map (primers over the variant track,
#'   sharing an x axis) beside the tumbler list from [plot_checks()].  For a
#'   `marker_doomed` target, a single context plot with the reason instead.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci)
#' variants <- tibble::tibble(
#'   chrom = loci$chrom, pos = loci$pos, ref = loci$ref, alt = loci$alt,
#'   consequence = "synonymous_variant"
#' )
#' plot_primers(assay, variants, loci)
#' }
#' @family plots
#' @export
plot_primers <- function(assay, variants, snp, genome = NULL, features = NULL) {
  if (is_kasp_doomed(assay)) {
    return(plot_doomed_target(assay, variants, snp, features))
  }
  win <- design_window(assay, snp)
  aa <- assay$assay[1, ]
  strand <- aa$strand %||% "+"
  plus <- identical(strand, "+")
  d <- primer_features(assay)
  letters <- primer_letters(assay)

  # alignment-level arrow: the amplicon the placement buys, on its own row
  amp_ok <- dplyr::between(aa$amplicon, KASP_LIMITS$amplicon[1], KASP_LIMITS$amplicon[2])
  amp <- tibble::tibble(
    role = factor("amplicon", levels = c(levels(d$role), "amplicon")),
    start = min(d$start), end = max(d$end), forward = plus,
    status = dplyr::if_else(amp_ok, "ok", "broken")
  )
  # one size for the whole map: primers and amplicon share the reference scale
  letter_size <- letter_size_for_span(c(amp$start, amp$end))
  # the amplicon's own sequence (plus strand) as letters on its row
  ti <- function(p) p - snp$pos + assay$template$snp_i
  amp_seq <- as.character(Biostrings::subseq(assay$template$seq, ti(amp$start), ti(amp$end)))
  letters <- dplyr::bind_rows(letters, tibble::tibble(
    role = "amplicon", x = seq.int(amp$start, amp$end),
    letter = strsplit(amp_seq, "")[[1L]], face = "plain", colour = "#292524"
  ))

  # the Tm as a small number just past the 3' end, and its deviation from the
  # target as a green/red tumbler next to it
  lab <- d
  lab$target <- dplyr::if_else(lab$role == "COM", COM_TM_TARGET, AS_TM_TARGET)
  lab$dtm <- lab$tm - lab$target
  lab$meas_x <- dplyr::if_else(lab$forward, lab$end + 1, lab$start - 1)
  lab$hj <- dplyr::if_else(lab$forward, 0, 1)
  lab$dtm_x <- dplyr::if_else(lab$forward, lab$meas_x + 4, lab$meas_x - 4)
  lab$dtm_status <- dplyr::if_else(abs(lab$dtm) <= 3, "ok", "broken")

  # pins: the 3' discriminating base of every oligo, and the injected
  # destabilising mismatch of the AS primers (if the placement carries one)
  pins <- tibble::tibble(role = as.character(d$role), x = d$three_x)
  asx <- d[!is.na(d$stab), ]
  if (nrow(asx)) {
    pins <- dplyr::bind_rows(pins, tibble::tibble(
      role = as.character(asx$role),
      x = dplyr::if_else(asx$forward, asx$end - (asx$stab - 1L), asx$start + (asx$stab - 1L))
    ))
  }

  span <- range(c(d$start, d$end, lab$meas_x, lab$dtm_x, amp$start, amp$end))

  # the annotated features on the same axis: a KASP marker should sit on a gene,
  # not in intergenic space, so if no (overlapping) annotation is given we say so
  feat <- normalise_features(features, span, chrom = snp$chrom)

  p_map <- ggplot2::ggplot() +
    gggenes::geom_gene_arrow(
      data = amp,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$role, forward = .data$forward),
      fill = "grey94", colour = "grey55", linewidth = 0.25
    ) +
    gggenes::geom_gene_arrow(
      data = d,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$role, forward = .data$forward),
      fill = "grey98", colour = "grey75", linewidth = 0.25
    ) +
    letters_layer(letters, size = letter_size) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_point(
      data = pins, ggplot2::aes(x = .data$x, y = .data$role),
      shape = 1, size = 3.4, colour = "grey35"
    ) +
    ggplot2::geom_text(
      data = lab,
      ggplot2::aes(x = .data$meas_x, y = .data$role, label = glue::glue("{round(.data$tm)}\u00b0C"), hjust = .data$hj),
      size = 2.6, colour = "black"
    ) +
    ggplot2::geom_point(
      data = lab, ggplot2::aes(x = .data$dtm_x, y = .data$role, fill = .data$dtm_status),
      shape = 21, size = 2.8, colour = "grey40"
    ) +
    ggplot2::geom_vline(xintercept = snp$pos, linetype = "dashed", colour = "grey45") +
    # feature-annotation row (baseline + boxes, or a note when there is none)
    ggplot2::geom_segment(
      data = tibble::tibble(x = span[1], xend = span[2]),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = "feature", yend = "feature"),
      colour = "grey72", linewidth = 0.3
    )
  if (is.null(feat)) {
    p_map <- p_map + ggplot2::annotate("text",
      x = mean(span), y = "feature", label = "no feature annotation",
      size = 2.1, colour = "grey55", fontface = "italic"
    )
  } else {
    p_map <- p_map +
      gggenes::geom_gene_arrow(
        data = feat,
        ggplot2::aes(xmin = .data$start, xmax = .data$end, y = "feature", forward = .data$forward),
        fill = "grey72", colour = "grey40"
      )
    if (any(!is.na(feat$name))) {
      lab_feat <- feat[!is.na(feat$name), ]
      # keep the label inside the visible window even for a feature that spans it
      lab_feat$x <- (pmax(lab_feat$start, span[1]) + pmin(lab_feat$end, span[2])) / 2
      p_map <- p_map + ggplot2::geom_text(
        data = lab_feat,
        ggplot2::aes(x = .data$x, y = "feature", label = .data$name),
        size = 2.1, colour = "black", fontface = "bold", vjust = 0.5
      )
    }
  }
  p_map <- p_map +
    ggplot2::scale_fill_manual(values = KASP_STATUS, guide = "none") +
    ggplot2::scale_y_discrete(limits = c("feature", "amplicon", "COM", "AS-alt", "AS-ref")) +
    ggplot2::coord_cartesian(xlim = span) +
    ggplot2::labs(x = NULL, y = NULL) +
    gggenes::theme_genes() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black")
    )

  vw <- variants_in_window(variants, snp$chrom, win)
  reg3p <- range(c(aa$win3p_start, aa$win3p_end))
  as_span <- c(aa$as_start, aa$as_end)
  com_span <- c(aa$com_start, aa$com_end)
  in_foot <- dplyr::between(vw$pos, as_span[1], as_span[2]) |
    dplyr::between(vw$pos, com_span[1], com_span[2])
  foot <- vw[in_foot, ]
  # the AS 3'-terminal window and the COM footprint: the two zones a variant
  # must not sit in.  A zone is "clean" when the expected dropout stays within
  # tolerance (DECIPHER), falling back to "no known variant" when the estimate is
  # unavailable.  Same green/red logic for both.
  tol <- KASP_LIMITS$footprint_dropout
  zone_status <- function(drop, n_var) {
    if (length(drop) && !is.na(drop)) {
      if (drop <= tol) "ok" else "broken"
    } else if (isTRUE(n_var == 0L)) {
      "ok"
    } else {
      "broken"
    }
  }
  win_status <- zone_status(aa[["as_dropout"]], aa$as_n_3p)
  com_status <- zone_status(aa[["com_dropout"]], aa$n_var)

  # the local variants: a thin strip on the same x axis.  When the assay carries
  # a linked haplotype, the target + its co-segregating bases are drawn as two
  # parallel tracks (REF vs ALT haplotype), so the linkage is visible.
  lk <- assay$selection$linked
  # when a full haplotype table is present the dedicated layered track carries
  # the linkage, so the two-row REF/ALT view is not drawn twice
  has_link <- !is.null(lk) && nrow(lk) > 0 && is.null(assay$haplotypes)
  if (has_link) {
    gp <- snp$pos - assay$template$snp_i + lk$offset
    hap <- dplyr::bind_rows(
      tibble::tibble(x = snp$pos, y = 1, base = snp$ref, kind = "ref_target"),
      tibble::tibble(x = snp$pos, y = 2, base = snp$alt, kind = "alt_target"),
      tibble::tibble(x = gp, y = 1, base = lk$ref, kind = "ref_linked"),
      tibble::tibble(x = gp, y = 2, base = lk$alt, kind = "alt_linked")
    )
    hap_col <- c(
      ref_target = KASP_STATUS[["neutral"]], alt_target = "#d7301f",
      ref_linked = "grey62", alt_linked = "#e6ab02"
    )
    hap_face <- c(ref_target = "bold", alt_target = "bold", ref_linked = "plain", alt_linked = "bold")
    ybreaks <- c(1, 2)
    ylabs <- c("REF hap", "ALT hap")
    zone <- c(0.5, 2.5)
    tall <- c(0.28, 2.72)
  } else {
    ybreaks <- 1
    ylabs <- "SNP"
    zone <- c(0.6, 1.4)
    tall <- c(0.22, 1.78)
  }
  ys <- ybreaks
  shared <- tibble::tibble(
    pos = rep(vw$pos, times = length(ys)),
    y = rep(ys, each = nrow(vw))
  )
  foot_pts <- tibble::tibble(
    pos = rep(foot$pos, times = length(ys)),
    y = rep(ys, each = nrow(foot))
  )
  p_var <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
      xmin = com_span[1] - 0.5, xmax = com_span[2] + 0.5,
      ymin = zone[1], ymax = zone[2], fill = KASP_STATUS[[com_status]], alpha = 0.22
    ) +
    # the AS 3' window is the discriminating side: drawn as a taller square so it
    # reads as the important zone, coloured by its own dropout
    ggplot2::annotate("rect",
      xmin = reg3p[1] - 0.5, xmax = reg3p[2] + 0.5,
      ymin = tall[1], ymax = tall[2], fill = KASP_STATUS[[win_status]], alpha = 0.30
    ) +
    ggplot2::geom_point(
      data = shared, ggplot2::aes(.data$pos, .data$y),
      colour = "grey70", size = 0.8
    ) +
    ggplot2::geom_point(data = foot_pts, ggplot2::aes(.data$pos, .data$y), shape = 1, colour = "black", size = 2.1)
  if (has_link) {
    p_var <- p_var +
      ggplot2::geom_point(
        data = hap[hap$kind %in% c("ref_target", "alt_target"), ],
        ggplot2::aes(.data$x, .data$y, fill = .data$kind), shape = 21, colour = "white", size = 2.6
      ) +
      ggplot2::geom_text(
        data = hap, ggplot2::aes(.data$x, .data$y, label = .data$base, colour = .data$kind, fontface = .data$kind),
        vjust = -0.9, size = 2.3, family = "mono", show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = hap_col, guide = "none") +
      ggplot2::scale_colour_manual(values = hap_col, guide = "none") +
      ggplot2::scale_discrete_manual("fontface", values = hap_face, guide = "none")
  } else {
    p_var <- p_var +
      variant_points(vw, 1) +
      ggplot2::geom_point(
        data = snp, ggplot2::aes(.data$pos, 1),
        shape = 21, fill = "firebrick", colour = "white", size = 2.4
      ) +
      ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL)
    # a multi-allelic target: name every ALT at the site, so the track reads as
    # one *site* with several alleles, not one allele
    aa <- snp_alleles(snp)
    if (length(aa) > 1L) {
      p_var <- p_var + ggplot2::annotate("text",
        x = snp$pos, y = 1, label = paste(snp$ref[1L], paste(aa, collapse = "/"), sep = ">"),
        vjust = 2.1, size = 2.2, family = "mono", colour = "black"
      )
    }
  }
  p_var <- p_var +
    ggplot2::scale_y_continuous(breaks = ybreaks, labels = ylabs, limits = c(zone[1] - 0.15, zone[2] + 0.15)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(xlim = span) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black")
    )

  # genome-wide ideogram: where the primers hit and where in-silico PCR amplifies
  p_ideo <- assay_ideogram(assay, snp, genome)
  p_hap <- haplotype_layers(assay, snp, span)

  # primer rows are packed tight; the feature row sits at the bottom
  parts <- list(p_map, p_var)
  heights <- c(2.3, 0.4)
  if (!is.null(p_hap)) {
    parts <- c(parts, list(p_hap))
    heights <- c(heights, 0.22 * length(unique(assay$haplotypes$haplotype)))
  }
  if (!is.null(p_ideo)) {
    parts <- c(parts, list(p_ideo))
    heights <- c(heights, 3)
  }
  panels <- patchwork::wrap_plots(parts, ncol = 1, heights = heights)

  # the tumbler needs room for its criterion labels and observed values; a
  # width below ~0.2 clips them into the map panel on the left
  patchwork::wrap_plots(panels, plot_checks(assay), widths = c(1, 0.22)) +
    patchwork::plot_annotation(
      title = "KASP assay: local map + acceptance checks",
      subtitle = glue::glue(
        "target {snp$chrom}:{snp$pos} ({snp_label(snp)}) on strand {strand}"
      )
    )
}

#' Format one acceptance criterion's observed value for the tumbler column
#'
#' Every value is monospaced and at most 4 characters wide, so the column lines
#' up: GC and discrimination as whole percents (`dd%`), Tm and the four dG terms
#' as `dd.d` / `-d.d`, and integer counts plain (compacted to `dK` / `dM` for
#' extreme values).
#' @keywords internal
#' @noRd
format_check_value <- function(name, value) {
  out <- if (!is.numeric(value)) {
    as.character(value)
  } else {
    dplyr::case_when(
      name %in% c("AS GC", "COM GC", "discrim.", "AS drop", "COM drop") ~ sprintf("%.0f%%", 100 * value),
      name %in% c("AS Tm", "COM Tm", "Tm", "dimer", "homodimer", "hairpin", "COM\u2265AS Tm") ~ sprintf("%.1f", value),
      abs(value) >= 1e6 ~ paste0(round(value / 1e6), "M"),
      abs(value) >= 1e3 ~ paste0(round(value / 1e3), "K"),
      TRUE ~ as.character(round(value))
    )
  }
  out[is.na(value)] <- "-"
  out
}

# One tumbler for KASP acceptance checks and for the per-class gate tables.
# `chk` needs `level`, `name`, `status` (ok/broken/neutral) and `value`.
#
# A neutral row (not evaluated, or the check does not apply to this class) shows
# `NA` rather than the prose in `detail`: a sentence such as "not evaluated: no
# population variant table was supplied" overruns the narrow tumbler column and
# hides the neighbouring values.  The full detail stays in the evaluation table.
tumbler_plot <- function(chk) {
  if (!is.null(chk$value)) {
    chk$value <- as.character(chk$value)
    neutral <- !is.na(chk$status) & chk$status == "neutral"
    chk$value[neutral] <- "NA"
  }
  chk$level <- droplevels(factor(chk$level, levels = unique(chk$level)))
  ggplot2::ggplot(chk, ggplot2::aes(x = 0.5, y = factor(.data$name, levels = rev(levels(factor(.data$name)))))) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$status), width = 0.3, height = 0.28) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$value),
      x = 0.92, hjust = 1, size = 2.3, colour = "black", family = "mono"
    ) +
    ggplot2::facet_grid(level ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_manual(values = KASP_STATUS, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0.32, 0.97), breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6.5, colour = "black"),
      strip.text.y = ggplot2::element_blank()
    )
}

#' Acceptance checks as a thin list of tumblers (feature vs alignment)
#'
#' Runs [marker_design_checks()] and shows one named tumbler per criterion - a
#' small green (ok) / red (broken) marker with its observed value and limit - in
#' a single narrow column, faceted by the level judged: the *site* level is the
#' target itself (KASP class, ALT count, linked variants), the *feature* level is
#' one oligo on its own, the *alignment* level is how the primers sit together on
#' the template.  Sits as a thin strip beside [plot_primers()].
#' @param assay assembled assay.
#' @param limits thresholds (see `KASP_LIMITS`).
#' @return a `ggplot` tumbler list.  A row that was not evaluated (or does not
#'   apply) is drawn neutral and shows `NA`; the prose explanation is kept in
#'   `assay$evaluation$detail`, not in the figure.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' assay <- design_hrm_marker(genome, loci)
#' plot_checks(assay)
#' @family plots
#' @export
plot_checks <- function(assay, limits = KASP_LIMITS) {
  chk <- assay$evaluation %||% marker_design_checks(assay, limits)
  if (!"status" %in% names(chk)) {
    chk$status <- dplyr::if_else(chk$pass, "pass", "fail")
  }
  chk <- chk |>
    dplyr::mutate(
      level = factor(.data$level, levels = c("site", "feature", "alignment", "specificity")),
      name = check_short(.data$criterion),
      status = dplyr::case_when(
        .data$status %in% c("unknown", "not_applicable") ~ "neutral",
        .data$status == "pass" ~ "ok",
        TRUE ~ "broken"
      ),
      value = dplyr::coalesce(.data$detail, format_check_value(.data$name, .data$observed))
    )
  # categorical checks draw a short code, never the prose in `detail`: the
  # site class is the whole detail, the Kwok mismatch class its first word
  first_word <- function(x) sub("^(\\S+).*$", "\\1", x)
  cls <- c(biallelic = "bial", degenerate = "degn", complementary = "comp", indel = "indl", none = "--")
  is_site <- chk$criterion == "target site class"
  chk$value[is_site] <- unname(cls[first_word(chk$detail[is_site])])
  sev <- c(severe = "sev", moderate = "mod", minimal = "min", context_dependent = "ctx")
  is_sev <- chk$criterion == "AS 3' mismatch base pair"
  chk$value[is_sev] <- unname(sev[first_word(chk$detail[is_sev])])
  # population facts are chemistry-independent: append ALT frequency and, when
  # available, haplotype coverage to the same check column
  pop <- marker_population_checks(assay)
  if (nrow(pop)) chk <- dplyr::bind_rows(chk, pop)
  tumbler_plot(chk)
}

#' Points of a constant-height arc (bump) spanning `[a, b]` above `y0`
#' @keywords internal
#' @noRd
bump_points <- function(a, b, y0 = 0.7, h = 1, n = 40L) {
  t <- seq(0, pi, length.out = n)
  tibble::tibble(
    x = (a + b) / 2 + (b - a) / 2 * cos(t),
    y = y0 + h * sin(t)
  )
}

#' Hairpins, inverted repeats and SSRs within the assay oligos
#'
#' The sequence-level structure of the three primers side by side: each oligo's
#' bases, its SSR spans underlined in red, and its inverted-repeat / hairpin stems
#' drawn as arcs above the sequence.  The visual twin of the thermodynamic
#' hairpin/SSR screen, so a primer that folds is seen, not merely scored.  Arcs
#' are constant-height bumps (`geom_path()` over computed points) rather than
#' `ggforce::geom_arc()` semicircles: a semicircle's height grows with the stem,
#' which would dominate a compact strip, while the horizontal span already shows
#' the stem length.
#' @param assay assembled assay.
#' @param min_armlength,max_looplength passed to [find_palindromes()].
#' @return a `ggplot` faceted by oligo.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci)
#' plot_oligo_structure(assay)
#' }
#' @family plots
#' @export
plot_oligo_structure <- function(assay, min_armlength = 4L, max_looplength = 6L) {
  p <- assay_primers(assay)
  oligos <- tibble::tibble(oligo = names(p), seq = unname(p), L = nchar(unname(p)))
  pals <- purrr::map(oligos$seq, function(s) {
    find_palindromes(s, min_armlength, max_looplength)[[1L]]
  })
  ssrs <- purrr::map(oligos$seq, repeat_spans)
  lab <- glue::glue(
    "{oligos$oligo}  {purrr::map_int(pals, nrow)} hairpin / {purrr::map_int(ssrs, nrow)} SSR"
  )
  panel_of <- function(oligo) lab[match(oligo, oligos$oligo)]

  # one row of letters per oligo, red where an SSR span sits
  letters <- purrr::pmap(oligos, function(oligo, seq, L) {
    colour <- rep("grey45", L)
    sp <- repeat_spans(seq)
    if (nrow(sp)) {
      idx <- unique(purrr::list_c(purrr::map(seq_len(nrow(sp)), function(j) sp[j, 1L]:sp[j, 2L])))
      colour[idx[idx >= 1L & idx <= L]] <- "#d7191c"
    }
    tibble::tibble(
      panel = panel_of(oligo), oligo = oligo, x = seq_len(L),
      letter = strsplit(seq, "")[[1L]], colour = colour
    )
  }) |>
    purrr::list_rbind()
  letters$panel <- factor(letters$panel, levels = lab)

  ssr <- purrr::imap(ssrs, function(sp, i) {
    if (!nrow(sp)) {
      return(NULL)
    }
    tibble::tibble(panel = lab[i], x = sp[, 1L], xend = sp[, 2L])
  }) |>
    purrr::list_rbind()

  arcs <- purrr::imap(pals, function(h, i) {
    if (!nrow(h)) {
      return(NULL)
    }
    purrr::pmap(h, function(start, end, width) {
      dplyr::mutate(bump_points(start, end), panel = lab[i], arc = paste(lab[i], start, end))
    }) |>
      purrr::list_rbind()
  }) |>
    purrr::list_rbind()

  gg <- ggplot2::ggplot(letters, ggplot2::aes(
    x = .data$x, y = 0, label = .data$letter, colour = .data$colour
  )) +
    ggplot2::geom_text(family = "mono", size = 2.0) +
    ggplot2::scale_colour_identity()
  if (nrow(ssr)) {
    gg <- gg + ggplot2::geom_segment(
      data = ssr,
      ggplot2::aes(x = .data$x - 0.45, xend = .data$xend + 0.45, y = -0.45, yend = -0.45),
      colour = "#d7191c", linewidth = 1.5
    )
  }
  if (nrow(arcs)) {
    gg <- gg + ggplot2::geom_path(
      data = arcs,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$arc),
      colour = "#4477aa", alpha = 0.75, linewidth = 0.4
    )
  }
  gg +
    ggplot2::facet_grid(panel ~ ., scales = "free_x") +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(-0.9, 2.2), breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0, size = 6.5)
    )
}

#' Self-complementarity dot-plot of one oligo
#'
#' Plots every position pair that can Watson-Crick pair inside the same oligo
#' (`i` pairs `j` when base `i` is the complement of base `j`).  A hairpin shows
#' as an off-diagonal block; the closer the block sits to the main diagonal, the
#' tighter its stem/loop.  Sequence-level, so it complements the thermodynamic
#' hairpin free energy with the raw pairing pattern.
#' @param seq one oligo, 5'->3'.
#' @return a `ggplot` dot-plot.
#' @examples
#' plot_self_complementarity("GAAGGTGACCAAGTTCATGCT")
#' @family plots
#' @export
plot_self_complementarity <- function(seq) {
  seq <- as.character(seq)[1L]
  b <- strsplit(seq, "")[[1L]]
  # i pairs j when base i is the complement of base j
  good <- outer(b, chartr("ACGT", "TGCA", b), "==")
  hit <- which(good, arr.ind = TRUE)
  n <- length(b)
  ggplot2::ggplot(tibble::tibble(x = hit[, 1L], y = hit[, 2L]), ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(width = 0.9, height = 0.9, fill = "#4477aa") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.2) +
    ggplot2::coord_equal(xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5)) +
    ggplot2::labs(x = "position", y = "position") +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Gene model with the variant track and a zoom on the design window
#'
#' The context of the target: the gene model and the 1001-Genomes variant track,
#' with `ggforce::facet_zoom()` adding a zoomed panel on the design window and
#' boxing the zoomed area on the full view - the "where is this SNP" picture that
#' sits above the local [plot_primers()] map.
#' @param gene gene record.
#' @param variants variant table.
#' @param snp target variant (one row).
#' @param assay assembled assay.
#' @return a `ggplot` object: the gene model and variant track, faceted by
#'   `ggforce::facet_zoom()` into the full view and a zoom on the design window.
#' @examples
#' genome <- read_genome_fasta(system.file("extdata", "example_region.fa", package = "tidymarkers"))
#' loci <- utils::read.csv(
#'   system.file("extdata", "example_region_targets.csv", package = "tidymarkers"),
#'   colClasses = "character"
#' )[, c("chrom", "pos", "ref", "alt")]
#' loci$pos <- as.integer(loci$pos)
#' \dontrun{
#' assay <- design_asq_marker(genome, loci)
#' # a gene record over this slice; ensembl_gene("AT4G05320") fetches the real one
#' # over the network
#' gene <- list(symbol = "AT4G05320", chrom = "Chr4", start = 2400L, end = 3200L, strand = "+")
#' variants <- tibble::tibble(
#'   chrom = loci$chrom, pos = loci$pos, ref = loci$ref, alt = loci$alt,
#'   consequence = "synonymous_variant"
#' )
#' plot_variant_landscape(gene, variants, loci, assay)
#' }
#' @family plots
#' @export
plot_variant_landscape <- function(gene, variants, snp, assay) {
  gene_df <- tibble::tibble(
    molecule = gene$symbol, start = gene$start, end = gene$end,
    strand = gene$strand
  )
  v <- variants |>
    dplyr::mutate(consequence = lump_top_n(.data$consequence, 6L))
  # label the target from the variant itself - a bare (chrom, pos, ref, alt)
  # record carries no `id`, so the plot must not depend on one; a multi-allelic
  # target shows every ALT
  extra_alt <- if ("alt2" %in% names(snp) && !is.na(snp$alt2[1L])) paste0("/", snp$alt2) else ""
  snp <- snp |>
    dplyr::mutate(label = glue::glue("{.data$chrom}:{.data$pos} {.data$ref}>{.data$alt}{extra_alt}"))
  win <- design_window(assay, snp)

  ggplot2::ggplot() +
    gggenes::geom_gene_arrow(
      data = gene_df,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = 1, forward = .data$strand > 0),
      fill = "grey80"
    ) +
    variant_points(v, 2) +
    ggplot2::geom_point(
      data = snp, ggplot2::aes(.data$pos, 2),
      shape = 21, fill = "firebrick", colour = "white", size = 3
    ) +
    ggrepel::geom_text_repel(
      data = snp, ggplot2::aes(.data$pos, 2, label = .data$label),
      size = 2.6, nudge_y = 0.6, min.segment.length = 0
    ) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "consequence") +
    ggplot2::scale_y_continuous(
      limits = c(0.5, 3), breaks = c(1, 2),
      labels = c(gene$symbol, "variants")
    ) +
    ggplot2::labs(
      x = glue::glue("{gene$chrom} position (bp)"), y = NULL,
      title = glue::glue("Variation across {gene$symbol} (1001 Genomes), zoom on the design window")
    ) +
    ggforce::facet_zoom(xlim = win, zoom.size = 1) +
    ggplot2::theme_minimal()
}
