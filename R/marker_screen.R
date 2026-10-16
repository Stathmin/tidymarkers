# =============================================================================
#  marker_screen.R  --  a cheap, vectorised "promising / doomed" target screen
# =============================================================================
#  Answers, without any thermodynamics, BLAST or genome-wide scan, "can a KASP
#  design exist at this SNP at all?" - so whole genes can be filtered before the
#  expensive steps (pick_target_snp()/assemble_kasp_assay()/rBLAST).
#
#  For every target SNP it counts how many of the canonical AS placements
#  (strand x length) pass the *compositional* tests that make a primer hopeless:
#  a GC desert, a homopolymer run, an SSR, a variant hotzone, or another variant
#  sitting in the 3' window.  The per-chromosome work (GC prefix sums, run/SSR
#  masks, one IRanges overlap pass) is done once; the per-SNP loop is trivial.
#
#  A linked variant (co-segregating with the target) is *not* a sin: pass it via
#  `linked` and the 3' window forgives it (see ROADMAP, "linked variants").
# =============================================================================

#' Cheap, vectorised screen of candidate target SNPs
#'
#' For each target SNP, counts how many of the canonical AS placements
#' (`strands x lengths`) clear the cheap compositional tests - GC within
#' `limits$gc`, no homopolymer run, no SSR, no variant hotzone, and a clean 3'
#' window - and returns a `promising`/`doomed` verdict with a `k/n` score.  No
#' thermodynamics, no BLAST, no genome-wide scan: milliseconds per gene, so a
#' gene scoring `0/n` can be skipped before the expensive steps.
#'
#' @param variants variant table (`chrom`, `pos`, `ref`, `alt`; or the Ensembl
#'   shape with `ref_genome`, from which `ref`/`alt` are derived).
#' @param genome a `marker_genome` or `DNAStringSet`.
#' @param lengths AS-primer lengths to try (each `strand` x `length` is one
#'   placement).
#' @param strands strands to try.
#' @param clean_3p size of the 3'-terminal window that must be variant-free.
#' @param hotzone most *other* variants tolerated under a placement before it is
#'   a hotzone.
#' @param linked optional two-column table (`target`, `linked`) of positions
#'   known to co-segregate with a target; a variant there is forgiven in the 3'
#'   window.
#' @param limits thresholds (see `KASP_LIMITS`).
#' @return one row per target SNP: `id`, `chrom`, `pos`, `ref`, `alt`, `class`
#'   (see [marker_site_class()]), `n_tried`, `n_ok`, `score` (`"7/8"`), `verdict`
#'   and `reasons`.
#' @examples
#' \dontrun{
#' screen_snp_targets(ensembl_variants("Chr4", 268901, 271503), ath_genome())
#' }
#' @family KASP design
#' @export
screen_snp_targets <- function(variants, genome, lengths = 20:27, strands = c("+", "-"),
                               clean_3p = 5L, hotzone = 8L, linked = NULL,
                               limits = KASP_LIMITS) {
  v <- tibble::as_tibble(variants)
  if (!all(c("chrom", "pos") %in% names(v))) {
    cli::cli_abort(gettext("`variants` needs `chrom` and `pos` columns."))
  }
  chroms <- unique(v$chrom)
  purrr::map(chroms, function(ch) {
    screen_one_chrom(
      v[v$chrom == ch, , drop = FALSE], ch, genome,
      lengths, strands, clean_3p, hotzone, linked, limits
    )
  }) |>
    purrr::list_rbind()
}

#' Screen the target SNPs of one chromosome
#' @keywords internal
#' @noRd
screen_one_chrom <- function(g, chrom, genome, lengths, strands, clean_3p, hotzone, linked, limits) {
  lmax <- max(lengths)
  lo <- max(1L, min(g$pos) - lmax - clean_3p)
  hi <- max(g$pos) + lmax + clean_3p
  s <- as.character(genome_seq(genome, chrom, lo, hi))
  n <- nchar(s)
  chars <- strsplit(s, "")[[1L]]

  # --- one pass over the sequence: GC prefix sums + run/SSR mask ---------------
  gc <- chars %in% c("G", "C")
  cum_gc <- c(0L, cumsum(as.integer(gc)))
  banned <- rep(FALSE, n)
  pats <- c(
    repeat_regex(1L, max(limits$max_run, 2L)), # homopolymers
    vapply(2:4, repeat_regex, character(1), min_copies = 3L) # SSRs
  )
  for (pat in pats) {
    m <- stringr::str_locate_all(s, pat)[[1L]]
    if (nrow(m)) {
      for (k in seq_len(nrow(m))) banned[m[k, 1L]:m[k, 2L]] <- TRUE
    }
  }
  cum_ban <- c(0L, cumsum(as.integer(banned)))

  # --- derive the allele set / REF for the targets ----------------------------
  alleles <- if ("ref_genome" %in% names(g)) g$ref_genome else NULL
  ref <- if ("ref" %in% names(g)) as.character(g$ref) else substring(s, g$pos - lo + 1L, g$pos - lo + 1L)
  alt_sets <- if (!is.null(alleles)) {
    lapply(seq_along(ref), function(i) {
      a <- setdiff(as.character(alleles[[i]]), ref[i])
      a[nzchar(a) & grepl("^[ACGTN]+$", a)]
    })
  } else if ("alt" %in% names(g)) {
    lapply(as.character(g$alt), function(a) a[nzchar(a) & grepl("^[ACGTN]+$", a)])
  } else {
    rep(list(character()), length(ref))
  }
  # a representative ALT for display: the first SNP if there is one
  alt <- vapply(alt_sets, function(a) {
    s <- a[nchar(a) == 1L & grepl("^[ACGT]$", a)]
    if (length(s)) s[[1L]] else if (length(a)) a[[1L]] else NA_character_
  }, character(1))
  # per-site KASP class (biallelic / degenerate / complementary / indel); a
  # multi-allelic site is designable when one degenerate ALT primer covers it
  cls <- purrr::map2(ref, alt_sets, marker_site_class)
  is_snp <- nchar(ref) == 1L & nchar(alt) == 1L
  disc <- is_snp & purrr::map_lgl(cls, "disc")

  # every variant on the chromosome is a subject for the variant-density tests
  vs_ir <- IRanges::IRanges(g$pos, g$pos)
  linked_ir <- function(p) {
    if (is.null(linked)) {
      return(NULL)
    }
    IRanges::IRanges(linked$linked[linked$target == p], width = 1L)
  }

  grid <- tidyr::crossing(strand = strands, length = lengths) # canonical placements
  n_tried <- nrow(grid)

  rows <- purrr::map(which(is_snp), function(i) {
    p <- g$pos[i]
    lt <- dplyr::if_else(grid$strand == "+", p - grid$length + 1L, p)
    rt <- dplyr::if_else(grid$strand == "+", p, p + grid$length - 1L)
    ia <- lt - lo + 1L
    ib <- rt - lo + 1L
    ok_range <- ia >= 1L & ib <= n
    ia_c <- pmax(ia, 1L)
    ib_c <- pmin(ib, n)
    gc_frac <- (cum_gc[ib_c + 1L] - cum_gc[ia_c]) / (ib_c - ia_c + 1L)
    gc_ok <- ok_range & gc_frac >= limits$gc[1L] & gc_frac <= limits$gc[2L]
    ban_ok <- ok_range & (cum_ban[ib_c + 1L] - cum_ban[ia_c]) == 0L
    # 3'-terminal window; the target's own base is excluded (it is the anchor)
    w3 <- IRanges::IRanges(
      dplyr::if_else(grid$strand == "+", p - clean_3p + 1L, p),
      dplyr::if_else(grid$strand == "+", p, p + clean_3p - 1L)
    )
    n3 <- IRanges::countOverlaps(w3, vs_ir) - 1L
    link <- linked_ir(p)
    if (!is.null(link)) n3 <- n3 - IRanges::countOverlaps(w3, link)
    var_ok <- n3 <= 0L
    hot <- IRanges::countOverlaps(IRanges::IRanges(lt, rt), vs_ir) - 1L
    hot_ok <- hot <= hotzone
    ok <- gc_ok & ban_ok & var_ok & hot_ok
    tibble::tibble(
      n_ok = sum(ok),
      reasons = paste(
        unique(c(
          if (any(!gc_ok)) "gc",
          if (any(!ban_ok)) "run/ssr",
          if (any(!var_ok)) "var_3p",
          if (any(!hot_ok)) "hotzone"
        )),
        collapse = ","
      )
    )
  }) |>
    purrr::list_rbind()

  if (nrow(rows) == 0L) {
    return(tibble::tibble())
  }
  out <- tibble::tibble(
    id = if ("id" %in% names(g)) as.character(g$id[which(is_snp)]) else NA_character_,
    chrom = chrom, pos = g$pos[which(is_snp)], ref = ref[which(is_snp)], alt = alt[which(is_snp)],
    class = purrr::map_chr(cls, "class")[which(is_snp)],
    n_tried = n_tried, n_ok = rows$n_ok,
    score = paste0(rows$n_ok, "/", n_tried),
    verdict = dplyr::if_else(rows$n_ok > 0L, "promising", "doomed"),
    reasons = rows$reasons
  )
  # a discriminable SNP has real placements; a complementary one cannot be primed
  out$n_tried[!disc[which(is_snp)]] <- 0L
  out$n_ok[!disc[which(is_snp)]] <- 0L
  out$score[!disc[which(is_snp)]] <- "0/0"
  out$verdict[!disc[which(is_snp)]] <- "doomed"
  out$reasons[!disc[which(is_snp)]] <- "not_discriminable"
  out
}

#' Plot a target screen along a gene
#'
#' A compact per-variant view of [screen_snp_targets()]: the `n_ok` placements
#' per SNP, coloured promising (green) / doomed (red), over the gene model.
#' @param screen result of [screen_snp_targets()].
#' @param gene optional gene record (drawn as a track).
#' @param threshold draw a line at this `n_ok` (default 1).
#' @return a `ggplot`.
#' @examples
#' \dontrun{
#' plot_snp_screen(screen_snp_targets(variants, ath_genome()), ensembl_gene("AT4G00650"))
#' }
#' @family KASP plots
#' @export
plot_snp_screen <- function(screen, gene = NULL, threshold = 1L) {
  if (nrow(screen) == 0L) cli::cli_abort(gettext("`screen` is empty."))
  p <- ggplot2::ggplot(screen, ggplot2::aes(.data$pos, .data$n_ok, colour = .data$verdict)) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$pos, yend = 0), linewidth = 0.3, alpha = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", colour = "grey55") +
    ggplot2::scale_colour_manual(values = c(promising = "#1b7837", doomed = "#d7301f"), name = NULL) +
    ggplot2::labs(
      x = glue::glue("{unique(screen$chrom)} position (bp)"),
      y = "AS placements that pass (n_ok)",
      title = "Target screen: is a KASP design possible here?"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")
  if (!is.null(gene)) {
    p <- p + ggplot2::annotate("rect",
      xmin = gene$start, xmax = gene$end, ymin = -0.6, ymax = -0.2,
      fill = "grey75", colour = "grey45"
    ) +
      ggplot2::annotate("text",
        x = (gene$start + gene$end) / 2, y = -0.4, label = gene$symbol,
        size = 2.2, colour = "black", fontface = "bold"
      )
  }
  p
}
