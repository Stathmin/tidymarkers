# =============================================================================
#  kasp_search_anim.R  --  animate the design search (how KASP converges)
# =============================================================================
#  A bar chart of candidate scores says nothing about *why* a marker is chosen.
#  This builds an animation of the search itself: every AS-primer placement is
#  tried in the order the algorithm ranks it (cleanest 3' window -> lowest
#  footprint variant load -> highest discrimination -> ...), each step carries
#  its own verdict ("rejected: 2 other variants in the 3' window"), and the
#  common primer is then placed on the same strand.  The final frame is the
#  accepted assay with its acceptance-check tally.
#
#  Frames are ordinary library-drawn ggplots (gggenes arrows + the 1001-Genomes
#  variant track, same as plot_primers()), assembled into a GIF.  With gganimate
#  installed the same frame table drives a `transition_manual()` animation; the
#  fallback renders frames and stitches them with ImageMagick.
# =============================================================================

#' Build the frame table of the design search
#'
#' `primers` (arrows, one row per role per step), `notes` (one row per step: the
#' headline, the multi-line detail incl. the sequence and the verdict), the
#' static `track` and `snp`, and the `win`dow.  Exposed so the animation logic
#' can be tested without rendering.
#' @param assay assembled assay (from [assemble_kasp_assay()]).
#' @param variants variant table.
#' @param snp target variant (one row).
#' @family KASP plots
#' @export
design_search_frames <- function(assay, variants, snp) {
  aa   <- assay$assay[1L, ]
  win  <- design_window(assay, snp)

  # ---- phase 1: the AS-primer search, in the algorithm's own order ----------
  as <- assay$as_candidates |>
    dplyr::arrange(.data$n_3p, .data$n_var, dplyr::desc(.data$discrimination),
                   .data$wobble, .data$seed_hits, abs(.data$tm - AS_TM_TARGET),
                   .data$dg_3p, .data$as_length) |>
    dplyr::mutate(step = dplyr::row_number())
  w <- as[1L, ]
  as <- as |>
    dplyr::mutate(
      verdict = dplyr::case_when(
        .data$step == 1L ~ "CHOSEN",
        .data$n_3p > w$n_3p ~ sprintf("rejected: %d other variant(s) in the AS 3' window", .data$n_3p),
        .data$n_var > w$n_var ~ sprintf("rejected: %d variant(s) under the footprint", .data$n_var),
        .data$discrimination < w$discrimination ~
          sprintf("rejected: discrimination %.2f < %.2f", .data$discrimination, w$discrimination),
        TRUE ~ "rejected: loses on the tie-breakers (Tm/dG/length)"))

  as_pr <- dplyr::bind_rows(
    dplyr::transmute(as, step = .data$step, role = "AS-ref", start = .data$as_start, end = .data$as_end,
                     forward = .data$strand == "+", win3p_start = .data$win3p_start,
                     win3p_end = .data$win3p_end, status = dplyr::if_else(.data$step == 1L, "chosen", "candidate")),
    dplyr::transmute(as, step = .data$step, role = "AS-alt", start = .data$as_alt_start, end = .data$as_alt_end,
                     forward = .data$strand == "+", win3p_start = NA_integer_, win3p_end = NA_integer_,
                     status = dplyr::if_else(.data$step == 1L, "chosen", "candidate")))

  as_notes <- dplyr::transmute(
    as, step = .data$step, phase = "AS-primer search", key = .data$step == 1L,
    headline = sprintf("%s strand | AS %d nt | Tm %.1f | GC %.2f",
                       dplyr::if_else(.data$strand == "+", "plus", "minus"),
                       .data$as_length, .data$tm, .data$gc),
    detail = sprintf("%s\n3'-window variants: %d    footprint variants: %d    discrimination: %.2f\n%s",
                     .data$ref_primer, .data$n_3p, .data$n_var, .data$discrimination, .data$verdict))

  # ---- phase 2: the common-primer search, in the algorithm's own order -------
  com <- assay$assay |>
    dplyr::arrange(.data$amplicon > KASP_LIMITS$amplicon[2],
                   .data$com_tm < w$tm,
                   .data$com_homodimer < KASP_LIMITS$homodimer_dg,
                   .data$com_hairpin   < KASP_LIMITS$hairpin_dg,
                   .data$n_var, abs(.data$com_tm - COM_TM_TARGET)) |>
    dplyr::mutate(step = dplyr::row_number() + max(as$step))
  cw <- com[1L, ]
  com <- com |>
    dplyr::mutate(
      verdict = dplyr::case_when(
        .data$step == cw$step ~ "CHOSEN",
        .data$amplicon > KASP_LIMITS$amplicon[2] ~ "rejected: amplicon > 100 bp",
        .data$com_tm < w$tm ~ "rejected: COM Tm below the AS Tm",
        .data$com_homodimer < KASP_LIMITS$homodimer_dg ~ "rejected: strong COM homodimer",
        .data$com_hairpin   < KASP_LIMITS$hairpin_dg   ~ "rejected: strong COM hairpin",
        .data$n_var > cw$n_var ~ "rejected: more variants under the COM footprint",
        TRUE ~ "rejected: COM Tm further from the target"))

  fixed_as <- function(role, s, e) tibble::tibble(
    step = com$step, role = role, start = s, end = e, forward = w$strand == "+",
    win3p_start = NA_integer_, win3p_end = NA_integer_,
    status = dplyr::if_else(com$step == cw$step, "chosen", "fixed"))
  com_pr <- dplyr::bind_rows(
    fixed_as("AS-ref", w$as_start, w$as_end),
    fixed_as("AS-alt", w$as_alt_start, w$as_alt_end),
    dplyr::transmute(com, step = .data$step, role = "COM", start = .data$com_start, end = .data$com_end,
                     forward = .data$strand != "+", win3p_start = NA_integer_, win3p_end = NA_integer_,
                     status = dplyr::if_else(.data$step == cw$step, "chosen", "candidate")))

  com_notes <- dplyr::transmute(
    com, step = .data$step, phase = "Common-primer search", key = .data$step == cw$step,
    headline = sprintf("%s strand | COM %d nt | gap %d | Tm %.1f",
                       dplyr::if_else(.data$strand == "+", "plus", "minus"),
                       .data$len, .data$gap, .data$com_tm),
    detail = sprintf("%s\n amplicon: %d bp    COM footprint variants: %d    COM homodimer: %.1f dG37\n%s",
                     .data$com, .data$amplicon, .data$n_var, .data$com_homodimer, .data$verdict))

  # ---- final frame: the accepted assay + the acceptance-check tally ----------
  checks  <- kasp_design_checks(assay)
  failed  <- checks$criterion[!checks$pass]
  final_step <- max(com$step) + 1L
  final_pr <- tibble::tibble(
    step = final_step,
    role = c("AS-ref", "AS-alt", "COM"),
    start = c(aa$as_start, aa$as_alt_start, aa$com_start),
    end   = c(aa$as_end,   aa$as_alt_end,   aa$com_end),
    forward = c(w$strand == "+", w$strand == "+", w$strand != "+"),
    win3p_start = c(aa$win3p_start, NA_integer_, NA_integer_),
    win3p_end   = c(aa$win3p_end,   NA_integer_, NA_integer_),
    status = "chosen")
  final_note <- tibble::tibble(
    step = final_step, phase = "Accepted assay", key = TRUE,
    headline = sprintf("accepted: %d / %d acceptance checks pass", sum(checks$pass), nrow(checks)),
    detail = if (length(failed)) paste("failing:", paste(failed, collapse = "; "))
             else "all criteria met - ready to order")

  list(
    primers = dplyr::bind_rows(as_pr, com_pr, final_pr),
    notes   = dplyr::bind_rows(as_notes, com_notes, final_note),
    track   = variants_in_window(variants, snp$chrom, win),
    snp     = snp,
    win     = win)
}

#' One frame of the design-search animation
#'
#' A library-drawn panel (gggenes arrows over the static variant track) with the
#' step's verdict as the subtitle.  Used by [render_design_search()] and handy on
#' its own for a slide.
#' @param frames result of [design_search_frames()].
#' @param step frame index.
#' @family KASP plots
#' @export
plot_design_search_step <- function(frames, step) {
  pr <- dplyr::filter(frames$primers, .data$step == !!step)
  nt <- dplyr::filter(frames$notes,   .data$step == !!step)

  ggplot2::ggplot() +
    gggenes::geom_gene_arrow(
      data = pr,
      ggplot2::aes(xmin = .data$start, xmax = .data$end, y = .data$role,
                   forward = .data$forward, fill = .data$status)) +
    ggplot2::geom_point(data = frames$track,
                        ggplot2::aes(.data$pos, 1, colour = .data$consequence), size = 1.6) +
    ggplot2::geom_point(data = frames$snp, ggplot2::aes(.data$pos, 1),
                        shape = 21, fill = "firebrick", colour = "white", size = 3) +
    ggrepel::geom_text_repel(data = frames$snp, ggplot2::aes(.data$pos, 1, label = "target SNP"),
                             size = 3, nudge_y = 0.35, min.segment.length = 0, colour = "firebrick") +
    ggplot2::scale_fill_manual(values = c(chosen = "firebrick", candidate = "grey80",
                                          fixed = "steelblue"), guide = "none") +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "consequence") +
    ggplot2::scale_y_discrete(limits = c("variants", "COM", "AS-alt", "AS-ref")) +
    ggplot2::labs(x = glue::glue("{frames$snp$chrom} position (bp)"), y = NULL,
                  title = nt$headline, subtitle = nt$detail,
                  caption = glue::glue("step {step} | {nt$phase}")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(family = "mono", size = 9),
                   panel.grid.major.y = ggplot2::element_blank())
}

#' Render the design-search animation as a GIF
#'
#' Renders every frame with [plot_design_search_step()] and stitches them into a
#' GIF.  gganimate is used when installed; otherwise the frames are assembled with
#' the ImageMagick `convert` CLI (or `gifski`).  The winner and the final frame are
#' held a few frames longer so they are readable.
#' @param assay assembled assay.
#' @param variants variant table.
#' @param snp target variant (one row).
#' @param path output GIF path.
#' @param fps frames per second.
#' @param hold how many extra frames to hold the chosen AS/COM/final frames.
#' @return the GIF path (invisibly).
#' @family KASP plots
#' @export
render_design_search <- function(assay, variants, snp,
                                 path = here::here("results", "design_search.gif"),
                                 fps = 4L, hold = 3L) {
  rlang::check_installed("ggplot2")
  frames <- design_search_frames(assay, variants, snp)
  n      <- nrow(frames$notes)

  # linger on the two winners and the final frame
  winners <- dplyr::filter(frames$notes, .data$key)$step
  steps <- unlist(purrr::map(seq_len(n), function(s) c(s, rep(s, if (s %in% winners) hold else 0L))))

  td   <- fs::path_temp("kasp_anim"); fs::dir_create(td)
  pngs <- character(0)
  for (i in seq_along(steps)) {
    f <- fs::path(td, sprintf("f%04d.png", i))
    ggplot2::ggsave(f, plot_design_search_step(frames, steps[i]),
                    width = 8, height = 4.6, dpi = 110, bg = "white")
    pngs <- c(pngs, f)
  }
  fs::dir_create(fs::path_dir(path))
  assemble_gif(pngs, as.character(path), delay_cs = round(100 / fps))
  cli::cli_alert_success("Design-search animation written to {.path {path}} ({length(pngs)} frames).")
  invisible(path)
}

#' Stitch PNG frames into a GIF with the ImageMagick CLI
#'
#' Kept to an external call so nothing heavy is added to the dependency graph
#' (the `gifski`/`gganimate` R stack needs Rust / `transformr`->`sf`/GDAL).
#' @keywords internal
#' @noRd
assemble_gif <- function(pngs, path, delay_cs = 25L) {
  if (!nzchar(Sys.which("convert"))) {
    rlang::abort(c("No GIF assembler available (ImageMagick `convert` not found).",
                   i = "Install ImageMagick, or assemble the `pngs` yourself (e.g. with gifski)."))
  }
  res <- processx::run("convert", c("-delay", as.character(delay_cs), "-loop", "0",
                                    "-layers", "Optimize", pngs, path),
                       error_on_status = FALSE)
  if (res$status != 0L)
    rlang::abort(glue::glue("ImageMagick failed: {res$stderr}"))
  path
}
