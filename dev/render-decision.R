# =============================================================================
#  dev/render-decision.R  --  (re)draw the README decision-tree image
# =============================================================================
#  The trees are generated from the schema in R/marker_decision.R, which itself
#  reads the engine's constants (KASP_GATES_*, KASP_RELAX_*, KASP_TIERS_*) and
#  the chemistry registry (marker_schemes()).  One image *and one text double*
#  per marker chemistry are written to the git-ignored smoke dir, and the
#  README-facing SNP and multiplex images are committed.  Run this after touching the selection policy or the
#  scheme registry and the diagrams stay in step with the code - nothing is drawn
#  by hand.
#
#  Vector first: the committed figures are .svg (svglite, via ggraph), so text
#  stays selectable and searchable and the image scales without banding.  A
#  .png twin is written next to each .svg for viewers without SVG support; the
#  .svg is the source of truth and the one the README embeds.
#
#  Run:  Rscript dev/render-decision.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

if (!requireNamespace("ggraph", quietly = TRUE)) {
  stop("ggraph is needed to draw the decision tree")
}

fs::dir_create("man/figures")
fs::dir_create("results/plot-smoke")

# the layouter, written with an explicit svglite device: `ggsave(device =
# <a device function>)` opens one SVG page per plot in a single file, so a loop
# over the schemes produced N plots sharing the first page's size.  The page is
# sized from the solved layout, not from a constant, so a scheme with more
# checks gets the space it needs and the committed figure keeps its aspect.
render_svg_plot <- function(file, x = NULL, scheme = NULL) {
  p <- if (is.null(x)) plot_decision_tree(scheme = scheme) else plot_decision_tree(x)
  g <- if (is.null(x)) marker_decision_graph(scheme = scheme) else x
  bb <- decision_bbox(decision_layout(if (inherits(g, "marker_assay")) g$decision %||% g else g))
  # headroom for the title, subtitle and legend the plot draws outside the data
  svglite::svglite(file, width = bb$xmax - bb$xmin + 1.2, height = bb$ymax - bb$ymin + 1.9, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(p)
  invisible(file)
}

# rsvg is optional: without it only the .svg ships (a missing device must not
# abort the render, the .svg is the source of truth).
render_png_twin <- function(svg_file, width = 1600) {
  png <- sub("[.]svg$", ".png", svg_file)
  if (!requireNamespace("rsvg", quietly = TRUE)) {
    if (file.exists(png)) unlink(png)
    return(invisible(NULL))
  }
  rsvg::rsvg_png(svg_file, width = width)
  invisible(png)
}

# one schema per marker chemistry: the allele-specific core is shared, the
# marker-agnostic primer subset is the bare "primer" scheme, and a duplex
# chemistry adds its compatibility stage.  The README embeds the SNP one (the
# package's core) and the multiplex one.
for (sc in marker_schemes()$scheme) {
  svg <- file.path("results/plot-smoke", paste0("decision-", sc, ".svg"))
  render_svg_plot(svg, scheme = sc)
  render_png_twin(svg)
  # the text double: same schema, greppable, next to every image
  writeLines(
    marker_decision_text(sc),
    file.path("results/plot-smoke", paste0("decision-", sc, ".txt"))
  )
  cat("wrote results/plot-smoke/decision-", sc, ".svg + .txt\n", sep = "")
}
# the committed, README-facing images (kept in step with the code by this script)
render_svg_plot("man/figures/decision-tree.svg", scheme = "snp")
cat("wrote man/figures/decision-tree.svg\n")
render_svg_plot("man/figures/decision-multiplex.svg", scheme = "duplex")
cat("wrote man/figures/decision-multiplex.svg\n")

# an annotated example (a real assay's traversal) for eyeballing, git-ignored.
# Skipped with TIDYKASP_SKIP_EXAMPLE=1 (the auto-render in dev/check.R), so the
# gate never needs the demo data.
if (!nzchar(Sys.getenv("TIDYKASP_SKIP_EXAMPLE"))) {
  tryCatch(
    {
      genome <- ath_genome()
      variants <- ensembl_variants("Chr4", 268901, 271503)
      snp <- tibble::tibble(chrom = "Chr4", pos = 269666L, ref = "G", alt = "A")
      assay <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = variants, flank = 150L))
      render_svg_plot("results/plot-smoke/decision-fri.svg", x = assay)
      writeLines(
        marker_decision_text(assay$scheme %||% "snp"),
        "results/plot-smoke/decision-fri.txt"
      )
      cat("wrote results/plot-smoke/decision-fri.svg + .txt\n")
    },
    error = function(e) cat("annotated example skipped:", conditionMessage(e), "\n")
  )
}
