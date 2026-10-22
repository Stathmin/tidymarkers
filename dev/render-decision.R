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

use_dot <- nzchar(Sys.which("dot"))
if (!use_dot) message("graphviz 'dot' not found; falling back to ggraph")

# Graphviz dot is the routing solver: one cluster per stage, edges routed
# around nodes instead of long straight lines to DOOMED.  Dot emits the .dot
# and a -Tsvg directly; without dot the ggplot layouter draws the same schema
# through svglite.  The committed figures use the ggplot layouter: it is the
# only one that can pin stage N and N+1 side by side (dot's ranking re-flows
# the columns, see `decision_dot()`).
render_svg <- function(file, dot_file, scheme = NULL, assay = NULL) {
  if (use_dot) {
    sc <- if (is.null(assay)) scheme else (assay$scheme %||% "snp")
    writeLines(decision_dot(scheme = sc, assay = assay), dot_file)
    status <- system2("dot", c("-Tsvg", "-o", shQuote(file), shQuote(dot_file)))
    if (!identical(as.integer(status), 0L)) stop("dot failed for ", dot_file)
  } else {
    p <- if (is.null(assay)) plot_decision_tree(scheme = scheme) else plot_decision_tree(assay)
    ggplot2::ggsave(file, p, width = 18, height = 14, bg = "white", device = svglite::svglite)
  }
  invisible(file)
}

# the ggplot layouter, written with an explicit svglite device: `ggsave(device =
# <a device function>)` opens one SVG page per plot in a single file, so a loop
# over the schemes produced N plots sharing the first page's size.
render_svg_plot <- function(file, x = NULL, scheme = NULL, width = 11, height = 8) {
  p <- if (is.null(x)) plot_decision_tree(scheme = scheme) else plot_decision_tree(x)
  svglite::svglite(file, width = width, height = height, bg = "white")
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
