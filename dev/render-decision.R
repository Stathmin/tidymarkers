# =============================================================================
#  dev/render-decision.R  --  (re)draw the README decision-tree images
# =============================================================================
#  The trees are generated from the schema in R/marker_decision.R, which itself
#  reads the engine's constants (KASP_GATES_*, KASP_RELAX_*, KASP_TIERS_*) and
#  the chemistry registry (marker_schemes()).  One image *and one text double*
#  per marker chemistry are written to the git-ignored smoke dir, and the
#  README-facing SNP and multiplex images are committed.  Run this after
#  touching the selection policy or the scheme registry and the diagrams stay in
#  step with the code - nothing is drawn by hand.
#
#  Vector only: the committed figures are .svg, written by decision_render()
#  directly (no graphics device).  Every coordinate in them came out of the
#  vendored layout engine under .tools/elk (ELK for placement, libavoid for
#  edge routing); this script only chooses which graph to draw.  Text stays
#  selectable and the image scales without banding.
#
#  Run:  Rscript dev/render-decision.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

if (!decision_solver_available()) {
  stop(
    "the vendored layout engine is not runnable: rebuild .tools/elk/node_modules ",
    "or set TIDYKASP_ELK_DIR (see .tools/elk/README.md)"
  )
}

fs::dir_create("man/figures")
fs::dir_create("results/plot-smoke")

# one schema per marker chemistry: the allele-specific core is shared, the
# marker-agnostic primer subset is the bare "primer" scheme, and a duplex
# chemistry adds its compatibility stage.  The README embeds the SNP one (the
# package's core) and the multiplex one.
for (sc in marker_schemes()$scheme) {
  svg <- file.path("results/plot-smoke", paste0("decision-", sc, ".svg"))
  plot_decision_tree(scheme = sc, out = svg)
  # the text double: same schema, greppable, next to every image
  writeLines(
    marker_decision_text(sc),
    file.path("results/plot-smoke", paste0("decision-", sc, ".txt"))
  )
  cat("wrote results/plot-smoke/decision-", sc, ".svg + .txt\n", sep = "")
}
# a raster twin of the committed figures, for a viewer that cannot open SVG.  The
# SVG stays the artefact; `rsvg` is optional and its absence never fails a render.
if (requireNamespace("rsvg", quietly = TRUE)) {
  for (f in c("decision-tree", "decision-multiplex")) {
    rsvg::rsvg_png(
      file.path("man/figures", paste0(f, ".svg")),
      file.path("results/plot-smoke", paste0(f, ".png")),
      width = 2400
    )
    cat("wrote results/plot-smoke/", f, ".png\n", sep = "")
  }
}
# the committed, README-facing images (kept in step with the code by this script)
plot_decision_tree(scheme = "snp", out = "man/figures/decision-tree.svg")
cat("wrote man/figures/decision-tree.svg\n")
plot_decision_tree(scheme = "duplex", out = "man/figures/decision-multiplex.svg")
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
      plot_decision_tree(assay, out = "results/plot-smoke/decision-fri.svg")
      writeLines(
        marker_decision_text(assay$scheme %||% "snp"),
        "results/plot-smoke/decision-fri.txt"
      )
      cat("wrote results/plot-smoke/decision-fri.svg + .txt\n")
    },
    error = function(e) cat("annotated example skipped:", conditionMessage(e), "\n")
  )
}
