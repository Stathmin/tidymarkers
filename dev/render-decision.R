# =============================================================================
#  dev/render-decision.R  --  (re)draw the README decision-tree image
# =============================================================================
#  The trees are generated from the schema in R/marker_decision.R, which itself
#  reads the engine's constants (KASP_GATES_*, KASP_RELAX_*, KASP_TIERS_*) and
#  the chemistry registry (marker_schemes()).  One image per marker chemistry is
#  written to the git-ignored smoke dir, and the README-facing SNP and multiplex
#  images are committed.  Run this after touching the selection policy or the
#  scheme registry and the diagrams stay in step with the code - nothing is drawn
#  by hand.
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

# one schema per marker chemistry: the allele-specific core is shared, the
# marker-agnostic primer subset is the bare "primer" scheme, and a duplex
# chemistry adds its compatibility stage.  The README embeds the SNP one (the
# package's core) and the multiplex one.
for (sc in marker_schemes()$scheme) {
  p <- plot_decision_tree(scheme = sc)
  ggplot2::ggsave(
    file.path("results/plot-smoke", paste0("decision-", sc, ".png")), p,
    width = 9, height = 13, dpi = 110, bg = "white"
  )
  cat("wrote results/plot-smoke/decision-", sc, ".png\n", sep = "")
}
# the committed, README-facing images (kept in step with the code by this script)
ggplot2::ggsave(
  "man/figures/decision-tree.png", plot_decision_tree(scheme = "snp"),
  width = 9, height = 13, dpi = 110, bg = "white"
)
cat("wrote man/figures/decision-tree.png\n")
ggplot2::ggsave(
  "man/figures/decision-multiplex.png", plot_decision_tree(scheme = "duplex"),
  width = 9, height = 13, dpi = 110, bg = "white"
)
cat("wrote man/figures/decision-multiplex.png\n")

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
      ggplot2::ggsave(
        "results/plot-smoke/decision-fri.png", plot_decision_tree(assay),
        width = 9, height = 13, dpi = 110, bg = "white"
      )
      cat("wrote results/plot-smoke/decision-fri.png\n")
    },
    error = function(e) cat("annotated example skipped:", conditionMessage(e), "\n")
  )
}
