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
# around nodes instead of long straight lines to DOOMED.
render_png <- function(file, dot_file, scheme = NULL, assay = NULL) {
  if (use_dot) {
    sc <- if (is.null(assay)) scheme else (assay$scheme %||% "snp")
    writeLines(decision_dot(scheme = sc, assay = assay), dot_file)
    status <- system2("dot", c("-Tpng", "-Gdpi=160", "-o", shQuote(file), shQuote(dot_file)))
    if (!identical(as.integer(status), 0L)) stop("dot failed for ", dot_file)
  } else {
    p <- if (is.null(assay)) plot_decision_tree(scheme = scheme) else plot_decision_tree(assay)
    ggplot2::ggsave(file, p, width = 18, height = 14, dpi = 200, bg = "white")
  }
}

# one schema per marker chemistry: the allele-specific core is shared, the
# marker-agnostic primer subset is the bare "primer" scheme, and a duplex
# chemistry adds its compatibility stage.  The README embeds the SNP one (the
# package's core) and the multiplex one.
for (sc in marker_schemes()$scheme) {
  png <- file.path("results/plot-smoke", paste0("decision-", sc, ".png"))
  render_png(png, sub("[.]png$", ".dot", png), scheme = sc)
  # the text double: same schema, greppable, next to every image
  writeLines(
    marker_decision_text(sc),
    file.path("results/plot-smoke", paste0("decision-", sc, ".txt"))
  )
  cat("wrote results/plot-smoke/decision-", sc, ".png + .txt + .dot\n", sep = "")
}
# the committed, README-facing images (kept in step with the code by this script)
render_png(
  "man/figures/decision-tree.png",
  "results/plot-smoke/decision-snp.dot", scheme = "snp"
)
cat("wrote man/figures/decision-tree.png\n")
render_png(
  "man/figures/decision-multiplex.png",
  "results/plot-smoke/decision-duplex.dot", scheme = "duplex"
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
      render_png(
        "results/plot-smoke/decision-fri.png",
        "results/plot-smoke/decision-fri.dot", assay = assay
      )
      writeLines(
        marker_decision_text(assay$scheme %||% "snp"),
        "results/plot-smoke/decision-fri.txt"
      )
      cat("wrote results/plot-smoke/decision-fri.png + .txt + .dot\n")
    },
    error = function(e) cat("annotated example skipped:", conditionMessage(e), "\n")
  )
}
