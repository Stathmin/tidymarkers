# =============================================================================
#  dev/render-decision.R  --  (re)draw the README decision-tree image
# =============================================================================
#  The tree is generated from the schema in R/kasp_decision.R, which itself reads
#  the engine's constants (KASP_GATES_*, KASP_RELAX_*, KASP_TIERS_*).  Run this
#  after touching the selection policy and the committed image the README embeds
#  stays in step with the code - no hand-drawn diagram to drift.
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
p <- plot_decision_tree() # the bare schema: the policy, not one assay's luck
ggplot2::ggsave("man/figures/decision-tree.png", p, width = 9, height = 13, dpi = 110, bg = "white")
cat("wrote man/figures/decision-tree.png\n")

# an annotated example (a real assay's traversal) for eyeballing, git-ignored
fs::dir_create("results/plot-smoke")
tryCatch(
  {
    genome <- ath_genome()
    gene <- ensembl_gene("AT4G00650")
    variants <- ensembl_variants(gene$chrom, gene$start, gene$end)
    snp <- pick_target_snp(variants, genome, gene)
    assay <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = variants, flank = 150L))
    ggplot2::ggsave(
      "results/plot-smoke/decision-fri.png", plot_decision_tree(assay),
      width = 9, height = 13, dpi = 110, bg = "white"
    )
    cat("wrote results/plot-smoke/decision-fri.png\n")
  },
  error = function(e) cat("annotated example skipped:", conditionMessage(e), "\n")
)
