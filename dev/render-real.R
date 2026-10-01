# =============================================================================
#  dev/render-real.R  --  render the plots on the real Arabidopsis demo
# =============================================================================
#  Runs the bundled TAIR10 / FRI (AT4G00650) / 1001-Genomes pipeline (data are
#  cached under kasp_cache_dir(), so no network is needed after the first run)
#  and writes the plots to results/plot-smoke/ for eyeballing.
#
#  Run:  Rscript dev/render-real.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

if (!nzchar(Sys.which("hybrid-min"))) stop("OligoArrayAux not on PATH - see data-raw/install-oligoarrayaux.R")

# --- the real demo: TAIR10 genome + FRI gene model + 1001-Genomes variants ----
genome <- ath_genome()
gene <- ensembl_gene("AT4G00650") # FRI (FRIGIDA)
variants <- ensembl_variants(gene$chrom, gene$start, gene$end)
snp <- pick_target_snp(variants, genome, gene)
message("target: ", snp$chrom, ":", snp$pos, " ", snp$ref, ">", snp$alt, " (", gene$symbol, ")")

assay <- assemble_kasp_assay(genome, snp, variants = variants, flank = 150L,
                             strands = c("+", "-"), tm_range = c(58, 64), clean_3p = 5L)

out <- fs::path(root, "results", "plot-smoke")
fs::dir_create(out)
save <- function(name, plot, width = 10, height = 7.5) {
  f <- fs::path(out, paste0(name, ".png"))
  tryCatch(
    {
      ggplot2::ggsave(f, plot, width = width, height = height, dpi = 110, bg = "white")
      cat("wrote", f, "\n")
    },
    error = function(e) cat("FAIL", name, ":", conditionMessage(e), "\n")
  )
}

save("primers_fri", plot_primers(assay, variants, snp, genome = genome, features = gene))
save("variants_fri", plot_variant_landscape(gene, variants, snp, assay), width = 10, height = 4.5)
save("checks_fri", plot_checks(assay), width = 3, height = 6)

cat("done\n")
