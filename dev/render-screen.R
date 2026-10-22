# =============================================================================
#  dev/render-screen.R  --  render the cheap target screen on the FRI demo
# =============================================================================
#  Screens every SNP in FRI (1001 Genomes) with screen_snp_targets() and writes
#  the per-variant "promising / doomed" plot for eyeballing.  Needs no
#  thermodynamics, so it runs without OligoArrayAux.
#
#  Run:  Rscript dev/render-screen.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

genome <- ath_genome()
gene <- ensembl_gene("AT4G00650")
variants <- ensembl_variants(gene$chrom, gene$start, gene$end)

t0 <- Sys.time()
screen <- screen_snp_targets(variants, genome)
dt <- as.numeric(Sys.time() - t0, units = "secs")
message(sprintf("%d target SNPs screened in %.2f s (%.1f ms/SNP)", nrow(screen), dt, 1000 * dt / nrow(screen)))
print(table(screen$verdict))

out <- fs::path(root, "vignettes", "figures")
fs::dir_create(out)
f <- fs::path(out, "screen_fri.png")
ggplot2::ggsave(f, plot_snp_screen(screen, gene), width = 10, height = 4, dpi = 110, bg = "white")
cat("wrote", f, "\n")

# the divergent ALT/REF test target: G>T at 269907 with a 1-bp deletion at 269904
print(as.data.frame(screen[screen$pos %in% c(269904L, 269907L), c("pos", "ref", "alt", "n_ok", "score", "verdict", "reasons")]))
