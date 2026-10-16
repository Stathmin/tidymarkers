# =============================================================================
#  dev/render-real.R  --  render the plots on the real Arabidopsis demo
# =============================================================================
#  Runs the bundled TAIR10 / FRI (AT4G00650) / 1001-Genomes pipeline (data are
#  cached under marker_cache_dir(), so no network is needed after the first run)
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
# NB: with HARD limits (no variant in the AS 3' window, real discrimination) the
# default pick_target_snp() choice (269745) is *doomed* - it has no designable AS
# placement.  pick_target_snp() does not yet know about designability, so the demo
# pins a locus that is designable (269666 G>A).  See ROADMAP (doomed target +
# next-candidate iteration).
snp <- tibble::tibble(chrom = "Chr4", pos = 269666L, ref = "G", alt = "A")
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

# --- a real co-segregating substitution --------------------------------------
# FRI 269043 G>A has a substitution in its '+' 3' window (269042) that, being
# *unlinked*, would trip the HARD "no variant in the AS 3' window" gate and doom
# the strand.  Declaring 269042 linked rides it with the ALT haplotype (baked into
# the ALT primer), so that gate no longer fires and the strand becomes designable.
snp_l <- tibble::tibble(chrom = "Chr4", pos = 269043L, ref = "G", alt = "A")
linked <- tibble::tibble(target = 269043L, linked = 269042L)
assay_l <- assemble_kasp_assay(genome, snp_l, variants = variants, flank = 150L,
                               strands = "+", linked = linked)
save("primers_fri_linked", plot_primers(assay_l, variants, snp_l, genome = genome, features = gene))

cat("done\n")
