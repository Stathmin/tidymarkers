# =============================================================================
#  dev/render-targets.R  --  the assay map for every bundled example target
# =============================================================================
#  The CSV of typed example targets (inst/extdata/example_targets.csv) is our
#  test set; this renders the same `plot_primers()` view for each locus, so the
#  sites can be eyeballed side by side (SNP / MNP / insertion / deletion /
#  multi-allelic).  A locus the strict gates doom is reported, not fatal.
#
#  Run:  Rscript dev/render-targets.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

if (!nzchar(Sys.which("hybrid-min"))) stop("OligoArrayAux not on PATH - see data-raw/install-oligoarrayaux.R")

genome <- ath_genome()
gene <- ensembl_gene("AT4G00650")
variants <- ensembl_variants(gene$chrom, gene$start, gene$end)
out <- fs::path(root, "results", "plot-smoke")
fs::dir_create(out)

ex <- kasp_example()
for (i in seq_len(nrow(ex))) {
  r <- ex[i, ]
  snp <- tibble::tibble(chrom = r$chrom, pos = as.integer(r$pos), ref = r$ref, alt = r$alt)
  tag <- gsub("[^A-Za-z0-9]+", "_", r$id)
  msg <- tryCatch(
    {
      a <- suppressWarnings(assemble_kasp_assay(genome, snp,
        variants = variants, flank = 150L, strands = c("+", "-")
      ))
      p <- plot_primers(a, variants, snp, genome = genome, features = gene)
      ggplot2::ggsave(fs::path(out, paste0("target_", tag, ".png")), p,
        width = 11, height = 9, dpi = 110, bg = "white"
      )
      glue::glue("ok ({r$type}) strand {a$strand} disc {round(a$efficiency$discrimination[1], 2)}")
    },
    error = function(e) glue::glue("doomed: {conditionMessage(e)}")
  )
  cat(sprintf("%-22s %-6s %s:%s %s>%s -> %s\n", r$id, r$type, r$chrom, r$pos, r$ref, r$alt, msg))
}
cat("wrote target_*.png into", out, "\n")
