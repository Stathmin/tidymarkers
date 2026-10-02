# =============================================================================
#  dev/render-targets.R  --  the assay map for every bundled example target
# =============================================================================
#  The CSV of typed example targets (inst/extdata/example_targets.csv) is our
#  test set; this renders the same `plot_primers()` view for each locus, so the
#  sites can be eyeballed side by side (SNP / MNP / insertion / deletion /
#  multi-allelic).  A locus the strict gates doom is reported, not fatal.
#
#  Loci are independent, so they are rendered across cores (fork on Unix); each
#  assay is assembled serially (`cores = 1`) so cores are not oversubscribed.
#  Set TIDYKASP_CORES to cap the worker count.
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

render_one <- function(i) {
  r <- ex[i, ]
  snp <- tibble::tibble(chrom = r$chrom, pos = as.integer(r$pos), ref = r$ref, alt = r$alt)
  tag <- gsub("[^A-Za-z0-9]+", "_", r$id)
  msg <- tryCatch(
    {
      a <- suppressWarnings(assemble_kasp_assay(genome, snp,
        variants = variants, flank = 150L, strands = c("+", "-"), cores = 1L
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
  invisible(NULL)
}

ncores <- min(
  as.integer(Sys.getenv("TIDYKASP_CORES", unset = "1")),
  nrow(ex)
)
if (ncores > 1L && .Platform$OS.type == "unix") {
  bp <- BiocParallel::MulticoreParam(ncores)
  invisible(BiocParallel::bplapply(seq_len(nrow(ex)), render_one, BPPARAM = bp))
} else {
  invisible(lapply(seq_len(nrow(ex)), render_one))
}
cat("wrote target_*.png into", out, "\n")
