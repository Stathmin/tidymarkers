# =============================================================================
#  dev/render-plots.R  --  render every plot_* on one real assay (smoke test)
# =============================================================================
#  Builds a small random assay and writes each plot to a temp dir, so a broken
#  geom / aes / library call shows up immediately (the tests only assert classes).
#  Needs OligoArrayAux on PATH (see data-raw/install-oligoarrayaux.R).
#
#  Run:  Rscript dev/render-plots.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

if (!nzchar(Sys.which("hybrid-min"))) stop("OligoArrayAux not on PATH - see data-raw/install-oligoarrayaux.R")

make_assay <- function(seed = 4L, n = 20000L, pos = 10000L) {
  rand <- function(len) paste(sample(c("A", "C", "G", "T"), len, TRUE), collapse = "")
  for (s in seed:(seed + 19L)) {
    withr::with_seed(s, genome <- Biostrings::DNAStringSet(c(Chr1 = rand(n), Chr2 = rand(n), Chr3 = rand(n))))
    ref <- as.character(Biostrings::subseq(genome[["Chr1"]], pos, pos))
    snp <- tibble::tibble(
      chrom = "Chr1", pos = pos, ref = ref,
      alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
    )
    vs <- tibble::tibble(chrom = "Chr1", pos = pos, consequence = "missense_variant")
    a <- try(suppressWarnings(assemble_kasp_assay(genome, snp, variants = vs, flank = 150L)), silent = TRUE)
    if (!inherits(a, "try-error") && !inherits(a, "kasp_doomed")) {
      return(list(assay = a, snp = snp, variants = vs, genome = genome))
    }
  }
  stop("no random genome assembled into an assay")
}

m <- make_assay()
out <- fs::path(root, "results", "plot-smoke") # git-ignored, easy to eyeball
fs::dir_create(out)

save <- function(name, plot) {
  f <- fs::path(out, paste0(name, ".png"))
  tryCatch(
    {
      ggplot2::ggsave(f, plot, width = 9.5, height = 7.5, dpi = 110, bg = "white")
      cat("wrote", f, "\n")
    },
    error = function(e) cat("FAIL", name, ":", conditionMessage(e), "\n")
  )
}

save("primers", plot_primers(m$assay, m$variants, m$snp, genome = m$genome))
save("checks", plot_checks(m$assay))
save("oligo_structure", plot_oligo_structure(m$assay))
save("self_complementarity", plot_self_complementarity(m$assay$assay$as_ref[1L]))
cat("done\n")
