# =============================================================================
#  dev/render-targets.R  --  the assay map for every bundled example target
# =============================================================================
#  The CSV of typed example targets (inst/extdata/example_targets.csv) is our
#  test set; this renders the same `plot_primers()` view for each locus, so the
#  sites can be eyeballed side by side (SNP / MNP / insertion / deletion /
#  multi-allelic).  The targets now span genes and chromosomes, so each locus
#  pulls its own neighbourhood variants and gene model (from the `gene` column);
#  a locus the strict gates doom is reported, not fatal.
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
out <- fs::path(root, "vignettes", "figures")
fs::dir_create(out)

ex <- marker_example()

render_one <- function(i) {
  r <- ex[i, ]
  # carry the second ALT of a multi-allelic site through to the assay/plot, so the
  # site class and the track read as multi-allelic (not just the first ALT)
  snp <- tibble::tibble(
    chrom = r$chrom, pos = as.integer(r$pos), ref = r$ref, alt = r$alt, alt2 = r$alt2
  )
  tag <- gsub("[^A-Za-z0-9]+", "_", r$id)
  msg <- tryCatch(
    {
      # the gene model and the neighbourhood variants are per-locus now
      gene <- tryCatch(ensembl_gene(r$gene), error = function(e) NULL)
      win <- if (!is.null(gene)) c(gene$start - 500L, gene$end + 500L) else c(snp$pos - 1000L, snp$pos + 1000L)
      variants <- ensembl_variants(r$chrom, max(1L, win[1L]), win[2L])
      # KASP first; when it dooms, the ladder falls through to CAPS/ARMS/HRM, so a
      # locus no KASP takes still gets a marker of some chemistry
      d <- design_marker(genome, snp, variants = variants, flank = 150L)
      # the whole ladder (every scheme's attempt + verdict) next to the image, so
      # the fallback is dissectable without re-running the design
      saveRDS(d, fs::path(out, paste0("target_", tag, ".rds")))
      p <- if (is.na(d$scheme)) {
        # nothing works: show the KASP attempt's context, so the reason is visible
        plot_primers(d$tried$marker[[1L]], variants, snp, genome = genome, features = gene)
      } else {
        switch(d$scheme,
          snp = ,
          indel = ,
          multi = plot_primers(d$marker, variants, snp, genome = genome, features = gene),
          caps = plot_caps_gel(d$marker),
          arms = plot_arms_gel(d$marker),
          hrm = plot_hrm_amplicon(d$marker)
        )
      }
      ggplot2::ggsave(fs::path(out, paste0("target_", tag, ".png")), p,
        width = 11, height = 9, dpi = 110, bg = "white"
      )
      if (is.na(d$scheme)) {
        glue::glue("NO marker ({r$type}); kasp: {d$tried$marker[[1L]]$reason}")
      } else {
        glue::glue("best = {d$scheme} ({r$type}); {paste(d$tried$scheme, d$tried$ok, sep = '=', collapse = ', ')}")
      }
    },
    error = function(e) glue::glue("error: {conditionMessage(e)}")
  )
  cat(sprintf("%-22s %-9s %s:%s %s>%s -> %s\n", r$id, r$type, r$chrom, r$pos, r$ref, r$alt, msg))
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
