# =============================================================================
#  dev/dissect-targets.R  --  why is each CSV target doomed (or not)?
# =============================================================================
#  Dumps every example target's design record to a temp dir and prints, for the
#  doomed ones, which gate binds (footprint variant count, 3' window, GC/Tm,
#  discrimination, structure).  Run: Rscript dev/dissect-targets.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1]])))) else getwd()
suppressMessages(pkgload::load_all(root, quiet = TRUE))

out <- fs::path("/tmp/kasp-dissect")
fs::dir_create(out)

genome <- ath_genome()
ex <- marker_example()

rows <- lapply(seq_len(nrow(ex)), function(i) {
  r <- ex[i, ]
  snp <- tibble::tibble(chrom = r$chrom, pos = as.integer(r$pos), ref = r$ref, alt = r$alt)
  tag <- gsub("[^A-Za-z0-9]+", "_", r$id)
  # per-locus neighbourhood now that the targets span genes and chromosomes
  gene <- tryCatch(ensembl_gene(r$gene), error = function(e) NULL)
  win <- if (!is.null(gene)) c(gene$start - 500L, gene$end + 500L) else c(snp$pos - 1000L, snp$pos + 1000L)
  variants <- ensembl_variants(r$chrom, max(1L, win[1L]), win[2L])
  a <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = variants, flank = 150L, cores = 1L))
  saveRDS(a, fs::path(out, paste0(tag, ".rds")))

  inwin <- sum(variants$pos[variants$chrom == r$chrom] >= snp$pos - 150L &
    variants$pos[variants$chrom == r$chrom] <= snp$pos + 150L) - 1L
  ac <- a$as_candidates
  has <- !is.null(ac) && nrow(ac) > 0L
  tibble::tibble(
    id = r$id, type = r$type, pos = snp$pos, win_var = inwin,
    verdict = if (inherits(a, "marker_doomed")) "doomed" else "ok",
    stage = if (inherits(a, "marker_doomed")) a$stage else "-",
    n_as = if (has) nrow(ac) else 0L,
    n_var0 = if (has) sum(ac$n_var == 0L) else NA_integer_,
    min_n_var = if (has) min(ac$n_var) else NA_integer_,
    n_3p0 = if (has) sum(ac$n_3p == 0L) else NA_integer_,
    max_disc = if (has) round(max(ac$discrimination, na.rm = TRUE), 3) else NA_real_,
    reason = if (inherits(a, "marker_doomed")) a$reason else NA_character_
  )
})
tab <- purrr::list_rbind(rows)
csv <- fs::path(root, "vignettes", "figures", "dissect.csv")
fs::dir_create(fs::path_dir(csv))
readr::write_csv(tab, csv)
cat(sprintf(
  "%-20s %-8s %7s %5s %5s %6s %7s %6s\n",
  "id", "type", "pos", "winv", "n_as", "n_v0", "min_v", "disc"
))
for (i in seq_len(nrow(tab))) {
  cat(sprintf(
    "%-20s %-8s %7d %5d %5d %6s %7s %6s\n",
    substr(tab$id[i], 1, 20), tab$type[i], tab$pos[i], tab$win_var[i], tab$n_as[i],
    format(tab$n_var0[i]), format(tab$min_n_var[i]), format(tab$max_disc[i])
  ))
}
cat("\ndumped", nrow(tab), "records into", out, "\n")
