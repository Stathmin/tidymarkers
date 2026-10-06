# =============================================================================
#  dev/haplotype-panel.R  --  2-3 haplotypes / 100 bp around each target
# =============================================================================
#  Builds a small haplotype panel per target (a SNP or an indel at a chosen
#  position, each with a population frequency) and checks that the designer
#  reads frequencies: a rare variant under a primer is reported but not fatal,
#  a common one costs more.  Dumps every record to /tmp/kasp-hap for dissection.
#
#  Run: Rscript dev/haplotype-panel.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1]])))) else getwd()
suppressMessages(pkgload::load_all(root, quiet = TRUE))
haplotype_variants <- getFromNamespace("haplotype_variants", "tidykasp")

out <- fs::path("/tmp/kasp-hap")
fs::dir_create(out)

# A clean target on a fixed synthetic genome (a real one is one long design).
W <- withr::with_seed(4L, paste(sample(c("A", "C", "G", "T"), 500L, TRUE), collapse = ""))
genome <- Biostrings::DNAStringSet(c(Chr1 = W))
ref <- as.character(Biostrings::subseq(genome[["Chr1"]], 250L, 250L))
snp <- tibble::tibble(chrom = "Chr1", pos = 250L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
base <- suppressWarnings(assemble_kasp_assay(genome, snp, flank = 150L, cores = 1L))
aa <- base$assay[1L, ]
cat(sprintf(
  "baseline: AS %d-%d, COM %d-%d (strand %s)\n",
  aa$as_start, aa$as_end, aa$com_start, aa$com_end, base$strand
))

# three haplotypes over the ~100 bp around the target: a REF haplotype and two
# ALT haplotypes, one rare variant sitting inside the AS footprint and one
# variant co-segregating with the ALT allele (rides both ALT haplotypes).
tgt <- snp$alt
hap_row <- function(haplotype, freq, pos, alt = NULL) {
  ref <- substring(W, pos, pos)
  if (is.null(alt) || alt == ref) alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  tibble::tibble(
    haplotype = haplotype, freq = freq, chrom = "Chr1", pos = pos,
    ref = ref, alt = alt
  )
}
haps <- dplyr::bind_rows(
  hap_row("H_ref", 0.70, 215L),
  hap_row("H_alt1", 0.20, 250L, tgt),
  hap_row("H_alt1", 0.20, 262L),
  hap_row("H_alt2", 0.10, 250L, tgt),
  hap_row("H_alt2", 0.10, 262L), # co-segregates with the ALT allele
  hap_row("H_alt2", 0.10, aa$as_start + 9L) # rare variant inside the AS primer
)

a <- suppressWarnings(assemble_kasp_assay(genome, snp, haplotypes = haps, flank = 150L, cores = 1L))
saveRDS(list(assay = a, haplotypes = haps), fs::path(out, "panel.rds"))
if (!inherits(a, "kasp_doomed")) {
  v <- haplotype_variants(haps, snp)$variants
  p <- plot_primers(a, v, snp, genome = genome)
  ggplot2::ggsave(fs::path(out, "panel.png"), p, width = 11, height = 9, dpi = 100, bg = "white")
}
cat("\nverdict:", if (inherits(a, "kasp_doomed")) paste("DOOMED:", a$reason) else "ok", "\n")
if (!inherits(a, "kasp_doomed")) {
  cat("AS footprint variants:", a$assay$as_n_var[1L], " dropout:", round(a$assay$as_foot_dropout[1L], 4L), "\n")
  gpos <- a$selection$linked$offset + snp$pos - 151L # template offset -> genomic
  cat("linked (co-segregating) positions:", paste(gpos, collapse = ", "), "\n")
  print(as.data.frame(a$haplotype_hits))
}
cat("\ndumped into", out, "\n")
