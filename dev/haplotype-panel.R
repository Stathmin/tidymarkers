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
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

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

# three haplotypes over the ~100 bp around the target, at rising frequencies
hap <- function(id, freq, pos, kind = c("snp", "indel")) {
  kind <- match.arg(kind)
  b <- substring(W, pos, pos)
  alt <- if (kind == "snp") setdiff(c("A", "C", "G", "T"), b)[1L] else paste0(b, "TG")
  tibble::tibble(
    haplotype = id, freq = freq, kind = kind, pos = pos,
    ref = b, alt = alt, pop_freq = freq, chrom = "Chr1", consequence = "missense_variant"
  )
}
haps <- dplyr::bind_rows(
  hap("H1_rare_snp_in_AS", 0.01, aa$as_start + 8L, "snp"),
  hap("H2_rare_indel_in_AS", 0.01, aa$as_start + 12L, "indel"),
  hap("H3_common_indel_in_AS", 0.50, aa$as_start + 12L, "indel")
)

rows <- lapply(seq_len(nrow(haps)), function(i) {
  h <- haps[i, ]
  variants <- dplyr::bind_rows(
    tibble::tibble(
      chrom = "Chr1", pos = 250L, consequence = "missense_variant"
    ),
    dplyr::select(h, "chrom", "pos", "ref", "alt", "pop_freq", "consequence")
  )
  a <- suppressWarnings(assemble_kasp_assay(genome, snp, variants = variants, flank = 150L, cores = 1L))
  saveRDS(list(assay = a, haplotype = h), fs::path(out, paste0(h$haplotype, ".rds")))
  if (!inherits(a, "kasp_doomed")) {
    p <- plot_primers(a, variants, snp, genome = genome)
    ggplot2::ggsave(fs::path(out, paste0(h$haplotype, ".png")), p,
      width = 11, height = 8, dpi = 100, bg = "white"
    )
  }
  tibble::tibble(
    haplotype = h$haplotype, freq = h$freq, kind = h$kind,
    verdict = if (inherits(a, "kasp_doomed")) "doomed" else "ok",
    reported = if (inherits(a, "kasp_doomed")) NA_integer_ else a$assay$as_n_var[1L] + a$assay$n_var[1L],
    foot_dropout = if (inherits(a, "kasp_doomed")) NA_real_ else round(a$assay$as_foot_dropout[1L], 4L),
    reason = if (inherits(a, "kasp_doomed")) a$reason else NA_character_
  )
})
tab <- purrr::list_rbind(rows)
print(as.data.frame(tab), row.names = FALSE)
cat("\ndumped", nrow(tab), "records into", out, "\n")
