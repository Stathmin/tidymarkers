# =============================================================================
#  quickstart.R  --  one SNP -> one assay, as fast as the cache allows
# =============================================================================
#  Run from the project root:   Rscript inst/scripts/quickstart.R
#  Needs: kaspexample installed, OligoArrayAux on PATH, TAIR10 already cached
#  (ath_genome() downloads it once on first use).
# =============================================================================

suppressPackageStartupMessages({
  library(kaspexample)
  library(here)
})

cli::cli_h1("kaspexample quickstart")

genome <- ath_genome()                                   # cached after first run
gene   <- ensembl_gene("AT4G00650")                      # FRI

snp    <- pick_target_snp(ensembl_variants(gene$chrom, gene$start, gene$end),
                          genome, gene)
cli::cli_alert_info("Target: {snp$chrom}:{snp$pos} {snp$ref}>{snp$alt}")

assay <- assemble_kasp_assay(genome, snp)
print(assay)                                             # kasp_assay print method

chk <- kasp_design_checks(assay)
if (passes_all(assay)) {
  cli::cli_alert_success("{sum(chk$pass)}/{nrow(chk)} acceptance checks passed")
} else {
  cli::cli_alert_warning("{sum(!chk$pass)} of {nrow(chk)} checks failed")
}

invisible(assay)
