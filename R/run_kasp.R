# =============================================================================
#  run_kasp.R  --  end-to-end KASP assay design for an Arabidopsis gene
# =============================================================================
#  Run from the project root:   Rscript R/run_kasp.R
#
#  Pipeline:
#    1. load the TAIR10 genome (Ensembl Plants FASTA) and the FRI gene model
#    2. pull real 1001-Genomes SNPs and pick a target variant
#    3. design AS + common primers, score thermodynamics (TmCalculator),
#       Taq-based discrimination (DECIPHER) and specificity (Biostrings/BLAST)
#    4. write results/ (csv, rds, png) for the report to consume
#
#  External tool required: OligoArrayAux (`hybrid-min`) for DECIPHER, e.g.
#     conda install -c bioconda oligoarrayaux
#  Optional: NCBI BLAST+ for the (slow) genome-wide BLAST step.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(Biostrings)
  library(ggmsa)
})

source(here::here("R", "kasp_functions.R"))
source(here::here("R", "kasp_plots.R"))

# Toggle the (slow) genome-wide BLAST step; the fast seed screen always runs.
RUN_BLAST <- FALSE

if (!has_oligoarrayaux()) {
  rlang::abort(c("OligoArrayAux ('hybrid-min') not found on PATH.",
                 i = "conda install -c bioconda oligoarrayaux",
                 i = "or build it from UNAFold and put it on PATH."))
}

# --- 1. genome + gene --------------------------------------------------------
genome  <- ath_genome()          # TAIR10 FASTA (Ensembl Plants), cached
gene_id <- "AT4G00650"           # FRI (FRIGIDA): the textbook flowering-time locus

gene <- ensembl_gene(gene_id)
print(gene)

# --- 2. target SNP -----------------------------------------------------------
variants <- ensembl_variants(gene$chrom, gene$start, gene$end)
message(nrow(variants), " variants in ", gene$symbol)

snp <- pick_target_snp(variants, genome, gene)
# Force a specific marker when you have one in mind, e.g. the nonsense allele:
#   snp <- pick_target_snp(variants, genome, gene, variant_id = "tmp_4_269194_A_T")
print(snp)

# --- 3. design + QC ----------------------------------------------------------
# Two designer knobs are exposed as real arguments:
#   strands  : search the allele-specific primer on the "+" and "-" strands
#   tm_range : accepted Tm window (C) for both AS and common primers
assay <- assemble_kasp_assay(genome, snp, variants = variants, flank = 150L,
                             strands = c("+", "-"), tm_range = c(58, 64),
                             clean_3p = 5L, run_blast = RUN_BLAST)

# --- 4. report the assay -----------------------------------------------------
cat("\n=== KASP assay for", gene$symbol, "(strand", assay$strand, ") ===\n")
assay$assay |>
  dplyr::select(rank, strand, as_ref, as_alt, COM = com, discrimination, as_gc,
                amplicon, n_var, as_n_var, as_n_3p) |>
  print()

cat("\nPredicted PCR efficiencies (DECIPHER, Taq-aware):\n")
print(assay$efficiency)

cat("\nPrimers (5'->3', tails in brackets for synthesis):\n")
cat("  AS-ref (FAM):", assay$tails[["AS_ref"]], "\n")
cat("  AS-alt (HEX):", assay$tails[["AS_alt"]], "\n")
cat("  COM        :", assay$tails[["COM"]], "\n\n")

cat("Thermodynamics (TmCalculator):\n");                    print(assay$thermo)
cat("\nDuplex dG37, kcal/mol (more negative = stronger dimer):\n")
print(assay$complementarity)
cat("\n3'-seed specificity (hits should be 1):\n");          print(assay$specificity)
if (!is.null(assay$blast)) {
  cat("\nTop BLAST hits:\n")
  print(utils::head(dplyr::arrange(assay$blast, .data$evalue), 10))
}

# --- 5. discrimination figure ------------------------------------------------
p_disc <- assay$as_candidates |>
  dplyr::slice_max(.data$discrimination, n = 15) |>
  dplyr::mutate(label = paste0("strand ", .data$strand, " / len ", .data$as_length,
                               " / mm ", .data$stab_at, " / ", .data$stab_base)) |>
  ggplot(aes(x = discrimination, y = stats::reorder(label, discrimination))) +
  geom_col(fill = "darkorange") +
  labs(x = "discrimination (on- minus off-target efficiency)", y = NULL,
       title = "Allele-specific discrimination of candidate AS primers") +
  theme_minimal()

# --- 6. persist --------------------------------------------------------------
fs::dir_create(here::here("results"))
readr::write_csv(assay$assay,          here::here("results", "kasp_assay.csv"))
readr::write_csv(assay$as_candidates,  here::here("results", "kasp_as_candidates.csv"))
readr::write_csv(assay$thermo,         here::here("results", "kasp_thermo.csv"))
readr::write_csv(assay$specificity,    here::here("results", "kasp_specificity.csv"))
readr::write_csv(assay$complementarity, here::here("results", "kasp_complementarity.csv"))
if (!is.null(assay$blast)) {
  readr::write_csv(assay$blast, here::here("results", "kasp_blast_hits.csv"))
}
saveRDS(assay, here::here("results", "kasp_assay.rds"))
ggsave(here::here("results", "kasp_discrimination.png"), p_disc, width = 7, height = 4)

# --- 7. library-driven views (see R/kasp_plots.R) ---------------------------
# ggmsa emits internal warnings (aes_ deprecation + a manual-scale note) that are
# not actionable here; silenciar only while rendering the ggmsa-based panels.
ggsave(here::here("results", "plot_primers.png"),    plot_primers(assay, variants, snp),           width = 8,   height = 4.5)
ggsave(here::here("results", "plot_variants.png"),    plot_variant_landscape(gene, variants, snp, assay), width = 8, height = 3.5)
suppressWarnings(ggsave(here::here("results", "plot_alignment.png"), plot_alignment(assay),  width = 8, height = 3))
suppressWarnings(ggsave(here::here("results", "plot_dimers.png"),    plot_dimer_qc(assay),  width = 9, height = 7))

# --- 8. oligo ordering sheets (IDT / LGC) -----------------------------------
primers <- export_primers(assay, gene)
cat("\nOrder sheet:\n"); print(primers)

message("\nDone. See results/ and knit report/kasp_report.Rmd for the narrative.")
