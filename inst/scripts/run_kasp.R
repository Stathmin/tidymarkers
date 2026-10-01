# =============================================================================
#  run_kasp.R  --  end-to-end KASP assay design for an Arabidopsis gene
# =============================================================================
#  Demo script (kept out of the package namespace, hence under inst/scripts).
#  Run from the project root:   Rscript inst/scripts/run_kasp.R
#
#  Pipeline:
#    1. load the TAIR10 genome (Ensembl Plants FASTA) and the FRI gene model
#    2. pull real 1001-Genomes SNPs and pick a target variant
#    3. design AS + common primers, score thermodynamics (TmCalculator),
#       Taq-based discrimination (DECIPHER) and specificity (Biostrings/BLAST)
#    4. write results/ (consolidated csv, rds, png) for the report to consume
#
#  External tool required: OligoArrayAux (`hybrid-min`) for DECIPHER, e.g.
#     conda install -c bioconda oligoarrayaux
#  Optional: NCBI BLAST+ for the (slow) genome-wide BLAST step; a cached
#  Rsubread FM-index is the index-based alternative and is built once.
# =============================================================================

suppressPackageStartupMessages({
  library(tidykasp)
  library(here)
  library(dplyr)
  library(ggplot2)
  library(Biostrings)
  library(ggmsa)               # geom_seqlogo needs it attached, not just namespaced
})

# Toggle the (slow) genome-wide BLAST step; the fast seed screen always runs.
RUN_BLAST <- FALSE
# Build/attach an Rsubread FM-index for the whole-primer genome scan (cached).
RUN_FM <- TRUE

if (!has_oligoarrayaux()) {
  rlang::abort(c("OligoArrayAux ('hybrid-min') not found on PATH.",
                 i = "conda install -c bioconda oligoarrayaux",
                 i = "or build it from UNAFold and put it on PATH."))
}

# --- 1. genome + gene --------------------------------------------------------
genome  <- ath_genome()          # TAIR10 FASTA (Ensembl Plants), cached
gene_id <- "AT4G00650"           # FRI (FRIGIDA): the textbook flowering-time locus

# An FM-index of the genome, built once and reused on later runs (Rsubread).
# This is the index-based alternative to BLAST for the whole-primer scan below;
# `kasp_genome(..., fm_index = )` would let `assemble_kasp_assay(run_blast = TRUE)`
# use it too.
fm_index <- if (RUN_FM) kasp_fm_index(genome, basename = fs::path(kasp_cache_dir(), "tair10_fm"))

gene <- ensembl_gene(gene_id)
print(gene)

# --- 2. target SNP -----------------------------------------------------------
variants <- ensembl_variants(gene$chrom, gene$start, gene$end)
message(nrow(variants), " variants in ", gene$symbol)

snp <- pick_target_snp(variants, genome, gene)
# Force a specific marker when you have one in mind, e.g. the nonsense allele:
#   snp <- pick_target_snp(variants, genome, gene, variant_id = "tmp_4_269194_A_T")
print(snp)

cat("\nEasyKASP direction hint (GC of each flank):\n")
print(strand_gc_hint(genome, snp$chrom, snp$pos))

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
if (!is.null(fm_index)) {
  cat("\nWhole-primer genome scan (FM-index via Rsubread; 1 = unique):\n")
  fm_primers <- c(AS_ref = assay$assay$as_ref[1L], AS_alt = assay$assay$as_alt[1L],
                  COM = assay$assay$com[1L])
  print(specificity_fm(fm_primers, fm_index, nthreads = 2L))
}
if (!is.null(assay$blast)) {
  cat("\nTop BLAST hits:\n")
  print(utils::head(dplyr::arrange(assay$blast, .data$evalue), 10))
}

cat("\nAcceptance checks (kasp_design_checks):\n")
checks <- kasp_design_checks(assay)
print(checks)
if (!all(checks$pass)) {
  rlang::warn(glue::glue("Assay fails {sum(!checks$pass)} acceptance criterion/criteria."))
}

# --- 5-6. persist (consolidated: design / primer QC / dimers + order subsets) -
fs::dir_create(here::here("results"))
saveRDS(assay, here::here("results", "kasp_assay.rds"))
write_kasp_results(assay, gene)   # kasp_design.csv, kasp_primer_qc.csv, kasp_dimer_dg.csv, order subsets

# --- 7. library-driven views (see R/kasp_plots.R) ---------------------------
# ggmsa emits internal warnings (aes_ deprecation + a manual-scale note) that are
# not actionable here; silence them only while rendering the ggmsa-based panels.
ggsave(here::here("results", "plot_primers.png"),    plot_primers(assay, variants, snp),           width = 8,   height = 4.5)
ggsave(here::here("results", "plot_variants.png"),    plot_variant_landscape(gene, variants, snp, assay), width = 8, height = 3.5)
suppressWarnings(ggsave(here::here("results", "plot_alignment.png"), plot_alignment(assay),  width = 8, height = 3))
suppressWarnings(ggsave(here::here("results", "plot_dimers.png"),    plot_dimer_qc(assay),  width = 9, height = 7))

# how the search converges: every AS placement (and then every COM placement) is
# tried, tagged with why it loses, until the accepted assay (a GIF)
render_design_search(assay, variants, snp, path = here::here("results", "design_search.gif"))

# --- 8. oligo ordering sheet (the IDT/LGC files are column subsets of it) ----
primers <- kasp_primer_table(assay, gene)
cat("\nOrder sheet:\n"); print(primers)

message("\nDone. See results/ and knit report/kasp_report.Rmd for the narrative.")
