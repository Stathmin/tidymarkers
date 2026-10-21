suppressPackageStartupMessages({
  library(tidymarkers)
  library(dplyr)
  library(ggplot2)
  library(Biostrings)
})

out <- "/root/development/Merovingen/data/scratch/reports/text/figures"
fs::dir_create(out)

## 1. decision graph: the KASP gate schema (SVG)
g <- plot_decision_tree(scheme = "snp")
ggsave(file.path(out, "kasp_decision.svg"), g, width = 13, height = 8,
       device = svglite::svglite)
cat("decision svg done\n")

## 2. Arabidopsis FRI demo: design, then primer map + variant landscape + checks
genome  <- ath_genome()
gene_id <- "AT4G00650"
gene <- ensembl_gene(gene_id)
variants <- ensembl_variants(gene$chrom, gene$start, gene$end)
snp <- pick_target_snp(variants, genome, gene)
assay <- assemble_kasp_assay(genome, snp, variants = variants, flank = 150L,
                             strands = c("+", "-"), tm_range = c(58, 64),
                             clean_3p = 5L, run_blast = FALSE)
cat("assay done, rank1 amplicon", assay$assay$amplicon[1], "\n")

ggsave(file.path(out, "kasp_primers.svg"),
       plot_primers(assay, variants, snp, genome = genome, features = gene),
       width = 11, height = 7.5, device = svglite::svglite)
ggsave(file.path(out, "kasp_variants.svg"),
       plot_variant_landscape(gene, variants, snp, assay),
       width = 9, height = 3.6, device = svglite::svglite)
ggsave(file.path(out, "kasp_checks.svg"), plot_checks(assay),
       width = 9, height = 4.2, device = svglite::svglite)

print(assay$assay |> dplyr::select(rank, as_ref, as_alt, com, discrimination))
print(marker_design_checks(assay))
