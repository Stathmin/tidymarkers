# =============================================================================
#  haplotypes.R  --  a transcript FASTA -> alignment -> haplotypes -> a KASP
# =============================================================================
#  Run from the project root:   Rscript inst/scripts/haplotypes.R
#
#  Takes a gene family's transcripts (here the bundled SAT2 seqdump: 11 Triticum
#  accessions), aligns them, treats one as the reference ("genome") and the rest
#  as haplotypes, then designs a KASP at a haplotype variant - all offline, no
#  genome download.  Writes an alignment view and the assay's primer plot.
#
#    read_transcripts()            a FASTA -> one sequence per accession
#    haplotypes_from_alignment()   align -> reference + haplotype table (+ MAF)
#    assemble_kasp_assay(haplotypes=)  design with the haplotype context
#    plot_alignment() / plot_primers()  the visuals
# =============================================================================

suppressPackageStartupMessages({
  library(tidymarkers)
  library(ggplot2)
})

out_dir <- "results/sat2"
fs::dir_create(out_dir)

# --- 1. FASTA -> alignment -> haplotype table --------------------------------
f <- system.file("extdata", "sat2_transcripts.fa", package = "tidymarkers")

# illustrative minor-allele frequencies for a few variants (in practice: your own
# VCF / panel table).  These override the haplotype-count pop_freq for those sites.
maf <- tibble::tibble(
  pos = c(356L, 422L, 431L),
  ref = c("T", "C", "A"),
  alt = c("C", "T", "G"),
  maf = c(0.31, 0.02, 0.18)
)

hp <- haplotypes_from_alignment(read_transcripts(f), chrom = "SAT2", maf = maf)
cat("reference:", hp$reference_name, "|", length(hp$reference), "bp\n")
cat("haplotypes:", length(unique(hp$haplotypes$haplotype)), "| variants:", nrow(hp$variants), "\n")
print(utils::head(dplyr::arrange(hp$variants, .data$pop_freq), 5L)) # the rare ones (MAF-weighted)

# the reference's GenBank annotation (JX896657.1): gene 1..1209, CDS 54..1067
features <- tibble::tibble(
  start = c(1L, 54L), end = c(1209L, 1067L),
  strand = "+", name = c("SAT2", "CDS")
)

# --- 2. the alignment view ----------------------------------------------------
ggsave(fs::path(out_dir, "alignment.png"), plot_alignment(hp, features = features),
  width = 10, height = 4.4, dpi = 110, bg = "white"
)
ggsave(fs::path(out_dir, "alignment_all.png"), plot_alignment(hp, type = "all", features = features),
  width = 12, height = 4.4, dpi = 110, bg = "white"
)

# --- 3. design a KASP at a haplotype variant ----------------------------------
# pick an interior substitution the haplotype set carries (pos 431 A>G designs)
snp <- tibble::tibble(chrom = "SAT2", pos = 431L, ref = "A", alt = "G")
assay <- assemble_kasp_assay(
  hp$genome, snp,
  variants = hp$variants, haplotypes = hp$haplotypes,
  flank = 150L
)
print(assay)

# --- 4. the assay's primer plot ----------------------------------------------
ggsave(fs::path(out_dir, "primers.png"),
  plot_primers(assay, hp$variants, snp, genome = hp$genome, features = features),
  width = 11, height = 8, dpi = 110, bg = "white"
)

message("Done. See ", out_dir, "/ (alignment.png, alignment_all.png, primers.png)")
