#' tidykasp: tidy design and quality control of KASP genotyping assays
#'
#' An opinionated, reproducible toolkit for designing and quality-controlling
#' Kompetitive Allele-Specific PCR (KASP) assays.  It is genome-agnostic: give it
#' a reference (a FASTA path, an indexed [kasp_genome()], a `DNAStringSet` or a
#' Bioconductor `BSgenome`) and variants (a VCF via [read_variants_vcf()], or the
#' tidy variant table the package uses), and it designs allele-specific and common
#' primers on either strand and scores them.  Everything is done in R:
#' thermodynamics come from `TmCalculator`, `DECIPHER` and OligoArrayAux, sequence
#' work from `Biostrings`/`Rsamtools`, annotation from `VariantAnnotation`/`TxDb`,
#' and the visuals from dedicated plotting libraries (`gggenes`,
#' `ggrepel`, `patchwork`, `ggraph`).
#'
#' An *Arabidopsis thaliana* setup (TAIR10 genome + 1001 Genomes variation via
#' Ensembl) is bundled as a reproducible demo ([ath_genome()], [ensembl_gene()],
#' [ensembl_variants()]); see `vignette("kasp-intro")` and
#' `vignette("custom-genome")` for using your own genome and variants.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"
