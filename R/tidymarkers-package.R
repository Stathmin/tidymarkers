#' tidymarkers: tidy design and quality control of genetic markers
#'
#' An opinionated, reproducible toolkit for designing and quality-controlling
#' genetic markers across many chemistries: KASP, CAPS/dCAPS, tetra-primer ARMS,
#' HRM, Sanger, SSR, InDel-size, TaqMan/dPCR, SNaPshot, LAMP/RPA and CRISPR-Dx.
#' It is genome-agnostic: give it a reference (a FASTA path, an indexed
#' [marker_genome()], a `DNAStringSet` or a Bioconductor `BSgenome`) and variants
#' (a VCF via [read_variants_vcf()], or the tidy variant table the package uses),
#' and it designs and scores markers for the chosen chemistry.  Everything is
#' done in R: thermodynamics come from `TmCalculator`, `DECIPHER` and
#' OligoArrayAux, sequence work from `Biostrings`/`Rsamtools`, annotation from
#' `VariantAnnotation`/`TxDb`, and the visuals from dedicated plotting libraries
#' (`gggenes`, `ggrepel`, `patchwork`, `ggraph`).
#'
#' An *Arabidopsis thaliana* setup (TAIR10 genome + 1001 Genomes variation via
#' Ensembl) is bundled as a reproducible demo ([ath_genome()], [ensembl_gene()],
#' [ensembl_variants()]); see `vignette("kasp-intro")` and
#' `vignette("custom-genome")` for using your own genome and variants.  Other
#' vignettes cover the remaining chemistries and the workflow:
#' `vignette("markers-probe")`, `vignette("markers-isothermal")`,
#' `vignette("markers-crispr-sequencing")`, `vignette("markers-gel")` and
#' `vignette("reading-results")`.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"
