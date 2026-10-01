# =============================================================================
#  kasp_functions.R  --  design + QC primitives for KASP assays
#  100% R; heavy thermodynamics is delegated to packages.
# =============================================================================
#
#  KASP = two allele-specific forward primers (FAM / HEX tails, identical except
#  for the 3'-terminal SNP base) plus one common reverse primer.  Taq cannot
#  extend a 3'-mismatched primer, so discrimination is competitive annealing.
#
#  Defaults follow LGC/KBioscience guidance and EasyKASP (Zhang, Yang & Wen,
#  BMC Bioinformatics 2025, 26:292) -- defaults, not dogma; validate in your lab.
#  The machine-checkable form is kasp_design_checks()/KASP_LIMITS; the full
#  design manual is README section 3 and vignette("kasp-intro").
# =============================================================================


# --- constants ---------------------------------------------------------------

# Standard KASP 5' universal tails (added at synthesis; not used in Tm).
KASP_TAIL_FAM <- "GAAGGTGACCAAGTTCATGCT"  # "FAM" channel
KASP_TAIL_HEX <- "GAAGGTCGGAGTCAACGGATT"  # "HEX/VIC" channel

# Tm conditions (mM) for TmCalculator
PCR_COND <- list(Na = 50, Mg = 1.5, dNTPs = 0.6)
NN_TABLE <- "DNA_NN_SantaLucia_2004"
# PCR conditions (molar) for DECIPHER::CalculateEfficiencyPCR
PCR_EFF  <- list(temp = 60, P = 4e-7, ions = 0.2)

# Ensembl Plants FASTA defline -> TAIR chromosome name; the only sequences kept.
ATH_CHROM_MAP <- c("1" = "Chr1", "2" = "Chr2", "3" = "Chr3", "4" = "Chr4",
                   "5" = "Chr5", "Mt" = "ChrM", "Pt" = "ChrC")

# Design targets / thresholds (named, so they are not magic numbers)
AS_TM_TARGET    <- 58    # preferred AS-primer Tm (C); kept below the common primer
COM_TM_TARGET   <- 60    # preferred common-primer Tm (C); KASP wants COM Tm >= AS Tm
HAIRPIN_DG_WARN <- -3    # flag hairpins with dG37 below this (kcal/mol)
SEED_K          <- 15L   # 3'-seed length for the fast specificity screen
SEED_MISMATCH   <- 1L

# Ensembl consequence wording -> design priority (lower = preferred)
# Sequence-Ontology term -> design priority (lower = preferred).  These are the
# canonical SO terms exactly as Ensembl/VEP spell them, so matching is by full
# token (not substring).  `SO_RANK_OTHER` sits one below the worst known term.
SO_TERM_RANK <- c(
  stop_gained              = 1L,
  stop_lost                = 1L,
  start_lost               = 1L,
  frameshift_variant       = 1L,
  splice_acceptor_variant  = 2L,
  splice_donor_variant     = 2L,
  missense_variant         = 3L,
  inframe_insertion        = 3L,
  inframe_deletion         = 3L,
  splice_region_variant    = 4L,
  synonymous_variant       = 5L,
  stop_retained_variant    = 5L,
  `5_prime_UTR_variant`    = 6L,
  `3_prime_UTR_variant`    = 6L,
  intron_variant           = 7L,
  upstream_gene_variant    = 8L,
  downstream_gene_variant  = 8L,
  intergenic_variant       = 9L
)
SO_RANK_OTHER <- max(SO_TERM_RANK) + 1L

# Ensembl VEP IMPACT ladder, used only when the variant table carries an `impact`
# column (the Ensembl `overlap/region` REST endpoint does not provide it).
IMPACT_RANK <- c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)

# KASP acceptance thresholds (consensus + EasyKASP 2025; the manual mirrors this)
KASP_LIMITS <- list(
  amplicon     = c(50, 100),    # bp  (EasyKASP: products < 100 bp)
  as_tm        = c(55, 65),     # C
  com_tm       = c(58, 63),     # C  (EasyKASP aims ~59.9-62.9)
  gc           = c(0.40, 0.60), # fraction
  max_run      = 4L,            # longest homopolymer tolerated is max_run - 1
  hairpin_dg   = -3,            # dG37 (kcal/mol); must be >= this
  duplex_dg    = -9,            # cross-dimer dG37; must be >= this
  homodimer_dg = -9,            # self-dimer dG37; must be >= this
  discrimination = 0.5          # Taq on- minus off-target efficiency
)

#' OligoArrayAux duplex / hairpin free energy
#'
#' The DNA engine DECIPHER already uses: `hybrid-min` (duplex of two strands) or
#' `hybrid-ss-min` (intramolecular folding of one strand).  Returns dG in kcal/mol
#' (more negative = more stable); `NA` if OligoArrayAux is not on `PATH`.  A
#' non-zero exit or an unparsable result is reported as a warning (not a silent
#' `NA`), and the run goes through `processx` so failures are explicit.
#'
#' @param seq1,seq2 oligos 5'->3'; with `seq2 = NULL` the hairpin (folding) dG of
#'   `seq1` alone is returned.
#' @param na sodium concentration (M).
#' @param temp temperature (C).
#' @family KASP thermodynamics
#' @export
oligo_dg <- function(seq1, seq2 = NULL, na = 0.05, temp = 37) {
  bin <- if (is.null(seq2)) "hybrid-ss-min" else "hybrid-min"
  exe <- Sys.which(bin)
  if (!nzchar(exe)) {
    rlang::warn(glue::glue("OligoArrayAux ('{bin}') is not on the PATH."))
    return(NA_real_)
  }
  args <- c("-n", "DNA", "-t", as.character(temp), "-T", as.character(temp),
            "-N", as.character(na), "-E", "-q")
  args <- if (is.null(seq2)) c(args, seq1) else c(args, seq1, seq2)
  # -q makes UNAFold print a single number; run via processx so a non-zero exit
  # or an unparsable result is an explicit failure, not a silent NA.
  run <- tryCatch(processx::run(exe, args, error_on_status = FALSE),
                  error = function(e) e)
  if (inherits(run, "error")) {
    rlang::warn(glue::glue("OligoArrayAux ('{bin}') could not be run: {conditionMessage(run)}"))
    return(NA_real_)
  }
  if (run$status != 0L) {
    rlang::warn(glue::glue("OligoArrayAux ('{bin}') exited with status {run$status}: ",
                           "{trimws(run$stderr)}"))
    return(NA_real_)
  }
  tok <- utils::tail(strsplit(trimws(run$stdout), "[[:space:]]+")[[1L]], 1L)
  dg  <- suppressWarnings(as.numeric(tok))
  if (length(dg) != 1L || !is.finite(dg)) {
    rlang::warn(glue::glue("OligoArrayAux ('{bin}') returned no parsable dG (got '{tok}')."))
    return(NA_real_)
  }
  dg
}


# --- small utilities ---------------------------------------------------------

#' Normalise chromosome names to "Chr1".."Chr5"/"ChrM"/"ChrC"
#' @keywords internal
#' @noRd
ath_chrom <- function(x) {
  key <- sub("^(chromosome|chr)", "", x, ignore.case = TRUE)
  out <- unname(ATH_CHROM_MAP[key])
  ifelse(is.na(out), paste0("Chr", key), out)
}

#' Read one genomic region as a DNAString
#' @keywords internal
#' @noRd
get_region <- function(genome, chrom, start, end) genome_seq(genome, chrom, start, end)


# --- primer geometry (single source of truth) --------------------------------
# All intervals are plus-strand, with `snp_pos` the SNP coordinate.  The SAME
# helpers are reused by design (in template coordinates: pass `snp_i`), by the
# assay assembly and by the plots, so the geometry cannot drift between them.

#' Plus-strand interval occupied by a primer
#'
#' `role = "AS"` = allele-specific primer (3'-terminal base at the SNP);
#' `role = "COM"` = the common primer, `gap` bases away on the other side.
#' @param snp_pos plus-strand coordinate of the primer's 3'-terminal base.
#' @param strand binding strand, `"+"` or `"-"`.
#' @param length primer length in nt.
#' @param gap distance from `snp_pos` to the near edge of a common primer.
#' @param role `"AS"` (allele-specific) or `"COM"` (common).
#' @family KASP geometry
#' @export
primer_interval <- function(snp_pos, strand, length, gap = 0L, role = c("AS", "COM")) {
  role <- match.arg(role)
  if (identical(role, "AS")) {
    if (identical(strand, "+")) IRanges::IRanges(snp_pos - length + 1L, snp_pos)
    else                        IRanges::IRanges(snp_pos, snp_pos + length - 1L)
  } else {
    if (identical(strand, "+")) IRanges::IRanges(snp_pos + gap, snp_pos + gap + length - 1L)
    else                        IRanges::IRanges(snp_pos - gap - length + 1L, snp_pos - gap)
  }
}

#' Plus-strand interval of the AS primer's N 3'-terminal bases, SNP excluded
#'
#' `window` counts the 3'-terminal bases including the SNP itself, so the
#' returned interval has `window - 1` positions.
#' @param snp_pos coordinate of the variant's anchored base.
#' @param strand binding strand.
#' @param window number of 3'-terminal bases counted (includes the variant).
#' @family KASP geometry
#' @export
as_3p_interval <- function(snp_pos, strand, window = 5L) {
  if (identical(strand, "+")) IRanges::IRanges(snp_pos - window + 1L, snp_pos - 1L)
  else                        IRanges::IRanges(snp_pos + 1L, snp_pos + window - 1L)
}

#' Plus-strand coordinate of the AS primer's 3'-terminal base for a variant
#'
#' A variant is a substitution block `[v_start, v_start + ref_len - 1]` (a SNP is
#' the `ref_len = 1` case).  On the plus strand the AS primer's 3' end sits on
#' the block's last base; on the minus strand (where the primer is the reverse
#' complement) it sits on the block's first base.  All primer geometry is then
#' expressed through this single anchor, so SNP and MNP share one code path.
#' @param v_start plus-strand start of the variant block.
#' @param ref_len length of the REF allele (block size; 1 for a SNP).
#' @param strand binding strand.
#' @family KASP geometry
#' @export
as_anchor <- function(v_start, ref_len, strand) {
  if (identical(strand, "+")) v_start + ref_len - 1L else v_start
}

#' Index range of the variant block inside an AS-primer binding window
#'
#' For `strand = "+"` the block is at the 3' end of the window; for `"-"`
#' (reverse-complementary primer) it is at the 5' end.
#' @keywords internal
#' @noRd
block_in_seg <- function(ref_len, primer_len, strand) {
  if (identical(strand, "+")) IRanges::IRanges(primer_len - ref_len + 1L, primer_len)
  else                        IRanges::IRanges(1L, ref_len)
}

#' Count variant positions inside an interval (IRanges)
#'
#' `exclude` drops positions (typically the target SNP) from the count.
#' @param interval an `IRanges` interval.
#' @param vpos integer vector of variant positions (or `NULL`).
#' @param exclude positions to drop (typically the target variant).
#' @family KASP geometry
#' @export
count_variants_in <- function(interval, vpos, exclude = NA_integer_) {
  if (is.null(vpos)) return(NA_integer_)
  v <- vpos[!is.na(vpos) & !(vpos %in% exclude)]
  if (!length(v)) return(0L)
  as.integer(length(IRanges::findOverlaps(IRanges::IRanges(v, v), interval)))
}

#' Per-user cache directory for downloads and indices
#'
#' Everything heavy that the package fetches or builds (Ensembl FASTA/JSON, BLAST
#' and FM indices) lives here — in the OS user cache, not in the project tree — so
#' the working directory stays clean and the cache is shared across projects.
#' @return the cache directory path (created if needed).
#' @family KASP input
#' @export
kasp_cache_dir <- function() {
  d <- tools::R_user_dir("tidykasp", which = "cache")
  fs::dir_create(d)
  d
}

#' Download (once) and load the TAIR10 genome as a DNAStringSet
#'
#' The FASTA comes from Ensembl Plants and is cached in the per-user cache
#' directory (`tools::R_user_dir("tidykasp", "cache")`), not in the project.
#' Sequence names are normalised to Chr1..Chr5 / ChrM / ChrC (scaffolds
#' dropped), so they line up with the Ensembl gene records used elsewhere.
#' @param cache cache directory for the downloaded FASTA.
#' @param release Ensembl Plants release number.
#' @section Caching:
#' The download (~35 MB) happens once and is reused; delete the cached
#' `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz` to refresh.
#' @family KASP input
#' @export
ath_genome <- function(cache = kasp_cache_dir(), release = 63L) {
  fs::dir_create(cache)
  fa <- fs::path(cache, "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
  if (!fs::file_exists(fa)) {
    url <- sprintf(paste0(
      "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-%d/fasta/",
      "arabidopsis_thaliana/dna/",
      "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"), release)
    cli::cli_alert_info("Downloading TAIR10 genome (~35 MB) from Ensembl Plants into {.path {cache}} ...")
    utils::download.file(url, fa, mode = "wb", quiet = TRUE)
  }
  g <- Biostrings::readDNAStringSet(fa)
  nm <- sub("\\s.*$", "", names(g))   # keep only the first token of the defline
  dropped <- nm[!nm %in% names(ATH_CHROM_MAP)]
  if (length(dropped)) {
    cli::cli_alert_info(
      "Dropping {length(dropped)} non-chromosome sequence(s), e.g. {.val {utils::head(dropped, 3L)}}.")
  }
  g <- g[nm %in% names(ATH_CHROM_MAP)]
  names(g) <- unname(ATH_CHROM_MAP[nm[nm %in% names(ATH_CHROM_MAP)]])
  need <- c(paste0("Chr", 1:5), "ChrM", "ChrC")
  missing <- setdiff(need, names(g))
  if (length(missing)) {
    rlang::abort(glue::glue(
      "ath_genome(): TAIR10 FASTA is missing expected chromosome(s): ",
      "{paste(missing, collapse = ', ')} - Ensembl deflines may have changed."))
  }
  g[need]
}


# =============================================================================
#  1. Sequence retrieval
# =============================================================================

#' Fetch a gene record from Ensembl Plants REST (cached on disk)
#' @param gene_id Ensembl gene id (e.g. `"AT4G00650"`).
#' @param species Ensembl species name.
#' @param cache cache directory.
#' @family KASP input
#' @export
ensembl_gene <- function(gene_id, species = "arabidopsis_thaliana",
                         cache = kasp_cache_dir()) {
  fs::dir_create(cache)
  path <- fs::path(cache, glue::glue("{species}_{gene_id}.json"))
  if (!fs::file_exists(path)) {
    resp <- httr2::request("https://rest.ensembl.org") |>
      httr2::req_url_path_append("lookup", "id", gene_id) |>
      httr2::req_url_query(expand = 1) |>
      httr2::req_headers(`Content-Type` = "application/json") |>
      httr2::req_user_agent("tidykasp") |>
      httr2::req_retry(max_tries = 5L, backoff = ~ 3) |>
      httr2::req_perform()
    readr::write_lines(httr2::resp_body_string(resp), path)
  }
  j <- jsonlite::fromJSON(readr::read_file(path), simplifyDataFrame = FALSE)
  need <- c("id", "seq_region_name", "start", "end", "strand")
  missing <- setdiff(need, names(j))
  if (length(missing)) {
    rlang::abort(glue::glue(
      "ensembl_gene(): record is missing field(s) {paste(missing, collapse = ', ')} ",
      "- the Ensembl REST schema may have changed."))
  }
  tibble::tibble(gene_id = j$id, symbol = j$display_name %||% j$id,
                 chrom = ath_chrom(j$seq_region_name),
                 start = as.integer(j$start), end = as.integer(j$end),
                 strand = as.integer(j$strand),
                 description = j$description %||% NA_character_)
}

#' Natural variants overlapping a region, from Ensembl Plants (1001 Genomes)
#'
#' Returns `id`, `chrom`, `pos`, `ref_genome` (the allele vector) and
#' `consequence` — the Ensembl/Sequence-Ontology term(s).  The `overlap/region`
#' endpoint does **not** provide a VEP `impact` column, so [pick_target_snp()]
#' ranks on the SO terms; a table you build yourself with an `impact` column is
#' also accepted.
#' @param chrom,start,end region to fetch (chromosome is normalised).
#' @param species,cache as in [ensembl_gene()].
#' @family KASP input
#' @export
ensembl_variants <- function(chrom, start, end,
                             species = "arabidopsis_thaliana",
                             cache = kasp_cache_dir()) {
  fs::dir_create(cache)
  region <- glue::glue("{sub('^Chr', '', chrom)}:{start}-{end}")
  path <- fs::path(cache, glue::glue("{species}_var_{region}.json"))
  if (!fs::file_exists(path)) {
    resp <- httr2::request("https://rest.ensembl.org") |>
      httr2::req_url_path_append("overlap", "region", species, region) |>
      httr2::req_url_query(feature = "variation") |>
      httr2::req_headers(`Content-Type` = "application/json") |>
      httr2::req_user_agent("tidykasp") |>
      httr2::req_retry(max_tries = 5L, backoff = ~ 3) |>
      httr2::req_perform()
    readr::write_lines(httr2::resp_body_string(resp), path)
  }
  df <- jsonlite::fromJSON(readr::read_file(path)) |> tibble::as_tibble()
  need <- c("id", "seq_region_name", "start", "alleles")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    rlang::abort(glue::glue(
      "ensembl_variants(): response is missing field(s) {paste(missing, collapse = ', ')} ",
      "- the Ensembl REST schema may have changed."))
  }
  if (!"consequence_type" %in% names(df)) df$consequence_type <- NA_character_
  df |>
    dplyr::transmute(id = .data$id, chrom = ath_chrom(.data$seq_region_name),
                     pos = as.integer(.data$start), ref_genome = .data$alleles,
                     consequence = .data$consequence_type)
}

#' Extract a DNA window centred on a SNP as a DNAString
#'
#' The window is copied into a standalone `DNAString`: `subseq()` of a
#' `DNAStringSet` keeps a reference to the genome-wide shared pool, which would
#' otherwise bloat anything we serialise (e.g. the saved assay).
#' @param genome DNAStringSet.
#' @param chrom chromosome name.
#' @param pos variant position.
#' @param flank window half-width (bp).
#' @family KASP input
#' @export
snp_template <- function(genome, chrom, pos, flank = 150L) {
  seq <- get_region(genome, chrom, pos - flank, pos + flank)
  seq <- Biostrings::DNAString(as.character(seq))
  list(seq = seq, snp_i = as.integer(flank) + 1L)
}

#' GC fraction of the flanks on each side of the SNP (EasyKASP direction hint)
#'
#' EasyKASP chooses the read direction by flank GC: it prefers the side whose
#' +/-`flank` bases sit closest to 45%.  Returned here as a diagnostic.
#' @param genome DNAStringSet.
#' @param chrom chromosome name.
#' @param pos variant position.
#' @param flank bases to look at on each side.
#' @family KASP design
#' @export
strand_gc_hint <- function(genome, chrom, pos, flank = 25L) {
  up <- as.character(get_region(genome, chrom, max(1L, pos - flank), pos - 1L))
  dn <- as.character(get_region(genome, chrom, pos + 1L, pos + flank))
  tibble::tibble(side = c("upstream", "downstream"), gc = c(gc_frac(up), gc_frac(dn))) |>
    dplyr::mutate(distance_to_45 = abs(.data$gc - 0.45))
}

#' Map Ensembl consequence/impact to a design priority (lower = preferred)
#'
#' If an Ensembl VEP `impact` is available it wins (the standard HIGH/MODERATE/
#' LOW/MODIFIER severity ladder in `IMPACT_RANK`).  Otherwise the consequence
#' string is split on commas and matched against the canonical SO terms in
#' `SO_TERM_RANK` by *full token* (not substring), taking the most severe; terms
#' with no known rank (and `NA`s) fall back to `SO_RANK_OTHER`.
#' @keywords internal
#' @noRd
so_rank <- function(consequence, impact = NULL) {
  term_rank <- purrr::map_int(consequence, function(x) {
    toks <- stringr::str_trim(stringr::str_split(x, ",")[[1L]])
    r    <- SO_TERM_RANK[toks]
    r    <- r[!is.na(r)]
    if (length(r)) min(unname(r)) else SO_RANK_OTHER
  })
  if (is.null(impact)) return(term_rank)
  impact_rank <- unname(IMPACT_RANK[match(toupper(impact), names(IMPACT_RANK))])
  dplyr::coalesce(impact_rank, term_rank)
}

#' 3'-seed sequence of a would-be AS primer anchored at a variant
#'
#' `NA` if the seed would run off the sequence.
#' @keywords internal
#' @noRd
seed_seq <- function(genome, chrom, pos, ref_len, strand, k = SEED_K) {
  n <- genome_length(genome, chrom)
  anchor <- as_anchor(pos, ref_len, strand)
  lo <- if (identical(strand, "+")) anchor - k + 1L else anchor
  hi <- if (identical(strand, "+")) anchor else anchor + k - 1L
  if (lo < 1L || hi > n) return(NA_character_)
  seg <- get_region(genome, chrom, lo, hi)
  if (identical(strand, "+")) as.character(seg)
  else as.character(Biostrings::reverseComplement(seg))
}

#' 3'-seed specificity of a would-be AS primer anchored at a variant
#'
#' The seed is the `k` 3'-terminal bases of the AS primer: `[anchor-k+1, anchor]`
#' on the plus strand, or the reverse complement of `[anchor, anchor+k-1]` on the
#' minus strand (see [as_anchor()]).  Counted over the whole genome, both
#' strands, like [specificity_seed()].  `NA` if the seed runs off the sequence.
#' @keywords internal
#' @noRd
seed_hits_at <- function(genome, chrom, pos, ref_len, strand,
                         k = SEED_K, max_mismatch = SEED_MISMATCH, BPPARAM = NULL) {
  s <- seed_seq(genome, chrom, pos, ref_len, strand, k)
  if (is.na(s)) return(NA_integer_)
  specificity_of(c(seed = s), genome, k = k, max_mismatch = max_mismatch,
                 BPPARAM = BPPARAM)$hits
}

#' Choose a designable target variant from a variant table
#'
#' Keeps variants we can actually design on: a biallelic substitution block whose
#' REF and ALT are ACGT-only strings of *equal* length (SNP or MNP).  Indels and
#' multi-allelic records are left for a dedicated engine and are not silently
#' forced into a substitution.
#'
#' @section Selection policy:
#' Among the designable substitutions the best one is picked by, in order:
#' 1. the cleanest AS 3' window (`n_3p`, best of the two orientations);
#' 2. consequence severity — the Ensembl VEP `impact` column if the table has it,
#'    otherwise the canonical SO terms in `SO_TERM_RANK` (`missense_variant`, ...);
#' 3. a transition (A>G / C>T) over a transversion;
#' 4. 3'-seed specificity (`seed_hits`, best of the two strands);
#' 5. the leftmost position.
#' SNPs with complementary alleles (A/T, C/G) are skipped — a 3'-terminal base
#' cannot tell them apart.  The genome-wide seed scan runs only on the best
#' `specificity_candidates` loci.
#'
#' @param variants variant table (Ensembl REST shape).
#' @param genome DNAStringSet.
#' @param gene one-row gene record.
#' @param variant_id force a specific variant by id.
#' @param clean_3p number of flanking bases of the AS primer that should be free
#'   of *other* variants.  Because the AS primer may be placed on either strand,
#'   the "clean" count is the best (smallest) of the two orientations: the chosen
#'   variant just has to be clean in *at least one* direction.
#' @param specificity_candidates how many top-ranked loci to check for 3'-seed
#'   specificity.  The seed scan is genome-wide, so only this many of the best
#'   candidates (by the other criteria) are checked; set `0` to skip it.
#' @param BPPARAM optional `BiocParallel` backend for the seed scan (split by
#'   chromosome, see [specificity_seed()]).
#' @return a one-row tibble for the chosen variant (`id`, `chrom`, `pos`,
#'   `consequence`, `ref`, `alt`, `ref_len`, `transition`, `n_3p`, `seed_hits`).
#' @family KASP design
#' @export
pick_target_snp <- function(variants, genome, gene, variant_id = NULL, clean_3p = 5L,
                            specificity_candidates = 10L, BPPARAM = NULL) {
  if (!is.null(variant_id)) variants <- dplyr::filter(variants, .data$id == variant_id)

  # REF length comes from the VCF (its first allele); the reference bases are
  # then read from the genome, which is the truth we design against.
  v <- variants |>
    dplyr::mutate(
      ref_len = purrr::map_int(.data$ref_genome, ~ nchar(.x[1L])),
      ref = purrr::pmap_chr(list(.data$chrom, .data$pos, .data$ref_len),
                            function(ch, p, n) as.character(get_region(genome, ch, p, p + n - 1L))),
      # exactly one non-reference allele of the SAME length, else NA (dropped
      # below): a silent "first of several" choice would pick an arbitrary ALT,
      # and unequal lengths are indels, handled separately.
      alt = purrr::map2_chr(.data$ref_genome, .data$ref, function(alleles, ref) {
        other <- setdiff(alleles, ref)
        if (length(other) == 1L && nchar(other) == nchar(ref)) other else NA_character_
      })
    ) |>
    dplyr::filter(!is.na(.data$alt),
                  stringr::str_detect(.data$ref, "^[ACGT]+$"),
                  stringr::str_detect(.data$alt, "^[ACGT]+$"),
                  # a SNP whose alleles are complementary (A/T, C/G) cannot be told
                  # apart by a 3'-terminal base (one allele's primer pairs the other
                  # allele's template), so KASP discrimination collapses - skip them
                  !(.data$ref_len == 1L &
                      paste0(pmin(.data$ref, .data$alt), pmax(.data$ref, .data$alt)) %in% c("AT", "CG")))
  if (nrow(v) == 0L) rlang::abort("No usable substitution variant (SNP/MNP) left after filtering.")

  # every variant of the gene (incl. indels): used to penalise targets whose AS
  # 3' neighbourhood is itself variable (breaks allele specificity across lines)
  allpos <- variants$pos[variants$chrom == gene$chrom]

  # Ensembl VEP IMPACT, when the table carries it, is the library-defined
  # severity ladder; otherwise so_rank() falls back to the SO terms.  (`v[[..]]`
  # is NULL when the column is absent.)
  impact <- v[["impact"]]

  ranked <- v |>
    dplyr::mutate(
      transition = .data$ref_len == 1L &
        paste0(.data$ref, .data$alt) %in% c("AG", "GA", "CT", "TC"),
      rank = so_rank(.data$consequence, impact),
      # best-case (either strand) count of other variants in the AS 3' window,
      # vectorised over the whole variant table (IRanges).  The windows sit
      # strictly outside the variant block, so the target cannot count itself.
      n_3p = pmin(
        IRanges::countOverlaps(
          IRanges::IRanges(.data$pos - clean_3p + 1L, .data$pos - 1L),
          IRanges::IRanges(allpos, allpos)),
        IRanges::countOverlaps(
          IRanges::IRanges(.data$pos + .data$ref_len, .data$pos + .data$ref_len + clean_3p - 2L),
          IRanges::IRanges(allpos, allpos)))
    ) |>
    dplyr::arrange(.data$n_3p, .data$rank, dplyr::desc(.data$transition), .data$pos)

  if (specificity_candidates <= 0L) {
    return(ranked |> dplyr::slice(1L) |> dplyr::select(-"rank", -"ref_genome"))
  }

  # among the best candidates, prefer the most specific AS 3'-seed (best of the
  # two strands).  The genome-wide scan is one [specificity_seed()] call over
  # both strands' seeds of the whole shortlist, so it can fan out by chromosome.
  short <- dplyr::slice_head(ranked, n = specificity_candidates) |>
    dplyr::mutate(candidate = dplyr::row_number())
  # one seed per candidate per strand; a stable key carries it through the scan
  seed_tbl <- tidyr::crossing(candidate = short$candidate, strand = c("+", "-")) |>
    dplyr::mutate(
      seed = purrr::pmap_chr(list(short$chrom[.data$candidate], short$pos[.data$candidate],
                                  short$ref_len[.data$candidate], .data$strand),
                             function(ch, p, rl, st) seed_seq(genome, ch, p, rl, st)),
      key = paste0(.data$candidate, .data$strand)) |>
    dplyr::filter(!is.na(.data$seed))
  if (nrow(seed_tbl) > 0L) {
    seed_tbl <- specificity_of(stats::setNames(seed_tbl$seed, seed_tbl$key), genome,
                               BPPARAM = BPPARAM) |>
      dplyr::select(key = "primer", "hits") |>
      dplyr::right_join(seed_tbl, by = "key")
  } else {
    seed_tbl$hits <- integer()
  }
  # seed specificity breaks ties among otherwise comparable loci (it does not
  # outrank the biological/chemical preferences above); best of the two strands
  short <- seed_tbl |>
    dplyr::summarise(seed_hits = min(.data$hits), .by = "candidate") |>
    dplyr::right_join(short, by = "candidate")

  short |>
    dplyr::arrange(.data$n_3p, .data$rank, dplyr::desc(.data$transition),
                   .data$seed_hits, .data$pos) |>
    dplyr::slice(1L) |>
    dplyr::select(-"rank", -"ref_genome", -"candidate")
}


# =============================================================================
#  2. Thermodynamics via packages
# =============================================================================

#' GC fraction and homopolymer-run test
#' @keywords internal
#' @noRd
gc_frac  <- function(seq) unname(Biostrings::letterFrequency(Biostrings::DNAStringSet(seq), "GC", as.prob = TRUE)[, 1])
has_run  <- function(seq, n = 4L) stringr::str_detect(seq, glue::glue("([ACGT])\\1{{{n - 1L},}}"))

#' Make a parallel worker able to dispatch on Biostrings objects
#'
#' A `SnowParam`/SOCK worker is a fresh R process.  It resolves the package's
#' closures, but unless `Biostrings` is loaded there, S4 dispatch on a serialized
#' `DNAStringSet` (e.g. `genome[[i]]`, `subseq()`) falls back to base and fails
#' with "this S4 class is not subsettable".  Loading the namespace in the worker
#' fixes the dispatch; fork-based backends ([BiocParallel::MulticoreParam()])
#' inherit loaded namespaces and are unaffected.  Called at the top of every
#' worker closure that touches the genome.
#' @keywords internal
#' @noRd
worker_biostrings <- function() requireNamespace("Biostrings", quietly = TRUE)

#' Tandem repeats (simple sequence repeats) in a sequence
#'
#' EasyKASP screens candidate primers for SSR repeat units, not just
#' homopolymers.  Returns, for each input, the repeated units of length `units`
#' occurring at least `min_copies` times in tandem.
#' @param seq character vector of DNA sequences.
#' @param units repeat-unit lengths to test.
#' @param min_copies tandem copies required for a hit.
#' @family KASP QC
#' @export
find_repeats <- function(seq, units = 2:4, min_copies = 3L) {
  purrr::map(seq, function(s) {
    hits <- character(0)
    for (u in units) {
      pat <- glue::glue("([ACGT]{{{u}}})\\1{{{min_copies - 1L},}}")
      m   <- stringr::str_match_all(s, pat)[[1L]]
      if (nrow(m)) hits <- c(hits, m[, 2L])
    }
    unique(hits)
  })
}

#' Does a sequence contain an SSR? (see [find_repeats()])
#' @param seq character vector of DNA sequences.
#' @param units repeat-unit lengths to test.
#' @param min_copies tandem copies required for a hit.
#' @family KASP QC
#' @export
has_repeat <- function(seq, units = 2:4, min_copies = 3L)
  lengths(find_repeats(seq, units, min_copies)) > 0L

#' Inverted repeats / palindromes in a primer (library screen)
#'
#' `Biostrings::findPalindromes()` finds reverse-complementary arms - the
#' structural signature of a hairpin.  This complements the thermodynamic
#' hairpin screen in [primer_thermo()] with a sequence-level flag that does not
#' need OligoArrayAux.  The defaults are deliberately strict (arms >= 4 bp, loop
#' <= 6 bp, no mismatches): a well-designed primer returns no hits, whereas a
#' repetitive or self-symmetric one returns several.
#'
#' @param seq character vector of DNA sequences.
#' @param min_armlength minimum arm length of a palindromic hit.
#' @param max_looplength maximum loop (spacer) length.
#' @param max_mismatch mismatches allowed in the arms.
#' @return a list, one tibble per input, with `start`, `end`, `width` of each
#'   palindromic view (empty tibble when none).
#' @family KASP QC
#' @export
find_palindromes <- function(seq, min_armlength = 4L, max_looplength = 6L, max_mismatch = 0L) {
  nm <- names(seq)
  hits <- purrr::map(as.character(seq), function(s) {
    h <- Biostrings::findPalindromes(Biostrings::DNAString(s),
                                     min.armlength = min_armlength,
                                     max.looplength = max_looplength,
                                     max.mismatch  = max_mismatch)
    tibble::tibble(start = as.integer(IRanges::start(h)),
                   end   = as.integer(IRanges::end(h)),
                   width = as.integer(IRanges::width(h)))
  })
  stats::setNames(hits, nm)
}

#' Does a primer contain an inverted repeat? (see [find_palindromes()])
#' @param seq character vector of DNA sequences.
#' @param min_armlength minimum arm length of a palindromic hit.
#' @param max_looplength maximum loop (spacer) length.
#' @param max_mismatch mismatches allowed in the arms.
#' @family KASP QC
#' @export
has_palindrome <- function(seq, min_armlength = 4L, max_looplength = 6L, max_mismatch = 0L)
  purrr::map_int(find_palindromes(seq, min_armlength, max_looplength, max_mismatch), nrow) > 0L

#' Nearest-neighbour Tm of oligos (perfect duplex), salt-corrected (TmCalculator)
#'
#' `res$gr` uses TmCalculator's `$` accessor for the result GRanges; `Tm` is a
#' metadata column on it.
#' @param seq oligo(s) 5'->3'.
#' @family KASP thermodynamics
#' @export
primer_tm <- function(seq) {
  gr  <- TmCalculator::to_genomic_ranges(as.character(seq))
  res <- TmCalculator::tm_calculate(gr, method = "tm_nn", nn_table = NN_TABLE,
                                    Na = PCR_COND$Na, Mg = PCR_COND$Mg,
                                    dNTPs = PCR_COND$dNTPs)
  as.numeric(S4Vectors::mcols(res$gr)$Tm)
}

#' Is OligoArrayAux (`hybrid-min`) available on the PATH?
#'
#' DECIPHER's thermodynamics runs through OligoArrayAux.
#'
#' @family KASP utilities
#' @export
has_oligoarrayaux <- function() nzchar(Sys.which("hybrid-min"))

#' Point OligoArrayAux at its data tables for the duration of `code`
#'
#' The thermodynamics engine (DECIPHER via OligoArrayAux) needs `UNAFOLDDAT` to
#' find its parameter tables.  This temporarily sets it, so you need not fiddle
#' with the environment by hand.
#'
#' @param path directory with the OligoArrayAux data tables.
#' @param code expression to evaluate (e.g. `assemble_kasp_assay(...)`).
#' @family KASP utilities
#' @export
with_unafloddat <- function(path, code) {
  if (!fs::dir_exists(path)) rlang::abort(glue::glue("Directory not found: {path}"))
  withr::with_envvar(c(UNAFOLDDAT = path), code)
}

#' PCR amplification efficiency of primer/target pairs (DECIPHER)
#'
#' Models hybridisation (incl. side reactions) times Taq elongation efficiency,
#' so a 3'-terminal mismatch is strongly penalised - the physical basis of KASP.
#'
#' @param primers character vector, 5'->3'.
#' @param targets character vector (same length) with the reverse-complement of
#'   the binding site (as in `?DECIPHER::CalculateEfficiencyPCR`).
#' @param temp annealing temperature (C).
#' @param P primer concentration (M).
#' @param ions monovalent + divalent salt (M).
#' @family KASP thermodynamics
#' @export
pcr_efficiency <- function(primers, targets, temp = PCR_EFF$temp,
                           P = PCR_EFF$P, ions = PCR_EFF$ions) {
  rlang::check_installed("DECIPHER")
  if (!has_oligoarrayaux()) {
    rlang::abort(c("OligoArrayAux ('hybrid-min') is required for PCR efficiencies.",
                   i = "conda install -c bioconda oligoarrayaux"))
  }
  as.numeric(DECIPHER::CalculateEfficiencyPCR(
    primer = Biostrings::DNAStringSet(as.character(primers)),
    target = Biostrings::DNAStringSet(as.character(targets)),
    temp = temp, P = P, ions = ions, taqEfficiency = TRUE))
}


# =============================================================================
#  3. Common (reverse) primer: a small library-scored candidate scan
# =============================================================================

#' Design candidate common primers
#'
#' Enumerates the partner primer and keeps the ones whose length, GC, Tm and
#' complexity look like a good KASP common primer.  `strand` is the orientation
#' of the *allele-specific* primer: the common primer must sit on the opposite
#' side of the SNP and point back at it, otherwise the amplicon would not close.
#'
#'   * `strand = "+"`: AS extends downstream, so the common primer is a reverse
#'     primer *downstream* of the SNP (sequence = reverse complement).
#'   * `strand = "-"`: AS extends upstream, so the common primer is a forward
#'     primer *upstream* of the SNP (sequence = the plus strand itself).
#'
#' Metrics come from packages (TmCalculator, Biostrings); enumeration/filtering
#' is ours.
#' @param window plus-strand template window (from [snp_template()]).
#' @param snp_i template coordinate of the variant anchor.
#' @param strand assay orientation.
#' @param lengths candidate primer lengths.
#' @param gaps distance from the anchor to the common primer.
#' @param target_tm preferred common-primer Tm.
#' @param min_gc,max_gc GC bounds.
#' @param tm_range accepted Tm window (C).
#' @param n number of candidates to return.
#' @family KASP design
#' @export
design_common_primers <- function(window, snp_i, strand = "+", lengths = 18:25,
                                  gaps = 18:80, target_tm = COM_TM_TARGET,
                                  min_gc = 0.40, max_gc = 0.60,
                                  tm_range = c(57, 64), n = 5L) {
  wlen <- length(window)
  plus <- identical(strand, "+")
  grid <- tidyr::crossing(L = lengths, gap = gaps)

  cand <- purrr::pmap(grid, function(L, gap) {
    iv <- primer_interval(snp_i, strand, L, gap = gap, role = "COM")  # template coords
    start <- IRanges::start(iv); end <- IRanges::end(iv)
    if (start < 1L || end > wlen) return(NULL)
    seq_plus <- Biostrings::subseq(window, start, end)
    com <- if (plus) as.character(Biostrings::reverseComplement(seq_plus))
           else      as.character(seq_plus)
    tibble::tibble(com = com, len = L, gap = gap, start = start, end = end)
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(com_gc = gc_frac(.data$com),
                  com_tm = primer_tm(.data$com),
                  run    = has_run(.data$com),
                  ssr    = has_repeat(.data$com)) |>
    dplyr::filter(!.data$run, !.data$ssr,
                  dplyr::between(.data$com_gc, min_gc, max_gc),
                  dplyr::between(.data$com_tm, tm_range[1], tm_range[2])) |>
    dplyr::arrange(abs(.data$com_tm - target_tm))

  if (nrow(cand) == 0L) {
    rlang::abort(c("No common primer passed the filters.",
                   i = "Widen `tm_range`/GC or increase `flank` in snp_template()."))
  }
  cand |> dplyr::slice_head(n = n) |> dplyr::mutate(rank = dplyr::row_number())
}


# =============================================================================
#  4. Allele-specific primers (the KASP-specific part)
# =============================================================================

#' Plus-strand binding window of an AS primer ending at `anchor`
#' @keywords internal
#' @noRd
as_window <- function(window, anchor, length, strand) {
  iv <- primer_interval(anchor, strand, length, role = "AS")
  Biostrings::subseq(window, IRanges::start(iv), IRanges::end(iv))
}

#' DECIPHER target (= strand bound by the primer) for a binding window
#' @keywords internal
#' @noRd
as_target <- function(seg, strand) {
  if (identical(strand, "+")) Biostrings::reverseComplement(seg) else seg
}

#' Build one allele-specific primer on a given strand
#'
#' `anchor` is the plus-strand coordinate of the primer's 3'-terminal base (see
#' [as_anchor()]).  The allele's block (1 nt for a SNP, several nt for an MNP) is
#' inserted at the 3' end of the binding window, so the primer ends on the
#' variant.  If `stab_at` is given, the base `stab_at` positions from the 3' end
#' is replaced by `stab_base` (the deliberate stabilising mismatch).
#'
#' `strand = "+"`: primer matches the plus strand, 3' base at the anchor.
#' `strand = "-"`: primer is the reverse complement, 3' base at the anchor.
#' @keywords internal
#' @noRd
make_as_primer <- function(window, anchor, length, allele, strand = "+",
                           stab_at = NULL, stab_base = NULL) {
  seg <- as_window(window, anchor, length, strand)       # DNAString, contains the variant
  blk <- block_in_seg(nchar(allele), length, strand)     # variant block inside `seg`
  p   <- Biostrings::replaceAt(seg, blk, allele)
  if (!identical(strand, "+")) p <- Biostrings::reverseComplement(p)
  if (!is.null(stab_at)) {
    i <- length - stab_at + 1L                            # 3'-relative -> 5'-index
    p <- Biostrings::replaceAt(p, IRanges::IRanges(i, i), stab_base)
  }
  as.character(p)
}

#' Grid-search allele-specific primers maximising Taq-based discrimination
#'
#' For every strand, length, secondary-mismatch position and mismatch base we
#' build the REF and ALT allele-anchored primers and predict the PCR efficiency of
#' each on both alleles (DECIPHER).  Discrimination = on-target minus off-target
#' efficiency, averaged over the two primers.  `ref`/`alt` are strings, so the
#' same engine handles SNP (`nchar == 1`) and substitution blocks / MNP (longer).
#'
#' @param v_start plus-strand start of the variant block (a single base for a SNP).
#' @param ref,alt the two alleles as ACGT strings of equal length.
#' @param strands which strands to consider for the AS primer; `"+"` = forward
#'   primer with the variant as its 3'-terminal block extending downstream, `"-"` =
#'   the mirror-image primer extending upstream.  Both are valid KASP layouts.
#' @examples
#' W <- Biostrings::DNAString(paste(rep("ACGT", 20L), collapse = ""))
#' \dontrun{
#' design_as_primers(W, 40L, ref = "T", alt = "A", lengths = 20:22)
#' }
#' @param window plus-strand template window.
#' @param lengths AS-primer lengths to try.
#' @param stab_ats positions (from the 3' end) for the stabilising mismatch;
#'   `NA` = none.
#' @param temp,P,ions as in [pcr_efficiency()].
#' @family KASP design
#' @export
design_as_primers <- function(window, v_start, ref, alt, strands = c("+", "-"),
                              lengths = 20:27, stab_ats = c(NA, 2L, 3L),
                              temp = PCR_EFF$temp, P = PCR_EFF$P, ions = PCR_EFF$ions) {
  ref_len <- nchar(ref)
  # candidate mismatch bases at `k` bases from the 3' end: any base other than
  # the one the (perfectly matching) primer already carries there.
  mismatch_bases <- function(L, k, strand) {
    base_p <- make_as_primer(window, as_anchor(v_start, ref_len, strand), L, ref,
                             strand = strand)
    cur    <- substr(base_p, L - k + 1L, L - k + 1L)
    setdiff(c("A", "C", "G", "T"), cur)
  }

  # enumerate configurations (strand x length x mismatch position x mismatch base)
  cfgs <- tidyr::crossing(strand = strands, as_length = lengths, stab_at = stab_ats) |>
    dplyr::mutate(stab_base = purrr::pmap(
      list(.data$as_length, .data$stab_at, .data$strand),
      function(L, k, st) if (is.na(k)) NA_character_ else mismatch_bases(L, k, st))) |>
    tidyr::unnest_longer("stab_base")

  # build all primer/target pairs and score them in ONE batched DECIPHER call
  pairs <- purrr::pmap(cfgs, function(strand, as_length, stab_at, stab_base) {
    L <- as_length
    k <- if (is.na(stab_at))   NULL else stab_at
    b <- if (is.na(stab_base)) NULL else stab_base
    anchor <- as_anchor(v_start, ref_len, strand)
    # binding window on the plus strand, and the same window carrying the ALT allele
    seg_ref <- as_window(window, anchor, L, strand)
    seg_alt <- Biostrings::replaceAt(seg_ref, block_in_seg(ref_len, L, strand), alt)
    ref_p <- make_as_primer(window, anchor, L, ref, strand = strand, stab_at = k, stab_base = b)
    alt_p <- make_as_primer(window, anchor, L, alt, strand = strand, stab_at = k, stab_base = b)
    t_ref <- as.character(as_target(seg_ref, strand))
    t_alt <- as.character(as_target(seg_alt, strand))
    # explicit on/off pairing: the "on" target carries the primer's own allele,
    # the "off" target carries the other one.  Written out (not c()-positional)
    # so a swap cannot silently invert the discrimination metric.
    specs <- tibble::tibble(
      which  = c("ref_on", "ref_off", "alt_on", "alt_off"),
      primer = c(ref_p,    ref_p,     alt_p,    alt_p),
      target = c(t_ref,    t_alt,     t_alt,    t_ref))
    dplyr::mutate(specs, strand = strand, as_length = L,
                  stab_at = stab_at, stab_base = stab_base)
  }) |>
    purrr::list_rbind()
  pairs$eff <- pcr_efficiency(pairs$primer, pairs$target, temp = temp, P = P, ions = ions)

  wide <- pairs |>
    dplyr::select("strand","as_length","stab_at","stab_base","which","eff") |>
    tidyr::pivot_wider(names_from = "which", values_from = "eff") |>
    dplyr::mutate(discrimination = ((.data$ref_on - .data$ref_off) +
                                    (.data$alt_on - .data$alt_off)) / 2)

  # attach sequences/metrics for the surviving configurations
  wide |>
    dplyr::rowwise() |>
    dplyr::mutate(
      ref_primer = make_as_primer(window, as_anchor(v_start, ref_len, .data$strand),
                                  .data$as_length, ref, strand = .data$strand,
                                  if (is.na(.data$stab_at)) NULL else .data$stab_at,
                                  if (is.na(.data$stab_base)) NULL else .data$stab_base),
      alt_primer = make_as_primer(window, as_anchor(v_start, ref_len, .data$strand),
                                  .data$as_length, alt, strand = .data$strand,
                                  if (is.na(.data$stab_at)) NULL else .data$stab_at,
                                  if (is.na(.data$stab_base)) NULL else .data$stab_base),
      # template-local footprint of the AS primer (genomic coords are added by
      # the caller, which knows the window offset)
      tmpl_start = IRanges::start(primer_interval(as_anchor(v_start, ref_len, .data$strand),
                                                  .data$strand, .data$as_length, role = "AS")),
      tmpl_end   = IRanges::end(primer_interval(as_anchor(v_start, ref_len, .data$strand),
                                                .data$strand, .data$as_length, role = "AS"))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(gc = gc_frac(.data$ref_primer), tm = primer_tm(.data$ref_primer),
                  dg_3p = terminal_dg(.data$ref_primer)) |>
    dplyr::select("strand","as_length","stab_at","stab_base","ref_primer","alt_primer",
                  "ref_on","ref_off","alt_on","alt_off","discrimination","tm","gc",
                  "dg_3p","tmpl_start","tmpl_end") |>
    dplyr::arrange(dplyr::desc(.data$discrimination))
}


# =============================================================================
#  4b. Indels: allele templates + junction-scan allele-specific primers
# =============================================================================

#' Plus-strand REF and ALT allele templates for an indel
#'
#' Both templates share the left flank (up to `snp_i - 1`).  The REF template is
#' the genomic window; the ALT template replaces the REF block with the ALT
#' allele.  Coordinates left of the block are identical, and to the right they
#' differ by `delta = nchar(alt) - nchar(ref)`.
#' @param ref_tmpl plus-strand REF template.
#' @param snp_i template coordinate of the first REF base.
#' @param ref,alt alleles as ACGT strings.
#' @family KASP design
#' @export
allele_templates <- function(ref_tmpl, snp_i, ref, alt) {
  ref_len <- nchar(ref)
  left  <- if (snp_i > 1L)                  Biostrings::subseq(ref_tmpl, 1L, snp_i - 1L)               else Biostrings::DNAString("")
  right <- if (snp_i + ref_len <= length(ref_tmpl)) Biostrings::subseq(ref_tmpl, snp_i + ref_len, length(ref_tmpl)) else Biostrings::DNAString("")
  list(ref = ref_tmpl,
       alt = Biostrings::xscat(left, Biostrings::DNAString(alt), right),
       delta = nchar(alt) - ref_len,
       block_ref = Biostrings::DNAString(ref),
       block_alt = Biostrings::DNAString(alt))
}

#' Allele-specific primer candidates for an insertion/deletion (junction scan)
#'
#' For indels the two alleles differ in *length*, so a fixed-offset substitution
#' no longer applies.  Instead each primer is built from its own allele template
#' and slid along the junction: every candidate is a perfect match to its own
#' allele, and its PCR efficiency against the OTHER allele is scored with
#' DECIPHER (the same Taq-aware model used elsewhere).  The placement that
#' actually discriminates is therefore chosen by the physics, not by a hand rule.
#'
#' @param ref_tmpl plus-strand template (small design window).
#' @param snp_i template coordinate of the first REF base.
#' @param ref,alt alleles as ACGT strings with `nchar(ref) != nchar(alt)`.
#' @param lengths AS-primer lengths to try.
#' @return tidy tibble of candidates: `allele` (`ref`/`alt`), `strand`,
#'   `as_length`, `as_start`/`as_end` (this allele template's coords), `primer`,
#'   `on_eff`, `off_eff`, `discrimination`, `tm`, `gc`.
#' @param strands strands to scan.
#' @param temp,P,ions as in [pcr_efficiency()].
#' @family KASP design
#' @export
design_as_primers_indel <- function(ref_tmpl, snp_i, ref, alt,
                                    strands = c("+", "-"), lengths = 20:27,
                                    temp = PCR_EFF$temp, P = PCR_EFF$P, ions = PCR_EFF$ions) {
  ref_len <- nchar(ref); alt_len <- nchar(alt)
  if (ref_len == alt_len) {
    rlang::abort("design_as_primers_indel() needs unequal alleles; use design_as_primers() for a substitution.")
  }
  at <- allele_templates(ref_tmpl, snp_i, ref, alt)

  # The OTHER allele's version of a window on THIS template: the window's block
  # overlap is replaced by the other allele, flanks kept.  Windows that miss the
  # block are identical on both alleles (no discrimination).
  substitute_block <- function(tmpl, other_block, block_len, s, e) {
    block_end <- snp_i + block_len - 1L
    if (s > block_end || e < snp_i) return(Biostrings::subseq(tmpl, s, e))
    left  <- if (s < snp_i)     Biostrings::subseq(tmpl, s, snp_i - 1L)    else Biostrings::DNAString("")
    right <- if (e > block_end) Biostrings::subseq(tmpl, block_end + 1L, e) else Biostrings::DNAString("")
    Biostrings::xscat(left, other_block, right)
  }

  build <- function(allele, tmpl, other_block, block_len, strand, L, e) {
    if (identical(strand, "+")) { s <- e - L + 1L; ee <- e } else { s <- e; ee <- e + L - 1L }
    if (s < 1L || ee > length(tmpl)) return(NULL)
    seg <- Biostrings::subseq(tmpl, s, ee)
    primer <- if (identical(strand, "+")) as.character(seg)
              else as.character(Biostrings::reverseComplement(seg))
    off_seg <- substitute_block(tmpl, other_block, block_len, s, ee)
    tibble::tibble(allele = allele, strand = strand, as_length = L,
                   as_start = s, as_end = ee, primer = primer,
                   on = as.character(as_target(seg, strand)),
                   off = as.character(as_target(off_seg, strand)))
  }

  cand <- purrr::pmap_dfr(
    tidyr::crossing(allele = c("ref", "alt"), strand = strands, as_length = lengths),
    function(allele, strand, as_length) {
      tmpl  <- if (allele == "ref") at$ref else at$alt
      oblk  <- if (allele == "ref") at$block_alt else at$block_ref
      blen  <- if (allele == "ref") ref_len else alt_len
      # 3'-end placements: from one base before the block to one past its end
      ends <- seq.int(max(1L, snp_i - 1L), min(length(tmpl), snp_i + blen))
      purrr::map_dfr(ends, function(e) build(allele, tmpl, oblk, blen, strand, as_length, e))
    })
  if (nrow(cand) == 0L) rlang::abort("No indel AS-primer placement fit inside the design window.")

  cand$on_eff  <- pcr_efficiency(cand$primer, cand$on,  temp = temp, P = P, ions = ions)
  cand$off_eff <- pcr_efficiency(cand$primer, cand$off, temp = temp, P = P, ions = ions)
  cand |>
    dplyr::mutate(discrimination = .data$on_eff - .data$off_eff,
                  tm = primer_tm(.data$primer), gc = gc_frac(.data$primer)) |>
    dplyr::select("allele","strand","as_length","as_start","as_end","primer",
                  "on_eff","off_eff","discrimination","tm","gc") |>
    dplyr::arrange(dplyr::desc(.data$discrimination))
}


# =============================================================================
#  5. Structural QC : Tm, hairpins, self- and cross-complementarity
# =============================================================================

#' Full structural profile of one primer
#' @param primer oligo(s) 5'->3'.
#' @family KASP thermodynamics
#' @export
primer_thermo <- function(primer) {
  s <- as.character(primer)
  tibble::tibble(
    primer       = s,
    len          = nchar(s),
    tm           = primer_tm(s),
    gc           = gc_frac(s),
    hairpin_dg   = oligo_dg(s),        # hybrid-ss-min, kcal/mol
    terminal_dg  = terminal_dg(s))     # 3'-terminal n-mer duplex, kcal/mol
}

#' Duplex free energy of a primer's 3'-terminal n-mer (3'-end stability)
#'
#' The last `n` bases are annealed to their reverse complement and scored with
#' OligoArrayAux (kcal/mol, 37 C); more negative = a more stable 3' anchor, the
#' "GC clamp" a polymerase needs to initiate.  Used as a design diagnostic.
#' @param primer oligo(s) 5'->3'.
#' @param n terminal window length.
#' @family KASP thermodynamics
#' @export
terminal_dg <- function(primer, n = 5L) {
  vapply(as.character(primer), function(s) {
    if (nchar(s) < n) return(NA_real_)
    mer <- substr(s, nchar(s) - n + 1L, nchar(s))
    oligo_dg(mer, as.character(Biostrings::reverseComplement(Biostrings::DNAString(mer))))
  }, numeric(1), USE.NAMES = FALSE)
}

#' Is the 3'-terminal allelic mismatch a weak "wobble" pair (G:T / A:C)?
#'
#' `on_allele`/`off_allele` are the plus-strand alleles.  The primer carries the
#' "on" base and faces the complement of the "off" base, so the mismatch set is
#' `{on, complement(off)}`; G:T and A:C are the least destabilising and give the
#' poorest allelic discrimination (LGC/KBioscience note).  Reported as a flag.
#' @param on_allele,off_allele plus-strand alleles (the primer carries `on`, the
#'   other allele is `off`).
#' @family KASP thermodynamics
#' @export
wobble_mismatch <- function(on_allele, off_allele) {
  faced <- as.character(Biostrings::complement(Biostrings::DNAStringSet(toupper(off_allele))))
  pair  <- purrr::map2_chr(toupper(on_allele), faced,
                           ~ paste(sort(c(.x, .y)), collapse = ""))
  pair %in% c("GT", "AC")
}

#' Duplex (dimer) free energy for every primer pair (OligoArrayAux, DNA)
#'
#' `dg` is the hybridisation free energy of the two strands (kcal/mol, 37 C);
#' more negative = stronger dimer.  The diagonal is the homodimer of each primer.
#' Returns the FULL matrix (both orders); symmetry is imposed here, once.  Each
#' pair is one OligoArrayAux subprocess, so for large primer sets pass a
#' `BPPARAM` to run the n(n+1)/2 pairs in parallel.
#' @param primers named character vector of oligos.
#' @param BPPARAM optional `BiocParallel` backend; pairs are split across workers.
#' @family KASP thermodynamics
#' @export
dimer_thermo <- function(primers, BPPARAM = NULL) {
  nm <- names(primers)
  # upper triangle (i <= j) is all DISTINCT pairs; the full grid is then a join
  grid <- tidyr::crossing(a = nm, b = nm) |>
    dplyr::mutate(lo = pmin(match(.data$a, nm), match(.data$b, nm)),
                  hi = pmax(match(.data$a, nm), match(.data$b, nm)))
  keys <- dplyr::distinct(grid, .data$lo, .data$hi)
  one  <- function(k) oligo_dg(primers[[keys$lo[k]]], primers[[keys$hi[k]]])
  keys$dg <- if (is.null(BPPARAM)) purrr::map_dbl(seq_len(nrow(keys)), one)
             else unlist(BiocParallel::bplapply(seq_len(nrow(keys)), one, BPPARAM = BPPARAM))
  grid |>
    dplyr::left_join(keys, by = c("lo", "hi")) |>
    dplyr::transmute(a = .data$a, b = .data$b, dg = .data$dg)
}


# =============================================================================
#  6. Genome-wide specificity (pure R, plus optional BLAST)
# =============================================================================

#' Count genome-wide near-matches of each primer's 3'-terminal k-mer
#'
#' For the exact case (`max_mismatch = 0`) counts come from `Biostrings::PDict`
#' (`vcountPDict`) - the precomputed-dictionary fast path, one call for all
#' primers.  With mismatches Biostrings has no counting index (`PDict()` builds an
#' `MTB_PDict`, but `vcountPDict()` does not support it), so the scan is done with
#' `vcountPattern` **parallelised over the genome's sequences** (and, per
#' sequence, once for all seeds), which is how this scales to large genomes.
#' Results are summed over both strands and over all genome sequences.
#' @param primers named character vector of oligos.
#' @param genome DNAStringSet.
#' @param k 3'-seed length.
#' @param max_mismatch mismatches allowed in the seed.
#' @param BPPARAM optional `BiocParallel` backend; the mismatch scan is split by
#'   genome sequence, so `MulticoreParam()` uses one worker per chromosome.
#' @family KASP specificity
#' @export
specificity_seed <- function(primers, genome, k = SEED_K, max_mismatch = SEED_MISMATCH,
                             BPPARAM = NULL) {
  if (inherits(genome, "kasp_genome")) {
    rlang::abort(c("specificity_seed() needs an in-memory DNAStringSet.",
                   i = "For an indexed genome use specificity_of()/specificity_blast()."))
  }
  seeds <- unname(substr(primers, nchar(primers) - k + 1L, nchar(primers)))
  if (max_mismatch == 0L) {
    hits <- exact_seed_counts(seeds, genome)
  } else {
    # vcountPattern takes a single pattern, so loop over seeds *per chromosome*;
    # parallelism is over the genome's sequences (one worker per chromosome)
    fwd <- Biostrings::DNAStringSet(seeds)
    rev <- Biostrings::reverseComplement(fwd)
    count_chr <- function(i) {
      worker_biostrings()
      chr <- genome[[i]]
      cf <- vapply(seq_along(fwd), function(j)
        Biostrings::countPattern(fwd[[j]], chr, max.mismatch = max_mismatch), numeric(1))
      cr <- vapply(seq_along(rev), function(j)
        Biostrings::countPattern(rev[[j]], chr, max.mismatch = max_mismatch), numeric(1))
      cf + cr
    }
    per_chr <- if (is.null(BPPARAM)) lapply(seq_along(genome), count_chr)
               else BiocParallel::bplapply(seq_along(genome), count_chr, BPPARAM = BPPARAM)
    hits <- Reduce(`+`, per_chr)
  }
  tibble::tibble(primer = names(primers), seed = seeds, hits = as.integer(hits))
}

#' Exact both-strand occurrence counts of seeds via a `PDict` dictionary
#' @keywords internal
#' @noRd
exact_seed_counts <- function(seeds, genome) {
  rc   <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seeds)))
  pats <- c(seeds, rc)
  keep <- !duplicated(pats)
  # vcountPDict returns one row per dictionary element, in the given order
  m    <- Biostrings::vcountPDict(Biostrings::PDict(pats[keep]), genome)
  cnt  <- stats::setNames(as.integer(rowSums(m)), pats[keep])
  as.integer(cnt[seeds] + cnt[rc])
}

#' In-silico PCR: every product a primer pair would amplify
#'
#' Uses `Biostrings::matchProbePair` (the library's own virtual PCR) to find
#' each product the forward/reverse pair generates, chromosome by chromosome.
#' A KASP test must amplify exactly one product per allele.
#'
#' @param F_primer forward primer (5'->3'), binds the plus strand.
#' @param R_primer reverse primer (5'->3').
#' @param genome DNAStringSet.
#' @param max_mismatch mismatches tolerated in the primer match.
#' @param BPPARAM optional `BiocParallel` backend; the virtual PCR is split by
#'   genome sequence (one worker per chromosome with `MulticoreParam()`).
#' @return tibble `chrom`, `start`, `end`, `width` (one row per product).
#' @family KASP specificity
#' @export
in_silico_pcr <- function(F_primer, R_primer, genome, max_mismatch = 0L, BPPARAM = NULL) {
  Fp <- Biostrings::DNAString(F_primer)
  Rp <- Biostrings::DNAString(R_primer)
  chroms <- genome_chroms(genome)
  scan_one <- function(nm) {
    worker_biostrings()
    # genome_seq() reads one sequence at a time, so an indexed `kasp_genome` stays
    # out of RAM; a DNAStringSet is simply subscripted
    pr <- Biostrings::matchProbePair(Fp, Rp, genome_seq(genome, nm),
                                     max.mismatch = max_mismatch, verbose = FALSE)
    if (is.null(pr) || length(pr) == 0L) return(NULL)
    tibble::tibble(chrom = nm, start = IRanges::start(pr), end = IRanges::end(pr),
                   width = IRanges::width(pr))
  }
  res <- if (is.null(BPPARAM)) lapply(chroms, scan_one)
         else BiocParallel::bplapply(chroms, scan_one, BPPARAM = BPPARAM)
  purrr::list_rbind(res)
}

#' Products predicted for each allele-specific primer against the common one
#'
#' On the reference genome the REF primer must give exactly one product; the
#' ALT primer (the other base / the indel) should give none.  Uses the
#' template-matched primers (no deliberate stabiliser), like [specificity_seed()].
#' @param assay assembled assay (from [assemble_kasp_assay()]).
#' @param genome DNAStringSet.
#' @param max_mismatch mismatches tolerated in the primer match.
#' @param BPPARAM optional `BiocParallel` backend (see [in_silico_pcr()]).
#' @family KASP specificity
#' @export
assay_amplicons <- function(assay, genome, max_mismatch = 0L, BPPARAM = NULL) {
  aa      <- assay$assay[1L, ]
  ref_len <- nchar(assay$snp$ref)
  # template-local anchor (assay$template$seq is the small design window)
  anchor  <- as_anchor(assay$template$snp_i, ref_len, aa$strand)
  seed <- c(AS_ref = make_as_primer(assay$template$seq, anchor,
                                    aa$as_len, assay$snp$ref, strand = aa$strand),
            AS_alt = make_as_primer(assay$template$seq, anchor,
                                    aa$as_len, assay$snp$alt, strand = aa$strand))
  purrr::imap_dfr(seed, function(as_p, nm) {
    fr <- if (identical(aa$strand, "+")) c(F = as_p, R = aa$com) else c(F = aa$com, R = as_p)
    prods <- in_silico_pcr(fr[["F"]], fr[["R"]], genome, max_mismatch = max_mismatch,
                           BPPARAM = BPPARAM)
    tibble::tibble(primer = nm, n_products = nrow(prods),
                   width = if (nrow(prods)) max(prods$width) else NA_integer_)
  })
}

#' BLAST every primer against a local genome database (rBLAST, optional)
#' @param primers named character vector of oligos.
#' @param db_path BLAST database path.
#' @family KASP specificity
#' @export
blast_primers <- function(primers, db_path) {
  if (is.null(db_path) || !fs::file_exists(db_path)) return(NULL)
  rlang::check_installed("rBLAST")
  query <- Biostrings::DNAStringSet(stats::setNames(primers, names(primers)))
  out <- tryCatch({
    bl <- rBLAST::blast(db_path, type = "blastn")
    stats::predict(bl, newdata = query,
                   BLAST_args = c("-task", "blastn-short", "-word_size", "7", "-dust", "no"))
  }, error = function(e) { rlang::warn(glue::glue("BLAST failed: {conditionMessage(e)}")); NULL })
  if (is.null(out) || nrow(out) == 0L) return(out)
  tibble::as_tibble(out)
}

#' Build (once) a BLAST nucleotide database from the TAIR10 genome
#' @param genome DNAStringSet.
#' @param path FASTA path to write (and index).
#' @family KASP specificity
#' @export
build_blast_db <- function(genome, path = fs::path(kasp_cache_dir(), "tair10.fa")) {
  fs::dir_create(fs::path_dir(path))
  if (!fs::file_exists(path)) Biostrings::writeXStringSet(genome, path)
  rBLAST::makeblastdb(path, dbtype = "nucl", verbose = FALSE)
  path
}


# =============================================================================
#  7. Assemble + score one complete assay
# =============================================================================

#' Assemble and score the best KASP assay for a target SNP
#'
#' When `variants` is supplied the selection is *variant-aware*: among primers
#' that pass the physico-chemical filters we minimise the number of other
#' variants falling under the primer footprint (the 3'-terminal 5 nt are handled
#' earlier by [pick_target_snp()]), then maximise discrimination.
#'
#' @section Selection policy:
#' The AS placement is chosen lexicographically by:
#' 1. cleanest 3' window (`n_3p`, best of the two orientations);
#' 2. fewest other variants under the primer footprint (`n_var`);
#' 3. highest discrimination (DECIPHER);
#' 4. not a weak G:T / A:C wobble mismatch;
#' 5. quieter 3'-seed genome-wide (`seed_hits`, best of the two strands) — a
#'    specificity tie-break computed only for the top `specificity_candidates`
#'    placements, mirroring [pick_target_snp()];
#' 6. then `|Tm - target|`, 3'-stability and length.
#' The common primer is then placed on the winning strand and ranked by
#' amplicon length, `COM Tm >= AS Tm`, dimers/hairpin, footprint load and Tm
#' proximity (see [design_common_primers()]).
#' @section KASP design pipeline:
#' Typical use is `ath_genome()` -> `ensembl_gene()`/`ensembl_variants()` ->
#' `pick_target_snp()` -> `assemble_kasp_assay()` -> `kasp_design_checks()`.
#' See `vignette("kasp-intro")`.
#'
#' @param genome DNAStringSet (from [ath_genome()]).
#' @param snp    one-row tibble (chrom, pos, ref, alt) from [pick_target_snp()].
#' @param variants variant table (for the same chromosome) to count overlaps.
#' @param strands which AS-primer orientations to search (see [design_as_primers()]).
#' @param tm_range accept AS primers whose Tm falls in this window (C).
#' @param clean_3p size of the "clean" AS 3' window used to count other variants.
#' @param max_footprint_variants hard cap on *other* variants under the AS primer
#'   footprint ("N-masking" as a constraint): `0` allows none, the default `Inf`
#'   only penalises them in the ranking.  Use [mask_positions()] to blank such
#'   positions in your own template if you need the literal `N`s.
#' @param min_disc minimum discrimination (on- minus off-target efficiency).
#' @param specificity_candidates how many leading AS placements to screen for
#'   3'-seed genome-wide uniqueness when breaking ties (mirrors
#'   [pick_target_snp()]; `0` disables the screen).
#' @param run_blast also build/query a genome BLAST db (needs NCBI BLAST+).
#' @param BPPARAM optional `BiocParallel` backend for the genome-wide steps
#'   (seed specificity and in-silico PCR; split by chromosome).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' v <- ensembl_variants("Chr4", 268901, 271503)
#' snp <- pick_target_snp(v, genome, ensembl_gene("AT4G00650"))
#' assemble_kasp_assay(genome, snp)
#' }
#' @section External tool:
#' The thermodynamics engine (DECIPHER via OligoArrayAux) needs `hybrid-min` on
#' the `PATH`; if the data tables are not found automatically, point them out:
#' `with_unafloddat("/path/to/share/oligoarrayaux", assemble_kasp_assay(...))`.
#' @param flank template half-width (bp).
#' @param min_gc,max_gc accepted AS GC window.
#' @family KASP design
#' @export
assemble_kasp_assay <- function(genome, snp, variants = NULL, flank = 150L,
                                strands = c("+", "-"),
                                min_gc = 0.35, max_gc = 0.70,
                                tm_range = c(56, 65), clean_3p = 5L,
                                max_footprint_variants = Inf,
                                min_disc = 0.5, specificity_candidates = 10L,
                                run_blast = FALSE, BPPARAM = NULL) {
  tmpl    <- snp_template(genome, snp$chrom, snp$pos, flank)
  vpos    <- if (is.null(variants)) NULL else variants$pos[variants$chrom == snp$chrom]
  ref_len <- nchar(snp$ref)
  alt_len <- nchar(snp$alt)
  is_indel <- ref_len != alt_len
  off     <- snp$pos - flank           # template coordinate -> genomic coordinate
  vblock  <- seq.int(snp$pos, snp$pos + max(ref_len, alt_len) - 1L)  # the target's own bases

  # all primer geometry comes from the shared helpers (plus-strand coordinates);
  # the REF AS primer's 3'-terminal base is the variant anchor for its strand
  anchor_of <- function(strand) as_anchor(snp$pos, ref_len, strand)
  com_iv    <- function(strand, L, gap) primer_interval(anchor_of(strand), strand, L, gap = gap, role = "COM")

  # candidate AS primers, unified to one schema for either variant kind:
  #   ref_primer/alt_primer, ref_on/ref_off/alt_on/alt_off, discrimination,
  #   tm, gc, as_length, strand, and genomic footprints as_* and as_alt_*
  if (!is_indel) {
    as_cand <- design_as_primers(tmpl$seq, tmpl$snp_i, snp$ref, snp$alt,
                                 strands = strands) |>
      dplyr::filter(dplyr::between(.data$gc, min_gc, max_gc),
                    dplyr::between(.data$tm, tm_range[1], tm_range[2]),
                    !has_run(.data$ref_primer), !has_repeat(.data$ref_primer)) |>
      dplyr::mutate(as_start     = .data$tmpl_start + off,
                    as_end       = .data$tmpl_end   + off,
                    as_alt_start = .data$as_start,
                    as_alt_end   = .data$as_end,
                    wobble       = wobble_mismatch(snp$ref, snp$alt))
  } else {
    ind <- design_as_primers_indel(tmpl$seq, tmpl$snp_i, snp$ref, snp$alt,
                                   strands = strands) |>
      dplyr::filter(dplyr::between(.data$gc, min_gc, max_gc),
                    dplyr::between(.data$tm, tm_range[1], tm_range[2]),
                    !has_run(.data$primer), !has_repeat(.data$primer))
    if (nrow(ind) == 0L) rlang::abort("No indel AS primer passed the GC/Tm filter.")
    # alt-template coordinate -> REF/genomic coordinate (block clamps to REF)
    alt_to_ref <- function(x) dplyr::case_when(
      x < tmpl$snp_i                    ~ x,
      x <= tmpl$snp_i + alt_len - 1L    ~ pmin(x, tmpl$snp_i + ref_len - 1L),
      TRUE                              ~ x - (alt_len - ref_len))
    top <- function(a) ind |>
      dplyr::filter(.data$allele == a) |>
      dplyr::group_by(.data$strand) |>
      dplyr::slice_max(.data$discrimination, n = 1L, with_ties = FALSE) |>
      dplyr::ungroup()
    ref_p <- top("ref"); alt_p <- top("alt")
    as_cand <- dplyr::inner_join(ref_p, alt_p, by = "strand", suffix = c("_ref", "_alt")) |>
      dplyr::transmute(
        strand = .data$strand,
        as_length    = .data$as_length_ref,
        ref_primer   = .data$primer_ref,
        alt_primer   = .data$primer_alt,
        ref_on       = .data$on_eff_ref,  ref_off  = .data$off_eff_ref,
        alt_on       = .data$on_eff_alt,  alt_off  = .data$off_eff_alt,
        discrimination = (.data$discrimination_ref + .data$discrimination_alt) / 2,
        tm = .data$tm_ref, gc = .data$gc_ref,
        dg_3p = terminal_dg(.data$primer_ref),
        wobble = NA,                                   # undefined for length-changing alleles
        as_start     = .data$as_start_ref + off,
        as_end       = .data$as_end_ref   + off,
        as_alt_start = alt_to_ref(.data$as_start_alt) + off,
        as_alt_end   = alt_to_ref(.data$as_end_alt)   + off)
  }

  # clean-3' window and footprint variant load come from the stored interval
  as_cand <- as_cand |>
    dplyr::mutate(
      win3p_start = dplyr::if_else(.data$strand == "+", .data$as_end - clean_3p + 1L, .data$as_start + 1L),
      win3p_end   = dplyr::if_else(.data$strand == "+", .data$as_end - 1L,           .data$as_start + clean_3p - 1L),
      n_3p  = purrr::map2_int(.data$win3p_start, .data$win3p_end, function(a, b)
                count_variants_in(IRanges::IRanges(a, b), vpos, vblock)),
      n_var = purrr::map2_int(.data$as_start, .data$as_end, function(a, b)
                count_variants_in(IRanges::IRanges(a, b), vpos, vblock)))

  # "N-masking" as a hard constraint: optionally drop candidates sitting on
  # other catalogued variants (their footprint would be ambiguous)
  if (is.finite(max_footprint_variants)) {
    as_cand <- dplyr::filter(as_cand, .data$n_var <= max_footprint_variants)
  }
  if (nrow(as_cand) == 0L) {
    rlang::abort(glue::glue(
      "No AS primer passed the GC/Tm filter (max_footprint_variants = {max_footprint_variants})."))
  }

  # prioritise: clean 3' window, then minimal footprint variant load, then
  # discrimination, then a quieter 3'-seed (genome-wide), then the tie-breakers.
  # Only the leading placements are seed-screened (like pick_target_snp()), so the
  # scan stays a tie-break and does not dominate the cost; an indexed genome
  # (no in-memory DNAStringSet) skips it.
  ordered <- as_cand |>
    dplyr::arrange(.data$n_3p, .data$n_var, dplyr::desc(.data$discrimination),
                   .data$wobble, abs(.data$tm - AS_TM_TARGET), .data$dg_3p, .data$as_length) |>
    dplyr::mutate(.as_row = dplyr::row_number())
  ordered$seed_hits <- NA_integer_
  if (specificity_candidates > 0L && !inherits(genome, "kasp_genome")) {
    short <- dplyr::slice_head(ordered, n = specificity_candidates)
    seeds <- stats::setNames(
      substr(short$ref_primer, nchar(short$ref_primer) - SEED_K + 1L, nchar(short$ref_primer)),
      short$.as_row)
    hits <- specificity_seed(seeds, genome, k = SEED_K, max_mismatch = SEED_MISMATCH,
                             BPPARAM = BPPARAM)$hits
    ordered$seed_hits[short$.as_row] <- hits
  }
  best_as <- ordered |>
    dplyr::arrange(.data$n_3p, .data$n_var, dplyr::desc(.data$discrimination),
                   .data$wobble, .data$seed_hits, abs(.data$tm - AS_TM_TARGET),
                   .data$dg_3p, .data$as_length) |>
    dplyr::slice_head(n = 1L)
  as_cand <- dplyr::select(ordered, - ".as_row")

  if (!is.na(best_as$n_3p) && best_as$n_3p > 0L) {
    rlang::warn(glue::glue(
      "Target SNP {snp$chrom}:{snp$pos} has {best_as$n_3p} other variant(s) in the ",
      "AS 3'-terminal {clean_3p} nt (on strand {best_as$strand}) - allele specificity ",
      "may not hold across accessions."))
  }
  if (!is.na(best_as$n_var) && best_as$n_var > 0L) {
    rlang::warn(glue::glue(
      "AS primer overlaps {best_as$n_var} other variant(s) in its footprint ",
      "(no candidate with fewer passed GC/Tm/discrimination filters)."))
  }
  if (best_as$discrimination < min_disc) {
    rlang::warn(glue::glue(
      "Best AS primer discrimination = {round(best_as$discrimination, 2)} ",
      "is below min_disc = {min_disc}; consider widening `strands` or `tm_range`."))
  }

  # COM: design on the SAME strand as the winning AS primer so the amplicon
  # closes.  design_common_primers() works in TEMPLATE coordinates, so its
  # anchor is the template-local variant anchor (not the genomic one).
  common <- design_common_primers(tmpl$seq, as_anchor(tmpl$snp_i, ref_len, best_as$strand),
                                  strand = best_as$strand, tm_range = tm_range) |>
    dplyr::mutate(com_homodimer = purrr::map_dbl(.data$com, ~ oligo_dg(.x, .x)),
                  com_hairpin   = purrr::map_dbl(.data$com, ~ oligo_dg(.x)))

  as_tm_mean <- mean(best_as$tm)

  assays <- common |>
    dplyr::mutate(
      strand = best_as$strand,
      com_start = purrr::pmap_int(list(.data$strand, .data$len, .data$gap),
                     function(st, L, g) IRanges::start(com_iv(st, L, g))),
      com_end   = purrr::pmap_int(list(.data$strand, .data$len, .data$gap),
                     function(st, L, g) IRanges::end(com_iv(st, L, g))),
      n_var  = purrr::map2_int(.data$gap, .data$len, function(g, L)
                count_variants_in(com_iv(best_as$strand, L, g), vpos, vblock)),
      as_ref = best_as$ref_primer, as_alt = best_as$alt_primer,
      discrimination = best_as$discrimination, as_gc = best_as$gc,
      as_len = best_as$as_length, as_n_var = best_as$n_var, as_n_3p = best_as$n_3p,
      # worst-case product over both alleles (indels differ by delta)
      amplicon = pmax(pmax(best_as$as_end, .data$com_end) - pmin(best_as$as_start, .data$com_start) + 1L,
                      pmax(best_as$as_alt_end, .data$com_end) - pmin(best_as$as_alt_start, .data$com_start) + 1L),
      # store the geometry once, for the plots to read instead of recomputing
      as_start    = best_as$as_start,     as_end     = best_as$as_end,
      as_alt_start = best_as$as_alt_start, as_alt_end = best_as$as_alt_end,
      win3p_start = best_as$win3p_start,  win3p_end  = best_as$win3p_end) |>
    # Robustness first: keep the amplicon <= 100 bp and the COM Tm >= AS Tm
    # (EasyKASP), then reject strong self-structure, then minimise footprint
    # variants, then get the COM Tm closest to target.  Logical terms sort
    # FALSE (= criterion met) before TRUE.
    dplyr::arrange(.data$amplicon > KASP_LIMITS$amplicon[2],
                   .data$com_tm < as_tm_mean,
                   .data$com_homodimer < KASP_LIMITS$homodimer_dg,
                   .data$com_hairpin   < KASP_LIMITS$hairpin_dg,
                   .data$n_var, abs(.data$com_tm - COM_TM_TARGET))

  best    <- assays |> dplyr::slice_head(n = 1L)
  primers <- c(AS_ref = best$as_ref, AS_alt = best$as_alt, COM = best$com)

  # specificity uses the template-matched allele for the AS primers (no
  # deliberate mismatch): that answers "is the target site unique?".
  # make_as_primer() works in TEMPLATE coordinates, hence the template anchor.
  tmpl_anchor <- as_anchor(tmpl$snp_i, ref_len, best$strand)
  seed_primers <- c(
    AS_ref = make_as_primer(tmpl$seq, tmpl_anchor, best$as_len, snp$ref, strand = best$strand),
    AS_alt = make_as_primer(tmpl$seq, tmpl_anchor, best$as_len, snp$alt, strand = best$strand),
    COM    = best$com)

  # in-silico PCR: the REF primer must give one product, the ALT primer none
  amplicons <- tryCatch(
    assay_amplicons(list(assay = best, template = tmpl, snp = snp), genome, BPPARAM = BPPARAM),
    error = function(e) { rlang::warn(glue::glue("in-silico PCR failed: {conditionMessage(e)}")); NULL })

  # whole-primer genome scan (only when requested): an attached FM-index is
  # queried first (index-based, multithreaded), otherwise the BLAST database
  blast <- NULL
  if (isTRUE(run_blast)) {
    blast <- if (inherits(genome, "kasp_genome") && !is.null(genome$fm_index)) {
      tryCatch(specificity_fm(primers, genome$fm_index),
               error = function(e) { rlang::warn(glue::glue("FM scan failed: {conditionMessage(e)}")); NULL })
    } else {
      db <- tryCatch(build_blast_db(genome), error = function(e) {
        rlang::warn(glue::glue("BLAST db unavailable: {conditionMessage(e)}")); NULL })
      blast_primers(primers, db)
    }
  }

  list(
    snp = snp, template = tmpl, strand = best$strand, assay = assays,
    as_candidates = as_cand,
    efficiency = best_as |> dplyr::select("ref_on","ref_off","alt_on","alt_off","discrimination"),
    thermo = purrr::imap(primers, ~ dplyr::mutate(primer_thermo(.x), name = .y)) |>
      purrr::list_rbind(),
    complementarity = dimer_thermo(primers),
    specificity = specificity_of(seed_primers, genome, BPPARAM = BPPARAM),
    amplicons = amplicons,
    blast = blast,
    tails = c(AS_ref = paste0(KASP_TAIL_FAM, best$as_ref),
              AS_alt = paste0(KASP_TAIL_HEX, best$as_alt),
              COM    = best$com)
  ) |> structure(class = "kasp_assay")
}

#' Is `x` a KASP assay object?
#'
#' @param x object to test.
#' @family KASP utilities
#' @export
is_kasp_assay <- function(x) inherits(x, "kasp_assay")

#' Do all acceptance checks pass?
#'
#' @param assay assembled assay.
#' @param limits thresholds; see [kasp_design_checks()].
#' @family KASP checks
#' @export
passes_all <- function(assay, limits = KASP_LIMITS) {
  all(kasp_design_checks(assay, limits)$pass)
}

#' Print a KASP assay
#'
#' Compact one-screen summary: target, the three oligos, amplicon size and
#' discrimination.
#'
#' @param x a [assemble_kasp_assay()] result.
#' @param ... unused.
#' @method print kasp_assay
#' @family KASP utilities
#' @export
print.kasp_assay <- function(x, ...) {
  a <- x$assay[1L, ]
  cat("<kasp_assay>",
      glue::glue("{x$snp$chrom}:{x$snp$pos} {x$snp$ref}>{x$snp$alt} (strand {x$strand})"), "\n")
  cat("  AS-ref ", a$as_ref, "\n")
  cat("  AS-alt ", a$as_alt, "\n")
  cat("  COM    ", a$com, "\n")
  cat(" ", glue::glue("amplicon {a$amplicon} bp | discrimination {round(x$efficiency$discrimination[1L], 3)}"), "\n")
  invisible(x)
}


# =============================================================================
#  7b. Batch design over many loci
# =============================================================================

#' Classify a variant by its REF/ALT strings
#' @param ref,alt allele strings.
#' @family KASP design
#' @export
variant_type <- function(ref, alt) {
  if (nchar(ref) == nchar(alt)) ifelse(nchar(ref) == 1L, "SNP", "MNP")
  else ifelse(nchar(alt) > nchar(ref), "insertion", "deletion")
}

#' Design KASP markers for a table of loci, with a per-locus verdict
#'
#' Runs [assemble_kasp_assay()] independently over every row of `loci`, so one
#' hard locus cannot abort the batch: errors are captured as `verdict = "fail"`,
#' warnings or failing acceptance checks give `verdict = "warn"`, and a fully
#' passing assay gives `verdict = "pass"`.  `loci` needs columns
#' `chrom`, `pos`, `ref`, `alt` (as in [read_variants_vcf()]).
#'
#' @param genome DNAStringSet.
#' @param loci tibble of loci (`chrom`, `pos`, `ref`, `alt`).
#' @param variants variant table for footprint accounting; defaults to `loci`.
#' @param BPPARAM optional `BiocParallel` backend; when given, loci are processed
#'   in parallel (a forked backend avoids copying the genome).
#' @param ... passed on to [assemble_kasp_assay()].
#' @return tidy tibble, one row per locus.
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' loci <- tibble::tibble(chrom = "Chr4", pos = 269090L, ref = "G", alt = "A")
#' design_kasp_markers(genome, loci)
#' }
#' @family KASP design
#' @export
design_kasp_markers <- function(genome, loci, variants = loci, BPPARAM = NULL, ...) {
  need <- c("chrom", "pos", "ref", "alt")
  missing <- setdiff(need, names(loci))
  if (length(missing)) {
    rlang::abort(glue::glue("`loci` is missing column(s): {paste(missing, collapse = ', ')}"))
  }
  one <- function(ch, pos, ref, alt) {
    snp  <- tibble::tibble(chrom = as.character(ch), pos = as.integer(pos),
                           ref = as.character(ref), alt = as.character(alt))
    type <- variant_type(snp$ref, snp$alt)
    warns <- character(0)
    res <- withCallingHandlers(
      tryCatch({
        a   <- assemble_kasp_assay(genome, snp, variants = variants, ...)
        chk <- kasp_design_checks(a)
        tibble::tibble(
          chrom = snp$chrom, pos = snp$pos, ref = snp$ref, alt = snp$alt, type = type,
          strand = a$strand, amplicon = a$assay$amplicon[1],
          discrimination = a$efficiency$discrimination[1],
          n_fail = sum(!chk$pass),
          verdict = if (all(chk$pass)) "pass" else "warn",
          note = paste(chk$criterion[!chk$pass], collapse = "; "))
      }, error = function(e) {
        tibble::tibble(chrom = snp$chrom, pos = snp$pos, ref = snp$ref, alt = snp$alt,
                       type = type, strand = NA_character_, amplicon = NA_integer_,
                       discrimination = NA_real_, n_fail = NA_integer_,
                       verdict = "fail", note = conditionMessage(e))
      }),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    # surface any assembly warnings on an otherwise built assay
    if (res$verdict != "fail" && length(warns)) {
      res$verdict <- if (res$verdict == "pass") "warn" else res$verdict
      res$note <- paste(c(res$note, warns), collapse = "; ")
    }
    res
  }

  if (is.null(BPPARAM)) {
    n  <- nrow(loci)
    id <- cli::cli_progress_bar("Designing KASP markers", total = n)
    res <- purrr::imap_dfr(seq_len(n), function(i, ...) {
      cli::cli_progress_update(id = id)
      one(loci$chrom[i], loci$pos[i], loci$ref[i], loci$alt[i])
    })
    cli::cli_progress_done(id = id)
    res
  } else {
    res <- BiocParallel::bplapply(seq_len(nrow(loci)), function(i) {
      worker_biostrings()
      one(loci$chrom[i], loci$pos[i], loci$ref[i], loci$alt[i])
    }, BPPARAM = BPPARAM)
    purrr::list_rbind(res)
  }
}


# =============================================================================
#  8. Reaction recipe (output / ordering)
# =============================================================================

#' Tidy KASP primer-mix + reaction recipe
#'
#' Standard KBioscience/LGC KASP chemistry: the two allele-specific primers are
#' used at half the concentration of the common primer, AS1:AS2:COM = 2:2:5 (the
#' usual working mix is 12/12/30 uM).  This returns the recipe as tidy data:
#'
#' * `stage = "primer mix"` - how to build the working primer mix from the
#'   orderable 100 uM oligo stocks;
#' * `stage = "reaction"` - the per-reaction setup (2x master mix, assay mix,
#'   DNA, water), scaled to `n_reactions` with a small `overfill` for pipetting.
#'
#' It is pure arithmetic, so it carries no chemistry assumptions beyond the
#' published LGC recipe; change `as_uM`/`com_uM` for other chemistries.
#'
#' @param n_reactions number of wells to prepare for.
#' @param mix_uL total volume of working primer mix to prepare.
#' @param stock_uM oligo stock concentration as ordered (IDT default 100 uM).
#' @param as_uM,com_uM working concentrations of each AS and the common primer.
#' @param reaction_uL final reaction volume per well.
#' @param master_mix_conc fold-concentrate of the supplied master mix.
#' @param dna_uL template DNA volume per reaction.
#' @param overfill fraction of extra volume to compensate pipetting loss.
#' @return tibble `stage`, `component`, `stock_uM`, `final_uM`, `per_rxn_uL`,
#'   `total_uL`.
#' @examples
#' kasp_mix_recipe(n_reactions = 96)
#' @family KASP output
#' @export
kasp_mix_recipe <- function(n_reactions = 96L, mix_uL = 100,
                            stock_uM = 100, as_uM = 12, com_uM = 30,
                            reaction_uL = 10, master_mix_conc = 2,
                            dna_uL = 1, overfill = 0.10) {
  ratios      <- c(AS1_FAM = as_uM, AS2_HEX = as_uM, COM = com_uM)
  stock_share <- ratios / stock_uM                       # fraction of the mix
  water_share <- 1 - sum(stock_share)
  if (water_share < 0) {
    rlang::abort(glue::glue(
      "Requested working mix ({as_uM}/{as_uM}/{com_uM} uM) needs more than ",
      "100% of {stock_uM} uM stock; raise `stock_uM` or lower the targets."))
  }

  primer_mix <- tibble::tibble(
    component = c(names(ratios), "nuclease-free water"),
    stock_uM  = c(rep(stock_uM, length(ratios)), NA_real_),
    final_uM  = c(unname(ratios), NA_real_),
    per_rxn_uL = NA_real_,
    total_uL  = c(unname(stock_share) * mix_uL, water_share * mix_uL),
    stage     = "primer mix")

  assay_per_rxn <- round(0.014 * reaction_uL, 3)
  master_per_rxn <- reaction_uL / master_mix_conc
  water_per_rxn  <- reaction_uL - master_per_rxn - assay_per_rxn - dna_uL
  if (water_per_rxn < 0) {
    rlang::abort(glue::glue(
      "Reaction volumes exceed {reaction_uL} uL; check `dna_uL`/`reaction_uL`."))
  }
  total_rxn <- ceiling(n_reactions * (1 + overfill))

  reaction <- tibble::tibble(
    component  = c(glue::glue("{master_mix_conc}x KASP master mix"), "primer mix",
                   "template DNA", "nuclease-free water"),
    stock_uM   = NA_real_,
    final_uM   = NA_real_,
    per_rxn_uL = c(master_per_rxn, assay_per_rxn, dna_uL, water_per_rxn),
    total_uL   = c(master_per_rxn, assay_per_rxn, dna_uL, water_per_rxn) * total_rxn,
    stage      = "reaction")

  dplyr::bind_rows(primer_mix, reaction)
}


# =============================================================================
#  9. Acceptance checks - "is this a good KASP?"
# =============================================================================

#' Longest homopolymer run in a sequence
#' @keywords internal
#' @noRd
longest_run <- function(seq) {
  max(nchar(stringr::str_extract_all(seq, "([ACGT])\\1*")[[1L]]))
}

#' Check an assembled assay against the KASP acceptance criteria
#'
#' Turns the prose design manual into a tidy pass/fail table so "good" and
#' "bad" assays are graded the same way every time (thresholds: `KASP_LIMITS`).
#' Reads only the assay object produced by [assemble_kasp_assay()], so it is
#' pure and cheap to test.
#'
#' @param assay result of [assemble_kasp_assay()].
#' @param limits list of thresholds; defaults to `KASP_LIMITS`.
#' @return tibble with `criterion`, `observed`, `limit` and `pass`.
#' @examples
#' \dontrun{
#' assay <- assemble_kasp_assay(genome, snp)
#' kasp_design_checks(assay)
#' }
#' @family KASP checks
#' @export
kasp_design_checks <- function(assay, limits = KASP_LIMITS) {
  aa    <- assay$assay[1L, ]
  th    <- assay$thermo
  d     <- assay$complementarity
  disc  <- assay$efficiency$discrimination[1L]

  as_tm <- mean(th$tm[th$name %in% c("AS_ref", "AS_alt")])
  com_tm <- aa$com_tm
  hair  <- min(th$hairpin_dg)
  cross <- min(d$dg[d$a != d$b])          # worst hetero-dimer
  homo  <- min(d$dg[d$a == d$b])          # worst self-dimer
  run   <- max(longest_run(aa$as_ref), longest_run(aa$com))

  # in-silico PCR (only if computed): REF primer -> 1 product, ALT primer -> 0
  amp        <- assay$amplicons
  amp_obs    <- if (is.null(amp)) NA_real_ else max(amp$n_products)
  amp_pass   <- if (is.null(amp)) FALSE else {
    ref <- amp$n_products[amp$primer == "AS_ref"]; alt <- amp$n_products[amp$primer == "AS_alt"]
    length(ref) == 1L && length(alt) == 1L && ref == 1L && alt == 0L
  }

  tibble::tibble(
    criterion = c("amplicon length (bp)", "AS Tm (C)", "COM Tm (C)",
                  "COM Tm >= AS Tm", "AS GC", "COM GC",
                  "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)",
                  "inverted repeat (arms>=4, loop<=6)", "worst hairpin dG37",
                  "worst cross-dimer dG37", "worst homodimer dG37",
                  "discrimination", "in-silico PCR (REF 1 / ALT 0)"),
    observed  = c(aa$amplicon, as_tm, com_tm, com_tm - as_tm, aa$as_gc, aa$com_gc,
                  run, sum(lengths(find_repeats(c(aa$as_ref, aa$com)))),
                  sum(purrr::map_int(find_palindromes(c(aa$as_ref, aa$as_alt, aa$com)), nrow)),
                  hair, cross, homo, disc, amp_obs),
    limit     = c(glue::glue("[{limits$amplicon[1]}, {limits$amplicon[2]}]"),
                  glue::glue("[{limits$as_tm[1]}, {limits$as_tm[2]}]"),
                  glue::glue("[{limits$com_tm[1]}, {limits$com_tm[2]}]"),
                  ">= 0",
                  glue::glue("[{limits$gc[1]}, {limits$gc[2]}]"),
                  glue::glue("[{limits$gc[1]}, {limits$gc[2]}]"),
                  glue::glue("< {limits$max_run}"),
                  "0",
                  "0",
                  glue::glue(">= {limits$hairpin_dg}"),
                  glue::glue(">= {limits$duplex_dg}"),
                  glue::glue(">= {limits$homodimer_dg}"),
                  glue::glue(">= {limits$discrimination}"),
                  "REF=1, ALT=0"),
    pass      = c(dplyr::between(aa$amplicon, limits$amplicon[1], limits$amplicon[2]),
                  dplyr::between(as_tm, limits$as_tm[1], limits$as_tm[2]),
                  dplyr::between(com_tm, limits$com_tm[1], limits$com_tm[2]),
                  com_tm >= as_tm,
                  dplyr::between(aa$as_gc, limits$gc[1], limits$gc[2]),
                  dplyr::between(aa$com_gc, limits$gc[1], limits$gc[2]),
                  run < limits$max_run,
                  !has_repeat(aa$as_ref) && !has_repeat(aa$com),
                  !any(has_palindrome(c(aa$as_ref, aa$as_alt, aa$com))),
                  hair >= limits$hairpin_dg,
                  cross >= limits$duplex_dg,
                  homo >= limits$homodimer_dg,
                  disc >= limits$discrimination,
                  amp_pass)
  )
}
