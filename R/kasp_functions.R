# =============================================================================
#  kasp_functions.R  --  design + QC primitives for KASP assays
#  Target organism: Arabidopsis thaliana (TAIR10)
#  100% R.  No Python.  Heavy thermodynamics is delegated to packages.
# =============================================================================
#
#  WHAT IS KASP (in one paragraph)
#  ---------------------------------------------------------------------------
#  Kompetitive Allele-Specific PCR (KASP, LGC Biosearch Technologies) is a
#  homogeneous, endpoint-fluorescence assay for scoring SNPs and small indels.
#  One assay = THREE primers:
#     * AS1 - allele-specific forward primer carrying a FAM  universal 5' tail
#     * AS2 - allele-specific forward primer carrying a HEX/VIC universal tail
#     * COM - one common (reverse) primer shared by both alleles
#  The two AS primers are identical except for their 3'-terminal base, which is
#  the SNP itself.  A FRET cassette + Taq polymerase read out which tail got
#  incorporated, so discrimination comes from *competitive annealing* of the two
#  AS primers + the inability of Taq to extend a 3'-mismatched primer.
#
#  DESIGN MANUAL (condensed, as of 2026)
#  ---------------------------------------------------------------------------
#  Follows LGC/KBioscience KASP guidance together with the peer-reviewed
#  EasyKASP study (Zhang, Yang & Wen, BMC Bioinformatics 2025, 26:292,
#  doi:10.1186/s12859-025-06322-x), which re-derived the optimal ranges from
#  500 validated KASP assays.  Defaults, not dogma - validate in your lab.
#  The machine-checkable form of this checklist is kasp_design_checks(), with
#  thresholds in KASP_LIMITS, so prose and code cannot drift apart.
#
#  Amplicon
#    * keep it short: 50-100 bp (EasyKASP: amplification products < 100 bp).
#
#  Common primer (COM)
#    * length 18-25 nt, GC 40-60%, no runs of >= 4 identical bases; low 3'-end
#      stability; no self/cross 3' dimers.
#    * its Tm should be HIGHER than the two allele-specific primers (EasyKASP:
#      the common primer typically has lower GC but more length and higher Tm).
#
#  Allele-specific primers (AS1/AS2)
#    * the SNP is the 3'-TERMINAL base - the whole point of the assay.
#    * target-specific part 20-28 nt, within 20-30 nt of the variant; Tm
#      (matching allele) ~ 55-65 C, below the common primer's Tm.
#    * prefer the strand whose +/-25 nt around the variant has GC ~ 45% and no
#      repeats (EasyKASP step 1) - the both-strand search below does this.
#    * add a DELIBERATE secondary mismatch 1-3 bases from the 3' end
#      (usually -2 or -3): it destabilises the non-target allele more than the
#      target one and is the main lever for discrimination.  We search for it
#      by maximising the on-target amplification efficiency (below).
#    * avoid the last ~5 bases covering another common SNP/indel, and 4+ base
#      homopolymer runs.
#
#  Discrimination metric
#    * we score the *Taq*-relevant quantity directly: the predicted PCR
#      efficiency (hybridisation + 3' elongation) of each AS primer on the
#      target vs the non-target allele (DECIPHER::CalculateEfficiencyPCR, which
#      models side reactions and Taq extension bias at 3' mismatches).
#      A good assay has high on-target efficiency and near-zero off-target.
#
#  Secondary structure / cross-complementarity
#    * reject hairpin-prone oligos (palindromes) and, above all, 3'-end dimers
#      between any primer pair (AS1-COM, AS2-COM, AS1-AS2, COM self).
#
#  Specificity
#    * the complete assay must amplify ONE product genome-wide (primers must be
#      single-copy; EasyKASP).  An off-target with an intact 3'-terminal match
#      is a killer - drop it.
#    * non-target variants in the flanking window matter; EasyKASP recommends
#      Sanger-checking the variant and its flanks for critical markers.
#
#  Validation
#    * always genotype three controls - homozygous REF, homozygous ALT and a
#      heterozygote - plus a no-template control.
#
#  ACCEPTANCE CHECKLIST (what a "good" assay must satisfy; KASP_LIMITS/checks)
#    * amplicon 50-100 bp; |AS Tm| in 55-65 C; COM Tm in 57-62 C and >= AS Tm
#    * GC 0.40-0.60; no homopolymer run >= 4 nt; hairpin dG37 >= -3 kcal/mol
#    * no cross-dimer or homodimer below -9 kcal/mol (37 C)
#    * Taq discrimination >= 0.5; the 3'-seed is unique genome-wide
#
#  ---------------------------------------------------------------------------
#  WHICH PACKAGE DOES WHAT (all R)
#    - DECIPHER       : PCR efficiency incl. Taq 3'-mismatch bias (the KASP
#                       discrimination metric).  Needs OligoArrayAux (`hybrid-min`).
#    - TmCalculator   : nearest-neighbour Tm with salt/DMSO correction (SantaLucia 2004 et al.)
#    - Biostrings     : sequence ops, GC (letterFrequency), local alignments
#                       (dimers), vcountPattern (specificity)
#    - rBLAST         : optional genome-wide BLAST (needs NCBI BLAST+ in PATH)
#    - httr2/jsonlite : Ensembl Plants REST (variation, gene lookup)
#  We hand-write only the KASP-specific bit (3'-anchored AS primers and the
#  candidate scan for the common primer) because no maintained R package does
#  KASP design natively (primeR was retired; openPrimeR needs MAFFT/ViennaRNA/MELTING).
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
SO_TERM_RANK  <- c(missense = 1L, splice = 2L, synonymous = 3L, UTR = 4L)
SO_RANK_OTHER <- 5L

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

# OligoArrayAux thermodynamics (the DNA engine DECIPHER already uses):
#   hybrid-min    -> duplex (dimer) free energy of two strands
#   hybrid-ss-min -> intramolecular folding (hairpin) free energy of one strand
# `na` is [Na+] in M; returns dG in kcal/mol (more negative = more stable).
oligo_dg <- function(seq1, seq2 = NULL, na = 0.05, temp = 37) {
  bin <- if (is.null(seq2)) "hybrid-ss-min" else "hybrid-min"
  exe <- Sys.which(bin)
  if (!nzchar(exe)) return(NA_real_)
  args <- c("-n", "DNA", "-t", as.character(temp), "-T", as.character(temp),
            "-N", as.character(na), "-E", "-q")
  args <- if (is.null(seq2)) c(args, seq1) else c(args, seq1, seq2)
  out <- suppressWarnings(system2(exe, args, stdout = TRUE, stderr = FALSE))
  out <- out[nzchar(trimws(out))]
  if (!length(out)) return(NA_real_)
  dg <- suppressWarnings(as.numeric(utils::tail(out, 1L)))
  if (length(dg) != 1L || !is.finite(dg)) {
    rlang::warn(glue::glue("OligoArrayAux ('{bin}') returned no parsable dG."))
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
get_region <- function(genome, chrom, start, end) {
  Biostrings::subseq(genome[[chrom]], start = start, end = end)
}


# --- primer geometry (single source of truth) --------------------------------
# All intervals are plus-strand, with `snp_pos` the SNP coordinate.  The SAME
# helpers are reused by design (in template coordinates: pass `snp_i`), by the
# assay assembly and by the plots, so the geometry cannot drift between them.

#' Plus-strand interval occupied by a primer
#'
#' `role = "AS"` = allele-specific primer (3'-terminal base at the SNP);
#' `role = "COM"` = the common primer, `gap` bases away on the other side.
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
#' @export
count_variants_in <- function(interval, vpos, exclude = NA_integer_) {
  if (is.null(vpos)) return(NA_integer_)
  v <- vpos[!is.na(vpos) & !(vpos %in% exclude)]
  if (!length(v)) return(0L)
  as.integer(length(IRanges::findOverlaps(IRanges::IRanges(v, v), interval)))
}

#' Download (once) and load the TAIR10 genome as a DNAStringSet
#'
#' The FASTA comes from Ensembl Plants and is cached under `data/cache/`.
#' Sequence names are normalised to Chr1..Chr5 / ChrM / ChrC (scaffolds
#' dropped), so they line up with the Ensembl gene records used elsewhere.
#' @export
ath_genome <- function(cache = here::here("data", "cache"), release = 63L) {
  fs::dir_create(cache)
  fa <- fs::path(cache, "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
  if (!fs::file_exists(fa)) {
    url <- sprintf(paste0(
      "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-%d/fasta/",
      "arabidopsis_thaliana/dna/",
      "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"), release)
    message("Downloading TAIR10 genome (~35 MB) from Ensembl Plants ...")
    utils::download.file(url, fa, mode = "wb", quiet = TRUE)
  }
  g <- Biostrings::readDNAStringSet(fa)
  nm <- sub("\\s.*$", "", names(g))   # keep only the first token of the defline
  dropped <- nm[!nm %in% names(ATH_CHROM_MAP)]
  if (length(dropped)) {
    rlang::inform(glue::glue(
      "ath_genome(): dropping {length(dropped)} non-chromosome sequence(s) ",
      "(e.g. {paste(utils::head(dropped, 3L), collapse = ', ')})."))
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
#' @export
ensembl_gene <- function(gene_id, species = "arabidopsis_thaliana",
                         cache = here::here("data", "cache")) {
  fs::dir_create(cache)
  path <- fs::path(cache, glue::glue("{species}_{gene_id}.json"))
  if (!fs::file_exists(path)) {
    resp <- httr2::request("https://rest.ensembl.org") |>
      httr2::req_url_path_append("lookup", "id", gene_id) |>
      httr2::req_url_query(expand = 1) |>
      httr2::req_headers(`Content-Type` = "application/json") |>
      httr2::req_user_agent("kasp-example (lecture demo)") |>
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
#' @export
ensembl_variants <- function(chrom, start, end,
                             species = "arabidopsis_thaliana",
                             cache = here::here("data", "cache")) {
  fs::dir_create(cache)
  region <- glue::glue("{sub('^Chr', '', chrom)}:{start}-{end}")
  path <- fs::path(cache, glue::glue("{species}_var_{region}.json"))
  if (!fs::file_exists(path)) {
    resp <- httr2::request("https://rest.ensembl.org") |>
      httr2::req_url_path_append("overlap", "region", species, region) |>
      httr2::req_url_query(feature = "variation") |>
      httr2::req_headers(`Content-Type` = "application/json") |>
      httr2::req_user_agent("kasp-example (lecture demo)") |>
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
#' @export
strand_gc_hint <- function(genome, chrom, pos, flank = 25L) {
  up <- as.character(get_region(genome, chrom, max(1L, pos - flank), pos - 1L))
  dn <- as.character(get_region(genome, chrom, pos + 1L, pos + flank))
  tibble::tibble(side = c("upstream", "downstream"), gc = c(gc_frac(up), gc_frac(dn))) |>
    dplyr::mutate(distance_to_45 = abs(.data$gc - 0.45))
}

#' Map Ensembl consequence strings to a design priority (lower = preferred)
#'
#' Uses the `SO_TERM_RANK` constant so the vocabulary lives in one place.
#' @keywords internal
#' @noRd
so_rank <- function(consequence) {
  vapply(consequence, function(x) {
    hit <- names(SO_TERM_RANK)[stringr::str_detect(x, names(SO_TERM_RANK))]
    if (length(hit)) unname(SO_TERM_RANK[[hit[1L]]]) else SO_RANK_OTHER
  }, integer(1), USE.NAMES = FALSE)
}

#' Choose a designable target variant from a variant table
#'
#' Keeps variants we can actually design on: a biallelic substitution block whose
#' REF and ALT are ACGT-only strings of *equal* length (SNP or MNP).  Indels and
#' multi-allelic records are left for a dedicated engine and are not silently
#' forced into a substitution.
#'
#' @param clean_3p number of flanking bases of the AS primer that should be free
#'   of *other* variants.  Because the AS primer may be placed on either strand,
#'   the "clean" count is the best (smallest) of the two orientations: the chosen
#'   variant just has to be clean in *at least one* direction.
#' @export
pick_target_snp <- function(variants, genome, gene, variant_id = NULL, clean_3p = 5L) {
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
                  stringr::str_detect(.data$alt, "^[ACGT]+$"))
  if (nrow(v) == 0L) rlang::abort("No usable substitution variant (SNP/MNP) left after filtering.")

  # every variant of the gene (incl. indels): used to penalise targets whose AS
  # 3' neighbourhood is itself variable (breaks allele specificity across lines)
  allpos <- variants$pos[variants$chrom == gene$chrom]

  v |>
    dplyr::mutate(
      transition = .data$ref_len == 1L &
        paste0(.data$ref, .data$alt) %in% c("AG", "GA", "CT", "TC"),
      rank = so_rank(.data$consequence),
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
    dplyr::arrange(.data$n_3p, .data$rank, dplyr::desc(.data$transition), .data$pos) |>
    dplyr::slice(1L) |>
    dplyr::select(-"rank", -"ref_genome")
}


# =============================================================================
#  2. Thermodynamics via packages
# =============================================================================

#' GC fraction and homopolymer-run test
#' @keywords internal
#' @noRd
gc_frac  <- function(seq) unname(Biostrings::letterFrequency(Biostrings::DNAStringSet(seq), "GC", as.prob = TRUE)[, 1])
has_run  <- function(seq, n = 4L) stringr::str_detect(seq, glue::glue("([ACGT])\\1{{{n - 1L},}}"))

#' Tandem repeats (simple sequence repeats) in a sequence
#'
#' EasyKASP screens candidate primers for SSR repeat units, not just
#' homopolymers.  Returns, for each input, the repeated units of length `units`
#' occurring at least `min_copies` times in tandem.
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
#' @export
has_repeat <- function(seq, units = 2:4, min_copies = 3L)
  lengths(find_repeats(seq, units, min_copies)) > 0L

#' Nearest-neighbour Tm of oligos (perfect duplex), salt-corrected (TmCalculator)
#'
#' `res$gr` uses TmCalculator's `$` accessor for the result GRanges; `Tm` is a
#' metadata column on it.
#' @export
primer_tm <- function(seq) {
  gr  <- TmCalculator::to_genomic_ranges(as.character(seq))
  res <- TmCalculator::tm_calculate(gr, method = "tm_nn", nn_table = NN_TABLE,
                                    Na = PCR_COND$Na, Mg = PCR_COND$Mg,
                                    dNTPs = PCR_COND$dNTPs)
  as.numeric(S4Vectors::mcols(res$gr)$Tm)
}

#' OligoArrayAux availability (needed by DECIPHER for thermodynamics)
#' @export
has_oligoarrayaux <- function() nzchar(Sys.which("hybrid-min"))

#' PCR amplification efficiency of primer/target pairs (DECIPHER)
#'
#' Models hybridisation (incl. side reactions) times Taq elongation efficiency,
#' so a 3'-terminal mismatch is strongly penalised - the physical basis of KASP.
#'
#' @param primers character vector, 5'->3'.
#' @param targets character vector (same length) with the reverse-complement of
#'   the binding site (as in `?DECIPHER::CalculateEfficiencyPCR`).
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
                                  if (is.na(.data$stab_base)) NULL else .data$stab_base)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(gc = gc_frac(.data$ref_primer), tm = primer_tm(.data$ref_primer)) |>
    dplyr::select("strand","as_length","stab_at","stab_base","ref_primer","alt_primer",
                  "ref_on","ref_off","alt_on","alt_off","discrimination","tm","gc") |>
    dplyr::arrange(dplyr::desc(.data$discrimination))
}


# =============================================================================
#  5. Structural QC : Tm, hairpins, self- and cross-complementarity
# =============================================================================

#' Full structural profile of one primer
#' @export
primer_thermo <- function(primer) {
  s <- as.character(primer)
  tibble::tibble(
    primer       = s,
    len          = nchar(s),
    tm           = primer_tm(s),
    gc           = gc_frac(s),
    hairpin_dg   = oligo_dg(s))        # hybrid-ss-min, kcal/mol
}

#' Duplex (dimer) free energy for every primer pair (OligoArrayAux, DNA)
#'
#' `dg` is the hybridisation free energy of the two strands (kcal/mol, 37 C);
#' more negative = stronger dimer.  The diagonal is the homodimer of each primer.
#' Returns the FULL matrix (both orders); symmetry is imposed here, once.
#' @export
dimer_thermo <- function(primers) {
  nm <- names(primers); n <- length(nm)
  m <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
  for (i in seq_len(n)) {
    m[i, i] <- oligo_dg(primers[[i]], primers[[i]])          # homodimer
    for (j in seq_len(i - 1L)) {                             # strictly lower triangle
      v <- oligo_dg(primers[[i]], primers[[j]])
      m[i, j] <- v; m[j, i] <- v
    }
  }
  d <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  tibble::tibble(a = as.character(d$Var1), b = as.character(d$Var2),
                 dg = as.numeric(d$Freq))
}


# =============================================================================
#  6. Genome-wide specificity (pure R, plus optional BLAST)
# =============================================================================

#' Count genome-wide near-matches of each primer's 3'-terminal k-mer
#'
#' Uses `Biostrings::vcountPattern` (vectorised over the genome sequences),
#' summing over both strands for each primer.
#' @export
specificity_seed <- function(primers, genome, k = SEED_K, max_mismatch = SEED_MISMATCH) {
  seeds <- unname(substr(primers, nchar(primers) - k + 1L, nchar(primers)))
  count_one <- function(s) {
    s  <- Biostrings::DNAString(s)
    rc <- Biostrings::reverseComplement(s)
    sum(Biostrings::vcountPattern(s,  genome, max.mismatch = max_mismatch)) +
      sum(Biostrings::vcountPattern(rc, genome, max.mismatch = max_mismatch))
  }
  tibble::tibble(primer = names(primers), seed = seeds,
                 hits = as.integer(purrr::map_int(seeds, count_one)))
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
#' @return tibble `chrom`, `start`, `end`, `width` (one row per product).
#' @export
in_silico_pcr <- function(F_primer, R_primer, genome, max_mismatch = 0L) {
  Fp <- Biostrings::DNAString(F_primer)
  Rp <- Biostrings::DNAString(R_primer)
  purrr::imap(genome, function(chr, nm) {
    pr <- Biostrings::matchProbePair(Fp, Rp, chr, max.mismatch = max_mismatch, verbose = FALSE)
    if (is.null(pr) || length(pr) == 0L) return(NULL)
    tibble::tibble(chrom = nm, start = IRanges::start(pr), end = IRanges::end(pr),
                   width = IRanges::width(pr))
  }) |>
    purrr::list_rbind()
}

#' Products predicted for each allele-specific primer against the common one
#'
#' On the reference genome the REF primer must give exactly one product; the
#' ALT primer (the other base / the indel) should give none.  Uses the
#' template-matched primers (no deliberate stabiliser), like [specificity_seed()].
#' @export
assay_amplicons <- function(assay, genome, max_mismatch = 0L) {
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
    prods <- in_silico_pcr(fr[["F"]], fr[["R"]], genome, max_mismatch = max_mismatch)
    tibble::tibble(primer = nm, n_products = nrow(prods),
                   width = if (nrow(prods)) max(prods$width) else NA_integer_)
  })
}

#' BLAST every primer against a local genome database (rBLAST, optional)
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
#' @export
build_blast_db <- function(genome, path = here::here("data", "cache", "tair10.fa")) {
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
#' @param genome DNAStringSet (from [ath_genome()]).
#' @param snp    one-row tibble (chrom, pos, ref, alt) from [pick_target_snp()].
#' @param variants variant table (for the same chromosome) to count overlaps.
#' @param strands which AS-primer orientations to search (see [design_as_primers()]).
#' @param tm_range accept AS primers whose Tm falls in this window (C).
#' @param clean_3p size of the "clean" AS 3' window used to count other variants.
#' @param min_disc minimum discrimination (on- minus off-target efficiency).
#' @param run_blast also build/query a genome BLAST db (needs NCBI BLAST+).
#' @export
assemble_kasp_assay <- function(genome, snp, variants = NULL, flank = 150L,
                                strands = c("+", "-"),
                                min_gc = 0.35, max_gc = 0.70,
                                tm_range = c(56, 65), clean_3p = 5L,
                                min_disc = 0.5, run_blast = FALSE) {
  tmpl    <- snp_template(genome, snp$chrom, snp$pos, flank)
  vpos    <- if (is.null(variants)) NULL else variants$pos[variants$chrom == snp$chrom]
  ref_len <- nchar(snp$ref)
  vblock  <- seq.int(snp$pos, snp$pos + ref_len - 1L)   # the target's own bases

  # all primer geometry comes from the shared helpers (plus-strand coordinates);
  # the AS primer's 3'-terminal base is the variant anchor for its strand
  anchor_of <- function(strand) as_anchor(snp$pos, ref_len, strand)
  as_iv    <- function(strand, L) primer_interval(anchor_of(strand), strand, L, role = "AS")
  com_iv   <- function(strand, L, gap) primer_interval(anchor_of(strand), strand, L, gap = gap, role = "COM")
  win3p_iv <- function(strand) as_3p_interval(anchor_of(strand), strand, clean_3p)

  as_cand <- design_as_primers(tmpl$seq, tmpl$snp_i, snp$ref, snp$alt,
                               strands = strands) |>
    dplyr::filter(dplyr::between(.data$gc, min_gc, max_gc),
                  dplyr::between(.data$tm, tm_range[1], tm_range[2]),
                  !has_run(.data$ref_primer), !has_repeat(.data$ref_primer)) |>
    dplyr::mutate(
      n_3p  = purrr::map_int(.data$strand, function(st)
                count_variants_in(win3p_iv(st), vpos, vblock)),
      n_var = purrr::map2_int(.data$strand, .data$as_length, function(st, L)
                count_variants_in(as_iv(st, L), vpos, vblock))
    )
  if (nrow(as_cand) == 0L) rlang::abort("No AS primer passed the GC/Tm filter.")

  # prioritise: clean 3' window, then minimal footprint variant load, then
  # discrimination, then Tm closest to the AS target
  best_as <- as_cand |>
    dplyr::arrange(.data$n_3p, .data$n_var, dplyr::desc(.data$discrimination),
                   abs(.data$tm - AS_TM_TARGET), .data$as_length) |>
    dplyr::slice_head(n = 1L)

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
      n_var  = purrr::map2_int(.data$gap, .data$len, function(g, L)
                count_variants_in(com_iv(best_as$strand, L, g), vpos, vblock)),
      as_ref = best_as$ref_primer, as_alt = best_as$alt_primer,
      discrimination = best_as$discrimination, as_gc = best_as$gc,
      as_len = best_as$as_length, as_n_var = best_as$n_var, as_n_3p = best_as$n_3p,
      amplicon = .data$gap + .data$len + .data$as_len - 1L,
      # store the geometry once, for the plots to read instead of recomputing
      as_start   = purrr::map2_int(.data$strand, .data$as_len,
                     function(st, L) IRanges::start(as_iv(st, L))),
      as_end     = purrr::map2_int(.data$strand, .data$as_len,
                     function(st, L) IRanges::end(as_iv(st, L))),
      com_start  = purrr::pmap_int(list(.data$strand, .data$len, .data$gap),
                     function(st, L, g) IRanges::start(com_iv(st, L, g))),
      com_end    = purrr::pmap_int(list(.data$strand, .data$len, .data$gap),
                     function(st, L, g) IRanges::end(com_iv(st, L, g))),
      win3p_start = purrr::map_int(.data$strand,
                      function(st) IRanges::start(win3p_iv(st))),
      win3p_end   = purrr::map_int(.data$strand,
                      function(st) IRanges::end(win3p_iv(st)))) |>
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
    assay_amplicons(list(assay = best, template = tmpl, snp = snp), genome),
    error = function(e) { rlang::warn(glue::glue("in-silico PCR failed: {conditionMessage(e)}")); NULL })

  blast <- NULL
  if (isTRUE(run_blast)) {
    db <- tryCatch(build_blast_db(genome), error = function(e) {
      rlang::warn(glue::glue("BLAST db unavailable: {conditionMessage(e)}")); NULL })
    blast <- blast_primers(primers, db)
  }

  list(
    snp = snp, template = tmpl, strand = best$strand, assay = assays,
    as_candidates = as_cand,
    efficiency = best_as |> dplyr::select("ref_on","ref_off","alt_on","alt_off","discrimination"),
    thermo = purrr::imap(primers, ~ dplyr::mutate(primer_thermo(.x), name = .y)) |>
      purrr::list_rbind(),
    complementarity = dimer_thermo(primers),
    specificity = specificity_seed(seed_primers, genome),
    amplicons = amplicons,
    blast = blast,
    tails = c(AS_ref = paste0(KASP_TAIL_FAM, best$as_ref),
              AS_alt = paste0(KASP_TAIL_HEX, best$as_alt),
              COM    = best$com)
  )
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
#' @param limits list of thresholds; defaults to [KASP_LIMITS].
#' @return tibble with `criterion`, `observed`, `limit` and `pass`.
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
                  "longest homopolymer (bp)", "SSR (2-4 nt, >=3 copies)", "worst hairpin dG37",
                  "worst cross-dimer dG37", "worst homodimer dG37",
                  "discrimination", "in-silico PCR (REF 1 / ALT 0)"),
    observed  = c(aa$amplicon, as_tm, com_tm, com_tm - as_tm, aa$as_gc, aa$com_gc,
                  run, sum(lengths(find_repeats(c(aa$as_ref, aa$com)))),
                  hair, cross, homo, disc, amp_obs),
    limit     = c(glue::glue("[{limits$amplicon[1]}, {limits$amplicon[2]}]"),
                  glue::glue("[{limits$as_tm[1]}, {limits$as_tm[2]}]"),
                  glue::glue("[{limits$com_tm[1]}, {limits$com_tm[2]}]"),
                  ">= 0",
                  glue::glue("[{limits$gc[1]}, {limits$gc[2]}]"),
                  glue::glue("[{limits$gc[1]}, {limits$gc[2]}]"),
                  glue::glue("< {limits$max_run}"),
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
                  hair >= limits$hairpin_dg,
                  cross >= limits$duplex_dg,
                  homo >= limits$homodimer_dg,
                  disc >= limits$discrimination,
                  amp_pass)
  )
}
