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
#  Values below follow LGC/KBioscience guidance and the standard primer-design
#  literature.  They are defaults, not dogma - validate in your lab.
#
#  Amplicon
#    * keep it short: 60-120 bp is the sweet spot (hard limits ~50-150 bp).
#
#  Common primer (COM)
#    * length 18-25 nt, Tm 57-62 C, GC 40-60%; low 3'-end stability;
#      no runs of >= 4 identical bases; no self/cross 3' dimers.
#
#  Allele-specific primers (AS1/AS2)
#    * the SNP is the 3'-TERMINAL base - the whole point of the assay.
#    * target-specific part 20-28 nt; Tm (matching allele) ~ 60-65 C.
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
#    * the complete assay must amplify ONE product genome-wide.  An off-target
#      with an intact 3'-terminal match is a killer - drop it.
#
#  Validation
#    * always genotype three controls - homozygous REF, homozygous ALT and a
#      heterozygote - plus a no-template control.
#
#  ---------------------------------------------------------------------------
#  WHICH PACKAGE DOES WHAT (all R)
#    - DECIPHER       : PCR efficiency incl. Taq 3'-mismatch bias (the KASP
#                       discrimination metric).  Needs OligoArrayAux (`hybrid-min`).
#    - TmCalculator   : nearest-neighbour Tm with salt/DMSO correction (SantaLucia 2004 et al.)
#    - Biostrings     : sequence ops, GC, palindromes (hairpins), local
#                       alignments (dimers), vmatchPattern (specificity), PDict
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
  suppressWarnings(as.numeric(utils::tail(out, 1L)))
}


# --- small utilities ---------------------------------------------------------

#' Normalise chromosome names to "Chr1".."Chr5"/"ChrM"/"ChrC"
#' @keywords internal
#' @noRd
ath_chrom <- function(x) {
  x <- sub("^chromosome", "", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  paste0("Chr", x)
}

#' Read one genomic region as a DNAString
#' @keywords internal
#' @noRd
get_region <- function(genome, chrom, start, end) {
  Biostrings::subseq(genome[[chrom]], start = start, end = end)
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
  nm[nm == "Mt"] <- "ChrM"; nm[nm == "Pt"] <- "ChrC"
  nm[nm %in% as.character(1:5)] <- paste0("Chr", nm[nm %in% as.character(1:5)])
  names(g) <- nm
  g[intersect(c(paste0("Chr", 1:5), "ChrM", "ChrC"), names(g))]
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
  jsonlite::fromJSON(readr::read_file(path)) |>
    tibble::as_tibble() |>
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

#' Choose a designable target SNP from a variant table
#' @export
#' @param clean_3p number of 3'-terminal bases of the AS primer that should be
#'   free of *other* variants.  Because the AS primer may be placed on either
#'   strand, the "clean" count is the best (smallest) of the two orientations:
#'   the chosen SNP just has to be clean in *at least one* direction.
pick_target_snp <- function(variants, genome, gene, variant_id = NULL, clean_3p = 5L) {
  v <- variants |>
    dplyr::filter(
      purrr::map_lgl(.data$ref_genome, ~ length(.x) == 2L),
      purrr::map_lgl(.data$ref_genome, ~ all(nchar(.x) == 1L)),
      purrr::map_lgl(.data$ref_genome, ~ !any(.x %in% c("-", "")))
    )
  if (!is.null(variant_id)) v <- dplyr::filter(v, .data$id == variant_id)
  if (nrow(v) == 0L) rlang::abort("No usable biallelic SNP left after filtering.")

  # every variant of the gene (incl. indels): used to penalise SNPs whose AS
  # 3' neighbourhood is itself variable (breaks allele specificity across lines)
  allpos <- variants$pos[variants$chrom == gene$chrom]

  v |>
    dplyr::mutate(
      ref = purrr::map2_chr(.data$chrom, .data$pos,
                            ~ as.character(get_region(genome, .x, .y, .y))),
      alt = purrr::map2_chr(.data$ref_genome, .data$ref,
                            ~ setdiff(.x, .y)[1] %||% NA_character_),
      transition = paste0(.data$ref, .data$alt) %in% c("AG", "GA", "CT", "TC"),
      rank = dplyr::case_when(
        stringr::str_detect(.data$consequence, "missense")   ~ 1L,
        stringr::str_detect(.data$consequence, "splice")     ~ 2L,
        stringr::str_detect(.data$consequence, "synonymous") ~ 3L,
        stringr::str_detect(.data$consequence, "UTR")        ~ 4L,
        TRUE                                                 ~ 5L
      ),
      # best-case (either strand) count of other variants in the AS 3' window:
      #   "+" orientation -> [p-clean_3p+1, p-1];  "-" orientation -> [p+1, p+clean_3p-1]
      n_3p = pmin(
        vapply(.data$pos, function(p)
          sum(allpos >= p - clean_3p + 1L & allpos <= p - 1L), integer(1)),
        vapply(.data$pos, function(p)
          sum(allpos >= p + 1L & allpos <= p + clean_3p - 1L), integer(1)))
    ) |>
    dplyr::filter(!is.na(.data$alt)) |>
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
gc_frac  <- function(seq) Biostrings::letterFrequency(Biostrings::DNAStringSet(seq), "GC", as.prob = TRUE)[, 1]
has_run  <- function(seq, n = 4L) stringr::str_detect(seq, sprintf("([ACGT])\\1{%d,}", n - 1L))

#' Nearest-neighbour Tm of oligos (perfect duplex), salt-corrected (TmCalculator)
#' @export
primer_tm <- function(seq) {
  gr  <- TmCalculator::to_genomic_ranges(as.character(seq))
  res <- TmCalculator::tm_calculate(gr, method = "tm_nn", nn_table = NN_TABLE,
                                    Na = PCR_COND$Na, Mg = PCR_COND$Mg,
                                    dNTPs = PCR_COND$dNTPs)
  as.numeric(res$gr$Tm)
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
                                  gaps = 18:80, target_tm = 60,
                                  min_gc = 0.40, max_gc = 0.60,
                                  tm_range = c(57, 64), n = 5L) {
  wlen <- length(window)
  grid <- expand.grid(L = lengths, gap = gaps, KEEP.OUT.ATTRS = FALSE)
  plus <- identical(strand, "+")

  cand <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    L <- grid$L[i]; gap <- grid$gap[i]
    # plus-strand template interval of the common primer
    if (plus) { start <- snp_i + gap;             end <- start + L - 1L }
    else      { end   <- snp_i - gap;             start <- end - L + 1L }
    if (start < 1L || end > wlen) return(NULL)
    seq_plus <- Biostrings::subseq(window, start, end)
    com <- if (plus) as.character(Biostrings::reverseComplement(seq_plus))
           else      as.character(seq_plus)
    tibble::tibble(com = com, len = L, gap = gap, start = start, end = end)
  }) |>
    dplyr::mutate(com_gc = gc_frac(.data$com),
                  com_tm = primer_tm(.data$com),
                  run    = has_run(.data$com)) |>
    dplyr::filter(!.data$run,
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

#' Plus-strand binding window of an AS primer of a given strand/length
#' @keywords internal
#' @noRd
as_window <- function(window, snp_i, length, strand) {
  if (identical(strand, "+")) Biostrings::subseq(window, snp_i - length + 1L, snp_i)
  else                        Biostrings::subseq(window, snp_i, snp_i + length - 1L)
}

#' DECIPHER target (= strand bound by the primer) for a binding window
#' @keywords internal
#' @noRd
as_target <- function(seg, strand) {
  if (identical(strand, "+")) Biostrings::reverseComplement(seg) else seg
}

#' Build one allele-specific primer on a given strand
#'
#' The 3'-terminal base is the SNP allele.  If `stab_at` is given, the base
#' `stab_at` positions from the 3' end is replaced by `stab_base` (the deliberate
#' stabilising mismatch).
#'
#' `strand = "+"`: SNP is the 3'-terminal base, the primer extends upstream.
#' `strand = "-"`: primer is the reverse complement (extends downstream).
#' @keywords internal
#' @noRd
make_as_primer <- function(window, snp_i, length, allele, strand = "+",
                           stab_at = NULL, stab_base = NULL) {
  seg <- as.character(as_window(window, snp_i, length, strand))
  p   <- strsplit(seg, "")[[1]]
  p[if (identical(strand, "+")) length(p) else 1L] <- allele
  s <- paste0(p, collapse = "")
  if (!identical(strand, "+"))
    s <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
  if (!is.null(stab_at)) {
    q <- strsplit(s, "")[[1]]
    q[length(q) - stab_at + 1L] <- stab_base
    s <- paste0(q, collapse = "")
  }
  s
}

#' Grid-search allele-specific primers maximising Taq-based discrimination
#'
#' For every strand, length, secondary-mismatch position and mismatch base we
#' build the REF and ALT 3'-anchored primers and predict the PCR efficiency of
#' each on both alleles (DECIPHER).  Discrimination = on-target minus off-target
#' efficiency, averaged over the two primers.
#'
#' @param strands which strands to consider for the AS primer; `"+"` = forward
#'   primer with the SNP as its 3'-terminal base extending downstream, `"-"` =
#'   the mirror-image primer extending upstream.  Both are valid KASP layouts.
#' @export
design_as_primers <- function(window, snp_i, ref, alt, strands = c("+", "-"),
                              lengths = 20:27, stab_ats = c(NA, 2L, 3L),
                              temp = PCR_EFF$temp, P = PCR_EFF$P, ions = PCR_EFF$ions) {
  # candidate mismatch bases at `k` bases from the 3' end: any base other than
  # the one the (perfectly matching) primer already carries there.
  mismatch_bases <- function(L, k, strand) {
    base_p <- make_as_primer(window, snp_i, L, ref, strand = strand)
    cur    <- substr(base_p, L - k + 1L, L - k + 1L)
    setdiff(c("A", "C", "G", "T"), cur)
  }

  # enumerate configurations (strand x length x mismatch position x mismatch base)
  cfgs <- purrr::map_dfr(strands, function(st) {
    purrr::map_dfr(lengths, function(L) {
      purrr::map_dfr(stab_ats, function(k) {
        bases <- if (is.na(k)) NA_character_ else mismatch_bases(L, k, st)
        tibble::tibble(as_length = L, stab_at = k, stab_base = bases, strand = st)
      })
    })
  })

  # build all primer/target pairs and score them in ONE batched DECIPHER call
  pairs <- purrr::pmap_dfr(cfgs, function(as_length, stab_at, stab_base, strand) {
    L <- as_length
    k <- if (is.na(stab_at))   NULL else stab_at
    b <- if (is.na(stab_base)) NULL else stab_base
    # binding window on the plus strand, and the same window carrying the ALT allele
    seg_ref <- as.character(as_window(window, snp_i, L, strand))
    seg_alt <- seg_ref
    pos_in_seg <- if (identical(strand, "+")) L else 1L
    substr(seg_alt, pos_in_seg, pos_in_seg) <- alt
    ref_p <- make_as_primer(window, snp_i, L, ref, strand = strand, stab_at = k, stab_base = b)
    alt_p <- make_as_primer(window, snp_i, L, alt, strand = strand, stab_at = k, stab_base = b)
    t_ref <- as.character(as_target(Biostrings::DNAString(seg_ref), strand))
    t_alt <- as.character(as_target(Biostrings::DNAString(seg_alt), strand))
    tibble::tibble(
      strand = strand, as_length = L, stab_at = stab_at, stab_base = stab_base,
      primer = c(ref_p, ref_p, alt_p, alt_p), which = c("ref_on","ref_off","alt_on","alt_off"),
      target = c(t_ref, t_alt, t_alt, t_ref))
  })
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
      ref_primer = make_as_primer(window, snp_i, .data$as_length, ref, strand = .data$strand,
                                  if (is.na(.data$stab_at)) NULL else .data$stab_at,
                                  if (is.na(.data$stab_base)) NULL else .data$stab_base),
      alt_primer = make_as_primer(window, snp_i, .data$as_length, alt, strand = .data$strand,
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
    primer     = s,
    len        = nchar(s),
    tm         = primer_tm(s),
    gc         = gc_frac(s),
    hairpin_dg = oligo_dg(s))          # hybrid-ss-min, kcal/mol
}

#' Duplex (dimer) free energy for every primer pair (OligoArrayAux, DNA)
#'
#' `dg` is the hybridisation free energy of the two strands (kcal/mol, 37 C);
#' more negative = stronger dimer.  The diagonal is the homodimer of each primer.
#' @export
dimer_thermo <- function(primers) {
  nm <- names(primers)
  purrr::map_dfr(seq_along(nm), function(i) {
    purrr::map_dfr(i:length(nm), function(j) {
      tibble::tibble(a = nm[i], b = nm[j],
                     dg = oligo_dg(primers[[i]], primers[[j]]))
    })
  })
}


# =============================================================================
#  6. Genome-wide specificity (pure R, plus optional BLAST)
# =============================================================================

#' Count genome-wide near-matches of each primer's 3'-terminal k-mer
#' @export
specificity_seed <- function(primers, genome, k = 15L, max_mismatch = 1L) {
  count_strand <- function(pat) {
    m <- Biostrings::vmatchPattern(Biostrings::DNAString(pat), genome,
                                   max.mismatch = max_mismatch, fixed = "subject")
    sum(S4Vectors::elementNROWS(m))
  }
  tibble::tibble(primer = names(primers),
                 seed = substr(primers, nchar(primers) - k + 1L, nchar(primers))) |>
    dplyr::mutate(hits = purrr::map_int(.data$seed, function(s)
      count_strand(s) + count_strand(as.character(
        Biostrings::reverseComplement(Biostrings::DNAString(s))))))
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
  tmpl   <- snp_template(genome, snp$chrom, snp$pos, flank)
  vpos   <- if (is.null(variants)) NULL else variants$pos[variants$chrom == snp$chrom]

  # number of *other* variants under a plus-strand genomic interval
  n_var_between <- function(from, to) {
    if (is.null(vpos)) return(NA_integer_)
    sum(vpos >= from & vpos <= to & vpos != snp$pos)
  }
  # plus-strand interval covered by an AS primer of a given strand/length
  as_footprint <- function(strand, L) {
    if (identical(strand, "+")) c(snp$pos - L + 1L, snp$pos - 1L)
    else                        c(snp$pos + 1L,     snp$pos + L - 1L)
  }
  # plus-strand interval of the AS 3'-terminal `clean_3p` bases
  as_3p_window <- function(strand) {
    if (identical(strand, "+")) c(snp$pos - clean_3p + 1L, snp$pos - 1L)
    else                        c(snp$pos + 1L,           snp$pos + clean_3p - 1L)
  }

  as_cand <- design_as_primers(tmpl$seq, tmpl$snp_i, snp$ref, snp$alt,
                               strands = strands) |>
    dplyr::filter(dplyr::between(.data$gc, min_gc, max_gc),
                  dplyr::between(.data$tm, tm_range[1], tm_range[2])) |>
    dplyr::mutate(
      n_3p  = purrr::map_int(.data$strand, function(st) {
        w <- as_3p_window(st); n_var_between(w[1], w[2]) }),
      n_var = purrr::map2_int(.data$strand, .data$as_length, function(st, L) {
        w <- as_footprint(st, L); n_var_between(w[1], w[2]) })
    )
  if (nrow(as_cand) == 0L) rlang::abort("No AS primer passed the GC/Tm filter.")

  # prioritise: clean 3' window, then minimal footprint variant load, then
  # discrimination, then Tm closest to 62 C
  best_as <- as_cand |>
    dplyr::arrange(.data$n_3p, .data$n_var, dplyr::desc(.data$discrimination),
                   abs(.data$tm - 62), .data$as_length) |>
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

  # COM: design on the SAME strand as the winning AS primer so the amplicon closes
  common <- design_common_primers(tmpl$seq, tmpl$snp_i, strand = best_as$strand,
                                  tm_range = tm_range)
  com_footprint <- function(gap, L) {
    if (identical(best_as$strand, "+")) c(snp$pos + gap,             snp$pos + gap + L - 1L)
    else                                c(snp$pos - gap - L + 1L,    snp$pos - gap)
  }

  assays <- common |>
    dplyr::mutate(
      strand = best_as$strand,
      n_var  = purrr::map2_int(.data$gap, .data$len, function(g, L) {
        w <- com_footprint(g, L); n_var_between(w[1], w[2]) }),
      as_ref = best_as$ref_primer, as_alt = best_as$alt_primer,
      discrimination = best_as$discrimination, as_gc = best_as$gc,
      as_len = best_as$as_length, as_n_var = best_as$n_var, as_n_3p = best_as$n_3p,
      amplicon = .data$gap + .data$len + .data$as_len - 1L) |>
    dplyr::arrange(.data$n_var, abs(.data$com_tm - 60))

  best    <- assays |> dplyr::slice_head(n = 1L)
  primers <- c(AS_ref = best$as_ref, AS_alt = best$as_alt, COM = best$com)

  # specificity uses the template-matched allele for the AS primers (no
  # deliberate mismatch): that answers "is the target site unique?".
  seed_primers <- c(
    AS_ref = make_as_primer(tmpl$seq, tmpl$snp_i, best$as_len, snp$ref, strand = best$strand),
    AS_alt = make_as_primer(tmpl$seq, tmpl$snp_i, best$as_len, snp$alt, strand = best$strand),
    COM    = best$com)

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
    thermo = purrr::imap_dfr(primers, ~ dplyr::mutate(primer_thermo(.x), name = .y)),
    complementarity = dimer_thermo(primers),
    specificity = specificity_seed(seed_primers, genome),
    blast = blast,
    tails = c(AS_ref = paste0(KASP_TAIL_FAM, best$as_ref),
              AS_alt = paste0(KASP_TAIL_HEX, best$as_alt),
              COM    = best$com)
  )
}
