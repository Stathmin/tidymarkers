# =============================================================================
#  Shared helpers for the scientific tests (ZW0 / ZW1 / ZW3 / Z7 / Z8).
# =============================================================================
#  Two kinds of fixture:
#
#   * a synthetic genome, generated deterministically, for design-level
#     invariants (never needs the network or a downloaded genome);
#   * a SYNTHETIC MARKER ASSAY built from REAL published primers, for the
#     check-level tests.  The primers are Makhoul 2020's Table S1; the surrounding
#     fields (Tm, GC, dropout) are computed by this package, so the fixture tests
#     the CHECKS, not the published assay outcome.
#
#  One exception: the Makhoul supplement fixture is the real workbook.

## --- verdict access that also understands a doomed design -------------------
# `marker_doomed` carries `reason`/`stage`, not `$evaluation`, so passing it to
# `evaluation_verdict()` returns "unknown" -- an artifact of the accessor, not a
# verdict.  Comparing through this helper is what keeps M2 from failing on that.
design_outcome <- function(m) {
  if (is.null(m)) {
    return("none")
  }
  if (inherits(m, "marker_doomed")) {
    return("doomed")
  }
  evaluation_verdict(m)
}

## --- sequence helpers -------------------------------------------------------

#' Reverse complement of a DNA string
rc <- function(x) {
  x <- as.character(x)
  chartr("ACGTacgt", "TGCAtgca", paste(rev(strsplit(x, "")[[1L]]), collapse = ""))
}

#' Write a single-line-per-sequence FASTA that `indexFa` accepts
#'
#' Rsamtools' faidx rejects a file whose lines differ in width, and a single
#' unwrapped line has only one width, so this is the safest shape for fixtures.
#' (Wrapping at 60 requires the sequence length to be a multiple of 60.)
write_fasta <- function(path, seqs, names = NULL) {
  names <- names %||% "Chr1"
  seqs <- as.character(seqs)
  names <- rep_len(as.character(names), length(seqs))
  writeLines(c(rbind(paste0(">", names), seqs)), path)
  invisible(path)
}

#' A reproducible synthetic genome with one chromosome
#'
#' @param seed RNG seed; the same seed always gives the same sequence.
#' @param n length in bp (exactly).
#' @return list with `path` (FASTA) and `seq`.
syn_genome <- function(seed = 1L, n = 400L) {
  set.seed(seed)
  # keep homopolymer runs short and repeats rare, so most loci are designable
  pool <- c("A", "C", "G", "T", "AC", "CG", "GT", "TA", "GCA", "TGC")
  s <- paste0(sample(pool, ceiling(n / 2) + 1L, replace = TRUE), collapse = "")
  if (nchar(s) < n) s <- paste0(s, strrep("ACG", n))
  s <- substr(s, 1L, n)
  # plain tempfile: the result outlives this call and its caller owns the file
  path <- tempfile("syntmp", fileext = ".fa")
  write_fasta(path, s)
  list(path = path, seq = s)
}

#' One-locus SNP tibble at `pos` with a transition (or caller-chosen ALT)
locus_at <- function(seq, pos, alt = NULL, chrom = "Chr1") {
  ref <- substr(seq, pos, pos)
  if (is.null(alt)) {
    alt <- if (ref %in% c("A", "G")) "C" else "T"
    if (alt == ref) alt <- if (ref == "C") "A" else "G"
  }
  tibble::tibble(
    chrom = chrom, pos = as.integer(pos),
    ref = ref, alt = alt, alt2 = NA_character_
  )
}

#' A genome with a BLAST index attached (needed by the specificity stage)
indexed_syn_genome <- function(gen = syn_genome()) {
  g0 <- marker_genome(gen$path)
  db <- marker_blast_index(g0)
  marker_genome(gen$path, blast_db = db)
}

## --- a marker_assay built from real published primers ------------------------

#' The Makhoul 2020 KASP primer table (13 real assay designs)
makhoul_primers <- function() {
  readr::read_csv(
    test_path("../../inst/extdata/makhoul_kasp_primers.csv"),
    show_col_types = FALSE
  )
}

#' A `marker_assay` whose primers are REAL, so `marker_design_checks()` can run
#'
#' Every field the KASP evaluator reads is derived here from the primer strings
#' (Tm, GC, GC, hairpin/dimer stand-ins) or set to a neutral value.  Nothing in
#' this object claims a wet-lab result.
make_assay <- function(as_ref, as_alt, com, ref = "T", alt = "C", amplicon = 49L,
                       pos = 100L, hairpin = -1, cross = -2, disc = 0.9) {
  tm <- primer_tm(c(as_ref, as_alt, com))
  structure(
    list(
      scheme = "snp",
      snp = tibble::tibble(
        chrom = "5B", pos = as.integer(pos),
        ref = ref, alt = alt, alt2 = NA_character_
      ),
      assay = tibble::tibble(
        start = 1L, end = as.integer(amplicon), len = as.integer(amplicon),
        as_ref = as_ref, as_alt = as_alt, as_alt_degen = FALSE,
        com = com, amplicon = as.integer(amplicon),
        as_gc = gc_frac(as_ref), com_gc = gc_frac(com),
        as_tm = tm[1L], com_tm = tm[3L],
        as_len = nchar(as_ref),
        n_var = 0L, w_var = 0L, as_n_var = 0L, as_n_3p = 0L,
        as_dropout = 0, as_foot_dropout = 0, com_dropout = 0,
        com_hairpin = hairpin, com_homodimer = hairpin,
        run = 2L, ssr = 0L, rank = 1L, strand = "+",
        as_start = 1L, as_end = nchar(as_ref),
        as_alt_start = 1L, as_alt_end = nchar(as_alt),
        com_start = 30L, com_end = 30L + nchar(com),
        win3p_start = 1L, win3p_end = 5L
      ),
      thermo = list(
        tm = tm, name = c("AS_ref", "AS_alt", "COM"),
        hairpin_dg = c(hairpin, hairpin, hairpin)
      ),
      complementarity = list(
        a = c("f", "f", "f", "h", "h", "c"),
        b = c("f", "h", "c", "h", "c", "c"),
        dg = c(-1, cross, -3, -1, cross, -1)
      ),
      efficiency = list(discrimination = disc),
      site = tibble::tibble(
        n_alt = 1L, alt_covered = 1L, disc = TRUE,
        class = "biallelic"
      )
    ),
    class = "marker_assay"
  )
}

#' A synthetic assay for Makhoul row `i`
makhoul_assay <- function(i = 1L) {
  p <- makhoul_primers()
  make_assay(p$fam[[i]], p$hex[[i]], p$common[[i]],
    ref = substr(gsub("[^ACGT]", "", p$snp_kasp[[i]]), 1L, 1L),
    alt = substr(gsub("[^ACGT]", "", p$snp_kasp[[i]]), 2L, 2L)
  )
}
