# =============================================================================
#  data-raw/dpcr_tm_jrc2025.R  --  Z0.2: a dPCR assay's own printed oligo Tm
# =============================================================================
#  JRC H5Nx ddPCR assays (Euro Surveill 2025;30(33):2500183,
#  10.2807/1560-7917.es.2025.30.33.2500183, PMC12372895).  Table 1 prints, for
#  two hydrolysis-probe assays, the sequence, length, GC % and Tm of every
#  primer and probe, plus the amplicon.  This is the first dPCR assay in this
#  package that prints an oligo Tm at all.
#
#  WHAT THIS SET IS, AND WHAT IT IS NOT
#
#    The Tm values are COMPUTED, not measured: the Methods say the oligos were
#    "designed ... [with] Primer3Plus" and Table 1 reports the design tool's
#    output.  No UV melting curve was run.  The buffer is NOT printed either:
#    the reactions use a vendor kit ("one-step ddRT-PCR Advanced Kit for
#    Probes", Bio-Rad; QIAcuity 4x One-Step Advanced Probe Master Mix), so
#    [Na+] and [Mg2+] are unknown numbers.  So this set CANNOT be used to
#    calibrate primer_tm() against a measured anchor -- it is landed as a
#    CONTRAST set, tier "computed, not measured", for Z0.2 only.
#
#    What it does show is that the textbook primer rule is not what a real
#    dPCR assay obeys: neither probe (60.9 C and 57.8 C) is within +8-10 C of
#    its primers (59.2/60.2 and 61.0/61.9 C), and the two land on OPPOSITE
#    sides -- JRC-HA's probe is 0.7 C above its hotter primer, JRC-MP's is
#    4.1 C below -- and the assays work (annealing optimised to 60 C by
#    gradient, 50-62 C).  It says nothing about the -16.86 C level offset.
#
#  ARITHMETIC ASSERTED BELOW is the paper's own, recomputed from its printed
#  sequences only: length, GC %, and amplicon size / primer ends.  The Tm
#  column is transcribed verbatim and NOT asserted against any model.
#
#  Run:  Rscript data-raw/dpcr_tm_jrc2025.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

doi <- "10.2807/1560-7917.es.2025.30.33.2500183"
pmcid <- "PMC12372895"
source_file <- "Euro Surveill 2025;30(33):2500183, Table 1"

# printed Table 1, verbatim (probe sequences stripped of fluorophore/quencher)
t1 <- data.frame(
  assay_id = c("jrc_ha", "jrc_ha", "jrc_ha", "jrc_mp", "jrc_mp", "jrc_mp"),
  primer_or_probe = c("primer_fwd", "primer_rev", "probe_fam", "primer_fwd", "primer_rev", "probe_vic"),
  sequence = c(
    "ACTGGGCTCAGAAATAGTCCTC",
    "TCCCCTGCTCATTGCTATGATG",
    "AGAGGGAGGATGGCAGGGAATGGTT",
    "TGTCTTTGCAGGGAAGAACACC",
    "ACGGTGAGCGTGAACACAAATC",
    "ATCTTGAGGCTCTCATGGAATGGCT"
  ),
  length_printed = c(22L, 22L, 25L, 22L, 22L, 25L),
  gc_printed = c(50.0, 50.0, 56.0, 50.0, 50.0, 48.0),
  tm_printed_C = c(59.2, 60.2, 60.9, 61.0, 61.9, 57.8),
  stringsAsFactors = FALSE
)
# printed amplicons (Table 1) and the two printed amplicon sizes
amplicon <- c(jrc_ha = "ACTGGGCTCAGAAATAGTCCTCTAAGAGAAAAGAGAAGAAAAAGAGGTCTGTTTGGGGCGATAGCAGGGTTTATAGAGGGAGGATGGCAGGGAATGGTTGATGGTTGGTATGGGTACCATCATAGCAATGAGCAGGGGA",
              jrc_mp = "TGTCTTTGCAGGGAAGAACACCGATCTTGAGGCTCTCATGGAATGGCTAAAGACAAGACCAATCCTGTCACCTCTGACTAAGGGGATTTTGGGATTTGTGTTCACGCTCACCGT")
size_printed <- c(jrc_ha = 139L, jrc_mp = 114L)

# --- the paper's own arithmetic, recomputed from its printed sequences ------
stopifnot(
  nrow(t1) == 6L,
  !anyDuplicated(paste(t1$assay_id, t1$primer_or_probe)),
  nchar(t1$sequence) == t1$length_printed,
  vapply(t1$sequence, function(s) all(strsplit(s, "")[[1L]] %in% c("A", "C", "G", "T")), logical(1)),
  # GC % as printed: round(100 * #GC / length, 1)
  round(
    vapply(t1$sequence, function(s) 100 * sum(strsplit(s, "")[[1L]] %in% c("G", "C")) / nchar(s), 0),
    1L
  ) == t1$gc_printed,
  # amplicon size as printed, and the Fwd/Rev primers really are its two ends
  nchar(amplicon) == size_printed,
  startsWith(amplicon[[1L]], t1$sequence[1L]), startsWith(amplicon[[2L]], t1$sequence[4L]),
  endsWith(amplicon[[1L]], paste(chartr("ACGT", "TGCA", rev(strsplit(t1$sequence[2L], "")[[1L]])), collapse = "")),
  endsWith(amplicon[[2L]], paste(chartr("ACGT", "TGCA", rev(strsplit(t1$sequence[5L], "")[[1L]])), collapse = "")),
  # each probe lies inside its amplicon
  grepl(t1$sequence[3L], amplicon[["jrc_ha"]], fixed = TRUE),
  grepl(t1$sequence[6L], amplicon[["jrc_mp"]], fixed = TRUE)
)

oligo_definition <- paste(
  "Tm of a Primer3Plus-designed oligo, ordered by length, GC and the Tm column",
  "of Table 1; oligos ran at 900 nM (primers) / 400 nM (probes) in a Bio-Rad",
  "one-step ddRT-PCR Advanced Kit on a QX600, annealing optimised 50-62 C, 60 C used;",
  "kit-dependent [Na+]/[Mg2+] not reported, so the Tm basis is unknown"
)
amplicon_definition <- paste(
  "Amplicon span of the assay, printed in Table 1 as a sequence and a size and",
  "re-checked here against the two primers; no Tm is measured on it, so the row",
  "carries the span only, to let a reader line a Tm up with the sequence it was",
  "designed against"
)

out <- data.frame(
  assay_id = t1$assay_id,
  primer_or_probe = t1$primer_or_probe,
  sequence = t1$sequence,
  tm_measured_C = t1$tm_printed_C,
  tm_method = "Primer3Plus (design tool output, printed in Table 1; not a UV melt)",
  salt_mM = NA_integer_,
  mg_mM = NA_integer_,
  amplicon_size_bp = size_printed[t1$assay_id],
  outcome_definition = oligo_definition,
  evidence_tier = "computed, not measured",
  source_file = source_file,
  doi = doi,
  pmcid = pmcid,
  stringsAsFactors = FALSE
)

# the amplicon itself, one row per assay: `sequence` holds the full span, the
# size column is repeated so every row of an assay carries its amplicon size
amplicon_rows <- data.frame(
  assay_id = names(amplicon),
  primer_or_probe = "amplicon",
  sequence = unname(amplicon),
  tm_measured_C = NA_real_,
  tm_method = NA_character_,
  salt_mM = NA_integer_,
  mg_mM = NA_integer_,
  amplicon_size_bp = size_printed,
  outcome_definition = amplicon_definition,
  evidence_tier = "computed, not measured",
  source_file = source_file,
  doi = doi,
  pmcid = pmcid,
  stringsAsFactors = FALSE
)
out <- rbind(out, amplicon_rows)
# `tm_method` is NA on the amplicon rows; the column is character, not factor
stopifnot(nrow(out) == 8L, !anyNA(out$sequence), !anyNA(out$amplicon_size_bp))

utils::write.csv(out, "inst/extdata/dpcr_tm_jrc2025.csv", row.names = FALSE)
message("wrote inst/extdata/dpcr_tm_jrc2025.csv: ", nrow(out), " rows")
