# =============================================================================
#  data-raw/readout_hrm.R  --  HRM melt: conditions -> Tm (diatom, PLoS ONE 2025)
# =============================================================================
#  Source (DOI is PROVENANCE; every number below is re-read from the fetched
#  full text, nothing is typed from memory):
#
#    Zhou Y, Ha S, Xu Y, Qin X, Ma Y, Lu J, Wang B, Cai J, Duan Z, Cong B, ...
#    Establishment of a simple prediction method for DNA melting temperature:
#    high-resolution melting curve analysis based on sequence characteristics.
#    PLoS ONE 2025.  DOI 10.1371/journal.pone.0321885   EPMC PMC12002481
#    Fetched: .tmp/a6/PMC12002481.xml (Europe PMC fullTextXML) -> plain text.
#
#  WHY THIS SET EXISTS
#  ROADMAP §V may only switch an HRM model on when a "conditions -> outcome"
#  table exists that the model must reproduce.  Here the SAME HRM run melts
#  20 amplicons (5 diatom species x 4 primer pairs) and the paper prints, for
#  every amplicon: the Sanger product sequence (Table 7), its GC% and length
#  (Table 8), its summed dH/dS (Table 8), and the measured Tm (Table 6).  So
#  one table constrains a nearest-neighbour Tm prediction AND is reproducible
#  from the printed sequences - which is exactly what a `simulated = TRUE` HRM
#  model must do.
#
#  WHAT IT DOES *NOT* CONSTRAIN
#   * Salt: the paper prints no Mg2+ and no dNTP/dsDNA concentration.  The
#     primers go into a proprietary "MonAmp 2x Taq mix Pro" and a 20x EvaGreen
#     dye; the reaction is 20 uL.  Mg is therefore NA in every row and this set
#     CANNOT fit a salt term - it constrains the sequence term only.
#   * dTm: these are ABSOLUTE Tm of different amplicons, not the dTm of an
#     HRM genotyping call.  The existing dTm sets (jsbbs2018, ijms2022,
#     gundry2008) are untouched and remain the discrimination anchor.
#   * The paper's Table 4 (psaA-2 Tm vs DNA input and cycle number, and its
#     supplementary S1 Table) is NOT pooled here: Table 4 prints psaA-2
#     Navicula at 82.38 in one place and 83.40 in another, and the two tables
#     disagree on psaA-2 Navicula/Skeletonema by ~1 C.  The inconsistency is
#     reported in data-raw/proposed/v-readout-conditions.md, not reconciled.
#
#  ASSERTIONS (the source's own arithmetic, checked rather than trusted)
#   * Table 8 GC% recomputed from the Table 7 sequence equals the printed GC%
#     to 0.01 for 22 of the 25 sequenced amplicons; the 3 that disagree, the
#     +1 nt between Table 7 length and Table 8 bp and the Dia18S rows with no
#     sequence are kept as columns rather than corrected.
#   * Every printed GC% lies in the paper's own stated 34.62-55.24% window.  Its
#     stated bp window ("60 to 152") is contradicted by its own Table 8 (59 and
#     154), so the bp check asserts the printed extremes instead.
#   * Table 9's "Difference" column is |measured - predicted| for every row, and
#     the printed "average difference 0.67 C" is the mean of |measured -
#     predicted| over the same eight cells.
#
#  Run:  Rscript data-raw/readout_hrm.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

doi <- "10.1371/journal.pone.0321885"
pmcid <- "PMC12002481"
source_file <- "fulltext.xml :: Table 6 (Tm), Table 7 (sequences), Table 8 (GC/bp/dH/dS), Table 9 (predicted vs measured)"

# --- Table 7: Sanger product sequences, one per (primer, species) ------------
seq_text <- c(
  UPA99 = paste0(
    "CCTCTTATGCTTTACACTCTAGATACGATTTCTACCCGTTCAGAGGAAACTTTGTACGTTTCC ",
    "CCTTTATGCCTTTATACTCTAGATACGATTTCTATCCGTTCAGAGGAAACCTTTGTACGTTTCC ",
    "TTATGCCTTTCACTCTAGATACGATTTCTACCCGTTCAGAGGATACCTTTGTACGTTTCC ",
    "CCTTTATGCCTTTACACTCTAGATACGATTTCTACCCGTTCAGAGGAAACTTTGTACGTTTCC ",
    "CTCTTATGCCTTTACACTCTAGATACGATTTCTACACGTTCAGAGGAAACCTTTGTACGTTTCC"
  ),
  `18S-V7` = paste0(
    "TGTTTTCCTTGTGAAGGCAGGCACTTCTTAGAGGGACAATTGGTGTTTAGCCAATGGAAGTGGGAGGCAATAACAGGTCTGTGATGCC ",
    "CCGCGGATGATCTTTCATTGGCGAGAGCTTCTTAGAGGGACGACCGTCTACTAGACGGAGGAAGATAGGGGCAATAACAGGTCTGTGATGCC ",
    "CGCGATGGCTTTTCATTGGCGAGGTCTTCTTAGAGGGACGTTCATTCTACAAGATGAAGGAAGATGGCGGCAATAACAGGTCTGTGATGCC ",
    "CGATGACTTTTCATTGGCGATGGTTTCTTAGAGGGACATGTAGTACAAAGCTACAGGAAGATTGGGGCAATAACAGGTCTGTGATGCC ",
    "GAGTGAGTTTTCACTGGGTGAAGCTTCTTAGAGGGACGTGCATTCTATTAGATGCAGGAAGATAGGGGCAATAACAGGTCTGTGATGCC"
  ),
  ND = paste0(
    "ATCGAGGTTTAAGGTGGTTTAGACTTCTTAAAAGATGATGAGAACATTAACTCTCAACCATTCATGCGTTGGAGAGAA ",
    "AATGAGGTTTAAGGTGGTTTAGACTTCTTAAAAGATGATGAAAACATTAACTCTCAACCATTCATGCGTTGGAGAGAA ",
    "ATGACGTTTAAAGGTGGTTTAGACTTCTTAAAGGATGATGAAAATATTAACTCTCAACCATTCATGCGTTGGAGAGAA ",
    "TTTGAGTTTAAAGGTGGTTTAGACTTCTTAAAAGACGATGAGAACATCAATTCTCAACCATTCATGCGTTGGAGAGAA ",
    "ATGAGGTTTAAAGGTGGTTTAGACTTCTTAAAAGATGATGAAAACATTAACTCTCAACCATTCATGCGTTGGAGAGAA"
  ),
  Primer1 = paste0(
    "TATTCAGCTCATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTGTGGCATTGGTTGTGGCGTCCATTGATTTGGTTTTGCTGCGATCGTGCCATCCTTGGGTGGATTTTGTGTGGCATTAGGTTGTCGTGCAGG ",
    "TATTCAGCTCATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATATCTGGTGGGAGCAATGGGTCTCGCGCTTAGCGCGGGTACTTCAGTTGACTCCGGCCATCCTCGGGGAGAGCCCGTCTGGCATTAGGTTGTCGTGCAGG ",
    "TCAGCTCATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCTGGCAGGAGTGACTGACCACAAACTCTGTTTGTGAGTTGTGTCATTCTGGCCATCCTTGGTAAGATCCTTTTTGGCATTAGGTTGTCGTGCAGG ",
    "TATTCAGCTCATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTGTGGTGCGACGGATCGGTCCGACCTTTGGTGGGTACTCGATGTTGTCGCGCCATCCTTGAGTGGTTCGTCCTGGCATTAGGTTGTCGTGCAGG ",
    "TATTCAGCTCATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTGTGGCTGTCGATAGCGGCCTGTCACTTAGTGTCAGTGTTTGCTGTCGTCGCCATCCTTGGGTGGAGCCTGTGTGGCATTAGGTTGTCGTGCAGG"
  ),
  `psaA-2` = paste0(
    "CGGTCCAGCTGTTATAAATCAGTCCATACAACTGTCCATGTTGCTGTAGAAGATTCACCTGCTACTGCTGCAGCAGCTTCTACTGGATCTACACCTGGTTGTGGTGTGATACGGA ",
    "CGCTCCAGCTGTTATAAATCTGTCCAAACAACAGTCCAAGTAGCTGTTGAAGATTCACCAGCTACTGCAGCTGCAGCTTCTACTGGATCTACACCTGGTTGTGGTGTGATACGGA ",
    "CGCTCCAGCTGTTATAAATCTGTCCATACAACAGTCCAAGTTGCAGTTGAAGATTCACCTGCTACAGCAGCAGCAGCTTCTACTGGATCTACACCTGGTTGTGGTGTGATACGGA ",
    "CGTTCCACGCTGTTATAAATCTGTCCAAACAACAGTCCATGTTGCTGTTGAAGATTCACCAGCTACAGCAGCGGCAGCTTCTACTGGATCAACACCTGGTTGTGGTGTGATACGGA ",
    "CGGTCCAGCTGTTATAAATCAGTCCATACAACTGTCCAAGTTGCTGTAGAAGATTCACCAGCTACTGCTGCTGCAGCTTCTACTGGATCTACACCTGGTTGTGGTGTGATACGGA"
  )
)
species <- c("Navicula sp.", "Odontella sp.", "Skeletonema sp.", "Chaetoceros sp.", "Nitzschia sp.")
seqs <- do.call(rbind, lapply(names(seq_text), function(p) {
  s <- strsplit(seq_text[[p]], " ", fixed = TRUE)[[1L]]
  stopifnot(length(s) == 5L)
  data.frame(
    primer = rep(p, 5L), species = species, sequence = s,
    stringsAsFactors = FALSE
  )
}))

# --- Table 8: GC%, bp, dH, dS ------------------------------------------------
# 18S-V7 and ND are not amplified for all five species?  No - the paper prints
# all 5 for those; Dia18S is the one with 3 species.  Table 8 rows, verbatim.
t8 <- data.frame(
  primer = c(rep("UPA99", 5), rep("18S-V7", 5), rep("ND", 5), rep("Primer1", 5), rep("psaA-2", 5), rep("Dia18S", 3)),
  species = c(species, species, species, species, species, "Navicula sp.", "Skeletonema sp.", "Nitzschia sp."),
  gc_pct = c(42.19, 39.01, 41.67, 42.19, 40.63, 47.73, 53.26, 49.45, 45.45, 47.19, 37.18, 34.62, 34.62, 35.90, 34.62, 44.30, 51.61, 44.00, 49.01, 48.03, 47.83, 47.83, 47.83, 48.28, 47.83, 51.02, 46.00, 55.24),
  bp = c(63, 63, 59, 63, 63, 87, 91, 90, 87, 88, 77, 77, 77, 77, 77, 148, 154, 149, 150, 151, 114, 114, 114, 115, 114, 97, 99, 104),
  dh_kcal_mol = c(-510.4, -506.8, -467.8, -510.3, -511.3, -706.9, -750.4, -743.7, -708.3, -715.2, -617, -614.4, -615.2, -619, -614.4, -1219.8, -1275.6, -1216.3, -1237.2, -1232.8, -932, -934.6, -934.6, -949.6, -932, -800.7, -800.9, -865.2),
  ds_cal_k_mol = c(-1390, -1385.8, -1272.8, -1388.9, -1395.3, -1898.3, -2010.6, -1997.9, -1912.9, -1914.6, -1681.8, -1678.5, -1681.5, -1688.5, -1678.5, -3283.5, -3283.5, -3283.1, -3318.3, -3318.3, -2507, -2512.9, -2512.8, -2549.2, -2507.1, -2135.7, -2153, -2303.2),
  stringsAsFactors = FALSE
)
# The GC window is the paper's own (34.62-55.24%).  Its stated bp range
# "60 to 152" is NOT reproduced by its own Table 8: UPA99/Skeletonema prints
# 59 bp and Primer1/Odontella prints 154 bp.  The assertion therefore checks the
# printed values against the range they actually occupy; the mismatch with the
# paper's prose is reported, not patched.
stopifnot(nrow(t8) == 28L, all(t8$gc_pct >= 34.62 & t8$gc_pct <= 55.24))
stopifnot(min(t8$bp) == 59L, max(t8$bp) == 154L)

# --- Table 6: measured Tm (C), 5 primers x 5 species, ND = not detected ------
t6 <- data.frame(
  primer = c(rep("UPA99", 5), rep("18S-V7", 5), rep("ND", 5), rep("Primer1", 5)),
  species = rep(species, 4L),
  # read off Table 6 in column order Navicula/Odontella/Skeletonema/Chaetoceros/Nitzschia
  tm_c = c(
    79.03, 78.37, 79.27, 79.17, 79.17,
    81.35, 82.40, 82.85, 81.30, 82.03,
    76.80, 76.42, 76.50, 77.50, 76.57,
    83.40, 84.72, 82.43, 84.42, 84.08
  ),
  stringsAsFactors = FALSE
)
stopifnot(!anyNA(t6$tm_c), nrow(t6) == 20L)

# --- Table 9: predicted vs measured Tm --------------------------------------
t9 <- data.frame(
  primer = c(rep("psaA-2", 5), rep("Dia18S", 3)),
  species = c(species, "Navicula sp.", "Skeletonema sp.", "Nitzschia sp."),
  predicted_c = c(81.68, 82.05, 82.20, 82.78, 81.97, 83.30, 81.50, 84.75),
  measured_c = c(82.38, 82.54, 82.56, 83.05, 82.34, 84.43, 82.91, 84.14),
  difference_c = c(0.70, 0.49, 0.36, 0.27, 0.37, 1.13, 1.41, 0.61),
  stringsAsFactors = FALSE
)

# --- ASSERTIONS -------------------------------------------------------------
# 1. Table 9's difference is |measured - predicted|, as printed: row 8
#    (Dia18S, Nitzschia) is 84.14 - 84.75 = -0.61 and the paper prints 0.61.
stopifnot(all(abs(abs(t9$measured_c - t9$predicted_c) - t9$difference_c) < 0.011))
# 2. The printed "average difference 0.67" is the mean of |measured - predicted|
#    (0.6675).  The signed mean is 0.515, so the paper means the absolute one.
stopifnot(abs(mean(abs(t9$measured_c - t9$predicted_c)) - 0.67) < 0.006)
# 3. GC% recomputed from the Table 7 sequence matches Table 8.
seqs$gc_from_seq <- 100 * vapply(
  strsplit(seqs$sequence, ""),
  function(b) sum(b %in% c("G", "C")) / length(b), numeric(1)
)
seqs$seq_len <- nchar(seqs$sequence)
m <- merge(seqs, t8[, c("primer", "species", "gc_pct", "bp")], by = c("primer", "species"), all = TRUE)
# 3. The overlap must be the 20 screening amplicons; psaA-2/Dia18S sequences are
#    in Table 7 too, so the merge covers all 28 rows.
# 3. Table 7 sequences every amplicon EXCEPT the three Dia18S cells, whose
#    sequence the paper never prints; the merge must be 25 overlapping cells
#    and must reproduce the 3 Dia18S rows as missing, never as filled.
stopifnot(nrow(m) == 28L, all(is.na(m$sequence) == (m$primer == "Dia18S")))
# Table 8 GC% agrees with the Table 7 sequence on 22 of 25 sequenced cells; the
# three that do not (UPA99 Navicula/Chaetoceros 42.19 vs 41.27, Primer1
# Skeletonema 44.00 vs 43.33) are the paper's own arithmetic slips and are
# recorded, not corrected.
m$gc_slip <- !is.na(m$gc_from_seq) & abs(m$gc_from_seq - m$gc_pct) >= 0.06
stopifnot(sum(m$gc_slip) == 3L)
# 4. Table 7 length vs Table 8 bp: 23 of the 25 sequenced cells are +1, and
#    the two UPA99 cells that match exactly (Navicula, Chaetoceros) are exactly
#    the two whose GC% disagrees.  Recorded, not corrected.
stopifnot(sum(m$seq_len - m$bp == 1L, na.rm = TRUE) == 23L)
stopifnot(sum(m$seq_len == m$bp, na.rm = TRUE) == 2L)
n_len_mismatch <- sum(m$seq_len != m$bp, na.rm = TRUE)

# --- OUTPUT -----------------------------------------------------------------
hrm <- rbind(
  merge(t6, t8, by = c("primer", "species"), sort = FALSE),
  merge(
    data.frame(primer = t9$primer, species = t9$species, tm_c = t9$measured_c, stringsAsFactors = FALSE),
    t8,
    by = c("primer", "species"), sort = FALSE
  )
)
hrm$seq_len <- m$seq_len[match(paste(hrm$primer, hrm$species), paste(m$primer, m$species))]
hrm$gc_from_seq <- round(m$gc_from_seq[match(paste(hrm$primer, hrm$species), paste(m$primer, m$species))], 2)
hrm$gc_slip <- m$gc_slip[match(paste(hrm$primer, hrm$species), paste(m$primer, m$species))]
hrm$seq_len_vs_table_bp <- ifelse(is.na(hrm$seq_len), NA_integer_, hrm$seq_len - hrm$bp)
hrm$primer <- factor(hrm$primer, levels = c("UPA99", "18S-V7", "ND", "Primer1", "psaA-2", "Dia18S"))
hrm <- hrm[order(hrm$primer, match(hrm$species, species)), ]
hrm$primer <- as.character(hrm$primer)

out <- data.frame(
  chemistry = "hrm_melt",
  assay = "HRM (EvaGreen, 60->95 C at 0.1 C/s)",
  primer_pair = hrm$primer,
  species = hrm$species,
  amplicon_bp = hrm$bp,
  gc_pct = hrm$gc_pct,
  seq_len = hrm$seq_len,
  gc_from_seq = hrm$gc_from_seq,
  gc_printed_vs_recomputed_slip = hrm$gc_slip,
  seq_len_vs_table_bp = hrm$seq_len_vs_table_bp,
  dh_kcal_mol = hrm$dh_kcal_mol,
  ds_cal_k_mol = hrm$ds_cal_k_mol,
  mg_mm = NA_real_,
  dna_ng_ul = NA_real_,
  tm_measured_c = hrm$tm_c,
  outcome_definition = paste0(
    "Outcome = the ABSOLUTE melting temperature Tm (C) of ONE amplicon printed by the paper under ONE HRM run: ",
    "same EvaGreen master mix, 35 cycles, ramp 60->95 C at 0.1 C/s, Table 6 for the four screening primer pairs and ",
    "the measured column of Table 9 for psaA-2 and Dia18S. This is NOT a dTm of a genotyping call - the amplicons differ ",
    "in sequence, so the set constrains a Tm-prediction model (nearest-neighbour dH/dS + salt), not a discrimination rule. ",
    "mg_mm and dna_ng_ul are NA because the paper prints neither: Mg2+ lives inside a proprietary 2x master mix, ",
    "so this set cannot fit a salt term. gc_pct/bp/dH/dS are the paper's Table 8; seq_len is nchar of the Table 7 product ",
    "and seq_len_vs_table_bp is its disagreement with Table 8 bp (never pooled away). Table 4 (Tm vs DNA input and cycle ",
    "number) is excluded: it contradicts Table 6 for psaA-2 Navicula."
  ),
  evidence_tier = "measured",
  source_file = source_file,
  doi = doi,
  pmcid = pmcid,
  stringsAsFactors = FALSE
)

dir.create(file.path(root, "inst", "extdata"), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, file.path(root, "inst", "extdata", "readout_hrm_diatom2025.csv"), row.names = FALSE)

cat(sprintf(
  "readout_hrm_diatom2025.csv: %d rows (%d primer pairs x species + Table 9 subset); %d sequences differ in length from Table 8 bp; Mg not printed (NA).\n",
  nrow(out), length(unique(out$primer_pair)), n_len_mismatch
))
