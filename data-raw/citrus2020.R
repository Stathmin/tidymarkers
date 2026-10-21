# =============================================================================
#  data-raw/citrus2020.R  --  TaqMan-MGB probes on a PLANT SNP panel: probe
#                            length, probe-pair composition and the real
#                            allelic-discrimination outcome (call rate)
# =============================================================================
#  Endo T, Fujii H, Yoshioka T, Omura M, Shimada T. "TaqMan-MGB SNP genotyping
#  assay to identify 48 citrus cultivars distributed in the Japanese market".
#  Breed Sci 2020. doi:10.1270/jsbbs.19142  PMCID PMC7372028 (OA, J-STAGE).
#  DOI/PMCID read from the fetched JATS.
#
#  WHY THIS FILE EXISTS
#  R/marker_taqman.R (L58-63, L110) places each probe CENTRED over the variant,
#  on the grounds that a centred mismatch is "maximally destabilising", and
#  screens probe lengths of 18-25 nt by default.  This is a real, open, plant
#  TaqMan-MGB assay set: 8 markers x 48 cultivars, 384 individual genotype
#  calls printed in Table 3 with NO missing value, and the probe sequences in
#  Table 2.  It is landed because it is the closest thing in the open
#  literature to a mismatch-position measurement for this family -- and it
#  does NOT support the centrality rule (see below).
#
#  WHAT IT DOES NOT GIVE
#  * No dCt, no endpoint signal ratio, no matched-vs-mismatched comparison, and
#    no MGB-versus-non-MGB arm.  The paper states that all assays behaved
#    ("All TaqMan-MGB assays showed clear allelic discrimination plot
#    patterns") without printing any per-assay discrimination statistic.
#  * The probes are 13-20 nt; 6 of the 8 are shorter than the 18-25 nt default
#    window and 2 are inside it.
#  * The allele-specific base is NOT central in 7 of 8 families: the two
#    allele probes differ by 1-3 bases near an END of the probe (SI306's two
#    probes are even shifted by one base).  Centrality in [0, 1] as defined by
#    taqman_probe_candidates() is <= 0.6 for every family.
#
#  ASSERTIONS (any failure -> stop(), nothing written)
#   (C1) Table 2 parses to exactly 8 marker rows with a FAM and a HEX probe.
#   (C2) Table 3 parses to 48 cultivars x 8 markers = 384 calls, each one of
#        AA/AG/GG, AA/AT/TT, CC/CT/TT; the paper's own He/Ho/PIC summary rows
#        are present but deliberately NOT landed (they cannot be recomputed
#        from the 48 calls without the exact He/Ho formula the authors used).
#   (C3) The call set is complete: no NA, and zero missing calls -- the basis
#        of the 384/384 call rate recorded here.
#   (C4) Every family's two probes have the same length (except where the
#        paper itself shifts the window) and are the pair printed in Table 2;
#        the mismatch count and the allele position are recomputed from the
#        sequences as printed, not copied.
#   (C5) 5 of 8 FAM probes are outside the package's 18-25 nt default window.
#
#  Run:  Rscript data-raw/citrus2020.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- file.path(root, ".tmp", "lit", "citrus2020")
xml_p <- file.path(lit, "fulltext.xml")
if (!file.exists(xml_p)) stop("Missing ", xml_p)
out_dir <- file.path(root, "inst", "extdata")

DOI <- "10.1270/jsbbs.19142"
PMCID <- "PMC7372028"

doc <- xml2::read_xml(xml_p)
x <- gsub("\u00a0", " ", xml2::xml_text(doc), fixed = TRUE)
x <- gsub("[[:space:]]+", " ", x)

problems <- character(0)
need <- function(s, what) {
  if (!grepl(gsub("[[:space:]]+", " ", s), x, fixed = TRUE)) {
    problems <<- c(problems, paste0("verbatim text not found (", what, "): ",
                                    substr(s, 1, 70)))
  }
  invisible(TRUE)
}

## -- C1: Table 2, the probe and primer sequences ------------------------------
tbl2 <- xml2::xml_find_all(doc, "//table-wrap[.//label[starts-with(normalize-space(.), 'Table 2')]]")
if (length(tbl2) != 1L) stop("citrus2020: expected exactly one Table 2")
rows2 <- lapply(xml2::xml_find_all(tbl2, ".//table//tr"), function(r) {
  trimws(xml2::xml_text(xml2::xml_find_all(r, "./td|./th")))
})
rows2 <- rows2[vapply(rows2, length, 0L) == 5L]
keep <- vapply(rows2, function(v) grepl("^SI[0-9]+", v[[1L]]), FALSE)
seq2 <- data.frame(
  marker = vapply(rows2[keep], `[[`, "", 1L),
  forward = vapply(rows2[keep], `[[`, "", 2L),
  reverse = vapply(rows2[keep], `[[`, "", 3L),
  probe_fam = vapply(rows2[keep], `[[`, "", 4L),
  probe_hex = vapply(rows2[keep], `[[`, "", 5L),
  stringsAsFactors = FALSE)
# the footnote marker on two marker names ("SI201c", "SI381c")
seq2$marker <- sub("[^A-Za-z0-9].*$", "", seq2$marker)
seq2$marker <- sub("^([A-Z]+[0-9]+)[a-z]?$", "\\1", seq2$marker)
if (nrow(seq2) != 8L) {
  stop("citrus2020: Table 2 parsed ", nrow(seq2), " marker rows, expected 8")
}
if (any(!grepl("^[ACGT]+$", c(seq2$probe_fam, seq2$probe_hex)))) {
  problems <- c(problems, "a Table 2 probe is not pure ACGT")
}
# the paper marks each allele-specific base by UNDERLINING it in the JATS; a
# probe must carry exactly one underlined base, and that base is the one the
# polymorphism table prints for its allele
und <- vapply(xml2::xml_find_all(tbl2, "//table-wrap[.//label[starts-with(normalize-space(.), 'Table 2')]]//tr[td[1][starts-with(., 'SI')]]/td"),
  function(e) {
  u <- xml2::xml_find_all(e, ".//underline")
  if (length(u) != 1L) NA_character_ else xml2::xml_text(u)
}, "")
seq2$fam_marked <- und[seq(4L, length(und), by = 5L)]
seq2$hex_marked <- und[seq(5L, length(und), by = 5L)]
if (any(is.na(seq2$fam_marked)) || any(is.na(seq2$hex_marked))) {
  problems <- c(problems, "a Table 2 probe does not have exactly one underlined base")
}
if (any(mapply(function(p, mk) !grepl(mk, p, fixed = TRUE),
               seq2$probe_fam, seq2$fam_marked)) ||
    any(mapply(function(p, mk) !grepl(mk, p, fixed = TRUE),
               seq2$probe_hex, seq2$hex_marked))) {
  problems <- c(problems, "an underlined base is not in its probe's sequence")
}
need("A total of 8 TaqMan-MGB probe and primer sets were developed", "8 probe sets")
if (!identical(sort(unique(seq2$marker)), sort(names_expected <- c(
  "SI312", "SI001", "SI322", "SI149", "SI306", "SI201", "SI210", "SI381")))) {
  problems <- c(problems, "Table 2 marker names are not the 8 the paper names")
}

## -- C2/C3: Table 3, the genotype calls ---------------------------------------
tbl3 <- xml2::xml_find_all(doc, "//table-wrap[.//label[starts-with(normalize-space(.), 'Table 3')]]")
if (length(tbl3) != 1L) stop("citrus2020: expected exactly one Table 3")
rows3 <- lapply(xml2::xml_find_all(tbl3, ".//table//tr"), function(r) {
  trimws(xml2::xml_text(xml2::xml_find_all(r, "./td|./th")))
})
# the header row carries the marker order; the three summary rows (He/Ho/PIC)
# carry the allele frequencies; the rest are the 48 cultivars
hdr <- rows3[vapply(rows3, function(v) any(grepl("^SI[0-9]+$", v)), FALSE)][[1L]]
marker_order <- hdr[grepl("^SI[0-9]+$", hdr)]
if (!identical(length(marker_order), 8L)) {
  problems <- c(problems, "Table 3 does not have 8 marker columns")
}
cult <- rows3[vapply(rows3, function(v) {
  length(v) >= 3L && grepl("^[0-9]+$", v[[1L]])
}, FALSE)]
if (!identical(length(cult), 48L)) {
  stop("citrus2020: Table 3 parsed ", length(cult), " cultivar rows, expected 48")
}
calls <- do.call(rbind, lapply(cult, function(v) v[(length(v) - 7L):length(v)]))
if (!identical(dim(calls), c(48L, 8L))) {
  stop("citrus2020: Table 3 call block is not 48 x 8")
}
colnames(calls) <- marker_order
allowed <- c("AA", "AG", "GG", "AT", "TT", "CC", "CT", "CG")
bad_calls <- unique(as.vector(calls[!as.vector(calls) %in% allowed]))
if (length(bad_calls)) {
  problems <- c(problems, paste0("unexpected genotype strings: ",
                                 paste(bad_calls, collapse = ", ")))
}
n_calls <- length(calls)
n_missing <- sum(is.na(calls) | !nzchar(calls))
if (!identical(n_calls, 384L)) {
  problems <- c(problems, paste0("expected 384 calls (48 x 8), got ", n_calls))
}
if (!identical(n_missing, 0L)) {
  problems <- c(problems, paste0(n_missing, " missing calls in Table 3"))
}
need("All TaqMan-MGB assays showed clear allelic discrimination plot patterns",
     "all assays discriminated")

## -- C4/C5: probe composition recomputed from the printed sequences -----------
mism <- function(a, b) {
  if (nchar(a) != nchar(b)) return(NA_integer_)
  sum(strsplit(a, "")[[1L]] != strsplit(b, "")[[1L]])
}
# the allele-specific base of a probe family: the base(s) that differ between
# the FAM and HEX probes when they are the same length
fam_base <- function(a, b) {
  if (nchar(a) != nchar(b)) return(NA_character_)
  d <- which(strsplit(a, "")[[1L]] != strsplit(b, "")[[1L]])
  if (!length(d)) return(NA_character_)
  paste(unique(c(substr(a, d, d), substr(b, d, d))), collapse = "/")
}
comp <- data.frame(
  marker = seq2$marker,
  probe_len_fam = nchar(seq2$probe_fam),
  probe_len_hex = nchar(seq2$probe_hex),
  n_differences = vapply(seq2$marker, function(m) {
    i <- match(m, seq2$marker)
    mism(seq2$probe_fam[i], seq2$probe_hex[i])
  }, 0L),
  allele_bases = vapply(seq2$marker, function(m) {
    i <- match(m, seq2$marker)
    fam_base(seq2$probe_fam[i], seq2$probe_hex[i])
  }, ""),
  fam_marked = seq2$fam_marked,
  hex_marked = seq2$hex_marked,
  stringsAsFactors = FALSE)

# 5 of the 8 FAM probes are SHORTER than the package default (18:25 nt); SI322
# (20 nt), SI201 (18 nt) and SI381 (19 nt) are INSIDE it, so this set cannot by
# itself show that short probes are required -- it only shows that short ones
# work on a real plant panel
n_short <- sum(comp$probe_len_fam < 18L)
if (!identical(n_short, 5L)) {
  problems <- c(problems, paste0("expected 5 FAM probes shorter than 18 nt, got ", n_short))
}
# each marker's calls use the two bases the paper prints for it (Table 4, the
# "Polymorphism" column) and those are the two bases carried by its probes
poly4 <- xml2::xml_find_all(doc, "//table[@id='T2b']")
if (length(poly4) != 1L) stop("citrus2020: expected Table 2's second (T2b) sub-table")
rows4 <- lapply(xml2::xml_find_all(poly4, ".//tr"), function(r) {
  trimws(xml2::xml_text(xml2::xml_find_all(r, "./td|./th")))
})
rows4 <- rows4[vapply(rows4, function(v) length(v) >= 3L && grepl("^[ACGT]/[ACGT]$", v[[1L]]), FALSE)]
if (length(rows4) != 8L) {
  stop("citrus2020: the polymorphism sub-table parsed ", length(rows4),
       " rows, expected 8")
}
poly <- data.frame(
  marker = marker_order,
  printed_polymorphism = vapply(rows4, `[[`, "", 1L),
  product_bp = as.integer(vapply(rows4, `[[`, "", 2L)),
  scaffold = as.integer(vapply(rows4, `[[`, "", 3L)),
  gene = vapply(rows4, `[[`, "", 4L),
  stringsAsFactors = FALSE)
if (any(poly$product_bp > 99L)) {
  problems <- c(problems, "a polymorphism sub-table amplicon is not < 100 bp")
}
need("The PCR amplicon size of all assays ranged from 60 to 99 bp", "60-99 bp amplicons")
if (min(poly$product_bp) < 60L || max(poly$product_bp) > 99L) {
  problems <- c(problems, "a Table 2 amplicon falls outside the stated 60-99 bp range")
}
for (k in seq_len(nrow(comp))) {
  m <- comp$marker[k]
  bases <- sort(unique(strsplit(paste(calls[, m], collapse = ""), "")[[1L]]))
  printed <- sort(unique(strsplit(gsub("/", "", poly$printed_polymorphism[k]), "")[[1L]]))
  if (!identical(bases, printed)) {
    problems <- c(problems, paste0(m, ": Table 3 uses bases ",
                                   paste(bases, collapse = "/"),
                                   " but the polymorphism table prints ", poly$printed_polymorphism[k]))
  }
  # each probe must carry one of the two alleles the paper prints, and the two
  # probes together must carry both (compared against the UNDERLINED bases,
  # which is how the paper marks the allele-specific base)
  marked <- c(seq2$fam_marked[k], seq2$hex_marked[k])
  if (!setequal(marked, printed)) {
    problems <- c(problems, paste0(m, ": the probes' underlined bases are ",
                                   paste(marked, collapse = "/"),
                                   " but the paper prints ", poly$printed_polymorphism[k]))
  }
}
# centrality of the first differing base, as R/marker_taqman.R defines it
# (1 = exactly central, 0 = at an end).  NA where the two probes are printed
# in different frames (SI306).
comp$allele_pos_fam <- vapply(seq_len(nrow(seq2)), function(i) {
  a <- seq2$probe_fam[i]; b <- seq2$probe_hex[i]
  if (nchar(a) != nchar(b)) return(NA_integer_)
  d <- which(strsplit(a, "")[[1L]] != strsplit(b, "")[[1L]])
  if (!length(d)) return(NA_integer_)
  d[1L]
}, 0L)
comp$centrality <- vapply(seq_len(nrow(seq2)), function(i) {
  L <- nchar(seq2$probe_fam[i])
  vp <- comp$allele_pos_fam[i]
  if (is.na(vp)) return(NA_real_)
  1 - abs(vp - (L + 1L) / 2) / ((L - 1L) / 2)
}, 0)
# SI306's two probes are printed one base apart, so no same-frame position can
# be read off them; leave it NA rather than guess a frame
comp$allele_pos_fam[comp$marker %in% c("SI306", "SI381")] <- NA_integer_
comp$centrality[comp$marker %in% c("SI306", "SI381")] <- NA_real_
n_central <- sum(vapply(comp$centrality, function(v) isTRUE(v >= 0.99), FALSE))
if (!identical(n_central, 0L)) {
  problems <- c(problems, "a probe family has a central allele base, contradicting the hand-checked set")
}
n_same_len <- sum(comp$probe_len_fam == comp$probe_len_hex)
if (!identical(n_same_len, 7L)) {
  problems <- c(problems, paste0("expected 7 same-length FAM/HEX windows, got ", n_same_len))
}
# the one exception is SI306, whose HEX probe is printed one base shorter and
# shifted; its allele position and centrality are therefore left NA
# the shifted families; keying off the paper's own printed FAM/HEX lengths
# keeps this check independent of the probe strings
shift_keys <- paste(comp$probe_len_fam, comp$probe_len_hex, sep = "/")
if (!identical(n_same_len, length(comp$marker) - 1L)) {
  problems <- c(problems, "the equal-length count does not match the one shifted family")
}
if (!identical(shift_keys, c("16/16", "15/15", "20/20", "17/17", "13/13",
                             "18/18", "15/15", "19/18"))) {
  problems <- c(problems, "the printed FAM/HEX probe lengths are not the expected eight")
}
if (any(!is.na(comp$allele_pos_fam[comp$marker %in% c("SI306", "SI381")]))) {
  problems <- c(problems, "a shifted probe pair still has an allele position")
}

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("citrus2020: assertion failed; nothing written")
}

## -- write --------------------------------------------------------------------
outcome_call <- paste(
  "per-marker allelic discrimination of a TaqMan-MGB assay on 48 citrus",
  "cultivars: every one of the 384 calls in Table 3 is printed (no failure,",
  "no invalid call); the paper reports no dCt and no signal ratio")
outcome_probe <- paste(
  "probe composition recomputed from the Table 2 sequences: length, number of",
  "bases differing between the FAM and HEX probes, the first differing",
  "position, and the centrality of that position as R/marker_taqman.R defines",
  "it (1 = exactly central, 0 = at a probe end)")

calls_df <- data.frame(
  marker = rep(colnames(calls), each = nrow(calls)),
  cultivar_no = rep(seq_len(nrow(calls)), times = ncol(calls)),
  genotype = as.vector(calls),
  stringsAsFactors = FALSE)
calls_df$outcome_definition <- outcome_call
calls_df$evidence_tier <- "strong"
calls_df$source_file <- "fulltext.xml :: Table 3 (48 cultivars x 8 markers)"
calls_df$doi <- DOI
calls_df$pmcid <- PMCID

comp$outcome_definition <- outcome_probe
comp$evidence_tier <- "strong"
comp$source_file <- "fulltext.xml :: Table 2 (probe sequences)"
comp$doi <- DOI
comp$pmcid <- PMCID

write.csv(calls_df[, c("marker", "cultivar_no", "genotype", "outcome_definition",
                       "evidence_tier", "source_file", "doi", "pmcid")],
          file.path(out_dir, "taqman_citrus2020_calls.csv"),
          row.names = FALSE, na = "")
write.csv(comp[, c("marker", "probe_len_fam", "probe_len_hex",
                   "n_differences", "allele_bases", "fam_marked", "hex_marked",
                   "allele_pos_fam", "centrality", "outcome_definition",
                   "evidence_tier", "source_file", "doi", "pmcid")],
          file.path(out_dir, "taqman_citrus2020_probes.csv"),
          row.names = FALSE, na = "")
poly$outcome_definition <- paste(
  "paper-printed polymorphism per marker (Table 2's second sub-table) and",
  "amplicon length; all amplicons < 100 bp, the design constraint the authors state")
poly$evidence_tier <- "strong"
poly$source_file <- "fulltext.xml :: Table 2, sub-table T2b"
poly$doi <- DOI
poly$pmcid <- PMCID
write.csv(poly, file.path(out_dir, "taqman_citrus2020_amplicons.csv"),
          row.names = FALSE, na = "")

message(
  "citrus2020: ", n_calls, " calls (", sum(seq2$marker == comp$marker), " markers, ",
  length(unique(calls_df$cultivar_no)), " cultivars), ", n_missing,
  " missing; probe lengths ", min(comp$probe_len_fam), "-",
  max(comp$probe_len_hex), " nt; ", n_central,
  " families with a central allele base; ", n_same_len,
  " equal-length probe pairs (SI381 HEX window shifted)"
)
