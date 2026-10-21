# =============================================================================
#  data-raw/grewal2020-supplement.R  --  Y3.2's two-axis homeolog rule, measured
# =============================================================================
#  Grewal S, Hubbart-Edwards S, Yang C, et al. "Rapid identification of
#  homozygosity and site of wild relative introgressions in wheat through
#  chromosome-specific KASP genotyping assays." Plant Biotechnol J 2020.
#  PMID 31465620, OA (CC-BY), PMCID PMC7006891.
#  DOI is read from tidymarkers.bib at run time so it cannot drift.
#
#  WHY THIS FILE EXISTS
#  Y3.2 claimed the AS primer's 3' end separates A/B/D. The reviewer objected:
#  the AS 3' base carries the ALLELE; subgenome specificity comes from the
#  COMMON primer's 3' end sitting on a homeolog-diagnostic base (PolyMarker).
#  Grewal's Fig 1b legend states exactly this, and Data S5 has 1000 real wheat
#  KASP assays to test our reading of it against.
#
#  WHAT DATA S5 CONTAINS (and the three traps it sets)
#  Every row names Primer_AlleleFAM and Primer_AlleleHEX, but the primers are
#  built from three different determinants:
#    418  two-allele     FAM 3' == allele_1,           HEX 3' == allele_2
#     91  wild-relative  FAM 3' == allele_1,           HEX 3' == complement(allele_1)
#    491  reverse        FAM 3' == complement(allele_1), HEX 3' == complement(allele_2)
#  Trap 1: a checker that only inspects the FAM primer calls all 491 reverse
#    assays "forward" and then blames the HEX primer for a base mismatch.
#  Trap 2: for transition SNPs (C/T, A/G) complement(allele) == allele, so the
#    reverse mode is invisible unless the PURE transversions are examined.
#  Trap 3: 352 of the 418 two-allele pairs have DIFFERENT lengths, because the
#    HEX allele primer carries an extra 5' tail. Asserting "the pair differs at
#    exactly one base" silently restricts you to 66 rows. The invariant that
#    survives is on the 3' BASE, not on the string edit distance.
#
#  THE ASSERTIONS
#  (1) Data S5 holds exactly 1000 assays.
#  (2) The three determinants are disjoint and together cover all 1000.
#  (3) Where the two AS primers have equal length (66 of them), they differ at
#      exactly one position, and it is the 3'-terminal one.
#  (4) NOT ASSERTED: that the common primer's 3' base encodes the subgenome.
#      Against the SNP alleles it coincides at chance rate (455/1000, chance
#      500/1000), so such an assertion cannot fail and is not a test. The
#      Fig. 1b rule is quoted verbatim, not verified.
#
#  Run:  Rscript data-raw/grewal2020-supplement.R [--dir .tmp/lit/grewal2020]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!exists("root")) {
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else getwd()
  setwd(root)
}
suppressPackageStartupMessages(library(xml2))
source(file.path(root, "data-raw", "lib", "xlsx-read.R"))

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "grewal2020")
s5 <- file.path(lit, "supplements", "PBI-18-743-s005.xlsx")
ft <- file.path(lit, "fulltext.xml")
for (f in c(s5, ft)) if (!file.exists(f)) {
  stop("missing ", f, "\n  fetch with: Rscript data-raw/fetch-sources.R")
}

DOI <- local({
  txt <- readLines(file.path(root, "inst", "extdata", "tidymarkers.bib"), warn = FALSE)
  i <- grep("^@[a-z]+\\{grewal2020,", txt)
  if (length(i) != 1L) stop("grewal2020 not in tidymarkers.bib; run fetch-sources.R")
  k <- grep("doi = \\{", txt[seq(i, length(txt))])[[1L]]
  sub(".*doi = \\{([^}]*)\\}.*", "\\1", txt[i - 1L + k])
})

## --- full text: pin the rule to the authors' words --------------------------
body <- xml_text(xml_find_all(read_xml(ft), "//body//p"))
rule <- grep("3\\W{0,2}end anchored to a base, which is unique to the target subgenome",
             body, value = TRUE)
if (!length(rule)) stop("anchored-common-primer sentence not found in the full text")

## --- Data S5 ---------------------------------------------------------------
rows <- xlsx_sheets(s5)[["Sheet1"]]
hdr <- which(vapply(rows, function(r) identical(xcell(r, "A"), "ID"), TRUE))
if (length(hdr) != 1L) stop("expected one header row starting with ID")
rec <- lapply(rows[seq(hdr[[1L]] + 1L, length(rows))], function(r) {
  id <- xcell(r, "A")
  if (is.na(id)) return(NULL)
  data.frame(id = id, fam = xcell(r, "B"), hex = xcell(r, "C"),
             common = xcell(r, "D"), allele_fam = xcell(r, "E"),
             allele_hex = xcell(r, "F"), seq = xcell(r, "G"),
             stringsAsFactors = FALSE)
})
kasp <- do.call(rbind, rec)
kasp <- kasp[!is.na(kasp$id), ]
if (nrow(kasp) != 1000L) {
  stop("Data S5 should hold 1000 KASP assays, parsed ", nrow(kasp))
}

end3 <- function(x) substring(x, nchar(x))
kasp$as_fam_3p <- end3(kasp$fam)
kasp$as_hex_3p <- end3(kasp$hex)
kasp$common_3p <- end3(kasp$common)

tag <- regmatches(kasp$seq, regexpr("\\[[ACGT]/[ACGT]\\]", kasp$seq))
if (length(tag) != nrow(kasp)) stop("some Data S5 rows lack an [X/Y] SNP tag")
kasp$snp_allele_1 <- substr(tag, 2L, 2L)
kasp$snp_allele_2 <- substr(tag, 4L, 4L)

rc1 <- function(x) chartr("ACGT", "TGCA", x)

## --- the three determinants, each defined by the bases ---------------------
is_wild <- kasp$as_fam_3p == kasp$snp_allele_1 &
  kasp$as_hex_3p == unname(rc1(kasp$snp_allele_1))
is_two_allele <- !is_wild & kasp$as_fam_3p == kasp$snp_allele_1 &
  kasp$as_hex_3p == kasp$snp_allele_2
is_reverse <- !is_wild & !is_two_allele &
  kasp$as_fam_3p == unname(rc1(kasp$snp_allele_1)) &
  kasp$as_hex_3p == unname(rc1(kasp$snp_allele_2))
if (!all(is_wild | is_two_allele | is_reverse)) {
  stop(sum(!(is_wild | is_two_allele | is_reverse)),
       " assays fit no known determinant; the axis reading is wrong")
}
if (any(is_two_allele & (is_wild | is_reverse))) stop("assay modes are not disjoint")
kasp$assay_mode <- ifelse(is_reverse, "reverse",
                     ifelse(is_wild, "wild-relative", "two-allele"))

## --- assertion (3): equal-length pairs differ only at the 3' base ----------
two <- kasp$assay_mode == "two-allele"
same_len <- nchar(kasp$fam) == nchar(kasp$hex)
check <- which(two & same_len)
n_diff <- vapply(check, function(i) {
  sum(strsplit(kasp$fam[i], "")[[1]] != strsplit(kasp$hex[i], "")[[1]])
}, integer(1))
dif_pos <- vapply(check, function(i) {
  a <- strsplit(kasp$fam[i], "")[[1]]
  b <- strsplit(kasp$hex[i], "")[[1]]
  which(a != b)
}, integer(1))
if (any(n_diff != 1L) || any(dif_pos != nchar(kasp$fam[check]))) {
  stop("equal-length two-allele pairs are not single 3'-terminal substitutions")
}

out <- data.frame(
  id = kasp$id,
  snp_allele_1 = kasp$snp_allele_1, snp_allele_2 = kasp$snp_allele_2,
  snp_class = ifelse(xor(kasp$snp_allele_1 %in% c("A", "G"),
                         kasp$snp_allele_2 %in% c("A", "G")),
                     "transversion", "transition"),
  assay_mode = kasp$assay_mode,
  as_fam_3p = kasp$as_fam_3p, as_hex_3p = kasp$as_hex_3p,
  fam_len = nchar(kasp$fam), hex_len = nchar(kasp$hex),
  common_3p = kasp$common_3p,
  doi = DOI,
  stringsAsFactors = FALSE
)
write.csv(out, file.path(root, "inst", "extdata", "grewal_kasp_axes.csv"),
          row.names = FALSE, na = "")

write.csv(data.frame(
  axis = c("allele", "subgenome"),
  carried_by = c("3'-terminal base of each allele-specific primer",
                 "3'-terminal base of the common primer (anchor)"),
  source = c("Grewal 2020 Data S5", "Grewal 2020 Fig. 1b legend"),
  verbatim = c(
    paste0(sum(is_two_allele), " two-allele (FAM 3'=allele_1, HEX 3'=allele_2), ",
           sum(is_wild), " wild-relative (HEX 3'=complement(allele_1)), ",
           sum(is_reverse), " reverse-strand (both complements)"),
    rule[[1L]]),
  doi = DOI
), file.path(root, "inst", "extdata", "kasp_homeolog_axes.csv"),
row.names = FALSE, na = "")

message(
  "grewal2020 [", DOI, "]: ", nrow(kasp), " assays = ", sum(is_two_allele),
  " two-allele + ", sum(is_wild), " wild-relative + ", sum(is_reverse),
  " reverse; ", length(check), " equal-length two-allele pairs, all single ",
  "3'-terminal substitutions; ", sum(!same_len & two),
  " two-allele pairs differ in length (5' tail) and are excluded on purpose"
)
