# =============================================================================
#  data-raw/hrmplant.R  --  HRM in plants: the substitution x context grid, the
#                           measured homozygote dTm panel, and what is actually
#                           open about the Liew 2004 numbers
# =============================================================================
#  READ FIRST. This script lands five tables and refuses two sources.
#
#  liew2004 (Clin Chem 2003;50:230-231? no: 2004, doi:10.1373/clinchem.2004.032136)
#  is CLOSED ACCESS (Europe PMC core: isOpenAccess N, inEPMC N, no PMCID).
#  Nothing here reads it and no paywall was bypassed. Its bands (0.8-1.4 C for
#  A::T<->G::C, <=0.4 C otherwise, ~4% of base-pair-neutral variants
#  unresolvable by nearest-neighbour symmetry) are landable ONLY as a SECONDARY
#  CITATION restated by an open paper, and every row of that CSV says so.
#
#  REFUSED 1: Reja et al., Clin Chem 2014, PMID 24718912. Closed access AND no
#  Europe PMC record (EXT_ID:24718912 returns an unrelated record). No file in
#  .tmp/lit, so nothing is landed from it.
#  REFUSED 2: crop HRM resolvability counts. The TITLE/ABSTRACT searches
#  logged in .tmp/sdd/report-hrmplant.md found no open crop paper reporting
#  per-marker resolvability against a truth set. Nothing landed.
#
#  WHAT IS OPEN (every number re-read from the named file):
#
#  A1  jsbbs2018 Supplemental Table 3 -> |ddG3bp| (kcal/mol) for both class I
#      bimorphisms (C/T and G/A) in all 16 flank contexts: 96 rows.
#  A2  jsbbs2018 Supplemental Table 1 -> the plant design frame: 70 numbered
#      marker rows (NOT 96 -- the sheet has blank rows removed here).
#  A3  jsbbs2018 Supplemental Tables 2, 4, 5 -> measured homozygote dTmobs for
#      38 + 14 + 18 marker rows, and the dTmobs-by-bimorphism summary.
#  A4  gundry2008 Tables 3 and 4 + prose -> 16 homozygote Tm rows from 8
#      amplicons, the 6-locus accuracy table (n sums to 1147), and the Liew
#      restatement quoted from the open text.
#  A5  liew2004 restated bands -> one CSV, evidence_tier "secondary_cited".
#
#  Deliberately NOT landed: the 1105-amplicon aggregate (not re-derivable from
#  the tables: sample x duplicate x target structure is not published);
#  BX647987 "resolution" (the paper reports a REGRESSION, 69 -> 64.3%);
#  Liew's per-substitution table (closed); Reja (no file); crop resolvability
#  (no qualifying paper).
#
#  Run:  Rscript data-raw/hrmplant.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

xml_j <- fs::path(root, ".tmp", "lit", "jsbbs2018", "jsbbs2018_fulltext.xml")
xlsx_j <- fs::path(root, ".tmp", "lit", "jsbbs2018", "supplements", "68_18048_2.xlsx")
xml_g <- fs::path(root, ".tmp", "lit", "gundry2008", "fulltext.xml")
for (p in c(xml_j, xlsx_j, xml_g)) if (!file.exists(p)) stop("Missing ", p)
out_dir <- fs::path(root, "inst", "extdata")

DOI_J <- "10.1270/jsbbs.18048"
DOI_G <- "10.1093/nar/gkn204"
DOI_L <- "10.1373/clinchem.2004.032136"

# -- helpers ------------------------------------------------------------------
norm <- function(x) {
  x <- gsub("\u2013|\u2014|\u2212", "-", x)   # en dash, em dash, minus
  x <- gsub("\u2032", "'", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}
paras_of <- function(doc) {
  vapply(xml2::xml_find_all(doc, "//*[local-name()='p']"),
         function(p) norm(xml2::xml_text(p)), character(1))
}
one_para <- function(paras, needle) {
  hit <- paras[grepl(needle, paras, fixed = TRUE)]
  if (length(hit) != 1L)
    stop("Expected exactly 1 paragraph containing '", needle, "', found ", length(hit))
  hit
}

source("data-raw/lib/xlsx-read.R")
sh <- xlsx_sheets(xlsx_j)
g <- function(r, k) {
  v <- r[[k]]
  if (is.null(v) || is.na(v) || !nzchar(trimws(v))) NA_character_ else trimws(v)
}
N <- function(r, k) suppressWarnings(as.numeric(g(r, k)))
bimorph_of <- function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_character_)
  m <- sort(c(substr(a, 2L, 2L), substr(b, 2L, 2L)))
  paste(m, collapse = "/")
}

# =============================================================================
# A1  jsbbs2018 Supplemental Table 3: Class I only (48 C/T + 48 G/A cells)
# (the sheet then repeats the same 4x4 grid for Class II G/T, Class III C/G and
#  Class IV A/T, also 96 rows; those are not base-pair neutral and are not used)
# =============================================================================
rows3 <- sh[["Supplemental_Table3"]]
# The sheet prints: Class I C/T (16), Class I G/A (16), then the SAME 4x4 grid
# repeated a second time with flipped 5'/3' labels (Class II G/T, III C/G and
# IV A/T interleaved). 96 unique (class, 5', 3') cells; the other 96 rows are
# the transposed duplicate. Keep the FIRST occurrence only and assert the
# duplicate to prove the reader is not double-counting.
cell3_all <- Filter(function(r) grepl("^Class", paste(g(r, "A"))) && !is.na(N(r, "J")), rows3)
key3 <- function(rs) vapply(rs, function(r) paste(g(r, "A"), g(r, "B"), g(r, "C"), g(r, "D")), character(1))
cell3 <- cell3_all[!duplicated(key3(cell3_all))]
if (length(cell3_all) - length(cell3) != 96L)
  stop("Supplemental Table 3: expected 96 transposed duplicate rows, got ",
       length(cell3_all) - length(cell3))
if (length(cell3) != 96L) stop("Supplemental Table 3 should have 96 Class I rows, got ", length(cell3))
grid <- tibble::tibble(
  snp_class = vapply(cell3, function(r) g(r, "A"), character(1)),
  nn5 = vapply(cell3, function(r) g(r, "B"), character(1)),
  bimorphism = vapply(cell3, function(r) g(r, "C"), character(1)),
  nn3 = vapply(cell3, function(r) g(r, "D"), character(1)),
  ref_3mer = vapply(cell3, function(r) g(r, "E"), character(1)),
  dg_ref = vapply(cell3, function(r) N(r, "F"), numeric(1)),
  alt_3mer = vapply(cell3, function(r) g(r, "H"), character(1)),
  dg_alt = vapply(cell3, function(r) N(r, "I"), numeric(1)),
  ddg_3bp = vapply(cell3, function(r) N(r, "J"), numeric(1)),
  best_combo = vapply(cell3, function(r) g(r, "K"), character(1)))
grid$ddg_recomputed <- grid$dg_alt - grid$dg_ref
grid$snp_class_num <- ifelse(grid$snp_class == "Class I", 1L, NA_integer_)
grid$evidence_tier <- "model"
grid$source_file <- "68_18048_2.xlsx :: Supplemental_Table3"
grid$doi <- DOI_J
grid$outcome_definition <- paste(
  "Outcome = |ddG3bp| (kcal/mol), a NEAREST-NEIGHBOUR FREE-ENERGY INDEX for one",
  "homozygote comparison: a 3-bp window (5' flank, SNP, 3' flank) scored by ref and",
  "alt triplet free energy. evidence_tier 'model': these are NN calculations from",
  "the paper, not measured Tm. This is NOT a melting temperature and NOT a dTm",
  "threshold. ddg_recomputed = |dg_alt - dg_ref| re-derives the printed column:",
  "the sheet carries the SIGN inconsistently (Class I prints dg_alt - dg_ref,",
  "Classes III/IV print dg_ref - dg_alt), so only the magnitude is additive.")
grid$ddg_recomputed <- abs(grid$dg_alt - grid$dg_ref)
bad <- abs(grid$ddg_recomputed - grid$ddg_3bp) > 1e-9
if (any(bad)) stop("Table 3 is not additive in row(s) ", paste(which(bad), collapse = ", "))
if (!all(paste0(grid$nn5, substr(grid$bimorphism, 1L, 1L), grid$nn3) == grid$ref_3mer))
  stop("Table 3 ref 3-mer != 5'flank + allele1 + 3'flank")
if (!all(paste0(grid$nn5, substr(grid$bimorphism, 3L, 3L), grid$nn3) == grid$alt_3mer))
  stop("Table 3 alt 3-mer != 5'flank + allele2 + 3'flank")
if (anyDuplicated(paste(grid$bimorphism, grid$nn5, grid$nn3)))
  stop("Table 3 has duplicated (bimorphism, 5', 3') cells")
# The 96 unique cells are 6 bimorphisms x 16, NOT the 2 that the Class I label
# suggests: the sheet prints all four substitution classes. Only C/T and G/A are
# base-pair neutral; C/G, G/T, A/T and C/A are not (kept as a column, not dropped).
tb <- table(grid$bimorphism)
if (length(tb) != 6L || any(tb != 16L))
  stop("Table 3 is not 6 bimorphisms x 16 cells, got ", paste(names(tb), tb, collapse = "/"))
if (!setequal(names(tb), c("C/T", "G/A", "C/G", "G/T", "A/T", "C/A")))
  stop("unexpected bimorphism labels: ", paste(names(tb), collapse = "/"))
bp_neutral <- c("C/T", "G/A")
if (sum(tb[bp_neutral]) != 32L) stop("base-pair-neutral share is not 32 of 96")
# every (bimorphism, 5') group must have all 4 3' contexts
grid1 <- grid[grid$snp_class == "Class I", ]
grp <- split(grid1, paste(grid1$bimorphism, grid1$nn5))
if (length(grp) != 8L || any(vapply(grp, nrow, integer(1)) != 4L))
  stop("Table 3 is not 8 groups of 4")
# winners: one per 5' group, 8 distinct 3-bp targets
win <- grid1[!is.na(grid1$best_combo), ]
if (nrow(win) != 8L || anyDuplicated(win$best_combo))
  stop("Table 3 should name 8 distinct winning 3-bp targets, got ", nrow(win))
for (i in seq_len(nrow(win))) {
  gr <- grid1[grid1$bimorphism == win$bimorphism[i] & grid1$nn5 == win$nn5[i], ]
  if (gr$best_combo[which.max(gr$ddg_3bp)] != win$best_combo[i])
    stop("Table 3: printed best combo is not the max |ddG| of its 5' group: ", win$best_combo[i])
}
# 3' context moves every C/T 5' group (the grid is not a pure 5' ranking)
ctrl <- vapply(split(grid1$ddg_3bp[grid1$bimorphism == "C/T"],
                     grid1$nn5[grid1$bimorphism == "C/T"]), function(v) diff(range(v)), numeric(1))
if (any(ctrl <= 0.1)) stop("Table 3: a C/T 5' group does not move with 3' context")
readr::write_csv(grid[, c("snp_class", "nn5", "bimorphism", "nn3", "ref_3mer", "dg_ref",
                          "alt_3mer", "dg_alt", "ddg_3bp", "ddg_recomputed", "best_combo",
                          "evidence_tier", "source_file", "doi", "outcome_definition")],
                 fs::path(out_dir, "hrm_jsbbs2018_context_grid.csv"))

# =============================================================================
# A2  jsbbs2018 Supplemental Table 1: the design frame (70 numbered rows)
# =============================================================================
rows1 <- sh[["Supplemental_Table1"]]
ids1 <- vapply(rows1, function(r) paste(g(r, "B")), character(1))
ids1 <- ids1[grepl("^[0-9]+$", ids1)]
if (length(ids1) != 70L) stop("Supplemental Table 1 should have 70 marker rows, got ", length(ids1))
# One contiguous numbering 1..80: Exp. 1 = 1:48 (10 drop out), Exp. 2 = 49:80.
if (!identical(sort(as.integer(ids1)), c(1:8, 10:21, 23:28, 30:32, 34L, 36:38, 40L, 42L,
                                         44:45, 48:80)))
  stop("Supplemental Table 1 ids are not 1:80 minus the ten absent Exp. 1 markers")
miss <- setdiff(as.character(1:48), ids1)
if (!identical(miss, c("9", "22", "29", "33", "35", "39", "41", "43", "46", "47")))
  stop("Supplemental Table 1: expected absent markers 1:48 minus {9,22,29,33,35,39,41,43,46,47}")
# The sheet carries NO M##_TF rows: every id is a plain integer, Exp. 2 is 49:80.
if (any(grepl("^M", ids1))) stop("Supplemental Table 1 unexpectedly contains M-prefixed ids")
if (sum(as.integer(ids1) >= 49L) != 32L || sum(as.integer(ids1) <= 48L) != 38L)
  stop("Supplemental Table 1: Exp. 2 is not 32 (49:80) with 38 remaining in Exp. 1")
exp1 <- ifelse(as.integer(ids1) >= 49L, "Exp. 2", "Exp. 1")
frame <- tibble::tibble(
  marker_id = ids1,
  experiment = exp1,
  n_designed = 70L,
  n_absent_from_1_48 = 10L,
  evidence_tier = "measured",
  source_file = "68_18048_2.xlsx :: Supplemental_Table1",
  doi = DOI_J,
  outcome_definition = paste(
    "Outcome = the set of SNP markers the soybean (Glycine max) study designed HRM",
    "primers for. 70 numbered rows: 38 in Exp. 1 (ids 1:48 less ten that never reach",
    "the sheet) and 32 in Exp. 2 (ids 49:80). There is no Exp. 3 and no M##_TF block.",
    "This is a DESIGN COUNT, not a success count: per-experiment denominators live in",
    "hrm_jsbbs2018_funnel.csv."))
readr::write_csv(frame, fs::path(out_dir, "hrm_jsbbs2018_table1_frame.csv"))

# =============================================================================
# A3  jsbbs2018 measured homozygote dTm: Tables 2, 4, 5 + by-bimorphism summary
# =============================================================================
t2 <- Filter(function(r) grepl("^[0-9]+$", paste(g(r, "A"))), sh[["Supplemental_Table2"]])
if (length(t2) != 38L) stop("Table 2: expected 38 markers, got ", length(t2))
m2 <- tibble::tibble(
  marker_id = vapply(t2, function(r) g(r, "A"), character(1)),
  product_bp = vapply(t2, function(r) N(r, "D"), numeric(1)),
  ref_triplet = vapply(t2, function(r) g(r, "E"), character(1)),
  alt_triplet = vapply(t2, function(r) g(r, "F"), character(1)),
  dtm_obs = vapply(t2, function(r) N(r, "I"), numeric(1)),
  ddg_3bp = vapply(t2, function(r) N(r, "J"), numeric(1)))
m2$bimorphism <- mapply(bimorph_of, m2$ref_triplet, m2$alt_triplet, USE.NAMES = FALSE)
if (any(is.na(m2$bimorphism))) stop("Table 2 triplet does not yield a bimorphism")

t4 <- Filter(function(r) grepl("^M[0-9]+_AF$", paste(g(r, "A"))), sh[["Supplemental_Table4"]])
if (length(t4) != 14L) stop("Table 4: expected 14 pairs, got ", length(t4))
m4 <- tibble::tibble(
  marker_nonnns = vapply(t4, function(r) g(r, "A"), character(1)),
  marker_nnns = vapply(t4, function(r) g(r, "B"), character(1)),
  product_bp = vapply(t4, function(r) N(r, "I"), numeric(1)),
  dtm_nonnns = vapply(t4, function(r) N(r, "T"), numeric(1)),
  dtm_nnns = vapply(t4, function(r) N(r, "U"), numeric(1)))
both4 <- !is.na(m4$dtm_nonnns) & !is.na(m4$dtm_nnns)
if (sum(both4) != 11L) stop("Table 4: expected 11 pairs with both dTm, got ", sum(both4))
if (!identical(vapply(t4, function(r) N(r, "V"), numeric(1)),
               ifelse(both4, m4$dtm_nnns - m4$dtm_nonnns, NA_real_)) %in% c(TRUE))
  stop("Table 4: delta column != NNNs - non-NNNs")

t5 <- Filter(function(r) grepl("^M[0-9]+_TF$", paste(g(r, "A"))), sh[["Supplemental_Table5"]])
if (length(t5) != 18L) stop("Table 5: expected 18 markers, got ", length(t5))
m5 <- tibble::tibble(
  marker_id = vapply(t5, function(r) g(r, "A"), character(1)),
  product_bp = vapply(t5, function(r) N(r, "F"), numeric(1)),
  ref_triplet = vapply(t5, function(r) g(r, "D"), character(1)),
  alt_triplet = vapply(t5, function(r) g(r, "E"), character(1)),
  dtm_obs = vapply(t5, function(r) N(r, "M"), numeric(1)))
m5$bimorphism <- mapply(bimorph_of, m5$ref_triplet, m5$alt_triplet, USE.NAMES = FALSE)

meas <- dplyr::bind_rows(
  tibble::tibble(experiment = "Exp. 1", marker_id = m2$marker_id,
                 product_bp = m2$product_bp, bimorphism = m2$bimorphism,
                 dtm_obs = m2$dtm_obs, ddg_3bp = m2$ddg_3bp),
  tibble::tibble(experiment = "Exp. 3", marker_id = m5$marker_id,
                 product_bp = m5$product_bp, bimorphism = m5$bimorphism,
                 dtm_obs = m5$dtm_obs, ddg_3bp = NA_real_))
meas$has_nnn <- experiment <- NULL
meas$has_nnn <- c(rep(FALSE, 38L), rep(TRUE, 18L))
if (nrow(meas) != 56L) stop("Expected 38 + 18 = 56 measured single-SNP rows")
if (any(meas$product_bp < 54 | meas$product_bp > 100))
  stop("measured products fall outside the paper's 54-100 bp window")

cl <- split(meas$dtm_obs, meas$bimorphism)
summ <- tibble::tibble(
  bimorphism = names(cl),
  n = vapply(cl, length, integer(1)),
  dtm_min = vapply(cl, min, numeric(1)),
  dtm_mean = vapply(cl, mean, numeric(1)),
  dtm_max = vapply(cl, max, numeric(1)))
if (sum(summ$n) != 56L) stop("bimorphism summary does not sum to 56")
if (identical(summ$bimorphism[which.max(summ$dtm_mean)], "C/T") == FALSE)
  stop("C/T is not the highest mean dTm class")
if (!identical(summ$bimorphism[which.min(summ$dtm_mean)], "C/G"))
  stop("C/G is not the lowest mean dTm class")
summ$largest_ddg_rank <- NA_integer_
ddg_by <- vapply(split(m2$ddg_3bp, m2$bimorphism), mean, numeric(1))
summ$ddg_mean_exp1 <- unname(ddg_by[summ$bimorphism])
if (names(ddg_by)[which.max(ddg_by)] != "C/T")
  stop("C/T is not the largest mean |ddG3bp| in Exp. 1, got ", names(ddg_by)[which.max(ddg_by)])
if (names(ddg_by)[which.min(ddg_by)] != "C/G")
  stop("C/G is not the smallest mean |ddG3bp| in Exp. 1, got ", names(ddg_by)[which.min(ddg_by)])
if (length(ddg_by) != 5L || names(sort(ddg_by)[1L]) != "C/G" ||
    names(sort(ddg_by, decreasing = TRUE)[2L]) != "A/G")
  stop("Exp. 1 class ordering is not C/T > A/G > A/C > G/T > C/G")
# G/T is 2nd-lowest on dTm but is NOT the 2nd on ddG -- the two indices disagree
if (names(sort(summ$dtm_mean)[2L]) != "G/T")
  stop("G/T is not the 2nd-lowest mean dTm class")
summ$evidence_tier <- "measured"
summ$source_file <- "68_18048_2.xlsx :: Supplemental_Table2 + Supplemental_Table5"
summ$doi <- DOI_J
summ$outcome_definition <- paste(
  "Outcome = mean measured homozygote dTm (C) within one bimorphism, pooled over the",
  "38 Exp. 1 markers (dual-labelled probe Tm) and 18 Exp. 3 markers (NNNs-HRM), all",
  "soybean Aokimame x Fukuyutaka and Oryza Fukuyutaka x Toyoshirome. 'measured'.",
  "G/T has the 2nd-lowest mean dTm but is NOT 2nd on mean |ddG3bp| (A/G is), so the",
  "free-energy index does not order the substitution classes observed.")
readr::write_csv(summ, fs::path(out_dir, "hrm_jsbbs2018_dtm_by_bimorphism.csv"))

# =============================================================================
# A4  gundry2008: Tables 3 and 4 + the prose that restates Liew
# =============================================================================
doc <- xml2::read_xml(xml_g)
txt <- xml2::xml_text(doc)
if (!grepl("10.1093/nar/gkn204", txt, fixed = TRUE))
  stop("gundry2008 fulltext.xml does not carry the DOI")
tabs <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
cap <- vapply(tabs, function(t) norm(xml2::xml_text(
  xml2::xml_find_first(t, ".//*[local-name()='caption']"))), character(1))
pick <- function(rx) {
  i <- which(startsWith(cap, rx))
  if (length(i) != 1L) stop("Expected one table whose caption starts '", rx, "'")
  lapply(xml2::xml_find_all(tabs[[i]], ".//*[local-name()='tr']"), function(r)
    norm(xml2::xml_text(xml2::xml_find_all(r,
      ".//*[local-name()='td' or local-name()='th']"))))
}

# Table 3: cell order on a gene's first row is
# gene, bp, genotype, symmetric, predicted, before(avg, sd, range), after(avg, sd, range), n
t3 <- pick("Tm values of small amplicon homozygotes")
t3_rows <- Filter(function(v) length(v) %in% c(10L, 12L) && grepl("^(Yes|No)$", v[[4L]]), t3)
if (length(t3_rows) != 16L) stop("Table 3: expected 16 homozygote rows, got ", length(t3_rows))
# The gene name AND amplicon size print only on the T/T row of each pair; the
# A/A row leaves both blank. Carry the pair forward, then assert the pairing.
gene3 <- character(length(t3_rows)); last <- ""; last_bp <- ""
for (i in seq_along(t3_rows)) {
  cell <- trimws(t3_rows[[i]][[1L]])
  if (nzchar(cell) && !grepl("^[0-9]", cell)) {
    if (i > 1L && !grepl("^[ACGT]/[ACGT]$", t3_rows[[i]][[3L]]))
      stop("Table 3: a named gene row does not carry a genotype in column 3")
    last <- cell
    last_bp <- trimws(t3_rows[[i]][[2L]])
    if (!grepl("^[0-9]+$", last_bp))
      stop("Table 3: row ", i, " names a gene but has no amplicon size in column 2")
  }
  gene3[i] <- last
}
# 6 genes over 8 amplicons: CHK2 is assayed at three amplicon lengths (42/50/86 bp).
bp3 <- character(length(t3_rows)); last_bp <- ""
for (i in seq_along(t3_rows)) {
  cell <- trimws(t3_rows[[i]][[2L]])
  if (grepl("^[0-9]+$", cell)) last_bp <- cell
  bp3[i] <- last_bp
}
if (length(unique(paste(gene3, bp3))) != 8L)
  stop("Table 3: not 8 distinct gene+bp amplicons, got ",
       length(unique(paste(gene3, bp3))))

tg <- table(gene3)
if (any(gene3 == "") || length(tg) != 6L || any(tg < 2L) ||
    !identical(unname(as.integer(tg[sort(names(tg))])), c(2L, 6L, 2L, 2L, 2L, 2L)))
  stop("Table 3: gene fill-down is not 6 genes over 8 amplicons, got ",
       paste(names(tg), tg, collapse = "/"))
pick_num <- function(v, from) {
  n <- v[grepl("^[0-9]+\\.?[0-9]*$", v)]
  suppressWarnings(as.numeric(n[from:from]))
}
t3d <- tibble::tibble(
  amplicon = paste(gene3, bp3, "bp"),
  gene = gene3,
  product_bp = as.integer(bp3),
  genotype = vapply(t3_rows, function(v) v[[3L]], character(1)),
  symmetric_nn = vapply(t3_rows, function(v) v[[4L]] == "Yes", logical(1)),
  tm_predicted = vapply(t3_rows, function(v) as.numeric(v[[5L]]), numeric(1)),
  tm_before = vapply(t3_rows, function(v) as.numeric(v[[6L]]), numeric(1)),
  # 1 gene,2 bp,3 genotype,4 symmetric,5 predicted,6-8 before (avg,sd,range),
  # 9-11 after (avg,sd,range),12 n
  tm_after = vapply(t3_rows, function(v) as.numeric(v[[9L]]), numeric(1)),
  n = vapply(t3_rows, function(v) as.integer(v[[12L]]), integer(1)))
if (nrow(t3d) != 16L || anyDuplicated(paste(t3d$amplicon, t3d$genotype)))
  stop("Table 3 is not 8 amplicons x 2 homozygotes")
t3d$dtm_before <- ave(t3d$tm_before, t3d$amplicon, FUN = function(v) abs(diff(v)))
t3d$dtm_after <- ave(t3d$tm_after, t3d$amplicon, FUN = function(v) abs(diff(v)))
t3d$dtm_predicted <- ave(t3d$tm_predicted, t3d$amplicon, FUN = function(v) diff(v))
if (any(abs(t3d$dtm_after - t3d$dtm_before) > 1.0))
  stop("Table 3: calibration moved a homozygote pair by more than 1.0 C")
# the two NN-symmetric amplicons are exactly the two where prediction says 0.00
if (!setequal(unique(t3d$amplicon[t3d$symmetric_nn]), c("CPS1 51 bp", "OTC 68 bp")) ||
    !all(t3d$dtm_predicted[t3d$symmetric_nn] == 0))
  stop("Table 3: symmetric-NN amplicons are not the zero-prediction ones")
# prose names the four non-symmetric predictions
p_pred <- one_para(paras_of(doc),
  "nearest neighbor analysis predicted Tm differences of 0.19")
if (!grepl("0.19\u00b0C, 0.18\u00b0C, 0.31\u00b0C and 0.04\u00b0C", xml2::xml_text(
  xml2::xml_find_all(doc, "//*[local-name()='p']")[grepl("predicted Tm differences of 0.19",
    vapply(xml2::xml_find_all(doc, "//*[local-name()='p']"),
           function(p) norm(xml2::xml_text(p)), character(1)))][[1L]]), fixed = TRUE))
  stop("gundry prose: predicted 0.19/0.18/0.31/0.04 not found")
# Table 3 prints SIX non-symmetric predictions; the prose quotes only four of them
# (0.19, 0.18, 0.31, 0.04) and omits the two CHK2 amplicons (0.24, 0.11).
pp <- sort(unique(round(abs(t3d$dtm_predicted[!t3d$symmetric_nn]), 2)))
if (!isTRUE(all.equal(pp, c(0.04, 0.11, 0.18, 0.19, 0.24, 0.31))))
  stop("Table 3 predicted differences are not 0.04/0.11/0.18/0.19/0.24/0.31, got ",
       paste(pp, collapse = "/"))
if (!all(c(0.04, 0.18, 0.19, 0.31) %in% pp))
  stop("the four predictions quoted in the prose are not a subset of the table")
t3d$evidence_tier <- "measured"
t3d$source_file <- "fulltext.xml :: Table 3"
t3d$doi <- DOI_G
t3d$outcome_definition <- paste(
  "Outcome = average Tm (C) of one homozygote class under LightScanner small-amplicon",
  "melting, before and after internal-oligonucleotide calibration; dtm_* = |difference|",
  "between the two homozygotes of one amplicon. 'measured'. tm_predicted is the paper's",
  "nearest-neighbour value (a MODEL column inside a measured table). dtm_after is the",
  "PHYSICAL separation actually achieved, which is what an HRM genotyping call needs.")
if (any(t3d$dtm_after < 0 | t3d$dtm_after > 1.0))
  stop("Table 3 dtm_after outside 0-1 C")
readr::write_csv(t3d[, c("amplicon", "gene", "product_bp", "genotype", "symmetric_nn",
                         "tm_predicted", "tm_before", "tm_after", "dtm_before",
                         "dtm_after", "dtm_predicted", "n", "evidence_tier",
                         "source_file", "doi", "outcome_definition")],
                 fs::path(out_dir, "hrm_gundry2008_amplicon_tm.csv"))

# Table 4: gene, bp, base change, NN symmetric, before %, after %, n
t4_rows <- Filter(function(v) length(v) == 7L && grepl("^[0-9.]+$", v[[5L]]),
                  pick("Homozygous genotyping accuracy"))
if (length(t4_rows) != 6L) stop("Table 4: expected 6 loci, got ", length(t4_rows))
acc <- tibble::tibble(
  gene = vapply(t4_rows, function(v) v[[1L]], character(1)),
  product_bp = vapply(t4_rows, function(v) as.integer(v[[2L]]), integer(1)),
  base_change = vapply(t4_rows, function(v) v[[3L]], character(1)),
  symmetric_nn = vapply(t4_rows, function(v) v[[4L]] == "Yes", logical(1)),
  accuracy_before_pct = vapply(t4_rows, function(v) as.numeric(v[[5L]]), numeric(1)),
  accuracy_after_pct = vapply(t4_rows, function(v) as.numeric(v[[6L]]), numeric(1)),
  n = vapply(t4_rows, function(v) as.integer(v[[7L]]), integer(1)))
if (sum(acc$n) != 1147L) stop("Table 4 n does not sum to 1147, got ", sum(acc$n))
if (sum(acc$accuracy_after_pct >= 99) != 5L)
  stop("Table 4: expected 5 of 6 loci at >= 99% after calibration")
if (which.min(acc$accuracy_after_pct - acc$accuracy_before_pct) != which(acc$gene == "BX647987"))
  stop("Table 4: the only locus that gets worse is not BX647987")
if (acc$accuracy_after_pct[acc$gene == "BX647987"] >= acc$accuracy_before_pct[acc$gene == "BX647987"])
  stop("Table 4: BX647987 accuracy did not regress")
# MSH2 is 49 bp in Table 3 but 48 bp in Table 4 -- a source discrepancy, kept
if (!(t3d$product_bp[t3d$gene == "MSH2"][1L] == 49L && acc$product_bp[acc$gene == "MSH2"] == 48L))
  stop("MSH2 amplicon length discrepancy (49 vs 48) is not reproduced")
acc$evidence_tier <- "measured"
acc$source_file <- "fulltext.xml :: Table 4"
acc$doi <- DOI_G
acc$outcome_definition <- paste(
  "Outcome = percent of samples whose measured Tm fell closer to the correct than to the",
  "wrong homozygote mean, with genotypes pre-called by probe/sequencing. 'measured',",
  "read from Table 4. Locus level, 42-296 samples each; not a per-SNP-class rate. Note",
  "the printed aggregate (91.7 -> 99.7% over 1105 amplicons) is NOT re-derived here:",
  "the sample x amplicon x target structure behind 1105 is not published in the file.")
readr::write_csv(acc[, c("gene", "product_bp", "base_change", "symmetric_nn",
                         "accuracy_before_pct", "accuracy_after_pct", "n",
                         "evidence_tier", "source_file", "doi", "outcome_definition")],
                 fs::path(out_dir, "hrm_gundry2008_homozygote_accuracy.csv"))

# =============================================================================
# A5  liew2004 restated bands -- SECONDARY CITATION ONLY
# =============================================================================
gp <- paras_of(doc)
# ONE paragraph carries every Liew-derived number. The tildes are U+223C in the
# XML, so the text is normalised before matching.
s_all <- one_para(gp, "resolves homozygous alleles in at least 84%")
s_all <- gsub("\u223c", "~", s_all, fixed = TRUE)
for (need in c("at least 84%", "0.8-1.4\u00b0C (10)", "~16% of single base variants",
               "0.4\u00b0C or less (10)", "4% of human single base variants")) {
  if (!grepl(need, s_all, fixed = TRUE))
    stop("gundry prose: Liew-derived phrase not found: ", need)
}
# The two references the prose points at must be the Liew papers.
refs <- vapply(xml2::xml_find_all(doc, "//*[local-name()='ref']"),
               function(r) norm(xml2::xml_text(r)), character(1))
lr <- refs[grepl("Liew", refs, fixed = TRUE)]
if (length(lr) < 2L) stop("gundry does not cite Liew twice (10) and (11)")
if (!any(grepl("10.1373/clinchem.2004.032136", lr, fixed = TRUE)))
  stop("gundry ref (10) does not carry the Liew 2004 DOI 10.1373/clinchem.2004.032136")
# Gundry cites Liew FOUR times; only ref (10) is the 2004 Clin Chem paper. Ref (11)
# is Palais 2005 (Anal Biochem), NOT a Liew NN-symmetry paper -- do not claim one.
if (length(lr) != 4L) stop("gundry should cite Liew 4 times, got ", length(lr))
if (sum(grepl("10.1373/clinchem.2004.032136", lr, fixed = TRUE)) != 1L)
  stop("gundry should carry the Liew 2004 DOI exactly once")
if (!any(grepl("10.1016/j.ab.2005.08.010", refs, fixed = TRUE)))
  stop("gundry ref (11) is not Palais 2005 (Anal Biochem)")

liew <- tibble::tibble(
  quantity = c("dtm_at_to_gc_lower_c", "dtm_at_to_gc_upper_c",
               "dtm_other_substitutions_upper_c",
               "nn_symmetric_unresolvable_pct_liew2004",
               "nn_symmetric_bp_neutral_pct_gundry",
               "homozygote_resolved_pct_gundry"),
  value = c(0.8, 1.4, 0.4, 4, 16, 84),
  unit = c("C", "C", "C", "pct", "pct", "pct"))
liew$evidence_tier <- "secondary_cited"
liew$source_file <- "fulltext.xml :: Gundry 2008, intro paragraph + ref (10)"
liew$doi <- DOI_G
liew$cited_by <- "Gundry 2008 prose, ref (10) = Liew 2004 DOI 10.1373/clinchem.2004.032136"
liew$outcome_definition <- paste(
  "NOT MEASURED HERE. These are the numbers Gundry 2008 attributes to Liew 2004, restated.",
  "Liew 2004 is closed access (no PMCID); it was not read. value = the printed bound:",
  "0.8 and 1.4 are the two ends of the A::T<->G::C homozygote Tm-difference range,",
  "0.4 is the upper bound for the remaining substitutions, 4 is Liew's symmetric",
  "unresolvable fraction, 16 is Gundry's own base-pair-neutral fraction, 84 is",
  "100 - Gundry's stated residual resolvable share. Any package text quoting these",
  "must label them 'cited, not verified'.")
if (nrow(liew) != 6L) stop("liew restatement should be 6 rows")
if (!all(liew$value[liew$quantity == "dtm_at_to_gc_lower_c"] <
         liew$value[liew$quantity == "dtm_at_to_gc_upper_c"]))
  stop("liew: 0.8 not below 1.4")
if (liew$value[liew$quantity == "dtm_other_substitutions_upper_c"] >=
    liew$value[liew$quantity == "dtm_at_to_gc_lower_c"])
  stop("liew: the 'other substitutions' bound (0.4) is not below the A::T<->G::C floor (0.8)")
if (abs(liew$value[liew$quantity == "homozygote_resolved_pct_gundry"] +
        liew$value[liew$quantity == "nn_symmetric_bp_neutral_pct_gundry"] - 100) > 0)
  stop("gundry: 84 and 16 do not sum to 100")
readr::write_csv(liew[, c("quantity", "value", "unit", "cited_by", "evidence_tier",
                          "source_file", "doi", "outcome_definition")],
                 fs::path(out_dir, "hrm_liew2004_restated_bands.csv"))

message("hrmplant: wrote ", length(list.files(out_dir, pattern = "^hrm_(jsbbs2018|gundry2008|liew2004)")),
        " CSVs to ", out_dir, "; liew2004 = secondary_cited, NOT verified")
