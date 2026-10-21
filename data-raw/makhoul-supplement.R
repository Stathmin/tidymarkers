# =============================================================================
#  data-raw/makhoul-supplement.R  --  the first supplement we actually opened
# =============================================================================
#  Makhoul et al. 2020 (TAG, 10.1007/s00122-020-03608-x, PMC7360542) is the only
#  source in the registry whose supplementary workbook we have downloaded and
#  parsed.  It is also the only KASP conversion study in the registry that
#  reports, per SNP, WHY an assay failed AND what the design strategy was.
#
#  Used by Z8 (gate-pass -> outcome) and Y3.2 (homeolog-specific design).
#
#  Fetch (a ZIP, not a tar, despite the endpoint name):
#    curl -sS -o <tmp>/makhoul.zip \
#      'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7360542/supplementaryFiles'
#  The per-SNP workbook is MOESM1_ESM.xlsx inside it.
#
#  Two things this script refuses to do:
#
#  1. It does not read cells that look 'merged'.  Table S4 leaves the
#     strategy columns (M, N, O) blank on continuation rows of a multi-hit
#     SNP block; a naive read counts 3 SNPs instead of 15.  We forward-fill
#     within a SNP block, then ASSERT the totals equal the numbers Table 3 of
#     the paper reports.  If the paper and the supplement ever disagree, this
#     script fails loudly rather than writing a silent average.
#
#  2. It does not translate the outcome vocabulary.  The authors distinguish
#     "Correct cluster", "Wrong cluster", "No cluster" and "Locus-specific
#     design not possible".  Those are readout states, not a binary, and the
#     CSV keeps the strings verbatim (plus a coarse ordinal for grouping).
#
#  Run:  Rscript data-raw/makhoul-supplement.R [path/to/MOESM1_ESM.xlsx]
#        (no installed xlsx reader is available in this toolchain, so the
#         workbook is parsed with the sheet XML directly)
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
xlsx <- if (length(args)) args[[1L]] else Sys.getenv("MAKHOUL_XLSX", "")
if (!nzchar(xlsx) || !file.exists(xlsx)) {
  stop("Pass the path to MOESM1_ESM.xlsx (see the header) or set MAKHOUL_XLSX.")
}

out_dir <- fs::path(root, "inst", "extdata")
DOI <- "10.1007/s00122-020-03608-x"

## --- minimal xlsx reader ------------------------------------------------------
# Shared strings + one worksheet, cell references resolved to A1 columns.  This
# is the whole of what we need and avoids adding a dependency for one file.
NS <- c(
  d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
  r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
)

xlsx_sheets <- function(path) {
  exdir <- tempfile("xlsx")
  unz <- utils::unzip(path, exdir = exdir)
  on.exit(unlink(exdir, recursive = TRUE), add = TRUE)
  sst_path <- file.path(exdir, "xl", "sharedStrings.xml")
  sst <- character()
  if (file.exists(sst_path)) {
    doc <- xml2::read_xml(sst_path)
    # an <si> can hold several <t> runs; concatenating them is what the text is
    sst <- vapply(xml2::xml_find_all(doc, ".//d1:si", NS), function(si) {
      paste(xml2::xml_text(xml2::xml_find_all(si, ".//d1:t", NS)), collapse = "")
    }, character(1))
  }
  wb <- xml2::read_xml(file.path(exdir, "xl", "workbook.xml"))
  rels <- xml2::read_xml(file.path(exdir, "xl", "_rels", "workbook.xml.rels"))
  # strip namespaces here: the relationships part has none, and `r:id` on a
  # stripped <sheet> reads as the plain attribute "id"
  wb <- xml2::xml_ns_strip(wb)
  rels <- xml2::xml_ns_strip(rels)
  rel_nodes <- xml2::xml_find_all(rels, ".//Relationship")
  rel <- stats::setNames(
    xml2::xml_attr(rel_nodes, "Target"),
    xml2::xml_attr(rel_nodes, "Id")
  )
  sheets <- xml2::xml_find_all(wb, ".//sheet")
  out <- list()
  for (sh in sheets) {
    nm <- xml2::xml_attr(sh, "name")
    rid <- xml2::xml_attr(sh, "id")
    path_sheet <- file.path(exdir, "xl", sub("^/?", "", rel[[rid]]))
    doc <- xml2::read_xml(path_sheet)
    rows <- xml2::xml_find_all(doc, ".//d1:row", NS)
    out[[nm]] <- lapply(rows, function(row) {
      cells <- xml2::xml_find_all(row, "./d1:c", NS)
      ref <- xml2::xml_attr(cells, "r")
      col <- sub("[0-9]+$", "", ref)
      type <- xml2::xml_attr(cells, "t")
      raw <- xml2::xml_child(cells, "d1:v", NS)
      v <- xml2::xml_text(raw)
      is_str <- !is.na(type) & type == "s"
      v[is_str] <- sst[as.integer(v[is_str]) + 1L]
      stats::setNames(as.list(v), col)
    })
  }
  out
}

sheets <- xlsx_sheets(xlsx)
stopifnot("S4" %in% names(sheets))   # Tab S4 = BLASTN + the three design strategies

## --- Table S4: one row per SNP, strategy outcomes resolved per block ---------
# NOT one row per physical spreadsheet row.  The authors left the three strategy
# columns empty on the row that names each SNP and filled them on the row for
# the 5B1 BLAST hit (which comes after the 5A/5D hits, sometimes two rows down).
# Reading the first row of each block loses 3 of the 15 SNPs and returns 3/6/8
# against the paper's 3/8/11; reading every row double-counts.  We therefore
# group by block, take the non-empty strategy values, and require exactly one.
get <- function(cells, k) if (is.null(cells[[k]])) NA_character_ else cells[[k]]

s4 <- sheets[["S4"]]
hap <- NA_character_
snp <- NULL
recs <- list()
for (row in s4) {
  if (!is.null(row[["A"]]) && grepl("^Hap-5B", row[["A"]])) hap <- row[["A"]]
  numbered <- !is.null(row[["B"]]) && grepl("^[0-9]+$", row[["B"]])
  named <- numbered && !is.null(row[["C"]]) && nzchar(row[["C"]])
  if (named) {
    snp <- list(
      hap = hap, snp_no = as.integer(row[["B"]]), snp = get(row, "C"),
      alleles_tg = get(row, "D"), n_hits = suppressWarnings(as.integer(row[["E"]])),
      best_hit = get(row, "F"),
      s1_flanks_only = character(), s2_single_ref = character(),
      s3_multi_sanger = character()
    )
    recs[[length(recs) + 1L]] <- snp
  }
  if (length(recs) && any(c("M", "N", "O") %in% names(row))) {
    last <- length(recs)
    for (k in c("M", "N", "O")) {
      v <- get(row, k)
      if (!is.na(v) && nzchar(v)) recs[[last]][[c(M = "s1_flanks_only", N = "s2_single_ref",
                                                 O = "s3_multi_sanger")[[k]]]] <-
        c(recs[[last]][[c(M = "s1_flanks_only", N = "s2_single_ref",
                          O = "s3_multi_sanger")[[k]]]], v)
    }
  }
}
# a block that reports two different strategy outcomes would be a supplement we
# misunderstand -- refuse rather than pick one.  A literal "-" is the authors'
# marker for "this approach was not attempted" (e.g. GENE-2890_482 lists it for
# its 5A1 and 5D1 hits); empty merged cells carry no value at all and are never
# read as an outcome.
distinct_one <- function(x) {
  u <- setdiff(unique(x), "-")
  if (length(u) > 1L) stop("conflicting outcomes within one SNP block: ", paste(u, collapse = " / "))
  if (!length(u)) NA_character_ else u
}
design <- purrr::map_dfr(recs, function(r) {
  r[c("s1_flanks_only", "s2_single_ref", "s3_multi_sanger")] <-
    lapply(r[c("s1_flanks_only", "s2_single_ref", "s3_multi_sanger")], distinct_one)
  tibble::as_tibble(r)
}) |>
  # the table's footnote rows sit under numbers too and carry a long note in the
  # SNP-name column; the real SNP ids are chip probe names, never sentences
  dplyr::filter(!grepl(" ", snp))

## --- the assertion that makes this a *verified* extract -----------------------
# Table 3 of the paper reports 3/15, 8/15 and 11/15.  We compute them from the
# supplement instead of typing them, and fail if they disagree.
count_ok <- function(x) sum(x == "Correct cluster", na.rm = TRUE)
observed <- c(
  `approach 1 (101 bp flanks only)` = count_ok(design$s1_flanks_only),
  `approach 2 (single reference genome)` = count_ok(design$s2_single_ref),
  `approach 3 (multiple assemblies + Sanger)` = count_ok(design$s3_multi_sanger)
)
expected <- c(3L, 8L, 11L)
if (!identical(unname(observed), expected) || nrow(design) != 15L) {
  stop(
    "Supplement and paper disagree: parsed ", nrow(design), " SNPs with correct-cluster ",
    "counts ", paste(observed, collapse = "/"), " against the paper's 3/8/11."
  )
}

readr::write_csv(
  dplyr::mutate(
    design,
    doi = DOI,
    # the readout state, verbatim from the workbook, is the outcome; the coarse
    # ordinal exists only for grouping and is NOT a probability
    outcome_ordinal = dplyr::case_when(
      s3_multi_sanger == "Correct cluster" ~ "usable",
      grepl("not possible", s3_multi_sanger) ~ "doomed",
      TRUE ~ "not usable"
    ),
    outcome_def = paste(
      "KASP cluster inspection: 'Correct cluster' = clusters separated in the",
      "expected positions; 'Wrong cluster' = amplified but clusters not assignable",
      "to the haplotypes; 'No cluster' = no amplification; 'Locus-specific design",
      "not possible' = the author judged no 5B-specific assay constructible"
    )
  ) |>
    dplyr::select(hap, snp_no, snp, alleles_tg, n_hits, best_hit,
                  s1_flanks_only, s2_single_ref, s3_multi_sanger,
                  outcome_ordinal, outcome_def, doi),
  fs::path(out_dir, "makhoul_design_strategies.csv")
)

## --- Table S1: real KASP primers with the allele in the AS 3' base -----------
# This is the Y3.2 artifact.  For each assay the two allele-specific primers
# differ in their 3'-terminal base (which carries the SNP) while the common
# primer is the one that can carry the subgenome difference.  We record both so
# the artifact can show the two axes separately instead of conflating them.
s1 <- sheets[["S1"]]
hap <- NA_character_
prim <- list()
for (row in s1) {
  if (!is.null(row[["A"]]) && grepl("^Hap-", row[["A"]])) hap <- row[["A"]]
  if (is.null(row[["D"]]) || is.na(row[["D"]])) next
  if (!grepl("^Hap[AB][0-9]", row[["D"]])) next
  prim[[length(prim) + 1L]] <- list(
    hap = hap, snp = get(row, "B"), assay = row[["D"]], snp_kasp = get(row, "E"),
    fam = get(row, "F"), fam_hap = get(row, "G"),
    hex = get(row, "H"), hex_hap = get(row, "I"),
    common = get(row, "J"), product_bp = get(row, "K")
  )
}
primers <- purrr::map_dfr(prim, tibble::as_tibble) |>
  dplyr::mutate(
    # the allele-bearing base: the last base of each allele-specific primer
    fam_3p = substr(fam, nchar(fam), nchar(fam)),
    hex_3p = substr(hex, nchar(hex), nchar(hex)),
    common_3p = substr(common, nchar(common), nchar(common)),
    # the allele is carried by the AS 3' base -- never by the common primer
    allele_in_as_3p = fam_3p != hex_3p,
    doi = DOI
  )
stopifnot(all(primers$allele_in_as_3p))

readr::write_csv(
  primers,
  fs::path(out_dir, "makhoul_kasp_primers.csv"),
  na = ""
)

message(
  "makhoul: ", nrow(design), " SNPs, correct-cluster ", paste(observed, collapse = "/"),
  " (paper 3/8/11); ", nrow(primers), " KASP primer sets"
)
