# =============================================================================
#  data-raw/vcf2caps2021.R  --  VCF2CAPS wet-lab verification, re-derived
# =============================================================================
#  Source (identifier resolution is PROVENANCE, not evidence):
#
#    Wesołowski M, Domnicz A, Augustynowicz J, et al.  VCF2CAPS - A
#    high-throughput CAPS marker design from VCF files and its test-use on a
#    genotyping-by-sequencing (GBS) dataset.  PLoS Computational Biology 2021;
#    17(5):e1008980.  DOI 10.1371/journal.pcbi.1008980  (PMID 34014924,
#    PMCID PMC8186816)
#
#  Fetched: .tmp/lit/vcf2caps2021/fulltext.xml (Europe PMC fullTextXML for
#  PMC8186816, 150142 bytes, contains <article>).  The DOI is confirmed from
#  the Crossref API and from the EPMC record; NOTE that a prior dispatch
#  carried "10.1371/journal.pcbi.1009001", which Crossref resolves to an
#  unrelated Group B Streptococcus paper - that DOI is WRONG and is not used.
#
#  ===========================================================================
#  THE REGISTRY PREMISE THIS CORRECTS
# ===========================================================================
#  `inst/extdata/scientific_backlog.csv` row vcf2caps2021 says the paper
#  "reports no wet-lab outcomes".  THAT IS FALSE.  It has a Results
#  subsection, "Experimental verification of the identified CAPS markers",
#  reporting 10 single-cut markers x 30 plants each (Table 1, Fig 8), and the
#  outcome of every marker is stated IN PROSE.
#
#  Exactly what is TABULATED vs ASSERTED, since this is the one CAPS number a
#  prior session flagged as unverified:
#    TABULATED - Table 1 only: marker name, primers, enzyme, and the 10/10/10
#      per-marker genotype design (30 plants per marker).  NOT a call table.
#    ASSERTED (prose) - the per-marker agreement between predicted and observed
#      genotype class; for 3_363 the count of five reclassified plants; for
#      3_300 and 6_022 merely "single supposed heterozygotes" (no count).
#    FIG 8 is a gel IMAGE and carries no table-wrap; this build has no
#      supplementary material (no <supplementary-material> element).
#  WHAT THIS FILE LANDS vs WHAT IT DOES NOT
#
#  * LANDS (counts the prose states unambiguously, verified by matching the
#    paper's own wording for each marker):
#      - n_markers_verified = 10, n_plants_each = 30 (both stated);
#      - the three observed genotype classes were 10/10/10 BY DESIGN, so the
#        design total must be 10 x 30 = 300 plant x marker observations;
#      - the split of the 10 markers into the prose's three wording groups
#        (exact match 2; "very close"/"similar" 4; "remaining two" 2) and the
#        exceptional group (2), where every predicted heterozygote read as a
#        non-cleavable homozygote.
#    Each group count is not a guess: the script locates the sentence that
#    names each marker and checks it carries the expected phrase.
#
#  * LANDS, as a count, ONE class-level correction that the prose states
#    explicitly: for marker 3_363, "the remaining five supposed non-cleavable
#    homozygotes were in fact heterozygotes" - i.e. 5 of the 10 plants predicted
#    non-cleavable were reclassified, and no other plant on that marker was.
#    It is asserted as 0 < 5 <= 10 (the paper gives no denominator table).
#
#  * DOES NOT LAND any per-marker "x/30 correct", any overall concordance, or
#    any error rate.  There is no per-plant call table; Fig 8 is a gel image
#    (no <table-wrap> exists for it) and the paper has no supplementary
#    material (<supplementary-material> absent from the XML).  Words like
#    "very close" and "looked like" cannot be turned into a number.  For 3_300
#    and 6_022 the prose says "single supposed heterozygotes" (plural noun,
#    unspecified count) so even the number of misread plants on those markers
#    is not stated.
#
#  Run:  Rscript data-raw/vcf2caps2021.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

xml_path <- fs::path(root, ".tmp", "lit", "vcf2caps2021", "fulltext.xml")
if (!file.exists(xml_path)) stop("Missing ", xml_path)
doc <- xml2::read_xml(xml_path)
out_dir <- fs::path(root, "inst", "extdata")
DOI <- "10.1371/journal.pcbi.1008980"
PMCID <- "PMC8186816"

# --- the paper's identifier must be the one we think it is -------------------
# Crossref holds the title/DOI pair; a wrong DOI cannot silently pass.
front <- xml2::xml_text(xml2::xml_find_first(doc, "//*[local-name()='article-title']"))
if (!grepl("VCF2CAPS", front, fixed = TRUE)) {
  stop("Fetched XML is not the VCF2CAPS paper: article-title = ", front)
}
xid <- xml2::xml_text(xml2::xml_find_all(doc, "//*[local-name()='article-id']"))
if (!DOI %in% xid) stop("XML does not carry DOI ", DOI, "; found: ", paste(xid, collapse = ", "))

## --- the Results subsection, as flat text --------------------------------
sec <- xml2::xml_find_first(doc,
  "//*[local-name()='sec'][*[local-name()='title' and contains(., 'Experimental verification of the identified')]]")
if (inherits(sec, "xml_missing")) stop("Results subsection 'Experimental verification of the identified CAPS markers' not found")
txt <- paste(xml2::xml_text(xml2::xml_find_all(sec, ".//text()")), collapse = " ")
txt <- gsub("[[:space:]]+", " ", txt)

## --- Table 1: the 10 markers taken forward to the bench -------------------
tw <- xml2::xml_find_all(doc, "//*[local-name()='table-wrap']")
cap <- vapply(tw, function(t) xml2::xml_text(xml2::xml_find_first(
  t, ".//*[local-name()='caption']")), character(1))
tbl1 <- tw[grepl("selected for experimental verification", cap)]
if (length(tbl1) != 1L) stop("Could not locate the verification Table 1")
trs <- xml2::xml_find_all(tbl1[[1L]], ".//*[local-name()='tr']")
cells <- lapply(trs, function(r)
  xml2::xml_text(xml2::xml_find_all(r, ".//*[local-name()='td' or local-name()='th']")))
data_rows <- Filter(function(v) length(v) == 6L && grepl("^[0-9]", v[1L]), cells)
if (length(data_rows) != 10L) stop("Expected 10 marker rows in verification Table 1, got ", length(data_rows))
table1 <- tibble::tibble(
  marker_original = vapply(data_rows, `[[`, character(1), 1L),
  marker_abbrev   = gsub("-", "_", gsub("\u2013", "_",
                    vapply(data_rows, `[[`, character(1), 2L))),
  chromosome      = vapply(data_rows, `[[`, character(1), 3L),
  forward_primer  = vapply(data_rows, `[[`, character(1), 4L),
  reverse_primer  = vapply(data_rows, `[[`, character(1), 5L),
  enzyme          = vapply(data_rows, `[[`, character(1), 6L)
)
if (anyDuplicated(table1$marker_abbrev)) stop("Duplicate marker ids in Table 1")

## --- the only two numbers the prose states outright -----------------------
# Every pattern below is the paper's own phrase, matched literally.
must <- c(
  n_markers         = "10 single-cut markers",
  n_plants          = "30 plants from population 398",
  design_classes    = "10 non-cleavable homozygotes, 10 heterozygotes and 10 cleavable",
  pcr_accordant     = "size of the obtained PCR product was accordant",
  fragment_accordant = "sizes of the derived restriction fragments",
  reclass_five      = "remaining five supposed non-cleavable homozygotes were in fact heterozygotes"
)
for (nm in names(must)) {
  if (!grepl(must[[nm]], txt, fixed = TRUE)) {
    stop("Prose no longer contains the phrase for '", nm, "': ", must[[nm]])
  }
}
n_markers_verified <- as.integer(regmatches(txt, regexpr("[0-9]+(?= single-cut markers)", txt, perl = TRUE)))
n_plants_each <- as.integer(regmatches(txt, regexpr("[0-9]+(?= plants from population)", txt, perl = TRUE)))
n_per_class <- as.integer(regmatches(must[["design_classes"]],
  gregexpr("[0-9]+", must[["design_classes"]], perl = TRUE))[[1L]])
if (n_markers_verified != 10L) stop("n_markers_verified parsed as ", n_markers_verified)
if (n_plants_each != 30L) stop("n_plants_each parsed as ", n_plants_each)
if (!identical(n_per_class, c(10L, 10L, 10L))) stop("design class counts are not 10/10/10")
# paper's own internal sum: the design split must add up to the tested n
if (sum(n_per_class) != n_plants_each) {
  stop("Design classes sum to ", sum(n_per_class), ", but the paper tests ", n_plants_each, " plants per marker")
}
# and Table 1 must carry exactly the markers the prose verified
if (nrow(table1) != n_markers_verified) {
  stop("Table 1 has ", nrow(table1), " markers but the prose says ", n_markers_verified)
}

## --- per-marker outcome groups, keyed by the paper's wording --------------
# Each entry carries a REGEX that must fire in the sentence naming that
# marker.  If the paper is revised and a marker moves between wording
# groups, the generator stops rather than silently re-labelling it.
sentences <- strsplit(txt, "(?<=[.!?]) (?=[A-Z])", perl = TRUE)[[1L]]
sentence_for <- function(id, window = 1L) {
  idx <- which(vapply(sentences, function(s) grepl(id, s, fixed = TRUE), logical(1)))
  idx <- idx[!grepl("DNA size marker", sentences[idx], fixed = TRUE)]
  if (length(idx) != 1L) stop("Expected exactly one sentence naming marker ", id, ", found ", length(idx))
  # `window = 2` reads the phrase from the marker's sentence plus the next one,
  # because the paper splits "...homozygotes. This marker showed that the
  # remaining five..." across a sentence boundary.
  paste(sentences[idx:min(idx + window - 1L, length(sentences))], collapse = " ")
}
groups <- list(
  exact_match = list(ids = c("1_028", "3_579"),
    re = "distribution of their restriction profiles exactly matched genotypic predictions"),
  weaker_cleavable_allele = list(ids = "3_351",
    re = "the band of the cleavable allele was weaker than that of the non-cleavable allele"),
  close_single_het = list(ids = c("3_714", "3_300", "6_022"),
    re = "the only discrepancy was shown for single supposed heterozygotes which marker-wise turned out to be homozygotes"),
  five_noncleavable_reclassified = list(ids = "3_363", window = 2L,
    re = "remaining five supposed non-cleavable homozygotes were in fact heterozygotes"),
  cleavable_hom_read_as_het_banding = list(ids = "1_017", window = 2L,
    re = "supposed cleavable homozygotes were similar to those of heterozygotes"),
  het_read_as_noncleavable = list(ids = c("2_525", "s_256"),
    re = "all supposed heterozygotes looked like the non-cleavable homozygotes")
)
for (g in names(groups)) {
  for (id in groups[[g]]$ids) {
    w <- if (is.null(groups[[g]]$window)) 1L else groups[[g]]$window
    s <- sentence_for(id, w)
    if (!grepl(groups[[g]]$re, s, fixed = TRUE)) {
      stop("Marker ", id, " no longer carries the '", g, "' wording; its sentence is: ", s)
    }
  }
}
group_of <- stats::setNames(
  rep(names(groups), vapply(groups, function(g) length(g$ids), integer(1))),
  unlist(lapply(groups, `[[`, "ids"), use.names = FALSE))
if (!setequal(names(group_of), table1$marker_abbrev)) {
  stop("Wording groups cover {", paste(names(group_of), collapse = ", "),
       "} but Table 1 is {", paste(table1$marker_abbrev, collapse = ", "), "}")
}
group_totals <- table(group_of)
EXPECTED_GROUP_N <- c(exact_match = 2L, weaker_cleavable_allele = 1L,
                      close_single_het = 3L, five_noncleavable_reclassified = 1L,
                      cleavable_hom_read_as_het_banding = 1L, het_read_as_noncleavable = 2L)
if (!identical(as.integer(group_totals[names(EXPECTED_GROUP_N)]), unname(EXPECTED_GROUP_N))) {
  stop("Wording-group sizes changed: ", paste(names(group_totals), group_totals, sep = "=", collapse = ", "))
}
if (sum(EXPECTED_GROUP_N) != n_markers_verified) {
  stop("Wording groups sum to ", sum(EXPECTED_GROUP_N), ", not ", n_markers_verified, " verified markers")
}
if (!any(sapply(groups, function(g) "3_300" %in% g$ids)) ||
    !any(sapply(groups, function(g) "6_022" %in% g$ids))) stop("3_300 / 6_022 missing from the groups")

## --- the one per-marker number the prose states ---------------------------
reclass <- table1 |>
  dplyr::mutate(
    predicted_noncleavable = dplyr::if_else(marker_abbrev == "3_363",
                                            as.integer(n_per_class[1L]), NA_integer_),
    observed_reclassified_as_het = dplyr::if_else(marker_abbrev == "3_363", 5L, NA_integer_),
    observed_supports_prediction = predicted_noncleavable - observed_reclassified_as_het
  )
# "five of ten" must sit inside the class the paper designed: 0 < 5 <= 10.
five <- reclass$observed_reclassified_as_het[reclass$marker_abbrev == "3_363"]
if (!identical(five, 5L)) stop("3_363 reclassified count is ", five, ", prose says five")
if (!(five > 0L && five <= n_per_class[1L])) {
  stop("Reclassified count ", five, " is outside 1..", n_per_class[1L])
}
# "no other plant on that marker was" - the prose names no other mismatch for
# 3_363; if it ever reports a second one, this script must be revisited.
if (grepl("3_363", txt, fixed = TRUE) &&
    grepl("other discrepancy", txt, fixed = TRUE)) stop("3_363 gained a second reported discrepancy")

OUTCOME_DEF <- paste(
  "VCF2CAPS wet-lab outcome = the genotype class read off a gel for a plant",
  "(non-cleavable homo / heterozygote / cleavable homo) compared with the class",
  "VCF2CAPS predicted. Reported per marker in prose; class agreement, not a",
  "fragment-size measurement. Counts exist only where the prose states them."
)
EVIDENCE <- "fulltext prose (PMC8186816); no per-plant call table, no supplement"

reclass <- reclass |>
  dplyr::mutate(
    outcome_group = unname(group_of[marker_abbrev]),
    outcome_definition = OUTCOME_DEF,
    n_plants_per_marker = as.integer(n_plants_each),
    n_plants_per_class_design = paste(n_per_class, collapse = "/"),
    n_plant_x_marker_observations = as.integer(n_markers_verified * n_plants_each),
    pdf_accordant_amplicon_and_fragments = TRUE,
    evidence_tier = "strong",
    evidence_basis = EVIDENCE,
    source_file = "fulltext.xml",
    fetch_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8186816/fullTextXML",
    pmcid = PMCID,
    doi = DOI,
    reported_verdict_verbatim = vapply(marker_abbrev, sentence_for, character(1), USE.NAMES = FALSE)
  ) |>
  dplyr::arrange(marker_abbrev)

# every reported reclassification must be a subset of its designed class
bad <- dplyr::filter(reclass, !is.na(observed_reclassified_as_het),
                     observed_reclassified_as_het > n_per_class[1L] | observed_supports_prediction < 0L)
if (nrow(bad)) stop("Per-marker counts are internally impossible: ", bad$marker_abbrev)

readr::write_csv(reclass, fs::path(out_dir, "vcf2caps2021_verification.csv"))
# The headline is a two-column key/value table with a real header.  Writing
# bare "key,value" lines (the previous form) makes any reader take the first
# row as the header and name the columns "n_markers_verified" and "X10".
headline <- data.frame(
  quantity = c("n_markers_verified", "n_plants_each", "n_plant_x_marker_observations",
               names(EXPECTED_GROUP_N), "marker_3_363_reclassified_as_het"),
  value = c(n_markers_verified, n_plants_each, n_markers_verified * n_plants_each,
            unname(EXPECTED_GROUP_N), five),
  stringsAsFactors = FALSE
)
# the group counts are a partition of the 10 markers, and the design counts of
# the 30 plants; both are asserted here so the file cannot drift from the table
if (sum(headline$value[headline$quantity %in% names(EXPECTED_GROUP_N)]) != n_markers_verified) {
  stop("Wording groups do not sum to the number of verified markers")
}
if (unname(n_per_class[1L]) * length(n_per_class) != n_plants_each) {
  stop("The designed class counts do not sum to the plants per marker")
}
readr::write_csv(headline, fs::path(out_dir, "vcf2caps2021_headline.csv"))

message("vcf2caps2021: 10 markers x 30 plants asserted from prose; wording groups ",
        paste(names(EXPECTED_GROUP_N), EXPECTED_GROUP_N, sep = "=", collapse = " "),
        "; 3_363 = 5 of 10 reclassified (assertion passed)")
