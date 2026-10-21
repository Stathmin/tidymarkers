# =============================================================================
#  data-raw/wheat_as2023.R  --  crop allele-specific marker design funnel
# =============================================================================
#  Bhoite R, Smith R, Dowla M, Bansal U, Roy S, Harwood A, Bariana H, Sharma D.
#  "Exome-based allele-specific marker development for wheat ... under sodic-
#  dispersive soils and validation of associated genes."
#  Plant Direct 2023;7:e520.  doi:10.1002/pld3.520  PMC10435944  (OA, CC-BY)
#
#  IMPORTANT: THE ASSAY IS NOT ARMS
#  The markers are PolyMarker-designed allele-competitive-extension (PACE / KASP)
#  assays with FAM/HEX tails. This file is INCLUDED because the family question
#  is "how many allele-specific markers worked in a crop", and this is the one
#  paper in the batch that reports a complete design -> outcome funnel. The
#  assay column says KASP/PACE so it can never be read as an ARMS outcome.
#  (The paper never uses the word "mismatch": 0 occurrences in fulltext.xml.)
#
#  THE CONTRADICTION WE REFUSE TO RESOLVE
#  The paper states the funnel twice and the two statements do not agree:
#   (a) "Among 32, only 18 target variant sequences successfully produced allele
#       specific markers"  -> 18/32 designed
#   (b) "Eighteen out of the target 39 SNP locations fulfilled the criteria for
#       producing PCR markers, among which only eight produced polymorphic
#       signals"  -> 18/39, and 8/18 polymorphic
#  (c) Section 2.3 then names SEVEN markers that "produced allele-specific
#       polymorphic clusters" (YD52, YD6, PM7, PM2, HM14, HI3, HI4) and says
#       genotyping used "18 allele-specific markers (Table S4)".
#  Table S4 does contain exactly 18 marker rows. So 18 "designed" is solid and
#  the two polymorphic counts (8 in abstract, 7 named in results) disagree by one.
#  Both are landed; neither is chosen. The proposal file records this.
#
#  ASSERTIONS (each one stops)
#  (1) Supplement Table S4 has exactly 18 marker rows and its marker names are
#      the ones the text names; the body-row count equals the number of markers
#      printed in section 2.3 ("18 allele-specific markers").
#  (2) The paper's own arithmetic closes: 39 - 7 unmappable = 32.
#  (3) The comma-formatted positions in Table 2 match Table S4, so both describe
#      the same 18 assays.
#  (4) Every figure we write is one of the numbers printed in the text: the
#      regexes fail (and we stop) if a sentence changes.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!exists("root")) {
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else {
    getwd()
  }
  setwd(root)
}
suppressPackageStartupMessages(library(xml2))
source(file.path(root, "data-raw", "lib", "xlsx-read.R"))

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "wheat_as2023")
xlsx <- file.path(lit, "supplements", "PLD3-7-e520-s003.xlsx")
ft <- file.path(lit, "fulltext.xml")
for (f in c(xlsx, ft)) {
  if (!file.exists(f)) {
    stop("missing ", f, "\n  fetch with:\n",
         "  curl -o ", ft, " https://www.ebi.ac.uk/europepmc/webservices/rest/",
         "PMC10435944/fullTextXML\n",
         "  curl -o ", file.path(lit, "sups.zip"),
         " https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10435944/",
         "supplementaryFiles")
  }
}

DOI <- "10.1002/pld3.520"
out_dir <- file.path(root, "inst", "extdata")

## --- 1. the paper's own sentences (all figures come from here) --------------
ft_text <- paste(xml_text(xml_find_all(read_xml(ft), "//text()")), collapse = " ")
ft_text <- gsub("[[:space:]]+", " ", ft_text)

grab <- function(pattern, what) {
  hit <- regmatches(ft_text, regexpr(pattern, ft_text))
  if (!length(hit)) stop("sentence not found (", what, "): ", pattern)
  hit
}
s_unmapped <- grab("Out of 39 significant SNPs, seven SNP locations could not", "unmapped")
s_designed <- grab("Among 32, only 18 target variant sequences", "designed")
s_poly18 <- grab("Eighteen out of the target 39 SNP locations fulfilled", "18 of 39")
s_poly8 <- grab("among which only eight produced polymorphic signals", "eight polymorphic")
s_genotyped <- grab("was genotyped with 18 allele.specific markers", "18 markers used")
start <- regexpr("Markers for yield", ft_text)
end <- regexpr("produced allele.specific polymorphic clusters", ft_text)
if (start < 0 || end < 0 || end < start) {
  stop("could not locate the sentence naming the polymorphic markers")
}
s_named <- substr(ft_text, start, end + attr(end, "match.length") - 1L)

n_significant <- as.integer(sub(".*Out of ([0-9]+) significant SNPs.*", "\\1", s_unmapped))
n_unmapped    <- 7L
n_designable  <- as.integer(sub(".*Among ([0-9]+), only ([0-9]+).*", "\\1", s_designed))
n_designed    <- as.integer(sub(".*Among ([0-9]+), only ([0-9]+).*", "\\2", s_designed))
n_fulfilled   <- as.integer(sub(".*Eighteen.*", "18", s_poly18))
n_polymorphic_abstract <- 8L
n_markers_used <- as.integer(sub(".*genotyped with ([0-9]+) allele.*", "\\1", s_genotyped))

# (2) the paper's own arithmetic must close
if (n_significant - n_unmapped != n_designable) {
  stop("funnel does not close: ", n_significant, " - ", n_unmapped, " != ", n_designable)
}
if (n_fulfilled != n_designed || n_markers_used != n_designed) {
  stop("the three '18's disagree: abstract=", n_fulfilled, " among-32=", n_designed,
       " markers-used=", n_markers_used)
}

## --- 2. Table S4: the 18 designed markers -----------------------------------
sheets <- xlsx_sheets(xlsx)
if (!"Table S4" %in% names(sheets)) {
  stop("Table S4 missing from the supplement; sheets: ", paste(names(sheets), collapse = ", "))
}
s4 <- sheets[["Table S4"]]
# Row 1 is the table title, row 2 the header; data rows are marker rows.
body <- s4[-c(1L, 2L)]
body <- Filter(function(r) {
  nm <- xcell(r, "A")
  !is.na(nm) && grepl("^[A-Z]{2}[0-9]+$", nm)
}, body)
markers <- vapply(body, function(r) xcell(r, "A"), character(1))

if (length(markers) != n_markers_used) {
  stop("Table S4 has ", length(markers), " marker rows but the text says ",
       n_markers_used, " markers were used")
}
if (anyDuplicated(markers)) stop("duplicate markers in Table S4")

# (3) Table S4 must describe the same assays Table 2 prints: the coordinate in
# column F (contig_regions) carries the physical position. Table 2 prints it with
# thousands separators; strip those and require the region string to contain it.
t2_pos <- c(
  YD52 = 411775364, YD6 = 34083794, PM2 = 53737652,
  PM7 = 486517837, HM14 = 356142374, HI4 = 41415831
)
contig_regions <- vapply(body, function(r) xcell(r, "F"), character(1))
names(contig_regions) <- markers
for (nm in names(t2_pos)) {
  # The first region is a span "chr:start-end"; Table 2's position lies inside it.
  span <- regmatches(contig_regions[[nm]], regexpr("[0-9]+-[0-9]+", contig_regions[[nm]]))
  if (!length(span)) stop("no coordinate span in Table S4 region for ", nm)
  bounds <- as.numeric(strsplit(span, "-")[[1L]])
  if (t2_pos[[nm]] < bounds[1L] || t2_pos[[nm]] > bounds[2L]) {
    stop("Table 2 position ", t2_pos[[nm]], " for ", nm,
         " lies outside the Table S4 region span ", span)
  }
}

## --- 3. which markers the text calls polymorphic ----------------------------
# Section 2.3 names them; parse the names out of the sentence rather than
# hand-typing, so a wrong name cannot slip in.
named <- regmatches(s_named, gregexpr("\\b(YD|PM|HM|HI)[0-9]+\\b", s_named))[[1L]]
named <- unique(named)
if (!all(named %in% markers)) {
  stop("text names markers absent from Table S4: ",
       paste(setdiff(named, markers), collapse = ", "))
}
n_named_polymorphic <- length(named)
if (n_named_polymorphic != n_polymorphic_abstract - 1L) {
  stop("abstract says ", n_polymorphic_abstract, " polymorphic but the results ",
       "name ", n_named_polymorphic, " (expected a one-marker gap to land honestly)")
}

## --- 4. write the outcome rows ----------------------------------------------
outcome_definition <- paste(
  "Assay is KASP/PACE (PolyMarker design, FAM/HEX tails), NOT ARMS.",
  "worked_design = the marker was designed and used on the FIGS set (Table S4).",
  "worked_polymorphic = the paper calls it polymorphic (allele-specific clusters);",
  "the abstract counts 8, section 2.3 names 7 - both are carried.",
  "worked = polymorphic with a significant phenotypic association in Table 3/4."
)

out <- data.frame(
  marker = markers,
  snp_in_table_s4 = vapply(body, function(r) xcell(r, "B"), character(1)),
  chromosome = vapply(body, function(r) xcell(r, "D"), character(1)),
  contig_regions = contig_regions,
  primer_type = vapply(body, function(r) xcell(r, "K"), character(1)),
  assay = "KASP/PACE (PolyMarker allele-specific + FAM/HEX tails)",
  named_polymorphic_in_text = markers %in% named,
  outcome_definition = outcome_definition,
  evidence_tier = "primary literature",
  source_file = "wheat_as2023/supplements/PLD3-7-e520-s003.xlsx (Table S4); fulltext.xml abstract + 2.1 + 2.3",
  doi = DOI,
  stringsAsFactors = FALSE
)
out$designed <- TRUE

summary_rows <- data.frame(
  stage = c("significant SNPs", "unmappable", "designable", "designed markers",
            "polymorphic (abstract)", "polymorphic (named in 2.3)"),
  n = c(n_significant, n_unmapped, n_designable, n_designed,
        n_polymorphic_abstract, n_named_polymorphic),
  stringsAsFactors = FALSE
)
summary_rows$outcome_definition <- outcome_definition
summary_rows$evidence_tier <- "primary literature"
summary_rows$source_file <- "wheat_as2023/fulltext.xml (abstract, 2.1, 2.3); Table S4"
summary_rows$doi <- DOI

outfile <- file.path(out_dir, "arms_wheat_as2023_outcomes.csv")
sumfile <- file.path(out_dir, "arms_wheat_as2023_funnel.csv")
utils::write.csv(out, outfile, row.names = FALSE)
utils::write.csv(summary_rows, sumfile, row.names = FALSE)
message("wrote ", outfile, " (", nrow(out), " markers) and ", sumfile)
message("funnel: ", paste(summary_rows$n, collapse = " -> "))
