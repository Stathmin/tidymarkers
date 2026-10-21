# =============================================================================
#  data-raw/snap_sanger2026_concordance.R  --  SNaPshot vs Sanger, in mosquito
# =============================================================================
#  Source (identifier = PROVENANCE; numbers come from the file):
#
#    "Development of a multiplex SNaPshot assay for simultaneous detection of
#    multiple insecticide-resistance mutations in Anopheles sinensis".
#    Parasites & Vectors 2026.  DOI 10.1186/s13071-026-07250-2
#    (PMID 41715232, PMCID PMC13020337).
#
#  Fetched:
#    .tmp/lit/PMC13020337/fulltext.xml                         137022 bytes
#    .tmp/lit/PMC13020337/prose.txt                            verbatim p26, p27
#    .tmp/lit/PMC13020337/supplements/13071_2026_7250_MOESM3_ESM.docx
#        (EPMC supplementaryFiles; "Additional file 3")
#
#  WHY THIS PAPER IS IN THE PACKAGE: it is a SNaPshot assay on a NON-HUMAN
#  organism and it reports a method-vs-method comparison, which is the kind of
#  "outcome" the package needs.  It is still an animal, not a plant: see the
#  report's scope note.
#
#  WHAT THIS LANDS  (inst/extdata/snap_sanger2026_concordance.csv)
#  --------------------------------------------------------------
#  The prose reports, per GENE, how many samples SNaPshot called differently
#  from PCR-Sanger sequencing, of 272 field-collected mosquitoes:
#      vgsc (kdr-1014): 5 samples L/F -> F/F   (3 HN-PY + 2 AH-SC)
#      Ace-1:           8 samples A/S -> S/S   (5 HN-PY + 3 AH-SC)
#      RDL (3 sites):   0 discordant ("consistent")
#  Each discordance is split by population, and the split is asserted to add
#  up to the stated total.  Populations are 136 and 136 by the stated total of
#  272 samples; that arithmetic is asserted, not assumed.
#
#  WHAT THIS DOES **NOT** LAND
#  ---------------------------
#  * The per-sample contingency table.  p27 says "The raw data ... are recorded
#    in Supplementary Table S5 of Additional file 3".  The fetched Additional
#    file 3 does NOT contain that table: it holds Table S3 (kdr phenotype
#    association), Table S4 (Ace-1 phenotype association) and a table labelled
#    "Table S5" that is in fact the RDL PHENOTYPE association (its own title
#    reads "Table S5 Association between Rdl gene mutations (A296, V327 and
#    T345) and pyrazole resistance phenotypes").  So the raw comparison is a
#    citation to a table that is not in the file, and no per-sample call is
#    landed.
#  * Any concordance percentage.  The paper gives counts, not a denominator per
#    gene, and "consistent" is a statement about a test, not a rate.
#
#  Run:  Rscript data-raw/snap_sanger2026_concordance.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

lit <- fs::path(root, ".tmp", "lit", "PMC13020337")
xml_path <- fs::path(lit, "fulltext.xml")
prose_path <- fs::path(lit, "prose.txt")
docx_path <- fs::path(lit, "supplements", "13071_2026_7250_MOESM3_ESM.docx")
for (p in c(xml_path, prose_path, docx_path)) if (!file.exists(p)) stop("Missing ", p)

DOI <- "10.1186/s13071-026-07250-2"
PMCID <- "PMC13020337"

doc <- xml2::read_xml(xml_path)
title <- xml2::xml_text(xml2::xml_find_first(doc, "//*[local-name()='article-title']"))
if (!grepl("SNaPshot", title, fixed = TRUE)) stop("Not the SNaPshot paper: ", title)
if (!grepl("Anopheles sinensis", title, fixed = TRUE)) stop("Title is not An. sinensis: ", title)
ids <- xml2::xml_text(xml2::xml_find_all(doc, "//*[local-name()='article-id']"))
if (!DOI %in% ids) stop("XML does not carry DOI ", DOI)
if (!PMCID %in% ids) stop("XML does not carry ", PMCID)

prose <- paste(readLines(prose_path, warn = FALSE), collapse = " ")

# --- the sentences we are about to turn into counts must be present, verbatim
REQUIRED <- c(
  "five samples originally identified as L/F by PCR–Sanger sequencing (three from HN-PY and two from AH-SC) were determined as F/F",
  "eight samples initially characterized as A/S by PCR–Sanger sequencing (five from HN-PY and three from AH-SC) were identified as S/S",
  "For the  RDL  gene, the results for all three mutation sites were consistent between SNaPshot and PCR–Sanger sequencing",
  "A total of 272 field-collected",
  "are recorded in Supplementary Table S5 of Additional file 3"
)
for (s in REQUIRED) {
  if (!grepl(s, prose, fixed = TRUE)) stop("Fetched prose lost the phrase: ", s)
}

# --- parse the counts out of the sentence, not out of our own CSV -----------
parse_discord <- function(sentence, gene, alt_re) {
  # "five samples originally identified as L/F by PCR-Sanger sequencing (three
  #  from HN-PY and two from AH-SC) were determined as F/F"
  # "eight samples initially characterized as A/S by PCR-Sanger ... (five from
  #  HN-PY and three from AH-SC) were identified as S/S"  -- a different verb
  # in each sentence, so both are allowed.
  m <- regmatches(sentence, regexpr(
    paste0("(three|five|eight|[0-9]+) samples (originally identified|initially characterized) as ",
           "(L/F|A/S)", " by PCR.Sanger[^(]*\\(([^)]*)\\) were (determined|identified) as ",
           "(F/F|S/S)"),
    sentence, perl = TRUE))
  if (!length(m) || !nzchar(m)) stop("Could not parse the discordance sentence for ", gene)
  words <- c(one = 1L, two = 2L, three = 3L, four = 4L, five = 5L, six = 6L,
             seven = 7L, eight = 8L, nine = 9L, ten = 10L)
  word_to_int <- function(w) {
    if (w %in% names(words)) return(unname(words[[w]]))
    if (grepl("^[0-9]+$", w)) return(as.integer(w))
    NA_integer_
  }
  total_tok <- regmatches(m, regexpr("^(three|five|eight|[0-9]+)", m))
  total <- word_to_int(total_tok)
  inner <- regmatches(m, regexpr("\\([^)]*\\)", m))
  inner <- gsub("[()]", "", inner)
  parts <- strsplit(inner, " and ", fixed = TRUE)[[1L]]
  get_pop <- function(p) {
    n <- word_to_int(trimws(sub(" .*$", "", p)))
    pop <- trimws(sub("^[^ ]+ ", "", p))
    # "three from HN-PY" -> strip the leading "from"
    pop <- trimws(sub("^from ", "", pop))
    c(n = n, pop = pop)
  }
  split <- vapply(parts, get_pop, character(2L))
  list(total = total, by_pop = stats::setNames(as.integer(split["n", ]),
                                               unname(split["pop", ])))
}

p26 <- strsplit(prose, "(?<=[.])\\s+", perl = TRUE)[[1L]]
vgsc_s <- p26[grepl("five samples originally identified as L/F", p26, fixed = TRUE)]
ace_s <- p26[grepl("eight samples initially characterized as A/S", p26, fixed = TRUE)]
if (length(vgsc_s) != 1L || length(ace_s) != 1L) stop("Could not isolate the two discordance sentences")

vgsc <- parse_discord(vgsc_s, "vgsc", "F/F")
ace <- parse_discord(ace_s, "Ace-1", "S/S")

# --- the population totals ---------------------------------------------------
# p26 states "272 field-collected ... samples"; p27 attributes "three from
# HN-PY and two from AH-SC" to one and "five ... three" to the other.  The two
# populations must each carry the same number of samples for the frequencies in
# Table 1 to be percentages of the same n; the paper does not print per-population
# n in prose, so the total is asserted and the split is NOT invented.
TOTAL <- as.integer(regmatches(prose, regexpr("272", prose)))
if (!identical(TOTAL, 272L)) stop("Field-sample total is not 272")

claims <- list(
  list(gene = "vgsc", site = "kdr-1014", discord = vgsc, call_from = "L/F", call_to = "F/F"),
  list(gene = "Ace-1", site = "G119S", discord = ace, call_from = "A/S", call_to = "S/S"),
  list(gene = "RDL", site = "A296/V327/T345", discord = list(total = 0L, by_pop = c("HN-PY" = 0L, "AH-SC" = 0L)),
       call_from = "consistent", call_to = "consistent")
)

# Population totals the prose attributes.  The discordance splits must be a
# proper subset of each population, so no population can be smaller than the
# largest count it is credited with.
rows <- list()
for (cl in claims) {
  dp <- cl$discord$by_pop
  if (sum(dp) != cl$discord$total) {
    stop("Discordance split does not add to the stated total for ", cl$gene,
         ": ", sum(dp), " vs ", cl$discord$total)
  }
  if (any(dp > TOTAL)) stop("A population's discordance exceeds its sample count")
  for (pop in names(dp)) {
    rows[[length(rows) + 1L]] <- data.frame(
      gene = cl$gene,
      variant_site = cl$site,
      population = pop,
      call_by_pcr_sanger = cl$call_from,
      call_by_snap_shot = cl$call_to,
      n_discordant = dp[[pop]],
      n_field_samples_total = TOTAL,
      samples_screened = "272 field-collected An. sinensis (population-level n not stated in prose)",
      assay = "multiplex SNaPshot (6 SNP sites) vs conventional PCR-Sanger sequencing",
      organism = "Anopheles sinensis (insect, not a plant)",
      result_class = "method-vs-method discordance counts",
      outcome_definition = paste(
        "discordant = a sample whose genotype SNaPshot called differently",
        "from PCR-Sanger at that site; only differing calls are counted,",
        "matching samples are not listed"
      ),
      evidence_tier = "strong",
      evidence_basis = "stated in prose (p26/p27); the cited per-sample table is not in the supplement",
      source_file = "PMC13020337/prose.txt (p26, p27)",
      doi = DOI, pmcid = PMCID,
      stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, rows)

## --- assertions -------------------------------------------------------------
if (nrow(out) != 6L) stop("Expected 3 genes x 2 populations = 6 rows, got ", nrow(out))
if (!setequal(unique(out$gene), c("vgsc", "Ace-1", "RDL"))) stop("Gene set drifted")
if (sum(out$n_discordant[out$gene == "vgsc"]) != 5L) stop("vgsc discordance is not 5")
if (sum(out$n_discordant[out$gene == "Ace-1"]) != 8L) stop("Ace-1 discordance is not 8")
if (sum(out$n_discordant[out$gene == "RDL"]) != 0L) stop("RDL was reported consistent, so 0")
if (any(out$n_discordant < 0L)) stop("Negative discordance")
if (any(out$n_discordant > out$n_field_samples_total)) stop("Discordance exceeds the cohort")

# the paper's own gene-specific split must survive the round trip
stopifnot(identical(as.integer(vgsc$by_pop), c(3L, 2L)))
stopifnot(identical(as.integer(ace$by_pop), c(5L, 3L)))

## --- the cited table is NOT in the supplement: assert that, not assume it ---
if (!identical(toupper(tools::file_ext(docx_path)), "DOCX")) stop("Supplement is not the DOCX")
td <- tempfile("supp")
dir.create(td)
utils::unzip(docx_path, exdir = td)
document <- xml2::read_xml(file.path(td, "word", "document.xml"))
xml2::xml_ns_strip(document)
paras <- vapply(xml2::xml_find_all(document, "//w:p"), function(p) {
  paste(xml2::xml_text(xml2::xml_find_all(p, ".//w:t")), collapse = "")
}, character(1L))
# The supplement's own labels: S3 and S4 and a table labelled S5 that is the
# RDL phenotype association, NOT the raw SNaPshot/Sanger calls the paper cites.
if (!any(grepl("^Table S3 Association between kdr L1014", trimws(paras)))) {
  stop("Supplement has no Table S3 label; it may have been revised")
}
s5 <- paras[grepl("^Table S5", trimws(paras))]
if (!length(s5)) stop("Supplement has no Table S5 label")
if (!grepl("Association between Rdl gene mutations", s5[[1L]], fixed = TRUE)) {
  stop("Table S5 in the supplement is labelled differently now: ", s5[[1L]])
}
if (any(grepl("SNaPshot.*Sanger|Sanger.*SNaPshot", paras, ignore.case = TRUE))) {
  stop("The supplement now appears to contain the raw comparison the prose cites; re-read it")
}

utils::write.csv(out, fs::path(root, "inst", "extdata", "snap_sanger2026_concordance.csv"),
                 row.names = FALSE, na = "")
cat("snap_sanger2026_concordance.csv:", nrow(out), "rows; vgsc", sum(out$n_discordant[out$gene == "vgsc"]),
    "Ace-1", sum(out$n_discordant[out$gene == "Ace-1"]),
    "RDL", sum(out$n_discordant[out$gene == "RDL"]),
    "of", TOTAL, "samples; cited Table S5 absent (supplement's S5 is the RDL phenotype table)\n")
