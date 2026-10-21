# =============================================================================
#  data-raw/sanger_floor.R  --  the Sanger detection floor, tiered by evidence
# =============================================================================
#  The package must never carry one "Sanger detects X%" number.  The floor moves
#  with the basecaller, the caller's mixing-ratio model, amplicon context and
#  whether the laboratory validated it.  This file collects what can actually be
#  read from a fetched file, tags each row with an evidence tier, and says which
#  rows are admissible as calibration targets.
#
#  Every number below is re-read from a fetched full text under .tmp/lit/:
#
#    PMC12188199  Improved Detection of HIV-1 Drug Resistance With Droplet
#                 Digital PCR.  (no DOI in the EPMC record)  -> "<20%"
#    PMC13410318  Ion Torrent Genexus ... Comparison with the GeneStudio S5.
#                 DOI 10.3390/ijms27146307 -> "below approximately 20%"
#    mosha2014    Mosha et al. 2014, BMC Infect Dis 14(Suppl 2):O1,
#                 DOI 10.1186/1471-2334-14-s2-o1 -> ABSTRACT ONLY (see below)
#
#  docs/ has the earlier notes; IMPORTANT negative results recorded here:
#
#   * liew2004 (10.1373/clinchem.2004.032136) is CLOSED ACCESS. Nothing was
#     fetched, so no Liew HRM number appears anywhere in this deliverable.
#   * qad.0000000000000267 (ultra-deep vs Sanger minority-mutation detection)
#     returned HTTP 500 from the EPMC fullTextXML endpoint; nothing fetched.
#   * PMC8621402 (10.3390/ijms222212527) was already landed by data-raw/
#     ijms2022.R from the same full text - the HRM panel file is that one.
#   * mosha2014's .tmp/lit/mosha2014/fulltext.xml is article-type="abstract";
#     its supplements.zip is 201 bytes and is NOT a zip.  The often-quoted
#     "8% clinical measurement" and a named vendor's "5%" claim are therefore
#     NOT landed here: they could not be read from a fetched file.  They are
#     listed in SOURCE-BLOCKED form in .tmp/sdd/report-readout.md instead.
#
#  ASSERTIONS
#   * Each of the two fetched articles must state its own floor sentence, and
#     the number must be extractable as 20.
#   * The two sentences must be WORDED differently (one "<20%", one "below
#     approximately 20%"), proving they are two restatements, not one string
#     copied twice; both are tier 'textbook' because neither author measured it.
#   * Every row's usable_as_calibration flag must agree with its tier.
#
#  Run:  Rscript data-raw/sanger_floor.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

paths <- c(
  hiv_ddpcr_review = fs::path(root, ".tmp", "lit", "PMC12188199", "fulltext.xml"),
  iontorrent = fs::path(root, ".tmp", "lit", "PMC13410318", "fulltext.xml"),
  mosha2014 = fs::path(root, ".tmp", "lit", "mosha2014", "fulltext.xml"))
miss <- names(paths)[!file.exists(paths)]
if (length(miss)) stop("Missing fetched file(s): ", paste(miss, collapse = ", "))

paras_of <- function(p) {
  d <- xml2::read_xml(p)
  gsub("[[:space:]]+", " ", xml2::xml_text(xml2::xml_find_all(d, "//*[local-name()='p']")))
}
hit_one <- function(ps, re) {
  h <- grep(re, ps, value = TRUE)
  if (length(h) < 1L) stop("Floor sentence not found for pattern: ", re)
  h[[1L]]
}

# --- measured-by-proxy statement: neither paper measured the floor ------------
p1 <- paras_of(paths[["hiv_ddpcr_review"]])
p2 <- paras_of(paths[["iontorrent"]])
s1 <- hit_one(p1, "Sanger sequencing, the current gold standard")
s2 <- hit_one(p2, "reliably detect variants that occur at intra-host frequencies")
if (!any(grepl("Sanger sequencing cannot detect minor variants constituting <20%", p1, fixed = TRUE)))
  stop("PMC12188199 no longer states the <20% floor verbatim")
if (!grepl("below approximately 20%", s2, fixed = TRUE))
  stop("PMC13410318 no longer states the 'below approximately 20%' floor")
if (identical(s1, s2)) stop("The two floor sentences are identical; expected two sources")
n1 <- as.numeric(sub(".*<([0-9.]+)%.*", "\\1", grep("constituting <[0-9.]+%", p1, value = TRUE)[[1L]]))
n2 <- as.numeric(sub(".*approximately ([0-9.]+)%.*", "\\1", s2))
if (!is.finite(n1) || !is.finite(n2)) stop("Could not parse the floor number")
if (n1 != 20 || n2 != 20) stop("Parsed floors are ", n1, "/", n2, ", expected 20/20")

# mosha2014 is an abstract; assert that, so nobody later mistakes it for full text.
m <- xml2::read_xml(paths[["mosha2014"]])
if (!identical(xml2::xml_attr(m, "article-type"), "abstract"))
  stop("mosha2014 fulltext.xml is no longer article-type='abstract'")

# --- rows --------------------------------------------------------------------
rows <- tibble::tibble(
  metric = c("detection_floor_pct", "detection_floor_pct"),
  value_pct = c(n1, n2),
  lower_bound = c(FALSE, FALSE),
  organism = c("HIV-1", "HIV-1"),
  method = c("population Sanger sequencing", "population Sanger sequencing"),
  basis = c("restated as the gold-standard limitation in the Introduction",
            "restated as the conventional-GRT limitation in the Introduction"),
  statement = c("Sanger sequencing cannot detect low-abundance variants constituting <20% of the virus population",
                "this method struggles to reliably detect variants that occur at intra-host frequencies below approximately 20%"),
  evidence_tier = c("textbook", "textbook"),
  usable_as_calibration = c(FALSE, FALSE),
  reason = c("not measured by these authors; no titration, no per-mix replicate series",
             "stated as an approximation ('approximately'); not measured here"),
  doi = c("NA", "10.3390/ijms27146307"),
  source_file = c(".tmp/lit/PMC12188199/fulltext.xml", ".tmp/lit/PMC13410318/fulltext.xml"))
if (any(rows$usable_as_calibration != (rows$evidence_tier == "measured")))
  stop("usable_as_calibration disagrees with evidence_tier")

out <- dplyr::mutate(rows, outcome_definition = paste0(
  "Outcome = the VARIANT ALLELE FRACTION below which a population (bulk) Sanger",
  "trace is said not to call a minor variant. A statement about calling a",
  "frequency, not about one base position's signal-to-noise; it presupposes a",
  "particular basecaller and mixing-ratio (peak-height) model, both unstated.",
  " Usable as a calibration target only when tier == 'measured'."))
readr::write_csv(out, fs::path(root, "inst", "extdata", "sanger_detection_floor.csv"))

message("sanger_floor: ", nrow(rows), " row(s), tiers ",
        paste(unique(rows$evidence_tier), collapse = "/"), "; both sources say ",
        n1, "% but neither measured it; mosha2014 is abstract-only -- assertions passed")
