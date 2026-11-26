# =============================================================================
#  data-raw/dnamate_operating_domain.R  --  the domain dnaMATE's Tm model is
#                                           valid over, and the salt term it uses
# =============================================================================
#  Panjkovich A, Norambuena T, Melo F. "dnaMATE: a consensus melting
#  temperature prediction server for short DNA sequences".  Nucleic Acids Res
#  2005; doi:10.1093/nar/gki379, PMC1160140.
#
#  WHY THIS FILE EXISTS
#  dnaMATE publishes no Tm table, so it cannot be cited for a measured value.
#  What it *does* publish is the boundary of its own model: the sequence length,
#  the monovalent-salt window, the accuracy claimed inside that window, the
#  salt-adjustment coefficient, the strand-count term and the three
#  nearest-neighbour parameter sets it averages.  Those are the numbers this
#  package's `primer_tm()` actually shares with it:
#    * `TM_SALT_METHOD` in R/marker_core.R is "Schildkraut2010", and the paper
#      states its salt factor is the Schildkraut-Lifson equation at 16.6 log[Na+];
#    * `NN_TABLE` is DNA_NN_SantaLucia_2004, one of the paper's three tables.
#  Pinning the domain keeps the comparison honest: outside 16-30 nt or 50-600 mM
#  na+, dnaMATE makes no accuracy claim, and neither may this package.
#
#  FETCH (verified route, 2026-11)
#    Europe PMC fullTextXML for PMC1160140 returns 500, so the page is the
#    source:  curl -sS -L -m 60 -o .tmp/lit/dnamate/page.html \
#      'https://pmc.ncbi.nlm.nih.gov/articles/PMC1160140/'
#
#  Run:  Rscript data-raw/dnamate_operating_domain.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

doi <- "10.1093/nar/gki379"
pmcid <- "PMC1160140"

page_path <- Sys.getenv(
  "DNAMATE_HTML",
  file.path(root, ".tmp", "lit", "dnamate", "page.html")
)
if (!file.exists(page_path)) {
  stop("missing ", page_path, "; see the FETCH block at the top of this script")
}

# Plain text of the page, used to gate every transcribed value.  xml2 decodes
# the numeric entities (&#xB0; for the degree sign) that a hand-written
# entity table would miss.
flat <- function(path) {
  x <- xml2::xml_text(xml2::read_html(path, encoding = "UTF-8"))
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  gsub("[[:space:]]+", " ", x)
}
page <- flat(page_path)

# The statements the extract is copied from, verbatim.  Each must occur once.
evidence <- c(
  length = "DNA sequences within the length range of 16\u201330 nt",
  batch = "up to 5000 DNA sequences can be submitted in a single run",
  salt_window = "in the range of 50\u2013600 mM monovalent salt concentration",
  accuracy = "89% of the melting temperature predictions will have an error or deviation of <5\u00b0C",
  gas_constant = "R is the gas constant (fixed at 1.987 cal/K mol)",
  strand_term = "Constant b adopts the value of 4 for non-self-complementary sequences or equal to 1 for duplexes of self-complementary strands",
  salt_rule = "The Schildkraut\u2013Lifson (8) equation is used as the salt adjustment factor, which corresponds to 16.6log",
  tables = "tables from Breslauer (5), SantaLucia (6) and Sugimoto (7)",
  consensus_error = "The Tm estimation error at the other regions where some consensus was observed is expected to be small (below 3\u20135\u00b0C)",
  black_zone = "The Tm estimations of oligonucleotides falling into the black regions of the consensus map by any of the methods could have a large error"
)
hits <- vapply(evidence, function(e) lengths(regmatches(page, gregexpr(e, page, fixed = TRUE))), integer(1))
if (any(hits == 0L)) {
  stop("evidence not found in ", page_path, ": ", paste(names(hits)[hits == 0L], collapse = ", "))
}

# `value` stays character: two rows are lists, not numbers.  Callers parse the
# numeric ones with as.numeric().
row <- function(item, value, unit, reported, section) {
  data.frame(
    item = item, value = as.character(value), unit = unit,
    reported = reported, section = section,
    doi = doi, pmcid = pmcid, stringsAsFactors = FALSE
  )
}

out <- rbind(
  row("sequence_length_min", 16L, "nt",
      "DNA sequences within the length range of 16\u201330 nt", "Abstract"),
  row("sequence_length_max", 30L, "nt",
      "DNA sequences within the length range of 16\u201330 nt", "Abstract"),
  row("monovalent_salt_min", 50L, "mM",
      "in the range of 50\u2013600 mM monovalent salt concentration", "Abstract"),
  row("monovalent_salt_max", 600L, "mM",
      "in the range of 50\u2013600 mM monovalent salt concentration", "Abstract"),
  row("accuracy_within_5C", 89L, "percent",
      "89% of the melting temperature predictions will have an error or deviation of <5\u00b0C from experimental data", "Abstract"),
  row("server_batch_limit", 5000L, "sequences",
      "up to 5000 DNA sequences can be submitted in a single run", "Abstract"),
  row("salt_adjustment_coefficient", 16.6, "log10[Na+] per decade",
      "The Schildkraut\u2013Lifson (8) equation is used as the salt adjustment factor, which corresponds to 16.6log[Na+]", "Methods"),
  row("gas_constant", 1.987, "cal/K/mol",
      "R is the gas constant (fixed at 1.987 cal/K mol)", "Methods"),
  row("b_nonself_complementary", 4L, "dimensionless",
      "Constant b adopts the value of 4 for non-self-complementary sequences", "Methods"),
  row("b_self_complementary_or_excess", 1L, "dimensionless",
      "or equal to 1 for duplexes of self-complementary strands or for duplexes when one of the strands is in significant excess", "Methods"),
  row("nn_tables", "Breslauer;SantaLucia;Sugimoto", "table names",
      "The thermodynamic data ... include the tables from Breslauer (5), SantaLucia (6) and Sugimoto (7)", "Methods"),
  row("consensus_zone_error_bound", 5L, "degC",
      "The Tm estimation error at the other regions where some consensus was observed is expected to be small (below 3\u20135\u00b0C)", "Results"),
  row("black_zone", "no consensus among any method; large error possible", "qualitative",
      "The Tm estimations of oligonucleotides falling into the black regions of the consensus map by any of the methods could have a large error", "Results")
)

stopifnot(
  nrow(out) == 13L,
  !anyDuplicated(out$item),
  all(nzchar(out$reported)),
  all(out$doi == doi),
  all(out$pmcid == pmcid),
  all(out$section %in% c("Abstract", "Methods", "Results")),
  as.numeric(out$value[out$item == "sequence_length_min"]) == 16,
  as.numeric(out$value[out$item == "sequence_length_max"]) == 30,
  as.numeric(out$value[out$item == "monovalent_salt_min"]) == 50,
  as.numeric(out$value[out$item == "monovalent_salt_max"]) == 600,
  as.numeric(out$value[out$item == "accuracy_within_5C"]) == 89,
  as.numeric(out$value[out$item == "server_batch_limit"]) == 5000,
  as.numeric(out$value[out$item == "salt_adjustment_coefficient"]) == 16.6,
  as.numeric(out$value[out$item == "gas_constant"]) == 1.987,
  as.numeric(out$value[out$item == "b_nonself_complementary"]) == 4,
  as.numeric(out$value[out$item == "b_self_complementary_or_excess"]) == 1,
  as.numeric(out$value[out$item == "consensus_zone_error_bound"]) == 5
)

readr::write_csv(out, "inst/extdata/dnamate_operating_domain.csv", na = "NA")
cat(
  "dnamate_operating_domain:", nrow(out), "rows;",
  sprintf("%s-%s nt, %s-%s mM na+; %s%% within 5 C over that window\n",
          out$value[out$item == "sequence_length_min"],
          out$value[out$item == "sequence_length_max"],
          out$value[out$item == "monovalent_salt_min"],
          out$value[out$item == "monovalent_salt_max"],
          out$value[out$item == "accuracy_within_5C"])
)
