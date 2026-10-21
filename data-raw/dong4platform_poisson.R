# =============================================================================
#  data-raw/dong4platform_poisson.R  --  Poisson PREDICTION vs MEASUREMENT
# =============================================================================
#  Source (DOI is PROVENANCE; every number is re-read from the fetched files):
#
#    Comparison of four digital PCR platforms for accurate quantification of DNA
#    copy number of a certified plasmid DNA reference material.  Sci Rep 2015.
#    DOI 10.1038/srep13174   EPMC PMC4548226
#    Fetched: .tmp/lit/PMC4548226/fulltext.xml   (Europe PMC fullTextXML)
#             .tmp/lit/PMC4548226/supplements/srep13174-s1.doc
#               (Europe PMC supplementaryFiles; Tables S5-S8 are the raw counts)
#
#  WHY THIS IS THE THEORY-vs-MEASUREMENT ANCHOR FOR THE POISSON MODEL
#  tidymarkers converts a positive-partition fraction to a concentration with the
#  Poisson model  lambda = -ln(1 - p).  Supplementary Tables S5-S8 print, for FOUR
#  platforms x 15 replicates, the RAW positive/negative partition counts AND the
#  instrument's own mean copies/partition.  So we can recompute lambda from the
#  counts and compare it with the printed value: the residual is the measured
#  deviation of the Poisson model from observation, at lambda ~ 1.4-1.6.
#
#  The .doc is Word 97.  There is no antiword/catdoc/LibreOffice on this machine
#  and pandoc cannot read .doc, so data-raw/lib/word97-text.R decodes the OLE2
#  container and the piece table.  That decoder was validated against an
#  independent Python implementation (same 18 343 characters) before use.
#
#  WHAT IS *NOT* LANDED
#   * Table 2's "Partition number" column (765 / 13800+/-464 / 64 / 1695000) is a
#     design figure, not a per-replicate observation, so it is not used.
#   * The linearity claim (r2 = 0.999, p. 18 of the article) is about a FIGURE;
#     the figure is an image and no number was read from it.
#
#  ASSERTIONS
#   * 4 tables x 15 replicates = 60 rows, and negative + positive must equal the
#     printed total in every row (the tables print all three, so it self-checks).
#   * The per-platform mean lambda recomputed from the counts must reproduce
#     Table 2's printed lambda (1.56 / 1.54 / 1.54 / 1.51) to 2 d.p.
#   * |lambda_poisson - printed| must stay under 0.01 for every replicate.
#
#  Run:  Rscript data-raw/dong4platform_poisson.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
source(fs::path(root, "data-raw", "lib", "word97-text.R"))

doc_path <- fs::path(root, ".tmp", "lit", "PMC4548226", "supplements", "srep13174-s1.doc")
xml_path <- fs::path(root, ".tmp", "lit", "PMC4548226", "fulltext.xml")
for (p in c(doc_path, xml_path)) if (!file.exists(p)) stop("Missing ", p)
doc <- xml2::read_xml(xml_path)
if (!grepl("10.1038/srep13174", xml2::xml_text(doc), fixed = TRUE)) {
  stop("fulltext.xml does not mention the DOI")
}

# --- Tables S5-S8 ------------------------------------------------------------
text <- wd_text(doc_path)
if (!grepl("Table S5", text, fixed = TRUE)) {
  stop("Word text extraction lost the Table S5 label; refusing to parse")
}
lines <- strsplit(text, "\r", fixed = TRUE)[[1L]]

edges <- grep("^Vials\x07", lines)
if (length(edges) != 4L) stop("expected 4 tables (S5-S8), found ", length(edges))
plat <- c("QuantStudio 12k", "BioMark", "QX100", "RainDrop")
labels <- c("Table S5", "Table S6", "Table S7", "Table S8")
rows <- list()
for (k in seq_along(edges)) {
  if (!grepl(labels[k], lines[edges[k] - 1L], fixed = TRUE)) {
    stop(labels[k], " label is not on the line above its table")
  }
  cells <- trimws(strsplit(lines[edges[k]], "\x07", fixed = TRUE)[[1L]])
  cells <- cells[nzchar(cells)]
  body <- cells[-(1:6)] # drop the 6-column header
  # Rows are [vial?] rep neg pos mean total; a vial label only appears at the
  # start of a new three-replicate block.
  j <- 1L
  cur <- "Vial-1"
  while (j + 4L <= length(body)) {
    if (!grepl("^[0-9]+$", body[j])) {
      cur <- body[j]
      j <- j + 1L
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      table = labels[k], platform = plat[k], vial = cur,
      replicate = as.integer(body[j]),
      negative = as.numeric(body[j + 1L]), positive = as.numeric(body[j + 2L]),
      reported_mean = as.numeric(body[j + 3L]),
      reported_total = as.numeric(body[j + 4L]), stringsAsFactors = FALSE
    )
    j <- j + 5L
  }
}
d <- do.call(rbind, Filter(Negate(is.null), rows))
if (nrow(d) != 60L) stop("expected 60 replicates, got ", nrow(d))
if (any(is.na(d$reported_total)) || any(is.na(d$negative))) {
  stop("a table row did not parse to six numeric cells")
}
if (!all(table(d$platform) == 15L)) stop("not 15 replicates per platform")

d$total <- d$negative + d$positive
if (any(d$total != d$reported_total)) {
  stop("negative + positive != printed total in ", sum(d$total != d$reported_total), " rows")
}
d$lambda_poisson <- -log(d$negative / d$total)
d$lambda_residual <- d$lambda_poisson - d$reported_mean
if (max(abs(d$lambda_residual)) > 0.01) {
  stop(
    "Poisson recomputation drifts from the printed mean by ",
    round(max(abs(d$lambda_residual)), 4)
  )
}

# Table 2 of the article prints lambda per platform to 2 d.p.
t2_lambda <- c("QuantStudio 12k" = 1.54, "BioMark" = 1.56, "QX100" = 1.54, "RainDrop" = 1.51)
per <- do.call(rbind, lapply(plat, function(g) {
  v <- d[d$platform == g, ]
  data.frame(
    platform = g, n = nrow(v),
    partitions_mean = round(mean(v$total)),
    lambda_reported = round(mean(v$reported_mean), 4),
    lambda_poisson = round(mean(v$lambda_poisson), 4),
    residual_mean = round(mean(v$lambda_residual), 3),
    residual_sd = round(sd(v$lambda_residual), 3),
    residual_max_abs = round(max(abs(v$lambda_residual)), 3),
    table2_lambda = unname(t2_lambda[[g]]),
    stringsAsFactors = FALSE
  )
}))
# Table 2 prints lambda to 2 d.p.; allow one unit in the last place.
# (RainDrop: mean of the 15 printed replicates is 1.5153, printed as 1.51).
if (any(abs(per$lambda_reported - per$table2_lambda) > 0.01)) {
  stop(
    "per-platform lambda does not match Table 2: ",
    paste(per$lambda_reported, "/", per$table2_lambda, collapse = " ")
  )
}
if (any(abs(per$lambda_poisson - per$table2_lambda) > 0.01)) {
  stop("per-platform lambda recomputed from counts does not match Table 2")
}

DEF <- paste(
  "Outcome = (a) the instrument's own mean copies per partition for ONE",
  "replicate, and (b) lambda_poisson = -ln(N_negative / N_total) from the SAME",
  "printed counts. lambda_residual = (b) - (a) is the Poisson-model deviation at",
  "that loading. This is a MEASUREMENT of model-vs-observation, not an assay",
  "detection limit."
)
d$outcome_definition <- DEF
d$evidence_tier <- "measured"
d$doi <- "10.1038/srep13174"
d$source_file <- "supplements/srep13174-s1.doc (Tables S5-S8)"
d <- d[, c(
  "table", "platform", "vial", "replicate", "negative", "positive",
  "total", "reported_total", "reported_mean", "lambda_poisson", "lambda_residual",
  "outcome_definition", "evidence_tier", "doi", "source_file"
)]
readr::write_csv(d, fs::path(root, "inst", "extdata", "dpcr_dong_poisson_residual.csv"))

per$outcome_definition <- paste(
  "Outcome = mean over 15 replicates of the Poisson-model prediction from the",
  "printed partition counts, compared with the instrument's printed lambda, on a",
  "certified plasmid reference material across four dPCR platforms."
)
per$evidence_tier <- "measured"
per$doi <- "10.1038/srep13174"
per$source_file <- "supplements/srep13174-s1.doc (Tables S5-S8) + Table 2"
readr::write_csv(per, fs::path(root, "inst", "extdata", "dpcr_dong_poisson_summary.csv"))

message(
  "dong4platform_poisson: ", nrow(d), " replicates over 4 platforms; max |residual| ",
  round(max(abs(d$lambda_residual)), 4), "; mean lambda ",
  paste(sprintf("%.2f", per$lambda_reported), collapse = "/"),
  " reproduces Table 2 -- assertions passed"
)
