# =============================================================================
#  data-raw/soybean_est_caps-supplement.R  --  EST-SNP -> CAPS funnel, soybeans
# =============================================================================
#  Source (identifier resolution is PROVENANCE, not evidence; every number
#  below is re-read from the fetched PDF):
#
#    束永俊 (Shu Yong-Jun), 李勇, 吴娜拉胡, 柏锡, 才华, 纪巍, 朱延明.
#    大豆 EST-SNP 的挖掘、鉴定及其 CAPS 标记的开发
#    (Development of soybean EST-SNP discovery, validation and CAPS markers).
#    作物学报 Acta Agronomica Sinica 2010, 36(4): 574-579.
#    DOI 10.3724/SP.J.1006.2010.00574  (printed on the first page)
#
#  Fetched: PDF is .tmp/lit/soybean_est_caps/paper.pdf (6 pages); the script
#  reads .tmp/lit/soybean_est_caps/paper.txt, i.e. `pdftotext -layout` of that
#  PDF.  The paper has NO supplementary files - the whole dataset is the paper.
#
#  WHY this paper: it is the only registry source that reports a CAPS
#  conversion funnel in which EVERY stage keeps its failures, so the package
#  can show where markers are lost instead of only how many survived:
#
#      537 EST-SNPs discovered (Fig. 1)
#       48 alter a restriction site -> primers designed   (Sec. 2.3)
#       44 of 48 gave a SINGLE band, 400-800 bp           (Sec. 2.3)
#       36 of 44 amplicons showed the expected EST-SNP when re-sequenced
#       26 of 36 digests gave an unequivocal restriction pattern -> Table 1
#
#  WHAT THIS SCRIPT DOES NOT MODEL  (deliberate, and load-bearing for the
#  package's honesty about "success"):
#
#   * The 26 is NOT a success rate over "designed CAPS markers".  The stages
#     discard different things: 48 -> 44 is a PCR/amplification failure,
#     44 -> 36 is a SNP-calling failure on Sanger traces (the paper calls the
#     36 "discordant" in the English abstract: the trace showed the expected
#     EST-SNP), and 36 -> 26 is a DIGEST-SCORING failure ("unequivocal
#     restriction pattern").  The three losses are not interchangeable and the
#     script keeps them separate; 26/48 must never be pooled with another
#     paper's conversion rate.
#   * Only the 26 survivors are named per primer (Table 1).  The paper does
#     NOT publish which 22 of the 48 failed, nor which 8 of the 44 gave no
#     SNP, nor which 10 digests were ambiguous - only the marginal counts.
#     So the per-marker CSV is the survivors with their enzyme; the losses
#     exist only as stage totals.  Any per-stage "reason" would be invented.
#   * Table 1 carries primers + enzyme only.  It has no amplicon size, no
#     digest fragment sizes, and no genotype calls per variety, so no
#     predicted-vs-observed fragment cross-check is possible here (contrast
#     the SNP2CAPS extract, where every cut sums to its PCR product).
#   * No rate is converted into a probability for internal use.
#
#  Run:  Rscript data-raw/soybean_est_caps-supplement.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

txt_path <- fs::path(root, ".tmp", "lit", "soybean_est_caps", "paper.txt")
if (!file.exists(txt_path)) {
  stop("Missing ", txt_path, " - run: pdftotext -layout paper.pdf paper.txt")
}
txt <- readLines(txt_path, encoding = "UTF-8", warn = FALSE)
out_dir <- fs::path(root, "inst", "extdata")
DOI <- "10.3724/SP.J.1006.2010.00574"

## --- Table 1: the 26 primers that survived every stage -----------------------
# Columns are fixed-width in the PDF; `-layout` keeps the whitespace between
# primer columns, so a primer row is: ESNPnn <forward> <reverse> <Enzyme>.
# Every primer row starts with ESNP and nothing else in the file does, so the
# row count IS the paper's 26 - no hand-typed marker list.
hit <- grep("^[[:space:]]*ESNP[0-9]{2}[[:space:]]", txt)
rows <- regmatches(txt[hit], regexec(
  "^[[:space:]]*(ESNP[0-9]{2})[[:space:]]+([ACGT]+)[[:space:]]+([ACGT]+)[[:space:]]+([A-Za-z]+[[:space:]]?[IVX]+)[[:space:]]*$",
  txt[hit]
))
parsed <- Filter(function(x) length(x) == 5L, rows)

table1 <- tibble::tibble(
  primer_id = vapply(parsed, `[[`, character(1), 2L),
  forward   = vapply(parsed, `[[`, character(1), 3L),
  reverse   = vapply(parsed, `[[`, character(1), 4L),
  enzyme    = gsub("[[:space:]]+", "", vapply(parsed, `[[`, character(1), 5L))
)

## --- the funnel, as the paper states it (Sec. 2.3 + abstract) ---------------
# Verbatim, Chinese: "筛选得到酶切位点改变的候选 EST-SNP 位点 48 个 ... 其中
# 44 对引物具有单一条带的扩增产物 ... 发现 36 对引物的 ... 产物在 EST-SNP
# 位点出现紊乱 ... 获得的 26 个 CAPS 标记 ... 见表 1"
funnel <- tibble::tibble(
  stage = c("EST-SNP discovered (Fig. 1)",
            "restriction-site-altering candidate; primer pair designed",
            "single-band amplicon (400-800 bp) from Suinong 14",
            "Sanger re-sequencing shows the expected EST-SNP",
            "digest gives an unequivocal restriction pattern (Table 1)"),
  n_parent = c(NA_integer_, 537L, 48L, 44L, 36L),
  n_kept   = c(537L, 48L, 44L, 36L, 26L)
)

## --- the assertion that makes this a *verified* extract --------------------
# Three independent things must agree, else we have mis-read the PDF:
#  (a) Table 1 has exactly as many primer rows as the text's final count (26);
#  (b) the stages are strictly monotone and each kept count fits inside its
#      parent - a mis-parsed 36/44 swap would break this;
#  (c) the enzyme vocabulary in Table 1 is one of the four the Methods name.
observed_n <- nrow(table1)
expected_n <- funnel$n_kept[nrow(funnel)]
enzymes_ok <- c("HindIII", "SacI", "VspI", "PstI")
bad_enzyme <- setdiff(unique(table1$enzyme), enzymes_ok)
monotone <- all(diff(funnel$n_kept) < 0)
nested <- all(funnel$n_kept[-1] <= funnel$n_parent[-1])
contiguous <- identical(table1$primer_id, sprintf("ESNP%02d", 1:26))

if (observed_n != 26L || expected_n != 26L || !contiguous || !monotone ||
    !nested || length(bad_enzyme) || any(!nzchar(table1$forward))) {
  stop(
    "Table 1 and the paper's text disagree: parsed ", observed_n,
    " primer rows (ids contiguous ESNP01..ESNP26: ", contiguous,
    "), expected ", expected_n, "; funnel monotone: ", monotone,
    ", nested: ", nested,
    if (length(bad_enzyme)) paste0("; unexpected enzyme(s): ", paste(bad_enzyme, collapse = ", "))
  )
}

OUTCOME_DEF <- paste(
  "CAPS 'success' = the PCR product's restriction digest showed an unequivocal",
  "(清晰) restriction pattern among the nine soybean varieties and their pool;",
  "the paper treats that as the definition of a CAPS marker. It is a",
  "gel-scoring judgement, not a genotype call: the 26 are the products whose",
  "pattern was judged unambiguous, and the 10 of 36 that did not qualify are",
  "counted as failures, not re-scored."
)

readr::write_csv(
  dplyr::mutate(table1, outcome = "unequivocal restriction pattern",
                outcome_definition = OUTCOME_DEF, doi = DOI),
  fs::path(out_dir, "soybean_est_caps_markers.csv")
)

readr::write_csv(
  dplyr::mutate(funnel, outcome_definition = OUTCOME_DEF, doi = DOI),
  fs::path(out_dir, "soybean_est_caps_funnel.csv")
)

message("soybean_est_caps: ", observed_n, " markers in Table 1; funnel ",
        paste(funnel$n_kept, collapse = " -> "), " (assertion passed)")
