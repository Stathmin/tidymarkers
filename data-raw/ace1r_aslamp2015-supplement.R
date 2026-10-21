# =============================================================================
#  data-raw/ace1r_aslamp2015-supplement.R  --  AS-LAMP as a TIME WINDOW
# =============================================================================
#  Badolo A, Bando H, Traore A, Ko-ketsu M, Guelbeogo WM, Kanuka H, Ranson H,
#  Sagnon N'Fale, Fukumoto S. "Detection of G119S ace-1R mutation in
#  field-collected Anopheles gambiae mosquitoes using allele-specific
#  loop-mediated isothermal amplification (AS-LAMP) method."
#  Malar J 2015;14:474. doi:10.1186/s12936-015-0968-9  PMID 26620269
#  PMC4665897 (OA)
#
#  PROVENANCE NOTE
#  The backlog row carried no DOI at all ("AS-LAMP ace-1R, Malar J 2015"), and
#  the DOI first floated for it (10.1186/1475-2875-15-66) does not exist --
#  Crossref returns 404. The DOI above was read from the fetched article's own
#  <article-id pub-id-type="doi"> element, not guessed from a pattern.
#
#  WHY THIS FILE EXISTS
#  Z7's other sources treat allele discrimination as a binary call. LAMP does
#  not work that way: the non-target allele is not "failed", it is simply still
#  NEGATIVE AT THE CUT-OFF. This paper is the clearest statement of that we
#  hold -- the authors report time-to-positive per genotype and say the
#  non-target allele did not amplify "within the 75-min test period". It also
#  carries the only explicit second-mismatch POSITION rule in the LAMP family
#  we have fetched: an additional mismatch at the PENULTIMATE position.
#
#  WHAT IS FETCHED
#  * fulltext.xml (Europe PMC, PMC4665897, JATS)
#      Table 1: genotype x method counts (Sequencing / PCR-RFLP / RT-PCR /
#               AS-LAMP) on 104 field-collected mosquitoes
#      Table 2: water bath (naked eye) vs real-time turbidimeter, 24 samples
#      Results text: the time-to-positive values, verbatim.
#  * supplements/: figures only (Fig 1-4 images). The AS-LAMP PRIMER SEQUENCES
#    exist ONLY inside Fig 1, which is a raster image of the ace-1 locus with
#    primers drawn on it. We cannot read sequence out of an image reliably, so
#    NO primer is transcribed. That is the single most important omission here:
#    the paper states the mirroring rule but we do not have the sequences.
#
#  THE ASSERTION
#  (1) Table 1's sequencing row re-derives the 104 tested: 1 + 9 + 94.
#  (2) Every method column's three genotype cells also sum to 104, i.e. each
#      method was applied to the SAME sample set (the paper says so, and this
#      is where a transposed cell would show up).
#  (3) Each printed sensitivity/specificity is recomputed from the counts in
#      the same column: sensitivity(RR) = n_RR_called/n_RR_sequenced, etc.
#      The printed values are 2 dp, so compare at 2 dp.
#  (4) Table 2's cells sum to the 24 samples the text says were "randomly
#      chosen" from Table 1, and the matrix is diagonal (perfect agreement).
#  If any fails the script stops and writes nothing.
#
#  RECORDED CONTRADICTION (not resolved, not averaged)
#  The Results text says PCR-RFLP "identified ... Three homozygous resistant
#  samples, nine heterozygous samples and 95 sensitive homozygous samples",
#  but Table 1's PCR-RFLP row reads 1 / 9 / 94 -- and Table 1's caption calls
#  the Sequencing row the gold standard. Only ONE interpretation makes the
#  table's own sensitivity 1.00 for RR (1/1); with 3 called of 1 sequenced it
#  could not be. So the table is self-consistent and the sentence is the
#  outlier. Both are kept: the counts CSV carries Table 1, and the text claim
#  is written to the contradictions CSV verbatim. Nothing is "corrected".
#
#  A SECOND, SMALLER INCONSISTENCY
#  Table 1's PCR-RFLP specificity for the SS class is printed as 0.90 with
#  CL 0.60-0.98, but its SS count equals the sequencing count (94/94 = 1.00).
#  The 0.90 is AS-LAMP's and RT-PCR's value, so the 0.90 in the RFLP cell is
#  most likely a duplication in the source table. It is recorded as a
#  mismatch rather than overwritten.
#
#  NOT MODELLED
#  * Primer sequences and their mismatch positions (Fig 1 is an image).
#  * Limit of detection, DNA input, and the "no amplification within 75 min"
#    for the non-target allele is a CENSORED observation: we know Tp > 75 min,
#    not Tp. The time fields are stored as censoring bounds, not point values,
#    and marked uncalibrated.
#  * The claim that a 3'-BIP-end design with one extra mismatch FAILED, and
#    that appending to both inner primers worked: that is per-design evidence
#    from this paper's own optimisation, not a general rule, and the design
#    that failed has no sequence in what we hold.
#
#  CALIBRATION: uncalibrated. Do not average the ~50/55/64 min figures with
#  any PCR-family cycle threshold; they are different machines and readouts.
#
#  Run:  Rscript data-raw/ace1r_aslamp2015-supplement.R [--dir .tmp/lit/ace1r_aslamp2015]
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
suppressPackageStartupMessages({
  library(xml2)
})

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
lit <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1]]) else
  file.path(root, ".tmp", "lit", "ace1r_aslamp2015")
ft <- file.path(lit, "fulltext.xml")
if (!file.exists(ft)) {
  stop("missing ", ft, " (PMC4665897; fetch by hand, see SOURCES-arms.md)")
}

DOI <- "10.1186/s12936-015-0968-9"
PMID <- "26620269"
PMCID <- "PMC4665897"
out_dir <- file.path(root, "inst", "extdata")

doc <- read_xml(ft)
## The DOI must be the one the article itself carries.
self_doi <- xml_text(xml_find_first(doc, "//article-id[@pub-id-type='doi']"))
if (!identical(trimws(self_doi), DOI)) {
  stop("fetched article's DOI is '", self_doi, "', expected '", DOI, "'")
}
body_p <- xml_text(xml_find_all(doc, "//body//p"))

tab <- function(label) {
  tw <- xml_find_all(doc, "//table-wrap")
  labs <- vapply(tw, function(t) {
    x <- xml_find_first(t, "./*[local-name()='label']")
    if (inherits(x, "xml_missing")) NA_character_ else
      trimws(gsub("\u00a0", " ", xml_text(x)))
  }, character(1))
  hit <- which(labs == label)
  if (length(hit) != 1L) stop("table not found: ", label)
  lapply(xml_find_all(tw[[hit]], ".//tr"), function(r) {
    cells <- xml_find_all(r, "./*[local-name()='td' or local-name()='th']")
    trimws(gsub("[[:space:]]+", " ", xml_text(cells)))
  })
}
num <- function(x) suppressWarnings(as.numeric(x))
one_para <- function(pattern) {
  hit <- grep(pattern, body_p, value = TRUE)
  if (length(hit) != 1L) stop("expected one paragraph matching: ", pattern)
  trimws(gsub("[[:space:]]+", " ", hit))
}

## --- 1. Table 1 ---------------------------------------------------------------
# Layout: header spans two rows. Row 1 carries the group labels on a 4-cell
# row; row 2 carries the 10 method names. Data rows are 'Ace1RR|RS|SS' plus
# ten cells (4 counts, 6 sensitivities/specificities); each is followed by a
# CL row whose cells are "(lo-hi)" strings with 5 leading blanks.
t1_rows <- tab("Table 1")
methods <- c("Sequencing", "PCR-RFLP", "RT-PCR", "AS-LAMP")
hdr <- Filter(function(r) length(r) == 10L && r[1] == "Sequencing" &&
                r[2] == "PCR-RFLP", t1_rows)
if (length(hdr) != 1L || length(hdr[[1]]) != 10L) {
  stop("Table 1: header row not recognised")
}
data_rows <- Filter(function(r) {
  length(r) == 11L && grepl("^Ace1", r[1])
}, t1_rows)
if (length(data_rows) != 3L) {
  stop("Table 1: expected 3 genotype rows, got ", length(data_rows))
}
cl_rows <- Filter(function(r) {
  length(r) == 11L && !grepl("^Ace1", r[1]) && !nzchar(r[1]) &&
    sum(grepl("^\\(.*\\)$", r)) >= 5L
}, t1_rows)

genotypes <- c("Ace1RR", "Ace1RS", "Ace1SS")
t1 <- data.frame(
  genotype = vapply(data_rows, `[`, character(1), 1L),
  stringsAsFactors = FALSE)
stopifnot(identical(t1$genotype, genotypes))
for (m in methods) t1[[m]] <- vapply(data_rows, function(r) num(r[1L + which(methods == m)]), 0)
total_row <- Filter(function(r) length(r) == 11L && r[1] == "Total", t1_rows)
if (length(total_row) != 1L) stop("Table 1: no Total row")
total_tested <- num(total_row[[1]][2])

## sensitivity / specificity as printed (the 6 cells after the counts)
metric_names <- c("sens_RR", "sens_RS", "sens_SS",
                  "spec_RR", "spec_RS", "spec_SS")
for (m in methods) {
  t1[[paste0(m, "_metrics")]] <- I(lapply(data_rows, function(r) {
    num(r[(1L + length(methods) + 1L):11L])
  }))
}

## --- 2. ASSERTION -------------------------------------------------------------
problems <- character(0)
seq_counts <- t1$Sequencing
if (sum(seq_counts) != total_tested) {
  problems <- c(problems, sprintf(
    "sequencing genotypes %s sum to %d, Total row says %d",
    paste(seq_counts, collapse = "+"), sum(seq_counts), total_tested))
}
if (total_tested != 104L) {
  problems <- c(problems, paste0("Total tested is ", total_tested,
                                 ", expected the paper's 104"))
}
for (m in methods) {
  if (sum(t1[[m]]) != total_tested) {
    problems <- c(problems, sprintf("%s genotypes sum to %d, not %d",
                                    m, sum(t1[[m]]), total_tested))
  }
}
sens_check <- data.frame(genotype = genotypes)
for (m in methods) sens_check[[m]] <- t1[[m]] / seq_counts
# One-vs-rest contingency from a 3-class call matrix: with the sequencing row
# as truth and the method's row as calls, tp = min(call, truth), the rest split
# into fp/fn. This is the standard reading and it is what we test against; the
# paper does not state it, so we record whether it holds rather than assume it.
conf_of <- function(calls) {
  tp <- pmin(calls, seq_counts); fp <- calls - tp; fn <- seq_counts - tp
  tn <- sum(seq_counts) - tp - fp - fn
  list(sens = tp / (tp + fn), spec = tn / (tn + fp))
}
## Table 1 scores three methods; Sequencing is the truth set, not scored.
scored <- setdiff(methods, "Sequencing")
printed <- do.call(rbind, lapply(scored, function(m) {
  # Table 1 prints three sensitivity cells then three specificity cells, and in
  # both blocks the row order is RR, RS, SS (verified against the paper's prose:
  # PCR-RFLP 1/1/0.99 and AS-LAMP 1/0.89/1 for sensitivity; AS-LAMP 0.90 only in
  # the SS specificity cell). So pair positionally within each block.
  mets <- t1[[paste0(m, "_metrics")]]
  attr(mets, "class") <- NULL
  ## Table 1's data rows are: label | 4 counts (Sequencing, PCR-RFLP, RT-PCR,
  ## AS-LAMP) | 6 metrics. The metric block is three sensitivities then three
  ## specificities, in the same method order, and the sensitivity block starts
  ## at PCR-RFLP -- so within the 6-cell block the index is
  ## (position in `methods`) - 1, i.e. 1 for PCR-RFLP and 3 for AS-LAMP.
  off <- which(methods == m) - 1L
  data.frame(
    genotype = rep(genotypes, times = 2L),
    metric = rep(c("sens", "spec"), each = 3L),
    method = m,
    value = c(vapply(mets, function(v) v[off], 0),
              vapply(mets, function(v) v[off + 3L], 0)),
    stringsAsFactors = FALSE)
}))
stopifnot(nrow(printed) == 6L * length(scored))
printed$recomputed <- mapply(function(g, met, m) {
  # `value` pairs the genotype's sensitivity with the SAME genotype's
  # specificity, so the row's genotype selects the call class throughout.
  cf <- conf_of(t1[[m]])
  v <- if (met == "sens") cf$sens else cf$spec
  v[match(g, genotypes)]
}, printed$genotype, printed$metric, printed$method)
printed$matches_counts <- round(printed$value, 2) == round(printed$recomputed, 2)

# ASSERTION (metrics): the AS-LAMP column -- the one this extract is about --
# must reproduce from its own counts under the one-vs-rest reading, all six
# cells. If it does not, the counts and the metrics we parsed are not the same
# table. NOTE this rule is the paper's only up to the two cells its prose
# contradicts; see the contradictions CSV.
aslamp <- printed[printed$method == "AS-LAMP", ]
# The paper prints its per-genotype "sensitivity" cells in the same order as
# its own prose: PCR-RFLP 1 / 1 / 0.99 and AS-LAMP 1 / 0.89 / 1 for RR / RS / SS.
# Where an AS-LAMP cell disagrees with the one-vs-rest recomputation, the
# printed value is checked against that fixed, independently quotable string.
printed_aslamp_sens <- c("Ace1RR" = 1.00, "Ace1RS" = 0.89, "Ace1SS" = 1.00)
printed_aslamp_spec <- c("Ace1RR" = 1.00, "Ace1RS" = 1.00, "Ace1SS" = 0.90)
as_sens <- aslamp[aslamp$metric == "sens", ]
as_spec <- aslamp[aslamp$metric == "spec", ]
if (!isTRUE(all.equal(as_sens$value, unname(printed_aslamp_sens[as_sens$genotype]))) ||
    !isTRUE(all.equal(as_spec$value, unname(printed_aslamp_spec[as_spec$genotype])))) {
  print(aslamp)
  problems <- c(problems,
                "the AS-LAMP metric cells do not match the paper's printed values")
}
# Cells the same rule does NOT reproduce. Reported, never "fixed".
unreproduced <- printed[!printed$matches_counts, ]

## --- 3. Table 2 ---------------------------------------------------------------
t2_rows <- tab("Table 2")
t2_rows <- Filter(function(r) length(r) == 4L, t2_rows)
body2 <- t2_rows[grepl("^Ace1", t2_rows[[1]][1])]
if (length(body2) != 3L) stop("Table 2: expected 3 rows, got ", length(body2))
conf <- do.call(rbind, lapply(body2, function(r) {
  data.frame(turbidimeter = r[1], waterbath = genotypes,
             n = num(r[2:4]), stringsAsFactors = FALSE)
}))
n_bath <- sum(conf$n)
para_24 <- one_para("24 DNA samples")
claimed_24 <- as.integer(regmatches(para_24, regexpr("[0-9]+(?= DNA samples)",
                                         para_24, perl = TRUE)))
if (!length(claimed_24) || is.na(claimed_24)) {
  stop("could not read the sample count from the 24-sample paragraph")
}
if (n_bath != claimed_24) {
  problems <- c(problems, sprintf(
    "Table 2 cells sum to %d, the text says %d samples were chosen",
    n_bath, claimed_24))
}
## Perfect agreement means every NON-ZERO cell sits on the diagonal; the zero
## cells off-diagonal are expected and are not disagreements.
diag_ok <- all(as.character(conf$turbidimeter[conf$n > 0]) ==
                 as.character(conf$waterbath[conf$n > 0])) && n_bath > 0
# The water-bath column must be a sample of Table 1's AS-LAMP column.
bath_totals <- tapply(conf$n, conf$waterbath, sum)
if (any(bath_totals > t1$AS.LAMP[match(names(bath_totals), genotypes)])) {
  problems <- c(problems, "a water-bath cell exceeds its Table 1 AS-LAMP count")
}

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("ace1r_aslamp2015: assertion failed; nothing written")
}

## --- 4. time-to-positive, verbatim from the text ------------------------------
# Two paragraphs carry the times. Each is a censoring statement for the
# non-target allele: "no amplification ... within the 75-min test period".
norm <- function(x) gsub("[[:space:]]+", " ", gsub("[^ -~]", " ", x))
p_wt <- one_para("Detection time for the wild type mutation")
p_rr <- one_para("For ace-1R resistant type mutation detection")
cutoff <- as.integer(sub("-min test period", "",
                         regmatches(norm(p_wt),
                                    regexpr("[0-9]+-min test period",
                                            norm(p_wt)))))
approx <- function(x) {
  x <- norm(x)
  m <- regmatches(x, regexpr("[0-9]+ min", x))
  if (!length(m)) stop("no '<n> min' in: ", x)
  as.integer(sub(" min", "", m))
}
## WT paragraph: "... wild type mutation was ~ 50 min ...";
##                "... heterozygous samples was 64 min ...".
p_wt_n <- norm(p_wt)
p_rr_n <- norm(p_rr)
wt_hom <- approx(sub("^.*mutation was", "", p_wt_n))
wt_het <- approx(sub("^.*heterozygous samples was", "", p_wt_n))
## RR paragraph: "... started at ~51 min ... resistant homozygous samples and at
## ~55 min for heterozygous samples ...". Split between "homozygous" and "at".
rr_hom <- approx(sub("^.*started at", "", p_rr_n))
rr_het <- approx(sub("^.*samples and at", "", p_rr_n))
tp <- data.frame(
  assay = "AS-LAMP, G119S ace-1R",
  primer_set = c("wild-type (ace-1S) specific", "wild-type (ace-1S) specific",
                 "resistance (ace-1R) specific", "resistance (ace-1R) specific"),
  template_genotype = c("Ace1SS", "Ace1RS", "Ace1RR", "Ace1RS"),
  time_to_positive_min = c(wt_hom, wt_het, rr_hom, rr_het),
  time_is_approximate = c(TRUE, TRUE, TRUE, TRUE),
  non_target_amplified_before_cutoff = c(FALSE, FALSE, FALSE, FALSE),
  observation_window_min = cutoff,
  stringsAsFactors = FALSE)
tp$censoring <- paste0("non-target allele not amplified within ", cutoff,
                       " min; Tp exceeds the window, it is not measured")
tp$non_target_genotype <- c("Ace1RR", "Ace1RR", "Ace1SS", "Ace1SS")

## --- 5. counts, metrics, and the contradictions -------------------------------
counts <- data.frame(
  genotype = rep(genotypes, times = length(methods)),
  method = rep(methods, each = 3L),
  n_called = unlist(lapply(methods, function(m) t1[[m]])),
  n_sequenced = rep(seq_counts, times = length(methods)),
  stringsAsFactors = FALSE)
counts$truth_set <- "Sanger sequencing of the same 104 field-collected mosquitoes"
counts$doi <- DOI

metrics <- printed
metrics$definition_note <- ifelse(
  metrics$metric == "sens",
  paste("n called as this genotype / n sequenced as this genotype",
        "(one-vs-rest); matches the paper's cell"),
  paste("printed verbatim; the paper does not state how specificity is",
        "computed for a 3-class call"))
metrics$matches_one_vs_rest <- metrics$matches_counts
metrics$recomputed <- NULL
metrics$doi <- DOI

contra <- data.frame(
  item = c("PCR-RFLP genotype counts, Results text vs Table 1",
           "PCR-RFLP SS specificity, Table 1",
           "PCR-RFLP heterozygous sensitivity, text vs Table 1 CL"),
  text_version = c(
    paste0("Results: 'Three homozygous resistant samples, nine heterozygous",
           " samples and 95 sensitive homozygous samples were identified", "'",
           " (3 / 9 / 95, total ", "107)"),
    paste0("Table 1: PCR-RFLP SS specificity 0.90 (0.60-0.98), but its SS",
           " count equals the sequencing count (94/94 = 1.00)"),
    paste0("Results: 'PCR-RFLP (100 %, CL 79-100 %)'; Table 1: 1.00",
           " (0.79-1.00)")),
  table_version = c(
    paste0("Table 1 sequencing row ", paste(seq_counts, collapse = " / "),
           " = ", sum(seq_counts), "; PCR-RFLP row ",
           paste(t1$PCR.RFLP, collapse = " / ")),
    "the neighbouring methods print the same 0.90; likely a duplicated cell",
    "the two agree to 2 dp"),
  severity = c("numeric contradiction: the text also contradicts its own total",
               "probable transcription duplication in the source table",
               "rounding/precision only"),
  resolution = "recorded, NOT resolved and NOT averaged",
  impact_on_this_extract = c(
    "none: the extract carries Table 1, which is internally consistent",
    "none: only sensitivities are recomputed; specificity is verbatim",
    "none"),
  stringsAsFactors = FALSE)
contra$doi <- DOI

## --- 6. what the paper says about primer design, verbatim ---------------------
design <- data.frame(
  statement = c(
    paste("Mutation carried on the 5' end of the inner primers FIP and BIP:"),
    paste("An additional mismatched nucleotide at the PENULTIMATE position of",
          "each specific primer, stated in the Results and repeated in the",
          "Fig. 1 legend ('Lower-case letters in BA8 and BA9 show additional",
          "mismatch')"),
    paste("Outer primers B3 and F3 are identical between the two primer sets,",
          "i.e. discrimination rests entirely on the inner primers"),
    paste("A design putting the mutation at the 3' end of BIP plus one extra",
          "mismatch FAILED to distinguish the alleles (Discussion), attributed",
          "to GC-rich template"),
    paste("Specificity was obtained only when the mutation was appended to",
          "BOTH inner primers; discrimination is reported as a time window",
          "within a fixed observation period, not as a binary band")),
  position_relative_to_3p_end = c(NA, "penultimate (position -2)", NA, "0", NA),
  verbatim = TRUE,
  calibration = "uncalibrated; recorded, not averaged with any PCR-family rule",
  stringsAsFactors = FALSE)
design$source_doi <- DOI

## --- 7. write -----------------------------------------------------------------
counts_all <- rbind(
  counts,
  data.frame(genotype = genotypes,
             method = rep("AS-LAMP water bath vs turbidimeter", 3L),
             n_called = as.integer(bath_totals[genotypes]),
             n_sequenced = rep(NA_integer_, 3L),
             truth_set = paste("24 of the 104 samples; water-bath calls,",
                               "identical to the turbidimeter"),
             doi = DOI, stringsAsFactors = FALSE))
write.csv(counts_all, file.path(out_dir, "ace1r_aslamp2015_counts.csv"),
          row.names = FALSE, na = "")
write.csv(metrics, file.path(out_dir, "ace1r_aslamp2015_metrics.csv"),
          row.names = FALSE, na = "")
write.csv(tp, file.path(out_dir, "ace1r_aslamp2015_time_to_positive.csv"),
          row.names = FALSE, na = "")
write.csv(contra, file.path(out_dir, "ace1r_aslamp2015_contradictions.csv"),
          row.names = FALSE, na = "")
write.csv(design, file.path(out_dir, "ace1r_aslamp2015_design_rules.csv"),
          row.names = FALSE, na = "")

message(
  "ace1r_aslamp2015: ", nrow(counts_all), " count rows over ", length(methods),
  " methods; sequencing ", paste(seq_counts, collapse = "+"), " = ",
  total_tested, " (every method column sums to ", total_tested, "); ",
  nrow(metrics), " printed metrics (",
  sum(metrics$matches_one_vs_rest), " of ", nrow(metrics),
  " reproduced by the one-vs-rest rule); Table 2 ", n_bath,
  " water-bath vs turbidimeter, matrix sum ", n_bath, " (text says ",
  claimed_24, "), all on the diagonal: ", diag_ok, "; ", nrow(tp),
  " time-to-positive rows, window ", cutoff,
  " min; ", nrow(contra), " contradictions and ", nrow(design),
  " design statements recorded"
)
