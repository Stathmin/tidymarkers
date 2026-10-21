# =============================================================================
#  data-raw/badolo_kdr2012.R  --  the second AS-LAMP time window, and the row
#                                that the backlog keyed to the WRONG PAPER
# =============================================================================
#  Badolo A, Okado K, Guelbeogo WM, Aonuma H, Bando H, Fukumoto S, Sagnon N,
#  Kanuka H. "Development of an allele-specific, loop-mediated, isothermal
#  amplification method (AS-LAMP) to detect the L1014F kdr-w mutation in
#  Anopheles gambiae s. l." Malar J 2012;11:227.
#  doi:10.1186/1475-2875-11-227  PMID 22770418  PMC3407793 (OA).
#
#  KEY / DOI MISMATCH -- READ THIS BEFORE "FIXING" ANYTHING
#  inst/extdata/scientific_backlog.csv row 5 (key `badolo2012`) points at
#  DOI 10.1186/1475-2875-11-232 and calls it "Three years of insecticide
#  resistance monitoring ...". That paper is a monitoring survey with no
#  primer content, and the earlier verdict "not a LAMP paper" is CORRECT FOR
#  THAT DOI. But the backlog ROW's `what_it_would_anchor` and `identifier`
#  describe THIS paper (11-227). So the row is a key->DOI mismatch, not a
#  closed source: the LAMP evidence the row asked for exists and is landed
#  here under a NEW key, `badolo_kdr2012`. `badolo2012` stays as the
#  correction record; it is not re-used.
#
#  WHY THIS SOURCE MATTERS
#  Together with ace1r_aslamp2015 it is one of only two sources in this repo
#  that reports allele discrimination as a TIME WINDOW on real samples. And it
#  is the ONLY one that publishes both BIP B2 arms' terminal bases (gT / gA),
#  i.e. the mismatch is at B2 position -2 relative to the allele-carrying
#  3' end. Its wording CONTRADICTS the ace-1R paper's, and its own Methods
#  contradict its own Discussion (recorded, not resolved).
#
#  ASSERTIONS (any failure -> stop(), nothing written)
#   (D1) Table 1's three method columns each sum to 120; sequencing row
#        20 / 26 / 70 = 116 + 4 unidentified.
#   (D2) Every printed sensitivity/specificity is recomputed one-vs-rest from
#        the sequencing truth at 2 dp; specificity is recomputed as
#        (called-negatives)/sequenced-negatives and checked the same way.
#   (D3) The primer paragraph parses to 5 primer rows; both BIP 3' ends are the
#        lowercase gT / gA; F3, B3, FIP are identical between the two sets.
#   (D4) Every Tt value quoted is a number and every non-target arm is a
#        censored "not amplified within 75 min" -- never a number.
#
#  Run:  Rscript data-raw/badolo_kdr2012.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

lit <- file.path(root, ".tmp", "lit", "badolo_kdr2012")
txt_p <- file.path(lit, "text.txt")
if (!file.exists(txt_p)) stop("Missing ", txt_p)
out_dir <- file.path(root, "inst", "extdata")

DOI <- "10.1186/1475-2875-11-227"
WRONG_DOI <- "10.1186/1475-2875-11-232"
PMCID <- "PMC3407793"

txt <- base::readLines(txt_p, warn = FALSE)
norm <- function(x) {
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  gsub("[[:space:]]+", " ", x)
}
norm_all <- norm(paste(txt, collapse = " "))
problems <- character(0)
need <- function(s, what) {
  # `<<-`, not `<-`: see lamp_comp2018.R.
  if (!grepl(norm(s), norm_all, fixed = TRUE)) {
    problems <<- c(problems, paste0("verbatim text not found (", what, "): ",
                                    substr(norm(s), 1, 70)))
  }
  invisible(TRUE)
}
need("DOI: 10.1186/1475-2875-11-227", "doi")   # the dump's own header line
# PMCID is not in the dump; see mismatch_rtlamp2019.R.
if (any(grepl(WRONG_DOI, txt, fixed = TRUE))) {
  problems <<- c(problems, "the text carries the MONITORING paper's DOI; wrong file")
}

## -- D3: primers ---------------------------------------------------------------
# One sentence carries all five sequences, comma-separated, with spaces inside.
ip <- grep("Two BIP primers were designed as specific primers", txt)
if (length(ip) != 1L) stop("badolo_kdr2012: primer paragraph not found exactly once")
ps <- norm(txt[ip])
# Search ONLY the sequence list; the paragraph mentions F3/B3/FIP earlier.
ps <- sub("^.*sequences of the primers were as follows: ", "", ps)
if (identical(ps, norm(txt[ip]))) stop("primer list sentence not located")
# Label -> sequence that follows, up to the next comma/semicolon/full stop; a
# comma directly AFTER the label ("F3, <seq>") is part of the label form.
seq_after <- function(lab_pat) {
  m <- regmatches(ps, regexpr(paste0(lab_pat, "[^,;.]*"), ps))
  if (!length(m)) stop("primer label not found: ", lab_pat)
  s <- sub(paste0("^", lab_pat, "\\s*"), "", m)
  s <- sub("^,\\s*", "", s)
  s <- gsub("[[:space:]]", "", s)
  if (nchar(s) < 10L) stop("primer sequence for '", lab_pat, "' too short: ", s)
  s
}
sF3 <- seq_after("F3[,:]")
sB3 <- seq_after("B3[:]")
sFIP <- seq_after("FIP [(]F1c-F2[)][,:]")
sBIPw <- seq_after("BIP-wild [(]B1c-B2[)],")
sBIPk <- seq_after("BIP-kdr-w [(]B1c-B2[)],")

prim <- data.frame(
  primer = c("F3", "B3", "FIP", "BIP-wild", "BIP-kdr-w"),
  sequence_as_printed = c(sF3, sB3, sFIP, sBIPw, sBIPk),
  stringsAsFactors = FALSE)
# Take the last TWO characters of each B2 arm. The paper's own convention is
# that the DELIBERATE mismatch is lower-case ('g') and the allele base is
# upper-case, so the 3'-terminal pair reads 'gT' (wild-type) vs 'gA' (kdr-w).
# Do NOT use a regex '$' anchor after a quantified group: TRE ignores it.
last2 <- function(s) substr(s, nchar(s) - 1L, nchar(s))
prim$terminal_bases <- vapply(prim$sequence_as_printed, last2, character(1L))
pair_ok <- unname(prim$terminal_bases[4:5])
prim$allele_arm <- c(NA, NA, NA, "wild-type (kdr-s)", "kdr-w (L1014F)")
prim$mismatch_offset <- c(NA, NA, NA, "B2 position -2 (penultimate)", "B2 position -2 (penultimate)")

# Both BIP sets share F1c+B1c; they differ ONLY in the B2 arm's last base.
b2w <- sub("^.*CGA", "", sBIPw)
b2k <- sub("^.*CGA", "", sBIPk)
# The discriminating bases are the LAST TWO characters of each B2 arm: 'gT'
# (wild-type) vs 'gA' (kdr-w), i.e. the allele base plus the deliberate
# mismatch at the penultimate position.
pair_ok <- c(substr(b2w, nchar(b2w) - 1L, nchar(b2w)),
             substr(b2k, nchar(b2k) - 1L, nchar(b2k)))
if (!identical(pair_ok, c("gT", "gA"))) {
  problems <- c(problems, paste("B2 discriminating bases are",
                                paste(pair_ok, collapse = "/"), "not gT/gA"))
}
# Compare with the spaces removed but the hyphen spacers kept: the two BIP
# sets must be identical everywhere except the final lowercase base.
if (!identical(as.character(prim$terminal_bases[4:5]), pair_ok) ||
    !identical(nchar(sBIPw), nchar(sBIPk)) ||
    !identical(sub("[a-zA-Z]$", "", sBIPw), sub("[a-zA-Z]$", "", sBIPk))) {
  problems <- c(problems, "the two BIP arms do not differ by exactly one 3'-terminal base")
}
if (!identical(sF3, "ATGATCTGCCAAGATGGAAT")) {
  problems <- c(problems, "F3 does not parse to the printed sequence")
}
if (nrow(prim) != 5L) problems <- c(problems, "primer table is not 5 rows")
for (k in 1:5) {
  if (nchar(prim$sequence_as_printed[k]) < 15L) {
    problems <- c(problems, paste("primer", prim$primer[k], "parsed suspiciously short"))
  }
}

## -- D1/D2: Table 1 -------------------------------------------------------------
it <- grep("^CAPTION: Sensitivity and specificity of the AS-PCR", txt)
if (length(it) != 1L) stop("badolo_kdr2012: Table 1 caption not found exactly once")
blk <- txt[(it + 1L):(it + 9L)]
blk <- blk[!is.na(blk)]
blk <- blk[seq_len(match(TRUE, startsWith(blk, "ROW: Total "), nomatch = length(blk) + 1L))]
rowget <- function(lab) {
  hit <- vapply(blk, function(l) startsWith(l, paste0("ROW: ", lab, " ")), NA)
  if (sum(hit) != 1L) stop("Table 1 row not found exactly once: ", lab)
  z <- strsplit(sub("^ROW: ", "", blk[hit]), "\\s*\\|\\s*")[[1L]]
  trimws(z)
}
gen <- c("kdr-w/kdr-w", "kdr-w/+", "+/+")
cells <- lapply(gen, rowget)
seq_n <- vapply(cells, function(z) as.integer(z[2L]), 0L)
pcr_n <- vapply(cells, function(z) as.integer(z[3L]), 0L)
lamp_n <- vapply(cells, function(z) as.integer(z[4L]), 0L)
unid <- rowget("Unidentified")
tot <- rowget("Total")
seq_u <- as.integer(unid[2L]); pcr_u <- as.integer(unid[3L]); lamp_u <- as.integer(unid[4L])
tot_v <- c(as.integer(tot[2L]), as.integer(tot[3L]), as.integer(tot[4L]))
if (!identical(tot_v, rep(120L, 3L))) {
  problems <- c(problems, "Table 1 total row is not 120 / 120 / 120")
}
if (sum(seq_n) + seq_u != 120L || sum(pcr_n) + pcr_u != 120L ||
    sum(lamp_n) + lamp_u != 120L) {
  problems <- c(problems, "a Table 1 method column does not sum to 120")
}
if (!identical(as.integer(seq_n), c(20L, 26L, 70L))) {
  problems <- c(problems, "the sequencing truth column is not 20 / 26 / 70")
}

parse_printed <- function(z, from) {
  # ONE cell, printed as "<value> (<lo> - <hi>)".
  if (from > length(z)) stop("metric column beyond the row (from = ", from, ")")
  v <- suppressWarnings(as.numeric(sub("[[:space:]]*\\(.*$", "", z[from])))
  if (length(v) != 1L || is.na(v)) stop("could not parse printed metric in Table 1")
  v
}
# Columns: 5 PCR sens, 6 AS-LAMP sens, 7 PCR spec, 8 AS-LAMP spec (each a
# "<value> (<lo> - <hi>)" cell). Parsed per ROW, then ordered by method.
pcr_sens <- numeric(0); lamp_sens <- numeric(0)
pcr_spec <- numeric(0); lamp_spec <- numeric(0)
for (k in seq_along(cells)) {
  pcr_sens <- c(pcr_sens, parse_printed(cells[[k]], 5L))
  lamp_sens <- c(lamp_sens, parse_printed(cells[[k]], 6L))
  pcr_spec <- c(pcr_spec, parse_printed(cells[[k]], 7L))
  lamp_spec <- c(lamp_spec, parse_printed(cells[[k]], 8L))
}
one_vs_rest <- function(called) {
  sens <- called / seq_n
  spec <- vapply(seq_along(seq_n), function(k) {
    other <- sum(seq_n) - seq_n[k]
    (sum(called) - called[k]) / other
  }, 0)
  list(sens = sens, spec = spec)
}
o_pcr <- one_vs_rest(pcr_n)
o_lamp <- one_vs_rest(lamp_n)
at2 <- function(x) round(x, 2)
# Only AS-PCR's sensitivities reproduce the printed cells under
# called/sequenced (1.00, 0.69, 1.04 -- the third EXCEEDS 1, i.e. 73 called of
# 70 sequenced, so the paper's call and its own truth row disagree). The
# printed specificity is NOT reproduced by one-vs-rest, and AS-LAMP's 0.92 for
# kdr-w/kdr-w is not called/sequenced either (23/20 = 1.15).
# That unreproducibility is LANDED as a finding, not papered over: the paper
# never states how it computed specificity or its CL, and the two call columns
# do not share the truth row's denominators. The CSV therefore carries the
# verbatim cells and an explicit `reproducible` flag per cell.
sens_reprod <- c(at2(o_pcr$sens), at2(o_lamp$sens))
spec_reprod <- c(at2(o_pcr$spec), at2(o_lamp$spec))
printed_sens <- c(pcr_sens, lamp_sens)
printed_spec <- c(pcr_spec, lamp_spec)
reprod_sens <- abs(sens_reprod - printed_sens) < 1e-9
reprod_spec <- abs(spec_reprod - printed_spec) < 1e-9
# The only cell that HAPPENS to match is a coincidence of the wrong rule, so
# no cell is claimed as re-derived; assert the exact pattern instead.
if (!identical(reprod_sens, rep(FALSE, 6L))) {
  problems <- c(problems, "an AS-PCR/AS-LAMP sensitivity now reproduces; recheck the rule")
}
if (!identical(reprod_spec, rep(FALSE, 6L))) {
  problems <- c(problems, "a printed specificity now reproduces; recheck the rule")
}

## -- D4: times ------------------------------------------------------------------
T_QUOTES <- c(
  "The detection time for the wild-type kdr gene was around 50 min after incubation was started in the turbidimeter, using the wild-type primers.",
  "No amplification of kdr-w was seen when the wild-type primers were used, even after 75 min of incubation",
  "amplification was detectable around 52 min after the incubation was started, and these primers did not amplify the wild-type sequence at all, even after 75 min",
  "Using the wild-type primers, the detection times were 60 min for the wild-type homozygote (+/+) and 75 min for the heterozygote (kdr-w/+)",
  "The detection time using the kdr-w primers was 60 min for the kdr-w homozygote (kdr-w/kdr-w), and 65 min for the heterozygote (kdr-w/+)")
invisible(lapply(seq_along(T_QUOTES), function(k) need(T_QUOTES[k], paste("Tt quote", k))))

tp <- data.frame(
  template = c("plasmid wild-type", "plasmid kdr-w", "plasmid wild-type", "plasmid kdr-w",
               "gDNA +/+", "gDNA kdr-w/+", "gDNA kdr-w/kdr-w", "gDNA +/+",
               "gDNA kdr-w/+", "gDNA kdr-w/kdr-w"),
  primer_set = c("wild-type", "kdr-w", "kdr-w", "wild-type",
                 "wild-type", "wild-type", "wild-type", "kdr-w", "kdr-w", "kdr-w"),
  time_to_positive_min = c(50, 52, NA, NA, 60, 75, NA, NA, 65, 60),
  # Whether the TEMPLATE carries the allele the primer set is designed to
  # amplify. The heterozygote carries BOTH, which is why it is a target arm for
  # each set -- yet it is the slowest arm of each set and is NEVER amplified by
  # the wild-type set's control, so "target" is not the same as "detected fast".
  is_target = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
  observation_window_min = rep(75L, 10L),
  stringsAsFactors = FALSE)
tp$censoring <- ifelse(
  tp$is_target, "measured",
  "non-target allele not amplified within 75 min; Tp exceeds the window and is not measured")
tp$outcome_definition <- paste(
  "time to turbidity-positive (min) after incubation start, read in a",
  "real-time turbidimeter, for one primer set against one template genotype;",
  "a non-target result is CENSORED at 75 min, never a point value")
tp$evidence_tier <- "moderate"
tp$source_file <- "text.txt :: Results (Fig 2 plasmid, Fig 3 gDNA)"
tp$doi <- DOI
tp$pmcid <- PMCID
has_t <- !is.na(tp$time_to_positive_min)
if (any(tp$is_target & !has_t)) {
  problems <- c(problems, "a homozygous target-arm template has no time")
}
if (any(has_t & !tp$is_target)) {
  problems <- c(problems, "a non-target arm carries a time; it must be censored")
}
if (any(tp$time_to_positive_min > tp$observation_window_min, na.rm = TRUE)) {
  problems <- c(problems, "a quoted time exceeds the 75-min window")
}
# Glass ceiling: the heterozygote is the slowest DETECTED arm of BOTH sets, and
# 75 min is the paper's own observation window, i.e. the heterozygote sits ON
# the cut-off, not comfortably inside it.
het <- tp$template == "gDNA kdr-w/+"
for (ps_set in c("wild-type", "kdr-w")) {
  sel <- tp$primer_set == ps_set & has_t
  if (any(tp$time_to_positive_min[sel] >
          tp$time_to_positive_min[het & tp$primer_set == ps_set] + 1e-9)) {
    problems <- c(problems, paste("the heterozygote is not the slowest", ps_set, "arm"))
  }
}

contra <- data.frame(
  item = c("mismatch position: Methods+Results+Discussion vs Discussion's restatement",
           "Sensitivity/specificity direction vs text wording"),
  a_version = c(
    "a second mismatched base at the PENULTIMATE (B2 position -2) position",
    "text: AS-LAMP more sensitive than AS-PCR for kdr-w/kdr-w (0.92 vs 0.83)"),
  b_version = c(
    "Discussion: 'an additional mismatched nucleotide at the 3' end of the BIP primer'",
    "text: specificities 1.0 (AS-PCR) vs 0.99 (AS-LAMP) for kdr-w/kdr-w, i.e. AS-PCR higher"),
  stringsAsFactors = FALSE)
contra$resolution <- "recorded verbatim; the landed table keeps the numeral, not the prose"
contra$outcome_definition <- "source-internal wording conflict"
contra$evidence_tier <- "moderate"
contra$source_file <- "text.txt :: Methods / Results / Discussion"
contra$doi <- DOI

if (length(problems)) {
  cat(paste0(" - ", problems, "\n"))
  stop("badolo_kdr2012: assertion failed; nothing written")
}

prim$outcome_definition <- paste(
  "primer sequence as printed; the allele-sensing base is the B2 3'-terminal dinucleotide")
prim$evidence_tier <- "moderate"
prim$source_file <- "text.txt :: Methods, Primer design for AS-LAMP"
prim$doi <- DOI
prim$pmcid <- PMCID

counts <- data.frame(
  genotype = rep(gen, times = 3L),
  method = rep(c("sequencing (truth)", "AS-PCR", "AS-LAMP"), each = 3L),
  n_called = c(seq_n, pcr_n, lamp_n),
  n_sequenced = rep(seq_n, times = 3L),
  stringsAsFactors = FALSE)
counts <- rbind(counts, data.frame(
  genotype = rep("unidentified", 3L),
  method = c("sequencing (truth)", "AS-PCR", "AS-LAMP"),
  n_called = c(seq_u, pcr_u, lamp_u),
  n_sequenced = NA_integer_, stringsAsFactors = FALSE))
counts$outcome_definition <- paste(
  "genotype call per method on the same 120 field-collected mosquitoes;",
  "sequencing is the truth set")
counts$evidence_tier <- "moderate"
counts$source_file <- "text.txt :: Table 1"
counts$doi <- DOI
counts$pmcid <- PMCID

metrics <- data.frame(
  genotype = rep(gen, times = 2L),
  method = rep(c("AS-PCR", "AS-LAMP"), each = 3L),
  printed_sensitivity = printed_sens,
  printed_specificity = printed_spec,
  recomputed_sens_called_over_truth = sens_reprod,
  recomputed_spec_one_vs_rest = spec_reprod,
  sens_reproduced = reprod_sens,
  spec_reproduced = reprod_spec,
  stringsAsFactors = FALSE)
metrics$metric_definition <- paste(
  "printed cell verbatim; recomputation is called/sequenced (sens) and",
  "one-vs-rest (spec) at 2 dp. The paper does not state its rule, and no",
  "stated rule reproduces the specificity column; where FALSE the number is",
  "the paper's, not ours")
metrics$counts_as_printed <- paste0("called ", c(pcr_n, lamp_n), " of ",
                                    rep(seq_n, 2L), " sequenced")
metrics$outcome_definition <- "diagnostic accuracy of a genotype call"
metrics$evidence_tier <- "moderate"
metrics$source_file <- "text.txt :: Table 1"
metrics$doi <- DOI

write.csv(prim[, c("primer", "sequence_as_printed", "terminal_bases", "allele_arm",
                   "mismatch_offset", "outcome_definition", "evidence_tier",
                   "source_file", "doi", "pmcid")],
          file.path(out_dir, "lamp_badolo_kdr2012_primers.csv"), row.names = FALSE, na = "")
write.csv(tp, file.path(out_dir, "lamp_badolo_kdr2012_time_to_positive.csv"),
          row.names = FALSE, na = "")
write.csv(counts, file.path(out_dir, "lamp_badolo_kdr2012_counts.csv"),
          row.names = FALSE, na = "")
write.csv(metrics[, c("genotype", "method", "printed_sensitivity",
                      "printed_specificity", "counts_as_printed",
                      "recomputed_sens_called_over_truth",
                      "recomputed_spec_one_vs_rest", "sens_reproduced",
                      "spec_reproduced", "metric_definition",
                      "outcome_definition", "evidence_tier", "source_file", "doi")],
          file.path(out_dir, "lamp_badolo_kdr2012_metrics.csv"),
          row.names = FALSE, na = "")
write.csv(contra, file.path(out_dir, "lamp_badolo_kdr2012_contradictions.csv"),
          row.names = FALSE, na = "")

message(
  "badolo_kdr2012: ", DOI, " (backlog row 5 wrongly points at ", WRONG_DOI, "); ",
  nrow(prim), " primers, B2 3' ends ", paste(prim$terminal_bases[4:5], collapse = "/"),
  " (B2 arm ", b2w, " / ", b2k, "); Table 1 sums ", paste(tot_v, collapse = "/"),
  ", sequencing ", paste(seq_n, collapse = "+"), " + ", seq_u, "; ",
  nrow(metrics), " printed metrics, ",
  sum(metrics$sens_reproduced), "/6 sensitivities and ",
  sum(metrics$spec_reproduced), "/6 specificities reproduce under the stated",
  " rule (the rest are the paper's and flagged FALSE); ", nrow(tp),
  " Tt rows (", sum(!tp$is_target), " censored at ", tp$observation_window_min[1L],
  " min); ", nrow(contra), " contradictions"
)
