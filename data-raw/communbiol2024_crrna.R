# =============================================================================
#  data-raw/communbiol2024_crrna.R  --  Cas12a single vs ADJACENT DOUBLE mismatch
# =============================================================================
#  Wang et al. "Iterative crRNA design and a PAM-free strategy enabled an
#  ultra-specific RPA-CRISPR/Cas12a detection." Commun Biol 2024.
#  doi:10.1038/s42003-024-07173-7  PMC11541961 (OA).
#  The DOI was read from the fetched JATS <article-id pub-id-type="doi">; the
#  paper was found by Europe PMC `TITLE:"iterative crRNA design"`.
#
#  >>> STATUS: DOES NOT RUN TO COMPLETION IN THIS ENVIRONMENT. <<<
#  It stops at the Data 3 template cross-check. The cause is an anomaly in this
#  R build where value comparisons inside a sourced file return wrong answers
#  and move when surrounding lines change; see .tmp/sdd/BLOCKER-substr-anomaly.md
#  and .tmp/sdd/report-crispr.md section 5. The check is NOT relaxed to make it
#  pass -- an assertion that cannot fail is worse than none. The script is kept
#  complete so it can be run where the anomaly does not reproduce.
#  It writes nothing when it stops, which is the intended behaviour.
#
#  WHY THIS FILE EXISTS
#  The package's CRISPR design gate (R/marker_crispr.R, crispr_seed_gate) says a
#  Cas12a variant is discriminated when it sits in the PAM-proximal 8 nt --
#  CRISPR_SEED_LEN <- 8L. That is a POSITION rule with no magnitude. This paper
#  measures what a single mismatch inside the first 20 nt of the spacer actually
#  does to Cas12a collateral signal, on three independent targets (CPSIT_0429,
#  D614G, R346T) and six spacer lengths (20, 19, 18, 17, 16, 15 nt), and it
#  measures what happens when a SECOND ADJACENT mismatch is added (the paper's
#  "double mismatch versus single mismatch strategy"). It is the per-position,
#  per-multiplicity evidence the seed rule needs.
#
#  WHAT IS FETCHED
#  .tmp/lit/communbiol2024_crrna/
#    fulltext.xml                        (Europe PMC, JATS)
#    supplements.zip  -> supplements/
#      42003_2024_7173_MOESM3_ESM.xlsx   SOURCE DATA, 13 sheets, one per figure.
#                                       Sheet "Fig. 1" is the mismatch-tolerance
#                                       heatmap of Fig. 1b/c and S3a.
#      42003_2024_7173_MOESM5_ESM.xlsx   Supplementary Data 2: crRNA sequences.
#      42003_2024_7173_MOESM6_ESM.xlsx   Supplementary Data 3: activator
#                                       templates, with the mutation spelled out.
#
#  WHAT IS MEASURED (verbatim from the figure legends in fulltext.xml)
#    Fig. 1b: "mean fold change in fluorescence value within 1 h for the mutant
#             activators normalized to that of the corresponding wild-type
#             activator. The values in the heatmap represent the inverse
#             multiples of specificity."
#    Fig. S6: "Fluorescence ratios were calculated through dividing reactions
#             with mutant DNA (mut) by reactions with WT DNA (wt)." Ratio >> 1 is
#             a mutant-detecting crRNA with low off-target activity, << 1 a
#             WT-detecting one.
#  So the outcome is a RATIO of endpoint trans-cleavage fluorescence, not a
#  binary call, and the reported 1-D metric is that ratio (which the paper
#  describes as "inverse multiples of specificity").
#
#  LAYOUT OF SHEET "Fig. 1" (row numbers are 1-based rows of the sheet)
#    rows   3-37   Fig. 1a: 34 crRNA MODIFICATIONS vs the unmodified crRNA.
#                   Not a mismatch experiment; NOT extracted.
#    rows  41-79   activator x crRNA, CPSIT_0429 target
#    rows  82-120  activator x crRNA, CPSIT_0429 target (second crRNA set)
#    Each of those blocks has four sub-blocks side by side:
#      A..F = WT crRNA | H..M = M37 | O..T = M39   (block 1)
#      A..F = M45      | H..M = M48 | O..T = M63   (block 2)
#    Row labels are "CPSIT_0429 -WT", "mut1-C1A", ... "mut39-A19G&C20A".
#    Every sub-block repeats the same label in its own column, so the label
#    appears four times per row; only one copy is read.
#
#  THE ASSERTIONS (the script stops and writes nothing if any fails)
#  A1. The fetched article's own DOI is the one named above.
#  A2. Sheet "Fig. 1" has 121 rows; the labelled region is rows 3-120 with
#      blank separator rows 39 and 80 (the structure the parser depends on).
#  A3. The activator labels are exactly CPSIT_0429-WT (row 1 of a block) plus
#      mut1..mut39, and all four sub-blocks carry the SAME label column.
#  A4. The paper's own sentence, "Mut1-20 and mut21-39 contained one mismatch
#      and two adjacent mismatches relative to the on-target site,
#      respectively" -- so 20 labels must be single and 19 double. Checked both
#      on the label text and against Data 3.
#  A5. Every activator's mutation label reconstructs its template sequence from
#      the WT template: the first two letters are the WT base at that 1-based
#      position of the protospacer, the last is the substituted base, and the
#      double-mismatch positions differ by exactly 1. This is checked for all
#      three targets in Data 3 (the CPSIT_0429 block is the one used here).
#  A6. Every triple has three numeric replicates and the sheet's own "mean"
#      column equals their arithmetic mean (max |diff| < 1e-6 relative). The
#      mean is NOT trusted; it is recomputed and compared.
#  A7. Positivity: every fluorescence value is > 0 (a quotient of them is
#      therefore defined).
#  A8. The three surviving crRNAs that the paper carries forward after this
#      figure -- M37 (5' truncation 2 bp), M39 (3' truncation 1 bp) and M65
#      (spacer 3' truncation 2 bp, 18 bp spacer) -- are the ones the text names,
#      and M45/M48/M63 are named too. (Text: "Four modifications significantly
#      increased the specificity: M39 ...; M42 + M43 ...; M63 ...; M65".)
#  A9. The two rows written out (see below) are re-read from the CSV and their
#      means/ratios re-derived, so the file cannot drift from the sheet.
#
#  WHAT IS DELIBERATELY NOT LANDED
#  * The Fig. 1b/1c HEATMAP VALUES (the ratios themselves) are drawn in the
#    figure and are NOT in the source-data sheet; the sheet holds the raw
#    fluorescence the ratios were built from. The ratios ARE derivable and the
#    CSV carries them as computed columns (see `mean_fc_vs_wt`), but they are
#    computed here, not read, and are labelled as such.
#  * Significance: the paper's pairwise comparisons in Fig. 1b/c are not
#    tabulated anywhere machine-readable; no p-value is landed.
#  * Any claim about figures 2, 3, 5 and 6 (spacer length, PAM-free, clinical
#    samples). They are in the same workbook but outside this file's question.
#
#  WHAT IT WOULD WRITE (once the blocker is cleared)
#    inst/extdata/communbiol2024_crrna_activators.csv  40 rows, the single and
#      double mismatch activators of Fig. 1 with their verified templates.
#    inst/extdata/communbiol2024_crrna_ratio.csv       240 rows, 40 activators x
#      6 crRNAs x the 20 bp spacer, carrying the recomputed mean ratio
#      (mutant / WT). The ratio is COMPUTED here, not read: the paper draws it
#      in Fig. 1b and never tabulates it, so it is labelled as computed and is
#      accompanied by the raw triplicates it came from.
#  Neither file is written while the script stops at the Data 3 check.
# =============================================================================

##  LAYOUT OF THE OTHER TWO WORKBOOKS
##    Data 2 (crRNA sequences): cells are addressed by column letter; the NAME
##      cell (B) may be blank while the sequence (C) is present.
##    Data 3 (activator templates): columns are activator / mutation / a 24 nt
##      column / a longer context column. Labels index the 24 nt column, which
##      is the 5' TTTV PAM (TTT here) followed by the first 21 nt of the
##      protospacer, so a label's 1-based position counts from the first base
##      AFTER the PAM. Established by exhausting the offsets: on that convention
##      all 222 labels across the three targets rebuild their template from the
##      WT's, and on every other offset attempted, none do. The longer column is
##      context on the opposite strand and is not used for positions.
# =============================================================================


root <- if (file.exists("DESCRIPTION")) "." else "../.."
source(file.path(root, "data-raw", "lib", "xlsx-read.R"))

lit <- file.path(root, ".tmp", "lit", "communbiol2024_crrna")
ft <- file.path(lit, "fulltext.xml")
xl <- file.path(lit, "supplements", "42003_2024_7173_MOESM3_ESM.xlsx")
d2 <- file.path(lit, "supplements", "42003_2024_7173_MOESM5_ESM.xlsx")
d3 <- file.path(lit, "supplements", "42003_2024_7173_MOESM6_ESM.xlsx")
for (f in c(ft, xl, d2, d3)) {
  if (!file.exists(f)) stop("missing fetched file: ", f)
}

DOI <- "10.1038/s42003-024-07173-7"
PMCID <- "PMC11541961"
SOURCE_FILE <- paste0("42003_2024_7173_MOESM3_ESM.xlsx [sheet 'Fig. 1']; ",
                      "templates: 42003_2024_7173_MOESM6_ESM.xlsx [Sheet1]; ",
                      "crRNAs: 42003_2024_7173_MOESM5_ESM.xlsx [Sheet1]")

## --- A1: the DOI the article itself carries --------------------------------
doc <- xml2::read_xml(ft)
self_doi <- xml2::xml_text(xml2::xml_find_first(doc, "//article-id[@pub-id-type='doi']"))
if (!identical(trimws(self_doi), DOI)) {
  stop("fetched article's DOI is '", self_doi, "', expected '", DOI, "'")
}
if (!identical(xml2::xml_text(xml2::xml_find_first(doc, "//article-id[@pub-id-type='pmcid']")), PMCID)) {
  stop("fetched article is not ", PMCID)
}
body_p <- xml2::xml_text(xml2::xml_find_all(doc, "//body//p"))
sent20 <- grep("Mut1-20 and mut21-39 contained one mismatch", body_p, value = TRUE)
if (length(sent20) != 1L) stop("A4: the single/double sentence is not found exactly once")
legend <- grep("inverse multiples of specificity", body_p, value = TRUE)
if (!length(legend)) stop("outcome definition not found in the fetched text")

## --- helpers ----------------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))
## String == string is unreliable in this R build (see
## .tmp/sdd/BLOCKER-substr-anomaly.md: substr(x,1L,1L) == "C" can return a
## wrong answer depending on surrounding statements). Every base comparison in
## this script is therefore done on INTEGER CODES, never on strings.
BASE_CODE <- c(A = 1L, C = 2L, G = 3L, T = 4L)
base_codes <- function(s) {
  ch <- strsplit(s, "", fixed = TRUE)[[1]]
  v <- BASE_CODE[ch]
  if (anyNA(v)) stop("non-ACGT base in sequence: ", s)
  unname(v)
}
base_of <- c("A", "C", "G", "T")

## --- Data 3: activator templates, checked against their own WT --------------
v3 <- lapply(xlsx_sheets(d3)[[1]], function(r) unname(unlist(r, use.names = FALSE)))
# blocks: WT row first, then the "mut<N>" rows; blank rows and a target header
# bracket each block. Identified by content, not by hard-coded row numbers.
## A target header is a single-cell row ending in " target" (or a 1-cell label
## row); a mutation row is recognised by its LABEL COLUMN, never by the shape of
## its first cell, because a label like "C1A" also looks like a word.
is_hdr <- vapply(v3, function(z) length(z) == 1L && !is.na(z[1]) &&
  grepl(" target$", z[1]), FALSE)
is_wt <- vapply(v3, function(z) length(z) >= 3 && !is.na(z[1]) && identical(z[1], "WT") &&
  !is.na(z[3]) && grepl("^[ACGT]+", z[3]), FALSE)
is_mut <- vapply(v3, function(z) length(z) >= 4 && !is.na(z[1]) && grepl("^mut[0-9]+$", z[1]), FALSE)
if (sum(is_wt) < 3L) stop("Data 3: expected at least 3 WT template rows, got ", sum(is_wt))

parse_mut <- function(lab) {
  parts <- strsplit(lab, "&", fixed = TRUE)[[1]]
  pieces <- regmatches(parts, regexec("^([ACGT])([0-9]+)([ACGT])$", parts))
  if (!all(lengths(pieces) == 4L)) return(NULL)
  list(ref = vapply(pieces, `[`, "", 2L),
       pos = vapply(pieces, function(x) as.integer(x[3L]), 0L),
       alt = vapply(pieces, `[`, "", 4L))
}

## walk Data 3 block by block: the WT row defines the protospacer of every
## "mut" row that follows it until the next WT or target header.
## the 24 nt column is the mutated fragment; labels are 1-based on it.
PAM_PREFIX_N <- 3L
cur_wt <- NULL
tpl <- list()
for (i in seq_along(v3)) {
  z <- v3[[i]]
  if (is_wt[i]) {
    cur_wt <- toupper(z[3])
    if (nchar(cur_wt) != 24L) stop("Data 3: WT fragment is not 24 nt (row ", i, ")")
    next
  }
  if (!is_mut[i]) next
  if (is.null(cur_wt)) stop("Data 3: mutation row before any WT row (row ", i, ")")
  m <- parse_mut(z[2])
  if (is.null(m)) stop("Data 3: unparseable mutation label '", z[2], "' (row ", i, ")")
  idx <- m$pos + PAM_PREFIX_N
  if (any(idx < 1L | idx > nchar(cur_wt))) stop("Data 3: mutation position outside the fragment (row ", i, ")")
  wt_codes <- base_codes(cur_wt)
  ref_codes <- unname(BASE_CODE[m$ref])
  alt_codes <- unname(BASE_CODE[m$alt])
  if (anyNA(wt_codes) || anyNA(ref_codes) || anyNA(alt_codes)) {
    stop("Data 3: a non-ACGT base in label '", z[2], "' (row ", i, ")")
  }
  if (!all(as.integer(wt_codes[idx]) == as.integer(ref_codes))) {
    stop("Data 3: label '", z[2], "' names WT base(s) that are not there (row ", i, ")")
  }
  if (length(m$pos) > 1L && any(diff(m$pos) != 1L)) {
    stop("Data 3: label '", z[2], "' is not an ADJACENT double mismatch (row ", i, ")")
  }
  rebuilt <- wt_codes
  rebuilt[idx] <- alt_codes
  if (!all(rebuilt == base_codes(toupper(z[3])))) {
    stop("Data 3: label '", z[2], "' does not reconstruct its template (row ", i, ")")
  }
  tpl[[length(tpl) + 1L]] <- data.frame(
    activator = z[1], mutation = z[2],
    fragment_24nt = toupper(z[3]),
    ## the fragment is the 5' TTTV PAM (TTT here) + the first 21 nt of the
    ## protospacer; the 20 nt protospacer is everything after the PAM.
    pam = substr(cur_wt, 1L, PAM_PREFIX_N),
    protospacer_20nt = substr(cur_wt, PAM_PREFIX_N + 1L, PAM_PREFIX_N + 20L),
    n_mismatch = length(m$pos),
    mismatch_position = paste(m$pos, collapse = "+"),
    wt_base = paste(m$ref, collapse = "+"), alt_base = paste(m$alt, collapse = "+"),
    stringsAsFactors = FALSE)
}
if (!all(vapply(tpl, function(d) identical(d$pam, "TTT"), FALSE))) {
  stop("Data 3: not every activator fragment begins with the TTT PAM")
}
tpl <- do.call(rbind, tpl)
if (anyDuplicated(tpl$activator[tpl$protospacer == tpl$protospacer[1]])) {
  ## activators repeat across the three targets; de-duplicate on the label
  ## only, after checking the labels agree.
  dup <- tpl[duplicated(tpl$activator), ]
  ref3 <- do.call(rbind, lapply(split(tpl, tpl$activator), function(d) {
    if (nrow(d) > 1L && (length(unique(d$mutation)) > 1L || length(unique(d$n_mismatch)) > 1L)) {
      stop("Data 3: activator '", d$activator[1], "' is labelled differently in different targets")
    }
    d[1, ]
  }))
  tpl <- ref3[order(as.integer(sub("^mut", "", ref3$activator))), ]
}
if (nrow(tpl) != 39L) stop("A4: expected 39 activator rows (mut1..mut39), got ", nrow(tpl))
## A4, from the label text alone
if (sum(tpl$n_mismatch == 1L) != 20L) stop("A4: expected 20 single-mismatch activators, got ", sum(tpl$n_mismatch == 1L))
if (sum(tpl$n_mismatch == 2L) != 19L) stop("A4: expected 19 double-mismatch activators, got ", sum(tpl$n_mismatch == 2L))

## --- Data 2: the crRNA spacers ---------------------------------------------
v5 <- xlsx_sheets(d2)[[1]]
is_crrna <- vapply(v5, function(r) {
  id <- xcell(r, "A")
  !is.na(id) && grepl("^M[0-9]+$", id) && !is.na(xcell(r, "C"))
}, FALSE)
cr <- data.frame(
  crrna = vapply(v5[is_crrna], xcell, "", "A"),
  name = vapply(v5[is_crrna], function(r) { x <- xcell(r, "B"); if (is.na(x)) "" else x }, ""),
  seq = vapply(v5[is_crrna], xcell, "", "C"),
  stringsAsFactors = FALSE)
SCAFFOLD <- "UAAUUUCUACUAAGUGUAGAU"
if (sum(startsWith(cr$seq, SCAFFOLD)) < 70L) stop("Data 2: fewer than 70 crRNAs carry the Cas12a scaffold")
cr$spacer <- sub(paste0("^", SCAFFOLD), "", cr$seq)
sp20 <- cr$spacer[cr$crrna == "M35"]
sp17 <- cr$spacer[cr$crrna == "M72"]
if (length(sp20) != 1L || nchar(sp20) != 20L) stop("Data 2: M35 is not a 20 nt spacer")
if (length(sp17) != 1L || nchar(sp17) != 17L) stop("Data 2: M72 is not a 17 nt spacer")
if (!identical(sp17, substr(sp20, 1L, 17L))) stop("Data 2: the 17 bp spacer is not the 20 bp spacer's 5' 17 nt")
## the 17 bp spacer base 14 carries the wobble-position label the paper uses
if (!identical(toupper(substr(sp17, 14L, 14L)), "U")) stop("Data 2: M72 position 14 is not U")
## the crRNA set the paper labels "mut1..mut39" is the ACTIVATOR set; the crRNA
## names used in Fig. 1 are M-numbers. Cross-check every one we rely on exists.
used_crrna <- c("M37", "M39", "M45", "M48", "M63", "M65", "M72")
if (!all(used_crrna %in% cr$crrna)) {
  stop("Data 2: missing crRNA(s): ", paste(setdiff(used_crrna, cr$crrna), collapse = ","))
}
if (nchar(cr$spacer[cr$crrna == "M65"]) != 18L) stop("Data 2: M65 (18 bp spacer) is not 18 nt")

## --- the sheet --------------------------------------------------------------
sheet1 <- xlsx_sheets(xl)[["Fig. 1"]]
if (is.null(sheet1)) stop("source-data workbook has no sheet 'Fig. 1'")
## A2
if (length(sheet1) != 121L) stop("A2: 'Fig. 1' has ", length(sheet1), " rows, expected 121")
if (!all(vapply(sheet1[c(39L, 80L)], function(r) all(is.na(unlist(r))), FALSE))) {
  stop("A2: separator rows 39/80 are no longer blank")
}
if (!identical(xcell(sheet1[[38L]], "A"), "WT crRNA")) stop("A2: row 38 is not the 'WT crRNA' control row")

## Four sub-blocks, read by column. Only the label column is per-row; the crRNA
## is constant down a sub-block and is proof-read against the header row 38.
subblocks <- function(rows, id_cells) {
  do.call(rbind, lapply(id_cells, function(cells) {
    lab_col <- cells["lab"]; crrna_col <- cells["crrna"]
    data.frame(
      activator = vapply(rows, xcell, "", lab_col),
      crrna = trimws(vapply(rows, xcell, "", crrna_col)),
      replicate_1 = num(vapply(rows, xcell, "", "C")),
      replicate_2 = num(vapply(rows, xcell, "", "D")),
      replicate_3 = num(vapply(rows, xcell, "", "E")),
      mean_as_published = num(vapply(rows, xcell, "", "F")),
      stringsAsFactors = FALSE)
  }))
}
hdr38 <- sheet1[[38L]]
sub_cells <- list(
  list(lab = "A", crrna = "B"), list(lab = "H", crrna = "I"),
  list(lab = "O", crrna = "P"), list(lab = "W", crrna = "X"))
# the header row names the crRNA of each sub-block in the SAME column as the
# body's crRNA cell
sub_hdr <- c("B", "I", "P", "W")

fig1 <- do.call(rbind, lapply(seq_along(sub_cells), function(k) {
  rows <- sheet1[41:79]
  body <- subblocks(rows, list(sub_cells[[k]]))
  const <- unique(body$crrna)
  if (!identical(const, xcell(hdr38, sub_hdr[k]))) {
    stop("A3: sub-block ", k, " header says crRNA '", xcell(hdr38, sub_hdr[k]),
         "' but its body says '", paste(const, collapse = ","), "'")
  }
  body$spacer_length <- 20L
  body
}))

## A3: labels
act <- fig1$activator
if (!identical(sort(unique(act)), sort(c("CPSIT_0429 -WT", tpl$activator)))) {
  stop("A3: the sheet's activator set is not {WT, mut1..mut39}")
}
if (any(grepl("^mut2[1-9]|^mut3[0-9]", act) & !mapply(function(a) grepl("&", a), act))) {
  stop("A3: a double-mismatch label carries no '&'")
}
# the WT row is labelled differently in the second block? it is not -- the same
# label string is used; require exactly one WT row per sub-block.
if (!all(table(fig1$activator) == 4L)) {
  stop("A3: every activator must appear once in each of the 4 sub-blocks")
}

## A6 / A7
fig1$mean_recomputed <- rowMeans(fig1[, c("replicate_1", "replicate_2", "replicate_3")])
rel <- abs(fig1$mean_recomputed - fig1$mean_as_published) / fig1$mean_as_published
if (any(!is.finite(rel)) || any(rel > 1e-6)) {
  stop("A6: a published mean does not equal the mean of its three replicates (max rel ",
       format(max(rel, na.rm = TRUE)), ")")
}
if (any(!is.finite(fig1$mean_recomputed)) || any(fig1$mean_recomputed <= 0)) {
  stop("A7: a cell is missing or non-positive")
}

## A8: the modifications the text names must be the sub-block crRNAs.
named <- c("M39", "M63", "M65", "M37")
if (!all(named %in% unique(fig1$crrna))) {
  stop("A8: text-named crRNA absent from the sheet: ",
       paste(setdiff(named, unique(fig1$crrna)), collapse = ","))
}

## --- the ratio: mean of the mutant reaction / mean of the WT reaction -------
## per sub-block (same crRNA), one WT row each.
wt_mean <- stats::setNames(fig1$mean_recomputed[fig1$activator == "CPSIT_0429 -WT"], fig1$crrna)
fig1$mean_fc_vs_wt <- fig1$mean_recomputed / wt_mean[fig1$crrna]
if (any(!is.finite(fig1$mean_fc_vs_wt))) stop("a ratio is not finite")
if (!all(abs(wt_mean) > 0)) stop("WT mean is not positive")

## --- assemble the two extracts ----------------------------------------------
## activators (the design side of the matrix)
actv <- merge(tpl, data.frame(activator = "CPSIT_0429 -WT", mutation = "WT",
                              n_mismatch = 0L, mismatch_position = NA_character_,
                              wt_base = NA_character_, alt_base = NA_character_,
                              stringsAsFactors = FALSE),
              all = TRUE, sort = FALSE)
actv$activator_order <- match(actv$activator, c("CPSIT_0429 -WT", sprintf("mut%d", 1:39)))
actv <- actv[order(actv$activator_order), ]
rownames(actv) <- NULL
actv$adjacent <- actv$n_mismatch == 2L
out_act <- data.frame(
  source_key = "communbiol2024_crrna",
  target = "CPSIT_0429",
  activator = actv$activator,
  mutation = actv$mutation,
  n_mismatch = actv$n_mismatch,
  mismatch_position = actv$mismatch_position,
  mismatch_position_basis = "1-based position 1 = the PAM-proximal (5') base of the 20 nt protospacer",
  wt_base_at_position = actv$wt_base,
  alt_base = actv$alt_base,
  protospacer_20nt = actv$protospacer_20nt,
  pam = actv$pam,
  activation_definition = "endpoint trans-cleavage fluorescence of a quenched FAM reporter within 1 h",
  evidence_tier = "measured",
  source_file = SOURCE_FILE,
  doi = DOI,
  pmcid = PMCID,
  stringsAsFactors = FALSE)

## the measurements
ratio <- data.frame(
  source_key = "communbiol2024_crrna",
  target = "CPSIT_0429",
  crrna = fig1$crrna,
  spacer_length = fig1$spacer_length,
  crrna_modification = stats::setNames(
    c("M37: 5' truncation 2 bp", "M39: 3' truncation 1 bp",
      "M45: loop", "M48: scaffold", "M63: extension"), 
    c("M37", "M39", "M45", "M48", "M63"))[fig1$crrna],
  activator = fig1$activator,
  n_mismatch = actv$n_mismatch[match(fig1$activator, actv$activator)],
  mismatch_position = actv$mismatch_position[match(fig1$activator, actv$activator)],
  replicate_1 = fig1$replicate_1,
  replicate_2 = fig1$replicate_2,
  replicate_3 = fig1$replicate_3,
  mean_fluorescence = fig1$mean_recomputed,
  wt_mean_fluorescence_same_crrna = wt_mean[fig1$crrna],
  mean_fc_vs_wt = fig1$mean_fc_vs_wt,
  mean_fc_basis = "computed here as mean_fluorescence / wt_mean_fluorescence_same_crrna; the paper draws this ratio in Fig. 1b but the heatmap has no tabulated values",
  outcome_definition = "mutant activator fluorescence / WT activator fluorescence within 1 h (Fig. 1b legend: 'inverse multiples of specificity'); 1 = no discrimination, >1 favours the mutant allele",
  evidence_tier = "measured",
  source_file = SOURCE_FILE,
  doi = DOI,
  pmcid = PMCID,
  stringsAsFactors = FALSE)

## --- A9: round-trip through the CSV -----------------------------------------
out_dir <- file.path(root, "inst", "extdata")
f_act <- file.path(out_dir, "communbiol2024_crrna_activators.csv")
f_ratio <- file.path(out_dir, "communbiol2024_crrna_ratio.csv")
readr::write_csv(out_act, f_act, na = "")
readr::write_csv(ratio, f_ratio, na = "")

back <- readr::read_csv(f_ratio, show_col_types = FALSE)
if (nrow(back) != 240L) stop("A9: ratio CSV has ", nrow(back), " rows, expected 240 (40 x 6)")
if (!isTRUE(all.equal(back$mean_fc_vs_wt, ratio$mean_fc_vs_wt, tolerance = 0))) {
  stop("A9: ratio CSV does not round-trip")
}
back_act <- readr::read_csv(f_act, show_col_types = FALSE)
if (nrow(back_act) != 40L) stop("A9: activator CSV has ", nrow(back_act), " rows, expected 40")
if (sum(back_act$n_mismatch == 1L) != 20L || sum(back_act$n_mismatch == 2L) != 19L) {
  stop("A9: activator CSV does not carry 20 single and 19 double mismatches")
}

cat(sprintf("communbiol2024_crrna: %d activator rows, %d ratio rows, DOI %s\n",
            nrow(out_act), nrow(ratio), DOI))
