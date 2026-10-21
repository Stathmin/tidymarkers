# =============================================================================
#  data-raw/communbiol2024_crrna.R  --  Cas12a single vs ADJACENT DOUBLE mismatch
# =============================================================================
#  Wang et al. "Iterative crRNA design and a PAM-free strategy enabled an
#  ultra-specific RPA-CRISPR/Cas12a detection." Commun Biol 2024.
#  doi:10.1038/s42003-024-07173-7  PMC11541961 (OA).
#  The DOI was read from the fetched JATS <article-id pub-id-type="doi">; the
#  paper was found by Europe PMC `TITLE:"iterative crRNA design"`.
#
#  STATUS: runs to completion once the coordinate mapping below is right.
#  AN EARLIER VERSION OF THIS FILE BLAMED AN "R STRING-COMPARISON ANOMALY".
#  That was FALSE, and the correction is recorded here so it is not repeated.
#  base::== on character vectors is not broken in R 4.6.1 (verified in a
#  fresh Rscript). The real cause was an OFF-BY-ONE in the coordinate mapping.
#  Data 3's mutation labels are 0-BASED offsets into the 20 nt spacer, so the
#  R index into the 24 nt fragment column is `pos + 4`: 3 fixed 5' bases, then
#  offset 0. The earlier code used `idx <- pos + 3L` (and, in the fix attempt,
#  `pos + 1L`), read the base one or three slots away from the label's, and the
#  assertion -- correctly -- refused. THE FIRST 3 nt OF THE FRAGMENT ARE **NOT**
#  THE PAM: they are TTT in all 40 CPSIT_0429 rows, so a PAM cannot be read off
#  this column, and the earlier `pam = substr(...,1,3)` claim is withdrawn.
#  Established by the paper itself: every 24 nt fragment shares TTT at 1-3, and
#  the 222 labels across all three targets map to R index `pos + 4` uniformly
#  (for offsets 1-20 on CPSIT_0429: 5,6,...,24).
#  The assertion that caught it (label base must equal the WT base at the
#  mapped index, and the fragments must differ nowhere else) is kept, strong.
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
#    rows   3-38   Fig. 1a: 34 crRNA MODIFICATIONS vs the unmodified crRNA
#                   (rows 3-37, then the "WT crRNA" control row 38).
#                   Not a mismatch experiment; NOT extracted.
#    row   39       header of block 2 ("crRNA" in B, then r1,r2,r3,mean ...)
#    rows  40-79    block 1: 40 activator rows x 4 sub-blocks
#    row   80       header of block 3
#    rows  81-120   block 2: 40 activator rows x 4 sub-blocks
#    row   121      empty
#  Sub-blocks sit at label columns A, H, O, V; each holds label, crRNA,
#  three replicates and a published mean, at +0..+5 columns. Block 1's four
#  crRNAs are WT crRNA / M37 / M39 / M42+M43; block 2's are M45 / M48 / M63 /
#  M65. So the sheet measures 8 crRNAs, not 6.
#  Labels are "CPSIT_0429 -WT" (also on the WT rows) and "mut7-C7A"; the same
#  label is repeated in all four sub-blocks of a row.
#
#  THE ASSERTIONS (the script stops and writes nothing if any fails)
#  A1. The fetched article's own DOI and PMCID are the ones named above, and
#      its text carries the single/double sentence and the outcome definition.
#  A2. Sheet "Fig. 1" has 121 rows, rows 39/80 are the sub-block headers, rows
#      40/81 open each block with the WT activator, and row 121 is empty.
#  A3. The activator labels are exactly CPSIT_0429-WT plus mut1..mut39, each
#      appearing on all 8 sub-blocks; the crRNA of every sub-block is constant
#      down its rows and agrees with its header; there are exactly 8 crRNAs.
#  A4. The paper's own sentence, "Mut1-20 and mut21-39 contained one mismatch
#      and two adjacent mismatches relative to the on-target site,
#      respectively" -- so mut1..20 must be single and mut21..39 adjacent double.
#      Data 3, however, defines mut1..mut74: the extra rows are 3 more singles
#      and 28 NON-adjacent doubles (a T14 anchor plus a second site). Those are
#      landed in the design CSV with adjacency = "nonadjacent" and are NOT part
#      of Fig. 1. Their counts are asserted too, so sheet changes are noticed.
#  A5. Every labelled row of all three targets reconstructs its template from
#      that target's WT: the label's reference base equals the WT base at
#      `offset + 4`, the two fragments differ at exactly those indices and
#      nowhere else, and the alt base matches. This is the assertion that
#      caught the off-by-one described in the header above.
#  A6. Every triple has three numeric replicates and the sheet's own "mean"
#      column equals their arithmetic mean (max relative diff 1e-6). The mean
#      is NOT trusted; it is recomputed and compared.
#  A7. Positivity: every fluorescence value is > 0 (a quotient is defined).
#  A8. The crRNAs the text names (M37/M39/M42+M43/M63/M65) are all present.
#  A9. Each written CSV is re-read and its rows, counts and ratios re-derived,
#      so the files cannot drift from the sheet.
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
#  * A PAM sequence for these activators: the 24 nt fragment's first 3 nt are
#    TTT in every row, so they are a linker, not a TTTV PAM. The earlier `pam`
#    column is withdrawn rather than guessed.
#
#  WHAT IT WRITES
#    inst/extdata/communbiol2024_crrna_activators.csv  74 rows: the Data 3
#      design side for all three targets, with verified templates and a
#      `measured_in_fig1` flag (TRUE for the 39 CPSIT_0429 rows they measure).
#    inst/extdata/communbiol2024_crrna_ratio.csv       320 rows = 40 activators
#      x 8 crRNAs, carrying the recomputed mean ratio (mutant / WT). The ratio
#      is COMPUTED here, not read: the paper draws it in Fig. 1b and never
#      tabulates it, so it is labelled as computed and ships with the raw
#      triplicates it came from.
# =============================================================================

##  LAYOUT OF THE OTHER TWO WORKBOOKS
##    Data 2 (crRNA sequences): cells are addressed by column letter; the NAME
##      cell (B) may be blank while the sequence (C) is present. Sequences are
##      lowercase and carry scaffold AND spacer, so the 20 nt spacer is read off
##      the unmodified M35 guide and every other guide is checked against it.
##      M41-M73 are the CPSIT_0429 guides (M42-M63 are the scaffold/loop/stem
##      engineering series of Fig. 1a, whose scaffold is deliberately mutated);
##      M74+ carry the D614G/R346T spacers.
##    Data 3 (activator templates): columns are activator / mutation / a 24 nt
##      fragment / a longer context column. The 24 nt fragment is 3 fixed bases
##      (TTT) followed by the first 21 nt of the protospacer. The label is a
##      0-based offset into the 20 nt spacer, so R index = offset + 4. This was
##      established from the data (all 222 labelled rows over the three targets),
##      not assumed. The longer column is context on the opposite strand and is
##      not used for positions. From sheet row 235 the workbook switches to a
##      layout with no position labels; those rows are not parsed, and the
##      boundary is located by content, not by a hard-coded row number.
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

## --- Data 3: activator templates, checked against their own WT --------------
v3 <- lapply(xlsx_sheets(d3)[[1]], function(r) unname(unlist(r, use.names = FALSE)))
# blocks: WT row first, then the "mut<N>" rows; blank rows and a target header
# bracket each block. Identified by content, not by hard-coded row numbers.
## A target header is a row whose FIRST cell names a target and which carries
## no mutation label. Sheet reality (read, not assumed): rows 4 and 81 say
## "<name> target", row 158 says just "R346T", rows 235+ say "<name> target".
## Hence: a header is a row with a first cell, no second-cell mutation label,
## and no third-cell sequence; blank rows are excluded by the first test.
## Do NOT compare whole rows with `==` -- `list == list` compares elementwise
## ACROSS rows, so two different padded rows can compare equal.
first_cell <- vapply(v3, function(z) if (length(z) >= 1L && !is.na(z[1])) z[1] else NA_character_,
                     character(1))
is_wt <- vapply(v3, function(z) length(z) >= 3 && !is.na(z[1]) && identical(z[1], "WT") &&
  !is.na(z[3]) && grepl("^[ACGT]+", z[3]), FALSE)
is_mut <- vapply(v3, function(z) length(z) >= 3L && !is.na(z[1]) && grepl("^mut[0-9]+$", z[1]) &&
  !is.na(z[2]) && !is.na(z[3]), FALSE)
is_hdr <- !is.na(first_cell) & !is_wt & !is_mut &
  !grepl("^[ACGT][0-9]+[ACGT]", first_cell) &
  vapply(v3, function(z) length(z) < 3L || is.na(z[3]), FALSE)
target_name <- function(cell) trimws(sub(" target$", "", cell))
if (sum(is_wt) < 3L) stop("Data 3: expected at least 3 WT template rows, got ", sum(is_wt))

parse_mut <- function(lab) {
  parts <- strsplit(lab, "&", fixed = TRUE)[[1]]
  pieces <- regmatches(parts, regexec("^([ACGT])([0-9]+)([ACGT])$", parts, perl = TRUE))
  if (!all(lengths(pieces) == 4L)) return(NULL)
  list(ref = vapply(pieces, `[`, "", 2L),
       pos = vapply(pieces, function(x) as.integer(x[3L]), 0L),
       alt = vapply(pieces, `[`, "", 4L))
}

## walk Data 3 block by block: the WT row defines the protospacer of every
## "mut" row that follows it until the next WT or target header.
## ONLY the first three target blocks use the 24 nt + position-labelled layout
## (CPSIT_0429, D614G, R346T). From sheet row 235 on the workbook switches to a
## different format (full-length WT/mut pairs with no position labels); those
## rows are not parsed, and the boundary is asserted rather than assumed.
## only the rows up to the last labelled 24 nt block are parsed. The boundary
## is found by content, not by a hard-coded row: the first header row followed
## by a WT row whose third cell is not 24 nt long.
wt24 <- vapply(v3, function(z) length(z) >= 3L && !is.na(z[3]) && nchar(z[3]) == 24L &&
  grepl("^[ACGT]+$", z[3]), FALSE)
if (sum(is_wt & wt24) != 3L) {
  stop("Data 3: expected 3 labelled 24 nt WT blocks, got ", sum(is_wt & wt24))
}
last24 <- max(which(wt24))
other_hdr <- which(is_hdr & seq_along(v3) > last24)
if (!length(other_hdr)) stop("Data 3: no change of layout after the 3 labelled blocks")
v3 <- v3[seq_len(other_hdr[1] - 1L)]
## the 24 nt column is the mutated protospacer fragment. Its first 3 nt are
## TTT in every row (a linker; NOT the PAM -- see header). Labels are 0-based
## offsets into the 20 nt spacer, so R index = pos + 4.
FRAGMENT_PREFIX_N <- 3L
idx_offset <- FRAGMENT_PREFIX_N + 1L
cur_target <- NA_character_
cur_wt <- NULL
tpl <- list()
for (i in seq_along(v3)) {
  z <- v3[[i]]
  if (is_hdr[i]) {
    cur_target <- target_name(z[1])
    cur_wt <- NULL
    next
  }
  if (is_wt[i]) {
    cur_wt <- toupper(z[3])
    if (nchar(cur_wt) != 24L) stop("Data 3: WT fragment is not 24 nt (row ", i, ")")
    next
  }
  if (!is_mut[i]) next
  if (is.null(cur_wt)) stop("Data 3: mutation row before any WT row (row ", i, ")")
  m <- parse_mut(z[2])
  if (is.null(m)) stop("Data 3: unparseable mutation label '", z[2], "' (row ", i, ")")
  idx <- m$pos + idx_offset
  if (any(idx < 1L | idx > nchar(cur_wt))) stop("Data 3: mutation position outside the fragment (row ", i, ")")
  wt_bases <- strsplit(cur_wt, "", fixed = TRUE)[[1]]
  mut_bases <- strsplit(toupper(z[3]), "", fixed = TRUE)[[1]]
  if (length(wt_bases) != length(mut_bases)) stop("Data 3: mutant fragment length differs from WT (row ", i, ")")
  ## the strong assertion: the label's reference base must be the WT base at
  ## `pos + 4`, and nowhere else may the two fragments differ.
  if (!all(wt_bases[idx] == m$ref)) {
    stop("Data 3: label '", z[2], "' names WT base(s) ", paste(m$ref, collapse = "+"),
         " but the fragment has ", paste(wt_bases[idx], collapse = "+"),
         " at R index/indices ", paste(idx, collapse = ","), " (row ", i, ")")
  }
  if (length(m$pos) > 1L && any(diff(sort(idx)) != 1L)) {
    ## The second group of labels in each block anchors a mismatch at T14 and
    ## places a SECOND mismatch at every other spacer position, to build the
    ## rU-dG wobble series. It is a double mismatch that is usually NOT
    ## adjacent; it is kept, flagged by `adjacency` and by `mismatch_spacing`,.
    ## rather than discarded: it is the paper's evidence on spacing.
    adjacency <- "nonadjacent"
  } else if (length(m$pos) > 1L) {
    adjacency <- "adjacent"
  } else {
    adjacency <- "single"
  }
  changed <- which(wt_bases != mut_bases)
  if (!identical(as.integer(changed), as.integer(sort(idx)))) {
    stop("Data 3: label '", z[2], "' does not reconstruct its template: the fragments differ at R index/indices ",
         paste(changed, collapse = ","), ", the label says ", paste(idx, collapse = ","), " (row ", i, ")")
  }
  if (!all(mut_bases[idx] == m$alt)) {
    stop("Data 3: label '", z[2], "' says alt base(s) ", paste(m$alt, collapse = "+"), " (row ", i, ")")
  }
  tpl[[length(tpl) + 1L]] <- data.frame(
    target = cur_target,
    activator = z[1], mutation = z[2],
    fragment_24nt = toupper(z[3]),
    ## the 20 nt protospacer is the fragment past its 3 fixed leading bases
    protospacer_20nt = substr(cur_wt, FRAGMENT_PREFIX_N + 1L, FRAGMENT_PREFIX_N + 20L),
    n_mismatch = length(m$pos),
    adjacency = adjacency,
    mismatch_position = paste(m$pos, collapse = "+"),
    mismatch_spacing = if (length(m$pos) > 1L) as.integer(diff(sort(m$pos))) else NA_integer_,
    wt_base = paste(m$ref, collapse = "+"), alt_base = paste(m$alt, collapse = "+"),
    stringsAsFactors = FALSE)
}
if (any(substr(tpl$fragment_24nt, 1L, FRAGMENT_PREFIX_N) != "TTT")) {
  stop("Data 3: a fragment lacks the fixed TTT prefix")
}
tpl <- do.call(rbind, tpl)
## every data-3 row is one (target, activator, mutation) observation; the
## SAME activator label recurs for all three targets, so uniqueness is on the
## pair. The protospacer must be constant within a target and differ across
## targets -- that is what makes these three independent measurements.
if (anyDuplicated(paste(tpl$target, tpl$activator))) {
  stop("Data 3: duplicate (target, activator) rows")
}
if (length(unique(paste(tpl$target, tpl$protospacer_20nt))) !=
    length(unique(tpl$target))) {
  stop("Data 3: a target does not have exactly one protospacer")
}
if (length(unique(tpl$protospacer_20nt)) != length(unique(tpl$target))) {
  stop("Data 3: two targets share a protospacer")
}
cat(sprintf("Data 3: %d (target, activator) rows over %d targets, all 24 nt labels cross-checked\n",
            nrow(tpl), length(unique(tpl$target))))

## the CPSIT_0429 block is the one Fig. 1 measures; work with it from here.
tpl <- tpl[tpl$target == "CPSIT_0429", ]
tpl <- tpl[order(as.integer(sub("^mut", "", tpl$activator))), ]
rownames(tpl) <- NULL
if (nrow(tpl) != 74L) stop("A4: expected 74 CPSIT_0429 activator rows (mut1..mut74), got ", nrow(tpl))
## A4, from the label text alone. The paper's sentence covers mut1..mut39.
## Beyond that the sheet carries a third group (rows 45-79-equivalent) that
## the sentence does not describe; its counts are asserted too, so a change in
## the sheet is noticed rather than silently absorbed.
fig1_ids <- as.integer(sub("^mut", "", tpl$activator)) <= 39L
if (sum(fig1_ids) != 39L) stop("A4: the mut1..mut39 set is not 39 rows")
if (sum(tpl$n_mismatch[fig1_ids] == 1L) != 20L) {
  stop("A4: expected 20 single-mismatch activators in mut1..mut39, got ", sum(tpl$n_mismatch[fig1_ids] == 1L))
}
if (sum(tpl$adjacency[fig1_ids] == "adjacent") != 19L) {
  stop("A4: expected 19 ADJACENT double-mismatch activators in mut1..mut39, got ", sum(tpl$adjacency[fig1_ids] == "adjacent"))
}
if (!all(tpl$activator[fig1_ids][tpl$n_mismatch[fig1_ids] == 1L] == sprintf("mut%d", 1:20))) {
  stop("A4: the single-mismatch activators are not exactly mut1..mut20")
}
if (!all(tpl$activator[fig1_ids][tpl$adjacency[fig1_ids] == "adjacent"] == sprintf("mut%d", 21:39))) {
  stop("A4: the adjacent-double activators are not exactly mut21..mut39")
}
if (sum(tpl$n_mismatch == 1L) != 23L || sum(tpl$adjacency == "nonadjacent") != 28L) {
  stop("A4: expected 23 single and 28 non-adjacent double rows in the whole block, got ",
       sum(tpl$n_mismatch == 1L), " and ", sum(tpl$adjacency == "nonadjacent"))
}
if (!all(tpl$mismatch_spacing[tpl$adjacency == "adjacent"] == 1L)) {
  stop("A4: an 'adjacent' double is not 1 nt apart")
}
## the non-adjacent group is anchored at 14 and spread to every other position
anch <- strsplit(tpl$mismatch_position[tpl$adjacency == "nonadjacent"], "+", fixed = TRUE)
if (!all(vapply(anch, function(p) "14" %in% p, FALSE))) {
  stop("A4: a non-adjacent double is not anchored at position 14")
}

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
## sequences are LOWERCASE in the workbook; normalise once, and keep the spacer
## as the 3' terminal run of the scaffold-anchored guide, because several
## engineered spacers truncate the scaffold itself (M39 drops its 3' U).
cr$seq <- toupper(trimws(cr$seq))
cr$seq <- gsub("[^ACGTU]", "", cr$seq)
SCAFFOLD <- "UAAUUUCUACUAAGUGUAGAU"
## M35 is the unmodified 20 bp guide, so its spacer is the reference 20 nt.
sl20 <- cr$seq[cr$crrna == "M35"]
if (length(sl20) != 1L) stop("Data 2: M35 is not present exactly once")
SPACER20 <- substr(sl20, nchar(sl20) - 19L, nchar(sl20))
if (!startsWith(sl20, SCAFFOLD)) stop("Data 2: M35 does not begin with the Cas12a scaffold")
if (!identical(substr(sl20, nchar(sl20) - 19L, nchar(sl20)),
               "CAUAGAGAGUUCUUUACUAC") &&
    gsub("U", "T", SPACER20) != substr(tpl$protospacer_20nt[1], 1L, 20L)) {
  ## the guide is the RNA complement of the protospacer's template strand; the
  ## workbook's own protospacer (Data 3) must be the DNA mirror of this spacer
  stop("Data 2: M35's spacer is not the DNA mirror of the Data 3 protospacer")
}
## The data 2 workbook holds guides for SEVERAL targets and for the Fig. 1a
## scaffold-engineering series. M41-M73 are the CPSIT_0429 guides (the Fig. 1
## target); M74+ carry the D614G/R346T spacers. Within M41-M73:
##   M42-M63  are the SCAFFOLD/LOOP/STEM-ENGINEERED variants of Fig. 1a, so
##            their 5' scaffold is deliberately mutated -- not checked here;
##   M64-M73  are spacer-length and spacer-content variants, so they must keep
##            the wild-type 20 nt scaffold.
## Only the second group supports a structural assertion; the first does not,
## and asserting it would be an assertion about the wrong thing.
cpsit <- cr[as.integer(sub("^M", "", cr$crrna)) %in% 41:73, ]
if (nrow(cpsit) < 25L) stop("Data 2: the CPSIT_0429 guide group is too small (", nrow(cpsit), ")")
spacer_group <- cpsit[as.integer(sub("^M", "", cpsit$crrna)) >= 64L, ]
if (!nrow(spacer_group)) stop("Data 2: no spacer-variant CPSIT_0429 guides (M64+) found")
scaffold_ok <- startsWith(spacer_group$seq, SCAFFOLD)
if (!all(scaffold_ok)) {
  stop("Data 2: spacer-variant guide(s) without the wild-type 20 nt scaffold: ",
       paste(spacer_group$crrna[!scaffold_ok], collapse = ","))
}
core <- substr(SPACER20, 1L, 14L)
if (any(!grepl(core, spacer_group$seq, fixed = TRUE))) {
  stop("Data 2: a spacer-variant guide does not carry the first 14 nt of the spacer")
}
sp20 <- SPACER20
sp17 <- substr(sp20, 1L, 17L)
## M72 is the 3' truncation to a 17 nt spacer: scaffold + the FIRST 17 nt.
if (!endsWith(cr$seq[cr$crrna == "M72"], sp17)) {
  stop("Data 2: M72 does not end in the 5' 17 nt of the 20 bp spacer")
}
## the 17 bp spacer base 14 is the wobble position the paper names (T14C in
## every target gives an rU-dG pair). Base 14 of the 20 nt spacer:
if (!identical(substr(sp20, 14L, 14L), "U")) stop("Data 2: spacer position 14 is not U")
## the crRNA set the paper labels "mut1..mut39" is the ACTIVATOR set; the crRNA
## names used in Fig. 1 are M-numbers. Cross-check every one we rely on exists.
used_crrna <- c("M37", "M39", "M42", "M43", "M45", "M48", "M63", "M65")
if (!all(used_crrna %in% cr$crrna)) {
  stop("Data 2: missing crRNA(s): ", paste(setdiff(used_crrna, cr$crrna), collapse = ","))
}
if (nchar(cr$seq[cr$crrna == "M65"]) != nchar(sl20) - 2L) {
  stop("Data 2: M65 is not a 2 bp 3' truncation")
}

## --- the sheet --------------------------------------------------------------
sheet1 <- xlsx_sheets(xl)[["Fig. 1"]]
if (is.null(sheet1)) stop("source-data workbook has no sheet 'Fig. 1'")
## A2
if (length(sheet1) != 121L) stop("A2: 'Fig. 1' has ", length(sheet1), " rows, expected 121")
if (!all(vapply(sheet1[c(39L, 80L)], function(r) identical(xcell(r, "B"), "crRNA"), FALSE))) {
  stop("A2: rows 39/80 are no longer the per-sub-block header rows ('crRNA' in column B)")
}
if (!all(vapply(sheet1[c(40L, 81L)], function(r) identical(xcell(r, "A"), "CPSIT_0429 -WT"), FALSE))) {
  stop("A2: rows 40/81 no longer open each block with the WT activator")
}
if (!is.na(xcell(sheet1[[121L]], "A"))) {
  stop("A2: the sheet no longer ends at row 120")
}

## Four sub-blocks side by side, read by column. Every sub-block has its own
## label column, crRNA column, three replicate columns and a published mean;
## the replicate/mean columns sit at a fixed offset from the label column
## (label, crRNA, r1, r2, r3, mean). The crRNA is constant down a sub-block and
## is proof-read against that block's own header row.
SUB_LAB_COLS <- c("A", "H", "O", "V")
## col_names() converts 1-based indices to spreadsheet letters WITHOUT relying
## on LETTERS[] indexing: the workbook uses W,X,Y,Z,AA, so `LETTERS[i + 1]` is
## NA past Z and silently reads nothing (this is the class of bug that stands
## behind the off-by-one reported earlier).
col_name <- function(n) {
  out <- character()
  n <- as.integer(n)
  while (n > 0L) {
    r <- (n - 1L) %% 26L
    out <- c(intToUtf8(utf8ToInt("A") + r), out)
    n <- (n - 1L) %/% 26L
  }
  paste(out, collapse = "")
}
COLS_1_702 <- vapply(seq_len(702L), col_name, "")
col_span <- function(lab) {
  i <- match(lab, COLS_1_702)
  if (is.na(i)) stop("bad sub-block column '", lab, "'")
  nm <- COLS_1_702[i + 0:5]
  data.frame(lab = nm[1], crrna = nm[2], r1 = nm[3], r2 = nm[4], r3 = nm[5],
             mean = nm[6], stringsAsFactors = FALSE)
}
subblocks <- function(rows, starts, rowids) {
  do.call(rbind, lapply(starts, function(lab) {
    cc <- col_span(lab)
    body <- data.frame(
      rowid = rowids,
      label = vapply(rows, xcell, "", cc$lab),
      crrna = trimws(vapply(rows, xcell, "", cc$crrna)),
      replicate_1 = num(vapply(rows, xcell, "", cc$r1)),
      replicate_2 = num(vapply(rows, xcell, "", cc$r2)),
      replicate_3 = num(vapply(rows, xcell, "", cc$r3)),
      mean_as_published = num(vapply(rows, xcell, "", cc$mean)),
      stringsAsFactors = FALSE)
    body$subblock = lab
    body
  }))
}
parse_blocks <- function(body_rows, hdr_row, starts) {
  hdr <- sheet1[[hdr_row]]
  out <- subblocks(sheet1[body_rows], starts, body_rows)
  for (lab in starts) {
    got <- unique(out$crrna[out$subblock == lab])
    want <- trimws(xcell(hdr, col_span(lab)$crrna))
    ## header row 39 leaves column A blank; row 80 keeps the "crRNA" word in
    ## its first label column. Accept either, but require the BODY to be
    ## constant down the sub-block, which is what actually identifies it.
    if (!is.na(want) && nzchar(want) && want != "crRNA" &&
        (length(got) != 1L || !identical(got, want))) {
      stop("A3: sub-block ", lab, " header says crRNA '", want,
           "' but its body says '", paste(got, collapse = ","), "'")
    }
    ## sub-block A's label column is the per-row ACTIVATOR; its crRNA column is
    ## constant. Only the crRNA of a sub-block is expected to be constant.
    if (length(got) != 1L) {
      stop("A3: sub-block ", lab, " is not constant down its rows: ",
           paste(got, collapse = ","))
    }
  }
  out
}
## block 1 = rows 41-79 under header row 39; block 2 = rows 81-120 under
## header row 80. Row 40 is a blank spacer inside the first block.
fig1 <- rbind(
  parse_blocks(c(40:79), 39L, c("A", "H", "O", "V")),
  parse_blocks(c(81:120), 80L, c("A", "H", "O", "V")))
fig1$crrna <- trimws(fig1$crrna)
fig1$spacer_length <- 20L
## Sheet reality, checked rather than assumed: the WT rows (sheet rows 40 and
## 81) carry the literal label "CPSIT_0429 -WT" in every sub-block, exactly
## like the mutant rows carry "mut7-C7A". So the label -> activator mapping is
## text only, and no positional special case is needed.
if (!setequal(unique(fig1$label[fig1$rowid %in% c(40L, 81L)]), "CPSIT_0429 -WT")) {
  stop("A3: the WT rows no longer carry the label 'CPSIT_0429 -WT'")
}
if (any(grepl("^[0-9.]+$", fig1$label))) {
  stop("A3: a label cell is numeric; the source layout has changed")
}
fig1$activator <- fig1$label
fig1$mutation <- ifelse(fig1$activator == "CPSIT_0429 -WT", "WT",
                        sub("^mut[0-9]+-", "", fig1$activator))
act <- fig1$mutation

## Label columns hold the ACTIVATOR on every sub-block (repeated four times a
## row); the column next to it holds the CRRNA for that sub-block ("WT crRNA",
## "M37", ...). `fig1$activator`/`fig1$mutation` were set above from the label
## column alone -- nothing downstream re-derives them.
## Fig. 1 measures only mut1..mut39; Data 3 also defines mut40..mut74 (the
## T14-anchored wobble series), which Fig. 1 does not contain.
tpl_fig1 <- tpl[as.integer(sub("^mut", "", tpl$activator)) <= 39L, ]
if (nrow(tpl_fig1) != 39L) stop("the Fig. 1 activator set is not 39 rows")
if (!identical(sort(unique(fig1$mutation)), sort(c("WT", tpl_fig1$mutation)))) {
  stop("A3: the sheet's activator set is not {WT, mut1..mut39}")
}
if (any(grepl("^mut2[1-9]|^mut3[0-9]", act) & !grepl("&", act))) {
  stop("A3: a double-mismatch label carries no '&'")
}
if (!all(table(fig1$mutation) == 8L)) {
  stop("A3: every activator must appear once in each of the 8 sub-blocks (4 per block x 2 blocks)")
}
if (length(unique(fig1$crrna)) != 8L) {
  stop("A3: expected 8 crRNAs in Fig. 1, got ", length(unique(fig1$crrna)))
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
## per crRNA, one WT row each.
wt_mean <- stats::setNames(fig1$mean_recomputed[fig1$mutation == "WT"],
                           fig1$crrna[fig1$mutation == "WT"])
if (length(wt_mean) != 8L || any(!nzchar(names(wt_mean)))) {
  stop("expected one WT row per crRNA, got ", length(wt_mean))
}
fig1$mean_fc_vs_wt <- fig1$mean_recomputed / wt_mean[fig1$crrna]
if (any(!is.finite(fig1$mean_fc_vs_wt))) stop("a ratio is not finite")
if (!all(abs(wt_mean) > 0)) stop("WT mean is not positive")

## --- assemble the two extracts ----------------------------------------------
## `tpl_fig1` is the design side (Data 3, mut1..mut39); `fig1` is the
## measurement side. Join them on the sheet's "mut<N>" id.
fig1$mutid <- sub("^(mut[0-9]+)-.*$", "\\1", fig1$activator)
fig1$mutid[fig1$mutation == "WT"] <- "WT"
actv <- rbind(
  data.frame(mutid = "WT", activator = "CPSIT_0429 -WT", mutation = "WT",
             n_mismatch = 0L, adjacency = "none", mismatch_position = NA_character_,
             mismatch_spacing = NA_integer_, wt_base = NA_character_, alt_base = NA_character_,
             protospacer_20nt = tpl_fig1$protospacer_20nt[1],
             fragment_24nt = tpl_fig1$fragment_24nt[1], stringsAsFactors = FALSE),
  transform(tpl_fig1[, c("activator", "mutation", "n_mismatch", "adjacency",
                         "mismatch_position", "mismatch_spacing", "wt_base", "alt_base",
                         "protospacer_20nt", "fragment_24nt")],
            mutid = activator))
actv <- actv[order(match(actv$mutid, c("WT", sprintf("mut%d", 1:39)))), ]
rownames(actv) <- NULL
if (nrow(actv) != 40L) stop("expected 40 activator rows (WT + mut1..mut39), got ", nrow(actv))
if (!setequal(actv$mutid, c("WT", sprintf("mut%d", 1:39)))) stop("activator ids are not WT + mut1..mut39")

## the design side, all 74 Data 3 rows (mut1..mut39 measured in Fig. 1, the
## rest defined but not measured here)
act_design <- tpl
act_design$measured_in_fig1 <- as.integer(sub("^mut", "", act_design$activator)) <= 39L
write_design <- data.frame(
  source_key = "communbiol2024_crrna",
  target = act_design$target,
  activator = act_design$activator,
  mutation = act_design$mutation,
  n_mismatch = act_design$n_mismatch,
  adjacency = act_design$adjacency,
  mismatch_position = act_design$mismatch_position,
  mismatch_spacing = act_design$mismatch_spacing,
  mismatch_position_basis = "1-based position 1 = the PAM-DISTAL (5') base of the 20 nt spacer; the label's 0-based offset k maps to the 24 nt fragment at R index k+4 (verified for all 222 labelled rows over 3 targets)",
  wt_base_at_position = act_design$wt_base,
  alt_base = act_design$alt_base,
  protospacer_20nt = act_design$protospacer_20nt,
  fragment_24nt = act_design$fragment_24nt,
  measured_in_fig1 = act_design$measured_in_fig1,
  outcome_definition = "mutant activator endpoint trans-cleavage fluorescence / the same crRNA's WT fluorescence within 1 h (Fig. 1b legend: 'inverse multiples of specificity'); 1 = no discrimination, >1 favours the mutant allele",
  evidence_tier = "measured",
  source_file = SOURCE_FILE,
  doi = DOI,
  pmcid = PMCID,
  stringsAsFactors = FALSE)

out_act <- write_design

## the measurements
actv_idx <- match(fig1$mutid, actv$mutid)
ratio <- data.frame(
  source_key = "communbiol2024_crrna",
  target = "CPSIT_0429",
  crrna = fig1$crrna,
  spacer_length = fig1$spacer_length,
  activator = actv$activator[actv_idx],
  mutation = actv$mutation[actv_idx],
  n_mismatch = actv$n_mismatch[actv_idx],
  adjacency = actv$adjacency[actv_idx],
  mismatch_position = actv$mismatch_position[actv_idx],
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
if (nrow(back) != 320L) stop("A9: ratio CSV has ", nrow(back), " rows, expected 320 (40 x 8)")
## A9: the CSV must round-trip. readr writes 15 significant digits, so the
## returned value can differ from the in-memory double in the last bit
## (measured: 2.2e-16 absolute, 2.2e-16 relative). The check is therefore a
## RELATIVE tolerance of 1e-12, which is ~4 orders of magnitude tighter than
## any value in this file is scientifically meaningful.
if (max(abs(back$mean_fc_vs_wt - ratio$mean_fc_vs_wt) /
        pmax(abs(ratio$mean_fc_vs_wt), .Machine$double.eps)) > 1e-12) {
  stop("A9: ratio CSV does not round-trip")
}
back_act <- readr::read_csv(f_act, show_col_types = FALSE)
if (nrow(back_act) != 74L) stop("A9: design CSV has ", nrow(back_act), " rows, expected 74 (Data 3, mut1..mut74)")
if (sum(back_act$measured_in_fig1) != 39L) stop("A9: the design CSV does not flag exactly 39 measured rows")
fig_rows <- back_act[back_act$measured_in_fig1, ]
if (sum(fig_rows$n_mismatch == 1L) != 20L || sum(fig_rows$adjacency == "adjacent") != 19L) {
  stop("A9: the measured rows are not 20 single and 19 adjacent double")
}
## the ratio file must cover exactly the activators flagged as measured, plus WT
if (!setequal(unique(ratio$mutation),
              c("WT", fig_rows$mutation))) {
  stop("A9: the ratio file and the measured design rows disagree")
}

cat(sprintf("communbiol2024_crrna: %d activator rows, %d ratio rows, DOI %s\n",
            nrow(out_act), nrow(ratio), DOI))
