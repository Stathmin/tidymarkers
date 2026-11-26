# =============================================================================
#  data-raw/ye2001_arms_primer_set.R  --  the Ye 2001 tetra-primer ARMS-PCR
#                                         set: 12 primers, 3 polymorphisms
# =============================================================================
#  S. Ye, S. Dhillon, X. Ke, A. R. Collins, I. N. M. Day, "An efficient
#  procedure for genotyping single nucleotide polymorphisms", Nucleic Acids
#  Research 29(17):e88 (2001), Table 2.  PMC55900.
#
#  WHY THIS FILE EXISTS
#  `mismatch_discrimination.csv` already carries the paper's *rule* (a strong
#  3'-terminal mismatch tolerates a weak second mismatch, and vice versa).
#  This table is the *worked example* behind that rule: twelve real primers,
#  each with the author's own Tm, the deliberate second mismatch two bases from
#  the 3' terminus in both inner primers, and the amplicon each pair makes.
#  It is the package's only source of an author-printed Tm for a primer whose
#  allele-specific base and second mismatch are both known, so the NN Tm model
#  can be scored against numbers it did not produce.
#
#  TABLE SHAPE (Table 2; 13 <tr>, 7 <td> each, no rowspan)
#    row 0   Genetic polymorphism | Primer sequence | Tm | Touchdown PCR
#            | Annealing temperature | Mg2+ | Amplicon size
#    rows 1-4  TNF -308G>A   (inner A, inner G, outer F, outer R)
#    rows 5-8  IL6 -174G>C   (inner G, inner C, outer F, outer R)
#    rows 9-12 AGTR1 1166A>C (inner A, inner C, outer F, outer R)
#  Only the first row of each assay carries the touchdown / annealing / Mg
#  cells; the other three leave them blank and inherit downwards.
#
#  DERIVATION RULES
#   * cells contain nested <hr> and <em>, so raw text arrives with newlines:
#     collapse all whitespace before matching.
#   * `-` and `>` are the en dash U+2013 and the rightwards arrow U+2192; the
#     prime in `5'` is U+2032.  All three are normalised or stripped.
#   * `68°C` -> 68, `2.5 mM` -> 2.5, `154 bp (A allele)` -> 154 + `A allele`.
#   * `323 bp (from two outer primers)` names the product `outer pair`.
#   * outer primers carry no allele in the paper: `allele` is NA for them.
#   * the paper's Tm column is carried verbatim, not recomputed.  Methods
#     attribute the values to the nearest-neighbour parameters of Breslauer
#     et al. (ref 13) with the formula of Rychlik & Rhoads (ref 14).
#
#  FETCH (verified route, 2026-11-26)
#    curl -sS -L -m 120 -o .tmp/lit/ye2001/page.html \
#      'https://pmc.ncbi.nlm.nih.gov/articles/PMC55900/'
#    (page.html is 140 989 bytes; the page holds exactly two <table>s and the
#    target is the one whose text contains "Genetic polymorphism")
#
#  Run:  Rscript data-raw/ye2001_arms_primer_set.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

page_path <- Sys.getenv(
  "YE2001_HTML",
  file.path(root, ".tmp", "lit", "ye2001", "page.html")
)
if (!file.exists(page_path)) {
  stop("missing ", page_path, "; see the FETCH block at the top of this script")
}

doc <- xml2::read_html(page_path)
tabs <- xml2::xml_find_all(doc, "//table")
if (length(tabs) != 2L) {
  stop("expected 2 <table>s on the page, got ", length(tabs))
}
hit <- tabs[vapply(
  tabs,
  function(t) grepl("Genetic polymorphism", xml2::xml_text(t), fixed = TRUE),
  logical(1)
)]
if (length(hit) != 1L) {
  stop("Table 2 not found exactly once; got ", length(hit))
}

flat <- function(x) {
  # blank cells are U+00A0 (no-break space); [[:space:]] does not match it
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", x))
}
norm_bases <- function(x) {
  x <- gsub("\u2013", "-", x, fixed = TRUE)
  gsub("\u2192", ">", x, fixed = TRUE)
}
rows <- lapply(
  xml2::xml_find_all(hit[[1L]], ".//tr"),
  function(r) norm_bases(flat(xml2::xml_text(
    xml2::xml_find_all(r, ".//td|.//th")
  )))
)
if (length(rows) != 13L) {
  stop("expected 13 rows (1 header + 12 primers), got ", length(rows))
}
if (!all(lengths(rows) == 7L)) {
  stop("every row must have 7 cells; got ", paste(lengths(rows), collapse = ","))
}

expected_header <- c(
  "Genetic polymorphism", "Primer sequence", "Tm", "Touchdown PCR",
  "Annealing temperature", "Mg2+", "Amplicon size"
)
if (!identical(rows[[1L]], expected_header)) {
  stop("unexpected header: ", paste(rows[[1L]], collapse = " | "))
}
data_rows <- rows[-1L]

role_of <- function(cell) {
  roles <- c(
    "Forward inner primer" = "forward_inner",
    "Reverse inner primer" = "reverse_inner",
    "Forward outer primer" = "forward_outer",
    "Reverse outer primer" = "reverse_outer"
  )
  hit_role <- names(roles)[startsWith(cell, names(roles))]
  if (length(hit_role) != 1L) {
    stop("unknown primer role in cell: ", cell)
  }
  unname(roles[[hit_role]])
}

seq_of <- function(cell) {
  m <- regmatches(cell, regexpr("5[\u2032']-[ACGT]+", cell))
  if (length(m) != 1L) {
    stop("no `5'-SEQ` in cell: ", cell)
  }
  sub("^5[\u2032']-", "", m)
}

allele_of <- function(cell) {
  m <- regmatches(cell, regexpr("\\([ACGT] allele\\)", cell))
  if (length(m) != 1L) {
    return(NA_character_)
  }
  sub("^\\(([ACGT]) allele\\)$", "\\1", m)
}

num_of <- function(cell) {
  m <- regmatches(cell, regexpr("[0-9.]+", cell))
  if (length(m) != 1L) {
    return(NA_real_)
  }
  as.numeric(m)
}

amplicon_of <- function(cell) {
  if (!nzchar(cell)) {
    return(list(bp = NA_real_, product = NA_character_))
  }
  m <- regmatches(cell, regexpr("^[0-9]+ bp \\((.+)\\)$", cell, perl = TRUE))
  if (length(m) != 1L) {
    stop("unparsable amplicon cell: ", cell)
  }
  label <- sub("^[0-9]+ bp \\((.+)\\)$", "\\1", cell)
  list(
    bp = as.numeric(sub(" .*$", "", cell)),
    # verbatim product label; "from two outer primers" is not an allele
    product = if (identical(label, "from two outer primers")) {
      "outer pair"
    } else {
      label
    }
  )
}

# the paper prints conditions on the first row of each assay only
# `td` is TRUE/FALSE/NA exactly as the printed cell reads
parse_touchdown <- function(cell) {
  if (identical(cell, "Yes")) TRUE else if (identical(cell, "No")) FALSE else NA
}

out <- vector("list", length(data_rows))
group <- NULL
cond <- list(touchdown = NA, annealing_c = NA_real_, mgcl2_mm = NA_real_)
for (i in seq_along(data_rows)) {
  v <- data_rows[[i]]
  if (nzchar(v[[1L]])) {
    group <- v[[1L]]
    cond <- list(
      touchdown = parse_touchdown(v[[4L]]),
      annealing_c = num_of(v[[5L]]),
      mgcl2_mm = num_of(v[[6L]])
    )
  }
  amp <- amplicon_of(v[[7L]])
  out[[i]] <- data.frame(
    polymorphism = group,
    primer_role = role_of(v[[2L]]),
    allele = allele_of(v[[2L]]),
    sequence = seq_of(v[[2L]]),
    tm_c = num_of(v[[3L]]),
    touchdown = cond$touchdown,
    annealing_c = cond$annealing_c,
    mgcl2_mm = cond$mgcl2_mm,
    amplicon_bp = amp$bp,
    product = amp$product,
    stringsAsFactors = FALSE
  )
}
out <- do.call(rbind, out)

expected_polymorphisms <- c(
  "TNF -308G>A", "TNF -308G>A", "TNF -308G>A", "TNF -308G>A",
  "IL6 -174G>C", "IL6 -174G>C", "IL6 -174G>C", "IL6 -174G>C",
  "AGTR1 1166A>C", "AGTR1 1166A>C", "AGTR1 1166A>C", "AGTR1 1166A>C"
)
expected_roles <- rep(
  c("forward_inner", "reverse_inner", "forward_outer", "reverse_outer"),
  3L
)

stopifnot(
  nrow(out) == 12L,
  ncol(out) == 10L,
  identical(out$polymorphism, expected_polymorphisms),
  identical(out$primer_role, expected_roles),
  !anyDuplicated(out$sequence),
  all(grepl("^[ACGT]{20,40}$", out$sequence)),
  # the paper's own Tm column, exactly these six values
  identical(sort(unique(out$tm_c)), c(62, 64, 65, 67, 68, 72)),
  # the paper's inner-primer design rule: 26 nt or longer
  min(nchar(out$sequence[grepl("inner", out$primer_role)])) >= 26L,
  # outer primers carry no allele; inner ones do
  all(is.na(out$allele[grepl("outer", out$primer_role)])),
  all(!is.na(out$allele[grepl("inner", out$primer_role)])),
  # per-assay conditions, in TNF/IL6/AGTR1 order
  identical(unique(out[c("polymorphism", "touchdown", "annealing_c", "mgcl2_mm")])$touchdown,
            c(TRUE, TRUE, FALSE)),
  identical(unique(out[c("polymorphism", "touchdown", "annealing_c", "mgcl2_mm")])$annealing_c,
            c(63, 63, 58)),
  identical(unique(out[c("polymorphism", "touchdown", "annealing_c", "mgcl2_mm")])$mgcl2_mm,
            c(2.5, 2.5, 3.5)),
  # per-assay amplicons in printed order (outer reverse has none)
  identical(out$amplicon_bp, c(154, 224, 323, NA, 205, 176, 326, NA, 251, 224, 421, NA)),
  identical(out$product, c(
    "A allele", "G allele", "outer pair", NA,
    "G allele", "C allele", "outer pair", NA,
    "A allele", "C allele", "outer pair", NA
  )),
  # the two inner amplicons of an assay must differ, or the assay cannot
  # genotype: that is the whole point of a tetra-primer design
  all(vapply(split(out$amplicon_bp, out$polymorphism), function(bp) {
    inner <- bp[1:2]
    !anyNA(inner) && length(unique(inner)) == 2L
  }, logical(1)))
)

readr::write_csv(out, "inst/extdata/ye2001_arms_primer_set.csv")
cat(
  "ye2001_arms_primer_set:", nrow(out), "primers;",
  length(unique(out$polymorphism)), "polymorphisms;",
  sprintf("Tm %.0f-%.0f, inner length %d-%d\n",
          min(out$tm_c), max(out$tm_c),
          min(nchar(out$sequence[grepl("inner", out$primer_role)])),
          max(nchar(out$sequence[grepl("inner", out$primer_role)])))
)
