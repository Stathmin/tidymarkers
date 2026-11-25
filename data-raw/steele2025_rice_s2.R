# =============================================================================
#  data-raw/steele2025_rice_s2.R  --  Y6: per-marker rice KASP call rates
# =============================================================================
#  Steele K, Quinton-Tulloch M, Vyas D, Witcombe J. "Thousands of trait-specific
#  KASP markers designed for diverse breeding applications in rice", G3
#  Genes|Genomes|Genetics 15(1):jkae251, 10.1093/g3journal/jkae251,
#  PMID 39486028, PMCID PMC11708223 (CC BY).
#
#  This is the only source in the registry that reaches a per-marker reading of
#  KASP call success on a real panel: 3,366 markers x 178 rice varieties.
#  Supplementary File S2 carries two sheets:
#
#    Table I  -- one row per marker: BU assay id, Shuhui498 position, and the
#                `Reference allele call summary` / `Alternative allele call
#                summary` columns.  A cell is "<homozygous genotype>: <n>"
#                (e.g. "G:G: 145", "TAAT:TAAT: 30", "-:-: 12", "?: 147").
#                IMPORTANT: these columns count HOMOZYGOUS CALLS ONLY -- het
#                cells and no-calls are simply absent from them, so
#                ref + alt != 178 and must never be asserted as such.  A whole
#                column prints "?: N" when that allele was never called.
#    Table J  -- 178 variety rows (6..183) x 3,366 marker columns (B..DYM);
#                rows 1-5 are the marker index, Indica position, working id,
#                alternative id, and the "Variety" header.  Cell vocabulary:
#                "A:A"/"A:B"/"-:-" are genotype calls, "?" is the no-call, and
#                the literal strings "Uncallable" / "Bad" are annotations, not
#                calls.
#
#  What the paper does NOT report and this file does NOT invent:
#    * Table I does not give a per-marker call total; the count is rebuilt from
#      Table J and cross-checked against Table I's homozygous tallies.
#    * The text says 3,366 markers "gave successful genotype calls in >90% of
#      samples" and that 635 designs (15.9 %) failed in all samples.  The
#      shipped matrix shows 91 marker columns at <= 90 % called and exactly one
#      all-no-call column.  We record the matrix and name the discrepancy in
#      ROADMAP Y6; we do not adjudicate the paper's sentence.
#    * The `KASP Alternative ID` row repairs 4 of the 7 id disagreements
#      between Table I and Table J; it is not an id map.  We align by column
#      index (verified: positions agree at all 3,366 indexes).
#
#  Route (the only one that worked, 2026-11-25):
#    curl -L --max-time 900 -o File_S2_G3-2024-405461.xlsx \
#      "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11708223/supplementaryFiles?includeInlineImage=true"
#  then unzip and keep File_S2_G3-2024-405461.xlsx in .tmp/steele/suppl/s2/
#  (or pass --dir=<dir containing it>).
#
#  Run:  Rscript data-raw/steele2025_rice_s2.R [--dir=PATH]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

source(file.path(root, "data-raw", "lib", "xlsx-read.R"))
suppressPackageStartupMessages(library(xml2))

args <- commandArgs(trailingOnly = TRUE)
dir_arg <- grep("^--dir=", args, value = TRUE)
dir <- if (length(dir_arg)) sub("^--dir=", "", dir_arg[[1L]]) else ".tmp/steele/suppl/s2"
xlsx <- file.path(dir, "File_S2_G3-2024-405461.xlsx")
if (!file.exists(xlsx)) {
  stop(
    "File S2 not found at ", xlsx, "\n",
    "  Download it with:\n",
    "  curl -L --max-time 900 -o File_S2_G3-2024-405461.xlsx \\\n",
    "    \"https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11708223/supplementaryFiles?includeInlineImage=true\"\n",
    "  then unzip and keep File_S2_G3-2024-405461.xlsx in ", dir, "\n",
    call. = FALSE
  )
}

sh <- xlsx_sheets(xlsx)
stopifnot(c("Table I", "Table J") %in% names(sh))

# --- Table I: per-marker allele labels -------------------------------------
# "<geno>: <n>" -> geno; a "?" geno means the column was never called.
parse_allele <- function(x) {
  if (is.na(x)) {
    return(list(geno = NA_character_, n = 0L))
  }
  geno <- sub(":[^:]*$", "", x)
  n <- as.integer(sub(".*:", "", x))
  if (!is.na(geno) && geno == "?") {
    return(list(geno = NA_character_, n = n))
  }
  list(geno = geno, n = n)
}

ti <- sh[["Table I"]]
n_marker <- length(ti) - 1L
assay_id <- character(n_marker)
position <- character(n_marker)
ref_allele <- character(n_marker)
alt_allele <- character(n_marker)
ref_i <- integer(n_marker)
alt_i <- integer(n_marker)
ref_q <- logical(n_marker)
alt_q <- logical(n_marker)
for (r in seq_len(n_marker)) {
  cells <- ti[[r + 1L]]
  assay_id[r] <- xcell(cells, "A")
  position[r] <- xcell(cells, "B")
  cr <- xcell(cells, "C")
  ca <- xcell(cells, "D")
  pr <- parse_allele(cr)
  pa <- parse_allele(ca)
  ref_allele[r] <- pr$geno
  alt_allele[r] <- pa$geno
  ref_i[r] <- pr$n
  alt_i[r] <- pa$n
  ref_q[r] <- !is.na(cr) && sub(":[^:]*$", "", cr) == "?"
  alt_q[r] <- !is.na(ca) && sub(":[^:]*$", "", ca) == "?"
}
q_markers <- ref_q | alt_q

# --- Table J: the 178 x 3,366 call matrix ----------------------------------
# Column letters -> 1-based index; Table J puts marker k in column k + 1.
col_index <- function(x) {
  vapply(strsplit(x, "", fixed = TRUE), function(ch) {
    Reduce(function(acc, cc) acc * 26L + match(cc, LETTERS), ch, 0L)
  }, integer(1))
}

tj <- sh[["Table J"]]
n_row <- length(tj)
n_variety <- n_row - 5L
if (n_variety != 178L) {
  stop("expected 178 variety rows in Table J, found ", n_variety, call. = FALSE)
}

j_number <- as.integer(unlist(tj[[1L]][-1L], use.names = FALSE))
j_pos <- unlist(tj[[2L]][-1L], use.names = FALSE)
j_work <- unlist(tj[[3L]][-1L], use.names = FALSE)
if (length(j_number) != n_marker || length(j_pos) != n_marker) {
  stop("Table J header does not have one column per Table I marker", call. = FALSE)
}

# Each variety row is a named vector whose first element is the variety label,
# so drop it positionally and map the rest onto marker 1..n_marker.
stopifnot(identical(col_index(names(tj[[6L]]))[-1L] - 1L, seq_len(n_marker)))

cnt <- matrix(0L,
  nrow = n_marker, ncol = 5L,
  dimnames = list(NULL, c("ref", "alt", "het", "missing", "other"))
)
for (r in seq_len(n_variety)) {
  cells <- tj[[r + 5L]]
  v <- rep(NA_character_, n_marker)
  v[col_index(names(cells))[-1L] - 1L] <- unlist(cells, use.names = FALSE)[-1L]
  miss <- !is.na(v) & v == "?"
  r_ref <- !is.na(v) & !is.na(ref_allele) & v == ref_allele
  r_alt <- !is.na(v) & !is.na(alt_allele) & v == alt_allele
  two <- !is.na(v) & grepl(":", v, fixed = TRUE)
  s1 <- sub(":.*$", "", v)
  s2 <- sub("^.*:", "", v)
  r_het <- two & !miss & !r_ref & !r_alt & s1 != "?" & s2 != "?" & s1 != s2
  cnt[, "ref"] <- cnt[, "ref"] + as.integer(r_ref)
  cnt[, "alt"] <- cnt[, "alt"] + as.integer(r_alt)
  cnt[, "het"] <- cnt[, "het"] + as.integer(r_het)
  cnt[, "missing"] <- cnt[, "missing"] + as.integer(miss)
  cnt[, "other"] <- cnt[, "other"] + as.integer(!is.na(v) & !miss & !r_ref & !r_alt & !r_het)
}

# Position strings are the paper's own allele-marker coordinate: Table I and
# Table J agree at every index.  Split once, from Table I.
chrom <- sub(":.*$", "", position)
pos <- as.integer(sub("^.*:", "", position))
stopifnot(!anyNA(pos), all(nzchar(chrom)))

n_called <- 178L - cnt[, "missing"]
call_rate <- n_called / 178
passes_90 <- n_called >= 161L

markers <- data.frame(
  assay_id = assay_id,
  working_id = j_work,
  chrom = chrom,
  pos = pos,
  ref_allele = ref_allele,
  alt_allele = alt_allele,
  ref_calls = cnt[, "ref"],
  alt_calls = cnt[, "alt"],
  het_calls = cnt[, "het"],
  missing_calls = cnt[, "missing"],
  other_calls = cnt[, "other"],
  n_called = n_called,
  call_rate = round(call_rate, 4L),
  passes_90 = passes_90,
  stringsAsFactors = FALSE
)

# --- the seven id disagreements, recorded rather than repaired -------------
id_mismatch_at <- c(506L, 1541L, 1853L, 1968L, 2385L, 2449L, 3329L)
id_mismatch <- markers$assay_id != markers$working_id
stopifnot(
  # shape
  n_marker == 3366L,
  n_variety == 178L,
  nrow(markers) == 3366L,
  !anyDuplicated(markers$assay_id),
  !anyDuplicated(paste(markers$chrom, markers$pos)),
  identical(j_number, seq_len(n_marker)),
  # the two sheets carry the same coordinate at every index
  identical(unname(j_pos), position),
  # Table J conserves the panel: five classes sum to 178 for every marker
  all(rowSums(cnt) == 178L),
  # measured totals (python/openpyxl and this reader agree)
  sum(cnt) == 178L * 3366L,
  sum(cnt[, "missing"]) == 9810L,
  sum(cnt[, "ref"]) == 417319L,
  sum(cnt[, "alt"]) == 168706L,
  sum(cnt[, "het"]) == 2004L,
  sum(cnt[, "other"]) == 1309L,
  sum(markers$n_called) == 178L * 3366L - 9810L,
  # Table I's homozygous tallies ARE the Table J exact-genotype tallies
  all(ref_i[!ref_q] == cnt[!ref_q, "ref"]),
  all(alt_i[!alt_q] == cnt[!alt_q, "alt"]),
  # ... and a "?: N" column carries that allele's no-call count
  all(ref_i[ref_q] == cnt[ref_q, "missing"]),
  all(alt_i[alt_q] == cnt[alt_q, "missing"]),
  sum(q_markers) == 51L,
  # empty allele columns: 6 reference, 50 alternative (45 print "?: N", 5 are
  # blank), and exactly one marker with no call printed for either allele
  sum(is.na(ref_allele)) == 6L,
  sum(is.na(alt_allele)) == 50L,
  sum(is.na(ref_allele) & is.na(alt_allele)) == 1L,
  markers$assay_id[is.na(ref_allele) & is.na(alt_allele)] == "RM214_SNP_nn_2",
  # "-:-" is a real homozygous deletion call, not a no-call
  sum(grepl("-", ref_allele, fixed = TRUE) | grepl("-", alt_allele, fixed = TRUE)) == 21L,
  # 11 markers carry a multi-base allele; "-:-" is a single-base deletion
  sum((!is.na(ref_allele) & nchar(ref_allele) > 3L) |
    (!is.na(alt_allele) & nchar(alt_allele) > 3L)) == 11L,
  # call-rate distribution, as measured
  sum(markers$n_called <= 160L) == 91L,
  sum(markers$passes_90) == 3275L,
  sum(markers$n_called == 0L) == 1L,
  markers$assay_id[markers$n_called == 0L] == "RM214_SNP_nn_2",
  # the 7 id disagreements are exactly the 7 measured, and nothing else
  identical(sum(id_mismatch), 7L),
  identical(which(id_mismatch), id_mismatch_at),
  identical(
    markers$assay_id[id_mismatch],
    c("id2001301", "id6001041", "id6016635", "gs_id1201_fn",
      "FGR8_SNP1", "id8007144", "id12007741")
  ),
  identical(
    markers$working_id[id_mismatch],
    c("RM211_SNP_nn_5", "RM19311_SNP_ff_1", "OsR498G0613266000_SNP_ff_1",
      "Xa8_NBSLRR_SNP3", "OS08G0424500_SNP_nn_3", "RM1345_SNP_nn_2",
      "OsR498G1221875800_SNP_fn_3")
  )
)

out <- "inst/extdata/steele2025_rice_kasp_calls.csv"
utils::write.csv(markers, out, row.names = FALSE)
cat(
  "steele2025_rice_s2:", nrow(markers), "markers x", n_variety, "varieties;",
  "no-call", sum(markers$missing_calls), "/", 178L * 3366L,
  sprintf("(%.2f %%)", 100 * sum(markers$missing_calls) / (178 * 3366)),
  ";", sum(markers$passes_90), "markers at >90 % called;",
  sum(markers$n_called == 0L), "all-no-call;",
  "id disagreements", sum(id_mismatch), "\n"
)
