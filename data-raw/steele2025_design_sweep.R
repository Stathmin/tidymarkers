# =============================================================================
#  data-raw/steele2025_design_sweep.R  --  ambiguity budget vs. assay outcomes
# =============================================================================
#  Steele K, Quinton-Tulloch M, Vyas D, Witcombe J. "Thousands of trait-specific
#  KASP markers designed for diverse breeding applications in rice", G3
#  15(1):jkae251, 10.1093/g3journal/jkae251, PMCID PMC11708223 (CC BY).
#
#  Table 3 ("Effect of reducing the cutoff for the maximum number of bases
#  permitted in the flanking sequences") is the paper's only sweep of a DESIGN
#  parameter against wet-lab outcomes: it re-runs the design search with the
#  flanking-ambiguity budget cut from 4 to 0 and reports how many designs
#  survive and how many of the tested ones work.
#
#  What the sweep says, and why it matters here: tightening the ambiguity
#  cutoff costs 84 % of the design space (1,392,359 -> 223,783) and moves the
#  median spacing between candidate designs 8.5x (68 -> 580 bp), yet the
#  successful-assay rate stays inside 4.2 points (83.76-87.92 %).  Design count
#  and wet-lab success are close to decoupled across this range.  That is a
#  check on any prior that treats "more candidate designs" as "better assays".
#
#  Values are transcribed from Table 3, not derived; the stopifnot block is the
#  transcription guard (the columns are monotone in the cutoff except for the
#  two outcome rates, which are not).
#
#  Run:  Rscript data-raw/steele2025_design_sweep.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

sweep <- data.frame(
  max_ambiguous_bases = c(4L, 3L, 2L, 1L, 0L),
  designs = c(1392359L, 1135654L, 845121L, 532572L, 223783L),
  median_spacing_bp = c(68L, 93L, 144L, 242L, 580L),
  markers_tested = c(3542L, 3121L, 2527L, 1787L, 828L),
  pct_producing_results = c(85.09, 86.48, 87.65, 89.20, 88.53),
  pct_successful = c(83.76, 85.16, 86.36, 87.92, 87.44),
  doi = "10.1093/g3journal/jkae251",
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(sweep) == 5L,
  identical(sweep$max_ambiguous_bases, 4:0),
  # design space and spacing are monotone in the cutoff
  all(diff(sweep$designs) < 0),
  all(diff(sweep$median_spacing_bp) > 0),
  # tightest cutoff keeps 16.1 % of the designs and loses 84 % of them
  sweep$designs[[5L]] / sweep$designs[[1L]] > 0.16,
  sweep$designs[[5L]] / sweep$designs[[1L]] < 0.17,
  # the outcome rates barely move: the whole point of the table
  min(sweep$pct_successful) >= 83.7,
  max(sweep$pct_successful) <= 87.95,
  max(sweep$pct_successful) - min(sweep$pct_successful) < 5,
  all(sweep$pct_producing_results >= sweep$pct_successful)
)

utils::write.csv(sweep, "inst/extdata/steele2025_design_sweep.csv", row.names = FALSE)
cat(
  "steele2025_design_sweep:", nrow(sweep), "cutoffs;",
  "designs", sweep$designs[[1L]], "->", sweep$designs[[nrow(sweep)]], ";",
  sprintf("success %.2f-%.2f %%\n", min(sweep$pct_successful), max(sweep$pct_successful))
)
