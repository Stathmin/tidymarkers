# =============================================================================
#  dev/as-tm-crosscheck.R  --  the allelic-discrimination layer vs Biopython Tm
# =============================================================================
#  Primer3 does not design allele-specific primers, so the AS layer has no
#  comparator in dev/primer3-crosscheck.R.  This script gives it one: a
#  nearest-neighbour Tm that CAN score a 3'-terminal mismatch.
#
#  WHAT IS HERE
#    * An 8-row panel: four 21-mer KASP-style AS primers (15-nt tail + 6-nt
#      allele region), each placed in two mismatch geometries -- the 3' terminal
#      base (geom 0) and the penultimate base (geom 1).  For each row the REF
#      template is the exact reverse complement of the primer and the ALT
#      template carries the one-base substitution the primer's ALT version
#      targets.
#    * Two independent numbers per row, from different models:
#        - `discrimination` and `dgap`: DECIPHER's PCR efficiency at the
#          package's conditions (4-term formula used by `design_asq_marker()`;
#          see R/marker_assemble.R).
#        - `dtm`: Biopython `Tm_NN` (Allawi & SantaLucia 1997 NN table) with
#          the 3'-terminal mismatch table (SantaLucia & Peyret 2001), matching
#          the same salt and strand concentrations.
#    * Written to results/as-tm-crosscheck.csv (git-ignored).
#
#  WHAT IS NOT HERE
#    * No claim that a Tm gap and an efficiency gap are the same quantity.
#      They are not: `dtm` is a thermodynamic penalty on one duplex, `dgap` is a
#      model of hybridisation plus Taq extension.  The script checks DIRECTION
#      (both disfavour the mismatch) and reports the rank correlation.  It does
#      not assert numeric agreement, and a magnitude would be a property of the
#      two models, not of a real assay.
#    * No orientation claim inside a primer pair.  Sign conventions differ
#      between engines; only `dtm > 0` / `dgap > 0` and the per-primer
#      terminal-vs-penultimate ordering are asserted downstream.
#
#  Requires: OligoArrayAux (for `pcr_efficiency()`) and python3 with Biopython.
#  The companion test skips with a named reason when either is absent.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

as_panel_row <- function(primer, pos_from_3p) {
  tgt <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(primer)))
  L <- nchar(primer)
  # Primer position from the 3' end: 0 = terminal base, 1 = penultimate.
  ppos <- L - pos_from_3p
  base <- substr(primer, ppos, ppos)
  comp <- c(A = "T", C = "G", G = "C", T = "A")[[base]]
  as_alt <- primer
  substr(as_alt, ppos, ppos) <- c(A = "C", C = "A", G = "A", T = "C")[[base]]
  # Target index pairs the primer's offset: 1-based from the target 5' end.
  idx <- pos_from_3p + 1L
  t_alt <- tgt
  substr(t_alt, idx, idx) <- c(A = "C", C = "A", G = "A", T = "C")[[comp]]
  tibble::tibble(
    primer = primer, as_alt = as_alt, t_ref = tgt, t_alt = t_alt,
    geom = as.integer(pos_from_3p)
  )
}

as_panel <- function() {
  primers <- c(
    "GAAGGTGACCAAGTTCATGCT", "CTGAGCGTCAGCTTCAGGAT",
    "AACCGTTCAGGTACCGATGA", "GTTCAGCCAGATTGCCGTAA"
  )
  dplyr::bind_rows(lapply(primers, as_panel_row, pos_from_3p = 0L)) |>
    dplyr::bind_rows(lapply(primers, as_panel_row, pos_from_3p = 1L))
}

panel <- as_panel()
panel$discrimination <- NA_real_
panel$dgap <- NA_real_
for (i in seq_len(nrow(panel))) {
  r <- panel[i, ]
  eff <- pcr_efficiency(
    c(r$primer, r$primer, r$as_alt, r$as_alt),
    c(r$t_ref, r$t_alt, r$t_alt, r$t_ref)
  )
  panel$discrimination[i] <- ((eff[[1L]] - eff[[2L]]) + (eff[[3L]] - eff[[4L]])) / 2
  panel$dgap[i] <- pcr_efficiency(r$primer, r$t_ref) - pcr_efficiency(r$primer, r$t_alt)
}

dir.create("results", showWarnings = FALSE)
panel_csv <- file.path("results", "as-tm-panel.csv")
out_csv <- file.path("results", "as-tm-crosscheck.csv")
utils::write.csv(panel, panel_csv, row.names = FALSE)

py <- file.path("dev", "as-tm-crosscheck.py")
if (!file.exists(py)) {
  stop("run this from the repo root: dev/as-tm-crosscheck.py not found")
}
if (!nzchar(Sys.which("python3"))) {
  stop("python3 not on PATH; the Biopython leg cannot run")
}
status <- system2("python3", c(shQuote(py), shQuote(panel_csv), shQuote(out_csv)))
if (status != 0L) {
  stop("python3 dev/as-tm-crosscheck.py exited ", status)
}

merged <- utils::read.csv(out_csv, stringsAsFactors = FALSE)
merged$primer_short <- substr(merged$primer, nchar(merged$primer) - 5L, nchar(merged$primer))
print(
  merged[, c("primer_short", "geom", "discrimination", "dgap", "tm_ref", "tm_alt", "dtm")],
  row.names = FALSE
)

rho <- stats::cor(merged$dtm, merged$discrimination, method = "spearman")
cat("\nDirection: dtm > 0 in", sum(merged$dtm > 0), "of", nrow(merged), "rows;")
cat(" dgap > 0 in", sum(merged$dgap > 0), "of", nrow(merged), "rows\n")
cat("Spearman rho(dtm, discrimination) =", round(rho, 3), "(reported, not asserted)\n")
cat("Wrote", out_csv, "\n")
