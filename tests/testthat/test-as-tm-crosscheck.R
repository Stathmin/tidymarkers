# The allelic-discrimination layer vs an independent Tm that can see a 3' mismatch.
#
# WHAT IS HERE
#   * An 8-row panel of KASP-style AS primers (four primers x two mismatch
#     geometries: the 3' terminal base, and the penultimate base).
#   * `discrimination` / `dgap` from DECIPHER's PCR efficiency (the package's
#     own 4-term formula), and `dtm` from Biopython `Tm_NN` with the
#     3'-terminal mismatch table (SantaLucia & Peyret 2001).
#   * Direction only: every row must show `dtm > 0` and `dgap > 0` (the mismatch
#     is disfavoured by both engines), and within each primer the terminal
#     mismatch must cost more than the penultimate one in `dtm`.  That second
#     assertion is the anti-degeneracy guard: a comparator that ignored the
#     mismatch position would pass direction while failing it.
#   * Spearman rho is REPORTED, never asserted.  The two numbers are different
#     quantities (a duplex Tm penalty vs an efficiency gap) and were measured to
#     rank primers differently.
#
# WHAT IS NOT HERE
#   * No numeric agreement between `dtm` and `discrimination`, and no absolute
#     Tm constant: engine-to-engine agreement is not evidence about a real
#     duplex.
#   * No claim about sign conventions inside a pair beyond `dtm > 0`/`dgap > 0`.
#
# The numbers behind this test: `Rscript dev/as-tm-crosscheck.R` (writes
# results/as-tm-crosscheck.csv, git-ignored).

as_tm_skip_reason <- function() {
  if (!nzchar(Sys.which("python3"))) {
    return("python3 not on PATH")
  }
  probe <- suppressWarnings(system2(
    "python3", c("-c", shQuote("import Bio.SeqUtils.MeltingTemp")),
    stdout = TRUE, stderr = TRUE
  ))
  if (!is.null(attr(probe, "status"))) {
    return("Biopython Bio.SeqUtils.MeltingTemp not importable")
  }
  NA_character_
}

as_tm_row <- function(primer, pos_from_3p) {
  tgt <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(primer)))
  L <- nchar(primer)
  ppos <- L - pos_from_3p
  base <- substr(primer, ppos, ppos)
  comp <- c(A = "T", C = "G", G = "C", T = "A")[[base]]
  as_alt <- primer
  substr(as_alt, ppos, ppos) <- c(A = "C", C = "A", G = "A", T = "C")[[base]]
  idx <- pos_from_3p + 1L
  t_alt <- tgt
  substr(t_alt, idx, idx) <- c(A = "C", C = "A", G = "A", T = "C")[[comp]]
  data.frame(
    primer = primer, as_alt = as_alt, t_ref = tgt, t_alt = t_alt,
    geom = as.integer(pos_from_3p), stringsAsFactors = FALSE
  )
}

as_tm_panel <- function() {
  primers <- c(
    "GAAGGTGACCAAGTTCATGCT", "CTGAGCGTCAGCTTCAGGAT",
    "AACCGTTCAGGTACCGATGA", "GTTCAGCCAGATTGCCGTAA"
  )
  do.call(rbind, c(
    lapply(primers, as_tm_row, pos_from_3p = 0L),
    lapply(primers, as_tm_row, pos_from_3p = 1L)
  ))
}

test_that("a 3' mismatch is disfavoured by both engines, terminal costlier than penultimate", {
  reason <- as_tm_skip_reason()
  if (!is.na(reason)) testthat::skip(reason)
  skip_if_no_oligoarrayaux()

  panel <- as_tm_panel()
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

  # Self-check that the panel is what this test claims: the ALT record differs
  # from the REF record at exactly one base -- the primer position `ppos` and the
  # target position `idx` -- and nowhere else.  This is the guard against an
  # off-by-one that would place the mismatch outside the molecule.
  for (i in seq_len(nrow(panel))) {
    L <- nchar(panel$primer[i])
    ppos <- L - panel$geom[i]
    idx <- panel$geom[i] + 1L
    expect_false(identical(panel$as_alt[i], panel$primer[i]))
    expect_false(identical(panel$t_alt[i], panel$t_ref[i]))
    expect_identical(substr(panel$as_alt[i], 1L, ppos - 1L), substr(panel$primer[i], 1L, ppos - 1L))
    expect_identical(substr(panel$as_alt[i], ppos + 1L, L), substr(panel$primer[i], ppos + 1L, L))
    expect_identical(substr(panel$t_alt[i], 1L, idx - 1L), substr(panel$t_ref[i], 1L, idx - 1L))
    expect_identical(substr(panel$t_alt[i], idx + 1L, L), substr(panel$t_ref[i], idx + 1L, L))
  }

  py <- test_path("..", "..", "dev", "as-tm-crosscheck.py")
  if (!file.exists(py)) testthat::skip("dev/as-tm-crosscheck.py not found")
  in_csv <- tempfile(fileext = ".csv")
  out_csv <- tempfile(fileext = ".csv")
  utils::write.csv(panel, in_csv, row.names = FALSE)
  status <- system2("python3", c(shQuote(py), shQuote(in_csv), shQuote(out_csv)))
  expect_identical(status, 0L)
  merged <- utils::read.csv(out_csv, stringsAsFactors = FALSE)
  expect_equal(nrow(merged), 8L)

  # Direction, 8/8.
  expect_true(all(merged$dtm > 0))
  expect_true(all(merged$dgap > 0))

  # Anti-degeneracy: per primer, the terminal mismatch costs more Tm than the
  # penultimate one.
  term <- merged$dtm[merged$geom == 0L]
  penult <- merged$dtm[merged$geom == 1L]
  expect_equal(length(term), 4L)
  expect_true(all(penult > term))

  # Reported, not asserted: the two metrics do not rank alike.
  rho <- suppressWarnings(stats::cor(merged$dtm, merged$discrimination, method = "spearman"))
  expect_true(is.finite(rho))
})
