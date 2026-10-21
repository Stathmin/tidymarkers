# ARMS/AS-PCR family, batch 2: sorghum ARMS mismatch positions and Cq outcomes.
#
# Zhu et al. 2021 (Foods 10:2218) is the only source in this batch whose
# mismatch POSITIONS can be verified from sequence rather than prose. These tests
# re-derive the positions from the sequences stored in the CSV, so the file is
# self-checking even without the PDF.

sorghum_csv <- function() {
  path <- test_path("../../inst/extdata/arms_sorghum2021_outcomes.csv")
  skip_if_not(file.exists(path), "sorghum outcomes CSV not built")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# Position of each mismatch relative to the 3' end (-1 = last base).
mismatch_positions <- function(seq, ref) {
  a <- strsplit(seq, "")[[1L]]
  b <- strsplit(ref, "")[[1L]]
  n <- length(a)
  if (n != length(b)) stop("length mismatch")
  idx <- which(a != b)
  as.integer(sort(-(n - idx + 1L)))
}

test_that("mismatch positions recomputed from the sequences match the CSV", {
  d <- sorghum_csv()
  ref <- d$sequence[d$primer_group == "wxc_f0"]
  expect_length(ref, 1L)
  for (i in seq_len(nrow(d))) {
    got <- mismatch_positions(d$sequence[i], ref)
    stated <- if (d$mismatch_positions[i] == "0") integer(0) else
      as.integer(strsplit(d$mismatch_positions[i], ";")[[1L]])
    expect_equal(got, stated, info = d$primer_group[i])
    expect_equal(length(got), d$n_mismatch[i], info = d$primer_group[i])
  }
})

test_that("the primer the authors call specific carries its mismatch only at -2", {
  d <- sorghum_csv()
  wxc4 <- d[d$primer_group == "wxc_f4", ]
  expect_equal(wxc4$mismatch_positions, "-2")
  expect_equal(wxc4$n_mismatch, 1L)
})

test_that("the primer called 'fourth' also carries a -1 mismatch", {
  d <- sorghum_csv()
  wxt4 <- d[d$primer_group == "wxt_f4", ]
  # The paper's prose says "fourth from the 3' end"; the sequence says -1 AND -4.
  expect_equal(wxt4$n_mismatch, 2L)
  expect_setequal(strsplit(wxt4$mismatch_positions, ";")[[1L]], c("-4", "-1"))
})

test_that("Cq outcomes exist only for the two scored primers and are consistent", {
  d <- sorghum_csv()
  scored <- d[d$outcome_available, ]
  expect_setequal(scored$primer_group, c("wxc_f4", "wxt_f4"))
  expect_false(any(is.na(scored$delta_ct)))
  # the mismatched primer must be later (larger Cq) than the matched one
  expect_true(all(scored$delta_ct > 0))
  expect_equal(scored$delta_ct, scored$ct_mismatched - scored$ct_matched)
})

test_that("unscored primers carry no outcome values", {
  d <- sorghum_csv()
  unscored <- d[!d$outcome_available, ]
  expect_true(all(is.na(unscored$delta_ct)))
  expect_true(all(is.na(unscored$ct_matched)))
})

test_that("the CSV still matches the fetched Table S1 when it is present", {
  lit <- test_path("../../.tmp/lit/sorghum_arms2021/supplements/s001",
                   "foods-1353494-supplementary.pdf")
  skip_if_not(file.exists(lit), "sorghum supplement not fetched here")
  tf <- tempfile(fileext = ".txt")
  system2("pdftotext", c("-layout", shQuote(lit), shQuote(tf)),
          stdout = TRUE, stderr = TRUE)
  skip_if_not(file.exists(tf), "pdftotext unavailable")
  lines <- readLines(tf, warn = FALSE)
  unlink(tf)
  pat <- "^[[:space:]]*([Ww]x[ct][[:space:]]*F[0-4])[[:space:]]+([ACGT]{20})[[:space:]]*$"
  hit <- Filter(length, regmatches(lines, regexec(pat, lines)))
  seqs <- stats::setNames(
    vapply(hit, `[[`, character(1), 3L),
    gsub("[[:space:]]+", "_", tolower(vapply(hit, `[[`, character(1), 2L))))
  d <- sorghum_csv()
  expect_equal(nrow(d), length(seqs))
  expect_equal(seqs[d$primer_group], d$sequence, ignore_attr = TRUE)
})
