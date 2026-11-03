# Banned wording in the vignettes (ROADMAP Z11.2): a term that reads as a fact
# but is context-dependent, or a label that drops its qualifier.  Text-only: the
# sources are scanned in a checkout, no package loading and no rendering.

vig_lines <- function() {
  dir <- testthat::test_path("..", "..", "vignettes")
  files <- sort(list.files(dir, pattern = "[.](Rmd|md|qmd)$", full.names = TRUE))
  # the vignette sources are not shipped with the installed package
  testthat::skip_if_not(length(files) > 0L, "vignette sources not present (installed check)")
  unlist(lapply(files, function(f) {
    x <- readLines(f, warn = FALSE)
    paste0(basename(f), ":", seq_along(x), ": ", x)
  }))
}

# the sentence is the unit: a line can hold several, and a qualifier in the next
# sentence must not excuse a fact in this one
sentences <- function(lines) {
  body <- sub("^[^:]*:[0-9]+: ", "", lines)
  s <- trimws(unlist(strsplit(body, "(?<=[.!?]) |[.]{3,}", perl = TRUE)))
  s[nzchar(s)]
}

# "Taq does not extend a mismatched 3' base" stated flat is the banned claim.
# The correct wording states the context: a ranking (Kwok), "some mismatches",
# "less reliably", "context-dependent", "well-extended", or a mechanism (MGB).
is_flat_extension <- function(s) {
  grepl("extend", s, ignore.case = TRUE) &
    grepl("mismatch|3'|3\u2032", s) &
    !grepl(
      paste0(
        "(context-dependent|well-extended|less reliably|readily|ranking|rather|",
        "MGB|minor.groove|suppresses|discriminat|efficienc|some )"
      ),
      s,
      ignore.case = TRUE
    )
}

test_that("no vignette calls a CRISPR readout dsDNase", {
  # the readouts are ssDNase (Cas12a, Chen 2018) / ssRNase (Cas13, Gootenberg 2017)
  expect_equal(grep("dsDNase", vig_lines(), ignore.case = TRUE, value = TRUE), character(0))
})

# Cas12a's target is dsDNA, but its collateral substrate is ssDNA: the activated
# nuclease is an ssDNase that cuts the ssDNA reporter, so a sentence that puts
# Cas12a's collateral cleavage on dsDNA has the substrate wrong.
is_cas12a_ds_collateral <- function(s) {
  grepl("Cas12a", s, ignore.case = TRUE) &
    grepl("collateral|cleav|DNase", s, ignore.case = TRUE) &
    grepl("dsDNA|double[ -]?stranded", s, ignore.case = TRUE)
}

test_that("Cas12a collateral cleavage is not claimed to act on dsDNA", {
  s <- sentences(vig_lines())
  expect_equal(s[is_cas12a_ds_collateral(s)], character(0))
})

test_that("PAM-free is never stated without its PFS qualifier", {
  hit <- grep("PAM-free", vig_lines(), ignore.case = TRUE, value = TRUE)
  # legitimate use: the bare label, named as one the package does not make, in a
  # sentence that also names the PFS
  ok <- grepl("PFS", hit) && grepl("(not|never)", hit)
  expect_equal(hit[!ok], character(0))
})

test_that("no vignette treats a 3' mismatch as non-extendable", {
  expect_equal(grep("non-?extendable", vig_lines(), ignore.case = TRUE, value = TRUE), character(0))
})

test_that("extension past a 3' mismatch is not stated as a fact", {
  s <- sentences(vig_lines())
  expect_equal(s[is_flat_extension(s)], character(0))
})
