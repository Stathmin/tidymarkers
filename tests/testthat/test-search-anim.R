# Design-search animation: the frame table must walk the search and say why.

# A tiny real assay on a random genome (skips without OligoArrayAux).  Random
# genomes occasionally yield no GC/Tm-passing window, so try a few seeds.
make_search_assay <- function(seed = 4L, n = 500L, pos = 250L, variants = NULL) {
  for (s in seed:(seed + 19L)) {
    W <- withr::with_seed(s, paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = ""))
    genome <- Biostrings::DNAStringSet(c(Chr1 = W))
    ref <- as.character(Biostrings::subseq(genome[["Chr1"]], pos, pos))
    snp <- tibble::tibble(
      chrom = "Chr1", pos = pos, ref = ref,
      alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
    )
    vs <- variants %||% tibble::tibble(chrom = "Chr1", pos = pos, consequence = "missense_variant")
    assay <- try(suppressWarnings(
      assemble_kasp_assay(genome, snp, variants = vs, flank = 150L)
    ), silent = TRUE)
    if (!inherits(assay, "try-error") && !inherits(assay, "kasp_doomed")) {
      return(list(assay = assay, snp = snp, variants = vs))
    }
  }
  testthat::skip("no random genome assembled into an assay")
}

test_that("design_search_frames walks the AS then COM search and ends accepted", {
  skip_if_no_oligoarrayaux()
  m <- make_search_assay()
  fr <- design_search_frames(m$assay, m$variants, m$snp)

  expect_true(all(c("primers", "notes", "track", "snp", "win") %in% names(fr)))
  expect_equal(fr$notes$step, seq_len(nrow(fr$notes))) # contiguous
  expect_equal(fr$notes$phase[1L], "AS-primer search")
  expect_equal(fr$notes$phase[nrow(fr$notes)], "Accepted assay")
  expect_true(fr$notes$key[1L] && fr$notes$key[nrow(fr$notes)]) # hold the winners

  # step 1 IS the chosen AS placement
  expect_equal(
    dplyr::filter(fr$primers, .data$step == 1L, .data$role == "AS-ref")$start,
    m$assay$assay$as_start[1L]
  )

  # every search step is tagged with a verdict (chosen or a reason)
  search <- fr$notes[grepl("search", fr$notes$phase), ]
  expect_true(all(grepl("CHOSEN|rejected", search$detail)))
})

test_that("the animation's winner is the algorithm's winner (minimal 3' window)", {
  skip_if_no_oligoarrayaux()
  m <- make_search_assay(seed = 7L, n = 600L, pos = 300L)
  fr <- design_search_frames(m$assay, m$variants, m$snp)
  w <- m$assay$as_candidates
  # the first frame is the placement with the cleanest 3' window
  expect_equal(m$assay$as_candidates$n_3p[1L], min(w$n_3p))
  expect_match(fr$notes$detail[1L], "CHOSEN")
  expect_match(fr$notes$headline[nrow(fr$notes)], "^accepted")
})

test_that("plot_design_search_step draws a panel for the first and last step", {
  skip_if_no_oligoarrayaux()
  m <- make_search_assay(seed = 11L, n = 400L, pos = 200L)
  fr <- design_search_frames(m$assay, m$variants, m$snp)
  expect_s3_class(plot_design_search_step(fr, 1L), "ggplot")
  expect_s3_class(plot_design_search_step(fr, nrow(fr$notes)), "ggplot")
})
