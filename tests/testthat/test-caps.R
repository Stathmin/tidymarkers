# CAPS / dCAPS: the enzyme scan and the flanking-pair design.

test_that("caps_fragments cuts at the enzyme site, or leaves one uncut band", {
  f <- caps_fragments("AAAAAGAATTCTTTTT") # one EcoRI site at 6 (cut after it)
  e <- f[f$enzyme == "EcoRI", ]
  expect_equal(e$n_sites, 1L)
  expect_equal(e$fragments[[1L]], c(6L, 10L))
  # an enzyme with no site in a short sequence leaves the molecule uncut
  expect_equal(caps_fragments("AAAAAAAAAA")$fragments[[1L]], 10L)
})

test_that("find_caps_enzymes flags the enzyme the variant turns on or off", {
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 400L, TRUE), collapse = ""))
  substr(w, 200L, 205L) <- "GAATTC" # a unique EcoRI site at 200
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  snp <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = "G", alt = "A") # G>A destroys it
  hits <- find_caps_enzymes(genome, snp, flank = 150L)
  eco <- hits[hits$enzyme == "EcoRI", ]
  expect_equal(nrow(eco), 1L)
  expect_true(eco$resolvable)
  expect_false(identical(eco$fragments_ref[[1L]], eco$fragments_alt[[1L]]))
})

test_that("design_caps_marker returns a classed marker whose pair spans the site", {
  withr::with_seed(2L, w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  substr(w, 300L, 305L) <- "GAATTC"
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "G", alt = "A")
  m <- design_caps_marker(genome, snp, flank = 250L)
  expect_s3_class(m, "kasp_caps")
  expect_true(m$amplicon$start < 300L && m$amplicon$end > 305L)
  expect_true(m$primers$product >= 120L)
})

test_that("the caps scheme is registered and yields a DAG with its own stage", {
  expect_true("caps" %in% kasp_schemes()$scheme)
  crit <- kasp_criteria(scheme = "caps")
  expect_true(all(paste0("caps:gate:", KASP_GATES_CAPS) %in% crit$id))
  # the marker-agnostic primer checks apply, the KASP 3'-terminal ones do not
  expect_true("AS GC" %in% crit$label)
  expect_false("in-silico PCR (REF 1 / ALT 0)" %in% crit$label)
  s <- kasp_decision_schema("caps")
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(ord[s$edges$to] > ord[s$edges$from]))
})
