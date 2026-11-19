# Tetra-primer ARMS-PCR: the outer pair, the inner allele-specific primers, the bands.

test_that("an outer flanking pair sits outside the region and codes a product", {
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  p <- tidymarkers:::design_flanking_pair(w, lo = 300L, hi = 300L)
  expect_true(p$start < 300L && p$end > 300L)
  expect_true(p$product >= 200L && p$product <= 800L)
})

test_that("design_arms_marker gives two resolvable, different-sized diagnostic bands", {
  skip_if_no_oligoarrayaux()
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 800L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  ref <- substr(w, 400L, 400L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 400L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_arms_marker(genome, snp, flank = 250L)
  expect_s3_class(m, "marker_arms")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  b <- m$bands
  expect_setequal(b$band, c("control", "REF", "ALT"))
  expect_gte(abs(b$size[b$band == "REF"] - b$size[b$band == "ALT"]), 20L)
  expect_true(all(b$size > 0))
  expect_true(inherits(plot(m), c("patchwork", "ggplot")))
})

test_that("design_arms_marker declines non-substitution alleles", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 40L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "A", alt = "AA")
  expect_message(suppressWarnings(design_arms_marker(genome, snp, flank = 40L)), "single-base")
  # §W: the refusal is also a DAG criterion, not only a message
  expect_true("family_single_base" %in% marker_criteria(scheme = "arms")$id)
})

test_that("plot_arms_gel draws the AA / AB / BB readout", {
  marker <- structure(
    list(
      snp = tibble::tibble(chrom = "Chr1", pos = 10L, ref = "A", alt = "G"),
      bands = tibble::tibble(
        band = c("control", "REF", "ALT"), size = c(400L, 200L, 250L)
      )
    ),
    class = "marker_arms"
  )
  p <- plot_arms_gel(marker)
  expect_s3_class(p, "ggplot")
  expect_true(all(p$data$col %in% c(GEL_DIAG, GEL_SHARED)))
  expect_false(any(p$data$col %in% unname(KASP_STATUS[c("neutral", "broken")])))
  expect_true(all(p$data$width > 0.2 & p$data$width <= 0.4))
})

test_that("the arms scheme is registered and yields a DAG with its own stage", {
  expect_true("arms" %in% marker_schemes()$scheme)
  crit <- marker_criteria(scheme = "arms")
  expect_true(all(paste0("arms:gate:", KASP_GATES_ARMS) %in% crit$id))
  s <- marker_decision_schema("arms")
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(ord[s$edges$to] > ord[s$edges$from]))
})

test_that("design_arms_marker refuses a signed mismatch position (Z7.1)", {
  withr::local_envvar(LANGUAGE = "en")
  # the position is a distance from the 3' end (3L); the policy constant is the
  # signed offset (-3).  A signed value used to clear the request in silence
  # deep inside design_as_primers(), after the designer had already run.
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 200L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 400L, ref = "A", alt = "T")
  expect_error(
    design_arms_marker(genome, snp, flank = 100L, stab_ats = ARMS_EXTRA_MISMATCH_POSITION),
    "distance from the 3' end"
  )
  expect_error(
    design_arms_marker(genome, snp, flank = 100L, stab_ats = c(-3L, 3L)),
    "distance from the 3' end"
  )
})

test_that("design_arms_marker tolerates a missing stab_ats entry (Z7.1)", {
  skip_if_no_oligoarrayaux()
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 800L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  ref <- substr(w, 400L, 400L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 400L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  # NA is "no extra mismatch", not a negative distance: it is ignored, not an abort
  expect_s3_class(
    design_arms_marker(genome, snp, flank = 250L, stab_ats = NA_integer_),
    "marker_arms"
  )
})

test_that("an inner primer on only one strand returns NULL, with a message", {
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 800L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  ref <- substr(w, 400L, 400L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 400L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  # the designer returns inner primers for the "+" strand only
  inner <- tibble::tibble(
    strand = "+",
    ref_primer = "AAAACCCCGGGGTTTTACGT", alt_primer = "AAAACCCCGGGGTTTTACGA",
    stab_at = NA_integer_, weak_terminal = FALSE,
    tmpl_start = 390L, tmpl_end = 405L, tm = 55, discrimination = 0.9
  )
  outer <- tibble::tibble(
    forward = "ACGTACGTACGTACGTACGT", reverse = "ACGTACGTACGTACGTACGT",
    start = 1L, end = 480L, product = 480L, tm_f = 55, tm_r = 55, dtm = 0
  )
  local_mocked_bindings(
    arms_outer_pairs = function(...) outer,
    design_as_primers = function(...) inner,
    .package = "tidymarkers"
  )
  res <- NULL
  expect_message(
    res <- design_arms_marker(genome, snp, flank = 250L),
    "An inner primer is missing on one strand."
  )
  expect_null(res)
})

test_that("the arms DAG obeys the contract", {
  expect_dag_invariants("arms")
})

test_that("arms declares its gate set, hardness, both orders and real floors", {
  crit <- marker_criteria(scheme = "arms")
  gates <- crit[crit$role == "gate", ]
  expect_setequal(gates$id, paste0("arms:gate:", KASP_GATES_ARMS))
  # application order (grave / cheap first) is the declared vector
  expect_equal(gates$id[order(gates$order)], paste0("arms:gate:", KASP_GATES_ARMS))
  # only "no allele-specific inner primer" is hard (C8)
  expect_equal(gates$id[gates$hard], "arms:gate:gate_arms_disc")
  expect_setequal(gates$id[!gates$hard], paste0("arms:gate:", KASP_RELAX_ARMS))
  # the ladder is least grave first, and the band-gap floor is *below* its bar -
  # the audit's "floor equals threshold, drawn soft, can never relax" bug
  expect_equal(KASP_RELAX_ARMS, c("gate_arms_outer", "gate_arms_bands"))
  expect_lt(KASP_GATE_FLOOR_ARMS[["gate_arms_bands"]], 0)
  expect_lt(KASP_ARMS_GAP_FLOOR, KASP_ARMS_WINDOW$band_gap)
  # limits are real values, never a `[lo, hi]` placeholder (C9)
  expect_false(any(grepl("\\[lo, ?hi\\]", gates$limit)))
})

test_that("design_arms_marker evaluates the arms gates and gives a resolvable gel", {
  skip_if_no_oligoarrayaux()
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 800L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  ref <- substr(w, 400L, 400L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 400L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_arms_marker(genome, snp, flank = 250L)
  # the gates are *evaluated*, not decorative (C7): the marker carries the outcome
  expect_setequal(names(m$selection), c("relaxed", "violated"))
  expect_true(all(unlist(m$selection) %in% KASP_GATES_ARMS))
  # the control band is the largest product (the inner primers are nested), and
  # both diagnostic bands clear the accepted separation from each other and it
  b <- m$bands
  ctrl <- b$size[b$band == "control"]
  ref_size <- b$size[b$band == "REF"]
  alt_size <- b$size[b$band == "ALT"]
  expect_gt(ctrl, max(ref_size, alt_size))
  # the band sizes are inclusive spans of the template coordinates (Z7)
  expect_equal(ref_size, m$outer$end - m$inner$ref$start + 1L)
  expect_equal(alt_size, m$inner$alt$end - m$outer$start + 1L)
  expect_gte(min(abs(ref_size - alt_size), ctrl - max(ref_size, alt_size)), KASP_ARMS_WINDOW$band_gap)
  expect_gte(m$outer$product, KASP_ARMS_WINDOW$control[1L])
  expect_lte(m$outer$product, KASP_ARMS_WINDOW$control[2L])
})

test_that("a below-floor design reaches the second pass instead of being doomed", {
  skip_if_no_oligoarrayaux()
  withr::with_seed(1L, w <- paste(sample(c("A", "C", "G", "T"), 800L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  ref <- substr(w, 400L, 400L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 400L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  # every inner primer is below the hard discrimination floor, so no candidate
  # clears the gates: the ladder hands the least-bad one to the second pass, which
  # accepts it (the scheme is tolerant), instead of dooming the target
  low <- tibble::tibble(
    strand = c("+", "-"),
    ref_primer = c("AAAACCCCGGGGTTTTACGT", "TTTTGGGGCCCCAATGCA"),
    alt_primer = c("AAAACCCCGGGGTTTTACGA", "TTTTGGGGCCCCAATGCG"),
    stab_at = NA_integer_, weak_terminal = FALSE,
    tmpl_start = c(390L, 410L), tmpl_end = c(405L, 425L), tm = 55,
    discrimination = 0.1
  )
  local_mocked_bindings(design_as_primers = function(...) low, .package = "tidymarkers")
  w <- capture_warnings(m <- design_arms_marker(genome, snp, flank = 250L))
  expect_length(w, 2L)
  expect_match(w[[1L]], "No tetra-primer set's outer pair gives a single in-silico PCR product", fixed = TRUE)
  expect_match(w[[2L]], "least-bad", fixed = TRUE)
  expect_s3_class(m, "marker_arms")
  expect_true("gate_arms_disc" %in% m$selection$violated)
})
