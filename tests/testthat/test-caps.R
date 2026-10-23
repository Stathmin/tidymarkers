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
  expect_s3_class(m, "marker_caps")
  # the pair spans the discriminating site, i.e. straddles the variant
  expect_true(m$amplicon$start < 300L && m$amplicon$end > 300L)
  expect_true(m$primers$product >= 120L)
  # the object carries the actual recognition site and its enzyme
  enz <- restriction_enzymes()
  sel <- enz[enz$enzyme == m$enzyme, ]
  expect_equal(m$site$enzyme, m$enzyme)
  expect_equal(m$site$motif, sel$site)
  expect_equal(m$site$cut, sel$cut)
  expect_gte(m$site$start, m$amplicon$start)
  expect_lte(m$site$end, m$amplicon$end)
  expect_equal(m$site$end - m$site$start + 1L, nchar(m$site$motif))
  # the map shows a real block, not a zero-width point, and names the enzyme
  spec <- marker_map_spec(m, "caps")
  site <- spec[spec$kind == "site", ]
  expect_equal(nrow(site), 1L)
  expect_equal(site$start, m$site$start)
  expect_equal(site$end, m$site$end)
  expect_match(site$label, m$enzyme, fixed = TRUE)
  expect_match(site$label, "\\^")
  expect_true(inherits(plot(m), c("patchwork", "ggplot")))
})

test_that("find_dcaps_enzymes finds a site a primer can complete", {
  withr::with_seed(3L, w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  substr(w, 300L, 305L) <- "GAGTTC" # one mismatch from EcoRI (G^AATTC) at offset 2
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  snp <- tibble::tibble(chrom = "Chr1", pos = 305L, ref = "C", alt = "A") # variant inside the site
  d <- find_dcaps_enzymes(genome, snp, flank = 200L)
  eco <- d[d$enzyme == "EcoRI", ]
  expect_true(nrow(eco) >= 1L)
  expect_equal(eco$allele[1L], "ref")
  expect_equal(eco$mm_offset[1L], 2L)
  expect_equal(eco$fix_base[1L], "A")
  expect_true(eco$resolvable[1L])
})

test_that("dCAPS scanning declines non-substitution alleles", {
  withr::local_envvar(LANGUAGE = "en")
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 20L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 20L, ref = "A", alt = "AA")
  expect_message(
    suppressWarnings(find_dcaps_enzymes(genome, snp, flank = 15L)),
    "single-base"
  )
})

test_that("design_caps_marker falls back to dCAPS when no natural site exists", {
  withr::with_seed(3L, w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  substr(w, 300L, 305L) <- "GAGTTC"
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  snp <- tibble::tibble(chrom = "Chr1", pos = 305L, ref = "C", alt = "A")
  # the fixture only admits a relaxed dCAPS call; the relaxation is the point
  m <- suppressWarnings(design_caps_marker(genome, snp, flank = 200L))
  expect_s3_class(m, "marker_caps")
  expect_true(isTRUE(m$dcaps))
  expect_true(m$enzyme %in% restriction_enzymes()$enzyme)
})

test_that("a resolvable CAPS call needs a site change, not just a length shift", {
  withr::with_seed(5L, w <- paste(sample(c("A", "C", "G", "T"), 300L, TRUE), collapse = ""))
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  ins <- tibble::tibble(chrom = "Chr1", pos = 150L, ref = "A", alt = "AA") # 1-bp insertion
  h <- find_caps_enzymes(genome, ins, flank = 100L)
  if (nrow(h)) expect_false(any(h$resolvable & !h$site_change))
})

test_that("plot_caps_gel draws a monochrome, Poisson-width two-lane readout", {
  withr::with_seed(2L, w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  substr(w, 300L, 305L) <- "GAATTC"
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "G", alt = "A")
  m <- design_caps_marker(genome, snp, flank = 250L)
  p <- plot_caps_gel(m)
  expect_s3_class(p, "ggplot")
  # allele identity is carried by lane + label, not by reserved blue/red
  expect_true(all(p$data$col %in% c("#44403C", "#D6D3D1")))
  expect_false(any(p$data$col %in% unname(KASP_STATUS[c("neutral", "broken")])))
  expect_true(all(p$data$width > 0.2 & p$data$width <= 0.4))
  expect_match(p$labels$title, m$enzyme, fixed = TRUE)
  expect_match(p$labels$title, "\\^")
})

test_that("the caps scheme is registered and yields a DAG with its own stage", {
  expect_true("caps" %in% marker_schemes()$scheme)
  crit <- marker_criteria(scheme = "caps")
  expect_true(all(paste0("caps:gate:", KASP_GATES_CAPS) %in% crit$id))
  # a flanking pair has no AS/COM oligos: the marker-agnostic primer checks apply,
  # the AS/COM-labelled ones do not (C2 regression)
  expect_false("AS GC" %in% crit$label)
  expect_false("AS Tm (C)" %in% crit$label)
  expect_false("COM GC" %in% crit$label)
  expect_false("COM Tm (C)" %in% crit$label)
  expect_false("in-silico PCR (REF 1 / ALT 0)" %in% crit$label)
  s <- marker_decision_schema("caps")
  ord <- stats::setNames(s$nodes$order, s$nodes$name)
  expect_true(all(ord[s$edges$to] > ord[s$edges$from]))
})

test_that("the caps DAG obeys the contract", {
  expect_dag_invariants("caps")
})

test_that("caps declares its gate set, hardness, both orders and real floors", {
  crit <- marker_criteria(scheme = "caps")
  gates <- crit[crit$role == "gate", ]
  expect_setequal(gates$id, paste0("caps:gate:", KASP_GATES_CAPS))
  # application order (grave / cheap first) is the declared vector
  expect_equal(gates$id[order(gates$order)], paste0("caps:gate:", KASP_GATES_CAPS))
  # only "no discriminating enzyme" is hard (C8); the rest form the ladder
  expect_equal(gates$id[gates$hard], "caps:gate:gate_caps_site")
  expect_setequal(gates$id[!gates$hard], paste0("caps:gate:", KASP_RELAX_CAPS))
  expect_true(all(KASP_RELAX_CAPS %in% KASP_GATES_CAPS))
  # the band gap is relaxable but *bounded*: a relaxed winner is never an
  # unresolvable digest (the audit's "relaxing gate_caps_resolvable" risk)
  expect_lt(KASP_GATE_FLOOR_CAPS[["gate_caps_resolvable"]], 0)
  expect_gt(KASP_CAPS_GAP_FLOOR, 0)
  expect_lt(KASP_CAPS_GAP_FLOOR, KASP_CAPS_GAP_MIN)
  # limits are real values, never a `[lo, hi]` placeholder (C9)
  expect_false(any(grepl("\\[lo, ?hi\\]", gates$limit)))
  expect_true(any(grepl("20 bp", gates$limit)))
  expect_true(any(grepl("120, 1200", gates$limit)))
})

test_that("design_caps_marker evaluates the caps gates and reports the product gel", {
  withr::with_seed(2L, w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = ""))
  substr(w, 300L, 305L) <- "GAATTC"
  genome <- Biostrings::DNAStringSet(c(Chr1 = w))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "G", alt = "A")
  m <- design_caps_marker(genome, snp, flank = 250L)
  # the gates are *evaluated*, not decorative (C7): the marker carries the outcome
  expect_setequal(names(m$selection), c("relaxed", "violated"))
  expect_true(all(unlist(m$selection) %in% KASP_GATES_CAPS))
  expect_false("gate_caps_resolvable" %in% m$selection$violated)
  # the reported digest is the *product*'s (the gel the bench would show): both
  # allele fragments add up to the amplicon, and the closest differing bands clear
  # the accepted 20 bp bar
  expect_equal(sum(m$fragments_ref), m$primers$product)
  expect_equal(sum(m$fragments_alt), m$primers$product)
  expect_gte(m$gap, KASP_CAPS_GAP_MIN)
})

## --- the CAPS extracts -------------------------------------------------------
# The enzyme catalogue and the three paper tables behind the CAPS family are read
# from inst/extdata, so what the scanner is handed is asserted, not assumed.

caps_ext <- function(name) {
  readr::read_csv(testthat::test_path("..", "..", "inst", "extdata", name),
    show_col_types = FALSE, na = c("", "NA")
  )
}

test_that("restriction_enzymes: one row per palindromic site, cut inside the site", {
  e <- caps_ext("restriction_enzymes.csv")
  expect_equal(nrow(e), 1236L)
  expect_false(any(duplicated(e$enzyme)))
  # an unnamed enzyme cannot be reported back to the user
  expect_false(any(is.na(e$enzyme) | !nzchar(e$enzyme)))
  # the scanner reads the plus strand alone, so an asymmetrically cutting site
  # would be modelled wrongly: every site must be its own reverse complement
  iupac <- c(
    A = "T", C = "G", G = "C", T = "A", R = "Y", Y = "R", S = "S", W = "W",
    K = "M", M = "K", B = "V", V = "B", D = "H", H = "D", N = "N"
  )
  rc <- vapply(strsplit(e$site, ""), function(b) {
    paste(rev(unname(iupac[match(b, names(iupac))])), collapse = "")
  }, "")
  expect_false(anyNA(rc)) # only IUPAC letters, no stray characters
  expect_true(all(rc == e$site))
  # a cut outside 0..nchar(site) cannot produce two fragments
  expect_true(all(nchar(e$site) >= 4L & nchar(e$site) <= 12L))
  expect_true(all(e$cut >= 0 & e$cut <= nchar(e$site)))
  expect_lte(max(e$cut), 8)
  # the five named in the documentation still carry their published motif/cut
  canon <- c(EcoRI = "GAATTC", HindIII = "AAGCTT", BamHI = "GGATCC", PstI = "CTGCAG", Sau3AI = "GATC")
  got <- e[match(names(canon), e$enzyme), ]
  expect_equal(got$site, unname(canon))
  expect_equal(got$cut, c(1, 1, 1, 5, 0))
})

test_that("liu2016_caps_snp_types: the pct is within family, and the families sum to 100", {
  d <- caps_ext("liu2016_caps_snp_types.csv")
  # the ten rows are three DENOMINATORS, not one: a percent of 244961 would make
  # every class a fraction of a percent, so the denominator is the family
  expect_equal(nrow(d), 10L)
  expect_setequal(unique(d$family), c("Transversions", "Transitions", "Indels"))
  # the table prints its shares with a % sign, so the column is text on disk
  expect_true(is.character(d$pct))
  d$pct_num <- as.numeric(sub("%$", "", d$pct))
  expect_false(anyNA(d$pct_num))
  fam <- split(d, d$family)
  expect_equal(sum(fam$Transversions$n), 72285L)
  expect_equal(sum(fam$Transitions$n), 153408L)
  expect_equal(sum(fam$Indels$n), 19268L)
  # a within-family share must add to 100 (the published rows are rounded)
  for (f in names(fam)) {
    expect_equal(sum(fam[[f]]$pct_num), 100, tolerance = 0.6, label = f)
  }
  # transversions and transversions-alt are four classes each, transitions two
  expect_equal(nrow(fam$Transversions), 4L)
  expect_equal(nrow(fam$Transitions), 2L)
  expect_equal(nrow(fam$Indels), 4L)
  # the total the shares are of
  expect_equal(sum(d$n), 244961L)
})

test_that("snap2caps_accession_groups: the predicted pattern is the one observed", {
  g <- caps_ext("snap2caps_accession_groups.csv")
  expect_equal(nrow(g), 29L)
  # fragment_sum_bp is the digest total, so it IS the PCR product
  expect_equal(g$fragment_sum_bp, g$pcr_product_bp)
  # every observed band list also adds up to the product (nothing was lost)
  obs_n <- lengths(strsplit(g$observed_fragments_bp, ","))
  obs_bp <- vapply(strsplit(g$observed_fragments_bp, ","), function(x) {
    sum(as.integer(trimws(x)))
  }, integer(1))
  expect_equal(obs_bp, g$pcr_product_bp)
  # a cut is predicted exactly when a site is predicted; the two are one fact
  expect_equal(g$predicted_cut, !is.na(g$predicted_sites_bp))
  # a predicted cut leaves more than one band; an uncut accession leaves one
  expect_true(all(obs_n[g$predicted_cut] >= 2L))
  expect_true(all(obs_n[!g$predicted_cut] == 1L))
  expect_equal(obs_bp[!g$predicted_cut], g$pcr_product_bp[!g$predicted_cut])
  # 29 marker x enzyme rows over 13 markers, and every one confirmed on the gel
  expect_equal(length(unique(g$marker)), 13L)
  expect_true(all(g$outcome == "predicted pattern confirmed on gel"))
})

test_that("soybean_est_caps_markers: one primer pair per marker, pure ACGT", {
  s <- caps_ext("soybean_est_caps_markers.csv")
  expect_equal(nrow(s), 26L)
  expect_false(any(duplicated(s$primer_id)))
  # a CAPS primer pair is a literal sequence: an ambiguity code has no Tm
  expect_true(all(grepl("^[ACGT]+$", s$forward)))
  expect_true(all(grepl("^[ACGT]+$", s$reverse)))
  expect_true(all(nzchar(s$primer_id)))
})
