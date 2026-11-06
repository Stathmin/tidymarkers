# ASQ / PACE: a real engine (Kalendar et al. 2022).
#
# ASQ reads the same competitive allele-specific extension as KASP but places the
# variant at the penultimate 3' base of the allele-specific primer and reports
# through universal probes against a single common quencher.  The reporter /
# plexity rules are panel/format properties and stay out of the gates (C7/C8);
# the geometry and ASQ's windows (100-200 bp, 62-64 C) are the engine's job.

test_that("asq is a real scheme: a designer, a ladder slot, a registered row", {
  expect_dag_invariants("asq")
  expect_true(exists("design_asq_marker", inherits = TRUE))
  expect_true("asq" %in% per_locus_schemes())
  expect_true("asq" %in% marker_schemes()$scheme)
})

test_that("the asq stage declares no gate of its own, because the shared stages compute it", {
  # the AS / pair stages already gate ASQ under asq_limits(); an extra
  # placement "gate" would be a gate no code path computes (C7)
  expect_length(KASP_GATES_ASQ, 0L)
  s <- marker_decision_schema("asq")
  expect_false(any(grepl("^asq:gate:", s$nodes$name)))
  expect_false(any(grepl("gate_asq_cassette|gate_asq_readout|gate_asq_asp_geometry", s$nodes$name)))
  txt <- marker_decision_text("asq")
  expect_true(any(grepl("ASQ primer geometry", txt, fixed = TRUE)))
  expect_true(any(grepl(marker_scheme("asq")$discrimination, txt, fixed = TRUE)))
})

test_that("the gate-less asq stage emits only its ranking tail, never a DOOMED edge (C13)", {
  s <- marker_decision_schema("asq")
  expect_false("asq:decision" %in% s$nodes$name)
  expect_false(any(grepl("^asq:", s$edges$from[s$edges$kind == "doomed"])))
  expect_true("asq:survivors" %in% s$nodes$name)
  expect_true(all(paste0("asq:tier:", names(KASP_TIERS_ASQ)) %in% s$nodes$name))
  ed <- paste(s$edges$from, s$edges$to, sep = " -> ")
  expect_true("pair:best -> asq:survivors" %in% ed)
  expect_length(KASP_RELAX_ASQ, 0L)
  expect_false(any(s$nodes$kind == "relax" & s$nodes$stage == "asq"))
})

test_that("asq keeps the KASP family's checks, because its core is KASP's (C2)", {
  chk <- marker_criteria(scheme = "asq")
  chk <- chk$label[chk$role == "check"]
  expect_true(all(KASP_CHECK_PRIMER %in% chk))
  expect_true(all(KASP_CHECK_AS %in% chk))
  expect_true(all(c("AS Tm (C)", "COM Tm (C)", "AS GC", "COM GC") %in% chk))
  expect_true(all(KASP_CHECK_SITE %in% chk))
})

test_that("design_asq_marker places the variant at the penultimate 3' base", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  a <- suppressWarnings(design_asq_marker(f$genome, f$snp, cores = 1L))
  expect_s3_class(a, "marker_assay")
  expect_identical(a$as_variant_offset, 1L)

  aa <- a$assay[1L, ]
  L <- nchar(aa$as_ref)
  k <- L - 1L
  # the two allele primers differ at the penultimate base and agree at the
  # terminus: that is the ASQ geometry (KASP would differ at the terminus)
  expect_false(substr(aa$as_ref, k, k) == substr(aa$as_alt, k, k))
  expect_identical(substr(aa$as_ref, L, L), substr(aa$as_alt, L, L))

  # the penultimate base is the site's REF / ALT allele on the primer's strand
  s_ref <- a$snp$ref
  s_alt <- a$snp$alt
  if (identical(aa$strand, "-")) {
    s_ref <- as.character(Biostrings::complement(Biostrings::DNAString(s_ref)))
    s_alt <- as.character(Biostrings::complement(Biostrings::DNAString(s_alt)))
  }
  expect_identical(substr(aa$as_ref, k, k), s_ref)
  expect_identical(substr(aa$as_alt, k, k), s_alt)

  # ASQ's own windows, carried on the assay and used by the acceptance checks
  expect_identical(a$limits$amplicon, ASQ_AMPLICON)
  expect_identical(a$limits$as_tm, ASQ_AS_TM_LIMITS)
  expect_true(dplyr::between(aa$amplicon, ASQ_AMPLICON[1L], ASQ_AMPLICON[2L]))
  chk <- a$evaluation
  amp <- chk[chk$criterion == "amplicon length (bp)", ]
  expect_true(amp$pass)
  # the ASQ AS-Tm rule is graded on the designed REF primer, not the REF/ALT
  # mean (whose ALT member carries the SNP mismatch)
  expect_true(chk$pass[chk$criterion == "AS Tm (C)"])
  # without a variant table some checks are NA, so the verdict is `unknown` -
  # never `fail`
  expect_false(identical(evaluation_verdict(a), "fail"))
})

test_that("asq takes the rescue ladder and the batch designer", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  loci <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = f$ref, alt = f$alt)
  out <- suppressWarnings(design_markers(f$genome, loci, scheme = "asq"))
  expect_true(out$ok)
})

test_that("marker_criteria shows ASQ's own windows, not KASP's", {
  chk <- marker_criteria(scheme = "asq")
  expect_true(any(grepl("100", chk$limit[chk$label == "amplicon length (bp)"], fixed = TRUE)))
  expect_true(any(grepl("62", chk$limit[chk$label == "AS Tm (C)"], fixed = TRUE)))
})

test_that("an indel is doomed rather than silently designed as KASP geometry", {
  f <- snp_fixture()
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = f$ref, alt = paste0(f$ref, "AAAA"))
  d <- suppressWarnings(design_asq_marker(f$genome, snp))
  expect_s3_class(d, "marker_doomed")
})
