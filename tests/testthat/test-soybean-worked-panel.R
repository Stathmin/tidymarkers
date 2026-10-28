# =============================================================================
#  The ordered soybean panel  --  a design benchmark, not an outcome table
# =============================================================================
#  Source: "Напоминалка по праймерам" (2026-09-10), 14 pages, Williams 82
#  Wm82.a4.v1 (GCF_000004515.6).  The document was compiled from the saved
#  selection tables; its header says every number was computed.  It lists 36
#  assays (15 Sanger amplicons, 10 KASP sets, 3 Dt2-flank KASP sets and 8
#  reserve Dt2-body amplicons) and 85 oligos.
#
#  The document prints a per-assay design, not a wet-lab call.  It never says
#  "this band was seen" or "this assay failed", so this file does NOT assert a
#  success rate.  What it does assert is the arithmetic that a real design
#  table must satisfy, and the Tm that the source pipeline printed:
#
#    * 36 assays / 85 oligos, 15 + 10 + 3 + 8 by scheme;
#    * every Sanger product = end - start + 1, and the 10 assays with printed
#      SNP offsets satisfy band = off_F + off_R + len(F) + len(R) - 1;
#    * every KASP product = gap + allele-primer body + common primer;
#    * the two allele-specific bodies differ only at their 3' base, which is
#      the stated allele or its complement (the primer may bind either strand);
#    * the printed Tm is reproduced by primer_tm(..., salt_method =
#      "SantaLucia1998-2") to 0.15 C.  That assertion lives in
#      test-tm-oracle.R, where the cross-model set is defined;
#    * every exact full-length genome count is computed from GCF_000004515.6,
#      not read from the prose: the 30 printed counts are reproduced, every
#      on-target oligo is unique on its own chromosome, and an ALT-allele primer
#      that matches elsewhere shows the off-target and its chromosome.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

assays <- function() ext("soybean_worked_panel_assays.csv")
primers <- function() ext("soybean_worked_panel_primers.csv")
gates <- function() ext("soybean_worked_panel_kasp_gates.csv")

test_that("soybean panel: the order list is 36 assays and 85 oligos", {
  a <- assays()
  p <- primers()
  expect_equal(nrow(a), 36L)
  expect_equal(nrow(p), 85L)
  expect_equal(
    as.integer(table(a$scheme)[c(
      "sanger_amplicon", "kasp", "kasp_dt2_flank", "sanger_dt2_body"
    )]),
    c(15L, 10L, 3L, 8L)
  )
  expect_equal(as.integer(table(a$queue)[c("1", "2", "3")]), c(15L, 13L, 8L))
  expect_false(anyDuplicated(a$assay_id) > 0)
  expect_true(all(a$assay_id %in% p$assay_id))
  expect_equal(
    as.integer(table(p$role)[c("forward", "reverse", "allele1", "allele2", "common")]),
    c(23L, 23L, 13L, 13L, 13L)
  )
  # an assay cannot be present in the joins with no primer, and a primer cannot
  # point at an assay that is not in the panel
  expect_setequal(unique(p$assay_id), a$assay_id)
  expect_false(any(is.na(a$locus) | !nzchar(a$locus)))
  expect_false(any(is.na(a$population) | !nzchar(a$population)))
  expect_false(any(is.na(a$target) | !nzchar(a$target)))
  expect_true(all(a$genome_build == "GCF_000004515.6 (Wm82.a4.v1)"))
  # the document is a design list: no outcome vocabulary may creep in
  expect_false(any(grepl("outcome|success|prob|calibrat|worked|failed", names(a))))
  expect_false(any(grepl("outcome|success|prob|calibrat|worked|failed", names(p))))
})

test_that("soybean panel: every oligo is a pure DNA sequence with the right tail", {
  p <- primers()
  expect_true(all(grepl("^[ACGT]+$", p$sequence)))
  expect_true(all(grepl("^[ACGT]+$", p$body)))
  expect_true(all(nchar(p$body) >= 18L & nchar(p$body) <= 30L))
  expect_equal(sort(unique(p$tail)), c("FAM", "HEX"))
  expect_equal(sum(p$tail == "FAM", na.rm = TRUE), 13L)
  expect_equal(sum(p$tail == "HEX", na.rm = TRUE), 13L)
  expect_equal(sum(is.na(p$tail)), 59L)
  expect_true(all(nchar(p$sequence[p$tail %in% "FAM"]) > nchar(p$body[p$tail %in% "FAM"])))
  expect_true(all(startsWith(
    p$sequence[p$tail %in% "FAM"], "GAAGGTGACCAAGTTCATGCT"
  )))
  expect_true(all(startsWith(
    p$sequence[p$tail %in% "HEX"], "GAAGGTCGGAGTCAACGGATT"
  )))
  # body == sequence minus the published tail, for every row
  body <- ifelse(p$tail %in% "FAM",
    sub("^GAAGGTGACCAAGTTCATGCT", "", p$sequence),
    ifelse(p$tail %in% "HEX",
      sub("^GAAGGTCGGAGTCAACGGATT", "", p$sequence),
      p$sequence
    )
  )
  expect_equal(p$body, body)
  expect_equal(length(unique(p$body)), nrow(p))
})

test_that("soybean panel: the allele-specific 3' base names the allele", {
  p <- primers()
  a <- assays()
  allele <- p[p$role %in% c("allele1", "allele2"), ]
  expect_equal(nrow(allele), 26L)
  expect_true(all(allele$strand %in% c("same", "reverse")))
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  end <- substr(allele$body, nchar(allele$body), nchar(allele$body))
  expect_true(all(end == allele$allele | end == comp[allele$allele]))
  expect_true(all(
    (allele$strand == "same" & end == allele$allele) |
      (allele$strand == "reverse" & end == comp[allele$allele])
  ))
  # the two bodies of one assay share everything except the 3' base
  for (id in unique(allele$assay_id)) {
    two <- allele[allele$assay_id == id, ]
    expect_equal(nrow(two), 2L, info = id)
    expect_equal(
      substr(two$body[[1L]], 1L, nchar(two$body[[1L]]) - 1L),
      substr(two$body[[2L]], 1L, nchar(two$body[[2L]]) - 1L),
      info = id
    )
    expect_false(
      substr(two$body[[1L]], nchar(two$body[[1L]]), nchar(two$body[[1L]])) ==
        substr(two$body[[2L]], nchar(two$body[[2L]]), nchar(two$body[[2L]])),
      info = id
    )
  }
  # and the assay's allele columns are the two bases carried by the primers
  expect_true(all(is.na(a$allele1) | grepl("^[ACGT]+$", a$allele1)))
  expect_true(all(is.na(a$allele2) | grepl("^[ACGT]+$", a$allele2)))
})

test_that("soybean panel: Sanger sizes close on coordinates and SNP offsets", {
  a <- assays()
  p <- primers()
  s <- a[a$scheme %in% c("sanger_amplicon", "sanger_dt2_body"), ]
  expect_equal(nrow(s), 23L)
  expect_true(all(s$product_bp == s$amplicon_end - s$amplicon_start + 1L))
  expect_true(all(s$amplicon_start < s$amplicon_end))

  # The 10 peak assays print the SNP's distance from each primer.  The document's
  # convention is band = off_F + off_R + len(F) + len(R) - 1; a one-base change
  # in either offset must fail this test.
  off <- a[!is.na(a$snp_offset_forward_bp), ]
  expect_equal(nrow(off), 10L)
  fw <- p[p$role == "forward", c("assay_id", "body")]
  rv <- p[p$role == "reverse", c("assay_id", "body")]
  off <- merge(off, fw, by = "assay_id")
  off <- merge(off, rv, by = "assay_id", suffixes = c("_f", "_r"))
  expect_equal(
    off$product_bp,
    off$snp_offset_forward_bp + off$snp_offset_reverse_bp +
      nchar(off$body_f) + nchar(off$body_r) - 1L
  )
  expect_true(all(off$snp_offset_forward_bp > 0))
  expect_true(all(off$snp_offset_reverse_bp > 0))
})

test_that("soybean panel: KASP product = gap + allele body + common", {
  a <- assays()
  p <- primers()
  k <- a[a$scheme %in% c("kasp", "kasp_dt2_flank"), ]
  expect_equal(nrow(k), 13L)
  a1 <- p[p$role == "allele1", c("assay_id", "body")]
  a2 <- p[p$role == "allele2", c("assay_id", "body")]
  com <- p[p$role == "common", c("assay_id", "body")]
  k <- merge(k, a1, by = "assay_id")
  k <- merge(k, a2, by = "assay_id", suffixes = c("_a1", "_a2"))
  k <- merge(k, com, by = "assay_id")
  expect_equal(
    k$product_bp,
    k$gap_bp + nchar(k$body_a1) + nchar(k$body)
  )
  expect_true(all(k$gap_bp > 0))
  expect_true(all(k$product_bp > k$gap_bp))
  expect_equal(nchar(k$body_a1), nchar(k$body_a2))
})

test_that("soybean panel: the GBS class counts keep their units", {
  a <- assays()
  # expected GBS class counts exist for 8 of the 10 KASP 2.1 sets; the branch
  # proxy and the joint model have no per-class count in the document
  expect_equal(sum(!is.na(a$gbs_a1a1)), 8L)
  gbs <- a[!is.na(a$gbs_a1a1), ]
  expect_true(all(gbs$gbs_a1a1 >= 0 & gbs$gbs_het >= 0 & gbs$gbs_a2a2 >= 0))
  expect_true(all(gbs$gbs_a1a1 + gbs$gbs_het + gbs$gbs_a2a2 > 0))
  expect_equal(sum(!is.na(a$het_fraction)), 3L)
  expect_true(all(a$het_fraction[!is.na(a$het_fraction)] >= 0))
  expect_true(all(a$het_fraction[!is.na(a$het_fraction)] <= 1))
  expect_true(all(a$het_n[!is.na(a$het_n)] > 0))
})

test_that("soybean panel: genome hits are computed against Wm82.a4.v1", {
  a <- assays()
  p <- primers()
  # genome_hits is an exact full-length occurrence count, both strands, over the
  # 20 chromosomes of GCF_000004515.6 - not a number copied from the document
  expect_false(any(is.na(p$genome_hits)))
  expect_true(all(p$genome_hits %in% 0:1))
  # every on-target oligo is unique: one exact match, on the assay's chromosome
  on <- merge(p[p$role != "allele2", ], a[, c("assay_id", "target")], by = "assay_id")
  expect_true(all(on$genome_hits == 1L))
  # readr infers the numeric-looking chromosome column as a number
  expect_equal(
    as.character(on$genome_hit_chrom),
    sub("^Gm0*", "", regmatches(on$target, regexpr("Gm[0-9]+", on$target)))
  )
  # an ALT-allele primer matches only where the reference happens to carry its
  # exact sequence; zero is the normal case
  alt <- p[p$role == "allele2", ]
  expect_true(all(alt$genome_hits %in% 0:1))
  expect_true(all(is.na(alt$genome_hit_chrom[alt$genome_hits == 0L])))
  expect_equal(sum(alt$genome_hits == 0L), 12L)
  # the ten KASP 2.1 sets print counts; the genome reproduces all 30 exactly
  printed <- p[!is.na(p$genome_hits_printed), ]
  expect_equal(nrow(printed), 30L)
  expect_equal(printed$genome_hits, printed$genome_hits_printed)
  # The AGP9 height row is resolved by the genome: the A2 body has no exact
  # match, so the structured line (A2 = 0) is right and the prose's "one match"
  # is not an exact one.
  agp9 <- a[a$assay_id == "KASP_AGP9_height_chr19_46.10Mb", ]
  agp9_a2 <- p[
    p$assay_id == "KASP_AGP9_height_chr19_46.10Mb" & p$role == "allele2",
  ]
  expect_equal(agp9_a2$genome_hits, 0L)
  expect_equal(agp9_a2$genome_hits_printed, 0L)
  expect_match(agp9$note, "genome", ignore.case = TRUE)
})

test_that("soybean panel: each exact hit carries a coordinate and a strand", {
  p <- primers()
  has <- p$genome_hits > 0L
  expect_true(all(is.na(p$genome_hit_coord[!has])))
  expect_true(all(is.na(p$genome_hit_strand[!has])))
  coords <- strsplit(as.character(p$genome_hit_coord[has]), ",", fixed = TRUE)
  strands <- strsplit(as.character(p$genome_hit_strand[has]), ",", fixed = TRUE)
  expect_true(all(lengths(coords) == p$genome_hits[has]))
  expect_true(all(lengths(coords) == lengths(strands)))
  expect_true(all(grepl("^[0-9]+:[0-9]+-[0-9]+$", unlist(coords))))
  expect_true(all(unlist(strands) %in% c("+", "-")))
  # the chromosome in the coordinate is the one the row already reports
  expect_equal(
    sub(":.*$", "", unlist(coords)),
    unlist(strsplit(as.character(p$genome_hit_chrom[has]), ",", fixed = TRUE))
  )
})

test_that("soybean panel: each Sanger pair has one exact product on the genome", {
  a <- assays()
  p <- primers()
  s <- a[a$scheme %in% c("sanger_amplicon", "sanger_dt2_body"), ]
  expect_equal(nrow(s), 23L)
  expect_true(all(s$genome_products == 1L))
  chr <- function(x) sub("^Gm0*", "", regmatches(x, regexpr("Gm[0-9]+", x)))
  expect_equal(
    as.character(s$genome_product_coord),
    paste0(chr(s$target), ":", s$amplicon_start, "-", s$amplicon_end)
  )
  # both primers have exactly one exact hit, on the product's chromosome; the
  # forward hit opens the product on the plus strand, the reverse one closes it
  # on the minus strand
  for (i in seq_len(nrow(s))) {
    id <- s$assay_id[[i]]
    f <- p[p$assay_id == id & p$role == "forward", ]
    r <- p[p$assay_id == id & p$role == "reverse", ]
    expect_equal(f$genome_hits, 1L, info = id)
    expect_equal(r$genome_hits, 1L, info = id)
    expect_equal(f$genome_hit_strand, "+", info = id)
    expect_equal(r$genome_hit_strand, "-", info = id)
    f_pos <- sub("^[0-9]+:", "", f$genome_hit_coord)
    r_pos <- sub("^[0-9]+:", "", r$genome_hit_coord)
    expect_equal(sub("-.*$", "", f_pos), as.character(s$amplicon_start[[i]]), info = id)
    expect_equal(sub("^.*-", "", r_pos), as.character(s$amplicon_end[[i]]), info = id)
  }
})

test_that("soybean panel: an ALT primer's off-target is placed on its chromosome", {
  p <- primers()
  # the Dt2 high-flank ALT primer's only exact match is on chr1, not its chr18
  # target; the check did not print a count for this set, but the genome shows it
  off <- p[p$assay_id == "KASP_DT2_flank_high_chr18_56.61Mb" & p$role == "allele2", ]
  expect_equal(off$genome_hits, 1L)
  expect_equal(as.character(off$genome_hit_chrom), "1")
})

test_that("soybean panel: the four redesigned sets are the four named", {
  a <- assays()
  expect_equal(sum(!is.na(a$status)), 4L)
  expect_setequal(
    a$assay_id[!is.na(a$status)],
    c(
      "KASP_AGP9_height_chr19_46.10Mb",
      "KASP_AGPO2_height_chr19_49.74Mb",
      "KASP_AGPO2_height2_chr19_51.07Mb",
      "KASP_JOINT_protein_chr15_11.28Mb"
    )
  )
  expect_true(all(a$status[!is.na(a$status)] == "redesigned after external check"))
  expect_true(all(nzchar(a$note[!is.na(a$status)])))
})

test_that("soybean panel: the 13 KASP sets are graded by the package gates", {
  g <- gates()
  a <- assays()
  k <- a[a$scheme %in% c("kasp", "kasp_dt2_flank"), ]
  n_check <- sum(marker_criteria(scheme = "snp")$role == "check")
  expect_equal(nrow(g), 13L * n_check)
  expect_setequal(unique(g$assay_id), k$assay_id)
  expect_true(all(table(g$assay_id) == n_check))
  expect_true(all(g$status %in% c("pass", "fail", "unknown")))
  # the table carries the package's evaluation columns, not doc prose
  expect_true(all(
    c("stage", "criterion", "observed", "limit", "status", "detail", "level") %in% names(g)
  ))
  # a pass always carries the number it was judged on (textual rows excepted)
  textual <- c("target site class", "AS 3' mismatch base pair")
  num <- g$status == "pass" & !g$criterion %in% textual
  expect_false(any(is.na(g$observed[num])))
})

test_that("soybean panel: the variant/linkage gates are unknown, not passed", {
  g <- gates()
  no_pop <- c(
    "linked variants skipped (n)", "design relies on linkage",
    "AS 3' window variants (n)", "AS footprint variants (n)",
    "COM footprint variants (n)", "AS 3' window dropout",
    "AS footprint dropout", "COM footprint dropout"
  )
  expect_setequal(unique(g$criterion[g$criterion %in% no_pop]), no_pop)
  expect_true(all(g$status[g$criterion %in% no_pop] == "unknown"))
  expect_true(all(is.na(g$observed[g$criterion %in% no_pop])))
})

test_that("soybean panel: the genome-level KASP gates pass", {
  g <- gates()
  off <- g[g$criterion == "off-target exact primer hits (n)", ]
  expect_equal(nrow(off), 13L)
  expect_true(all(off$status == "pass"))
  expect_true(all(off$observed <= 1))
  pcr <- g[g$criterion == "in-silico PCR (REF 1 / ALT 0)", ]
  expect_equal(nrow(pcr), 13L)
  expect_true(all(pcr$status == "pass"))
  # the ALT primer with one exact match on chr1 makes no product with COM
  p <- primers()
  off_target <- p[
    p$assay_id == "KASP_DT2_flank_high_chr18_56.61Mb" & p$role == "allele2",
  ]
  expect_equal(off_target$genome_hits, 1L)
  expect_equal(as.character(off_target$genome_hit_chrom), "1")
})

test_that("soybean panel: the KASP gate failures are the measured ones", {
  g <- gates()
  failed <- function(crit) sort(g$assay_id[g$criterion == crit & g$status == "fail"])
  expect_equal(failed("amplicon length (bp)"), sort(c(
    "KASP_AGP10_height_chr18_58.02Mb",
    "KASP_AGPO2_nodes_chr19_49.65Mb",
    "KASP_AGPO8_protein_chr15_10.93Mb",
    "KASP_DT2_near_chr18_56.14Mb"
  )))
  expect_equal(failed("discrimination"), sort(c(
    "KASP_AGP9_height_chr19_46.10Mb",
    "KASP_AGPO2_height_chr19_49.74Mb",
    "KASP_AGPO2_nodes_chr19_49.65Mb",
    "KASP_AGPO2_seed_chr19_48.03Mb",
    "KASP_AGPO8_protein_chr15_10.93Mb"
  )))
  # Tm passes under the salt model that reproduces the printed Tm; the COM >= AS
  # rule is the one most sets do not meet
  expect_true(all(g$status[g$criterion %in% c("AS Tm (C)", "COM Tm (C)")] == "pass"))
  expect_equal(sum(g$criterion == "COM Tm >= AS Tm" & g$status == "fail"), 11L)
})

test_that("soybean panel: the KASP gate observations are computed", {
  g <- gates()
  a <- assays()
  p <- primers()
  # the amplicon is the genome span, equal to the printed product (asserted by
  # the generator); it is not read from the table
  amp <- g[g$criterion == "amplicon length (bp)", ]
  expect_equal(amp$observed, a$product_bp[match(amp$assay_id, a$assay_id)])
  # the off-target observation is the worst exact-hit count of the three primers
  off <- g[g$criterion == "off-target exact primer hits (n)", ]
  worst <- vapply(
    off$assay_id,
    function(id) max(p$genome_hits[p$assay_id == id]),
    numeric(1L)
  )
  expect_equal(unname(off$observed), unname(worst))
  # discrimination is a DECIPHER efficiency difference in [0, 1]
  disc <- g[g$criterion == "discrimination", ]
  expect_true(all(disc$observed >= 0 & disc$observed <= 1))
})

test_that("soybean panel: the generator reproduces the committed CSVs", {
  skip_if_not_installed("pkgload")
  root <- normalizePath(test_path("../.."), mustWork = FALSE)
  src <- file.path(root, "data-raw", "soybean_worked_panel.R")
  txt <- file.path(root, ".tmp", "lit", "soybean_worked_panel", "paper.txt")
  skip_if_not(file.exists(src), "the generator is not present")
  skip_if_not(file.exists(txt), "the fetched source text is not present")
  # the generator now computes genome_hits from the assembly; CI has no genome
  genome <- file.path(
    tools::R_user_dir("tidymarkers", "cache"), "gm_v4",
    "GCF_000004515.6_Glycine_max_v4.0_genomic.fna"
  )
  skip_if_not(file.exists(paste0(genome, ".fai")), "the Wm82.a4.v1 genome is not cached")
  csv <- file.path(
    root, "inst", "extdata",
    c(
      "soybean_worked_panel_assays.csv",
      "soybean_worked_panel_primers.csv",
      "soybean_worked_panel_kasp_gates.csv"
    )
  )
  skip_if_not(all(file.exists(csv)), "the committed CSVs are not present")
  before <- lapply(csv, readLines, warn = FALSE)
  old <- getwd()
  setwd(root)
  on.exit(setwd(old), add = TRUE)
  # any stop() from an assertion in the generator is the failure, and
  # write_checked() refuses to rewrite a CSV that disagrees with the table
  expect_no_error(
    suppressMessages(
      utils::capture.output(source(src, local = new.env()), type = "output")
    )
  )
  expect_identical(lapply(csv, readLines, warn = FALSE), before)
})
