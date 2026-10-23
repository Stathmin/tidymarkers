# tidymarkers against Primer3 2.6.1 -- the shared layer, not the marker.
#
# WHAT IS HERE
#   * Thermodynamics on the SAME oligos with the SAME salt and temperature:
#     `primer_tm()` vs `oligotm` (Primer3's SantaLucia 1998 NN table and its
#     default SantaLucia 1998 salt correction), and `hairpin_dg_thermo()` /
#     `duplex_dg_thermo()` (OligoArrayAux) vs `ntthal`.  The assertions are
#     bands on the measured offset and on the rank correlation, never "this
#     oligo melts at X C": an absolute number here would be a property of this
#     machine, not of the package.
#   * A parser self-check on `ntthal`: dG must equal dH - T dS, so the numbers
#     compared above cannot be a column read from the wrong field.
#   * Under TIDYMARKERS_FULL_TESTS=1 only: Primer3 and `design_kasp_markers()`
#     on three real Arabidopsis SNPs, with Primer3 constrained to the package's
#     own amplicon window (50-100 bp) so a length disagreement cannot pose as a
#     quality disagreement.
#
# WHAT IS NOT HERE
#   * No claim that Primer3 designs a KASP assay.  It does not design
#     allele-specific primers; the design test asks only whether the generic
#     "any primer pair" layer agrees on the same loci.
#   * No absolute Tm or dG constants, and no "agreement means the duplex melts
#     there": engine-to-engine agreement is not evidence about a real duplex.
#   * Nothing fetched from the network; every input is local.
#
# The numbers behind the bands: `Rscript dev/primer3-crosscheck.R` (writes
# results/primer3-thermo.csv and results/primer3-design.csv, both git-ignored).

pc_have <- function(bins) {
  vapply(bins, function(b) nzchar(Sys.which(b)), logical(1))
}

pc_skip_if_missing <- function(bins) {
  missing <- bins[!pc_have(bins)]
  if (length(missing)) {
    testthat::skip(paste0(
      "Primer3 binaries not on PATH: ", paste(missing, collapse = ", ")
    ))
  }
}

# 24 deterministic 20-mers (the seed of test-tm-oracle.R) plus the edge probes
# that helper defines, so the band is exercised at both GC extremes.
pc_set <- function() {
  set.seed(20040301L)
  r <- vapply(seq_len(24L), function(i) {
    paste(sample(c("A", "C", "G", "T"), 20L, replace = TRUE), collapse = "")
  }, character(1))
  s <- c(r, tm_edge_oligos())
  rlang::set_names(unname(s), c(paste0("oligo", seq_along(r)), names(tm_edge_oligos())))
}

# `oligotm` at the package's conditions.  -dm is omitted on purpose: oligotm
# 2.6.1 returns -273.15 when it is given an explicit `-dm 0`, while its own
# default (also 0 % DMSO) is the documented behaviour.
pc_oligotm <- function(seqs, cond = reaction_conditions()) {
  vapply(unname(seqs), function(s) {
    args <- c(
      "-mv", as.character(cond$Na), "-dv", as.character(cond$Mg),
      "-n", as.character(cond$dNTPs), "-d", "50"
    )
    out <- suppressWarnings(system2("oligotm", c(args, s), stdout = TRUE, stderr = TRUE))
    as.numeric(trimws(out[[1L]]))
  }, numeric(1))
}

# ntthal's answering lines: dG/dH/dS in cal/mol, t in C.
pc_ntthal <- function(a, s1, s2 = NULL, cond = reaction_conditions()) {
  args <- c(
    "-a", a, "-mv", as.character(cond$Na), "-dv", as.character(cond$Mg),
    "-n", as.character(cond$dNTPs), "-d", "50", "-t", as.character(cond$temperature),
    "-s1", s1, if (!is.null(s2)) c("-s2", s2)
  )
  paste(suppressWarnings(system2("ntthal", args, stdout = TRUE, stderr = TRUE)), collapse = "\n")
}

pc_field <- function(txt, key) {
  m <- regmatches(txt, regexpr(paste0(key, " = [^ \t\n]+"), txt))
  if (!length(m)) NA_real_ else as.numeric(sub(paste0(key, " = "), "", m[[1L]]))
}

# "No secondary structure could be calculated" means dG 0 (nothing folds), the
# same convention the package now applies to UNAFold's ~+1000 sentinel.
pc_dg <- function(txt) {
  if (grepl("No secondary structure", txt)) {
    return(0)
  }
  pc_field(txt, "dG") / 1000
}

# --- thermodynamics on the same oligos ---------------------------------------

test_that("oligotm and primer_tm agree at matched salt, in rank and in level", {
  pc_skip_if_missing("oligotm")
  seqs <- pc_set()
  pkg <- unname(primer_tm(unname(seqs), conditions = reaction_conditions()))
  p3 <- pc_oligotm(seqs)
  keep <- is.finite(pkg) & is.finite(p3) & pkg > 20 & pkg < 90
  expect_gte(sum(keep), 24L)
  d <- p3[keep] - pkg[keep]

  # the two engines must rank the panel alike before any level is discussed
  expect_gt(stats::cor(pkg[keep], p3[keep]), 0.98)
  # bands, not constants: measured mean ~+3.9 C, sd ~0.25 C, max ~4.7 C
  expect_lt(abs(mean(d)), 10)
  expect_lt(stats::sd(d), 2.5)
  expect_lt(max(abs(d)), 15)
})

test_that("ntthal and the OligoArrayAux wrapper agree on hairpins and dimers", {
  pc_skip_if_missing("ntthal")
  skip_if_no_oligoarrayaux()
  seqs <- pc_set()

  # parser self-check: a perfectly matched duplex must satisfy dG = dH - T dS.
  # ntthal prints each field rounded, so the identity holds to ~1 cal/mol, not
  # to the last digit.
  txt <- pc_ntthal("ANY", "AAAAAAAAAA", "TTTTTTTTTT")
  t_k <- 273.15 + reaction_conditions()$temperature
  identity_err <- abs(pc_field(txt, "dG") - (pc_field(txt, "dH") - t_k * pc_field(txt, "dS")))
  expect_lt(identity_err, 1)

  p3 <- vapply(unname(seqs), function(s) pc_dg(pc_ntthal("HAIRPIN", s)), numeric(1))
  pkg <- unname(hairpin_dg_thermo(unname(seqs), conditions = reaction_conditions()))
  keep <- is.finite(p3) & is.finite(pkg)
  expect_gte(sum(keep), 24L)
  d <- p3[keep] - pkg[keep]

  # no hairpin folds in either engine on a homopolymer (the sentinel case)
  expect_identical(hairpin_dg_thermo("AAAAAAAAAAAAAAAAAAAA"), 0)
  expect_gt(stats::cor(pkg[keep], p3[keep]), 0.5)
  # measured: mean -0.84, sd 0.76 kcal/mol -- one side's "no structure" is 0
  expect_lt(abs(mean(d)), 3)
  expect_lt(stats::sd(d), 2.5)

  a <- "ACGTACGTACGTACGTACGT"
  b <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(a)))
  d_dim <- pc_dg(pc_ntthal("ANY", a, b)) - unname(duplex_dg_thermo(a, b, conditions = reaction_conditions()))
  expect_lt(abs(d_dim), 5)
})

# --- design on the same loci (slow; opt-in) ----------------------------------

# The 300 bp window, SNP at 1-based 150, Primer3 constrained to the package's
# amplicon window so the two designers are asked the same question.
pc_primer3 <- function(seq) {
  inp <- c(
    paste0("SEQUENCE_TEMPLATE=", seq),
    "SEQUENCE_TARGET=150,1",
    "PRIMER_TASK=generic",
    "PRIMER_PICK_LEFT_PRIMER=1", "PRIMER_PICK_RIGHT_PRIMER=1",
    "PRIMER_PICK_INTERNAL_OLIGO=0",
    "PRIMER_OPT_SIZE=22", "PRIMER_MIN_SIZE=18", "PRIMER_MAX_SIZE=25",
    "PRIMER_OPT_TM=60", "PRIMER_MIN_TM=57", "PRIMER_MAX_TM=63",
    "PRIMER_PRODUCT_SIZE_RANGE=50-100",
    "PRIMER_NUM_RETURN=5",
    "="
  )
  out <- suppressWarnings(system2("primer3_core", input = inp, stdout = TRUE, stderr = FALSE))
  get1 <- function(k) {
    hit <- grep(paste0("^", k, "="), out, value = TRUE)
    if (!length(hit)) NA_character_ else sub("^[^=]*=", "", hit[[1L]])
  }
  n <- suppressWarnings(as.integer(get1("PRIMER_PAIR_NUM_RETURNED")))
  list(
    n_pairs = if (is.na(n)) 0L else n,
    product = suppressWarnings(as.integer(get1("PRIMER_PAIR_0_PRODUCT_SIZE")))
  )
}

pc_panel <- function(n = 3L) {
  genome <- ath_genome()
  json <- fs::path(marker_cache_dir(), "arabidopsis_thaliana_var_4:1500000-1560000.json")
  if (!fs::file_exists(json)) {
    return(NULL)
  }
  recs <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  snv <- purrr::keep(recs, function(r) {
    !is.null(r$start) && !is.null(r$end) && r$end == r$start &&
      length(r$alleles) == 2L && all(nchar(unlist(r$alleles)) == 1L) &&
      all(toupper(unlist(r$alleles)) %in% c("A", "C", "G", "T"))
  })
  ok <- purrr::keep(snv, function(r) {
    ch <- paste0("Chr", r$seq_region_name)
    ch %in% names(genome) &&
      toupper(as.character(Biostrings::subseq(genome[[ch]], r$start, r$start))) ==
        toupper(r$alleles[[1L]])
  })
  if (length(ok) < n) {
    return(NULL)
  }
  panel <- purrr::map_dfr(ok[seq_len(n)], function(r) {
    al <- toupper(unlist(r$alleles))
    tibble::tibble(
      chrom = paste0("Chr", r$seq_region_name), pos = as.integer(r$start),
      ref = al[[1L]], alt = al[[2L]]
    )
  })
  panel$seq <- vapply(seq_len(nrow(panel)), function(i) {
    as.character(Biostrings::subseq(genome[[panel$chrom[[i]]]], panel$pos[[i]] - 149L, panel$pos[[i]] + 150L))
  }, character(1))
  list(genome = genome, panel = panel)
}

test_that("design_kasp_markers and Primer3 both solve the same three loci", {
  pc_skip_if_missing(c("primer3_core", "oligotm", "ntthal"))
  skip_if_not_full()
  skip_if_no_oligoarrayaux()
  dat <- pc_panel(3L)
  skip_if_not(is.null(dat), "the cached variation panel is not available")

  p3 <- lapply(dat$panel$seq, pc_primer3)
  # Primer3 must actually find pairs on this windows, or the comparison is empty
  expect_true(all(vapply(p3, function(x) x$n_pairs > 0L, logical(1L))))

  design <- suppressWarnings(design_kasp_markers(dat$genome, dat$panel[c("chrom", "pos", "ref", "alt")]))
  expect_equal(nrow(design), nrow(dat$panel))
  ok <- design$verdict != "fail"
  expect_gte(sum(ok), 1L)
  # every assay the package returns sits inside the window it declares
  expect_true(all(dplyr::between(
    design$amplicon[ok], KASP_LIMITS$amplicon[1L], KASP_LIMITS$amplicon[2L]
  )))

  # on the loci both solved, the products are comparable (not identical: the
  # two optimise different objectives)
  both <- ok & !is.na(vapply(p3, function(x) x$product, integer(1L)))
  if (any(both)) {
    expect_true(all(abs(design$amplicon[both] - vapply(p3, function(x) x$product, integer(1L))[both]) < 60L))
  }
})
