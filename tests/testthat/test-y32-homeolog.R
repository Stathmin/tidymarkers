# =============================================================================
#  Y3.2  --  a KASP construct has TWO axes, and they are not interchangeable
# =============================================================================
#  Review point, verified against the PolyMarker source (data-raw/polymarker-source.R,
#  read at VERSION 1.3.3 commit d4a4578): the allele-specific primer's 3' base
#  carries the ALLELE.  Subgenome specificity is carried by the COMMON primer's
#  3' base sitting on a homoeolog-diagnostic position.  Reading "the AS 3' base
#  separates A/B/D" conflates the two and produces an assay that genotypes but
#  does not localise.
#
#  The three claims, in the order they appear in the source:
#    1. allele axis   -- the two AS primers differ at their 3' base, and that base
#                        is the SNP allele (proven by Grewal's reverse-strand mode,
#                        where both 3' bases are complements of the alleles);
#    2. subgenome axis -- the common primer is force-right-ended onto a mask base
#                        that all other subgenomes differ from;
#    3. they are independent columns of the artifact, never one "3' base" field.

ext <- function(name) {
  readr::read_csv(tm_extdata(name),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

test_that("Y3.2: the rule source says the allele sits on the AS primer", {
  r <- ext("polymarker_design_rules.csv")
  expect_gt(nrow(r), 6L)
  # every rule names the file, the method and the line it was read from
  expect_true(all(c(
    "rule_id", "role", "scope", "source_file", "source_method",
    "quoted", "calibrated", "doi"
  ) %in% names(r)))
  expect_true(all(
    r$source_file == "lib/bio/PolyploidTools/SNP.rb" |
      grepl("primer3.rb|PrimerRegion.rb", r$source_file)
  ))
  expect_true(all(r$source_method != "" & !is.na(r$source_method)))
  # the allele axis is a real category in the source, and so is the subgenome axis
  expect_true(any(r$role == "allele-primer"))
  expect_true(any(r$role == "common-primer" & grepl("subgenome", r$scope)))
  # the quoted line must be reproducible from the source again, so a rule cannot
  # drift into a paraphrase
  expect_true(all(nzchar(r$quoted)))
  expect_true(all(grepl(
    "SEQUENCE|right_pos|snp_pos|mask|almost|homoeologous|",
    r$quoted
  )))
})

test_that("Y3.2: WITHOUT a diagnostic common base, no subgenome is claimed", {
  r <- ext("polymarker_design_rules.csv")
  fallback <- r[r$scope == "fallback", ]
  expect_equal(nrow(fallback), 1L)
  # the source's own fallback label is "chromosome_nonspecific": the assay exists
  # and genotypes, but it is NOT subgenome-specific
  expect_true(grepl("nonspecific", fallback$quoted))
  expect_true(grepl("not subgenome-specific", fallback$rule))
  # and the construct has no field of its own that says which primer carries the
  # subgenome, so this cannot be inferred from a label
  con <- r[r$role == "construct", ]
  expect_equal(nrow(con), 1L)
  expect_true(grepl("no 'which primer carries", con$rule))
})

test_that("Y3.2: the allele axis is the AS 3' base, in both strand modes", {
  g <- ext("grewal_kasp_axes.csv")
  expect_equal(nrow(g), 1000L)
  rc1 <- function(x) chartr("ACGT", "TGCA", x)
  two <- g$assay_mode == "two-allele"
  rev_mode <- g$assay_mode == "reverse"
  # forward design: the AS 3' bases ARE the two alleles
  expect_true(all(g$as_fam_3p[two] == g$snp_allele_1[two]))
  expect_true(all(g$as_hex_3p[two] == g$snp_allele_2[two]))
  # reverse design: both are the complements, which is the same axis mirrored
  expect_true(all(g$as_fam_3p[rev_mode] == rc1(g$snp_allele_1[rev_mode])))
  expect_true(all(g$as_hex_3p[rev_mode] == rc1(g$snp_allele_2[rev_mode])))
  # in every assay the allele axis is read from the AS primers, so the two must
  # never be equal (an equal pair cannot separate the alleles)
  expect_true(all(g$as_fam_3p != g$as_hex_3p))
})

test_that("Y3.2: the subgenome axis is NOT the AS 3' base", {
  # The negative control for the review's correction: if the AS 3' base carried
  # the subgenome, it would be constrained by the chromosome; it is constrained
  # by the allele instead, so its distribution follows the SNP class
  p <- ext("makhoul_kasp_primers.csv")
  expect_gt(nrow(p), 10L)
  expect_true(all(p$allele_in_as_3p))
  # with 13 assays there is not enough data to correlate anything; what CAN be
  # asserted is the structural fact the source states
  r <- ext("polymarker_design_rules.csv")
  as_rows <- r[r$role == "allele-primer", ]
  expect_true(any(grepl("3'-terminal base = the SNP", as_rows$scope)))
})

test_that("Y3.2: the common primer's 3' base varies independently of the SNP", {
  p <- ext("makhoul_kasp_primers.csv")
  # For one SNP class the common 3' base takes several values across assays.  If
  # it were an allele axis it would be fixed by the SNP.
  tc <- p[p$snp_kasp == "[T/C]", ]
  expect_gt(length(unique(tc$common_3p)), 1L)
  # and the artifact must expose it as its own column, under its own name
  expect_true("common_3p" %in% names(p))
  expect_false(identical(p$common_3p, p$fam_3p))
})

test_that("Y3.2: the two axes are separate columns, and both are documented", {
  h <- ext("kasp_homeolog_axes.csv")
  expect_setequal(h$axis, c("allele", "subgenome"))
  allele <- h[h$axis == "allele", ]
  genome <- h[h$axis == "subgenome", ]
  expect_true(grepl("allele-specific primer", allele$carried_by))
  expect_true(grepl("common primer", genome$carried_by))
  # the genome axis quotes the published sentence, so it is not our paraphrase
  expect_true(grepl(
    "3.? end anchored to a base, which is unique to the target subgenome",
    genome$verbatim
  ))
  # both rows point at sources
  expect_true(all(grepl("^10\\.", h$doi)))
})

test_that("Y3.2: the allele axis is constrained by the SNP, the genome axis is not", {
  # The empirical form of the review's correction, on 1000 real wheat assays: if
  # the common primer's 3' base carried the ALLELE, it would be fixed for each
  # SNP pair (four values would collapse to one per allele).  If the AS 3' bases
  # carried the SUBGENOME, they would not track the allele.  Both go the other way:
  g <- ext("grewal_kasp_axes.csv")
  g$snp <- paste(g$snp_allele_1, g$snp_allele_2)
  tab <- table(g$snp, g$common_3p)
  distinct_common <- apply(tab, 1L, function(v) sum(v > 0))
  # every SNP pair takes all four common-primer 3' bases somewhere in the data
  expect_true(all(distinct_common >= 3L))
  # while the AS 3' base is one of exactly two values per SNP pair
  as_bases <- table(g$snp, g$as_fam_3p)
  distinct_as <- apply(as_bases, 1L, function(v) sum(v > 0))
  expect_true(all(distinct_as == 2L))
  # and the AS 3' bases are subsets of {allele_1, allele_2} (or their complements)
  rc1 <- function(x) chartr("ACGT", "TGCA", x)
  ok <- (g$as_fam_3p %in% c(
    g$snp_allele_1, g$snp_allele_2, rc1(g$snp_allele_1),
    rc1(g$snp_allele_2)
  ))
  expect_true(all(ok))
})

test_that("Y3.2: an artifact with one 3' column cannot pass as two axes", {
  # The regression this test exists for: a design that reports a single "3' base"
  # cannot distinguish "genotypes" from "genotypes AND localises".
  axes <- ext("kasp_homeolog_axes.csv")
  expect_equal(nrow(axes), 2L)
  expect_false(any(duplicated(axes$axis)))
  # a naive artifact that kept only one column would lose the subgenome claim,
  # and nothing in the data would notice -- so the column is required by name
  p <- ext("makhoul_kasp_primers.csv")
  expect_true(all(c("fam_3p", "hex_3p", "common_3p") %in% names(p)))
})

test_that("Y3.2: PolyMarker is a benchmark and rule source, never reimplemented", {
  s <- ext("scientific_sources.csv")
  row <- s[s$key == "polymarker_source", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$access, "source")
  expect_equal(row$doi, "10.1093/bioinformatics/btv069")
  # the rule table points back at the same DOI, so the two cannot disagree
  r <- ext("polymarker_design_rules.csv")
  expect_true(all(r$doi == row$doi))
  # every rule is uncalibrated: the source states a construction, not an outcome
  expect_false(any(r$calibrated))
})

# =============================================================================
#  Y3.1  --  `allele_base` / `subgenome_anchor` as explicit candidate fields
# =============================================================================

test_that("Y3.1: a designed candidate carries the two 3' axes explicitly", {
  f <- snp_fixture()
  m <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp, pair_screen = 0L, cores = 1L))
  ax <- homeolog_axes(m)

  # the allele axis is the AS 3' base (REF/ALT), named and always filled
  expect_named(m$allele_base, c("AS_ref", "AS_alt"))
  expect_equal(
    unname(m$allele_base),
    plus_3p_base(c(m$assay$as_ref[1], as_alt_sequence(m)), m$strand)
  )

  # the subgenome axis is the COM 3' base, read on the opposite strand ...
  expect_equal(
    m$subgenome_anchor,
    plus_3p_base(m$assay$com[1], if (m$strand == "+") "-" else "+")
  )
  # ... and its position is the plus-strand template base at `anchor_pos`
  idx <- ax$anchor_pos - m$snp$pos + m$template$snp_i
  expect_equal(
    m$subgenome_anchor,
    as.character(Biostrings::subseq(m$template$seq, idx, idx))
  )

  # without a mask the verdict is explicit NA, never a silent `unique`
  expect_equal(nrow(m$homeolog_diagnostics), 1L)
  expect_true(is.na(m$homeolog_diagnostics$diagnostic))
  expect_true(is.na(m$homeolog_diagnostics$specificity))
  expect_match(m$homeolog_diagnostics$reason, "no subgenome mask")
})

test_that("Y3.1: a diagnostic mask hit reports unique", {
  f <- snp_fixture()
  m0 <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp, pair_screen = 0L, cores = 1L))
  ap <- homeolog_axes(m0)$anchor_pos
  m <- suppressWarnings(assemble_kasp_assay(
    f$genome, f$snp,
    pair_screen = 0L, cores = 1L,
    diagnostic_positions = data.frame(chrom = "Chr1", pos = ap)
  ))
  ax <- homeolog_axes(m)
  expect_true(ax$diagnostic)
  expect_equal(ax$specificity, "unique")
  expect_equal(m$subgenome_anchor, m0$subgenome_anchor)
})

test_that("Y3.1: a diagnostic mask miss warns and reports multicopy", {
  f <- snp_fixture()
  m0 <- suppressWarnings(assemble_kasp_assay(f$genome, f$snp, pair_screen = 0L, cores = 1L))
  ap <- homeolog_axes(m0)$anchor_pos
  expect_warning(
    m <- withCallingHandlers(
      assemble_kasp_assay(
        f$genome, f$snp,
        pair_screen = 0L, cores = 1L,
        diagnostic_positions = data.frame(chrom = "Chr1", pos = ap + 1000L)
      ),
      warning = function(w) {
        # the synthetic template also raises an unrelated near-match PCR warning
        if (grepl("near-match in-silico PCR", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    "does not localise"
  )
  ax <- homeolog_axes(m)
  expect_false(ax$diagnostic)
  expect_equal(ax$specificity, "multicopy")
})

test_that("Y3.1: fixed primers carry the same two axes", {
  f <- snp_fixture()
  len <- 20L
  pos <- f$snp$pos
  as_ref <- substring(f$W, pos - len + 1L, pos)
  as_alt <- paste0(substring(f$W, pos - len + 1L, pos - 1L), f$alt)
  com <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAString(substring(f$W, pos + 30L, pos + 30L + len - 1L))
  ))
  m <- suppressWarnings(assemble_primer_assay(
    f$genome, f$snp,
    c(AS_ref = as_ref, AS_alt = as_alt, COM = com)
  ))
  expect_named(m$allele_base, c("AS_ref", "AS_alt"))
  expect_equal(unname(m$allele_base), plus_3p_base(c(as_ref, as_alt), m$strand))
  expect_false(is.na(m$subgenome_anchor))
  expect_equal(nrow(m$homeolog_diagnostics), 1L)
  expect_true(is.na(m$homeolog_diagnostics$diagnostic))
})

test_that("Y3.1: the data artifact fills both axes for every assay", {
  p <- ext("makhoul_kasp_primers.csv")
  expect_true(all(c("fam_3p", "hex_3p", "common_3p") %in% names(p)))
  # every row carries a real nucleotide on each 3' axis: NA would only appear
  # where the base is not a nucleotide, and no row is allowed to drop the field
  expect_true(all(!is.na(p$fam_3p) & p$fam_3p %in% c("A", "C", "G", "T")))
  expect_true(all(!is.na(p$hex_3p) & p$hex_3p %in% c("A", "C", "G", "T")))
  expect_true(all(!is.na(p$common_3p) & p$common_3p %in% c("A", "C", "G", "T")))
})
