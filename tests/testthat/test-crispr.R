# CRISPR-Dx (SHERLOCK / DETECTR): the guide design and the classed marker.

test_that("design_crispr_marker builds a Cas12a guide with a TTTV PAM and an amplicon", {
  skip_if_no_oligoarrayaux()
  W <- bseq(600L, 7L)
  substr(W, 293L, 296L) <- "TTTA" # 5'-TTTV PAM immediately 5' of the protospacer
  substr(W, 297L, 316L) <- "ACGTGCTAGCTAGCATGCA" # GC-balanced 20 nt spacer
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "T", alt = "A")
  expect_warning(
    m <- design_crispr_marker(genome, snp, flank = 200L, cas = "cas12a"),
    "guide/PAM conservation is not checked"
  )
  expect_s3_class(m, "marker_crispr")
  expect_equal(m$cas, "cas12a")
  expect_match(m$pam, "^TTT[ACG]$")
  expect_true(is.na(m$pfs)) # Cas12a uses a DNA PAM, not a Cas13 PFS
  expect_true(is.na(m$pfs_side))
  expect_true(m$gates[["gate_crispr_pam"]])
  expect_gte(nchar(m$guide), 20L)
  expect_lte(nchar(m$guide), 23L)
  expect_lt(m$seed_offset, 8L) # the variant sits in the PAM-proximal seed
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  expect_gte(m$amplicon$product, 80L)
  expect_lte(m$amplicon$product, 250L)
  expect_output(print(m), "marker_crispr")
})

test_that("design_crispr_marker builds a Cas13a guide with a 3' non-G PFS", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  m <- suppressWarnings(design_crispr_marker(f$genome, f$snp, flank = 200L, cas = "cas13"))
  expect_s3_class(m, "marker_crispr")
  expect_equal(m$cas, "cas13")
  expect_true(is.na(m$pam)) # no DNA PAM: Cas13 targets RNA
  expect_equal(m$pfs_side, "3'")
  expect_false(is.na(m$pfs))
  expect_false(identical(m$pfs, "G")) # the 3' PFS avoids G
  expect_true(m$gates[["gate_crispr_pam"]]) # the PFS licensing gate is evaluated
  expect_false(is.na(m$amplicon$product)) # ... and the marker still designs
  # LwaCas13a uses ~28-nt spacers: the default window is 26-28 nt.
  expect_gte(nchar(m$guide), 26L)
  expect_lte(nchar(m$guide), 28L)
  expect_false("gate_crispr_spacer" %in% m$relaxed)
  expect_output(print(m), "PFS")
})

test_that("cas13 and cas13a are the same Cas13a rule, while cas13b is explicit", {
  expect_equal(crispr_default_spacer("cas12a"), 20:23)
  expect_equal(crispr_default_spacer("cas13"), 26:28)
  expect_equal(crispr_default_spacer("cas13a"), 26:28)
  expect_equal(crispr_default_spacer("cas13b"), 28:30)
  expect_error(crispr_default_spacer("cas9"))

  s <- strrep("ACGT", 40L)
  alias <- crispr_guide_candidates(s, 80L, "cas13", 28L, c(0.35, 0.65), "+")
  explicit <- crispr_guide_candidates(s, 80L, "cas13a", 28L, c(0.35, 0.65), "+")
  expect_equal(alias, explicit)
  expect_true(all(alias$pfs_side == "3'"))
})

test_that("Cas13a rejects a G at the 3' flank of the plus-strand protospacer", {
  s <- strrep("ACGT", 40L)
  # len 28, offset 10: plus-strand protospacer 70..97, 3' flank at 98.
  substr(s, 98L, 98L) <- "G"
  cand <- crispr_guide_candidates(s, 80L, "cas13a", 28L, c(0.35, 0.65), "+")
  expect_false(any(cand$start == 70L)) # the G-PFS candidate is filtered out
  expect_true(all(cand$pfs != "G"))
  expect_true(all(is.na(cand$pam)))

  substr(s, 98L, 98L) <- "A"
  cand <- crispr_guide_candidates(s, 80L, "cas13a", 28L, c(0.35, 0.65), "+")
  expect_true(any(cand$start == 70L))
  expect_equal(cand$pfs[cand$start == 70L], "A")
  expect_true(all(cand$pfs_side == "3'"))
})

test_that("Cas13a reads the 3' PFS from the guide strand on the minus strand", {
  s <- strrep("ACGT", 40L)
  # len 28, offset 10 on "-": plus-strand protospacer 63..90, so the 3' flank of
  # the target RNA is complement(s[62]); C on the plus strand gives G in the RNA.
  substr(s, 62L, 62L) <- "C"
  cand <- crispr_guide_candidates(s, 80L, "cas13a", 28L, c(0.35, 0.65), "-")
  expect_false(any(cand$start == 63L))

  substr(s, 62L, 62L) <- "A" # complement T: an allowed 3' PFS
  cand <- crispr_guide_candidates(s, 80L, "cas13a", 28L, c(0.35, 0.65), "-")
  expect_true(any(cand$start == 63L))
  expect_equal(cand$pfs[cand$start == 63L], "T")
})

test_that("Cas13b enforces its 5' D-PFS (A/U/G, i.e. not C) on both strands", {
  s <- strrep("ACGT", 40L)
  # plus: len 28, offset 10 -> 70..97, the 5' flank is s[69].
  substr(s, 69L, 69L) <- "C"
  cand <- crispr_guide_candidates(s, 80L, "cas13b", 28L, c(0.35, 0.65), "+")
  expect_false(any(cand$start == 70L))
  substr(s, 69L, 69L) <- "A"
  cand <- crispr_guide_candidates(s, 80L, "cas13b", 28L, c(0.35, 0.65), "+")
  expect_true(any(cand$start == 70L))
  expect_equal(cand$pfs[cand$start == 70L], "A")
  expect_true(all(cand$pfs_side == "5'"))
  expect_true(all(cand$pfs != "C"))

  # minus: 63..90, the 5' flank of the target RNA is complement(s[91]).
  substr(s, 91L, 91L) <- "G" # complement C -> rejected
  cand <- crispr_guide_candidates(s, 80L, "cas13b", 28L, c(0.35, 0.65), "-")
  expect_false(any(cand$start == 63L))
  substr(s, 91L, 91L) <- "A" # complement T -> accepted
  cand <- crispr_guide_candidates(s, 80L, "cas13b", 28L, c(0.35, 0.65), "-")
  expect_true(any(cand$start == 63L))
  expect_equal(cand$pfs[cand$start == 63L], "T")
})

test_that("Cas13b warns that its 3' NAN/NNA PFS is not enforced", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  expect_warning(
    expect_warning(
      m <- design_crispr_marker(f$genome, f$snp, flank = 200L, cas = "cas13b"),
      "3' PFS"
    ),
    "guide/PAM conservation is not checked"
  )
  expect_s3_class(m, "marker_crispr")
  expect_equal(m$cas, "cas13b")
  expect_equal(m$pfs_side, "5'")
  expect_false(identical(m$pfs, "C"))
  expect_gte(nchar(m$guide), 28L)
  expect_lte(nchar(m$guide), 30L)
  expect_false("gate_crispr_spacer" %in% m$relaxed)
})

test_that("design_crispr_marker refuses a spacer length it would have truncated", {
  W <- bseq(600L, 7L)
  substr(W, 293L, 296L) <- "TTTA"
  substr(W, 297L, 316L) <- "ACGTGCTAGCTAGCATGCA"
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "T", alt = "A")

  # 20.9 used to become 20, so the caller got a guide one base shorter
  for (bad in list(c(20.9, 21.5), 20.5, 0L, -1L, NA, Inf, c(20L, NA))) {
    expect_error(
      design_crispr_marker(genome, snp, spacer_lengths = bad),
      "positive whole numbers of nucleotides"
    )
  }
  # the same whole series still designs a guide of a named length
  skip_if_no_oligoarrayaux()
  expect_warning(
    m <- design_crispr_marker(genome, snp, flank = 200L, cas = "cas12a", spacer_lengths = 20:23),
    "guide/PAM conservation is not checked"
  )
  expect_true(nchar(m$guide) %in% 20:23)
})

test_that("design_crispr_marker informs and returns NULL when no guide passes", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("AT", 300L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "C")
  expect_message(m <- design_crispr_marker(genome, snp, flank = 200L, cas = "cas12a"), "guide")
  expect_null(m)
})

test_that("crispr_guide_candidates honours the strand restriction", {
  W <- bseq(300L, 3L)
  args <- list(seq = W, snp_i = 150L, cas = "cas13", spacer_lengths = 22:28, gc_range = c(0.35, 0.65))
  both <- do.call(crispr_guide_candidates, args)
  plus <- do.call(crispr_guide_candidates, c(args, list(strands = "+")))
  minus <- do.call(crispr_guide_candidates, c(args, list(strands = "-")))
  expect_setequal(both$strand, c("+", "-"))
  expect_true(all(plus$strand == "+"))
  expect_true(all(minus$strand == "-"))
  expect_equal(nrow(plus) + nrow(minus), nrow(both))
})

test_that("the Cas12a readout is documented as ssDNase on an ssDNA reporter", {
  def <- scheme_def_crispr()
  expect_match(def$discrimination, "ssDNase")
  expect_match(def$discrimination, "ssDNA reporter")
  expect_match(def$readout, "ssDNA")
  expect_false(grepl("dsDNA", def$discrimination))
  expect_false(grepl("dsDNA", def$readout))
})

test_that("scheme_def_crispr declares its prefixed gates, tiers and lookups", {
  def <- scheme_def_crispr()
  expect_equal(def$scheme, "crispr")
  expect_equal(def$family, "crispr")
  expect_true(all(grepl("^gate_crispr_", def$gates)))
  expect_setequal(def$gates, names(def$gate_short))
  expect_setequal(def$gates, names(def$gate_limit))
  expect_true(all(def$relax %in% def$gates))
  expect_setequal(names(def$tiers), c("guide", "seed", "specificity", "amplicon"))
  expect_true(all(names(def$gate_floor_label) %in% def$relax))
  expect_true(all(names(def$gate_hard_label) %in% setdiff(def$gates, def$relax)))
})

test_that("the crispr scheme describes Cas13 by its PFS, not as simply PAM-free", {
  def <- scheme_def_crispr()
  expect_match(def$gate_short[["gate_crispr_pam"]], "PFS")
  expect_match(def$gate_limit[["gate_crispr_pam"]], "3' non-G")
  expect_match(def$gate_limit[["gate_crispr_pam"]], "5' D")
  expect_match(def$gate_limit[["gate_crispr_spacer"]], "26-28")
  expect_false(grepl("PAM-free", def$gate_short[["gate_crispr_pam"]]))
})

test_that("the CRISPR seed gate re-checks the guide's discriminating region", {
  # Cas12a: the PAM-proximal 8 nt are the seed
  expect_true(crispr_seed_gate("cas12a", 20L, 7L))
  expect_false(crispr_seed_gate("cas12a", 20L, 8L))
  expect_false(crispr_seed_gate("cas12a", 20L, NA_integer_))
  # Cas13: the variant must be in the central window, not at an end
  expect_true(crispr_seed_gate("cas13", 28L, 14L))
  expect_false(crispr_seed_gate("cas13", 28L, 0L))
})

test_that("guide specificity is a hard gate", {
  # an off-target protospacer could license the reporter on the wrong template,
  # so this gate is never relaxed (unlike the amplicon / GC / spacer gates).
  expect_true(
    "gate_crispr_specificity" %in% setdiff(scheme_def_crispr()$gates, scheme_def_crispr()$relax)
  )
})

test_that("design_crispr_marker supports only the RPA pre-amplification", {
  # `lamp` was accepted but only reused the flanking-pair designer, which is not
  # a LAMP reaction (no FIP/BIP): the choice is gone until a LAMP engine exists.
  expect_false(any(grepl("lamp", deparse(formals(design_crispr_marker)$preamp))))
  expect_identical(formals(design_crispr_marker)$preamp, "rpa")
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("AT", 300L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "A", alt = "C")
  expect_error(design_crispr_marker(genome, snp, preamp = "lamp"), "rpa")
})

test_that("the crispr DAG obeys the contract", expect_dag_invariants("crispr"))

test_that("design_crispr_marker ranks its gated guides and returns the first as the pick", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  m <- suppressWarnings(design_crispr_marker(f$genome, f$snp, flank = 200L, cas = "cas13"))
  cd <- m$candidates
  expect_s3_class(cd, "tbl_df")
  expect_gt(nrow(cd), 1L)
  expect_true(all(c("spacer", "off_target", "chosen", "not_chosen") %in% names(cd)))
  # one chosen, and it is the first row: the pick is the best-first candidate
  expect_equal(sum(cd$chosen), 1L)
  expect_true(cd$chosen[1L])
  expect_equal(cd$spacer[1L], m$guide)
  expect_equal(m$off_target, cd$off_target[1L])
  # ordering is the enumerator's (central then GC-closest), not re-sorted here
  tmpl <- tidymarkers:::snp_template(f$genome, f$snp$chrom, f$snp$pos, 200L)
  fresh <- crispr_guide_candidates(
    tmpl$seq, tmpl$snp_i, "cas13", crispr_default_spacer("cas13"), c(0.35, 0.65)
  )
  hits <- vapply(fresh$spacer, crispr_spacer_hits, integer(1L), genome = f$genome)
  usable <- fresh[hits - 1L <= 0L, ]
  expect_equal(cd$start, usable$start)
  expect_equal(cd$strand, usable$strand)
  # the specificity gate is hard: no returned candidate carries an off-target
  expect_true(all(cd$off_target[cd$chosen] <= 0L))
})

test_that("an off-target protospacer is excluded, or explains why a guide was not chosen", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  W <- f$W
  # a second exact copy of a guide spacer elsewhere in the template makes that
  # guide non-specific; it must never come back as the chosen guide
  tmpl <- tidymarkers:::snp_template(f$genome, f$snp$chrom, f$snp$pos, 200L)
  probes <- crispr_guide_candidates(
    tmpl$seq, tmpl$snp_i, "cas13", crispr_default_spacer("cas13"), c(0.35, 0.65)
  )
  sp <- probes$spacer[nrow(probes)]
  copy_at <- max(probes$end) + 20L
  if (copy_at + nchar(sp) <= nchar(W)) {
    substr(W, copy_at, copy_at + nchar(sp) - 1L) <- sp
  }
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  m <- suppressWarnings(design_crispr_marker(genome, f$snp, flank = 200L, cas = "cas13"))
  skip_if(is.null(m), "no guide survived the duplicated spacer")
  expect_false(identical(m$guide, sp))
  expect_lte(m$off_target, 0L)
  expect_true(all(m$candidates$off_target[m$candidates$chosen] <= 0L))
  # any candidate that was not chosen carries a named reason
  rejected <- m$candidates[!m$candidates$chosen, ]
  expect_true(all(!is.na(rejected$not_chosen)))
  expect_true(all(
    rejected$not_chosen %in%
      c("lower_ranked", "no_preamp_pair", "off_target_protospacer", "specificity_unscanned")
  ))
  ot_flagged <- which(m$candidates$not_chosen == "off_target_protospacer")
  expect_true(all(m$candidates$off_target[ot_flagged] > 0L))
})

test_that("the CRISPR readout stays unsimulated: no reporter kinetics model", {
  skip_if_no_oligoarrayaux()
  f <- snp_fixture()
  m <- suppressWarnings(design_crispr_marker(f$genome, f$snp, flank = 200L, cas = "cas13"))
  # no guide -> activity dataset exists, so a predicted signal would be invented
  expect_false(m$simulated)
  expect_equal(m$kinetics, "unknown")
  expect_null(m$readout)
  r <- marker_readout(m)
  expect_false(r$simulated)
  expect_match(r$note, "kinetics not modelled")
  expect_false(crispr_seed_evidence()$calibrated)
  expect_equal(crispr_seed_evidence()$magnitude, "unknown - no model predicts the ratio for a target")
})
