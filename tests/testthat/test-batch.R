# Batch design: one row per locus with a verdict; a bad locus must not abort.

test_that("variant_type labels SNV, MNP and indels", {
  expect_equal(variant_type("A", "G"), "SNP")
  expect_equal(variant_type("AC", "GT"), "MNP")
  expect_equal(variant_type("A", "AGG"), "insertion")
  expect_equal(variant_type("AGG", "A"), "deletion")
})

test_that("design_kasp_markers returns a tidy per-locus verdict without aborting", {
  skip_if_no_oligoarrayaux()
  W <- Biostrings::DNAString(bseq(500, 4))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  cand <- seq(160L, 340L, by = 1L)
  clean <- cand[vapply(cand, function(p) {
    !has_run(as.character(Biostrings::subseq(W, p - 35L, p + 35L)))
  }, logical(1))]
  p1 <- clean[1L]
  p2 <- clean[which(abs(clean - p1) > 60L)[1L]]
  ref1 <- as.character(Biostrings::subseq(W, p1, p1))
  alt1 <- setdiff(c("A", "C", "G", "T"), ref1)[1L]
  ref2 <- as.character(Biostrings::subseq(W, p2, p2))
  loci <- tibble::tibble(
    chrom = c("Chr1", "Chr1", "Chr1"),
    pos   = c(p1, p2, 5L), # last one is un-designable
    ref   = c(ref1, ref2, "A"),
    alt   = c(alt1, paste0(ref2, "GG"), "G")
  )
  res <- design_kasp_markers(genome, loci, flank = 150L)
  expect_equal(nrow(res), 3L)
  expect_true(all(c("chrom", "pos", "ref", "alt", "type", "verdict", "note") %in% names(res)))
  expect_equal(res$type, c("SNP", "insertion", "SNP"))
  # the un-designable edge locus is captured as a failure, not thrown ...
  expect_equal(res$verdict[3L], "fail")
  expect_true(is.na(res$strand[3L]))
  expect_true(nzchar(res$note[3L]))
  # ... and the verdict wiring is consistent: an assay exists iff the locus did
  # not fail, with finite coordinates and discrimination
  built <- res$verdict != "fail"
  expect_gte(sum(built), 1L)
  expect_false(any(is.na(res$strand[built])))
  expect_true(all(is.finite(res$discrimination[built])))
  expect_false(any(is.na(res$amplicon[built])))
})

test_that("a BiocParallel backend gives the same batch verdicts", {
  skip_if_no_oligoarrayaux()
  skip_if_not_installed("BiocParallel")
  W <- Biostrings::DNAString(bseq(400, 4))
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- as.character(Biostrings::subseq(W, 210L, 210L))
  loci <- tibble::tibble(
    chrom = "Chr1", pos = 210L, ref = ref,
    alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  serial <- suppressWarnings(design_kasp_markers(genome, loci, flank = 150L))
  skip_on_os("windows")
  skip_if_not_full()
  # an explicit SerialParam must agree with the default, and a SOCK pool with both
  par <- suppressWarnings(design_kasp_markers(genome, loci,
    flank = 150L,
    BPPARAM = BiocParallel::SerialParam()
  ))
  expect_equal(par, serial)
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  sock <- suppressWarnings(design_kasp_markers(genome, loci, flank = 150L, BPPARAM = bp))
  expect_equal(sock, serial)
  expect_equal(sock, par)
})

# -- design_marker_batch() --------------------------------------------------
# The decision graph a design carries is a `tbl_graph` whose print carries a
# per-build address, so two identical builds are not `identical()`.  Compare the
# graph by its vertices and edges, then `all.equal()` the rest.
batch_payload <- function(x) {
  if (inherits(x, "tbl_graph")) {
    return(list(
      v = as.data.frame(tidygraph::as_tibble(x, active = "nodes")),
      e = as.data.frame(tidygraph::as_tibble(x, active = "edges"))
    ))
  }
  if (is.list(x)) {
    out <- lapply(x, batch_payload)
    names(out) <- names(x)
    return(out)
  }
  x
}

# two SNPs that both design (260 was checked to be designable in this fixture),
# so the equivalence test is not vacuously comparing two doomed loci
batch_loci <- function() {
  fx <- snp_fixture(600L, seed = 2L, pos = 300L)
  tibble::tibble(
    chrom = c("Chr1", "Chr1"), pos = c(300L, 260L),
    ref = c(fx$ref, substr(fx$W, 260L, 260L)),
    alt = c(fx$alt, setdiff(c("A", "C", "G", "T"), substr(fx$W, 260L, 260L))[1L])
  )
}

# one fixture and one batch, reused by the tests below: a design is the
# expensive unit here, so each locus is designed once per comparison.  The
# batch carries an out-of-range locus in the middle so the same run also proves
# that a failing locus is captured and the batch goes on.
batch_cache <- new.env(parent = emptyenv())

batch_shared <- function() {
  if (is.null(batch_cache$batch)) {
    fx <- snp_fixture(600L, seed = 2L, pos = 300L)
    loci <- batch_loci()
    bad <- tibble::tibble(chrom = "Chr1", pos = 5L, ref = "A", alt = "G")
    batch_cache$fx <- fx
    batch_cache$loci <- loci
    # the good rows are distinguishable, so the ordering proves the resume
    batch_cache$mixed <- rbind(loci[1L, ], bad, loci[2L, ])
    batch_cache$batch <- design_marker_batch(fx$genome, batch_cache$mixed, design = "snp", flank = 100L)
  }
  batch_cache
}

test_that("design_marker_batch runs over many loci and keeps the design", {
  skip_if_no_oligoarrayaux()
  res <- batch_shared()$batch
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 3L)
  expect_true(all(c("chrom", "pos", "ref", "alt", "status", "scheme", "message", "elapsed", "result") %in% names(res)))
  expect_equal(res$status, c("ok", "fail", "ok"))
  expect_equal(res$scheme, rep("snp", 3L))
  expect_true(all(is.finite(res$elapsed)) && all(res$elapsed >= 0))
  # the design object itself comes back, not a summary of it
  expect_true(all(vapply(res$result[c(1L, 3L)], is_target_designable, logical(1))))
  expect_s3_class(res$result[[1L]], "marker_assay")
  expect_true(nzchar(res$result[[1L]]$assay$as_ref[1L]))
})

test_that("a locus that errors is captured, not thrown, and the batch goes on", {
  skip_if_no_oligoarrayaux()
  res <- batch_shared()$batch
  # the batch returns a row per locus rather than stopping at the error
  expect_equal(res$status, c("ok", "fail", "ok"))
  expect_true(nzchar(res$message[2L]))
  expect_null(res$result[[2L]])
  # the loci on both sides of the bad one still produced designs
  expect_true(is_target_designable(res$result[[1L]]))
  expect_true(is_target_designable(res$result[[3L]]))
})

test_that("design_marker_batch equals the same calls run one by one", {
  skip_if_no_oligoarrayaux()
  shared <- batch_shared()
  fx <- shared$fx
  batch <- shared$batch
  # the two good rows, run by hand on the canonical row the batch builds
  one_by_one <- lapply(c(1L, 3L), function(i) {
    locus <- tidymarkers:::locus_row(shared$mixed[i, , drop = FALSE])
    suppressWarnings(marker_design_at(fx$genome, locus, "snp", flank = 100L, cores = 1L))
  })
  expect_equal(batch_payload(batch$result[c(1L, 3L)]), batch_payload(one_by_one))
})

test_that("the default ladder matches a hand-run ladder", {
  skip_if_no_oligoarrayaux()
  shared <- batch_shared()
  locus <- shared$loci[1L, , drop = FALSE]
  ladder <- design_marker_batch(shared$fx$genome, locus, flank = 100L)
  # the batch hands the designer a canonical row (`alt2`/`strand` filled in), so
  # the hand run must start from the same row to be the same design
  hand <- suppressWarnings(
    design_auto_marker(shared$fx$genome, tidymarkers:::locus_row(locus), flank = 100L, cores = 1L)
  )
  # a ladder answer is an answer, so the row is always ok and carries the scheme
  expect_equal(ladder$status, "ok")
  expect_equal(ladder$scheme, hand$scheme)
  expect_true(nzchar(ladder$scheme))
  expect_equal(batch_payload(ladder$result[[1L]]), batch_payload(hand))
})

# -- one locus shape for every entry point ----------------------------------
# The multi-allelic row of the shipped panel: the one row whose `alt2` decides a
# design, so a path that drops it cannot pass by luck.
multi_locus <- function() {
  ex <- marker_example()
  ex[!is.na(ex$alt2), c("chrom", "pos", "ref", "alt", "alt2")][1L, ]
}

test_that("the batch and design_markers() design a multi-allelic locus the same", {
  skip_if_no_oligoarrayaux()
  skip_if_cold_demo_cache("genome")
  genome <- ath_genome()
  locus <- multi_locus()
  batch <- design_marker_batch(genome, locus, design = "snp", flank = 150L)
  ref <- design_markers(genome, locus, scheme = "snp", flank = 150L)
  expect_equal(batch$status, "ok")
  expect_equal(ref$ok, TRUE)
  # the second ALT reaches the design as the degenerate site it is, in both
  expect_equal(batch$result[[1L]]$site$class[1L], "degenerate")
  expect_equal(batch$result[[1L]]$site$n_alt[1L], 2L)
  expect_setequal(marker_snp_alts(locus)$alt_pass, c(locus$alt, locus$alt2))
  expect_equal(ref$scheme, "snp")
})

test_that("a locus table without alt2/strand is the same as one that names them", {
  skip_if_no_oligoarrayaux()
  skip_if_cold_demo_cache("genome")
  genome <- ath_genome()
  bare <- multi_locus()[, c("chrom", "pos", "ref", "alt")]
  named <- bare
  named$alt2 <- NA_character_
  named$strand <- "+"
  expect_identical(design_markers(genome, bare, flank = 150L), design_markers(genome, named, flank = 150L))
  # and the batch sees the same row either way
  a <- design_marker_batch(genome, bare, flank = 150L)$result[[1L]]
  b <- design_marker_batch(genome, named, flank = 150L)$result[[1L]]
  expect_equal(batch_payload(a), batch_payload(b))
})

test_that("design_kasp_markers' verdict schema does not see alt2/strand", {
  skip_if_no_oligoarrayaux()
  fx <- snp_fixture()
  res <- suppressWarnings(design_kasp_markers(fx$genome, fx$snp, flank = 150L))
  # the verdict table is the one HEAD returns: same columns in the same order,
  # one row, and the same values -- the canonical row adds `alt2`/`strand` for
  # the designer, but neither reaches this output
  expect_identical(
    names(res),
    c(
      "chrom", "pos", "ref", "alt", "type", "strand", "amplicon",
      "discrimination", "n_fail", "verdict", "note"
    )
  )
  expect_equal(nrow(res), 1L)
  expect_equal(res$chrom, fx$snp$chrom)
  expect_equal(res$pos, fx$snp$pos)
  expect_equal(res$type, "SNP")
  expect_equal(res$strand, "+")
  expect_equal(res$verdict, "warn")
  expect_equal(res$n_fail, 1L)
})
