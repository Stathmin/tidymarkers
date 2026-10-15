# TaqMan 5'-nuclease allelic discrimination: one flanking pair + two probes.

test_that("design_taqman_marker returns a short pair with two allele-specific probes", {
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = ref, alt = alt)

  m <- design_taqman_marker(genome, snp, flank = 150L)
  expect_s3_class(m, "kasp_taqman")
  expect_true(all(c("a", "b", "kind", "run", "pass") %in% names(m$compatibility)))
  # a short amplicon straddling the variant
  expect_lte(m$amplicon$product, 150L)
  expect_lt(m$amplicon$start, 300L)
  expect_gt(m$amplicon$end, 300L)
  # exactly one probe per allele, different fluorophores
  expect_equal(nrow(m$probes), 2L)
  expect_setequal(m$probes$allele, c("REF", "ALT"))
  expect_equal(length(unique(m$probes$fluorophore)), 2L)
  # the two probes are the same window, differing only at the variant base
  expect_equal(m$probes$start[1L], m$probes$start[2L])
  expect_equal(m$probes$end[1L], m$probes$end[2L])
  # the REF probe carries the reference base at the variant position
  vp <- m$probes$variant_pos[1L]
  expect_equal(substr(m$probes$probe[m$probes$allele == "REF"], vp, vp), ref)
  expect_equal(substr(m$probes$probe[m$probes$allele == "ALT"], vp, vp), alt)
  # whole probes are ACGT, 5'-G free, and inside the amplicon
  expect_true(all(grepl("^[ACGT]+$", m$probes$probe)))
  expect_true(all(substr(m$probes$probe, 1L, 1L) != "G"))
  expect_true(all(m$probes$start >= m$amplicon$start & m$probes$end <= m$amplicon$end))
  expect_output(print(m), "kasp_taqman")
  # the method itself (generic dispatch needs the package NAMESPACE, regenerated
  # by roxygen once the scheme is wired in)
  expect_output(print.kasp_taqman(m), "<kasp_taqman>")
  expect_output(print.kasp_taqman(m), "FAM")
  expect_output(print.kasp_taqman(m), "HEX")
})

test_that("probes are centred over the variant (maximal destabilisation)", {
  W <- bseq(600L, 2L)
  genome <- Biostrings::DNAStringSet(c(Chr1 = W))
  ref <- substr(W, 300L, 300L)
  snp <- tibble::tibble(
    chrom = "Chr1", pos = 300L, ref = ref, alt = setdiff(c("A", "C", "G", "T"), ref)[1L]
  )
  m <- design_taqman_marker(genome, snp, flank = 150L)
  expect_gte(m$probes$centrality[1L], 0.5)
  # centrality 1 = exactly central; the variant must sit near, not at, a probe end
  center <- (nchar(m$probes$probe[1L]) + 1L) / 2
  expect_lte(abs(m$probes$variant_pos[1L] - center), 2)
})

test_that("design_taqman_marker refuses non-substitutions with a message", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = bseq(400L, 3L)))
  snp <- tibble::tibble(chrom = "Chr1", pos = 200L, ref = "AC", alt = "GT")
  expect_null(suppressMessages(design_taqman_marker(genome, snp)))
})

test_that("the taqman scheme definition has the gates, tiers and labels wired", {
  d <- scheme_def_taqman()
  expect_equal(d$scheme, "taqman")
  expect_equal(d$family, "taqman")
  expect_equal(d$loci, 1L)
  expect_named(d, c(
    "scheme", "family", "label", "discrimination", "readout", "zygosity",
    "amplicon", "loci", "verdict", "gates", "relax", "tiers",
    "gate_short", "gate_limit", "gate_floor_label", "gate_hard_label", "stage_label"
  ))
  expect_true(all(grepl("^gate_taqman_", d$gates)))
  expect_true(all(d$relax %in% d$gates))
  expect_setequal(names(d$gate_short), d$gates)
  expect_setequal(names(d$gate_limit), d$gates)
  expect_true(all(names(d$gate_floor_label) %in% d$relax))
  expect_true(all(names(d$gate_hard_label) %in% setdiff(d$gates, d$relax)))
  expect_equal(unname(d$tiers[1L]), "probe discrimination (variant centrality)")
})
