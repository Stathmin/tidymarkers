# =============================================================================
#  Z8.5  --  analytical performance is reported one quantity at a time
# =============================================================================
#  The claim this file guards is narrow and load-bearing: "detectable" is only
#  meaningful with the quantity it names.  A dPCR assay reports a limit of blank,
#  a limit of detection and a limit of quantification against different
#  denominators, and folding them into one "sensitivity" makes the claim
#  unfalsifiable.  So the report names each quantity on its own line, the
#  artefact's decision row names them, and the sourced JAK2 table keeps them as
#  three rows with three definitions.
# =============================================================================

read_ext <- function(name) {
  readr::read_csv(tm_extdata(name), show_col_types = FALSE, na = c("", "NA"))
}

test_that("Z8.5 the dPCR report lists each analytical quantity separately", {
  txt <- marker_decision_text("dpcr")
  h <- which(grepl("^analytical performance", txt))
  expect_length(h, 1L)
  expect_true(grepl("never one 'sensitivity'", txt[[h]], fixed = TRUE))

  # one line per quantity, quantity name right-padded then its own definition
  block <- txt[(h + 1L):(h + 6L)]
  q <- sub("^  ([A-Za-z_]+) +.*$", "\\1", block)
  expect_setequal(q, c(
    "LoB", "LoD", "LoQ", "MAF", "false_positive_rate", "partition_count"
  ))

  # the three dPCR limits are present and do not share a definition
  defs <- sub("^  [A-Za-z_]+ +", "", block[match(c("LoB", "LoD", "LoQ"), q)])
  expect_equal(length(unique(defs)), 3L)
  expect_false(any(grepl("sensitivity", tolower(defs))))
})

test_that("Z8.5 the JAK2 limits land in separate columns, not one sensitivity", {
  l <- read_ext("jak2v617f_ddpcr_limits.csv")
  expect_setequal(l$quantity, c("LoB", "LoD", "LoQ"))
  # no merged column may exist in the source table
  expect_false("sensitivity" %in% names(l))

  # the reported value and definition are keyed by quantity, so a pivot cannot
  # silently overwrite one limit with another
  value <- stats::setNames(l$value_vaf_pct, l$quantity)
  definition <- stats::setNames(l$definition, l$quantity)
  wide <- tibble::as_tibble(as.list(value))
  expect_setequal(names(wide), c("LoB", "LoD", "LoQ"))
  expect_equal(length(unique(definition)), 3L)

  # LoB < LoD <= LoQ; a reversal means the table was misread
  v <- value[c("LoB", "LoD", "LoQ")]
  expect_equal(v, sort(v))
})

test_that("Z8.5 the artefact decision names the quantities per scheme", {
  genome <- Biostrings::DNAStringSet(c(Chr1 = strrep("ACGT", 40L)))
  plan <- marker_plan(genome) +
    marker_targets(chrom = "Chr1", pos = 20L, ref = "T", alt = "G") +
    marker_tech(schemes = "dpcr")
  tg <- plan$targets
  solved <- structure(
    list(
      report = tibble::tibble(
        label = tg$label, chrom = tg$chrom, pos = tg$pos, ref = tg$ref,
        alt = tg$alt, type = tg$type, scheme = NA_character_,
        cost = NA_real_, ok = FALSE, tried = "dpcr", marker = list(NULL)
      ),
      attempts = tibble::tibble(
        label = tg$label, scheme = "dpcr", cost = 1, ok = FALSE
      ),
      plan = plan, schemes = "dpcr",
      cache = new.env(parent = emptyenv()), knobs = list()
    ),
    class = "marker_solved"
  )
  dec <- artefact_decisions(solved)
  expect_true("analytical" %in% names(dec))
  expect_match(dec$analytical, "LoB")
  expect_match(dec$analytical, "LoD")
  expect_match(dec$analytical, "LoQ")
  expect_false(grepl("sensitivity", dec$analytical))

  # a scheme that reports no analytical limits says so with NA, not with a
  # merged column
  expect_null(scheme_analytical("snp"))
})

test_that("Z7.6 a deliberate mismatch records the late-cycle readout limit", {
  with_mm <- list(inner = list(ref = list(stab_at = 3L), alt = list(stab_at = NA)))
  plain <- list(inner = list(ref = list(stab_at = NA_integer_), alt = list(stab_at = NA)))
  solved_with <- function(m) list(report = tibble::tibble(marker = list(m)))

  lim <- tidymarkers:::artefact_limitations(solved_with(with_mm))
  expect_length(lim, 1L)
  expect_match(lim, "late-cycle discrimination", fixed = TRUE)
  # the statement names the mechanism, so it is a limitation and not a token
  expect_match(lim, "amplicon", fixed = TRUE)

  # a design without the mismatch claims nothing, and an empty marker is safe
  expect_length(tidymarkers:::artefact_limitations(solved_with(plain)), 0L)
  expect_length(tidymarkers:::artefact_limitations(solved_with(NULL)), 0L)

  # every assembled marker carries its decision tree under `$decision`, a
  # tbl_graph; `[[` on a graph is a vertex lookup and aborts, so the walk must
  # stop at the graph instead of descending into it
  skip_if_not_installed("tidygraph")
  g <- tidygraph::tbl_graph(edges = tibble::tibble(from = 1L, to = 2L))
  marker <- list(inner = list(ref = list(stab_at = 3L)), decision = g)
  expect_length(tidymarkers:::artefact_limitations(solved_with(marker)), 1L)
})

test_that("Z8.5 a non-dPCR scheme declares only the quantity it reports", {
  aq <- scheme_analytical("sanger")
  expect_equal(names(aq), "MAF")

  txt <- marker_decision_text("sanger")
  h <- which(grepl("^analytical performance", txt))
  expect_length(h, 1L)
  block <- txt[(h + 1L):(h + length(aq))]
  expect_equal(trimws(sub("^\\s*(\\S+).*$", "\\1", block)), "MAF")
  expect_false(any(grepl("sensitivity", tolower(block))))
  # the scheme's declared floor is the number the quantity refers to
  expect_equal(SANGER_DETECTION$floor, 0.15)
  expect_false(SANGER_DETECTION$calibrated)
})
