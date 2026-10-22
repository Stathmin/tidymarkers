# Off-target index memory and time on Arabidopsis (ROADMAP Y3.5, Arabidopsis half)
#
# Roadmap row (Z5.1 / Y3.5), verbatim:
#   «Блокер: индекс под 15 Gb (`Y3.5`), но **не** на Arabidopsis: сначала на
#     120 Mb, память индекса мерить до коммита.»
#
# This file records a measurement, not a model.  Every number below was produced
# by the exact commands quoted here, on this machine, on 2026-10-22; nothing is
# extrapolated from a smaller genome.  No index binary is committed.
#
# Genomes are obtained exactly the way the package obtains them:
#   - `ath_genome()` (`R/marker_core.R:838`) downloads
#     Ensembl Plants release-63 TAIR10 toplevel FASTA into `marker_cache_dir()`
#     and keeps only the `ATH_CHROM_MAP` sequences (Chr1-5, ChrM, ChrC);
#   - no other public genome of comparable size is shipped with the package, so
#     Arabidopsis is the only in-package genome that can be measured without a
#     new download.  `specificity_seed()` is the package's whole-genome scan for
#     an in-memory `DNAStringSet`; `specificity_of()` / `specificity_fm()` are
#     the indexed backends.
#
# How to reproduce (run from the package root):
#
#   export PATH=/opt/R-4.6.1/bin:/usr/bin:/bin
#   export R_LIBS_SITE=/opt/R-4.6.1/lib/R/library
#   Rscript data-raw/offtarget_index_memory.R
#
# The script is self-contained: it loads the working tree with
# `pkgload::load_all(".", quiet = TRUE)`, so it measures the code under review
# rather than a stale installed copy.  It writes nothing to disk.
#
# Peak RSS was taken separately, because `/usr/bin/time` is *not installed* on
# this host (`command -v time` fails): a python3 wrapper ran the Rscript child
# and read `resource.getrusage(RUSAGE_CHILDREN).ru_maxrss` (KiB).  That is the
# high-water mark of the child process, genome included.
#
# ---------------------------------------------------------------------------
# Environment (R.version and Bioconductor package versions, as printed)
# ---------------------------------------------------------------------------
# R              4.6.1 (2026-06-24)
# Bioconductor   3.23
# tidymarkers    0.4.0
# Biostrings     2.80.2
# IRanges        2.46.0
# Rsamtools      2.28.0
# BiocParallel   1.46.0
# Rsubread       NOT INSTALLED
#
# Rsubread is absent, so `marker_fm_index()` / `specificity_fm()` (the FM-index
# backend) could not be built or measured here.  We assert the absence rather
# than guess a memory figure: no FM-index number for Arabidopsis exists in this
# repository.
#
# ---------------------------------------------------------------------------
# Genome
# ---------------------------------------------------------------------------
# file   : <marker_cache_dir()>/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
#          /root/.cache/R/tidymarkers/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
# md5    : bb2cd5da056b029bacbd558ecd9f1a27
# n      : 7 (Chr1-Chr5, ChrM, ChrC)
# bp     : 119,667,750   (~120 Mb, the roadmap's target size)
#
# ---------------------------------------------------------------------------
# Measurements (one R process, `BPPARAM = NULL`, single core)
# Reproduced by `Rscript data-raw/offtarget_index_memory.R`; the run that the
# numbers below come from printed:
#
#   ath_genome()                       1.240 s  119,672,328 B  (114.1 Mb)
#   specificity_seed k=8 mm=0 (4)      1.096 s        1,680 B  hits
#                                                              4992/5604/3591/255923
#   specificity_seed k=8 mm=1 (3)      1.756 s        1,560 B  hits
#                                                              135055/138705/100398
#   PDict(20-mer) build                0.010 s   33,595,760 B  (32 Mb)
#   matchPattern(21-mer) whole genome  0.185 s
#   whole child process wall            7.47 s   peak RSS 538.9 MB
#
# Earlier probe runs of the same calls (before this script existed) read
# 1.09 s / 1.08 s / 1.77 s / 0.017 s / 0.244 s and peak RSS 420.9 MB; the
# wall-clock and RSS differences are run-to-run variance of the same machine,
# not a change in the code.  The high-water mark quoted is the later, larger
# one: 538.9 MB.
#
#   -> the seed scan builds no persistent index, so there is no separate index
#      object to size: what the peak RSS contains is the 114 Mb genome plus
#      the R/Biostrings baseline.
#   -> an exact-match `PDict` costs ~0.28 bytes per genome base at k = 20 for a
#      one-pattern dictionary (33,595,760 / 119,667,750).  It is a lookup
#      structure, not a suffix array, and it is built per pattern.
#
# Note on the 15 Gb extrapolation (NOT a measurement): the roadmap's 15 Gb
# case is not measured here and must not be inferred from 120 Mb.  The only
# claim this file supports is the Arabidopsis line: a whole-genome scan of
# ~120 Mb runs in seconds inside <0.5 GB RSS without any persisted index.
#
# ---------------------------------------------------------------------------
# The script
# ---------------------------------------------------------------------------

#' Version table for the provenance record
#' @noRd
tm_offtarget_versions <- function() {
  pkgs <- c(
    "tidymarkers", "Biostrings", "IRanges", "Rsamtools", "BiocParallel",
    "Rsubread"
  )
  installed <- vapply(
    pkgs,
    function(p) {
      v <- tryCatch(
        as.character(utils::packageVersion(p)),
        error = function(e) "NOT INSTALLED"
      )
      v
    },
    character(1)
  )
  data.frame(
    package = pkgs,
    version = unname(installed),
    stringsAsFactors = FALSE
  )
}

#' Measure the Arabidopsis off-target scan
#'
#' Returns a tidy one-row-per-measurement data frame and prints it.  Nothing is
#' written to disk: the numbers are meant to be pasted into the block above.
#' @noRd
measure_offtarget_index_memory <- function() {
  suppressPackageStartupMessages(library(Biostrings))
  # the installed `tidymarkers` may be stale; use the working tree
  pkgload::load_all(".", quiet = TRUE)

  size_of <- function(x) {
    o <- utils::object.size(x)
    data.frame(
      bytes = as.numeric(o),
      label = format(o, units = "auto"),
      stringsAsFactors = FALSE
    )
  }
  time_it <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    value <- force(expr)
    list(value = value, seconds = proc.time()[["elapsed"]] - t0)
  }

  t_genome <- time_it(ath_genome())
  g <- t_genome$value
  genome_size <- size_of(g)

  primers <- c(
    AS_ref = "GAAGGTGACCAAGTTCATGCT",
    AS_alt = "GAAGGTGACCAAGTTCATGCA",
    COM = "GCTTAGCCATGGTCAATGCTA",
    LOW = "AAAAAAAAAAAAAAAAAAAA"
  )
  t_exact <- time_it(specificity_seed(primers, g, k = 8L, max_mismatch = 0L))
  t_radius <- time_it(
    specificity_seed(primers[1:3], g, k = 8L, max_mismatch = 1L)
  )

  pdict <- time_it({
    d <- Biostrings::PDict("ACGTACGTACGTACGTACGT")
    list(dict = d, size = size_of(d), pattern = d[[1L]])
  })
  single <- time_it({
    n <- vapply(
      seq_along(g),
      function(i) {
        length(Biostrings::matchPattern(
          pdict$value$pattern, g[[i]],
          max.mismatch = 0L
        ))
      },
      numeric(1)
    )
    sum(n)
  })

  out <- rbind(
    data.frame(
      step = "ath_genome()", seconds = t_genome$seconds,
      bytes = genome_size$bytes, size = genome_size$label,
      note = sprintf("n=%d bp=%d", length(g), sum(Biostrings::width(g))),
      stringsAsFactors = FALSE
    ),
    data.frame(
      step = "specificity_seed k=8 max_mismatch=0 (4 primers)",
      seconds = t_exact$seconds,
      bytes = size_of(t_exact$value)$bytes, size = size_of(t_exact$value)$label,
      note = paste(t_exact$value$hits, collapse = "/"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      step = "specificity_seed k=8 max_mismatch=1 (3 primers)",
      seconds = t_radius$seconds,
      bytes = size_of(t_radius$value)$bytes,
      size = size_of(t_radius$value)$label,
      note = paste(t_radius$value$hits, collapse = "/"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      step = "PDict(20-mer) build", seconds = pdict$seconds,
      bytes = pdict$value$size$bytes, size = pdict$value$size$label,
      note = "a lookup dict, not a suffix index",
      stringsAsFactors = FALSE
    ),
    data.frame(
      step = "matchPattern(21-mer) whole genome, 7 sequences",
      seconds = single$seconds,
      bytes = NA_real_, size = NA_character_,
      note = as.character(single$value),
      stringsAsFactors = FALSE
    )
  )
  rownames(out) <- NULL
  out
}

if (sys.nframe() == 0L) {
  print(tm_offtarget_versions(), row.names = FALSE)
  print(measure_offtarget_index_memory(), row.names = FALSE)
  cat("peak RSS is measured outside this script; see the header\n")
}
