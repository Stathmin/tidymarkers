# =============================================================================
#  data-raw/mismatch-panel.R  --  build mismatch_panel.csv with BLAST as oracle
# =============================================================================
#  The point of this script is that the observation is made by an engine that is
#  NOT the one under test.  DECIPHER estimates a priming-efficiency margin; whether
#  BLAST can seed a full-length hit cannot be decided by DECIPHER, so the 0/1 flag
#  is genuinely external, and DECIPHER's margin is scored against it downstream.
#
#  What this panel is, precisely:
#
#    * The template is the exact reverse complement of the primer locus, so the
#      primer's 3' base pair with the template's first base.  Only then does a
#      lost 3' base mean "no full-length hit".
#    * blastn -task blastn-short seeds from the 5' end and stops extending at a
#      mismatch edge, so `qend == nchar(primer) && mismatch == 0` is the
#      full-length test.  There is no 1-base overhang artefact to correct for;
#      BLAST simply shortens the hit when the 3' base is lost.
#    * A perfect 21-mer is NOT reported as qend 21 in every context: blastn-short
#      left-extends a perfect hit by one base, so a perfect match on the template
#      comes back as qstart 0, qend 20.  The rule is therefore "the match covers
#      the whole primer AND extends past it", not "qend equals the length".
#    * Each primer is reported against its OWN template in a single-subject
#      search.  batching many subjects makes BLAST return the best-scoring one,
#      which silently confounds the base under test.
#
#  Run:  Rscript data-raw/mismatch-panel.R   (needs `blastn` on PATH)
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

blastn <- Sys.which("blastn")
if (!nzchar(blastn)) {
  stop("blastn not found on PATH; BLAST+ is required to rebuild mismatch_panel.csv")
}

# the 21-mer from the package's shared KASP fixture family (3' base A)
primer <- "CATGGACCTAGCATTGACGTA"
tgt <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(primer)))
stopifnot(substr(primer, 1L, 1L) ==
  as.character(Biostrings::complement(Biostrings::DNAString(
    substr(tgt, nchar(tgt), nchar(tgt))
  ))))

setb <- function(s, p, b) {
  s <- strsplit(s, "", fixed = TRUE)[[1L]]
  s[p] <- b
  paste(s, collapse = "")
}

wd <- c(
  AA = "TT", AC = "GT", AG = "CT", AT = "AT", CA = "TG", CC = "GG", CG = "CG",
  CT = "AG", GA = "TC", GC = "GC", GG = "CC", GT = "AC", TA = "TA", TC = "GA",
  TG = "CA", TT = "AA"
)

# one query against one subject, with the hit covering the whole primer
blast_full_length <- function(query, subject, n = nchar(query)) {
  tmp <- withr::local_tempdir()
  qf <- fs::path(tmp, "q.fa")
  sf <- fs::path(tmp, "s.fa")
  readr::write_lines(c(">q", query), qf)
  readr::write_lines(c(">s", subject), sf)
  out <- suppressWarnings(system2(blastn, c(
    "-task", "blastn-short", "-query", shQuote(qf), "-subject", shQuote(sf),
    "-word_size", "7", "-outfmt", "6 qstart qend length mismatch",
    "-evalue", "10"
  ), stdout = TRUE, stderr = FALSE))
  if (!length(out) || !nzchar(out[[1L]])) {
    return(list(hit = 0L, detail = NA_character_))
  }
  # best hit = the longest alignment, then the fewest mismatches
  f <- do.call(rbind, strsplit(out, "\t", fixed = TRUE))
  len <- as.integer(f[, 3L])
  mm <- as.integer(f[, 4L])
  k <- order(-len, mm)[[1L]]
  qs <- as.integer(f[k, 1L])
  qe <- as.integer(f[k, 2L])
  # covers the whole primer, either exactly or overhanging it by one base
  hit <- as.integer((qe - qs + 1L) >= n && mm[k] == 0L)
  list(hit = hit, detail = paste(qs, qe, len[k], mm[k], sep = "/"))
}

# every primer 3' base against every template 3' base, on the mirrored template
rows <- purrr::map_dfr(c("A", "C", "G", "T"), function(p3) {
  purrr::map_dfr(c("A", "C", "G", "T"), function(t3) {
    pr <- setb(primer, nchar(primer), p3)
    tg <- setb(tgt, 1L, as.character(Biostrings::complement(
      Biostrings::DNAString(t3)
    )))
    b <- blast_full_length(pr, tg)
    m <- pcr_efficiency(pr, tg, cores = 1L)
    tibble::tibble(
      primer_3p = p3,
      template_3p = t3,
      mismatch = p3 != t3,
      blast_full_length = b$hit,
      blast_detail = b$detail,
      decipher_efficiency = round(m, 4L),
      engine_blast = "blastn blastn-short w7",
      engine_decipher = "DECIPHER CalculateEfficiencyPCR"
    )
  })
})

stopifnot(
  # the perfect matches must be full-length hits, or the test is meaningless
  all(rows$blast_full_length[!rows$mismatch] == 1L)
)

out <- fs::path(root, "inst", "extdata", "mismatch_panel.csv")
readr::write_csv(rows, out)
cat(
  "wrote", out, ":", nrow(rows), "rows;",
  sum(rows$blast_full_length), "full-length primer hits;",
  sum(rows$decipher_efficiency[!rows$mismatch] > 0.5), "DECIPHER efficient matches\n"
)
