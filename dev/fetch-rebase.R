# =============================================================================
#  dev/fetch-rebase.R  --  regenerate the bundled enzyme catalogue from REBASE
# =============================================================================
#  The restriction-enzyme table (inst/extdata/restriction_enzymes.csv) is *data*
#  and must be sourced, not typed from memory.  This downloads the REBASE
#  "bionet" flat file (the authoritative catalogue; NE Bionet format, one line
#  per recognition sequence with `^` marking the top-strand cut) and rebuilds the
#  CSV deterministically.
#
#  Kept enzymes: Type II recognition sequences that are (a) single-line,
#  (b) palindromic (site == reverse-complement, so a plus-strand match is
#  unambiguous and the cut is symmetric - our caps_fragments() model), (c) IUPAC
#  DNA of length 4..12, and (d) annotated with a cut.  Methyltransferases (names
#  or aliases starting with `M.`) are dropped.
#
#  Run:  Rscript dev/fetch-rebase.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

url <- "http://rebase.neb.com/rebase/link_bionet"
tmp <- withr::local_tempfile()
utils::download.file(url, tmp, method = "libcurl", quiet = TRUE)
lines <- readr::read_lines(tmp)
# the version line looks like:  REBASE version 610        bionet.610
version <- sub("^.*bionet\\.([0-9]+).*$", "\\1", grep("bionet\\.[0-9]+", lines, value = TRUE)[1L])

# entries: NAME [(alias)] SITE-with-carets
rx <- "^([^ ]+)\\s+(\\([^)]*\\)\\s+)?([ACGTRYSWKMBDHVN]+\\^?[ACGTRYSWKMBDHVN]*|[\\^]?[ACGTRYSWKMBDHVN]*\\^?[ACGTRYSWKMBDHVN]*) *$"
hit <- regmatches(lines, regexec(rx, lines))
rows <- purrr::map_dfr(hit, function(g) {
  if (length(g) < 4L) {
    return(tibble::tibble())
  }
  name <- g[[2L]]
  alias <- g[[3L]]
  site <- g[[4L]]
  tibble::tibble(name = name, alias = trimws(alias), site = site)
})

clean <- rows |>
  dplyr::filter(
    !grepl("^M\\.", .data$name), !grepl("\\(M\\.", .data$alias),
    stringr::str_count(.data$site, "\\^") == 1L
  ) |>
  dplyr::mutate(
    cut = stringr::str_locate(.data$site, "\\^")[, 1L] - 1L,
    site = stringr::str_remove(.data$site, "\\^")
  ) |>
  dplyr::filter(
    nchar(.data$site) >= 4L, nchar(.data$site) <= 12L,
    grepl("^[ACGTRYSWKMBDHVN]+$", .data$site)
  ) |>
  # a plus-strand match is only unambiguous for palindromic sites
  dplyr::mutate(palindrome = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(.data$site))) == .data$site) |>
  dplyr::filter(.data$palindrome) |>
  dplyr::distinct(.data$name, .keep_all = TRUE) |>
  dplyr::arrange(.data$name) |>
  dplyr::transmute(enzyme = .data$name, site = .data$site, cut = .data$cut)

out <- fs::path(root, "inst", "extdata", "restriction_enzymes.csv")
readr::write_csv(clean, out)
cat("REBASE bionet", version, "->", nrow(clean), "palindromic type-II enzymes\n")
cat("wrote", out, "\n")
