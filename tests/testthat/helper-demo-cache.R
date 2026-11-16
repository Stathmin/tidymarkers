# The Atlas/TAIR10 demo data is fetched on first use by its exported getters
# (`ath_genome()`, `ensembl_variants()`, `ensembl_gene()`).  A test must never
# make that fetch: a cold cache means a network round-trip inside the suite,
# which is slow at best and an ERROR on a runner that cannot reach Ensembl or
# EBI.  Use the cache when it is warm; skip (never error) when it is cold.

# The three fixtures the suite reads.  Names are the cache file names the
# getters write, so a rename in `R/` fails these tests loudly instead of
# silently skipping them.
demo_cache_files <- c(
  genome = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
  variants = "arabidopsis_thaliana_var_4:1500000-1560000.json",
  gene = "arabidopsis_thaliana_AT4G00650.json"
)

#' Is the cached demo fixture present (so no download will happen)?
#' @param what any of `"genome"`, `"variants"`, `"gene"`.
#' @noRd
demo_cache_warm <- function(what = "genome") {
  stopifnot(all(what %in% names(demo_cache_files)))
  files <- demo_cache_files[what]
  # `fs::file_exists()` drops names; keep them so the skip message can name
  # the missing fixture instead of printing TRUE/FALSE.
  stats::setNames(fs::file_exists(fs::path(marker_cache_dir(), files)), what)
}

#' Skip when the demo fixture would have to be downloaded
#' @noRd
skip_if_cold_demo_cache <- function(what = "genome") {
  warm <- demo_cache_warm(what)
  cold <- names(warm)[!warm]
  if (length(cold)) {
    testthat::skip(paste0(
      "cold demo cache: ", paste(cold, collapse = ", "),
      " (warm it once with Rscript -e 'tidymarkers::ath_genome()')"
    ))
  }
}
