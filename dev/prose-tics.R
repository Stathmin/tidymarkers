# =============================================================================
#  dev/prose-tics.R  --  report LLM wordiness in docs and code comments
# =============================================================================
#  A report, not a gate. `--check` fails only on the "filler" tier: phrases that
#  carry no information and are a known LLM tell. The "contrastive" tier
#  ("X, not Y"; "rather than"; "not just") is this repository's deliberate
#  teaching style -- hundreds of hits -- so it is printed on request and never
#  gated. Dev-only: `dev/` is in .Rbuildignore.
#
#  Scans: DESCRIPTION, README.md, man/*.Rd, vignettes/*.Rmd and the roxygen
#  (`#'`) lines of R/*.R. Code lines are not prose and are skipped.
#
#  Run:  Rscript dev/prose-tics.R [--check] [--contrastive]
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) {
  dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
} else {
  getwd()
}
setwd(root)

args <- commandArgs(trailingOnly = TRUE)
check <- "--check" %in% args
show_contrastive <- "--contrastive" %in% args

# Two tiers: `filler` is gated by --check, `contrastive` is reported only.
tics <- list(
  filler = c(
    "it(')?s worth noting", "it is (important|worth) (to note|noting)",
    "\\bdelve", "\\btapestry", "\\bseamless", "\\bpivotal", "\\bcrucial",
    "\\bfurthermore", "\\bmoreover", "\\bthat said", "\\bleverag(e|ing)\\b"
  ),
  contrastive = c(
    "\\b(is|are|was|were)\\b[^.]{0,90},[[:space:]]+not\\b",
    "\\bnot\\b[^.]{0,60}\\bbut\\b", "\\brather than\\b",
    "\\bnot just\\b", "\\bnot only\\b", "\\badditionally\\b"
  )
)

prose_lines <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (endsWith(path, ".R")) {
    keep <- grepl("^\\s*#'", lines)
    return(list(lines[keep], which(keep)))
  }
  list(lines, seq_along(lines))
}

files <- c(
  if (file.exists("DESCRIPTION")) "DESCRIPTION",
  if (file.exists("README.md")) "README.md",
  list.files("man", pattern = "[.]Rd$", full.names = TRUE),
  list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE),
  list.files("R", pattern = "[.]R$", full.names = TRUE)
)

scan <- function(patterns, tier) {
  hits <- list()
  for (f in files) {
    p <- prose_lines(f)
    text <- p[[1L]]
    lines <- p[[2L]]
    for (pat in patterns) {
      m <- grepl(pat, text, ignore.case = TRUE, perl = TRUE)
      for (i in which(m)) {
        hits[[length(hits) + 1L]] <- data.frame(
          file = f, line = lines[[i]], tier = tier,
          text = trimws(substr(text[[i]], 1L, 160L)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(hits)) do.call(rbind, hits) else NULL
}

report <- function(hits, name) {
  if (is.null(hits)) {
    message("  ", name, ": 0")
    return(invisible())
  }
  message(
    "  ", name, ": ", nrow(hits), " (", length(unique(hits$file)), " files)"
  )
  for (i in seq_len(nrow(hits))) {
    message("    ", hits$file[i], ":", hits$line[i], "  ", hits$text[i])
  }
}

filler <- scan(tics$filler, "filler")
message("prose tics")
report(filler, "filler (gated by --check)")
if (show_contrastive) {
  report(scan(tics$contrastive, "contrastive"), "contrastive (never gated)")
} else {
  message("  contrastive: not scanned (pass --contrastive to list)")
}

if (check && !is.null(filler)) {
  stop(nrow(filler), " filler phrase(s) found", call. = FALSE)
}
