## Regenerate data-raw/proposed/crisprbase_registry.csv from the crisprBase
## source distribution's own data/*.rda files (see report-crispr.md section 3).
## The package is NOT installed here, so the S4 classes are re-declared minimally
## with the same slot names and the objects are loaded into that hierarchy.
suppressMessages({library(Biostrings); library(S4Vectors)})
Nuclease <- setClass("Nuclease", representation(nucleaseName = "character", targetType = "character",
  motifs = "ANY", weights = "numeric", metadata = "list"))
CrisprNuclease <- setClass("CrisprNuclease", contains = "Nuclease",
  representation(pam_side = "character", spacer_gap = "integer", spacer_length = "integer"))

d <- ".tmp/probe/crisprbase/crisprBase-devel/data"
if (!dir.exists(d)) stop("missing ", d, " -- re-fetch the crisprBase tarball (see report)")
stopifnot(file.exists(".tmp/probe/crisprbase/crisprBase-devel/DESCRIPTION"))
ver <- sub("^Version:[[:space:]]*", "",
  grep("^Version:", readLines(".tmp/probe/crisprbase/crisprBase-devel/DESCRIPTION"), value = TRUE))

keys <- c("AsCas12a", "CasRx", "Csm", "enAsCas12a", "MAD7", "SaCas9", "SpCas9", "SpG")
rows <- list()
for (n in keys) {
  f <- file.path(d, paste0(n, ".rda"))
  e <- new.env(); load(f, envir = e); o <- get(n, e)
  mo <- as.character(o@motifs)
  ## A PFS model exists only if the motif list constrains the flanking site.
  ## The two RNA-targeting nucleases ship the catch-all motif "N", which as an
  ## IUPAC code means any base: a well-formed object that restricts nothing.
  ## That is NOT a PFS model, and it is recorded as "none (catch-all N)".
  pfs <- if (!identical(o@targetType, "DNA")) {
    if (all(mo == "N")) "none (catch-all N)" else "present"
  } else "not applicable (DNA target)"
  ## No seed slot exists on the class at all; verified from the class definition.
  rows[[n]] <- data.frame(
    key = n,
    class = class(o)[1],
    target = o@targetType,
    motif_side = o@pam_side,
    spacer_len = o@spacer_length,
    n_motifs = length(mo),
    motif_letters = paste(sort(unique(unlist(strsplit(mo, "")))), collapse = ""),
    weights_are_binary = all(o@weights == 1),
    pfs_model = pfs,
    seed_model = "none (class has no seed slot)",
    metadata = paste(names(o@metadata), collapse = ";"),
    crisprbase_version = ver,
    rda_md5 = unname(tools::md5sum(f)),
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows); rownames(out) <- NULL
## the class definition must genuinely lack a seed/PFS slot, for itself AND for
## the Nuclease parent it contains. Parsed from the actual `slots = c(...)` block.
all_slots <- c()
for (cls_file in c("CrisprNuclease-class.R", "Nuclease-class.R")) {
  p <- file.path(".tmp/probe/crisprbase/crisprBase-devel/R", cls_file)
  if (!file.exists(p)) stop("missing class source: ", p)
  cls <- readLines(p)
  inside <- function(x) grepl("^[[:space:]]+[a-zA-Z_]+ = \"[a-zA-Z]+\"", x)
  all_slots <- c(all_slots, sub("^[[:space:]]+([a-zA-Z_]+) = .*", "\\1", Filter(inside, cls)))
  all_slots <- c(all_slots, sub("^#\\' @slot ([a-zA-Z_]+) .*", "\\1",
    grep("^#\\' @slot ", cls, value = TRUE)))
}
all_slots <- sort(unique(all_slots))
if (any(grepl("seed|pfs", all_slots, ignore.case = TRUE))) {
  stop("crisprBase now HAS a seed/PFS slot; the report's table is stale: ",
       paste(all_slots, collapse = ","))
}
out <- cbind(out, class_slots = paste(all_slots, collapse = ","))
stopifnot(nrow(out) == 8L, !"LwaCas13a" %in% out$key, "CasRx" %in% out$key)
dir.create("data-raw/proposed", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(out, "data-raw/proposed/crisprbase_registry.csv")
cat("wrote", nrow(out), "rows, crisprBase", ver, "\n")
