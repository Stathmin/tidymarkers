# A stable key for the per-file whole-design caches below.  It includes the
# genome, the target and every argument that changes the design, so a cached
# result can never be reused for a different fixture.
design_key <- function(genome, snp, args) {
  nm <- names(args)
  if (is.null(nm)) {
    nm <- rep("", length(args))
  }
  vals <- vapply(args, function(x) paste(x, collapse = ","), character(1))
  paste(
    paste(as.character(genome), collapse = ""),
    paste(unlist(snp), collapse = ","),
    paste(nm, vals, sep = "=", collapse = "|"),
    sep = "|"
  )
}
