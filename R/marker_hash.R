# =============================================================================
#  marker_hash.R  --  one content-hash seam
# =============================================================================
#  `digest` stays the implementation; this is the only place that names the
#  algorithm.  A manifest, an index shard and an artefact therefore hash the
#  same way, and switching the digest is one edit, not a tree-wide grep.
#
#  A path and a value are different inputs: `hash_file()` reads the bytes,
#  `hash_object()` serialises the R object.  Guessing wrong hashes nothing.
# =============================================================================

#' SHA-256 of an R object
#'
#' `digest::digest()` serialises the object; a character value is data here,
#' never a path (use `hash_file()` for bytes).
#' @param x any R object to serialise and hash.
#' @return a length-one lowercase hex string, `NA` when hashing fails.
#' @keywords internal
#' @noRd
hash_object <- function(x) {
  tryCatch(
    as.character(digest::digest(x, algo = "sha256")),
    error = function(e) NA_character_
  )
}

#' SHA-256 of a file's bytes
#' @param path a file path.
#' @return a length-one lowercase hex string, `NA` when the file is missing or
#'   unreadable.
#' @keywords internal
#' @noRd
hash_file <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NA_character_)
  }
  tryCatch(
    as.character(digest::digest(file = path, algo = "sha256")),
    error = function(e) NA_character_
  )
}
