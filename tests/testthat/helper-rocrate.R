# Structural readers for an RO-Crate metadata document.  Shared by the crate
# tests: the JSON-LD is read as nested lists (`simplifyVector = FALSE`) because
# `@graph` entries do not share a shape.

rocrate_meta <- function(crate) {
  jsonlite::fromJSON(
    readr::read_file(fs::path(crate, "ro-crate-metadata.json")),
    simplifyVector = FALSE
  )
}

rocrate_entry <- function(meta, id) {
  for (entry in meta[["@graph"]]) {
    if (identical(as.character(entry[["@id"]]), id)) {
      return(entry)
    }
  }
  NULL
}
