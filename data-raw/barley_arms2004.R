# =============================================================================
#  data-raw/barley_arms2004.R  --  BLOCKED source, recorded as such
# =============================================================================
#  Hayden MJ, Tabone TL, Mather DE. "Development and assessment of simple PCR
#  markers for SNP loci in barley."  Genome 2004;47:414-420.
#  doi:10.1139/g03-130
#
#  WHAT THIS SCRIPT IS FOR
#  The dispatch asked for "Genotyping single nucleotide polymorphisms in barley
#  by tetra-primer ARMS-PCR" at this DOI. The DOI actually resolves to the title
#  above (verified through Europe PMC EXT_ID 15060595 and Crossref); there is no
#  separate tetra-primer-ARMS barley paper at g03-130.
#
#  WHY NOTHING IS EXTRACTED
#  Europe PMC: no PMCID, isOpenAccess = false. Unpaywall: is_oa = false,
#  oa_status = "closed". The publisher (cdnsciencepub.com) answers every PDF and
#  full-text URL with HTTP 403 and a Cloudflare interstitial ("Just a moment...").
#  Under the worker rules a closed paper is recorded as blocked; abstracts and
#  vendor pages may never become outcomes. So there is NO numbers file for this
#  key, and this script's only job is to make that absence explicit and testable.
#
#  THE ASSERTION
#  The script re-checks the two facts that justify "blocked": that no full text
#  was downloaded into .tmp/lit/barley_arms2004, and that the DOI's metadata is
#  closed access. If someone later fetches the paper, the script stops and says
#  so - at which point a real extract should replace it. A blocked record must
#  not quietly outlive the block.
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!exists("root")) {
  root <- if (length(script)) {
    dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]]))))
  } else {
    getwd()
  }
  setwd(root)
}

DOI <- "10.1139/g03-130"
lit <- file.path(root, ".tmp", "lit", "barley_arms2004")

# A real fetch would be a full text plus a supplement. An empty or absent dir is
# the expected state; anything else means the block has lifted.
if (dir.exists(lit)) {
  files <- list.files(lit, recursive = TRUE)
  substantial <- files[file.size(file.path(lit, files)) > 20000]
  if (length(substantial)) {
    stop("barley_arms2004 now looks fetched (", paste(substantial, collapse = ", "),
         "). Re-check access and replace this blocked record with a real extract ",
         "rather than leaving it in place.")
  }
}

# The DOI must still be the paper we think it is. If the title changes, the DOI
# was wrong, and the reason for the block ("no tetra-primer ARMS barley paper at
# this DOI") would itself be wrong.
meta <- tryCatch({
  url <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=",
                "DOI:%2210.1139/g03-130%22&format=json&resultType=core")
  jsonlite::fromJSON(url, simplifyVector = FALSE)$resultList$result[[1L]]
}, error = function(e) NULL)

if (is.null(meta)) meta <- list()
or_null <- function(x) if (is.null(x) || !length(x)) "" else x
if (length(meta)) {
  title <- or_null(meta$title)
  if (!grepl("barley", title, ignore.case = TRUE)) {
    stop("DOI ", DOI, " now resolves to \"", title,
         "\" which is not the barley paper; the blocked record is stale")
  }
  journal <- or_null(meta$journalInfo$journal$title)
  year <- or_null(meta$pubYear)
  access <- if (isTRUE(meta$isOpenAccess)) "open" else "closed"
  if (access != "closed") {
    stop("DOI ", DOI, " is now open access; extract the paper instead")
  }
  message("blocked: ", DOI, " (", journal, " ", year, ") - ", title)
  message("  access: ", access, "; no PMCID. No numbers landed.")
}

# The row that would be proposed to the registry, with no outcomes attached.
blocked <- data.frame(
  key = "barley_arms2004",
  doi = DOI,
  citation = "Hayden MJ, Tabone TL, Mather DE. Development and assessment of simple PCR markers for SNP loci in barley",
  venue = "Genome",
  year = 2004L,
  status = "blocked",
  reason = paste(
    "closed access: no PMCID, isOpenAccess=false, Unpaywall oa_status=closed;",
    "publisher CDN Science returns HTTP 403 + Cloudflare challenge for PDF and",
    "full-text URLs. No paywall bypass attempted, so no per-SNP outcomes and no",
    "deliberate-mismatch positions were extracted."
  ),
  what_it_would_anchor = "barley SNP PCR marker design and success rates (not extracted)",
  outcome_definition = "NONE - blocked source, no numbers landed",
  evidence_tier = "blocked",
  source_file = NA_character_,
  stringsAsFactors = FALSE
)

outfile <- file.path(root, "inst", "extdata", "arms_barley_arms2004_blocked.csv")
utils::write.csv(blocked, outfile, row.names = FALSE, na = "")
message("wrote ", outfile, " (status = blocked, 0 numbers)")
