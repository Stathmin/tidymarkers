# =============================================================================
#  data-raw/polymarker-source.R  --  PolyMarker's primer-placement rules,
#  transcribed from the source with the enclosing method as the locator
# =============================================================================
#  PolyMarker (Ramirez-Gonzalez, Uauy & Caccamo; Bioinformatics 2015,
#  10.1093/bioinformatics/btv069) is the published wheat KASP design pipeline.
#  The paper is closed access, so the *source* is the evidence here: every row
#  carries the method (and, where useful, the line) it was read from.  Treat it
#  as a design RULE SOURCE and a benchmark, not as something to reimplement.
#
#  Fetch:
#    git clone --depth 1 \
#      https://github.com/Informatics-John-Innes-Centre/bio-polymarker.git \
#      .tmp/lit/polymarker/repo
#  Read at VERSION 1.3.3, commit d4a4578 (2025-07-25).  The script re-checks the
#  file digests, because a probe that moves upstream silently turns signed-off
#  rules into guesses.  If a digest fails, the assertion below names the file and
#  the fix is to re-read the method, not to update the hash.
#
#  Run:  Rscript data-raw/polymarker-source.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

repo <- Sys.getenv("POLYMARKER_REPO", file.path(root, ".tmp", "lit", "polymarker", "repo"))
if (!dir.exists(repo)) {
  stop("PolyMarker source not found at ", repo, "; see the header for the clone command")
}

read <- function(rel) readLines(file.path(repo, rel), warn = FALSE)
# NB: these lists are combined with c(), which DROPS duplicate names.  The
# source list uses `lines`, the row list uses `quoted`, so no name collides.
prov <- function(rel) {
  list(
    lines = read(rel), file = rel,
    digest = unname(tools::md5sum(file.path(repo, rel)))
  )
}

snp_rb <- prov("lib/bio/PolyploidTools/SNP.rb")
engine_rb <- prov("lib/bio/db/primer3.rb")
tags <- c(
  "lib/bio/PolyploidTools/SNP.rb" = "5fd84bf5f6c0b89a69fee7ebfdbe16ff",
  "lib/bio/PolyploidTools/ExonContainer.rb" = "630eb8a5b77e2b920b4a45e087f00f59",
  "lib/bio/db/primer3.rb" = "8abcc3ddff8330b1f8b9f51226776561"
)
for (nm in names(tags)) {
  got <- unname(tools::md5sum(file.path(repo, nm)))
  if (!identical(got, unname(tags[[nm]]))) {
    stop(
      nm, " changed (md5 ", got, ", read at ", tags[[nm]], ").\n",
      "  The rules below were read from a specific revision; re-read the method",
      "  and re-sign the row, do not just update the digest."
    )
  }
}

has <- function(p, pattern, method) {
  n <- grep(pattern, p$lines, fixed = TRUE)
  if (length(n) == 0L) {
    stop("no line matches '", pattern, "' in ", p$file)
  }
  list(file = p$file, method = method, line = n[[1L]],
       quoted = p$lines[[n[[1L]]]])
}

# --- the rows ---------------------------------------------------------------
# role/scope: which part of the ASSAY construct a rule governs.
#   allele-primer  -> the two allele-specific primers
#   common-primer  -> the shared primer
#   construct      -> which primers exist at all
rules <- list()

# 1. The allele-specific primers are position-locked onto the SNP.  The SNP's
#    index is the 3' end forced on the template, and the allele base is written
#    in before the design runs.
rules$as_pinned <- c(has(snp_rb, "seq_original[pr.snp_pos] = self.original",
                         "primer_3_all_strings"),
                     list(role = "allele-primer", scope = "3'-terminal base = the SNP",
                          rule = "the SNP base index is forced to the primer's 3' end, and the template is designed once with the reference base and once with the alternate base; the pair IS the allele-specific pair"))

# 2. pr.snp_pos is the SNP's index in the primer region, and both the forward
#    and the reverse-complement template are derived from it.
rules$snp_index <- c(has(snp_rb, "when mask[i] == '&'",
                         "primer_region"),
                     list(role = "allele-primer", scope = "SNP position",
                          rule = "the SNP's index in the primer region is recorded first (mask '&', or ':' when the SNP is itself homoeologous); every allele primer is built around that index"))

# 3. The common primer is where subgenome discrimination can come from: it is
#    forced onto a position the mask marks chromosome-specific.
rules$common_on_mask <- c(
  has(snp_rb, ":right_pos => pos", "primer_3_all_strings"),
  list(role = "common-primer",
       scope = "subgenome",
       rule = "for each mask position marked chromosome-specific, the common primer is forced to end (SEQUENCE_FORCE_RIGHT_END) on that position")
)
# The four forced-common classes are the two specificity levels x exon/intron;
# this row is the count, so a class silently dropped upstream fails here.
rules$common_on_mask_classes <- c(
  has(snp_rb, 'chromosome_nonspecific all', "primer_3_all_strings"),
  list(role = "common-primer", scope = "fallback",
       rule = "with no chromosome-specific position the common primer is designed with no force (chromosome_nonspecific) -- i.e. the assay is not subgenome-specific")
)
rules$common_avoid_homoeo <- c(
  has(snp_rb, "return nil if check_str != check_str.downcase", "return_primer_3_string"),
  list(role = "common-primer", scope = "guard",
       rule = "optionally refuse a common primer whose downstream window contains a homoeologous difference (README: 'not tested')")
)

# 4. What "chromosome specific" MEANS: a mask base is upper-cased only when the
#    other subgenomes differ from the target at that base and all genomes are
#    represented; lower case is the weaker 'almost' class.
rules$mask_definition <- c(has(snp_rb, 'masked_snps[i] = masked_snps[i].upcase if different == expected_snps and from_group == genomes_count',
                               "mask_aligned_chromosomal_snp"),
                           list(role = "common-primer", scope = "subgenome",
                                rule = "a position is chromosome-specific when every other aligned subgenome differs from the target AND all genomes are represented at it"))
rules$mask_almost <- c(
  has(snp_rb, "almost_chromosome_specific << position_in_region", "primer_region"),
  list(role = "common-primer", scope = "subgenome",
       rule = "lower-case mask positions (not every subgenome differs) are still offered as 'almost' specific common primers, and are scored lower")
)
rules$scores <- c(has(engine_rb, "@scores[:chromosome_semispecific] = 100", "initialize"),
                  list(role = "selection", scope = "ranking",
                       rule = "final pick ranks chromosome_specific (1000) > semi-specific (100) > non-specific (0), with exon designs above intron (50)"))
rules$homoeo_label <- c(has(snp_rb, '@snp_type = "homoeologous"', "primer_3_all_strings"),
                        list(role = "label", scope = "bookkeeping",
                             rule = "the SNP itself is recorded as homoeologous or non-homoeologous (a base differs between subgenomes at the marker); this is a LABEL on the SNP, separate from which primer carries subgenome specificity"))
rules$no_homoeo_axis <- c(has(engine_rb, "def homoeologous?", "Primer3Record#parse_header"),
                          list(role = "construct", scope = "explicit",
                               rule = "the construct has no 'which primer carries the subgenome' field: subgenome specificity is implied by the common primer being forced, never stated as a rule"))

# A rule-source row is only worth keeping if the quoted line exists verbatim.
for (nm in names(rules)) {
  r <- rules[[nm]]
  if (!is.character(r$quoted) || length(r$quoted) != 1L || is.na(r$quoted)) {
    stop("row '", nm, "' has no quoted line")
  }
  if (!any(grepl(r$quoted, read(r$file), fixed = TRUE))) {
    stop("row '", nm, "' quotes a line that is not in ", r$file)
  }
}
stopifnot(all(vapply(rules, function(r) r$role, character(1)) %in%
                c("allele-primer", "common-primer", "selection", "label", "construct")))

out <- do.call(rbind, lapply(names(rules), function(nm) {
  r <- rules[[nm]]
  data.frame(
    rule_id = nm,
    role = r$role,
    scope = r$scope,
    rule = r$rule,
    source_file = r$file,
    source_method = r$method,
    source_line = if (!is.null(r$line)) r$line else NA_integer_,
    quoted = r$quoted,
    relation_to_y3_2 = "subgenome discrimination is a COMMON-primer axis, not an AS-primer axis",
    source_type = "software source",
    evidence = "primary literature",
    calibrated = FALSE,
    doi = "10.1093/bioinformatics/btv069",
    version = "1.3.3",
    commit = "d4a4578",
    stringsAsFactors = FALSE
  )
}))
out <- out[order(out$rule_id), ]
rownames(out) <- NULL

## --- assertions on the SOURCE, not on our own output -------------------------
# The two claims Y3.2 hinges on must both be present: allele primers are forced
# onto the SNP, and the common primer is the one forced onto the subgenome-
# diagnostic base.  Count them, so a refactor that drops one fails here.
snp_txt <- read("lib/bio/PolyploidTools/SNP.rb")
stopifnot(
  nrow(out) == length(rules),
  sum(grepl("SEQUENCE_FORCE_LEFT_END", snp_txt, fixed = TRUE)) == 2L,
  sum(grepl("SEQUENCE_FORCE_RIGHT_END", snp_txt, fixed = TRUE)) == 1L,
  # the force flags are used on the allele template, and on the common template
  # only for the chromosome_specific classes (four classes: two specificity
  # levels x exon/intron)
  sum(grepl(":right_pos => pos", snp_txt, fixed = TRUE)) == 4L,
  grepl("chromosome_nonspecific all", snp_txt, fixed = TRUE) |> any(),
  # the common primer is never given a right_pos in the non-specific fallback
  !any(grepl("chromosome_nonspecific.*right_pos", snp_txt))
)

# PrimerRegion must declare the fields the two axes rest on.
pr <- paste(read("lib/bio/PolyploidTools/PrimerRegion.rb"), collapse = " ")
for (field in c("snp_pos", "chromosome_specific", "almost_chromosome_specific",
                "homoeologous")) {
  if (!grepl(field, pr, fixed = TRUE)) stop("PrimerRegion lacks field ", field)
}

## --- cross-check against the two published KASP rule sets we hold ----------
# Grewal 2020 (1000 wheat assays): the allele is in the AS 3' base even in the
# reverse-strand mode, so it is a separate axis from whatever the common primer
# does.  Makhoul: the common 3' base VARIES between assays for the same SNP
# class, which is what a subgenome anchor looks like (and the opposite of an
# allele axis).
g <- utils::read.csv(file.path(root, "inst", "extdata", "grewal_kasp_axes.csv"),
                     stringsAsFactors = FALSE)
m <- utils::read.csv(file.path(root, "inst", "extdata", "makhoul_kasp_primers.csv"),
                     stringsAsFactors = FALSE)

rc1 <- function(x) chartr("ACGT", "TGCA", x)
reverse_mode <- g$assay_mode == "reverse"
two_allele <- g$assay_mode == "two-allele"
stopifnot(
  # allele axis: AS 3' bases differ exactly when the design is two-allele
  all(g$as_fam_3p[two_allele] != g$as_hex_3p[two_allele]),
  all(g$as_fam_3p[reverse_mode] == unname(rc1(g$snp_allele_1[reverse_mode]))),
  all(g$as_hex_3p[reverse_mode] == unname(rc1(g$snp_allele_2[reverse_mode]))),
  # genome axis: in Makhoul the common 3' base is not tied to the SNP class
  length(unique(m$common_3p[m$snp_kasp == "[T/C]"])) > 1L,
  # and it is independent of the allele pair
  !identical(m$common_3p, m$fam_3p)
)

## --- write -------------------------------------------------------------------
out_dir <- fs::path(root, "inst", "extdata")
utils::write.csv(out, fs::path(out_dir, "polymarker_design_rules.csv"),
                 row.names = FALSE, na = "")

bib <- readLines(fs::path(out_dir, "tidymarkers.bib"), warn = FALSE)
if (!any(grepl("^@[a-zA-Z]+\\{polymarker_source,", bib))) {
  stop("polymarker_design_rules.csv cites polymarker_source but the bib has no entry")
}

message(
  "polymarker: ", nrow(out), " rules from source; subgenome axis on the common primer (",
  sum(out$role == "common-primer"), " rules), allele axis on the AS primer (",
  sum(out$role == "allele-primer"), " rules)"
)
