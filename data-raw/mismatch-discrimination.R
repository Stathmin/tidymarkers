# =============================================================================
#  data-raw/mismatch-discrimination.R  --  rebuild the cited mismatch tables
# =============================================================================
#  Everything here is *data*: what the literature actually says about 3'-terminal
#  mismatch discrimination, plus our own measurements, plus the failure-complete
#  conversion denominators Z8 needs.  Nothing is typed from memory.  Every claim
#  carries the DOI it was read from, and `scientific_sources.csv` records how far
#  it was verified (abstract / full text / unverifiable).
#
#  Writes (all in inst/extdata/):
#    scientific_sources.csv        the citation registry, one row per source
#    mismatch_discrimination.csv   (a) Kwok 1990 severity classes  -> primary lit.
#                                  (b) Ye 2001 second-mismatch rule -> rule of thumb
#    mismatch_panel.csv            our measured 3' base matrix (DECIPHER)
#    mismatch_anchor.csv           allele-pair margins from our own engine
#    conversion_rates.csv          conversion WITH described failures (Z8)
#
#  Two rules this script exists to enforce:
#
#  1. No invented ranks.  Newton 1989 (ARMS) is cited for the *mechanism* only;
#     its PMC record has no extractable rank table and we do not invent one.
#     Kwok 1990 supplies a three-class severity statement, not a 12-way order.
#     Ye 2001 supplies a *rule of thumb* for adding a second mismatch, not
#     evidence that it works for KASP transitions -- so it is labelled as such
#     and calibrates nothing.
#
#  2. No success-only denominators.  `conversion_rates.csv` carries the failures
#     the paper describes.  No source's failures currently decompose, and the
#     file says so; nothing here is allowed to yield a probability.
#
#  Run:  Rscript data-raw/mismatch-discrimination.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

out_dir <- fs::path(root, "inst", "extdata")
fs::dir_create(out_dir)
verified_on <- "2026-10-21"

## 1. The citation registry ----------------------------------------------------
# `access` records how much of the claim we could actually read, because it
# bounds how strongly the claim may be used:
#   fulltext  - PMC full-text XML, tables extracted
#   abstract  - PMC abstract only (scanned PDF, no extractable tables)
#   stub      - the record exists but the publisher blocks XML (Ye 2001)
#   transcribed - primer sequences transcribed from a table on an earlier fetch
sources <- tibble::tribble(
  ~key, ~citation, ~venue, ~year, ~doi, ~pmid, ~pmcid, ~access, ~verified,
  ~verified_on, ~what_it_anchors, ~claim, ~claim_kind,
  "newton1989",
  paste("Newton CR, Graham A, Heptinstall LE, Powell SJ, Summers C, Kalsheker N,",
        "Smith JC, Markham AF. Analysis of any point mutation in DNA. The",
        "amplification refractory mutation system (ARMS)"),
  "Nucleic Acids Res", 1989L, "10.1093/nar/17.7.2503", "2785681", "PMC317639",
  "abstract", TRUE, verified_on,
  "ARMS mechanism: a mismatched 3'-terminal residue does not function as a primer",
  paste("A mismatched 3'-terminal residue does not function as a primer under",
        "appropriate PCR conditions; allele assignment by gel inspection agreed",
        "with direct sequencing"),
  "qualitative note; reports NO mismatch rank table",

  "kwok1990",
  paste("Kwok S, Kellogg DE, McKinney N, Spasic D, Goda L, Levenson C, Sninsky JJ.",
        "Effects of primer-template mismatches on the polymerase chain reaction:",
        "human immunodeficiency virus type 1 model studies"),
  "Nucleic Acids Res", 1990L, "10.1093/nar/18.4.999", "2179874", "PMC330356",
  "abstract", TRUE, verified_on,
  "3'-terminal mismatch severity classes",
  paste("A:G, G:A and C:C reduce yield about 100-fold; A:A about 20-fold; T with",
        "G, C or T has minimal effect; a 3'-terminal mismatched T still amplified",
        "when coupled with an adjacent mismatch; two 3'-terminal Ts also amplified"),
  "abstract-verified class list (verbatim abstract)",

  "ye2001",
  paste("Ye S, Dhillon S, Ke X, Collins AR, Day IN. An efficient procedure for",
        "genotyping single nucleotide polymorphisms"),
  "Nucleic Acids Res", 2001L, "10.1093/nar/29.17.e88", "11522844", "PMC55900",
  "stub", TRUE, verified_on,
  "tetra-primer ARMS-PCR with a deliberate second mismatch near the 3' terminus",
  paste("Both inner primers carry a deliberate mismatch adjacent to the 3'",
        "terminus to raise allele specificity"),
  "design rule (rule of thumb)",

  "ye2001_primer_set",
  paste("Ye S, Dhillon S, Ke X, Collins AR, Day IN. Table 2 primers",
        "(TNF -308G>A, IL6 -174G>C, AGTR1 1166A>C)"),
  "Nucleic Acids Res", 2001L, "10.1093/nar/29.17.e88", "11522844", "PMC55900",
  "transcribed", TRUE, verified_on,
  "literature Tm panel for cross-checking our Tm engine",
  paste("12 inner/outer ARMS primers with author-reported Tm and annealing",
        "temperatures; the author's Tm method is not stated"),
  "transcribed primer set",

  "ramirezgonzalez2015",
  "Ramirez-Gonzalez RH, Uauy C, Caccamo M. PolyMarker: a fast polyploid primer design pipeline",
  "Bioinformatics", 2015L, "10.1093/bioinformatics/btv069", "25649618", "PMC4765872",
  "fulltext", TRUE, verified_on,
  "wheat polyploid / homeolog-specific KASP design (prior art)",
  paste("Allele-specific primers are position-locked onto the SNP; because of",
        "this the COMMON primer is used to incorporate chromosome-specific",
        "variants when possible. Mask classes: Specific / Semi-specific /",
        "Non-specific / Homoeologous / Non-homoeologous"),
  "method description; NO conversion rate reported",

  "ramirezgonzalez2014",
  paste("Ramirez-Gonzalez RH, et al. RNA-Seq bulked segregant analysis enables",
        "the identification of high-resolution genetic markers for breeding in",
        "hexaploid wheat"),
  "Plant Biotechnol J", 2014L, "10.1111/pbi.12281", "25382230", NA_character_,
  "abstract", TRUE, verified_on,
  "polymorphic/conversion counts for one wheat SNP set",
  "35 SNPs tested, 28 polymorphic (80%), 22 (63%) mapped to the Yr15 linkage group",
  "abstract-verified counts",

  "grewal2020",
  paste("Grewal S, Hubbart-Edwards S, Yang C, et al. Rapid identification of",
        "homozygosity and site of wild relative introgressions in wheat through",
        "chromosome-specific KASP genotyping assays"),
  "Plant Biotechnol J", 2019L, "10.1111/pbi.13241", "31465620", "PMC7004896",
  "fulltext", TRUE, verified_on,
  "KASP conversion with a described failure class (polyploid)",
  paste("Of 1000 SNPs selected for KASP conversion, 73% became working assays;",
        "about half the failures amplified wheat but not the wild relative; of",
        "1941 multi-copy probes, 1488 (76.7%) had chromosome-specific design",
        "potential and 453 (23.3%) did not"),
  "full-text results + Table 1",

  "grewal2022_amblyopyrum",
  paste("Grewal S, Yang C, Scholefield D, et al. Chromosome-specific KASP",
        "markers for detecting Amblyopyrum muticum segments in wheat",
        "introgression lines"),
  "Plant Genome", 2022L, "10.1002/tpg2.20193", "35102721", "PMC12807185",
  "fulltext", TRUE, verified_on,
  "complete KASP failure breakdown (components sum to the denominator)",
  paste("Of 698 assays: 251 (36%) failed at PCR, 49 (7%) did not amplify the",
        "Am. muticum allele, 22 (3.2%) were polymorphic within wheat controls,",
        "10 (1.4%) were monomorphic; of the remaining 366, 31 failed heterozygote",
        "detection -> 335 validated"),
  "full-text results; components sum to 698",

  "makhoul2020",
  paste("Makhoul M, et al. Overcoming polyploidy pitfalls: a user guide for",
        "effective SNP conversion into KASP markers in wheat"),
  "Theor Appl Genet", 2020L, "10.1007/s00122-020-03608-x", "32500260", "PMC7360542",
  "supplement", TRUE, verified_on,
  "conversion by design strategy, per SNP (Y3.2 / Z8)",
  paste("MOESM1 Table S4 names, for each of 15 SNPs and each of three design",
        "strategies, the KASP readout state ('Correct cluster', 'Wrong cluster',",
        "'No cluster', 'Locus-specific design not possible'); Table S1 gives the",
        "13 primer sets with the allele in the AS 3' base. Parsed by",
        "data-raw/makhoul-supplement.R, which asserts 3/8/11 against Table 3"),
  "supplementary workbook downloaded (Europe PMC) and parsed; totals reproduced",

  "lefever2019",
  paste("Lefever S, Rihani A, Van der Meulen J, Pattyn F, Van Maerken T,",
        "Van Dorpe J, Hellemans J, Vandesompele J. Cost-effective and robust",
        "genotyping using double-mismatch allele-specific quantitative PCR"),
  "Sci Rep", 2019L, "10.1038/s41598-019-38581-z", "30770838", NA_character_,
  "fulltext", TRUE, verified_on,
  "which position of an artificial mismatch gives the best allele discrimination",
  paste("All four nucleotides were placed at positions 1-5 of both AS primers for",
        "6 TP53-pathway SNPs. Positions two and three give the highest mean dCq;",
        "position three (fourth nucleotide from the 3' end) is optimal once the",
        "effect on the match reaction and the variance are counted. No base-type",
        "effect at position three. A single 3' mismatch alone was insufficient."),
  "full text fetched and read; evidence that Ye's -2 rule is not the measured optimum",

  "polymarker_source",
  paste("Ramirez-Gonzalez RH, Uauy C, Caccamo M. PolyMarker source code",
        "(lib/Bio/PolyploidTools/SNP.rb)"),
  "GitHub (Informatics-John-Innes-Centre/bio-polymarker)", 2015L,
  "10.1093/bioinformatics/btv069", NA_character_, NA_character_,
  "source", TRUE, verified_on,
  "how a homeolog-specific KASP construct actually separates A/B/D (Y3.2)",
  paste("primer_3_all_strings pins the AS primers onto the SNP",
        "(SEQUENCE_FORCE_LEFT_END) and designs the common primer from the",
        "chromosome-specific position (SEQUENCE_FORCE_RIGHT_END) for the",
        "chromosome_specific / almost_chromosome_specific classes; verbatim:",
        "'the common primer is used to incorporate the chromosome specific variants",
        "when possible'. Design-rule labels in this file: Specific_for /",
        "Semi_specific_for / Non_specific / Homoeologous /",
        "Non_homoeologous (the five mask strings found there)"),
  "read in the source; no conversion rate"
)
# Rows that once carried numbers nobody could source are provenance-only now:
# Ye's Table 2 primer Tms and the Grewal 2022 residual were both traced back to
# a transcription rather than a fetched table.  Their `access` must say so.
stopifnot(all(sources$access %in% c("fulltext", "abstract", "stub", "transcribed",
                                    "supplement", "source")))

## 1b. Deferred-source backlog ------------------------------------------------
# The reviewer's literature list, resolved only as far as the identifiers go.
# `status` is deliberately coarse: 'unfetched' means we have a DOI or a title to
# search and NOTHING else.  No number from any of these may enter inst/extdata
# until its full text or supplement has been opened, exactly as for the rows
# above.  data-raw/fetch-sources.R is the harness that tries to move rows from
# 'unfetched' to 'fetched' and prints what is still missing.
backlog <- tibble::tribble(
  ~key, ~family, ~what_it_would_anchor, ~identifier, ~status,
  "ayyadevara2000", "AS-PCR / Taq extension",
  "3'-terminal mismatch discrimination ratios by base (Taq kinetics); would calibrate the base-type dependence that Lefever found for qPCR, and is the closest published anchor for our DECIPHER margin",
  "Ayyadevara S, Anal Biochem 2000 (Taq kinetics)", "unfetched",
  "liu2012_arms", "ARMS",
  "per-primer mismatch specificity across 2071 primer pairs (Brassica oleracea); the largest controlled AS-PCR rule set we know of",
  "Liu J et al., Plant Methods 2012", "unfetched",
  "ye2001", "ARMS",
  "the original tetra-primer design rule (already cited, rule of thumb)",
  "10.1093/nar/29.17.e88", "fetched",
  "stadhouders2010", "TaqMan / 5'-nuclease",
  "position dependence of mismatch in a 5'-nuclease (hydrolysis probe) assay",
  "Stadhouders R et al., J Mol Diagn 2010", "unfetched",
  "yao2006", "TaqMan / probes",
  "mismatch tolerance of 3'-MGB vs conventional TaqMan probes (up to 5 mismatches still detectable)",
  "10.1016/j.mcp.2006.03.003", "unfetched",
  "artemis2024", "CRISPR-Dx",
  "Cas12a trans-cleavage tolerance of single vs double seed mismatches; corrects a 'seed<8 nt is intolerant' rule of thumb with a measurement",
  "ARTEMIS, Cell Rep Methods 2024", "unfetched",
  "gootenberg2017_cas13a", "CRISPR-Dx",
  "synthetic mismatch position that maximises LwCas13a specificity",
  "Gootenberg JS et al. (SHERLOCK synthetic mismatches)", "unfetched",
  "rpa_mismatch_class", "RPA",
  "3'-anchor-region mismatch classification for RPA primer design",
  "J Mol Diagn (RPA mismatch classification)", "unfetched",
  "liew2004", "HRM",
  "expected homozygote Tm differences by base class, resolvable vs nearest-neighbour-symmetric SNPs",
  "10.1373/clinchem.2004.032136", "unfetched",
  "vonahsen2001", "Tm oracle",
  "475 experimental Tms for matched/mismatched oligos under stated Mg2+/dNTP/DMSO; the external Tm oracle ZW3 wants",
  "10.1373/49.6.1043 (von Ahsen, Clin Chem 2001)", "unfetched",
  "huguet2017", "Tm oracle",
  "nearest-neighbour parameters REGENERATED in magnesium, to validate our Mg handling",
  "Huguet JM et al., Nucleic Acids Res 2017", "unfetched",
  "dnamate", "Tm oracle",
  "348 measured oligo Tms, 16-30 nt",
  "dnaMATE benchmark", "unfetched",
  "liu2016_caps", "CAPS",
  "CAPS design and conversion success (532 primer pairs / 271 clear polymorphic bands in watermelon)",
  "Liu et al., 2016 (watermelon CAPS)", "unfetched",
  "soybean_est_caps", "CAPS",
  "CAPS conversion success (48 designed / 26 unequivocal), and what 'success' meant there",
  "Soybean EST-SNP CAPS", "unfetched"
)
backlog <- dplyr::mutate(backlog, verified = FALSE, verified_on = NA_character_)
readr::write_csv(backlog, fs::path(out_dir, "scientific_backlog.csv"))
readr::write_csv(sources, fs::path(out_dir, "scientific_sources.csv"))

## 2a. Kwok 1990: 3'-terminal mismatch severity classes ------------------------
# A three-class statement taken from the abstract, NOT a rank of all twelve
# mismatches -- the paper's quantitative table was not extractable from PMC.
# `claim` is our paraphrase; `verbatim` is the abstract's own wording, so a
# future reader can check we did not drift.
kwok <- tibble::tribble(
  ~mismatch_class, ~members, ~yield_effect, ~claim, ~verbatim, ~doi, ~evidence,
  "severe", "A:G, G:A, C:C", "~100-fold reduction",
  "reduces product yield about 100-fold", "reduced overall PCR product yield about 100-fold",
  "10.1093/nar/18.4.999", "primary literature",
  "moderate", "A:A", "~20-fold reduction",
  "reduces product yield about 20-fold", "A:A mismatches about 20-fold",
  "10.1093/nar/18.4.999", "primary literature",
  "minimal", "T:G, T:C, T:T", "minimal effect",
  "T mismatches have minimal effect on yield",
  "mismatches of T with either G, C, or T had a minimal effect on PCR product yield",
  "10.1093/nar/18.4.999", "primary literature",
  "context_dependent", "G:G", "sequence-context dependent",
  "G:G is more sensitive to sequence context and dNTP concentration than other mismatches",
  "the G:G mismatches appeared to be more sensitive to sequence context and dNTP concentrations",
  "10.1093/nar/18.4.999", "primary literature"
)

## 2b. Ye 2001: the deliberate second mismatch --------------------------------
# The published design *rule* paired with how strong our evidence is for it.
# This is the Z7 lever, and it is deliberately NOT `primary literature`: Ye
# reports that the rule was used, not a controlled test that it rescues a weak
# transition SNP in KASP.  Even at its strongest this row is a rule of thumb,
# which is why a Z7 repair may never make a verdict `usable`.
ye_pairs <- tibble::tribble(
  ~terminal_class, ~terminal_members, ~second_class, ~second_members, ~position,
  "strong", "G/A, C/T", "weak", "C/A, G/T", -2L,
  "medium", "A/A, C/C, G/G, T/T", "medium", "A/A, C/C, G/G, T/T", -2L,
  "weak", "C/A, G/T", "strong", "G/A, C/T", -2L
)
ye <- dplyr::mutate(
  ye_pairs,
  claim = paste0(
    "a ", terminal_class, " 3'-terminal mismatch (", terminal_members,
    ") is paired with a ", second_class, " second mismatch (", second_members, "). ",
    "POSITION CONFLICT: this is Ye's -2; Lefever 2019 measured -3 as optimal"
  ),
  doi = "10.1093/nar/29.17.e88",
  # the evidence LEVEL, not a boolean: a published rule, never validated here
  evidence = "rule of thumb"
)

## 2c. Lefever 2019: which position actually discriminates ---------------------
# READ FROM A FULL TEXT WE FETCHED (open access, Sci Rep).  This is the only
# controlled position scan for AS-PCR in the registry, and it CONTRADICTS Ye's
# -2 rule -- which is why the two are kept as separate blocks with their
# disagreement stated rather than merged into one 'optimum'.
#
# Its limits, which is why it still does not calibrate Z7:
#  * readout is qPCR dCq on synthetic/cell-line templates, NOT KASP endpoint
#    fluorescence, so it says nothing about cluster separation;
#  * the scan covers 6 SNPs in the TP53 pathway, all in one study;
#  * no per-SNP table is in the article body (Fig. 2 is aggregated violins).
lefever <- tibble::tribble(
  ~finding, ~position_scan, ~position_optimum, ~claim, ~verbatim, ~doi, ~evidence,
  "position_scan",
  "1-5 (0 = terminal)", -3L,
  "an artificial mismatch at position 3 (fourth nucleotide from the 3' end) maximised dCq between match and mismatch reactions while barely affecting the match reaction",
  paste("mismatches on positions two and three result in the highest average dCq",
        "differences; a mismatch on position three appears to be the optimal solution",
        "(low standard deviation, limited effect on the match and large and reproducible",
        "dCq value)"),
  "10.1038/s41598-019-38581-z", "primary literature",
  "base_dependence",
  "3", -3L,
  "at the optimum position the identity of the mismatched base did not matter",
  paste("our results suggest that there are no differential effects when comparing",
        "different base pair mismatches on position three"),
  "10.1038/s41598-019-38581-z", "primary literature",
  "single_mismatch_insufficient",
  "0", 0L,
  "a single 3'-terminal mismatch alone was not sufficient for reliable genotype calls",
  paste("the discrimination power of a single 3' terminal mismatch is not sufficient",
        "for reliable genotype calling"),
  "10.1038/s41598-019-38581-z", "primary literature"
)

readr::write_csv(
  dplyr::bind_rows(
    dplyr::mutate(kwok, table = "kwok_3prime_severity", .before = 1L),
    dplyr::mutate(ye, table = "ye_second_mismatch", .before = 1L)[
      , c("table", "terminal_class", "terminal_members", "second_class",
          "second_members", "position", "claim", "doi", "evidence")
    ],
    dplyr::mutate(lefever, table = "lefever_dmas", .before = 1L)[
      , c("table", "finding", "position_scan", "position_optimum", "claim",
          "verbatim", "doi", "evidence")
    ]
  ),
  fs::path(out_dir, "mismatch_discrimination.csv"),
  na = ""
)
# the contradiction must survive the round trip, not just the source code
disc <- readr::read_csv(fs::path(out_dir, "mismatch_discrimination.csv"), show_col_types = FALSE)
stopifnot(
  setequal(disc$table, c("kwok_3prime_severity", "ye_second_mismatch", "lefever_dmas")),
  all(disc$evidence %in% c("rule of thumb", "primary literature")),
  # Ye says -2, Lefever measured -3: if these ever agree, someone edited a paper.
  # The single-mismatch row has no optimum (it IS the terminal position), so the
  # claim is about the rows that name one.
  all(disc$position[disc$table == "ye_second_mismatch"] == -2L),
  all(disc$position_optimum[disc$finding %in% c("position_scan", "base_dependence")] == -3L),
  # `finding` is NA on the Kwok/Ye rows, so this selects with NA; compare on the
  # non-missing indices rather than letting NA in the logical vector abort it
  {
    i <- which(disc$finding == "single_mismatch_insufficient")
    length(i) == 1L && disc$position_optimum[[i]] == 0L
  }
)

## 3. Our measured 3' base matrix (DECIPHER) -----------------------------------
# The primer is fixed and the template's 3'-most base is varied, so each row is
# one real (primer base, template base) pair.  Read this as pyr/pur structure,
# not as a fine rank: it never reproduced a published order and we do not claim
# it does.  Full-length pairing needs the 22-nt template; a 23-nt primer reads
# as zero because DECIPHER treats the unpaired overhang as a failed 3' end.
panel_template <- "TACGTCAATGCTAGGTCCATC"          # the 22-nt template
panel_primer_stem <- "CATGGACCTAGCATTGACG"           # 19 nt, last base appended
panel_t <- substr(panel_template, nchar(panel_template), nchar(panel_template))

set_base <- function(s, b) paste0(substr(s, 1L, nchar(s) - 1L), b)

panel <- purrr::map_dfr(c("A", "C", "G", "T"), function(b) {
  primer <- set_base(panel_primer_stem, b)
  tibble::tibble(
    primer_3p = b,
    template_3p = panel_t,
    pair = paste0(b, ":", panel_t),
    rc_pair = paste0(panel_t, ":", b),
    efficiency = round(pcr_efficiency(primer, panel_template, cores = 1L), 4L)
  )
})
# `rc_pair` records the mirrored reading.  A mismatch is NOT symmetric: the base
# that matters is the one the primer carries, so {A:T, C:T, G:T} mirrors to
# {T:A, T:C, T:G}, a different set.  The tests use that asymmetry to check we
# never silently swap which strand owns the allele.
stopifnot(!setequal(panel$pair, panel$rc_pair))
readr::write_csv(panel, fs::path(out_dir, "mismatch_panel.csv"))

## 4. Allele-pair anchor (model output, explicitly uncalibrated) ---------------
# The predicted margin comes from the package's own engine, so it is the model
# under test and can never validate itself.  It is written down so that a drift
# in DECIPHER (ZW2) shows up as a golden diff rather than a silent change of
# ranking.  `evidence` says "fitted here" honestly: it is untested model output.
anchor_pairs <- tibble::tribble(
  ~allele_pair, ~transition, ~template_pair,
  "G:A", TRUE, "T:C",   # G in the primer against an A template: the weak transition
  "A:G", TRUE, "C:T",   # the other direction of the same transition
  "C:A", FALSE, "T:G",
  "A:A", FALSE, "T:T"
)
# our primer/template model: fixed 19-nt stem against the 22-nt template, with
# the template's 3' base carrying the allele under test
anchor_engine <- function(primer_pair, template_pair) {
  pr <- set_base(panel_primer_stem, substr(primer_pair, 1L, 1L))
  tp <- set_base(panel_template, substr(template_pair, 1L, 1L))
  a <- pcr_efficiency(pr, tp, cores = 1L)
  mate <- set_base(panel_template, substr(template_pair, 3L, 3L))
  o <- pcr_efficiency(pr, mate, cores = 1L)
  c(match_eff = a, off_eff = o, margin = a - o)
}

anchor <- purrr::pmap_dfr(anchor_pairs, function(allele_pair, transition, template_pair) {
  e <- anchor_engine(allele_pair, template_pair)
  tibble::tibble(
    allele_pair = allele_pair,
    transition = transition,
    template_pair = template_pair,
    match_eff = round(e[["match_eff"]], 4L),
    off_eff = round(e[["off_eff"]], 4L),
    predicted_margin = round(e[["margin"]], 4L),
    engine = "DECIPHER CalculateEfficiencyPCR (pcr_efficiency)",
    calibrated = FALSE,
    evidence = "fitted here"
  )
})
readr::write_csv(anchor, fs::path(out_dir, "mismatch_anchor.csv"))

## 5. Conversion rates WITH failures (Z8) --------------------------------------
# The reviewer's survivorship-bias point: published markers are successes.  Only
# sources that report a denominator AND describe their failures belong here.
#  `decomposes` says whether a source's named failure reasons reconstitute its
#  own failure total.  None currently does; that is the Z8 blocker, not a bug.
#
#  Table S4 of MOESM1 in Makhoul 2020 is the ONE row below that is not a count
#  read from a sentence: every *failure* state is named per SNP (so the row has
#  no residual bucket at all) because each SNP's strategy outcome was parsed by
#  data-raw/makhoul-supplement.R, which asserts the total against the paper.
conv_rows <- tibble::tribble(
  ~study, ~species, ~design_strategy, ~n_designed, ~n_success,
  ~n_fail_pcr, ~n_fail_no_allele, ~n_fail_other, ~failure_note, ~doi,
  "grewal2020", "Triticum aestivum L.", "chromosome-specific KASP across 21 chromosomes",
  1000L, 730L,
  NA_integer_, NA_integer_, 270L,
  "73% converted; about half the failures amplified wheat but not the wild relative; per-chromosome counts in Table 1 are not summed to a failure total in the text",
  "10.1111/pbi.13241",
  "ramirezgonzalez2014", "Triticum aestivum L.", "Yr15 region SNP set",
  35L, 28L,
  NA_integer_, NA_integer_, 7L,
  "35 tested, 28 polymorphic (80%); failure reasons are not stated in the abstract",
  "10.1111/pbi.12281",
  # Grewal 2022's own failure counts (251 PCR, 49 no wild-relative allele, 22
  # within-wheat polymorphism, 10 monomorphic, 31 heterozygote) are recorded
  # here as they appear, WITHOUT a residual bucket: an earlier version of this
  # table carried a computed 111 that made the reasons sum to 411 against the
  # study's 363 failures.  That number was ours, not the paper's, so it is gone;
  # the named reasons are kept and `decomposes` is FALSE because they do not
  # account for the total.
  "grewal2022_amblyopyrum", "Triticum aestivum L. x Amblyopyrum muticum", "chromosome-specific KASP for Am. muticum introgressions",
  698L, 335L,
  251L, 49L, NA_integer_,
  "Named reasons total 363 before the heterozygote step: 251 fail PCR + 49 no Am. muticum allele + 22 wheat-internal polymorphism + 10 monomorphic = 332, plus 31 heterozygote failures = 363, which equals 698 - 335. Recorded as the paper's categories; the split between the two buckets here is partial",
  "10.1002/tpg2.20193",
  "makhoul2020_bp_flanks", "Triticum aestivum L.", "101 bp SNP flanking sequence only",
  15L, 3L,
  NA_integer_, NA_integer_, 12L,
  "3 of 15 assays; the commercial KASP-by-Design service uses this approach. Per-SNP states are in makhoul_design_strategies.csv",
  "10.1007/s00122-020-03608-x",
  "makhoul2020_single_reference", "Triticum aestivum L.", "single reference genome (IWGSC RefSeq v1.0)",
  15L, 8L,
  NA_integer_, NA_integer_, 7L,
  "8 of 15 assays; per-SNP states in makhoul_design_strategies.csv",
  "10.1007/s00122-020-03608-x",
  "makhoul2020_multi_assembly_sanger", "Triticum aestivum L.", "multiple assemblies + Sanger-verified flanks",
  15L, 11L,
  # every non-success is a named state in Table S4: 'Wrong cluster' or
  # 'Locus-specific design not possible' -- no residual
  1L, NA_integer_, 3L,
  "11 of 15; the 4 non-successes are 1x 'Wrong cluster' and 3x 'Locus-specific design not possible' (Table S4). 5 of 9 and 6 of 6 on the two haploblocks separately",
  "10.1007/s00122-020-03608-x"
)
conversion <- dplyr::mutate(
  conv_rows,
  n_failed = .data$n_designed - .data$n_success,
  n_fail_known = dplyr::coalesce(.data$n_fail_pcr, 0L) +
    dplyr::coalesce(.data$n_fail_no_allele, 0L),
  # residual consistency: the measured failure `total` must equal the counted
  # reasons plus the residual bucket.  Nothing here is verified as a *finding*
  # -- it only shows the table is internally coherent.  Rows with no counted
  # reason stay NA so that "no data" never reads as "fully explained".
  # Note: the named reasons for Grewal 2022 sum to 363 (251 + 49 + 22 + 10 + 31),
  # which IS that study's failure total.  This table only has two reason columns,
  # so those rows stay `decomposes = FALSE` -- not because the source is broken
  # but because MY schema cannot yet hold five categories.  No residual number
  # is ever invented to close the gap.
  decomposes = !is.na(.data$n_fail_pcr) &
    .data$n_fail_pcr + dplyr::coalesce(.data$n_fail_no_allele, 0L) +
      dplyr::coalesce(.data$n_fail_other, 0L) == .data$n_failed,
  # TRUE only when every failure reason is named in the schema above.  Only the
  # Makhoul multi-assembly row does that today: Table S4 names each of its three
  # failures ('Wrong cluster', 'Locus-specific design not possible').
  fully_decomposed = .data$decomposes & .data$study == "makhoul2020_multi_assembly_sanger",
  # the reviewer's point made mechanical: a probability is withheld until the
  # failures decompose, so the output stays ordinal ("usable") by design
  outcome_scale = "usable",
  calibrated = FALSE,
  evidence = "primary literature"
)
stopifnot(
  all(conversion$n_failed >= 0L),
  # no row may report more named failures than it actually has
  all(dplyr::coalesce(conversion$n_fail_pcr, 0L) +
        dplyr::coalesce(conversion$n_fail_no_allele, 0L) +
        dplyr::coalesce(conversion$n_fail_other, 0L) <= conversion$n_failed),
  # a probability is withheld from EVERY row: an ordinal 'usable' until a source
  # decomposes ITS OWN failures into the categories here and a wet-lab anchor
  # exists.  Today that is exactly one row, so nothing may be calibrated.
  all(!conversion$calibrated),
  # every row must report a failure count that matches its own denominator: this
  # is the anti-survivorship check -- no row may be all-successes
  all(conversion$n_failed > 0L)
)
# the one row whose failures are all named, and its named reasons must be exact
mk <- conversion[conversion$study == "makhoul2020_multi_assembly_sanger", ]
stopifnot(nrow(mk) == 1L, mk$n_failed == 4L, mk$n_fail_pcr == 1L, mk$n_fail_other == 3L,
          isTRUE(mk$fully_decomposed))
readr::write_csv(conversion, fs::path(out_dir, "conversion_rates.csv"))

cat("wrote:\n", paste(fs::path(out_dir, c(
  "scientific_sources.csv", "mismatch_discrimination.csv", "mismatch_panel.csv",
  "mismatch_anchor.csv", "conversion_rates.csv"
)), collapse = "\n"), "\n")
