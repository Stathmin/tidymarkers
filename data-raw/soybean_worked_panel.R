# =============================================================================
#  data-raw/soybean_worked_panel.R  --  the ordered soybean assay panel
# =============================================================================
#  Source: "Напоминалка по праймерам" (2026-09-10), 14 PDF pages, compiled from
#  the saved selection tables for the soybean (Glycine max) projects.  The
#  header of that document states that every number in it was computed, not
#  typed by hand.  Genome build: Williams 82 Wm82.a4.v1 (GCF_000004515.6).
#
#  What this panel is, and what it is NOT
#  --------------------------------------
#  It is a real design panel: 36 assays and 85 oligos that were selected for
#  ordering (15 Sanger amplicons, 10 KASP sets on 24 plants, 3 KASP sets on the
#  Dt2 flanks, 8 reserve Dt2 body amplicons).  Every assay carries the numbers
#  the document prints: primer sequences, Tm, amplicon/product size, the KASP
#  gap, the printed hairpin Tm and, where present, the GBS class counts.  It also
#  prints exact-genome-hit counts for the ten KASP 2.1 sets; those are kept as
#  `genome_hits_printed`, while `genome_hits` is recomputed from the assembly
#  itself (bullet (d) below), not read from the sentence.
#
#  It is NOT a wet-lab outcome table.  The document does not print a per-assay
#  call, a band score, or a success/failure verdict.  The operator reports that
#  the panel worked, but that report is not a measured table and is not turned
#  into a rate here.  The extract is therefore a DESIGN benchmark: it lets the
#  package re-derive the arithmetic and the Tm of a panel that was actually
#  ordered, without claiming that a design equals a working assay.
#
#  The four assertions that make this an extract rather than a transcription:
#
#    (a) every Sanger product size equals end - start + 1, and the 10 peak
#        assays also satisfy band = offset_F + offset_R + len(F) + len(R) - 1;
#    (b) every KASP product satisfies product = gap + len(AS body) + len(COM),
#        and the two allele-specific bodies differ only at their 3' base, which
#        is the stated allele;
#    (c) the printed Tm is reproduced by primer_tm() with the SantaLucia (1998)
#        salt correction to better than 0.15 C on all 85 oligos.  With the
#        package default (Schildkraut 2010) the same set reads ~3.7 C low, so
#        the panel pins the salt model, not merely the nearest-neighbour table;
#    (d) `genome_hits` is the exact full-length occurrence count of each primer
#        body in GCF_000004515.6, both strands, computed here with the package's
#        indexed reader and exact-hit table.  Every non-ALT oligo has exactly
#        one match (its target locus) and `genome_hit_coord`/`genome_hit_strand`
#        give where it sits; the ALT-allele-specific primers have zero where the
#        reference carries the other allele.  The document's 30 printed counts
#        are asserted equal to that computation.
#    (e) each of the 23 Sanger primer pairs amplifies exactly one exact product,
#        on its own chromosome and at the printed coordinates, confirmed through
#        in_silico_pcr(); genome-wide uniqueness of both primers rules out any
#        other exact product.
#    (f) each of the 13 KASP sets is graded against the package's own acceptance
#        gates by marker_design_checks(): the amplicon span, in-silico PCR and
#        off-target counts come from the assembly, the AS on/off windows give
#        the discrimination, and the structure/Tm gates use the panel sequences.
#        The panel carries no population variant table, so the footprint and
#        linkage checks are reported `unknown`, never passed.  The gate table is
#        `soybean_worked_panel_kasp_gates.csv`.
#
#  The PDF is .tmp/lit/soybean_worked_panel/paper.pdf; this script reads the
#  `pdftotext -layout` twin.  Recreate it with:
#    pdftotext -layout paper.pdf paper.txt
#
#  Run:  Rscript data-raw/soybean_worked_panel.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

txt_path <- fs::path(root, ".tmp", "lit", "soybean_worked_panel", "paper.txt")
if (!file.exists(txt_path)) {
  stop(
    "Missing ", txt_path,
    " - run: pdftotext -layout paper.pdf paper.txt"
  )
}

## --- the assembly, for the exact hits ----------------------------------------
# `genome_hits` is not copied from the document.  The exact full-length
# occurrence of every primer body, with its chromosome and coordinates, is
# computed against Wm82.a4.v1 with the package's own indexed reader
# (`marker_genome()` / `genome_seq()`) and exact-hit table (`exact_hit_table()`),
# one chromosome at a time.  The 946 MB FASTA is not in the repository: it is
# cached under `marker_cache_dir()` (like `ath_genome()`) and indexed once.
# Without the cache this script downloads it; the test that sources this script
# skips instead of downloading.
soy_cache <- fs::path(marker_cache_dir(), "gm_v4")
soy_assembly <- "GCF_000004515.6_Glycine_max_v4.0"
soy_fasta <- fs::path(soy_cache, paste0(soy_assembly, "_genomic.fna"))
soy_report <- fs::path(soy_cache, "assembly_report.txt")

fetch_soybean_genome <- function() {
  fs::dir_create(soy_cache)
  base <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/", soy_assembly
  )
  gz <- paste0(soy_fasta, ".gz")
  if (!fs::file_exists(soy_fasta) && !fs::file_exists(gz)) {
    cli::cli_alert_info(gettext("Downloading Wm82.a4.v1 genome (~283 MB) ..."))
    utils::download.file(
      paste0(base, "/", soy_assembly, "_genomic.fna.gz"), gz,
      mode = "wb", quiet = TRUE
    )
  }
  if (!fs::file_exists(soy_report)) {
    utils::download.file(
      paste0(base, "/assembly_report.txt"), soy_report,
      quiet = TRUE
    )
  }
  # fasta_index() decompresses the .gz once and builds the .fai
  fasta_index(if (fs::file_exists(soy_fasta)) soy_fasta else gz)
}

# Chromosome number -> RefSeq accession, straight from the NCBI assembly report.
# Only the assembled molecules are kept; the mitochondrion and chloroplast do
# not contain any of these primers (checked: adding all 284 records leaves every
# count unchanged), so the count runs on the 20 chromosomes alone.
soybean_chrom_map <- function(report) {
  tab <- utils::read.delim(report,
    comment.char = "#", header = FALSE, stringsAsFactors = FALSE
  )
  chr <- tab[
    tab[[2L]] == "assembled-molecule" & grepl("^[0-9]+$", tab[[3L]]), ,
    drop = FALSE
  ]
  stats::setNames(as.character(as.integer(chr[[3L]])), chr[[7L]])
}

soybean_genome <- function() {
  if (!fs::file_exists(paste0(soy_fasta, ".fai"))) {
    fetch_soybean_genome()
  }
  if (!fs::file_exists(soy_report)) {
    stop("Missing ", soy_report, " (NCBI assembly report)")
  }
  marker_genome(soy_fasta, chrom_map = soybean_chrom_map(soy_report))
}

# "Gm18:..." or "AGP1187691 Gm18:..." -> "18"
chrom_of_target <- function(x) {
  parts <- strsplit(x, "Gm", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) {
    return(NA_character_)
  }
  sub("^0*", "", regmatches(parts[[2L]], regexpr("^[0-9]+", parts[[2L]])))
}

lines <- readLines(txt_path, encoding = "UTF-8", warn = FALSE)
# PDF Candy emits a form feed where a page break falls inside a heading.
lines <- gsub("\f", "", lines)
trimmed <- trimws(lines)
heads <- grep("^(AMP|KASP)_[A-Za-z0-9_.]+$", trimmed)
if (!length(heads)) stop("no assay headings found in ", txt_path)

blocks <- lapply(seq_along(heads), function(i) {
  end <- if (i < length(heads)) heads[[i + 1L]] - 1L else length(lines)
  lines[heads[[i]]:end]
})
names(blocks) <- trimmed[heads]

flat <- function(b) paste(b, collapse = "\n")
# field patterns may run across a line break; collapse the whitespace once
flatw <- function(b) gsub("[[:space:]]+", " ", flat(b))

# Read one labelled field without letting `(.+)` swallow the rest of the block:
# labels sit either on the same line as their value or on a line of their own.
field_line <- function(b, labels) {
  pat <- paste0("^(", paste(labels, collapse = "|"), ")([[:space:]]|$)")
  i <- which(grepl(pat, b))
  if (!length(i)) {
    return(NULL)
  }
  i <- i[[1L]]
  rest <- trimws(sub(paste0("^(", paste(labels, collapse = "|"), ")[[:space:]]*"), "", b[[i]]))
  if (nzchar(rest)) {
    return(rest)
  }
  j <- i + 1L
  while (j <= length(b) && !nzchar(trimws(b[[j]]))) j <- j + 1L
  if (j <= length(b)) trimws(b[[j]]) else NULL
}

caps <- function(b, pattern, what) {
  x <- regmatches(flatw(b), regexec(pattern, flatw(b), perl = TRUE))[[1L]]
  if (!length(x)) stop("could not read ", what, " from ", b[[1L]])
  x[-1L]
}

opt_caps <- function(b, pattern) {
  x <- regmatches(flatw(b), regexec(pattern, flatw(b), perl = TRUE))[[1L]]
  if (!length(x)) NULL else x[-1L]
}

ru_num <- function(x) as.numeric(gsub(",", ".", gsub("[[:space:]]", "", gsub("\u2212", "-", x))))
ru_int <- function(x) as.integer(gsub("[[:space:]]", "", gsub("\u2212", "-", x)))

# The document is Russian; the extract is English.  The mapping is explicit so
# a reader can compare it with the source labels.
pop_en <- c(
  "СК Веда" = "SK Veda", "Клондайк" = "Klondike",
  "СК Сарита" = "SK Sarita", "Юнка" = "Yunka",
  "СК Аврора" = "SK Aurora", "Парус" = "Parus",
  "СК Дока" = "SK Doka", "РАТС-02" = "RATS-02"
)
trait_en <- c(
  "высота растения" = "plant height",
  "тип роста" = "growth habit",
  "число ветвей" = "branch number",
  "число узлов" = "node number",
  "масса семян с растения" = "seed mass per plant",
  "содержание белка" = "protein content"
)
translate <- function(x, map) {
  if (is.null(x) || is.na(x)) {
    return(NA_character_)
  }
  for (k in names(map)) x <- gsub(k, map[[k]], x, fixed = TRUE)
  x <- gsub("[\u00d7]", "x", x)
  gsub("[\u2014]", "-", x)
}
translate_locus <- function(x) {
  if (is.null(x) || is.na(x)) {
    return(NA_character_)
  }
  x <- sub("^покрывает[[:space:]]+", "", x)
  x <- translate(x, trait_en)
  x <- gsub("второй пик белка объединённой модели B (статья 1)",
    "second protein peak, joint model B (article 1)",
    x,
    fixed = TRUE
  )
  x <- gsub("прокси пикового индела, r2=0,989",
    "proxy for the peak indel, r2 = 0.989",
    x,
    fixed = TRUE
  )
  x <- gsub("промотор Dt2, окно", "Dt2 promoter, window", x, fixed = TRUE)
  x <- gsub("от старта трансляции", "from the translation start", x, fixed = TRUE)
  x <- gsub("п\\. н\\.", "bp", x)
  x <- gsub("экзон 1 Dt2", "Dt2 exon 1", x, fixed = TRUE)
  x <- gsub("несинонимичная замена G/A (Ser → Asn)",
    "non-synonymous G/A (Ser -> Asn)", x,
    fixed = TRUE
  )
  x <- gsub("[\u2013]", "-", x)
  gsub(";", ",", x, fixed = TRUE)
}

# The four sets the document marks as redesigned after an external check.  The
# notes are the document's own reasons, translated.
NOTES <- c(
  KASP_AGP9_height_chr19_46.10Mb = paste(
    "redesigned after external check; hairpin with tail was 55.6 C, now 47.6 C.",
    "The table prints A2 = 0 exact genome hits, while the prose says the A2",
    "primer has one unrelated exact match; the contradiction is kept, not resolved"
  ),
  KASP_AGPO2_height_chr19_49.74Mb = paste(
    "redesigned after external check; a perfectly linked SNP (r2 = 1.0) sat",
    "under the allele primer, so the matrices saw two mismatches; HEX Tm 62.7 C",
    "is the remaining concession"
  ),
  KASP_AGPO2_height2_chr19_51.07Mb = paste(
    "redesigned; a polymorphic site sat on the common primer's 3' base, which",
    "blocks extension, so the common primer was shifted"
  ),
  KASP_JOINT_protein_chr15_11.28Mb = paste(
    "redesigned; hairpin with tail is 51.4 C, above the 50 C working threshold;",
    "an amplicon on the same point is kept in reserve"
  ),
  KASP_AGP10_branch_chr6_13.303311Mb = paste(
    "proxy for the peak one-base indel (the indel itself is not usable for an",
    "allele-specific assay); r2 = 0.989 with the peak, 19.8 kb away"
  ),
  KASP_DT2_near_chr18_56.14Mb = paste(
    "GBS calls in this block are unreliable: raw heterozygote fraction 0.353,",
    "smoothed 0.864, while markers just outside the block read 0.47-0.52"
  ),
  AMP_DT2_EXON1_SerAsn = paste(
    "non-synonymous G/A (Ser -> Asn); below the association threshold in the",
    "source, so it is part of a haplotype rather than a stand-alone candidate"
  )
)

FAM_TAIL <- "GAAGGTGACCAAGTTCATGCT"
HEX_TAIL <- "GAAGGTCGGAGTCAACGGATT"
DNA_COMP <- c(A = "T", C = "G", G = "C", T = "A")
primer_strand <- function(body, allele) {
  end <- substr(body, nchar(body), nchar(body))
  if (identical(end, allele)) "same" else if (identical(end, DNA_COMP[[allele]])) "reverse" else NA_character_
}

parse_sanger <- function(id, b) {
  fwd <- caps(b, "прямой(?: праймер)? ([ACGT]{15,})", paste(id, "forward"))[[1L]]
  rev <- caps(b, "обратный(?: праймер)? ([ACGT]{15,})", paste(id, "reverse"))[[1L]]
  tm <- caps(b, "Tm ([0-9]+,[0-9]+) / ([0-9]+,[0-9]+)", paste(id, "Tm"))
  band <- caps(
    b, "(?:ожидаемая(?: полоса)?|полоса) ([0-9][0-9 ]*) п\\. н\\.",
    paste(id, "band size")
  )[[1L]]
  amp <- caps(
    b, "ампликон (?:Gm[0-9]+:)?([0-9][0-9 ]*)[[:space:]]*[\u2013-][[:space:]]*([0-9][0-9 ]*)",
    paste(id, "amplicon interval")
  )
  off <- opt_caps(
    b, "([0-9][0-9 ]*) п\\. н\\. от прямого, ([0-9][0-9 ]*) п\\. н\\. от обратного"
  )
  point <- opt_caps(b, "точка (Gm[0-9]+):([0-9 ]+), аллели ([A-Z/]+)")
  what <- field_line(b, c("что(?: покрывает)?"))
  alleles <- if (!is.null(what)) regmatches(what, regexpr("[ACGT]/[ACGT]", what)) else character()
  proxy <- id == "AMP_AGP10_branch_chr6_proxy"
  chrom <- if (proxy) "Gm06" else "Gm18"
  list(
    id = id,
    scheme = "sanger_amplicon",
    queue = 1L,
    locus = if (proxy) {
      "AGP10_chr6_13.3Mb, branch number (proxy amplicon)"
    } else if (!is.null(what)) {
      translate_locus(what)
    } else if (!is.null(point)) {
      translate_locus(field_line(b, "локус"))
    } else {
      NA_character_
    },
    population = if (!is.null(point)) {
      translate(field_line(b, "популяция"), pop_en)
    } else {
      "AGP10 (SK Sarita x Yunka)"
    },
    target = if (!is.null(point)) {
      paste0(point[[1L]], ":", gsub(" ", "", point[[2L]]))
    } else {
      paste0(chrom, ":", ru_int(amp[[1L]]), "-", ru_int(amp[[2L]]))
    },
    allele1 = if (!is.null(point)) {
      strsplit(point[[3L]], "/", fixed = TRUE)[[1L]][[1L]]
    } else if (length(alleles)) {
      strsplit(alleles, "/", fixed = TRUE)[[1L]][[1L]]
    } else {
      NA_character_
    },
    allele2 = if (!is.null(point)) {
      strsplit(point[[3L]], "/", fixed = TRUE)[[1L]][[2L]]
    } else if (length(alleles)) {
      strsplit(alleles, "/", fixed = TRUE)[[1L]][[2L]]
    } else {
      NA_character_
    },
    product_bp = ru_int(band),
    gap_bp = NA_integer_,
    amplicon_start = ru_int(amp[[1L]]),
    amplicon_end = ru_int(amp[[2L]]),
    snp_offset_forward_bp = if (is.null(off)) NA_integer_ else ru_int(off[[1L]]),
    snp_offset_reverse_bp = if (is.null(off)) NA_integer_ else ru_int(off[[2L]]),
    gbs_a1a1 = NA_integer_, gbs_het = NA_integer_, gbs_a2a2 = NA_integer_,
    het_fraction = NA_real_, het_n = NA_integer_,
    status = NA_character_,
    note = if (id %in% names(NOTES)) NOTES[[id]] else NA_character_,
    tm_forward_c = ru_num(tm[[1L]]),
    tm_reverse_c = ru_num(tm[[2L]]),
    forward = fwd,
    reverse = rev
  )
}

parse_kasp <- function(id, b) {
  allele_lines <- grep("^(FAM|HEX) = аллель", b, value = TRUE)
  if (length(allele_lines) != 2L) stop(id, ": expected FAM and HEX lines")
  parsed <- regmatches(
    allele_lines,
    regexec("^(FAM|HEX)[[:space:]]*=[[:space:]]*аллель[[:space:]]+([ACGT])[[:space:]]+([ACGT]{15,})[[:space:]]*$", allele_lines, perl = TRUE)
  )
  if (any(lengths(parsed) != 4L)) stop(id, ": could not parse an allele primer line")
  fam <- parsed[[which(vapply(parsed, `[[`, character(1L), 2L) == "FAM")]]
  hex <- parsed[[which(vapply(parsed, `[[`, character(1L), 2L) == "HEX")]]
  common <- caps(b, "общий праймер ([ACGT]{15,})", paste(id, "common primer"))[[1L]]
  tm <- caps(
    b, "Tm тел ([0-9]+,[0-9]+) / ([0-9]+,[0-9]+) °C, общий ([0-9]+,[0-9]+)",
    paste(id, "Tm")
  )
  pg <- caps(
    b, "продукт ([0-9][0-9 ]*) п\\. н\\., зазор(?: до общего праймера)? ([0-9][0-9 ]*)",
    paste(id, "product and gap")
  )
  hairpin <- caps(b, "шпилька с хвостом ([0-9]+,[0-9]+) / ([0-9]+,[0-9]+)", paste(id, "hairpin"))
  hits <- opt_caps(b, "совпадений в(?: геноме)? A1 ([0-9]+), A2 ([0-9]+), общий ([0-9]+)")
  gbs <- opt_caps(b, "классы GBS ([0-9]+) : ([0-9]+) : ([0-9]+)")
  het <- opt_caps(b, "доля гетерозигот в GBS ([0-9]+,[0-9]+) при n = ([0-9]+)")
  marker <- opt_caps(b, "маркер (AGP[0-9]+), (Gm[0-9]+:[0-9 ]+)")
  point <- opt_caps(b, "точка (Gm[0-9]+):([0-9 ]+)")
  dist <- opt_caps(b, "положение ([\u2212-]?[0-9]+,[0-9]+) кб")
  status <- opt_caps(b, "статус (.+)")

  if (!is.null(marker)) {
    target <- paste(marker[[1L]], gsub(" ", "", marker[[2L]]))
    locus <- if (id == "KASP_DT2_flank_low_chr18_55.26Mb") {
      "Dt2 flank, 650.6 kb below the gene"
    } else if (id == "KASP_DT2_near_chr18_56.14Mb") {
      "Dt2 flank, 217.2 kb from the gene"
    } else {
      "Dt2 flank, 690.9 kb above the gene"
    }
    population <- "AGP10 (SK Sarita x Yunka)"
  } else {
    target <- paste0(point[[1L]], ":", gsub(" ", "", point[[2L]]))
    locus <- translate_locus(field_line(b, "локус"))
    if (id == "KASP_JOINT_protein_chr15_11.28Mb") {
      locus <- "JOINT_protein_chr15_11.28Mb, second protein peak, joint model B (article 1)"
    }
    population <- translate(field_line(b, "популяция"), pop_en)
  }

  list(
    id = id,
    scheme = if (!is.null(marker)) "kasp_dt2_flank" else "kasp",
    queue = 2L,
    locus = locus,
    population = population,
    target = target,
    allele1 = fam[[3L]],
    allele2 = hex[[3L]],
    product_bp = ru_int(pg[[1L]]),
    gap_bp = ru_int(pg[[2L]]),
    amplicon_start = NA_integer_, amplicon_end = NA_integer_,
    snp_offset_forward_bp = NA_integer_, snp_offset_reverse_bp = NA_integer_,
    gbs_a1a1 = if (is.null(gbs)) NA_integer_ else ru_int(gbs[[1L]]),
    gbs_het = if (is.null(gbs)) NA_integer_ else ru_int(gbs[[2L]]),
    gbs_a2a2 = if (is.null(gbs)) NA_integer_ else ru_int(gbs[[3L]]),
    het_fraction = if (is.null(het)) NA_real_ else ru_num(het[[1L]]),
    het_n = if (is.null(het)) NA_integer_ else ru_int(het[[2L]]),
    status = if (is.null(status)) NA_character_ else "redesigned after external check",
    note = if (id %in% names(NOTES)) NOTES[[id]] else NA_character_,
    tm_allele1_c = ru_num(tm[[1L]]),
    tm_allele2_c = ru_num(tm[[2L]]),
    tm_common_c = ru_num(tm[[3L]]),
    hairpin_allele1_c = ru_num(hairpin[[1L]]),
    hairpin_allele2_c = ru_num(hairpin[[2L]]),
    genome_hits_allele1 = if (is.null(hits)) NA_integer_ else ru_int(hits[[1L]]),
    genome_hits_allele2 = if (is.null(hits)) NA_integer_ else ru_int(hits[[2L]]),
    genome_hits_common = if (is.null(hits)) NA_integer_ else ru_int(hits[[3L]]),
    allele1_sequence = fam[[4L]],
    allele2_sequence = hex[[4L]],
    common = common,
    dt2_distance_kb = if (is.null(dist)) NA_real_ else ru_num(dist[[1L]])
  )
}

parsed <- lapply(names(blocks), function(id) {
  if (startsWith(id, "AMP_")) parse_sanger(id, blocks[[id]]) else parse_kasp(id, blocks[[id]])
})

## --- the reserve Dt2 body table --------------------------------------------
# The PDF's last table is split across two page boundaries, so it is transcribed
# as eight rows rather than regexed.  Every row is still asserted below: the
# printed band size must equal end - start + 1, and the Tm must reproduce.
dt2_body <- tibble::tribble(
  ~id, ~covers, ~forward, ~reverse, ~product_bp, ~start, ~end, ~tm_forward_c, ~tm_reverse_c,
  "DT2_B1_A1", "exons 1,2", "ACATAACCGTACGGAGCGTC", "TGAACGAGTAGGAACCTGCG", 697L, 55909494L, 55910190L, 59.9, 59.8,
  "DT2_B1_A2", "exons 1,2", "TTCAGCTGAAGCGGATCGAG", "ACCAAGAAACCACACGCTGA", 808L, 55910007L, 55910814L, 60.2, 60.1,
  "DT2_B2_A1", "exon 3", "ATATTCCGCTCAGCATAGGC", "TCCACATACTTACCCGGTAGT", 489L, 55914536L, 55915024L, 57.6, 57.6,
  "DT2_B3_A1", "exons 4,5", "CCCGATTGGTCTGATAGCCA", "TGAAGAAGCACTCGTCACCT", 744L, 55915463L, 55916206L, 59.2, 59.0,
  "DT2_B3_A2", "exons 4,5", "TGACAGGGCTAGTTACGTGC", "CCGTTGTGTGGACCTTCAGT", 956L, 55916005L, 55916960L, 59.8, 60.2,
  "DT2_B4_A1", "exons 6,7,8", "CATCTCTGAAGCGCATCCGA", "TCACCTGCCAATAAGCTGCA", 1069L, 55916160L, 55917228L, 60.2, 60.0,
  "DT2_B5_A1", "exon 9", "TGTAGCTCGAATAAGATTACCTCGA", "TAGCACCAACCATCACGTGT", 415L, 55917722L, 55918136L, 59.2, 59.6,
  "DT2_B5_A2", "exon 9", "TGCGCTGCATGTCTGACTAG", "GAAATCCTGCGGTTTCGAGC", 523L, 55917911L, 55918433L, 60.5, 59.9
)

dt2_assays <- tibble::tibble(
  assay_id = dt2_body$id,
  scheme = "sanger_dt2_body",
  queue = 3L,
  locus = paste("Dt2 body,", dt2_body$covers),
  population = "AGP10 (SK Sarita x Yunka)",
  target = paste0("Gm18:", dt2_body$start, "-", dt2_body$end),
  allele1 = NA_character_, allele2 = NA_character_,
  product_bp = dt2_body$product_bp,
  gap_bp = NA_integer_,
  amplicon_start = dt2_body$start, amplicon_end = dt2_body$end,
  snp_offset_forward_bp = NA_integer_, snp_offset_reverse_bp = NA_integer_,
  gbs_a1a1 = NA_integer_, gbs_het = NA_integer_, gbs_a2a2 = NA_integer_,
  het_fraction = NA_real_, het_n = NA_integer_,
  status = NA_character_,
  note = ifelse(dt2_body$id == "DT2_B4_A1",
    "covers exons 6-8, including the six intronic haplotype SNPs", NA_character_
  ),
  tm_forward_c = dt2_body$tm_forward_c, tm_reverse_c = dt2_body$tm_reverse_c,
  forward = dt2_body$forward, reverse = dt2_body$reverse
)

## --- assemble the two tables -----------------------------------------------
assays <- dplyr::bind_rows(
  lapply(parsed, function(p) {
    tibble::tibble(
      assay_id = p$id, scheme = p$scheme, queue = p$queue,
      locus = p$locus, population = p$population, target = p$target,
      allele1 = p$allele1, allele2 = p$allele2,
      product_bp = p$product_bp, gap_bp = p$gap_bp,
      amplicon_start = p$amplicon_start, amplicon_end = p$amplicon_end,
      snp_offset_forward_bp = p$snp_offset_forward_bp,
      snp_offset_reverse_bp = p$snp_offset_reverse_bp,
      gbs_a1a1 = p$gbs_a1a1, gbs_het = p$gbs_het, gbs_a2a2 = p$gbs_a2a2,
      het_fraction = p$het_fraction, het_n = p$het_n,
      status = p$status, note = p$note,
      tm_forward_c = if (!is.null(p$tm_forward_c)) p$tm_forward_c else NA_real_,
      tm_reverse_c = if (!is.null(p$tm_reverse_c)) p$tm_reverse_c else NA_real_,
      dt2_distance_kb = if (!is.null(p$dt2_distance_kb)) p$dt2_distance_kb else NA_real_,
      genome_build = "GCF_000004515.6 (Wm82.a4.v1)"
    )
  }),
  dt2_assays |>
    dplyr::mutate(
      dt2_distance_kb = NA_real_,
      genome_build = "GCF_000004515.6 (Wm82.a4.v1)"
    )
) |>
  dplyr::arrange(.data$queue, .data$scheme, .data$assay_id)

primers <- dplyr::bind_rows(
  lapply(parsed, function(p) {
    if (p$scheme == "sanger_amplicon") {
      tibble::tibble(
        assay_id = p$id,
        role = c("forward", "reverse"),
        allele = NA_character_, strand = NA_character_, tail = NA_character_,
        sequence = c(p$forward, p$reverse),
        body = c(p$forward, p$reverse),
        tm_c = c(p$tm_forward_c, p$tm_reverse_c),
        hairpin_c = NA_real_,
        genome_hits = NA_integer_, genome_hits_printed = NA_integer_,
        note = NA_character_
      )
    } else {
      tibble::tibble(
        assay_id = p$id,
        role = c("allele1", "allele2", "common"),
        allele = c(p$allele1, p$allele2, NA_character_),
        strand = c(
          primer_strand(sub(paste0("^", FAM_TAIL), "", p$allele1_sequence), p$allele1),
          primer_strand(sub(paste0("^", HEX_TAIL), "", p$allele2_sequence), p$allele2),
          NA_character_
        ),
        tail = c("FAM", "HEX", NA_character_),
        sequence = c(p$allele1_sequence, p$allele2_sequence, p$common),
        body = c(
          sub(paste0("^", FAM_TAIL), "", p$allele1_sequence),
          sub(paste0("^", HEX_TAIL), "", p$allele2_sequence),
          p$common
        ),
        tm_c = c(p$tm_allele1_c, p$tm_allele2_c, p$tm_common_c),
        hairpin_c = c(p$hairpin_allele1_c, p$hairpin_allele2_c, NA_real_),
        genome_hits = NA_integer_,
        genome_hits_printed = c(
          p$genome_hits_allele1, p$genome_hits_allele2, p$genome_hits_common
        ),
        note = ifelse(p$scheme == "kasp_dt2_flank", NA_character_, NA_character_)
      )
    }
  }),
  dt2_body |>
    dplyr::mutate(
      role = "forward", allele = NA_character_, strand = NA_character_, tail = NA_character_,
      body = .data$forward, tm_c = .data$tm_forward_c,
      hairpin_c = NA_real_, genome_hits = NA_integer_,
      genome_hits_printed = NA_integer_, note = NA_character_
    ) |>
    dplyr::transmute(
      assay_id = .data$id, role, allele, strand, tail, sequence = .data$forward,
      body, tm_c, hairpin_c, genome_hits, genome_hits_printed, note
    ),
  dt2_body |>
    dplyr::mutate(
      role = "reverse", allele = NA_character_, strand = NA_character_, tail = NA_character_,
      body = .data$reverse, tm_c = .data$tm_reverse_c,
      hairpin_c = NA_real_, genome_hits = NA_integer_,
      genome_hits_printed = NA_integer_, note = NA_character_
    ) |>
    dplyr::transmute(
      assay_id = .data$id, role, allele, strand, tail, sequence = .data$reverse,
      body, tm_c, hairpin_c, genome_hits, genome_hits_printed, note
    )
) |>
  dplyr::arrange(.data$assay_id, factor(.data$role,
    levels = c("forward", "reverse", "allele1", "allele2", "common")
  ))

## --- the exact hits, from the genome (count, chromosome, coordinates) --------
soy <- soybean_genome()
stopifnot(length(genome_chroms(soy)) == 20L)
# One row per exact full-length occurrence, both strands, over the 20
# chromosomes.  `exact_hit_table()` streams the indexed FASTA one chromosome at
# a time and searches a precomputed dictionary per oligo width.  The primer
# column is the row index in `primers`, so the hits join back without a string
# key (two rows may legitimately carry the same body sequence).
hits <- exact_hit_table(
  rlang::set_names(primers$body, as.character(seq_len(nrow(primers)))),
  soy
)
hit_idx <- as.integer(hits$primer)
primers$genome_hits <- tabulate(hit_idx, nbins = nrow(primers))
hits_by_primer <- split(hits, factor(hit_idx, levels = seq_len(nrow(primers))))

hits_field <- function(z, f) {
  if (is.null(z) || !nrow(z)) {
    return(NA_character_)
  }
  f(z)
}
# Where each hit sits, as comma-separated chromosome / coordinate / strand
# lists (NA if none).  For an ALT-allele primer a match elsewhere is exactly
# the off-target we measure.
primers$genome_hit_chrom <- vapply(
  hits_by_primer, hits_field, character(1L),
  f = function(z) paste(unique(z$chrom), collapse = ",")
)
primers$genome_hit_coord <- vapply(
  hits_by_primer, hits_field, character(1L),
  f = function(z) paste0(z$chrom, ":", z$start, "-", z$end, collapse = ",")
)
primers$genome_hit_strand <- vapply(
  hits_by_primer, hits_field, character(1L),
  f = function(z) paste(z$strand, collapse = ",")
)

# Chromosome each assay targets, used below to check that an on-target
# primer's single exact match is really at the intended locus (an ALT-allele
# primer may legitimately sit elsewhere, which is the off-target we measure).
target_chr <- vapply(
  assays$target[match(primers$assay_id, assays$assay_id)],
  chrom_of_target,
  character(1L)
)

## --- assertions -------------------------------------------------------------
if (anyDuplicated(assays$assay_id)) {
  print(assays$assay_id[duplicated(assays$assay_id)])
}
stopifnot(
  nrow(assays) == 36L,
  nrow(primers) == 85L,
  sum(assays$scheme == "sanger_amplicon") == 15L,
  sum(assays$scheme == "kasp") == 10L,
  sum(assays$scheme == "kasp_dt2_flank") == 3L,
  sum(assays$scheme == "sanger_dt2_body") == 8L,
  sum(assays$queue == 1L) == 15L,
  sum(assays$queue == 2L) == 13L,
  sum(assays$queue == 3L) == 8L,
  !anyDuplicated(assays$assay_id),
  !anyDuplicated(paste(primers$assay_id, primers$role)),
  all(grepl("^[ACGT]+$", primers$sequence)),
  all(grepl("^[ACGT]+$", primers$body)),
  all(primers$body == primers$sequence | primers$tail %in% c("FAM", "HEX")),
  # every oligo is unique in the genome: at most one exact full-length match.
  # An ALT-allele primer wants zero (its 3' base is not in the reference); the
  # ten KASP 2.1 sets print 0 for both allele primers' ALT side and the genome
  # reproduces that.  A Dt2-flank allele primer whose allele is the reference
  # base legitimately has one, so only the non-allele roles are pinned to 1.
  all(primers$genome_hits %in% 0:1),
  all(primers$genome_hits[primers$role != "allele2"] == 1L),
  all(is.na(primers$genome_hit_coord) == (primers$genome_hits == 0L)),
  all(is.na(primers$genome_hit_strand) == (primers$genome_hits == 0L)),
  all(primers$genome_hits[
    primers$assay_id == "KASP_AGP9_height_chr19_46.10Mb" &
      primers$role == "allele2"
  ] == 0L),
  # a primer that is meant to be on-target lands on the chromosome the assay
  # targets; an ALT-allele primer may legitimately hit elsewhere (an off-target)
  all(vapply(
    which(primers$genome_hits > 0L & primers$role != "allele2"),
    function(i) identical(primers$genome_hit_chrom[[i]], target_chr[[i]]),
    logical(1L)
  )),
  # the document's printed counts are a claim about this assembly
  all(primers$genome_hits[!is.na(primers$genome_hits_printed)] ==
    primers$genome_hits_printed[!is.na(primers$genome_hits_printed)])
)

# Sanger: the printed interval and the printed band size are the same quantity.
sanger <- assays[assays$scheme %in% c("sanger_amplicon", "sanger_dt2_body"), ]
stopifnot(all(sanger$product_bp == sanger$amplicon_end - sanger$amplicon_start + 1L))

# The 10 peak assays also print the SNP's distance from each primer.  The
# printed convention is: band = off_F + off_R + len(F) + len(R) - 1.  The
# eleventh block (the branch proxy) has no SNP offsets and is excluded.
off <- assays[!is.na(assays$snp_offset_forward_bp), ]
off <- off |>
  dplyr::left_join(
    primers[primers$role == "forward", c("assay_id", "body")],
    by = "assay_id"
  ) |>
  dplyr::left_join(
    primers[primers$role == "reverse", c("assay_id", "body")],
    by = "assay_id", suffix = c("_f", "_r")
  )
stopifnot(nrow(off) == 10L)
stopifnot(all(
  off$product_bp ==
    off$snp_offset_forward_bp + off$snp_offset_reverse_bp +
      nchar(off$body_f) + nchar(off$body_r) - 1L
))

# KASP: product = gap + allele body + common primer, and the allele bodies
# differ only at the 3' base that names the allele.
kasp_ids <- assays$assay_id[assays$scheme %in% c("kasp", "kasp_dt2_flank")]
for (id in kasp_ids) {
  p <- primers[primers$assay_id == id, ]
  a1 <- p[p$role == "allele1", ]
  a2 <- p[p$role == "allele2", ]
  com <- p[p$role == "common", ]
  asy <- assays[assays$assay_id == id, ]
  stopifnot(
    nrow(a1) == 1L, nrow(a2) == 1L, nrow(com) == 1L,
    startsWith(a1$sequence, FAM_TAIL),
    startsWith(a2$sequence, HEX_TAIL),
    substr(a1$body, nchar(a1$body), nchar(a1$body)) %in% c(a1$allele, DNA_COMP[[a1$allele]]),
    substr(a2$body, nchar(a2$body), nchar(a2$body)) %in% c(a2$allele, DNA_COMP[[a2$allele]]),
    a1$strand == primer_strand(a1$body, a1$allele),
    a2$strand == primer_strand(a2$body, a2$allele),
    substr(a1$body, 1L, nchar(a1$body) - 1L) ==
      substr(a2$body, 1L, nchar(a2$body) - 1L),
    asy$product_bp == asy$gap_bp + nchar(a1$body) + nchar(com$body)
  )
}

## --- the 13 KASP sets, graded against the package's own KASP gates -----------
# The quality gates are not read from the document.  Each of the 13 sets is
# assembled into a minimal `marker_assay` from the panel's primers plus the
# assembly, and `marker_design_checks()` - the same criterion registry the
# designer uses - grades it.  The AS strand and the REF/ALT orientation come
# from the allele-1 primer's exact hit and the reference base at the target; the
# on/off targets are the AS binding window carrying REF vs ALT, so
# `pcr_efficiency()` computes the discrimination the designer would have seen;
# in-silico PCR is the exact REF/ALT product count.
#
# Two honest exclusions:
#  * the panel carries no population variant table, so the footprint-count and
#    -dropout checks and the two linkage checks are reported `unknown`, never
#    passed;
#  * Tm is read with the `SantaLucia1998-2` salt correction, the model that
#    reproduces the printed Tm (bullet (c)).  The package default (Schildkraut
#    2010) reads ~3.7 C low and would fail every COM Tm as a model offset, not
#    as a design property.
kasp_sets <- assays[assays$scheme %in% c("kasp", "kasp_dt2_flank"), ]
stopifnot(nrow(kasp_sets) == 13L)
target_pos <- function(x) as.integer(sub("^.*:", "", x))
kasp_chr <- list()
kasp_chrom_seq <- function(ch) {
  if (is.null(kasp_chr[[ch]])) kasp_chr[[ch]] <<- genome_seq(soy, ch)
  Biostrings::DNAString(kasp_chr[[ch]])
}
kasp_base <- function(ch, pos) {
  as.character(Biostrings::subseq(kasp_chrom_seq(ch), pos, pos))
}
kasp_seg <- function(ch, start, end) {
  as.character(Biostrings::subseq(kasp_chrom_seq(ch), start, end))
}
# checks that need a population variant table the panel does not carry
kasp_no_pop <- c(
  "linked variants skipped (n)", "design relies on linkage",
  "AS 3' window variants (n)", "AS footprint variants (n)",
  "COM footprint variants (n)", "AS 3' window dropout",
  "AS footprint dropout", "COM footprint dropout"
)
kasp_gates <- vector("list", nrow(kasp_sets))
for (i in seq_len(nrow(kasp_sets))) {
  id <- kasp_sets$assay_id[[i]]
  p <- primers[primers$assay_id == id, ]
  a1 <- p[p$role == "allele1", ]
  a2 <- p[p$role == "allele2", ]
  com <- p[p$role == "common", ]
  stopifnot(nrow(a1) == 1L, nrow(a2) == 1L, nrow(com) == 1L)
  h1 <- hits[
    hits$primer == as.character(which(primers$assay_id == id & primers$role == "allele1")),
  ]
  h2 <- hits[
    hits$primer == as.character(which(primers$assay_id == id & primers$role == "allele2")),
  ]
  hc <- hits[
    hits$primer == as.character(which(primers$assay_id == id & primers$role == "common")),
  ]
  stopifnot(nrow(h1) == 1L, nrow(hc) == 1L)
  chrom <- chrom_of_target(kasp_sets$target[[i]])
  pos <- target_pos(kasp_sets$target[[i]])
  strand <- h1$strand[[1L]]
  anchor <- if (strand == "+") h1$end[[1L]] else h1$start[[1L]]
  stopifnot(anchor == pos)
  ref <- kasp_base(chrom, pos)
  # a primer's 3' base pairs with the complementary template strand, so the
  # panel's `allele` is already the plus-strand base
  stopifnot(ref == a1$allele, a2$allele != ref)
  alt <- a2$allele
  seg <- kasp_seg(chrom, h1$start[[1L]], h1$end[[1L]])
  anchor_seg <- if (strand == "+") nchar(seg) else 1L
  stopifnot(substr(seg, anchor_seg, anchor_seg) == ref)
  seg_alt <- seg
  substr(seg_alt, anchor_seg, anchor_seg) <- alt
  t_ref <- as.character(as_target(Biostrings::DNAString(seg), strand))
  t_alt <- as.character(as_target(Biostrings::DNAString(seg_alt), strand))
  eff <- pcr_efficiency(
    c(a1$body, a1$body, a2$body, a2$body),
    c(t_ref, t_alt, t_alt, t_ref)
  )
  # the REF primer and the COM primer sit on one chromosome, opposite strands;
  # the span between them is the panel's printed product (asserted)
  span <- max(h1$end[[1L]], hc$end[[1L]]) - min(h1$start[[1L]], hc$start[[1L]]) + 1L
  stopifnot(span == kasp_sets$product_bp[[i]])
  # the ALT primer gives an exact product only if it and COM share a chromosome
  alt_products <- if (nrow(h2)) {
    sum(vapply(seq_len(nrow(h2)), function(j) {
      any(hc$chrom == h2$chrom[[j]] & hc$strand != h2$strand[[j]] &
        abs(hc$start - h2$start[[j]]) < 1000L)
    }, logical(1L)))
  } else {
    0L
  }
  stopifnot(alt_products == 0L)
  primers3 <- c(AS_ref = a1$body, AS_alt = a2$body, COM = com$body)
  conditions <- scheme_conditions("snp")
  thermo <- purrr::imap(primers3, function(s, nm) {
    dplyr::mutate(
      primer_thermo(s), name = nm,
      tm = primer_tm(s, salt_method = "SantaLucia1998-2")
    )
  }) |>
    purrr::list_rbind()
  as_gc <- gc_frac(a1$body)
  com_gc <- gc_frac(com$body)
  a1_body <- a1$body
  a2_body <- a2$body
  com_body <- com$body
  assay_row <- tibble::tibble(
    amplicon = span,
    com_tm = thermo$tm[thermo$name == "COM"],
    as_ref = a1_body, as_alt = a2_body, com = com_body,
    as_gc = as_gc, com_gc = com_gc,
    as_n_3p = NA_integer_, as_dropout = NA_real_,
    as_n_var = NA_integer_, as_foot_dropout = NA_real_,
    n_var = NA_integer_, com_dropout = NA_real_
  )
  ideogram_hits <- dplyr::bind_rows(
    if (nrow(h1)) tibble::tibble(primer = "AS_ref", chrom = h1$chrom, pos = h1$start),
    if (nrow(h2)) tibble::tibble(primer = "AS_alt", chrom = h2$chrom, pos = h2$start),
    if (nrow(hc)) tibble::tibble(primer = "COM", chrom = hc$chrom, pos = hc$start)
  )
  assay <- list(
    scheme = "snp",
    snp = tibble::tibble(chrom = chrom, pos = pos, ref = ref, alt = alt),
    assay = assay_row,
    efficiency = tibble::tibble(
      ref_on = eff[[1L]], ref_off = eff[[2L]],
      alt_on = eff[[3L]], alt_off = eff[[4L]],
      discrimination = ((eff[[1L]] - eff[[2L]]) + (eff[[3L]] - eff[[4L]])) / 2
    ),
    thermo = thermo,
    complementarity = dimer_thermo(primers3, conditions = conditions),
    amplicons = tibble::tibble(
      primer = c("AS_ref", "AS_alt"),
      n_products = c(1L, alt_products),
      width = c(span, NA_integer_)
    ),
    ideogram = list(hits = ideogram_hits),
    site = marker_site_class(ref, alt)
  )
  assay$site$alt_covered <- 1L
  assay <- structure(assay, class = "marker_assay")
  ev <- marker_design_checks(assay)
  # the panel carries no population variant table: report those checks unknown
  unknown_idx <- ev$criterion %in% kasp_no_pop
  ev$observed[unknown_idx] <- NA_real_
  ev$pass[unknown_idx] <- NA
  ev$status[unknown_idx] <- "unknown"
  ev$detail[unknown_idx] <- "not evaluated: the panel carries no population variant table"
  ev$detail[ev$criterion == "in-silico PCR (REF 1 / ALT 0)"] <-
    "REF pair gives one exact product; ALT pair gives none"
  ev$detail[ev$criterion == "off-target exact primer hits (n)"] <-
    "every primer has at most one exact genome hit"
  ev$detail[ev$criterion == "discrimination"] <-
    "DECIPHER Taq efficiency over the AS binding window (REF vs ALT)"
  ev$assay_id <- id
  kasp_gates[[i]] <- ev
}
kasp_gates <- dplyr::bind_rows(kasp_gates) |>
  dplyr::mutate(observed = round(.data$observed, 3)) |>
  dplyr::transmute(
    assay_id, stage, criterion, observed, limit, status, pass, detail, level
  )
stopifnot(
  nrow(kasp_gates) == 13L * sum(marker_criteria(scheme = "snp")$role == "check"),
  setequal(unique(kasp_gates$assay_id), kasp_sets$assay_id),
  all(kasp_gates$status %in% c("pass", "fail", "unknown")),
  # the variant/linkage checks are unknown, never silently passed
  all(kasp_gates$status[kasp_gates$criterion %in% kasp_no_pop] == "unknown"),
  # the genome-level gates are evaluated
  all(kasp_gates$status[kasp_gates$criterion == "off-target exact primer hits (n)"] == "pass"),
  all(kasp_gates$status[kasp_gates$criterion == "in-silico PCR (REF 1 / ALT 0)"] == "pass")
)

## --- the 23 Sanger products, from the genome ---------------------------------
# A forward primer's exact hit is on the plus strand and a reverse primer's on
# the minus strand, so a product needs both on one chromosome with the forward
# upstream.  Every Sanger primer is already asserted unique genome-wide above,
# hence at most one exact product can exist anywhere.  in_silico_pcr() then
# confirms that product on its target chromosome; scanning the whole genome
# again per pair would repeat the uniqueness result at 20x the cost.
sanger_idx <- which(assays$scheme %in% c("sanger_amplicon", "sanger_dt2_body"))
stopifnot(length(sanger_idx) == 23L)
genome_products <- rep(NA_integer_, nrow(assays))
genome_product_coord <- rep(NA_character_, nrow(assays))
chr_cache <- list()
for (i in sanger_idx) {
  id <- assays$assay_id[[i]]
  f_i <- which(primers$assay_id == id & primers$role == "forward")
  r_i <- which(primers$assay_id == id & primers$role == "reverse")
  stopifnot(length(f_i) == 1L, length(r_i) == 1L)
  chrom <- chrom_of_target(assays$target[[i]])
  stopifnot(!is.na(chrom))
  fh <- hits[hits$primer == as.character(f_i), ]
  rh <- hits[hits$primer == as.character(r_i), ]
  stopifnot(
    nrow(fh) == 1L, nrow(rh) == 1L,
    fh$chrom == chrom, rh$chrom == chrom,
    fh$strand == "+", rh$strand == "-"
  )
  start <- fh$start[[1L]]
  end <- rh$end[[1L]]
  stopifnot(
    start == assays$amplicon_start[[i]],
    end == assays$amplicon_end[[i]],
    end - start + 1L == assays$product_bp[[i]]
  )
  # the independent confirmation: the package's in-silico PCR on that chromosome
  if (is.null(chr_cache[[chrom]])) {
    chr_cache[[chrom]] <- Biostrings::DNAStringSet(
      rlang::set_names(list(genome_seq(soy, chrom)), chrom)
    )
  }
  prod <- suppressWarnings(in_silico_pcr(
    primers$body[[f_i]], primers$body[[r_i]], chr_cache[[chrom]],
    max_mismatch = 0L
  ))
  prod <- prod[prod$n_mismatch == 0L, ]
  stopifnot(
    nrow(prod) == 1L,
    prod$chrom == chrom, prod$start == start, prod$end == end,
    prod$width == assays$product_bp[[i]]
  )
  genome_products[[i]] <- 1L
  genome_product_coord[[i]] <- paste0(chrom, ":", start, "-", end)
}
assays$genome_products <- genome_products
assays$genome_product_coord <- genome_product_coord
stopifnot(all(assays$genome_products[sanger_idx] == 1L))

# Tm: the panel's printed values pin the salt correction.  The default
# (Schildkraut 2010) is a level shift, not a per-oligo error; the matched model
# must reproduce every printed Tm to better than 0.15 C.
tm_calc <- primer_tm(primers$body, salt_method = "SantaLucia1998-2")
tm_delta <- tm_calc - primers$tm_c
stopifnot(max(abs(tm_delta)) <= 0.15)
tm_default <- primer_tm(primers$body)
stopifnot(mean(tm_default - primers$tm_c) < -2)

# Exactly the four sets the document marks as redesigned.
stopifnot(
  sum(!is.na(assays$status)) == 4L,
  setequal(
    assays$assay_id[!is.na(assays$status)],
    c(
      "KASP_AGP9_height_chr19_46.10Mb",
      "KASP_AGPO2_height_chr19_49.74Mb",
      "KASP_AGPO2_height2_chr19_51.07Mb",
      "KASP_JOINT_protein_chr15_11.28Mb"
    )
  )
)

assays_out <- assays |>
  dplyr::select(
    "assay_id", "scheme", "queue", "locus", "population", "target",
    "allele1", "allele2", "product_bp", "gap_bp",
    "amplicon_start", "amplicon_end",
    "snp_offset_forward_bp", "snp_offset_reverse_bp",
    "gbs_a1a1", "gbs_het", "gbs_a2a2", "het_fraction", "het_n",
    "status", "note", "dt2_distance_kb",
    "genome_products", "genome_product_coord", "genome_build"
  )
primers_out <- primers |>
  dplyr::select(
    "assay_id", "role", "allele", "strand", "tail", "sequence", "body",
    "tm_c", "hairpin_c", "genome_hits", "genome_hits_printed",
    "genome_hit_chrom", "genome_hit_coord", "genome_hit_strand", "note"
  )

write_checked <- function(tbl, path) {
  written <- utils::capture.output(
    readr::write_csv(tbl, stdout(), na = ""),
    type = "output"
  )
  if (file.exists(path)) {
    on_disk <- readLines(path, warn = FALSE)
    if (!identical(on_disk, written)) {
      stop("the CSV on disk disagrees with the recomputed table: ", path)
    }
  }
  readr::write_csv(tbl, path, na = "")
  invisible(path)
}

out_dir <- fs::path(root, "inst", "extdata")
write_checked(assays_out, fs::path(out_dir, "soybean_worked_panel_assays.csv"))
write_checked(primers_out, fs::path(out_dir, "soybean_worked_panel_primers.csv"))
write_checked(kasp_gates, fs::path(out_dir, "soybean_worked_panel_kasp_gates.csv"))

message(
  "soybean_worked_panel: ", nrow(assays_out), " assays, ", nrow(primers_out),
  " oligos, ", nrow(kasp_gates), " KASP gate rows; ",
  "Tm matched to ", sprintf("%.3f", max(abs(tm_delta))), " C; ",
  "default-salt offset ", sprintf("%.2f", mean(tm_default - primers$tm_c)), " C"
)
