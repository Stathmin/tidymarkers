# =============================================================================
#  dev/primer3-crosscheck.R  --  tidymarkers against Primer3 2.6.1, end to end
# =============================================================================
#  Primer3 is the reference implementation the package's thermodynamics comes
#  from (the same SantaLucia nearest-neighbour family) and the field-standard
#  primer designer.  This script measures how far apart the two actually are,
#  on one fixed machine, and writes the numbers down.  Nothing here is a
#  tolerance: it reports, it does not "correct".
#
#  PART 1 - thermodynamics, same oligo set, same conditions
#    * `primer_tm()` (TmCalculator, SantaLucia 2004 NN) vs `oligotm` at the
#      Primer3 2.6.1 defaults it runs with when no flags are given: `-tp 1`
#      (SantaLucia 1998 NN table) and `-sc 1` (SantaLucia 1998 salt
#      correction).  `primer_tm()`'s default salt method is instead
#      "Schildkraut2010" (16.6*log10[Na+]), which is `oligotm -sc 0`
#      (Schildkraut & Lifson 1965), so the residual offset is the NN TABLE
#      (2004 vs 1998) PLUS the salt-formula difference and the two cannot be
#      separated on this panel alone.  A second column calls the package at
#      Mg = 0 with the same Na: the shift between the two columns is what the
#      package's divalent-cation term is worth here (~8.5 C), so the ~4 C
#      residual at matched salt is not a missing Mg correction.
#      The salt method is left at each engine's own default on purpose: the
#      question is how far apart the two tools are as shipped.
#    * `oligo_dg_vec()` / `duplex_dg_thermo()` (OligoArrayAux hybrid-ss-min /
#      hybrid-min, 37 C, -N = Na/1000) vs `ntthal` (-a HAIRPIN / -a ANY) at the
#      same -mv/-dv/-n/-d/-t.  These are DIFFERENT salt models (OligoArrayAux
#      takes one monovalent molarity; ntthal derives a [Na+]-equivalent from
#      -mv/-dv/-n).  The gap is reported as an offset and never used to adjust
#      a package number.
#    Conditions matched to `reaction_conditions()`: Na 50 mM, Mg 1.5 mM,
#    dNTPs 0.6 mM, DMSO 0 %, 37 C, plus oligotm -d 50 (nM strand).
#    Oligo panel: the 24 pseudo-random 20-mers of tests/testthat/test-tm-oracle.R
#    (set.seed(20040301L)), the 8 edge probes of helper-tm-tolerance.R
#    (`tm_edge_oligos()`) and the 2 universal KASP tails: 34 in total.  The
#    hetero-dimer set is 4 synthetic pairs, both engines on the same pairs.
#
#  PART 2 - design on 20 real Arabidopsis SNPs
#  Panel: the first records of the cached Ensembl variation JSON
#    (marker_cache_dir()/arabidopsis_thaliana_var_4:1500000-1560000.json) that
#    are SNVs, have exactly two single-base alleles, and whose first allele
#    equals the TAIR10 base in `ath_genome()`; ALT = the other allele.
#    The loci are thinned to >= 1 kb apart because the whole 60 kb cache region
#    is dense (~10 bp between SNVs): `design_kasp_markers()` passes `loci` as
#    its own `variants` table for footprint accounting, so clustered loci would
#    count each other as footprint variants and every assay would look doomed.
#    Spacing, not a change to `variants`, is the honest fix - a genome-wide
#    locus table has the same property.
#    Settings for both designers (a KASP amplicon, so short):
#      Primer3 2.6.1, boulder-IO, per locus, run TWICE:
#        SEQUENCE_TEMPLATE          300 bp, SNP-149 .. SNP+150
#        SEQUENCE_TARGET            SNP at 1-based 150, width 1
#        PRIMER_OPT_SIZE 22 / MIN 18 / MAX 25
#        PRIMER_OPT_TM 60 / MIN 57 / MAX 63
#        PRIMER_PRODUCT_SIZE_RANGE  60-150  ("wide", Primer3's usual range)
#                           and  50-100  ("kasp", the package's own window)
#        PRIMER_NUM_RETURN 5, PRIMER_PICK_ANYWAY 0, everything else default
#      The head-to-head gate table uses the 50-100 run, so "the amplicon is
#      longer than tidymarkers allows" cannot be the explanation for a
#      disagreement.  The wide run is reported only as "how often Primer3 would
#      have solved the locus with a product the package rejects on length".
#    tidymarkers `design_kasp_markers(genome, loci)` with its own defaults:
#      flank 150 (the same 300 bp window), KASP_LIMITS: amplicon 50-100 bp,
#      GC 0.40-0.60, AS Tm 55-65 C, COM Tm 58-63 C, max_run 4, hairpin dG >= -3,
#      cross/homodimer dG >= -9, discrimination >= 0.5.
#    The comparison is restricted to the "any primer pair" question: does a
#    tidymarkers hard gate reject a pair Primer3 accepted?  Primer3 does not
#    design allele-specific primers and this script does not claim it does.
#
#  PART 3 - tests/testthat/test-primer3-crosscheck.R asserts only what this
#  run supports (Tm offset band, rank correlation, the dG offset band, the
#  design-agreement floors).  Re-run this script to see them move.
#
#  Run:  Rscript dev/primer3-crosscheck.R
#  Needs: oligotm, ntthal, primer3_core on PATH (Primer3 2.6.1), and the
#  OligoArrayAux build in .tools/.  Results -> results/ (git-ignored).
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

needed <- c("oligotm", "ntthal", "primer3_core")
have <- Sys.which(needed)
if (any(!nzchar(have))) {
  stop("missing binary/binaries: ", paste(needed[!nzchar(have)], collapse = ", "))
}
cat(sprintf("binaries: %s\n", paste(needed, unname(have), sep = " = ", collapse = ", ")))

suppressMessages(pkgload::load_all(root, quiet = TRUE))
res_dir <- fs::path(root, "results")
fs::dir_create(res_dir)

COND <- reaction_conditions() # Na 50 mM, Mg 1.5 mM, dNTPs 0.6 mM, DMSO 0, 37 C
if (!has_oligoarrayaux()) stop("OligoArrayAux (hybrid-min) not found; the dG half cannot run")

fmt <- function(x, d = 2L) formatC(x, format = "f", digits = d)
reg <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) {
    return(c(r = NA_real_, rho = NA_real_))
  }
  c(r = stats::cor(x[ok], y[ok]), rho = stats::cor(x[ok], y[ok], method = "spearman"))
}
dg_stats <- function(pkg, p3, label) {
  keep <- is.finite(pkg) & is.finite(p3)
  pkg <- pkg[keep]
  p3 <- p3[keep]
  dd <- p3 - pkg
  n <- length(dd)
  list(
    n = n, label = label,
    pkg_mean = mean(pkg), p3_mean = mean(p3),
    delta_mean = mean(dd), delta_sd = stats::sd(dd),
    delta_min = min(dd), delta_max = max(dd),
    r = if (n >= 3L) stats::cor(pkg, p3) else NA_real_,
    rho = if (n >= 3L) stats::cor(pkg, p3, method = "spearman") else NA_real_,
    n_offset_lt_2 = sum(abs(dd) < 2)
  )
}

# --- PART 1: the panel ------------------------------------------------------
# The 24-mers of the oracle test, verbatim (same seed, same call).
set.seed(20040301L)
random20 <- vapply(seq_len(24L), function(i) {
  paste(sample(c("A", "C", "G", "T"), 20L, replace = TRUE), collapse = "")
}, character(1))

# The edge probes, copied by value from tests/testthat/helper-tm-tolerance.R
# (`tm_edge_oligos()`), because this script does not source test helpers.
edge <- c(
  at_rich = "AATATTACCGTCGATACTAT",
  gc_rich = "GCGCGCGCGGGCGCGGCGCG",
  homo_a = "GATGATAAAAAAAAATGATC",
  homo_g = "CATCATGGGGGGGGATCATG",
  three_prime_a = "ATTAAACAACACACAAAAAA",
  three_prime_g = "GGGGGGGGATGCATGCATCA",
  homo5 = "CAAGCAGCAAAAACAGACAAC",
  poly_a = "AAAAAAAAAAAAAAAAAAAA"
)

# The two universal KASP tails (R/marker_core.R).
kasptails <- c(kasp_tail_fam = KASP_TAIL_FAM, kasp_tail_hex = KASP_TAIL_HEX)

# --- the engine wrappers ----------------------------------------------------
# oligotm with the package's conditions and Primer3's defaults for table/salt.
oligotm_tm <- function(seqs, na = COND$Na, mg = COND$Mg, dntp = COND$dNTPs,
                       dmso = COND$DMSO) {
  # `-dm` is only passed when DMSO is actually present: oligotm 2.6.1 returns
  # -273.15 when it is given an explicit `-dm 0`, while its own default (also
  # 0 %) works.  The package default is DMSO = 0, so the run stays on the
  # documented default rather than the broken flag.
  dmso_flag <- if (is.finite(dmso) && dmso > 0) c("-dm", as.character(dmso)) else character(0)
  vapply(seqs, function(s) {
    out <- suppressWarnings(system2("oligotm", c(
      "-mv", as.character(na), "-dv", as.character(mg),
      "-n", as.character(dntp), "-d", "50", dmso_flag, s
    ), stdout = TRUE, stderr = TRUE))
    as.numeric(trimws(out[[1L]]))
  }, numeric(1))
}

ntthal_lines <- function(a, s1, s2 = NULL, na = COND$Na, mg = COND$Mg,
                         dntp = COND$dNTPs, temp = COND$temperature) {
  out <- suppressWarnings(system2("ntthal", c(
    "-a", a, "-mv", as.character(na), "-dv", as.character(mg),
    "-n", as.character(dntp), "-d", "50", "-t", as.character(temp),
    "-s1", s1, if (!is.null(s2)) c("-s2", s2)
  ), stdout = TRUE, stderr = TRUE))
  paste(out, collapse = "\n")
}

first_field <- function(txt, key) {
  m <- regmatches(txt, regexpr(paste0(key, " = [^ \t\n]+"), txt))
  if (!length(m)) NA_real_ else as.numeric(sub(paste0(key, " = "), "", m[[1L]]))
}

# ntthal prints "No secondary structure" for an oligo it cannot fold -> 0.0
# kcal/mol (no structure), counted separately.
parse_ntthal_dg <- function(txt) {
  if (grepl("No secondary structure", txt)) {
    return(c(dg = 0, tm = NA_real_, none = 1))
  }
  c(dg = first_field(txt, "dG") / 1000, tm = first_field(txt, "t"), none = 0)
}
ntthal_hairpin <- function(s) parse_ntthal_dg(ntthal_lines("HAIRPIN", s))

# Guard: prove the dG parser.  A perfectly matched decamer duplex has dG = dH -
# T*dS by definition, so the parsed dG must equal that identity.
self_check_ntthal <- function() {
  txt <- ntthal_lines("ANY", "AAAAAAAAAA", "TTTTTTTTTT", na = 50, mg = 0, dntp = 0)
  t <- 273.15 + COND$temperature
  c(raw = first_field(txt, "dG") / 1000,
    from_enthalpy = (first_field(txt, "dH") - t * first_field(txt, "dS")) / 1000)
}

build_rows <- function(panel, kind) {
  nm <- names(panel)
  hp <- t(vapply(unname(panel), ntthal_hairpin, numeric(3L)))
  tibble::tibble(
    oligo = nm, seq = unname(panel), kind = kind, len = nchar(unname(panel)),
    pkg_tm = rlang::set_names(primer_tm(unname(panel), conditions = COND), nm),
    # the same calculation with no Mg: the shift against pkg_tm is what the
    # package's divalent term is worth, so the residual at matched salt is not
    # a missing Mg correction
    pkg_tm_nomg = rlang::set_names(primer_tm(unname(panel), conditions = reaction_conditions(
      Na = COND$Na, Mg = 0, dNTPs = 0, DMSO = 0
    )), nm),
    p3_tm = rlang::set_names(oligotm_tm(unname(panel)), nm),
    pkg_hairpin = rlang::set_names(hairpin_dg_thermo(unname(panel), conditions = COND), nm),
    p3_hairpin = hp[, "dg"],
    p3_hairpin_tm = hp[, "tm"],
    p3_hairpin_none = hp[, "none"] > 0
  )
}

# hybrid-min vs ntthal -a ANY, matched to the package's duplex_dg_thermo() call
# (na = COND$Na / 1000 molar, T = COND$temperature).
hetero_rows <- function(nm, a, b) {
  tibble::tibble(
    pair = nm, a = a, b = b,
    pkg = unname(duplex_dg_thermo(a, b, conditions = COND)),
    p3 = vapply(seq_along(a), function(i) {
      parse_ntthal_dg(ntthal_lines("ANY", a[[i]], b[[i]]))[["dg"]]
    }, numeric(1))
  )
}

hetero_synth <- list(
  # the universal tail against its reverse complement (the strongest self-pair
  # a tailed KASP primer can form)
  tail_vs_tail_rc = c(KASP_TAIL_FAM, KASP_TAIL_HEX[1:1]),
  two_random = c(random20[[1L]], random20[[2L]]),
  random_vs_rc = c(random20[[3L]],
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(random20[[3L]])))),
  edge_pair = c(edge[["homo_g"]], edge[["homo_a"]])
)
hetero_synth[[1L]] <- c(KASP_TAIL_FAM,
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(KASP_TAIL_FAM))))

synth_dimer <- hetero_rows(
  names(hetero_synth),
  vapply(hetero_synth, `[[`, "", 1L), vapply(hetero_synth, `[[`, "", 2L)
)

selfcheck <- self_check_ntthal()

report1 <- function(tab, hetero, label) {
  bad <- tab$oligo[is.finite(tab$pkg_tm) & (tab$pkg_tm > 90 | tab$pkg_tm < 20)]
  u <- tab[!tab$oligo %in% bad, ]
  ok <- is.finite(u$pkg_tm) & is.finite(u$p3_tm)
  u <- u[ok, ]
  delta <- u$p3_tm - u$pkg_tm
  rr <- reg(u$pkg_tm, u$p3_tm)
  d_na <- u$p3_tm - u$pkg_tm_nomg
  rr_na <- reg(u$pkg_tm_nomg, u$p3_tm)
  hp <- dg_stats(tab$pkg_hairpin, tab$p3_hairpin, "hairpin")
  dim_ <- dg_stats(hetero$pkg, hetero$p3, "heterodimer")

  cat(sprintf("\n===== PART 1 (%s) =====\n", label))
  cat(sprintf("panel: %d oligos, length %d-%d nt; Tm stats on n = %d\n",
    nrow(tab), min(tab$len), max(tab$len), nrow(u)))
  if (length(bad)) cat("  excluded (non-physiological pkg Tm):", paste(bad, collapse = ", "), "\n")
  cat(sprintf(
    "Tm  oligotm - primer_tm  : mean %s C, sd %s C, range [%s, %s], max|d| %s; r = %s, rho = %s\n",
    fmt(mean(delta)), fmt(stats::sd(delta)), fmt(min(delta)), fmt(max(delta)), fmt(max(abs(delta))),
    fmt(rr[["r"]], 3L), fmt(rr[["rho"]], 3L)
  ))
  cat(sprintf(
    "    same with pkg at Mg = 0: mean %s C, sd %s C, max|d| %s, r = %s\n",
    fmt(mean(d_na)), fmt(stats::sd(d_na)), fmt(max(abs(d_na))), fmt(rr_na[["r"]], 3L)
  ))
  cat(sprintf(
    "hairpin dG ntthal - hybrid-ss-min: n = %s, mean %s, sd %s, range [%s, %s]; rho = %s; %s/%s within 2\n",
    as.integer(hp[["n"]]), fmt(hp[["delta_mean"]]), fmt(hp[["delta_sd"]]),
    fmt(hp[["delta_min"]]), fmt(hp[["delta_max"]]), fmt(hp[["rho"]], 3L),
    as.integer(hp[["n_offset_lt_2"]]), as.integer(hp[["n"]])
  ))
  cat(sprintf(
    "  ntthal found no hairpin on %d/%d (0.0); pkg < -0.5 on %d of those\n",
    sum(tab$p3_hairpin_none), nrow(tab),
    sum(tab$p3_hairpin_none & tab$pkg_hairpin < -0.5)
  ))
  cat(sprintf(
    "  package reports dG = 0 (no foldable structure) on %d/%d\n",
    sum(tab$pkg_hairpin == 0), nrow(tab)
  ))
  cat(sprintf(
    "heterodimer dG ntthal - hybrid-min: n = %s pairs, mean offset %s, sd %s, range [%s, %s], rho = %s\n",
    as.integer(dim_[["n"]]), fmt(dim_[["delta_mean"]]), fmt(dim_[["delta_sd"]]),
    fmt(dim_[["delta_min"]]), fmt(dim_[["delta_max"]]), fmt(dim_[["rho"]], 3L)
  ))
  cat(sprintf(
    "  parser self-check (A10/T10 duplex): dG field %s vs dH - T dS %s kcal/mol\n",
    fmt(selfcheck[["raw"]], 3L), fmt(selfcheck[["from_enthalpy"]], 3L)
  ))
  invisible(list(
    n_tm = nrow(u), tm_mean = mean(delta), tm_sd = stats::sd(delta),
    tm_min = min(delta), tm_max = max(delta), tm_absmax = max(abs(delta)),
    r = rr[["r"]], rho = rr[["rho"]],
    tm_mean_nomg = mean(d_na), tm_sd_nomg = stats::sd(d_na),
    hp = hp, dimer = dim_, excluded = bad
  ))
}

random20_named <- rlang::set_names(random20, paste0("oligo", seq_along(random20)))
tab_panel <- build_rows(c(random20_named, edge, kasptails), "panel") |>
  dplyr::mutate(kind = dplyr::case_when(
    grepl("^oligo", .data$oligo) ~ "random20",
    .data$oligo %in% names(edge) ~ "edge",
    TRUE ~ "kasp_tail"
  ))
stat_panel <- report1(tab_panel, synth_dimer, "panel only")

# =============================================================================
#  PART 2 - design: 20 real SNPs, Primer3 vs design_kasp_markers()
# =============================================================================
cat("\n===== PART 2: design on a shared SNP panel =====\n")

genome <- ath_genome()
json <- fs::path(marker_cache_dir(), "arabidopsis_thaliana_var_4:1500000-1560000.json")
if (!fs::file_exists(json)) stop("missing variation cache: ", json)
recs <- jsonlite::fromJSON(json, simplifyVector = FALSE)

snv <- purrr::keep(recs, function(r) {
  !is.null(r$start) && !is.null(r$end) && r$end == r$start &&
    length(r$alleles) == 2L && all(nchar(unlist(r$alleles)) == 1L) &&
    all(toupper(unlist(r$alleles)) %in% c("A", "C", "G", "T"))
})
panel_rec <- purrr::keep(snv, function(r) {
  ch <- paste0("Chr", r$seq_region_name)
  if (!ch %in% names(genome)) return(FALSE)
  base <- as.character(Biostrings::subseq(genome[[ch]], r$start, r$start))
  toupper(base) == toupper(r$alleles[[1L]])
})
cat(sprintf("cache: %d records, %d SNVs, %d with REF == TAIR10\n",
  length(recs), length(snv), length(panel_rec)))

N_LOCI <- as.integer(Sys.getenv("P3_N_LOCI", "20"))
# thin to >= 1 kb apart so no locus falls inside another's 300 bp window
MIN_GAP <- 1000L
kept <- list()
last_pos <- -Inf
for (r in panel_rec) {
  if (length(kept) >= N_LOCI) break
  if (r$start - last_pos >= MIN_GAP) {
    kept[[length(kept) + 1L]] <- r
    last_pos <- r$start
  }
}
panel_rec <- kept
N_LOCI <- length(panel_rec)
panel <- purrr::map_dfr(panel_rec, function(r) {
  al <- toupper(unlist(r$alleles))
  tibble::tibble(
    id = as.character(r$id), chrom = paste0("Chr", r$seq_region_name),
    pos = as.integer(r$start), ref = al[[1L]], alt = al[[2L]]
  )
})
panel$seq <- vapply(seq_len(nrow(panel)), function(i) {
  as.character(Biostrings::subseq(genome[[panel$chrom[[i]]]], panel$pos[[i]] - 149L, panel$pos[[i]] + 150L))
}, character(1))
cat(sprintf("panel: %d loci on %s\n", nrow(panel),
  paste(names(table(panel$chrom)), as.integer(table(panel$chrom)), sep = "=", collapse = ", ")))

# --- Primer3 on each 300 bp window -----------------------------------------
FLANK <- 149L # SNP at 1-based 150 of the 300 bp window
primer3_one <- function(seq, snp_i, id, size_range = "60-150") {
  inp <- c(
    paste0("SEQUENCE_ID=", id),
    paste0("SEQUENCE_TEMPLATE=", seq),
    paste0("SEQUENCE_TARGET=", snp_i, ",1"),
    "PRIMER_TASK=generic",
    "PRIMER_PICK_LEFT_PRIMER=1", "PRIMER_PICK_RIGHT_PRIMER=1",
    "PRIMER_PICK_INTERNAL_OLIGO=0",
    "PRIMER_OPT_SIZE=22", "PRIMER_MIN_SIZE=18", "PRIMER_MAX_SIZE=25",
    "PRIMER_OPT_TM=60", "PRIMER_MIN_TM=57", "PRIMER_MAX_TM=63",
    paste0("PRIMER_PRODUCT_SIZE_RANGE=", size_range),
    "PRIMER_NUM_RETURN=5",
    "="
  )
  out <- suppressWarnings(system2("primer3_core", input = inp, stdout = TRUE, stderr = FALSE))
  kv <- grep("=", out, value = TRUE)
  g <- function(k) {
    hit <- grep(paste0("^", k, "="), kv, value = TRUE)
    if (!length(hit)) NA_character_ else sub("^[^=]*=", "", hit[[1L]])
  }
  n <- suppressWarnings(as.integer(g("PRIMER_PAIR_NUM_RETURNED")))
  n <- if (is.na(n)) 0L else n
  pairs <- if (n > 0L) {
    purrr::map_dfr(seq_len(n), function(k) {
      j <- k - 1L
      tibble::tibble(
        rank = j,
        left = g(paste0("PRIMER_LEFT_", j, "_SEQUENCE")),
        right = g(paste0("PRIMER_RIGHT_", j, "_SEQUENCE")),
        left_pos = as.integer(sub(",.*", "", g(paste0("PRIMER_LEFT_", j, "=")))),
        right_pos = as.integer(sub(",.*", "", g(paste0("PRIMER_RIGHT_", j, "=")))),
        product = as.integer(g(paste0("PRIMER_PAIR_", j, "_PRODUCT_SIZE"))),
        left_tm = as.numeric(g(paste0("PRIMER_LEFT_", j, "_TM"))),
        right_tm = as.numeric(g(paste0("PRIMER_RIGHT_", j, "_TM")))
      )
    })
  } else {
    tibble::tibble(rank = integer(0), left = character(0), right = character(0),
      left_pos = integer(0), right_pos = integer(0), product = integer(0),
      left_tm = numeric(0), right_tm = numeric(0))
  }
  list(error = g("PRIMER_ERROR"), explain = g("PRIMER_PAIR_EXPLAIN"), pairs = pairs)
}

t_p3 <- numeric(nrow(panel))
p3_rows <- vector("list", nrow(panel))
p3_rows_k <- vector("list", nrow(panel))
for (i in seq_len(nrow(panel))) {
  t0 <- proc.time()[["elapsed"]]
  p3_rows[[i]] <- primer3_one(panel$seq[[i]], FLANK + 1L, panel$id[[i]], "60-150")
  p3_rows_k[[i]] <- primer3_one(panel$seq[[i]], FLANK + 1L, panel$id[[i]], "50-100")
  t_p3[[i]] <- proc.time()[["elapsed"]] - t0
}
cat(sprintf("primer3_core: %d loci x 2 runs; >= 1 pair: wide 60-150 %d, kasp 50-100 %d; %s s total (%s s/locus); errors: %s\n",
  nrow(panel),
  sum(vapply(p3_rows, function(x) nrow(x$pairs) > 0L, logical(1L))),
  sum(vapply(p3_rows_k, function(x) nrow(x$pairs) > 0L, logical(1L))),
  fmt(sum(t_p3), 1L), fmt(mean(t_p3), 2L),
  paste(stats::na.omit(vapply(p3_rows, function(x) x$error %||% NA_character_, "")), collapse = " | ")))

# --- tidymarkers on the same loci ------------------------------------------
loci <- panel[c("chrom", "pos", "ref", "alt")]
n_loci <- nrow(loci)
cat(sprintf("design_kasp_markers(): %d loci, single core\n", n_loci))
t_tm0 <- proc.time()[["elapsed"]]
tm_design <- suppressWarnings(design_kasp_markers(genome, loci))
t_tm_design <- proc.time()[["elapsed"]] - t_tm0
cat(sprintf("tidymarkers: %s s total, %s s per locus (n = %d); verdicts: %s\n",
  fmt(t_tm_design, 1L), fmt(t_tm_design / n_loci, 1L), n_loci,
  paste(names(table(tm_design$verdict)), as.integer(table(tm_design$verdict)), sep = "=", collapse = ", ")))

# The assays (cheap: design_kasp_markers() is the expensive pass; this repeats
# it once per locus to keep the primer sequences, which the table does not carry)
tm_assays <- vector("list", n_loci)
tm_checks <- vector("list", n_loci)
for (i in seq_len(n_loci)) {
  snp <- loci[i, ]
  a <- suppressWarnings(tryCatch(assemble_kasp_assay(genome, snp, variants = loci), error = function(e) NULL))
  ok_assay <- !is.null(a) && !is_kasp_doomed(a)
  # `[<-` with list(NULL), not `[[<- NULL`: the latter DELETES the slot and
  # renumbers every later locus, which silently mis-joined the tables
  tm_assays[i] <- list(if (ok_assay) a else NULL)
  tm_checks[i] <- list(if (ok_assay) {
    tryCatch(marker_design_checks(a), error = function(e) NULL)
  } else {
    NULL
  })
}

# --- does a tidymarkers hard gate reject what Primer3 accepted? -------------
# Restricted to the "any primer pair" gates: 3'-end placement off the SNP,
# amplicon length, Tm window, GC, homopolymer, hairpin and hetero-dimer.
gate_tm <- KASP_LIMITS$as_tm
gate_amp <- KASP_LIMITS$amplicon
gate_gc <- KASP_LIMITS$gc * 100
gate_run <- as.integer(KASP_LIMITS$max_run)
gc_frac_of <- function(s) 100 * nchar(gsub("[^GC]", "", s)) / nchar(s)
longest_run_of <- function(s) max(rle(strsplit(s, "")[[1L]])$lengths)

cmp <- purrr::map_dfr(seq_len(nrow(panel)), function(i) {
  # the gate table uses the KASP-window run (50-100 bp), so a length
  # disagreement cannot masquerade as a quality disagreement
  pr <- p3_rows_k[[i]]$pairs
  pr_wide <- p3_rows[[i]]$pairs
  a <- tm_assays[[i]]
  tm_amp <- if (is.null(a)) NA_integer_ else as.integer(a$assay$amplicon[[1L]])
  tm_as <- if (is.null(a)) NA_character_ else a$assay$as_ref[[1L]]
  tm_com <- if (is.null(a)) NA_character_ else a$assay$com[[1L]]
  base <- tibble::tibble(
    id = panel$id[[i]], chrom = panel$chrom[[i]], pos = panel$pos[[i]],
    ref = panel$ref[[i]], alt = panel$alt[[i]],
    tm_verdict = tm_design$verdict[[i]], tm_amplicon = tm_amp,
    tm_strand = tm_design$strand[[i]], tm_discrimination = tm_design$discrimination[[i]],
    tm_n_fail = tm_design$n_fail[[i]], tm_note = tm_design$note[[i]],
    tm_as_ref = tm_as, tm_com = tm_com,
    p3_product_wide = if (nrow(pr_wide)) pr_wide$product[[1L]] else NA_integer_
  )
  if (!nrow(pr)) {
    return(dplyr::bind_cols(base, tibble::tibble(
      p3_n_pairs = 0L, p3_left = NA_character_, p3_right = NA_character_,
      p3_product = NA_integer_, p3_left_tm = NA_real_, p3_right_tm = NA_real_,
      p3_gc_ok = NA, p3_tm_ok = NA, p3_amp_ok = NA, p3_run_ok = NA,
      p3_3p_ok = NA, p3_hairpin_ok = NA, p3_dimer_ok = NA,
      p3_worst_gate = NA_character_
    )))
  }
  # the SNP base in Primer3's 1-based window coordinates
  snp_i <- FLANK + 1L
  # 3' terminal base of each primer: 1-based template position of the last
  # base added (left: start + len - 1; right: start - len + 1, since a right
  # primer is given in genomic coordinates of its 5' end)
  left_3p <- pr$left_pos + nchar(pr$left) - 1L
  right_3p <- pr$right_pos - nchar(pr$right) + 1L
  # per pair gates; the locus passes a gate if ANY returned pair clears it
  pair_gate <- function(k) {
    primes <- c(pr$left[[k]], pr$right[[k]])
    h <- unname(hairpin_dg_thermo(primes, conditions = COND))
    dd <- unname(duplex_dg_thermo(pr$left[[k]], pr$right[[k]], conditions = COND))
    c(
      gc = all(dplyr::between(gc_frac_of(primes), gate_gc[1L], gate_gc[2L])),
      tm = all(dplyr::between(c(pr$left_tm[[k]], pr$right_tm[[k]]), gate_tm[1L], gate_tm[2L])),
      amp = dplyr::between(pr$product[[k]], gate_amp[1L], gate_amp[2L]),
      run = all(vapply(primes, longest_run_of, integer(1L)) < gate_run),
      p3p = !snp_i %in% c(left_3p[[k]], right_3p[[k]]),
      hairpin = all(h >= KASP_LIMITS$hairpin_dg),
      dimer = all(dd >= KASP_LIMITS$duplex_dg)
    )
  }
  pg <- t(vapply(seq_len(nrow(pr)), pair_gate, logical(7L)))
  colnames(pg) <- c("gc", "tm", "amp", "run", "p3p", "hairpin", "dimer")
  any_ok <- apply(pg, 2L, any)
  # the hardest gate for the accepted pairs: the gate the FEWEST pairs clear
  worst <- names(which.min(colSums(pg)))
  dplyr::bind_cols(base, tibble::tibble(
    p3_n_pairs = nrow(pr),
    p3_left = pr$left[[1L]], p3_right = pr$right[[1L]],
    p3_product = pr$product[[1L]],
    p3_left_tm = pr$left_tm[[1L]], p3_right_tm = pr$right_tm[[1L]],
    p3_left_3p_dist_snp = snp_i - left_3p[[1L]],
    p3_right_3p_dist_snp = snp_i - right_3p[[1L]],
    p3_gc_ok = any_ok[["gc"]], p3_tm_ok = any_ok[["tm"]], p3_amp_ok = any_ok[["amp"]],
    p3_run_ok = any_ok[["run"]], p3_3p_ok = any_ok[["p3p"]],
    p3_hairpin_ok = any_ok[["hairpin"]], p3_dimer_ok = any_ok[["dimer"]],
    p3_pairs_clearing_all = sum(rowSums(pg) == ncol(pg)),
    p3_worst_gate = worst
  ))
})

any_pair_gates <- c("p3_gc_ok", "p3_tm_ok", "p3_amp_ok", "p3_run_ok", "p3_3p_ok",
  "p3_hairpin_ok", "p3_dimer_ok")
cmp$p3_any_pair_pass <- ifelse(
  cmp$p3_n_pairs == 0L, NA,
  rowSums(as.matrix(cmp[any_pair_gates])) == length(any_pair_gates)
)
cmp$has_assay <- !is.na(cmp$tm_amplicon)

n_p3_pairs <- sum(cmp$p3_n_pairs)
n_p3_loci <- sum(cmp$p3_n_pairs > 0L)
n_pass_loci <- sum(cmp$p3_any_pair_pass %in% TRUE)
tract <- cmp[cmp$p3_any_pair_pass %in% TRUE, ]

# --- amplicon: the same locus, two designers --------------------------------
# Only loci where both sides produced something are comparable.  The sizes are
# not expected to agree: Primer3 optimises a 60-150 bp product, tidymarkers is
# gated to 50-100 bp, and the two pick different primer pairs on purpose.
both <- cmp[cmp$has_assay & is.finite(cmp$p3_product), ]
amp_d <- both$tm_amplicon - both$p3_product
amp_rr <- reg(both$tm_amplicon, both$p3_product)
cat(sprintf("\n--- amplicon: tidymarkers vs Primer3 product (n = %d loci with both) ---\n", nrow(both)))
if (nrow(both)) {
  cat(sprintf("  tidymarkers %d-%d bp (median %s); Primer3 %d-%d bp (median %s)\n",
    min(both$tm_amplicon), max(both$tm_amplicon), fmt(stats::median(both$tm_amplicon), 1L),
    min(both$p3_product), max(both$p3_product), fmt(stats::median(both$p3_product), 1L)))
  cat(sprintf(
    "  difference (tidymarkers - Primer3): median %s bp, range [%d, %d]; within 25 bp: %d/%d; rho = %s\n",
    fmt(stats::median(amp_d), 1L), min(amp_d), max(amp_d),
    sum(abs(amp_d) <= 25L), length(amp_d), fmt(amp_rr[["rho"]], 3L)))
}

# --- the shared question: does a tidymarkers gate reject an accepted pair? --
# `tract` holds the loci where at least one Primer3 pair clears all seven
# "any primer pair" gates.  The interesting cell is `!has_assay`: Primer3 says
# a usable pair exists, tidymarkers returned no KASP assay.
rejected <- tract[!tract$has_assay, ]
tm_only <- cmp[cmp$has_assay & !(cmp$p3_any_pair_pass %in% TRUE), ]
n_no_pair <- sum(cmp$p3_n_pairs == 0L)

cat(sprintf("\n--- Primer3 accepted, tidymarkers equivalent (n = %d loci on the panel) ---\n", nrow(cmp)))
cat(sprintf("  Primer3 returned >= 1 pair        : %d/%d (%d loci with none)\n",
  n_p3_loci, nrow(cmp), n_no_pair))
cat(sprintf("  ...that pair cleared all 7 gates  : %d/%d\n", nrow(tract), nrow(cmp)))
cat(sprintf("  ...and tidymarkers also designed  : %s  <- agreement floor\n",
  if (nrow(tract)) sprintf("%d/%d", sum(tract$has_assay), nrow(tract)) else "n/a (no such locus)"))
n_wide_long <- sum(cmp$p3_product_wide > gate_amp[2L], na.rm = TRUE)
cat(sprintf("  context: Primer3 on 60-150 bp would have used a product > %d bp at %d/%d loci\n",
  as.integer(gate_amp[2L]), n_wide_long, nrow(cmp)))
cat(sprintf("  Primer3 pair rejected by a gate   : %d/%d\n",
  sum(cmp$p3_any_pair_pass %in% FALSE), nrow(cmp)))
cat(sprintf("  tidymarkers designed, no complete Primer3 pair: %d\n", nrow(tm_only)))

if (nrow(rejected)) {
  # the note is a concatenation of criteria and warnings, with cli's ANSI codes
  clean_note <- function(x) {
    if (is.na(x)) return("(no note)")
    x <- gsub("\033\\[[0-9;]*m", "", x)
    x <- gsub("\\s+", " ", x)
    if (nchar(x) > 160L) paste0(substr(x, 1L, 157L), "...") else x
  }
  cat("\n  where they disagree (locus -> tidymarkers reason):\n")
  for (i in seq_len(nrow(rejected))) {
    cat(sprintf("    %s:%d %s>%s  v=%s  n_fail=%s  %s\n",
      rejected$chrom[[i]], rejected$pos[[i]], rejected$ref[[i]], rejected$alt[[i]],
      rejected$tm_verdict[[i]],
      ifelse(is.na(rejected$tm_n_fail[[i]]), "-", rejected$tm_n_fail[[i]]),
      clean_note(rejected$tm_note[[i]])))
  }
}

if (nrow(tm_only)) {
  cat("\n  tidymarkers designed where Primer3 had no complete pair; the gate that\n")
  cat("  barred the MOST Primer3 pairs at those loci: ")
  wg <- sort(table(tm_only$p3_worst_gate), decreasing = TRUE)
  cat(paste(names(wg), as.integer(wg), sep = "=", collapse = "; "), "\n")
  cat(sprintf("  (of %d loci, %d had no Primer3 pair at all)\n",
    nrow(tm_only), sum(is.na(tm_only$p3_worst_gate))))
}

# --- verdicts and where the assays sat inside the gates ---------------------
vd <- table(cmp$tm_verdict)
cat(sprintf("\n--- tidymarkers verdicts (n = %d) ---\n", nrow(cmp)))
cat(sprintf("  %s\n", paste(names(vd), as.integer(vd), sep = "=", collapse = ", ")))
with_a <- cmp[cmp$has_assay, ]
if (nrow(with_a)) {
  cat(sprintf(
    "  amplicon %d-%d bp, discrimination %s-%s, n_fail %s; strand %s\n",
    min(with_a$tm_amplicon), max(with_a$tm_amplicon),
    fmt(min(with_a$tm_discrimination), 2L), fmt(max(with_a$tm_discrimination), 2L),
    paste(sort(intersect(unique(with_a$tm_n_fail), with_a$tm_n_fail)), collapse = "/"),
    paste(ifelse(names(table(with_a$tm_strand)) == "-", "minus", "plus"),
      as.integer(table(with_a$tm_strand)), sep = "=", collapse = ", ")))
  cat("  criteria named on the non-passing assays (from marker_design_checks()): ")
  failed <- purrr::map_dfr(tm_checks, function(chk) {
    if (is.null(chk)) return(tibble::tibble(criterion = character()))
    tibble::tibble(criterion = chk$criterion[!is.na(chk$pass) & !chk$pass])
  })
  if (nrow(failed)) {
    ct <- sort(table(failed$criterion), decreasing = TRUE)
    cat(paste(names(ct), as.integer(ct), sep = "=", collapse = "; "), "\n")
  } else {
    cat("(none)\n")
  }
}

# --- write the numbers down ------------------------------------------------
thermo_out <- dplyr::bind_rows(
  tibble::tibble(
    item = tab_panel$oligo, seq_a = tab_panel$seq, seq_b = NA_character_,
    kind = tab_panel$kind, measure = "tm_c",
    pkg = tab_panel$pkg_tm, primer3 = tab_panel$p3_tm, note = NA_character_
  ),
  tibble::tibble(
    item = tab_panel$oligo, seq_a = tab_panel$seq, seq_b = NA_character_,
    kind = tab_panel$kind, measure = "tm_c_pkg_at_mg0",
    pkg = tab_panel$pkg_tm_nomg, primer3 = tab_panel$p3_tm, note = NA_character_
  ),
  tibble::tibble(
    item = tab_panel$oligo, seq_a = tab_panel$seq, seq_b = NA_character_,
    kind = tab_panel$kind, measure = "hairpin_dg",
    pkg = tab_panel$pkg_hairpin, primer3 = tab_panel$p3_hairpin,
    note = ifelse(tab_panel$p3_hairpin_none, "ntthal: no structure (0.0)", NA_character_)
  ),
  tibble::tibble(
    item = synth_dimer$pair, seq_a = synth_dimer$a, seq_b = synth_dimer$b,
    kind = "synth_pair", measure = "heterodimer_dg",
    pkg = synth_dimer$pkg, primer3 = synth_dimer$p3, note = NA_character_
  )
) |>
  dplyr::mutate(delta = .data$primer3 - .data$pkg) |>
  dplyr::arrange(.data$measure, .data$item)

utils::write.csv(thermo_out, file.path(res_dir, "primer3-thermo.csv"), row.names = FALSE)

cmp_out <- cmp |>
  dplyr::mutate(
    tm_minus_p3_product = .data$tm_amplicon - .data$p3_product,
    p3_any_pair_pass = .data$p3_any_pair_pass
  ) |>
  dplyr::arrange(.data$chrom, .data$pos)
utils::write.csv(cmp_out, file.path(res_dir, "primer3-design.csv"), row.names = FALSE)

# --- summary ----------------------------------------------------------------
cat("\n===== SUMMARY (the numbers the test may assert) =====\n")
smry <- c(
  n_oligo = nrow(tab_panel),
  tm_n = stat_panel$n_tm,
  tm_offset_mean = stat_panel$tm_mean, tm_offset_sd = stat_panel$tm_sd,
  tm_offset_absmax = stat_panel$tm_absmax,
  tm_r = stat_panel$r, tm_rho = stat_panel$rho,
  tm_offset_mean_pkg_mg0 = stat_panel$tm_mean_nomg,
  hairpin_n = as.numeric(stat_panel$hp[["n"]]),
  hairpin_offset_mean = stat_panel$hp[["delta_mean"]],
  hairpin_offset_sd = stat_panel$hp[["delta_sd"]],
  hairpin_rho = stat_panel$hp[["rho"]],
  dimer_n = as.numeric(stat_panel$dimer[["n"]]),
  dimer_offset_mean = stat_panel$dimer[["delta_mean"]],
  dimer_offset_sd = stat_panel$dimer[["delta_sd"]],
  n_loci = nrow(cmp),
  p3_loci_with_pair = n_p3_loci,
  p3_wide_long_product = n_wide_long,
  p3_loci_pair_clears_gates = nrow(tract),
  agreed = sum(tract$has_assay),
  tm_assays = sum(cmp$has_assay),
  amp_n = nrow(both),
  amp_median_diff = if (nrow(both)) stats::median(amp_d) else NA_real_,
  amp_rho = amp_rr[["rho"]],
  amp_within_25 = if (nrow(both)) sum(abs(amp_d) <= 25L) else 0L
)
print(round(smry, 3L))
cat(sprintf("\nwrote %s and %s\n",
  file.path("results", "primer3-thermo.csv"), file.path("results", "primer3-design.csv")))
