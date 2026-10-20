# =============================================================================
#  dev/render-gallery.R  --  one marker figure per scheme, for the vignettes
# =============================================================================
#  Every figure is produced by the package's own designer and `plot(marker)`;
#  nothing is drawn by hand.  The genome and a variant table are passed so the
#  same bricks as the KASP view are present: sequence letters, the variant
#  track and the genome-wide off-site panel.  PNG is for quick viewing; SVG is
#  the copyable version (sequence letters stay text, not paths).
#
#  Run:  Rscript dev/render-gallery.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

if (!nzchar(Sys.which("hybrid-min"))) {
  stop("OligoArrayAux not on PATH - see data-raw/install-oligoarrayaux.R")
}

out <- fs::path(root, "vignettes", "figures")
fs::dir_create(out)

bseq <- function(n, seed) {
  withr::with_seed(seed, paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = ""))
}

snp_fixture <- function(n = 600L, seed = 2L, pos = 300L) {
  W <- bseq(n, seed)
  ref <- substr(W, pos, pos)
  alt <- setdiff(c("A", "C", "G", "T"), ref)[1L]
  list(
    W = W, pos = pos,
    genome = Biostrings::DNAStringSet(c(Chr1 = W)),
    snp = tibble::tibble(chrom = "Chr1", pos = as.integer(pos), ref = ref, alt = alt)
  )
}

variants_of <- function(f) {
  tibble::tibble(
    chrom = "Chr1", pos = f$snp$pos, ref = f$snp$ref, alt = f$snp$alt,
    consequence = "missense_variant"
  )
}

save <- function(name, p, width = 12, height = 9) {
  f_png <- fs::path(out, paste0(name, ".png"))
  f_svg <- fs::path(out, paste0(name, ".svg"))
  ggplot2::ggsave(f_png, p, width = width, height = height, dpi = 130, bg = "white")
  # SVG keeps the sequence text as text, so a reader can copy letters and small
  # labels stay crisp at any zoom.
  svglite::svglite(f_svg, width = width, height = height, bg = "white")
  tryCatch(print(p), finally = grDevices::dev.off())
  cat("wrote", f_png, "and", f_svg, "\n")
}

save_marker <- function(name, m, variants = NULL, genome = NULL, width = 12, height = 9) {
  if (is.null(m) || inherits(m, "marker_doomed")) {
    cat("SKIP", name, ": doomed/empty\n")
    return(invisible(FALSE))
  }
  ok <- tryCatch(
    {
      # the feature track shares x with the sequence/map; give every scheme one
      features <- NULL
      pos <- if (!is.null(variants) && nrow(variants)) {
        range(variants$pos)
      } else if (!is.null(m$snp) && nrow(m$snp)) {
        range(m$snp$pos)
      } else {
        NULL
      }
      if (!is.null(pos)) {
        features <- tibble::tibble(
          start = pos[1L] - 120L, end = pos[2L] + 120L,
          strand = "+", name = "target"
        )
      }
      p <- if (is.null(variants) && is.null(genome)) {
        plot(m)
      } else {
        plot(m, variants = variants, genome = genome, features = features)
      }
      save(name, p, width, height)
      # QC the bricks, not just the file: primers, letters, off-site panel, checks
      sch <- m$scheme %||% sub("^marker_", "", class(m)[1L])
      pr <- marker_primers(m)
      pr <- pr[!is.na(pr) & nzchar(pr)]
      io <- tryCatch(marker_ideogram(m, genome), error = function(e) NULL)
      chk <- tryCatch(
        if (inherits(m, "marker_assay")) plot_checks(m) else plot_marker_checks(m, sch, genome),
        error = function(e) NULL
      )
      cat(
        "  QC ", name, ": primers=", length(pr), " nt=", sum(nchar(pr)),
        " offsite=", !is.null(io),
        " checks=", if (is.null(chk)) 0L else nrow(chk$data),
        " ok=", if (is.null(chk)) 0L else sum(chk$data$status == "ok", na.rm = TRUE),
        " broken=", if (is.null(chk)) 0L else sum(chk$data$status == "broken", na.rm = TRUE),
        "\n", sep = ""
      )
      TRUE
    },
    error = function(e) {
      cat("PLOT FAIL", name, ":", conditionMessage(e), "\n")
      FALSE
    }
  )
  invisible(ok)
}

build <- function(name, expr, variants = NULL, genome = NULL, width = 12, height = 9) {
  m <- tryCatch(suppressWarnings(eval(expr)), error = function(e) e)
  if (inherits(m, "error")) {
    cat("SKIP", name, ":", conditionMessage(m), "\n")
    return(invisible(FALSE))
  }
  save_marker(name, m, variants, genome, width, height)
}

# --- KASP family -------------------------------------------------------------
f <- snp_fixture()
a <- suppressWarnings(assemble_kasp_assay(
  f$genome, f$snp, variants = variants_of(f), flank = 150L, cores = 1L
))
save_marker("snp", a, variants_of(f), f$genome)

# insertion KASP: use the known-designable fixture from test-indel.R
Wi <- bseq(400L, 4L)
cand <- seq(160L, 240L, by = 1L)
clean <- cand[vapply(cand, function(p) {
  !has_run(as.character(Biostrings::subseq(Biostrings::DNAString(Wi), p - 35L, p + 35L)))
}, logical(1))]
gi <- Biostrings::DNAStringSet(c(Chr1 = Wi))
for (p in clean) {
  ref_i <- substr(Wi, p, p)
  ins <- tibble::tibble(chrom = "Chr1", pos = p, ref = ref_i, alt = paste0(ref_i, "GGG"))
  ai <- tryCatch(
    suppressWarnings(assemble_kasp_assay(gi, ins, variants = NULL, flank = 150L, cores = 1L)),
    error = function(e) NULL
  )
  if (inherits(ai, "marker_assay")) break
}
if (exists("ai") && inherits(ai, "marker_assay")) {
  save_marker(
    "indel", ai,
    tibble::tibble(chrom = "Chr1", pos = ai$snp$pos, ref = ai$snp$ref, alt = ai$snp$alt,
      consequence = "intron_variant"),
    gi
  )
} else {
  cat("SKIP indel: no designable insertion\n")
}

# multiallelic KASP (degenerate ALT primer)
fm <- snp_fixture(seed = 2L, pos = 200L)
others <- setdiff(c("A", "C", "G", "T"), fm$snp$ref)
fm$snp$alt2 <- others[2L]
fm$snp$alt <- others[1L]
am <- suppressWarnings(assemble_kasp_assay(fm$genome, fm$snp, flank = 150L, cores = 1L))
save_marker("multi", am, variants_of(fm), fm$genome)

# duplex KASP (two loci)
fd <- snp_fixture(n = 800L, seed = 1L, pos = 250L)
bases <- c(substr(fd$W, 250L, 250L), substr(fd$W, 550L, 550L))
duplex_loci <- tibble::tibble(
  chrom = "Chr1", pos = c(250L, 550L), ref = bases,
  alt = vapply(bases, function(b) setdiff(c("A", "C", "G", "T"), b)[1L], character(1))
)
duplex_vars <- dplyr::mutate(duplex_loci, consequence = "missense_variant")
build("duplex", quote(design_duplex(fd$genome, duplex_loci, flank = 150L, cores = 1L)),
  duplex_vars, fd$genome, 12, 10
)

# haplotype-aware lax KASP on the bundled SAT2 transcripts
f_sat <- system.file("extdata", "sat2_transcripts.fa", package = "tidymarkers")
if (nzchar(f_sat)) {
  hp <- haplotypes_from_alignment(read_transcripts(f_sat), chrom = "SAT2")
  ch <- marker_chunk(hp, from = 281L, to = 581L)
  lax_snp <- tibble::tibble(chrom = "SAT2", pos = 151L, ref = "A", alt = "G")
  build(
    "lax",
    quote(design_lax_marker(
      chunk_genome(ch), lax_snp, flank = 150L, cores = 1L,
      haplotypes = chunk_haplotypes(ch)
    )),
    chunk_variants(ch), chunk_genome(ch), 12, 10
  )
}

# ASQ has no separate engine: the KASP primer anatomy is what the scheme describes
fa <- snp_fixture()
a_asq <- suppressWarnings(assemble_kasp_assay(
  fa$genome, fa$snp, variants = variants_of(fa), flank = 150L, cores = 1L,
  scheme = "asq"
))
save_marker("asq", a_asq, variants_of(fa), fa$genome)

# the marker-agnostic primer core
fp <- snp_fixture()
a_primer <- suppressWarnings(assemble_kasp_assay(
  fp$genome, fp$snp, variants = variants_of(fp), flank = 150L, cores = 1L,
  scheme = "primer"
))
save_marker("primer", a_primer, variants_of(fp), fp$genome)

# --- gel / melt family -------------------------------------------------------
fc <- snp_fixture()
Wc <- fc$W
substr(Wc, 300L, 305L) <- "GAATTC"
fc$genome <- Biostrings::DNAStringSet(c(Chr1 = Wc))
fc$snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "G", alt = "A")
build("caps", quote(design_caps_marker(fc$genome, fc$snp, flank = 250L)),
  variants_of(fc), fc$genome
)

Warms <- bseq(800L, 1L)
farms <- list(
  W = Warms,
  genome = Biostrings::DNAStringSet(c(Chr1 = Warms)),
  snp = tibble::tibble(
    chrom = "Chr1", pos = 400L, ref = substr(Warms, 400L, 400L),
    alt = setdiff(c("A", "C", "G", "T"), substr(Warms, 400L, 400L))[1L]
  )
)
build("arms", quote(design_arms_marker(farms$genome, farms$snp, flank = 250L)),
  variants_of(farms), farms$genome
)

fh <- snp_fixture()
build("hrm", quote(design_hrm_marker(fh$genome, fh$snp, flank = 150L)),
  variants_of(fh), fh$genome
)

Ws <- bseq(600L, 7L)
at <- 296L
ssr <- "ACACACACAC"
Ws <- paste0(substr(Ws, 1L, at - 1L), ssr, substr(Ws, at + nchar(ssr), 600L))
fs <- list(
  genome = Biostrings::DNAStringSet(c(Chr1 = Ws)),
  snp = tibble::tibble(chrom = "Chr1", pos = at + nchar(ssr) %/% 2L, ref = "A", alt = "G")
)
build("ssr", quote(design_ssr_marker(fs$genome, fs$snp, flank = 250L)),
  NULL, fs$genome
)

fis <- snp_fixture(seed = 2L)
ref_is <- substr(fis$W, 300L, 300L)
fis$snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = ref_is, alt = paste0(ref_is, "TTTTTTTTTT"))
build("indel_size", quote(design_indel_size_marker(fis$genome, fis$snp, flank = 250L, min_gap = 5L)),
  variants_of(fis), fis$genome
)

# --- probe family ------------------------------------------------------------
ft <- snp_fixture()
build("taqman", quote(design_taqman_marker(ft$genome, ft$snp, flank = 150L)),
  variants_of(ft), ft$genome
)
fd2 <- snp_fixture(seed = 1L)
build("dpcr", quote(design_dpcr_marker(fd2$genome, fd2$snp, flank = 150L)),
  variants_of(fd2), fd2$genome
)
fsn <- snp_fixture(n = 800L, seed = 1L, pos = 400L)
build("snapshot", quote(design_snapshot_marker(fsn$genome, fsn$snp, flank = 250L)),
  variants_of(fsn), fsn$genome, 12, 10
)
fca <- snp_fixture()
build("castpcr", quote(design_castpcr_marker(fca$genome, fca$snp, flank = 150L)),
  variants_of(fca), fca$genome
)

# --- isothermal --------------------------------------------------------------
fl <- snp_fixture()
build("lamp", quote(design_lamp_marker(fl$genome, fl$snp, flank = 150L)),
  variants_of(fl), fl$genome
)
fr <- snp_fixture()
build("rpa", quote(design_rpa_marker(fr$genome, fr$snp, flank = 200L, design_probe = TRUE)),
  variants_of(fr), fr$genome
)

# --- CRISPR / Sanger ---------------------------------------------------------
Wcr <- bseq(600L, 7L)
substr(Wcr, 293L, 296L) <- "TTTA"
substr(Wcr, 297L, 316L) <- "ACGTGCTAGCTAGCATGCA"
fcr <- list(
  genome = Biostrings::DNAStringSet(c(Chr1 = Wcr)),
  snp = tibble::tibble(chrom = "Chr1", pos = 300L, ref = "T", alt = "A")
)
build("crispr", quote(design_crispr_marker(fcr$genome, fcr$snp, flank = 200L, cas = "cas12a")),
  variants_of(fcr), fcr$genome, 12, 10
)

Wsa <- bseq(1200L, 2L)
fsa <- list(
  genome = Biostrings::DNAStringSet(c(Chr1 = Wsa)),
  snp = tibble::tibble(
    chrom = "Chr1", pos = 600L, ref = substr(Wsa, 600L, 600L),
    alt = setdiff(c("A", "C", "G", "T"), substr(Wsa, 600L, 600L))[1L]
  )
)
build("sanger", quote(design_sanger_marker(fsa$genome, fsa$snp, flank = 500L)),
  variants_of(fsa), fsa$genome, 12, 10
)

cat("\ndone: ", length(fs::dir_ls(out, glob = "*.svg")), " figures\n", sep = "")
