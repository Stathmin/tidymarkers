# =============================================================================
#  dev/profile-gif.R  --  where does the design-search GIF spend its time?
# =============================================================================
#  Builds a small random assay and times the phases of render_design_search():
#  drawing the frames and stitching the GIF, with and without the ragg device
#  and across cores.  The numbers say whether the frames or the assembler
#  dominate, and what the parallel/ragg win actually is.
#
#  Needs OligoArrayAux on PATH (see data-raw/install-oligoarrayaux.R).
#  Run:  Rscript dev/profile-gif.R
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)
suppressMessages(pkgload::load_all(root, quiet = TRUE))

if (!nzchar(Sys.which("hybrid-min"))) {
  stop("OligoArrayAux not on PATH - see data-raw/install-oligoarrayaux.R")
}

# --- a small real assay (same recipe as test-search-anim.R) -----------------
make_assay <- function(seed = 4L, n = 500L, pos = 250L) {
  for (s in seed:(seed + 19L)) {
    set.seed(s)
    genome <- Biostrings::DNAStringSet(c(Chr1 = paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")))
    ref <- as.character(Biostrings::subseq(genome[["Chr1"]], pos, pos))
    snp <- tibble::tibble(chrom = "Chr1", pos = pos, ref = ref,
                          alt = setdiff(c("A", "C", "G", "T"), ref)[1L])
    vs <- tibble::tibble(chrom = "Chr1", pos = pos, consequence = "missense_variant")
    a <- try(suppressWarnings(assemble_kasp_assay(genome, snp, variants = vs, flank = 150L)), silent = TRUE)
    if (!inherits(a, "try-error")) return(list(assay = a, snp = snp, variants = vs))
  }
  stop("no random genome assembled into an assay")
}

m     <- make_assay()
frames <- design_search_frames(m$assay, m$variants, m$snp)
n      <- nrow(frames$notes)
cat(sprintf("frames: %d steps\n\n", n))

td <- fs::path_temp("kasp_anim"); fs::dir_create(td)
pngs <- fs::path(td, sprintf("f%04d.png", seq_len(n)))

draw <- function(device = NULL) {
  for (i in seq_len(n)) {
    ggplot2::ggsave(pngs[i], plot_design_search_step(frames, i),
                    width = 8, height = 4.6, dpi = 110, bg = "white", device = device)
  }
}


agg <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
tm  <- function(label, expr) {
  t <- system.time(expr)[["elapsed"]]
  cat(sprintf("  %-26s %6.2f s\n", label, t))
}
cat("frame drawing:\n")
tm("default device, serial", draw())
tm("ragg, serial", draw(agg))

cat("\nper-frame cost breakdown (first step):\n")
p <- plot_design_search_step(frames, 1L)
tm("build plot", p <- plot_design_search_step(frames, 1L))
tm("ggplot_build", ggplot2::ggplot_build(p))
tm("ggsave render", ggplot2::ggsave(fs::path(td, "one.png"), p,
  width = 8, height = 4.6, dpi = 110, bg = "white"
))

cat("\nassembly (ImageMagick):\n")
tm("assemble_gif", tidykasp:::assemble_gif(pngs, fs::path(td, "out.gif"), delay_cs = 25L))

cat("\nrender_design_search end-to-end:\n")
tm("serial (default)", render_design_search(m$assay, m$variants, m$snp,
  path = fs::path(td, "serial.gif"), fps = 4L, hold = 3L, cores = 1L
))
tm(sprintf("parallel (%d cores)", min(parallel::detectCores(), 40L)), render_design_search(m$assay, m$variants, m$snp,
  path = fs::path(td, "parallel.gif"), fps = 4L, hold = 3L, cores = parallel::detectCores()
))
