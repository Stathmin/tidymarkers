# =============================================================================
#  pipeline.R  --  the layered workflow: plan -> solve -> tubes -> export
# =============================================================================
#  Run from the project root:   Rscript inst/scripts/pipeline.R
#
#  This is the "how a user programs with tidymarkers" script.  It builds a *plan* in
#  layers (the grammar), solves it, checks which results can share a tube, then
#  plots and exports.
#
#    target layer :  marker_plan(genome) + marker_targets(...) + marker_context(...)
#    tech layer   :  + marker_tech(preset, schemes, costs, exclude)
#    solver layer :  marker_solve(plan)
#    tube layer   :  marker_tubes(solved)
#    viz/export   :  marker_plot(solved), marker_export(solved, dir)
#
#  Needs OligoArrayAux (tidymarkers::marker_install_tools()); TAIR10 is cached on the
#  first ath_genome() call.
# =============================================================================

suppressPackageStartupMessages({
  library(tidymarkers)
  library(ggplot2)
})

out_dir <- "results/plan"
fs::dir_create(out_dir)

# -----------------------------------------------------------------------------
# 1. TARGET LAYER - what genome, which loci, what context
# -----------------------------------------------------------------------------
genome <- ath_genome() # TAIR10 (cached)
ex <- marker_example() # the curated real panel

# Per-locus context: the neighbourhood variants of each target (its gene +/- 500
# bp, or +/- 1 kb when the gene is unknown) - the same idiom as
# dev/render-targets.R.  A design should know what else is polymorphic where it
# works, so a primer under a real variant is scored for it.
locus_variants <- function(chrom, pos, gene_id) {
  g <- tryCatch(ensembl_gene(gene_id), error = function(e) NULL)
  win <- if (!is.null(g)) c(g$start - 500L, g$end + 500L) else c(pos - 1000L, pos + 1000L)
  ensembl_variants(chrom, max(1L, win[1L]), win[2L])
}
context <- purrr::map(seq_len(nrow(ex)), function(i) {
  locus_variants(ex$chrom[i], as.integer(ex$pos[i]), ex$gene[i])
}) |>
  purrr::list_rbind() |>
  dplyr::bind_rows(locus_variants("Chr4", 269116L, "AT4G00650")) # the extra FRI target

plan <- marker_plan(genome) +
  # the bundled targets (their `id` labels the report, `alt2` keeps multi-allelic)
  marker_targets(ex) +
  # one ad-hoc target, labelled for the report
  marker_targets(
    chrom = "Chr4", pos = 269116L, ref = "G", alt = "A",
    type = "snp", label = "FRI G>A"
  ) +
  # the design should know what else is polymorphic there
  marker_context(variants = context)

print(plan) # cheap: nothing has been designed yet

# -----------------------------------------------------------------------------
# 2. TECH LAYER - which chemistries are on hand, and how costly each is
# -----------------------------------------------------------------------------
# Presets: "kasp" (KASP only), "gel" (KASP + fragment sizes), "lab" (everything
# bench-standard), "field" (isothermal / lateral-flow).  Or list schemes yourself;
# a custom `design_<name>_marker()` + `scheme_def_<name>()` is discovered too.
plan <- plan + marker_tech(
  preset = "lab",
  # cost schema: the solver walks cheapest-first and STOPS - a KASP win means the
  # fragment analysis is never run for that target.  Say "KASP is 10x easier":
  costs = marker_costs(snp = 1, caps = 10),
  exclude = "snapshot" # not set up in our lab
)

# -----------------------------------------------------------------------------
# 3. SOLVER LAYER - run it (cache-aware, progress, tidy report)
# -----------------------------------------------------------------------------
cache <- new.env() # reused across solves
solved <- marker_solve(plan, flank = 150L, cores = 2L, cache = cache)
print(solved) # target -> chemistry -> cost

# the per-scheme attempts: what was tried, what won, what was never needed
print(solved$attempts)

# solving the same plan again is free (identical designs come from the cache)
invisible(marker_solve(plan, flank = 150L, cores = 2L, cache = cache))

# -----------------------------------------------------------------------------
# 4. TUBE LAYER - which results can be multiplexed
# -----------------------------------------------------------------------------
tubes <- marker_tubes(solved)
print(tubes) # clashing primers flagged

# -----------------------------------------------------------------------------
# 5. VISUALIZATION + EXPORT
# -----------------------------------------------------------------------------
ggsave(fs::path(out_dir, "plan.png"), marker_plot(solved), width = 9, height = 5, dpi = 110)
marker_export(solved, dir = out_dir)

# the per-target marker objects are on the result, so nothing is recomputed:
first_ok <- which(solved$report$ok)[1L]
if (length(first_ok)) {
  cat("\nFirst solved marker (", solved$report$scheme[first_ok], "):\n", sep = "")
  print(solved$report$marker[[first_ok]])
}

message("Done. See ", out_dir, "/")

# tidy the BiocParallel forks (an otherwise cosmetic exit-time message)
if (requireNamespace("BiocParallel", quietly = TRUE)) try(BiocParallel::bpstop(), silent = TRUE)
