# =============================================================================
#  dev/targets/_targets.R  --  the freshness holder for generated artifacts
# =============================================================================
#  `targets` decides *when* to redraw, not *what*: every target calls the same
#  script dev/render-all.R calls, so the engine stays the only author.  The
#  payoff is content hashing -- a figure or the manifest is regenerated only
#  when something it actually reads changed, and `tar_outdated()` answers
#  "would a fresh clone still see this file as current?".
#
#  Do not run this file directly; dev/targets.R owns the store path and mode.
#  The heavy gallery panels are deliberately not here: they need OligoArrayAux
#  and would make every gate run depend on it (see dev/render-all.R --render).
#
#  Z8.2 harness: `panel_outcomes` joins the per-assay wet-lab label
#  (soybean_worked_panel_assays.csv: `outcome` + `basis`) to the per-assay gate
#  verdicts (soybean_worked_panel_kasp_gates.csv) and writes
#  inst/extdata/z8_gate_outcomes.csv.  It is a table, not a model: no rate and
#  no probability is computed while the failure class is this thin.
# =============================================================================

library(targets)
tar_option_set(packages = character())

# targets inherits the caller's working directory and dev/targets.R runs from
# the package root, so relative paths here mean what they mean everywhere else
if (!file.exists("DESCRIPTION")) stop("run from the package root (use dev/targets.R)")

# Z8.4: the source-level split and the gate-rate check live in one file, so the
# tests source the harness rather than re-declaring it.  No file artifact: the
# split is a check, not a generated table.
source("dev/z8-harness.R")

# every file the renderer reads: the engine, its data, and the script itself
# z8_gate_outcomes.csv is the harness's own output, so reading it would make
# the harness depend on its last result
engine_inputs <- function() {
  sort(c(
    list.files("R", pattern = "[.]R$", full.names = TRUE),
    setdiff(
      list.files("inst/extdata", pattern = "[.]csv$", full.names = TRUE),
      "inst/extdata/z8_gate_outcomes.csv"
    ),
    "dev/render-decision.R"
  ))
}

render_decision <- function(engine) {
  status <- system2(file.path(R.home("bin"), "Rscript"), "dev/render-decision.R")
  if (!identical(as.integer(status), 0L)) stop("dev/render-decision.R failed (", status, ")")
  sort(c(
    Sys.glob("man/figures/*.svg"),
    Sys.glob("vignettes/figures/decision-*.svg"),
    Sys.glob("vignettes/figures/decision-*.txt")
  ))
}

# the manifest is a pure function of the tree, so refreshing it here is the same
# work `dev/artifacts.R --write` does -- no second hashing contract
refresh_manifest <- function(figures) {
  source("dev/artifacts.R", local = TRUE)
  artifacts <- refresh_hashes(read_manifest())
  write_manifest(artifacts)
  write_projection(artifacts)
  c(ARTIFACT_MANIFEST, ARTIFACT_PROJECTION)
}

## --- Z8.2: gate verdict -> wet-lab outcome -----------------------------------

# The vocabulary is closed so a typo cannot become a seventh class.
z8_outcomes <- c(
  "worked", "failed", "ambiguous", "no_amplification", "no_polymorphism",
  "not_attempted"
)

# Assay IDs are the panel's own, so the join needs no new key.  A failed gate
# is a DESIGN verdict; only the wet-lab `outcome` may say `failed`, and only
# with the observation quoted in `basis`.
panel_outcomes <- function(assays, gates) {
  stopifnot(
    !anyDuplicated(assays$assay_id),
    all(assays$outcome %in% z8_outcomes),
    all(is.na(assays$outcome) == is.na(assays$basis)),
    all(!grepl("^failed$", assays$outcome) | nzchar(assays$basis)),
    all(gates$status %in% c("pass", "fail", "unknown"))
  )
  g <- table(factor(gates$assay_id, assays$assay_id), gates$status)
  out <- tibble::tibble(
    source = "operator_panel_soybean",
    row = assays$assay_id,
    scheme = assays$scheme,
    outcome = assays$outcome,
    basis = assays$basis,
    n_gate_pass = as.integer(g[, "pass"]),
    n_gate_fail = as.integer(g[, "fail"]),
    n_gate_unknown = as.integer(g[, "unknown"]),
    gate_status = ifelse(
      !as.character(assays$assay_id) %in% gates$assay_id, "not_evaluated",
      ifelse(g[, "fail"] > 0L, "fail", "pass")
    )
  )
  stopifnot(!anyDuplicated(out[c("source", "row")]))
  out[order(out$row), ]
}

# Byte-identical across runs and machines: the body is written from the
# recomputed table, so a re-run is a no-op and a diff means the data changed.
write_z8_outcomes <- function(tbl, path = "inst/extdata/z8_gate_outcomes.csv") {
  con <- file(path, open = "wb")
  on.exit(close(con))
  readr::write_csv(tbl, con, na = "")
  path
}

list(
  tar_target(engine, engine_inputs(), format = "file"),
  tar_target(figures, render_decision(engine), format = "file"),
  tar_target(manifest, refresh_manifest(figures), format = "file"),
  tar_target(
    panel_assays_file,
    "inst/extdata/soybean_worked_panel_assays.csv",
    format = "file"
  ),
  tar_target(
    panel_gates_file,
    "inst/extdata/soybean_worked_panel_kasp_gates.csv",
    format = "file"
  ),
  # The paths are file targets so editing either CSV invalidates the join;
  # reading the literal path would let the harness go stale silently.
  tar_target(
    panel_assays,
    readr::read_csv(panel_assays_file,
      show_col_types = FALSE, na = c("", "NA")
    ),
    format = "rds"
  ),
  tar_target(
    panel_gates,
    readr::read_csv(panel_gates_file,
      show_col_types = FALSE
    ),
    format = "rds"
  ),
  tar_target(
    z8_outcomes_file,
    write_z8_outcomes(panel_outcomes(panel_assays, panel_gates)),
    format = "file"
  )
)
