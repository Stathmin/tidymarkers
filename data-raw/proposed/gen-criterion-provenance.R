# Generates inst/extdata/criterion_provenance.csv (ROADMAP 223 + 1083, ZW4.6).
# Ids and sources are read from criterion_spec(); nothing here is hand-typed.
# Classes mirror the four reason classes in zw4.6-sources.md. Run from the
# package root with the R 4.6.1 env from AGENTS.md.
pkgload::load_all(".")

spec <- tidymarkers:::criterion_spec()
stopifnot(identical(anyDuplicated(spec$id), 0L))

# the chemistry's own key already sits on the assay row; these rows only restate
# the geometry that chemistry defines
inherited_geometry <- c(
  "family_crispr_pam", "family_lamp_fip_alt", "family_rpa_thf_off",
  "family_probe_extension", "family_probe_spans", "family_rpa_probe_spans",
  "family_rpa_variant_under", "family_crispr_seed", "family_crispr_guide_spans",
  "family_sanger_primers_inside", "family_sanger_pair_defined",
  "family_lamp_inner_inside", "family_lamp_fip_variant",
  "rpa:tier:tier3", "snapshot:gate:gate_snapshot_placement",
  "lamp:gate:gate_lamp_variant", "lamp:tier:tier1"
)

# constants of a scheme that no source in the collection prints
scheme_constant <- c(
  "family_caps_band_gap", "caps:gate:gate_caps_resolvable", "caps:tier:tier1",
  "arms:gate:gate_arms_bands", "arms:tier:tier2", "family_sanger_read_distance",
  "sanger:gate:gate_sanger_dist", "ssr", "family_lamp_inner_tm_balance",
  "family_size_repeat_unit", "family_size_repeat_copies", "family_size_delta",
  "max_run", "inverted_repeat", "family_size_fragment",
  "rpa:tier:tier2", "indel_size:gate:gate_indel_size_gap"
)

# thresholds applied to a number the package itself simulates
simulator_threshold <- c(
  "discrimination", "as_dropout", "as_foot_dropout", "com_dropout",
  "offtarget_hits", "in_silico_pcr", "family_crispr_offtarget",
  "ssr:gate:gate_ssr_offtarget", "family_probe_asp_discrimination",
  "taqman:gate:gate_taqman_discrimination", "taqman:tier:tier1",
  "duplex:gate:gate_duplex_targets", "duplex:gate:gate_duplex_complement",
  "duplex:gate:gate_duplex_selfcomp", "duplex:tier:tier1",
  "lamp:gate:gate_lamp_cross_dimer", "pair:tier:tier1",
  "lax:gate:gate_lax_cross_priming"
)

class_of <- function(id, label) {
  if (id %in% inherited_geometry) {
    return("inherited_geometry")
  }
  if (id %in% scheme_constant) {
    return("scheme_constant")
  }
  if (id %in% simulator_threshold) {
    return("simulator_threshold")
  }
  "conditions_unstated"
}

# one sentence per class (a few exceptions for the geometry rows, which name the
# concrete chemistry rather than the buffer)
geometry_reason <- function(id) {
  switch(id,
    family_crispr_pam = "PAM TTTV is the definition of the Cas12a chemistry; the assay row already cites its key.",
    "the bar is the definition of the scheme's own chemistry, whose key already sits on the assay row."
  )
}

reason_for <- function(id, label, cls) {
  switch(cls,
    conditions_unstated = sprintf(
      "no cited source states the %s bar under these reaction conditions.", label
    ),
    simulator_threshold = sprintf(
      "the %s tolerance thresholds a number the package simulates, so no wet-lab paper can state it.", label
    ),
    inherited_geometry = geometry_reason(id),
    scheme_constant = sprintf(
      "the %s constant is our own scheme default; no text in the collection prints this number.", label
    )
  )
}

named_reason <- function(src) {
  ifelse(grepl("^doi:", src), "inline snapshot rule states the bar.", "source states the bar.")
}

out <- data.frame(
  id = spec$id,
  source = spec$source,
  reason_class = character(nrow(spec)),
  reason = character(nrow(spec)),
  stringsAsFactors = FALSE
)
named <- out$source != "package assumption"
out$reason_class[named] <- "named_source"
out$reason[named] <- named_reason(out$source[named])
cl <- vapply(
  out$id[!named],
  function(id) class_of(id, spec$label[spec$id == id]),
  character(1)
)
lab <- spec$label[!named]
out$reason_class[!named] <- unname(cl)
out$reason[!named] <- mapply(reason_for, out$id[!named], lab, unname(cl), USE.NAMES = FALSE)

out <- out[order(out$id, method = "radix"), ]
stopifnot(all(nzchar(out$reason)), !anyNA(out$reason_class))
utils::write.csv(out, "inst/extdata/criterion_provenance.csv", row.names = FALSE, quote = TRUE)
print(table(out$reason_class))
