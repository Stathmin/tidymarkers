# =============================================================================
#  marker_result_registry.R  --  one S3 method per result generic
# =============================================================================
#  `plot()` (and, as more domains land, `print()` and `marker_map_spec()`)
#  dispatch through a named list `class -> function` instead of one S3 method
#  per class.  `marker_registry_lookup()` walks the whole `class()` vector, so a
#  wrapper such as `marker_explain`, which prepends its own class, still reaches
#  the underlying result's entry.

#' Plot methods by result class
#' @noRd
marker_plot_registry <- function() {
  list(
    marker_assay = marker_plot_assay,
    marker_doomed = marker_plot_doomed,
    marker_caps = marker_plot_caps,
    marker_arms = marker_plot_arms,
    marker_hrm = marker_plot_hrm,
    marker_taqman = marker_plot_taqman,
    marker_castpcr = marker_plot_castpcr,
    marker_dpcr = marker_plot_dpcr,
    marker_snapshot = marker_plot_snapshot,
    marker_lamp = marker_plot_lamp,
    marker_rpa = marker_plot_rpa,
    marker_crispr = marker_plot_crispr,
    marker_sanger = marker_plot_sanger,
    marker_ssr = marker_plot_ssr,
    marker_indel_size = marker_plot_indel_size,
    marker_duplex = marker_plot_duplex,
    marker_lax = marker_plot_lax
  )
}

#' Find the registered method for an object
#' @noRd
marker_registry_lookup <- function(registry, x, what) {
  for (cls in class(x)) {
    fn <- registry[[cls]]
    if (is.function(fn)) {
      return(fn)
    }
  }
  abort_msg(
    "No {what} formatter is registered for this object.",
    x = "Classes: {.cls {class(x)}}.",
    i = "Add an entry to marker_{what}_registry().",
    class = "tidymarkers_registry_error"
  )
}

#' @rdname plot-marker
#' @export
plot.marker <- function(x, variants = NULL, genome = NULL, features = NULL, ...) {
  # same formals as the retired methods, so a positional `plot(m, variants)`
  # still binds `variants` instead of being swallowed by the generic's `y`
  marker_registry_lookup(marker_plot_registry(), x, "plot")(
    x, variants, genome, features, ...
  )
}
