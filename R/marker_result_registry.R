# =============================================================================
#  marker_result_registry.R  --  one S3 method per result generic
# =============================================================================
#  `plot()`, `print()` and `marker_map_spec()` dispatch through a named list
#  `class -> function` instead of one S3 method per class.
#  `marker_registry_lookup()` walks the whole `class()` vector, so a wrapper
#  such as `marker_explain`, which prepends its own class, still reaches the
#  underlying result's entry.

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

#' Print methods by result class
#' @noRd
marker_print_registry <- function() {
  list(
    marker_assay = marker_print_assay,
    marker_doomed = marker_print_doomed,
    marker_caps = marker_print_caps,
    marker_arms = marker_print_arms,
    marker_hrm = marker_print_hrm,
    marker_taqman = marker_print_taqman,
    marker_castpcr = marker_print_castpcr,
    marker_dpcr = marker_print_dpcr,
    marker_snapshot = marker_print_snapshot,
    marker_lamp = marker_print_lamp,
    marker_rpa = marker_print_rpa,
    marker_crispr = marker_print_crispr,
    marker_sanger = marker_print_sanger,
    marker_ssr = marker_print_ssr,
    marker_indel_size = marker_print_indel_size,
    marker_duplex = marker_print_duplex,
    marker_lax = marker_print_lax
  )
}

#' Map-spec methods by result class
#' @noRd
marker_map_spec_registry <- function() {
  list(
    marker_assay = marker_map_assay,
    marker_caps = marker_map_caps,
    marker_arms = marker_map_arms,
    marker_hrm = marker_map_hrm,
    marker_taqman = marker_map_taqman,
    marker_castpcr = marker_map_castpcr,
    marker_dpcr = marker_map_dpcr,
    marker_snapshot = marker_map_snapshot,
    marker_lamp = marker_map_lamp,
    marker_rpa = marker_map_rpa,
    marker_crispr = marker_map_crispr,
    marker_sanger = marker_map_sanger,
    marker_ssr = marker_map_ssr,
    marker_indel_size = marker_map_indel_size,
    marker_duplex = marker_map_duplex
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

#' @export
print.marker <- function(x, ...) {
  # marker_explain prepends a headline and then prints the wrapped result
  if (inherits(x, "marker_explain")) {
    cat(attr(x, "headline"), "\n", sep = "")
  }
  marker_registry_lookup(marker_print_registry(), x, "print")(x, ...)
}

#' @export
marker_map_spec.marker <- function(marker, scheme = NULL, ...) {
  marker_registry_lookup(marker_map_spec_registry(), marker, "map_spec")(
    marker, scheme, ...
  )
}
