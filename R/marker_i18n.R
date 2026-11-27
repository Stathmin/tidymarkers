# =============================================================================
#  i18n: one chokepoint for translated user messages
# =============================================================================
#  Every user-facing abort goes through `abort_msg()`, every warning through
#  `warn_msg()`.  Both translate the message and each bullet exactly once,
#  keep cli's names (so `i`, `x`, `v` bullets render as before) and forward
#  `class`, `call` and `parent` to the corresponding `cli::cli_*()` call.
#  `po/README.md` registers them as custom translation functions so
#  `potools::po_extract()` still sees the literals.

# Translate the message and its bullets once, keeping the bullet names.
translate_msg <- function(message, bullets) {
  msg <- gettext(message)
  if (length(bullets)) {
    nm <- names(bullets)
    if (is.null(nm)) nm <- rep("", length(bullets))
    vals <- vapply(bullets, function(x) as.character(gettext(x)), character(1))
    msg <- c(msg, stats::setNames(vals, nm))
  }
  msg
}

#' Translate and raise a user-facing abort
#'
#' `message` is one literal and `...` are named cli bullets, all translated
#' once via `gettext()`.  `{...}` is interpolated in the calling frame.
#' @param message a single message literal.
#' @param ... named bullet literals.
#' @param class condition class(es) passed to [cli::cli_abort()].
#' @param call the call shown to the user.
#' @param parent an optional parent condition to chain.
#' @keywords internal
#' @noRd
abort_msg <- function(message, ..., class = NULL,
                      call = rlang::caller_env(), parent = NULL) {
  cli::cli_abort(translate_msg(message, list(...)),
    class = class, call = call,
    parent = parent, .envir = parent.frame()
  )
}

#' Translate and raise a user-facing warning
#'
#' Same contract as `abort_msg()`, but raises a warning instead of an error.
#' @param message a single message literal.
#' @param ... named bullet literals.
#' @keywords internal
#' @noRd
warn_msg <- function(message, ...) {
  cli::cli_warn(translate_msg(message, list(...)), .envir = parent.frame())
}
