# =============================================================================
#  i18n: one chokepoint for translated abort messages
# =============================================================================
#  Every user-facing abort goes through `abort_msg()`.  It translates the
#  message and each bullet exactly once, keeps cli's names (so `i`, `x`, `v`
#  bullets render as before) and forwards `class`, `call` and `parent` to
#  `cli::cli_abort()`.  `po/README.md` registers it as a custom translation
#  function so `potools::po_extract()` still sees the literals.

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
  env <- parent.frame()
  bullets <- list(...)
  msg <- gettext(message)
  if (length(bullets)) {
    nm <- names(bullets)
    if (is.null(nm)) nm <- rep("", length(bullets))
    vals <- vapply(bullets, function(x) as.character(gettext(x)), character(1))
    msg <- c(msg, stats::setNames(vals, nm))
  }
  cli::cli_abort(msg,
    class = class, call = call,
    parent = parent, .envir = env
  )
}
