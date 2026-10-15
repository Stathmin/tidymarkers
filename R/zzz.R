# =============================================================================
#  zzz.R  --  package session state
# =============================================================================
#  A holder for one-shot user hints (e.g. "OligoArrayAux missing"), so the
#  package never writes a `options()` entry or a global variable to say "already
#  told them".  A package-namespace environment is the standard place for such
#  state (memoise caches use the same idea).
# =============================================================================

.kasp_state <- new.env(parent = emptyenv())
