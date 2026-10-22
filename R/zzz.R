# =============================================================================
#  zzz.R  --  package session state
# =============================================================================
#  A holder for one-shot user hints (e.g. "OligoArrayAux missing"), so the
#  package never writes a `options()` entry or a global variable to say "already
#  told them".  A package-namespace environment is the standard place for such
#  state (memoise caches use the same idea).
# =============================================================================

.marker_state <- new.env(parent = emptyenv())

# Genome-derived, pattern-independent facts about a sequence: currently the runs
# of unambiguous bases used by `exact_pattern_hits()`.  Keyed by a content
# signature, so a changed sequence never reuses a stale entry.
.marker_state$blocks <- new.env(parent = emptyenv())
