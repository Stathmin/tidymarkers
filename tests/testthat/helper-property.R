# =============================================================================
#  helper-property.R  --  generators for the property tests (ZY11)
# =============================================================================
#  Only cheap, pure inputs are generated here: a primer-length ACGT string and
#  the scheme ids the registry actually declares.  Nothing in this file builds a
#  genome, an assay or a decision graph.
# =============================================================================

# A random ACGT string of n_nt bases (n_nt in 1..20 keeps every profile cheap).
prop_dna <- function(n_nt) {
  paste(sample(c("A", "C", "G", "T"), n_nt, replace = TRUE), collapse = "")
}

# The registered scheme ids, as the package itself lists them.
prop_schemes <- function() marker_schemes()$scheme
