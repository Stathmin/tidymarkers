# =============================================================================
#  kasp_asq.R  --  ASQ (allele-specific qPCR): a universal-probe readout of KASP
# =============================================================================
#  ASQ reads the same competitive allele-specific extension as KASP, but its
#  distinguishing features are *not* sequence properties of a candidate:
#   * the variant sits at the penultimate (or antepenultimate) 3' base of the
#     allele-specific primer, not the terminal base KASP uses; and
#   * the reporter is a set of universal probes (UP1-UP4, FAM/HEX/Cy3/Cy5) read
#     against a *single* common quencher (Uni-Q), i.e. a master-mix / instrument
#     (panel) property, up to 4-plex.
#
#  Design rules (researched, not invented).  Kalendar R, Shustov AV, Akhmetollayev
#  I, Kairov U (2022) "Designing Allele-Specific Competitive-Extension PCR-Based
#  Assays for High-Throughput Genotyping and Gene Characterization", Front Mol
#  Biosci 9:773956, doi:10.3389/fmolb.2022.773956 (open access, PMC8921500):
#   * "each allele-specific primer is designed such that the penultimate base
#     before the primer's 3' end base is positioned at the SNP site" - the tool
#     "automatically places a variable base at the second (penultimate) position
#     of the primer's 3' terminus", or the third, chosen thermodynamically;
#   * the FRET cassette becomes "two or more universal probes (UP-1-4) with
#     corresponding tags and different fluorescent dyes at the 5' end and a single
#     common UP with a quencher at the 3' end (Uni-Q)" - "4 in ASQ, 2 in KASP";
#   * the allele-specific primers "were selected to obtain a Tm close to 63 C
#     (range 62-64 C)"; the universal probes are designed higher (66-68 C) and the
#     quencher lower (53-55 C), the signal being read at a second 55 C step;
#   * "amplicons ranged in size from 100 to 200 bases" and encompassed the SNP.
#
#  Documentation-only (C13).  This scheme is *described*, not *designed* here:
#   (1) the reporter/plexity rules are panel and format properties - they are not
#       computable from a candidate, so they cannot be selection gates (C7/C8).
#       The two gates this stage used to declare (`gate_asq_cassette`,
#       `gate_asq_readout`) were exactly that mistake and are gone.
#   (2) the one candidate-level rule, the penultimate-3' placement, is not
#       expressible by the shared AS-primer engine: [design_as_primers()] anchors
#       the *terminal* base on the variant (KASP geometry); `stab_at` only adds a
#       deliberate destabilising mismatch, it does not move the variant.  ASQ's
#       100-200 bp amplicon window and 62-64 C AS-primer window would also
#       contradict the KASP acceptance layer (`KASP_LIMITS`), which this build
#       grades ASQ with.
#  So `asq` stays out of `per_locus_schemes()` and `kasp_design_at()` (no
#  `design_asq_marker()` engine).  With no engine no *candidate* is ever scored on
#  this stage, so it declares **no gate at all** (C7: a gate no code path computes
#  is decorative).  The decision schema therefore renders the ASQ stage as
#  documentation - its declared ranking tier, no gate / decision / DOOMED node -
#  which needs `kasp_decision_schema()` to tolerate `gates = character(0)`
#  (R/kasp_decision.R).
# =============================================================================

# ---- the ASQ selection stage (documentation-only) ---------------------------
# No gate: none is computable from a candidate without an ASQ engine.  The one
# rule the method adds - the variant at the penultimate (or antepenultimate) 3'
# base (Kalendar 2022) - is documented in the header, not gated.  No relaxation
# ladder either (C6's relaxation order is empty): there is nothing to relax.
KASP_GATES_ASQ <- character(0)
KASP_RELAX_ASQ <- character(0)
# The one ranking dimension ASQ *would* add: the AS primer Tm sits near its 63 C
# optimum (Kalendar 2022: 62-64 C).  It is a margin, not a restatement of a gate
# (C14); with no engine it is declared for documentation only.
KASP_TIERS_ASQ <- c(tier1 = "AS primer Tm (toward the 63 C ASQ optimum)")

# This scheme contributes no gate, so it has no entry in the shared gate-label
# lookups (assembled in R/kasp_select.R); the empty vectors keep those
# concatenations valid without a dangling label for a gate that no longer exists.
KASP_GATE_LABEL_ASQ <- character(0)
KASP_GATE_SHORT_ASQ <- character(0)
KASP_GATE_LIMIT_ASQ <- character(0)
KASP_GATE_HARD_LABEL_ASQ <- character(0)
