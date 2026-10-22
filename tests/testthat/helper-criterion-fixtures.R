# The falsifiable-fixture ledger (ZW1): every criterion the registry can publish
# must be able to FAIL, and the evidence for that failure lives in one place.
#
# This helper is the ONE source of truth for that ledger.  It is auto-loaded by
# testthat for `test-falsifiable.R` and sourced by `dev/check-criteria.R`, so the
# pull-request gate and the test apply the same rule to the same table (a new
# criterion id fails both, not just the test).
#
# The ledger is a claim about REAL tests in this repo, so it is not trusted: the
# test greps each cited file for the quoted `test_that()` header and fails if the
# citation no longer exists.

# A criterion id that is exempt from the fail-fixture proof, keyed by the id,
# valued by WHY the in-memory suite cannot drive it to `fail` (a genome, a
# BLAST database or a population the suite does not carry).  Every value must be
# a non-empty reason: adding a key is the one place this suite admits it does
# not test a criterion.  Empty is the healthy state.
zw1_exemptions <- function() {
  reason <- paste(
    "Kwok 1990's classes are a RANKING of extension kinetics, so no class is a",
    "fail; an unclassed substitution pair is reported `unknown`, not `fail`,",
    "so the in-memory suite cannot drive this criterion to `fail`"
  )
  stats::setNames(reason, "as_mismatch_type")
}

# One `id -> test` fixture row.  `file` is a test file under tests/testthat and
# `test` its exact `test_that()` header (including the leading `::`).  The short
# constructor keeps the table readable: the citations are long, the rows are not.
zw1_fixture <- function(id, file, test) {
  tibble::tibble(id = id, test = paste0("tests/testthat/", file, test))
}

# The real test that makes each criterion FAIL: one row per non-exempt id.
zw1_fail_fixtures <- function() {
  kasp <- "test-kasp-checks.R"
  dec <- "test-decision.R"
  caps <- "test-caps.R"
  lamp <- "test-lamp.R"
  rpa <- "test-rpa.R"
  crispr <- "test-crispr.R"
  spec <- "test-specificity.R"
  taq <- "test-taqman.R"
  snap <- "test-snapshot.R"
  san <- "test-sanger.R"
  ssr <- "test-ssr.R"
  tech <- "test-technology.R"
  dup <- "test-duplex.R"
  arms <- "test-arms.R"
  hrm <- "test-hrm.R"
  cast <- "test-castpcr.R"
  dpcr <- "test-dpcr.R"
  isz <- "test-indel-size.R"
  lax <- "test-lax.R"
  out <- list(
    zw1_fixture(
      "site_class", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "site_alt_alleles", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "site_linked_skipped", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "site_alt_covered", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "site_linked_dependent", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "amplicon", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "as_tm", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "com_tm", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "com_tm_ge_as_tm", kasp,
      "::COM Tm must not fall below the AS primers (EasyKASP)"
    ),
    zw1_fixture(
      "as_gc", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "com_gc", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "max_run", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "ssr", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "inverted_repeat", kasp,
      "::the inverted-repeat screen fails a self-symmetric primer"
    ),
    zw1_fixture(
      "hairpin_dg", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "cross_dimer_dg", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "homodimer_dg", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "as_n_3p", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "as_n_var", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "com_n_var", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "as_dropout", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "as_foot_dropout", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "com_dropout", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    # `as_mismatch_type` is EXEMPT, not fixture-backed: Kwok's classes rank
    # kinetics, so no class is a fail (see `zw1_exemptions()`).
    zw1_fixture(
      "discrimination", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "in_silico_pcr", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "offtarget_hits", kasp,
      "::each broken criterion flips the verdict to fail"
    ),
    zw1_fixture(
      "as:gate:gate_discrimination", dec,
      "::an assay's outcome is stamped on the gates"
    ),
    zw1_fixture(
      "as:gate:gate_wobble", dec,
      "::an assay's outcome is stamped on the gates"
    ),
    zw1_fixture(
      "as:gate:gate_as_footprint", dec,
      "::an assay's outcome is stamped on the gates"
    ),
    zw1_fixture(
      "as:gate:gate_as_window", dec,
      "::an assay's outcome is stamped on the gates"
    ),
    zw1_fixture(
      "pair:gate:gate_structure", dec,
      "::an assay's outcome is stamped on the gates"
    ),
    zw1_fixture(
      "pair:gate:gate_tm", dec,
      "::an assay's outcome is stamped on the gates"
    ),
    zw1_fixture(
      "pair:gate:gate_amplicon", dec,
      "::an assay's outcome is stamped on the gates"
    ),
    zw1_fixture(
      "pair:gate:gate_com_footprint", dec,
      "::an assay's outcome is stamped on the gates"
    ),
    zw1_fixture(
      "caps:gate:gate_caps_site", dec,
      "::a discovered scheme's own gate nodes are stamped too (C13)"
    ),
    zw1_fixture(
      "caps:gate:gate_caps_resolvable", dec,
      "::a discovered scheme's own gate nodes are stamped too (C13)"
    ),
    zw1_fixture(
      "caps:gate:gate_caps_amplicon", dec,
      "::a discovered scheme's own gate nodes are stamped too (C13)"
    ),
    zw1_fixture(
      "caps:gate:gate_caps_tm", dec,
      "::a discovered scheme's own gate nodes are stamped too (C13)"
    ),
    zw1_fixture(
      "family_caps_site_defined", caps,
      "::find_caps_enzymes flags the enzyme the variant turns on or off"
    ),
    zw1_fixture(
      "family_caps_site_separates", caps,
      "::a resolvable CAPS call needs a site change, not just a length shift"
    ),
    zw1_fixture(
      "family_caps_band_gap", caps,
      "::a resolvable CAPS call needs a site change, not just a length shift"
    ),
    zw1_fixture(
      "family_caps_digest_spans", caps,
      "::find_caps_enzymes flags the enzyme the variant turns on or off"
    ),
    zw1_fixture(
      "family_caps_site_conserved", tech,
      "::a monomorphic site fails, whatever the marker's schema"
    ),
    zw1_fixture(
      "family_lamp_fip_alt", lamp,
      "::B3 is a minus-strand oligo and FIP carries the ALT base at its 3' end"
    ),
    zw1_fixture(
      "family_lamp_inner_tm_balance", lamp,
      "::lamp_geometry names a failed loop / 3'-end geometry"
    ),
    zw1_fixture(
      "family_lamp_inner_inside", lamp,
      "::lamp_geometry flags an F2-B2 span above the amplicon ceiling"
    ),
    zw1_fixture(
      "family_lamp_fip_variant", lamp,
      "::the LAMP variant gate re-checks FIP placement (not a decorative TRUE)"
    ),
    zw1_fixture(
      "family_rpa_primer_length", rpa,
      "::the RPA variant gate checks the actual placement"
    ),
    zw1_fixture(
      "family_rpa_thf_off", rpa,
      "::the probe gate records the relaxed THF/FQ rung honestly (C6/C7)"
    ),
    zw1_fixture(
      "family_rpa_probe_spans", rpa,
      "::the exo probe lies inside the amplicon between the primer 3' ends"
    ),
    zw1_fixture(
      "family_rpa_variant_under", rpa,
      "::design_rpa_marker declines non-substitution alleles"
    ),
    zw1_fixture(
      "family_crispr_pam", crispr,
      "::design_crispr_marker builds a Cas13a guide with a 3' non-G PFS"
    ),
    zw1_fixture(
      "family_crispr_seed", crispr,
      "::the CRISPR seed gate re-checks the guide's discriminating region"
    ),
    zw1_fixture(
      "family_crispr_guide_spans", crispr,
      "::design_crispr_marker informs and returns NULL when no guide passes"
    ),
    zw1_fixture(
      "family_crispr_conserved", tech,
      "::a monomorphic site fails, whatever the marker's schema"
    ),
    zw1_fixture(
      "family_crispr_offtarget", spec,
      "::a verdict fails when the products are not specific"
    ),
    zw1_fixture(
      "family_probe_spans", taq,
      "::probes lie strictly between the primer 3' ends, not merely in the amplicon"
    ),
    zw1_fixture(
      "family_probe_tm_delta", taq,
      "::the mismatch dG gate is the real discrimination, not only centrality"
    ),
    zw1_fixture(
      "family_probe_asp_discrimination", spec,
      "::a verdict fails when the products are not specific"
    ),
    zw1_fixture(
      "family_probe_extension", snap,
      "::the extension primer's 3' end sits one base 5' of the variant"
    ),
    zw1_fixture(
      "family_sanger_read_distance", san,
      "::a variant beyond the usable read is warned about, not silently paired"
    ),
    zw1_fixture(
      "family_sanger_primers_inside", san,
      "::a read primer falls back to the flanking primer when no nested one qualifies"
    ),
    zw1_fixture(
      "family_sanger_pair_defined", san,
      "::design_sanger_marker declines when no pair keeps the variant off the read start"
    ),
    zw1_fixture(
      "family_size_repeat_unit", ssr,
      "::design_ssr_marker declines when no repeat is in the window"
    ),
    zw1_fixture(
      "family_size_repeat_copies", ssr,
      "::the polymorphism gate is really evaluated and relaxes to the detection floor (C7)"
    ),
    zw1_fixture(
      "family_size_delta", isz,
      "::a sub-resolvable size difference is declined with a message"
    ),
    zw1_fixture(
      "family_size_fragment", ssr,
      "::the polymorphism gate is really evaluated and relaxes to the detection floor (C7)"
    ),
    zw1_fixture(
      "duplex:gate:gate_duplex_targets", dup,
      "::a hard duplex gate dooms (C4); a soft one relaxes within a bounded floor"
    ),
    zw1_fixture(
      "duplex:gate:gate_duplex_complement", dup,
      "::duplex_compatibility flags inter-locus complementarity and the topology"
    ),
    zw1_fixture(
      "duplex:gate:gate_duplex_tm", dup,
      "::the duplex stage computes its Tm at the duplex scheme's conditions (Z3)"
    ),
    zw1_fixture(
      "duplex:gate:gate_duplex_selfcomp", dup,
      "::the 3' fold-back fails gate_duplex_selfcomp inside the five-base cap"
    ),
    zw1_fixture(
      "arms:gate:gate_arms_disc", arms,
      "::design_arms_marker evaluates the arms gates and gives a resolvable gel"
    ),
    zw1_fixture(
      "arms:gate:gate_arms_bands", arms,
      "::design_arms_marker evaluates the arms gates and gives a resolvable gel"
    ),
    zw1_fixture(
      "arms:gate:gate_arms_outer", arms,
      "::an outer flanking pair sits outside the region and codes a product"
    ),
    zw1_fixture(
      "hrm:gate:gate_hrm_amplicon", hrm,
      "::the amplicon gate is hard and the Tm balance is the ladder's only rung"
    ),
    zw1_fixture(
      "hrm:gate:gate_hrm_tm", hrm,
      "::the HRM ladder widens the Tm window and reports it"
    ),
    zw1_fixture(
      "castpcr:gate:gate_castpcr_as3", cast,
      "::design_castpcr_marker detects the MUTANT allele, with a WT blocker and a locus probe"
    ),
    zw1_fixture(
      "castpcr:gate:gate_castpcr_mismatch", cast,
      "::the castpcr DAG gates the mutant anchor, the WT blocker and the probe ordering"
    ),
    zw1_fixture(
      "castpcr:gate:gate_castpcr_blocker", cast,
      "::the castpcr DAG gates the mutant anchor, the WT blocker and the probe ordering"
    ),
    zw1_fixture(
      "castpcr:gate:gate_castpcr_as_tm", cast,
      "::castPCR thermo is evaluated at its own reaction conditions (Z3)"
    ),
    zw1_fixture(
      "castpcr:gate:gate_castpcr_probe_tm", cast,
      "::castPCR thermo is evaluated at its own reaction conditions (Z3)"
    ),
    zw1_fixture(
      "castpcr:gate:gate_castpcr_amplicon", cast,
      "::castPCR thermo is evaluated at its own reaction conditions (Z3)"
    ),
    zw1_fixture(
      "crispr:gate:gate_crispr_pam", crispr,
      "::Cas13a rejects a G at the 3' flank of the plus-strand protospacer"
    ),
    zw1_fixture(
      "crispr:gate:gate_crispr_seed", crispr,
      "::the CRISPR seed gate re-checks the guide's discriminating region"
    ),
    zw1_fixture(
      "crispr:gate:gate_crispr_spacer", crispr,
      "::design_crispr_marker informs and returns NULL when no guide passes"
    ),
    zw1_fixture(
      "crispr:gate:gate_crispr_specificity", crispr,
      "::guide specificity is a hard gate"
    ),
    zw1_fixture(
      "crispr:gate:gate_crispr_gc", crispr,
      "::design_crispr_marker informs and returns NULL when no guide passes"
    ),
    zw1_fixture(
      "crispr:gate:gate_crispr_amplicon", crispr,
      "::design_crispr_marker builds a Cas12a guide with a TTTV PAM and an amplicon"
    ),
    zw1_fixture(
      "dpcr:gate:gate_dpcr_probe_comp", dpcr,
      "::every declared dPCR gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "dpcr:gate:gate_dpcr_probe_specific", dpcr,
      "::an ALT probe with a second perfect genomic site fails the dPCR specificity gate"
    ),
    zw1_fixture(
      "dpcr:gate:gate_dpcr_probe_structure", dpcr,
      "::every declared dPCR gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "dpcr:gate:gate_dpcr_probe_tm", dpcr,
      "::the dPCR DAG shows the gates in application order and the ladder least grave first"
    ),
    zw1_fixture(
      "dpcr:gate:gate_dpcr_primer_tm", dpcr,
      "::every declared dPCR gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "dpcr:gate:gate_dpcr_amplicon", dpcr,
      "::the amplicon soft gate relaxes to its floor instead of silently failing (C6/C7)"
    ),
    zw1_fixture(
      "indel_size:gate:gate_indel_size_gap", isz,
      "::the resolvable-gap gate is evaluated and relaxes to the capillary floor (C7)"
    ),
    zw1_fixture(
      "indel_size:gate:gate_indel_size_amplicon", isz,
      "::a sub-resolvable size difference is declined with a message"
    ),
    zw1_fixture(
      "indel_size:gate:gate_indel_size_tm", isz,
      "::scheme_def_indel_size declares the InDel-size gates, tiers and labels"
    ),
    zw1_fixture(
      "indel_size:gate:gate_indel_size_structure", isz,
      "::scheme_def_indel_size declares the InDel-size gates, tiers and labels"
    ),
    zw1_fixture(
      "indel_size:gate:gate_indel_size_unique", isz,
      "::a duplicated locus dooms the design through the hard uniqueness gate (C7)"
    ),
    zw1_fixture(
      "lamp:gate:gate_lamp_variant", lamp,
      "::the LAMP variant gate re-checks FIP placement (not a decorative TRUE)"
    ),
    zw1_fixture(
      "lamp:gate:gate_lamp_inner_tm", lamp,
      "::lamp_geometry names a failed loop / 3'-end geometry"
    ),
    zw1_fixture(
      "lamp:gate:gate_lamp_outer_tm", lamp,
      "::design_lamp_marker surfaces failed inner-primer geometry"
    ),
    zw1_fixture(
      "lamp:gate:gate_lamp_structure", lamp,
      "::design_lamp_marker surfaces failed inner-primer geometry"
    ),
    zw1_fixture(
      "lamp:gate:gate_lamp_cross_dimer", lamp,
      "::LAMP cross-dimer dG is evaluated at the LAMP temperature, not 37 C (Z3)"
    ),
    zw1_fixture(
      "lamp:gate:gate_lamp_amplicon", lamp,
      "::lamp_geometry flags an F2-B2 span above the amplicon ceiling"
    ),
    zw1_fixture(
      "lax:gate:gate_lax_coverage", lax,
      "::the lax stage dooms on its hard gate and on an exhausted floor"
    ),
    zw1_fixture(
      "lax:gate:gate_lax_alt_specific", lax,
      "::lax_ref_leakage catches a linked allele the wild type also carries"
    ),
    zw1_fixture(
      "lax:gate:gate_lax_cross_priming", lax,
      "::design_lax_marker evaluates the lax gates and records them (C7)"
    ),
    zw1_fixture(
      "rpa:gate:gate_rpa_variant", rpa,
      "::the RPA variant gate checks the actual placement"
    ),
    zw1_fixture(
      "rpa:gate:gate_rpa_primer_len", rpa,
      "::design_rpa_marker returns long primers and a short amplicon"
    ),
    zw1_fixture(
      "rpa:gate:gate_rpa_structure", rpa,
      "::RPA dimers are screened at the RPA temperature, not 37 C (Z3)"
    ),
    zw1_fixture(
      "rpa:gate:gate_rpa_probe", rpa,
      "::the probe gate records the relaxed THF/FQ rung honestly (C6/C7)"
    ),
    zw1_fixture(
      "rpa:gate:gate_rpa_amplicon", rpa,
      "::design_rpa_marker can add an exo probe that covers the variant"
    ),
    zw1_fixture(
      "sanger:gate:gate_sanger_dist", san,
      "::a variant beyond the usable read is warned about, not silently paired"
    ),
    zw1_fixture(
      "sanger:gate:gate_sanger_amplicon", san,
      "::the amplicon gate is really evaluated and relaxes to its floor (C7)"
    ),
    zw1_fixture(
      "sanger:gate:gate_sanger_tm", san,
      "::scheme_def_sanger describes the chemistry and its gates"
    ),
    zw1_fixture(
      "sanger:gate:gate_sanger_primer", san,
      "::a read primer falls back to the flanking primer when no nested one qualifies"
    ),
    zw1_fixture(
      "sanger:gate:gate_sanger_unique", san,
      "::a duplicated locus dooms the design through the hard uniqueness gate (C7)"
    ),
    zw1_fixture(
      "snapshot:gate:gate_snapshot_placement", snap,
      "::the extension primer's 3' end sits one base 5' of the variant"
    ),
    zw1_fixture(
      "snapshot:gate:gate_snapshot_comp", snap,
      "::scheme_def_snapshot describes the SNaPshot scheme and its gates"
    ),
    zw1_fixture(
      "snapshot:gate:gate_snapshot_unique", snap,
      "::every declared SNaPshot gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "snapshot:gate:gate_snapshot_structure", snap,
      "::every declared SNaPshot gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "snapshot:gate:gate_snapshot_tm", snap,
      "::the extension-Tm soft gate relaxes to its floor, never past it (C8)"
    ),
    zw1_fixture(
      "snapshot:gate:gate_snapshot_amplicon", snap,
      "::the flanking pair gives a short amplicon enclosing the variant"
    ),
    zw1_fixture(
      "ssr:gate:gate_ssr_flanks", ssr,
      "::a duplicated locus dooms the design through the hard flank-uniqueness gate (C7)"
    ),
    zw1_fixture(
      "ssr:gate:gate_ssr_polymorphic", ssr,
      "::the polymorphism gate is really evaluated and relaxes to the detection floor (C7)"
    ),
    zw1_fixture(
      "ssr:gate:gate_ssr_amplicon", ssr,
      "::design_ssr_marker declines when no repeat is in the window"
    ),
    zw1_fixture(
      "ssr:gate:gate_ssr_tm", ssr,
      "::the ssr scheme definition carries its gates, tiers and labels"
    ),
    zw1_fixture(
      "ssr:gate:gate_ssr_offtarget", ssr,
      "::design_ssr_marker returns a clean flanking pair outside the repeat"
    ),
    zw1_fixture(
      "taqman:gate:gate_taqman_probe_comp", taq,
      "::every declared TaqMan gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "taqman:gate:gate_taqman_probe_specific", taq,
      "::an ALT probe with a second perfect genomic site fails the specificity gate"
    ),
    zw1_fixture(
      "taqman:gate:gate_taqman_probe_structure", taq,
      "::every declared TaqMan gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "taqman:gate:gate_taqman_probe_tm", taq,
      "::every declared TaqMan gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "taqman:gate:gate_taqman_primer_tm", taq,
      "::every declared TaqMan gate is evaluated by the designer (C7)"
    ),
    zw1_fixture(
      "taqman:gate:gate_taqman_discrimination", taq,
      "::the mismatch dG gate is the real discrimination, not only centrality"
    ),
    zw1_fixture(
      "taqman:gate:gate_taqman_centrality", taq,
      "::the centrality soft gate relaxes to its interior floor, never past it (C8)"
    ),
    zw1_fixture(
      "taqman:gate:gate_taqman_amplicon", taq,
      "::the TaqMan DAG shows the gates in application order, ordering only as a tier"
    )
  )
  do.call(rbind, out)
}

# The reachability citations: `unknown` (missing data is not a pass) and
# `not_applicable` (another marker class) must both be reachable and tested.
zw1_status_citations <- function() {
  tibble::tibble(
    status = c("unknown", "not_applicable"),
    test = c(
      "tests/testthat/test-kasp-checks.R::check statuses keep pass/fail, and missing data is unknown not pass",
      "tests/testthat/test-kasp-checks.R::family checks are not_applicable to another marker class"
    )
  )
}

# The registry ids this ledger must cover: every `gate` and `check` row, walked
# (not a hand-kept list), so a new criterion appears here with no edit.  `tier`
# rows are ranking weights, not pass/fail criteria.
zw1_registered_ids <- function(limits = KASP_LIMITS) {
  spec <- criterion_spec(limits)
  sort(unique(spec$id[spec$role %in% c("gate", "check")]))
}

# The ledger: what the registry publishes (`registered`) and the evidence the
# repo offers for it (`fixtures`, `exempt`).  Nothing derived is stored here:
# coverage and staleness are computed from these three fields at CALL time
# (`zw1_ledger_coverage()`), so editing `$fixtures` or `$exempt` can never leave
# a cached answer behind that quietly accepts a dropped row.
zw1_ledger <- function(limits = KASP_LIMITS) {
  list(
    registered = zw1_registered_ids(limits),
    fixtures = zw1_fail_fixtures(),
    exempt = zw1_exemptions()
  )
}

# Derive coverage from what the ledger says right now.  `registered` is a
# parameter (defaulting to the ledger's own id list) so a ledger carrying its
# own set of published ids is still judged against those ids.
zw1_ledger_coverage <- function(ledger, registered = ledger$registered) {
  cited <- c(ledger$fixtures$id, names(ledger$exempt))
  list(
    uncovered = sort(setdiff(registered, cited)),
    # an exemption of a criterion that no longer exists is stale, not evidence
    stale_exemptions = sort(setdiff(names(ledger$exempt), registered)),
    # a fixture row for a criterion that no longer exists is stale too
    stale_fixtures = sort(setdiff(ledger$fixtures$id, registered))
  )
}

# Step 1 of the rule: the ledger is well-formed on its own terms.
#   - every registration is exempt with a non-empty reason, or backed by exactly
#     one non-empty `tests/testthat/<file>.R::<test>` string
#   - exemptions and fixtures name registered ids (a stale row fails)
zw1_ledger_problems <- function(ledger, registered = ledger$registered) {
  problems <- character(0)
  coverage <- zw1_ledger_coverage(ledger, registered)
  ex <- ledger$exempt
  if (length(ex) && (is.null(names(ex)) || anyNA(names(ex)) || any(!nzchar(names(ex))))) {
    problems <- c(problems, "an exemption has no criterion id (a named vector is required)")
  }
  blank <- names(ex)[is.na(ex) | !nzchar(trimws(ex))]
  if (length(blank)) {
    problems <- c(problems, paste0("an exemption has an empty reason: ", paste(blank, collapse = ", ")))
  }
  if (anyDuplicated(names(ex))) {
    problems <- c(problems, paste0(
      "a criterion is exempt twice: ", paste(unique(names(ex)[duplicated(names(ex))]), collapse = ", ")
    ))
  }
  if (length(coverage$stale_exemptions)) {
    problems <- c(problems, paste0(
      "a stale exemption names an id the registry no longer has: ",
      paste(coverage$stale_exemptions, collapse = ", ")
    ))
  }
  fx <- ledger$fixtures
  bad <- fx$id[is.na(fx$test) | !nzchar(fx$test) | !grepl("^tests/testthat/.+[.]R::.+$", fx$test)]
  if (length(bad)) {
    problems <- c(problems, paste0(
      "a fixture is not a `tests/testthat/<file>.R::<test>` citation: ", paste(bad, collapse = ", ")
    ))
  }
  if (anyDuplicated(fx$id)) {
    problems <- c(problems, paste0(
      "a criterion has more than one fixture row: ", paste(unique(fx$id[duplicated(fx$id)]), collapse = ", ")
    ))
  }
  if (length(coverage$stale_fixtures)) {
    problems <- c(problems, paste0(
      "a stale fixture names an id the registry no longer has: ",
      paste(coverage$stale_fixtures, collapse = ", ")
    ))
  }
  both <- intersect(names(ex), fx$id)
  if (length(both)) {
    problems <- c(problems, paste0(
      "a criterion is both exempt and fixture-backed: ", paste(both, collapse = ", ")
    ))
  }
  problems
}

# Step 2 of the rule: every cited test really exists in the cited file.  A
# renamed or deleted test is a failing claim, not a silently satisfied one.  The
# package root is found upwards, so this works from `tests/testthat/` (where
# testthat runs) and from the root (`dev/check-criteria.R`).
zw1_package_root <- function() {
  path <- normalizePath(".", winslash = "/")
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      return(normalizePath(".", winslash = "/"))
    }
    path <- parent
  }
}

zw1_citation_exists <- function(citation) {
  parts <- strsplit(citation, "::", fixed = TRUE)[[1L]]
  test <- paste(parts[-1L], collapse = "::")
  file <- file.path(zw1_package_root(), parts[[1L]])
  if (!file.exists(file)) {
    return(FALSE)
  }
  lines <- trimws(readLines(file, warn = FALSE))
  headers <- regmatches(lines, regexpr('^test_that[(]".*"', lines))
  headers <- gsub('^test_that[(]"|"$', "", headers)
  test %in% headers
}

# Step 3 of the rule: no registered id is left uncovered.  Coverage is derived
# from the ledger as passed in, never read from a stored field, so a removed
# fixture row is reported here.  The message names every uncovered criterion
# with its label and says exactly how to satisfy it.
zw1_uncovered_problems <- function(ledger, registered = ledger$registered, limits = KASP_LIMITS) {
  uncovered <- zw1_ledger_coverage(ledger, registered)$uncovered
  if (!length(uncovered)) {
    return(character(0))
  }
  spec <- criterion_spec(limits)
  label <- spec$label[match(uncovered, spec$id)]
  paste0(
    "criterion ", uncovered, " (", label, ") has no fail fixture: ",
    "add a fail fixture to tests/testthat/test-falsifiable.R or an exemption with a reason"
  )
}

# The whole rule, as one vector of problems (empty means the ledger is honest).
# `dev/check-criteria.R` and `test-falsifiable.R` both call this, so the PR gate
# and the test cannot disagree about what "covered" means.
zw1_ledger_violations <- function(ledger = zw1_ledger(), registered = ledger$registered) {
  cited <- unique(c(ledger$fixtures$test, zw1_status_citations()$test))
  missing <- cited[!vapply(cited, zw1_citation_exists, logical(1))]
  c(
    zw1_ledger_problems(ledger, registered),
    zw1_uncovered_problems(ledger, registered),
    if (length(missing)) {
      paste0(
        "a cited test does not exist (renamed or removed): ", paste(missing, collapse = "; "),
        " - a fixture must name a real `test_that()` header"
      )
    }
  )
}
