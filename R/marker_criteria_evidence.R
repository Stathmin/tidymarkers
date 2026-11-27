# =============================================================================
#  marker_criteria_evidence.R  --  tier targets, sources, evidence, asymmetry
# =============================================================================
#  Split from marker_criteria.R (roadmap ZT T7); the registry declarations stay
#  in marker_criteria.R, the table/vctrs/spec build in marker_criteria_table.R.
# =============================================================================
criterion_tier_targets <- function() {
  c(
    "as:tier:tier1" = "mutant/WT discrimination margin (computed in-silico)",
    "as:tier:tier2" = "AS primer Tm 58-60 C (package bar 58-63 C)",
    "as:tier:tier3" = "3'-terminal mismatch weakness (Kwok rank)",
    "as:tier:tier4" = "footprint load (simulated dropout)",
    "pair:tier:tier1" = "hairpin + homodimer dG >= -3 kcal/mol",
    "pair:tier:tier2" = "common-primer Tm >= AS Tm, ideally +2 C",
    "pair:tier:tier3" = "amplicon 50-100 bp",
    "pair:tier:tier4" = "footprint load (simulated dropout)",
    "pair:tier:tier5" = "weakest COM 3'-terminal 5-mer dG (package rule)",
    "duplex:tier:tier1" = "<= 5 self-complementary bases",
    "duplex:tier:tier2" = "shared amplicon <= 150 bp",
    "asq:tier:tier1" = "AS primer Tm 63 C (range 62-64 C)",
    "caps:tier:tier1" = ">= 10 bp band gap (gel resolvability)",
    "caps:tier:tier2" = "primer pair dTm <= 2 C, Tm 58-60 C",
    "caps:tier:tier3" = "CAPS product 200-1000 bp",
    "arms:tier:tier1" = "inner-primer discrimination margin (Kwok rank)",
    "arms:tier:tier2" = ">= 10 bp band separation",
    "arms:tier:tier3" = "primer Tm spread <= 2 C across the tube",
    "hrm:tier:tier1" = "variant at the centre of the amplified sequence",
    "hrm:tier:tier2" = "primer pair dTm <= 2 C",
    "castpcr:tier:tier1" = "ASP mutant/WT discrimination margin (Kwok rank)",
    "castpcr:tier:tier2" = "AS/probe Tm 58-60 C, dTm <= 2 C",
    "castpcr:tier:tier3" = "amplicon 80-120 nt",
    "crispr:tier:guide" = "Cas12a TTTV PAM; spacer 19-24 nt; GC 30-70 %",
    "crispr:tier:seed" = "variant inside the PAM-proximal seed window",
    "crispr:tier:specificity" = "guide off-target count (computed)",
    "crispr:tier:amplicon" = "pre-amplification product window (scheme default)",
    "dpcr:tier:tier1" = "probe Tm 5-10 C above the primer Tm",
    "dpcr:tier:tier2" = "probe hairpin/homodimer dG >= -3 kcal/mol",
    "dpcr:tier:tier3" = "amplicon 70-150 bp",
    "indel_size:tier:tier1" = "primer pair dTm <= 2 C, Tm 58-60 C",
    "indel_size:tier:tier2" = "amplicon wide enough to resolve the InDel (scheme default)",
    "indel_size:tier:tier3" = "<= 2 G/C in the last 5 nt; no T at the 3' end",
    "lamp:tier:tier1" = "variant at the FIP 3' end / BIP 5' end (mismatch position)",
    "lamp:tier:tier2" = "F2/B2 arm Tm 60-65 C",
    "lamp:tier:tier3" = "no secondary structure in the primers",
    "lamp:tier:tier4" = "target 130-260 nt, at most 300 nt",
    "lax:tier:tier1" = "resolvable ALT haplotype frequency (computed)",
    "rpa:tier:tier1" = "3'-anchor mismatch class (classification)",
    "rpa:tier:tier2" = "primer length 30-35 nt",
    "rpa:tier:tier3" = "THF probe: >= 30 nt from the 5' end, >= 15 nt from the 3' end",
    "rpa:tier:tier4" = "amplicon 100-200 bp",
    "sanger:tier:tier1" = "variant central in the amplicon",
    "sanger:tier:tier2" = "primer pair dTm <= 2 C",
    "sanger:tier:tier3" = "amplicon length (scheme default)",
    "snapshot:tier:tier1" = "extension-primer Tm 53-60 C",
    "snapshot:tier:tier2" = "extension primer hairpin/homodimer dG >= -3 kcal/mol",
    "snapshot:tier:tier3" = "amplicon near the ideal (scheme default)",
    "ssr:tier:tier1" = "primer pair dTm <= 2 C, Tm 58-60 C",
    "ssr:tier:tier2" = "amplicon 100-300 bp",
    "ssr:tier:tier3" = "<= 2 G/C in the last 5 nt; no T at the 3' end",
    "taqman:tier:tier1" = "probe-variant mismatch dG margin (package rule)",
    "taqman:tier:tier2" = "probe Tm 68-70 C, 8-10 C above the primers",
    "taqman:tier:tier3" = "amplicon 70-150 bp"
  )
}

# The rung behind each raised row.  `criterion_sources()` says WHERE the number
# is printed; this says HOW STRONGLY it is fixed, and only ever raises a row that
# has a source (a rung with no citation would be provenance without a page).
# Rows absent here keep the constructor's `rule_of_thumb`.
criterion_evidence_sources <- function() {
  c(
    as_gc = "primary_literature",
    com_gc = "primary_literature",
    "pair:gate:gate_amplicon" = "primary_literature",
    "hrm:gate:gate_hrm_amplicon" = "primary_literature",
    # the number is a vendor (TwistDx/Piepenburg) default reproduced in a paper
    max_run = "vendor_guideline",
    family_crispr_pam = "primary_literature",
    # ZW4.6 round 4 expansion: a peer-reviewed paper or protocol prints the
    # same quantity as the bar (see each source's `access_detail` for the
    # printed number and its offset from the package `bounds`).
    as_tm = "primary_literature",
    com_tm = "primary_literature",
    com_tm_ge_as_tm = "primary_literature",
    "pair:gate:gate_tm" = "primary_literature",
    "duplex:gate:gate_duplex_targets" = "primary_literature",
    "lamp:gate:gate_lamp_inner_tm" = "primary_literature",
    "lamp:gate:gate_lamp_amplicon" = "primary_literature",
    "lamp:gate:gate_lamp_structure" = "primary_literature",
    "lamp:gate:gate_lamp_variant" = "primary_literature",
    "crispr:gate:gate_crispr_spacer" = "primary_literature",
    "crispr:gate:gate_crispr_gc" = "primary_literature",
    "dpcr:gate:gate_dpcr_probe_comp" = "primary_literature",
    "dpcr:gate:gate_dpcr_primer_tm" = "primary_literature",
    "dpcr:gate:gate_dpcr_amplicon" = "primary_literature",
    "dpcr:gate:gate_dpcr_probe_tm" = "primary_literature",
    "taqman:gate:gate_taqman_probe_comp" = "primary_literature",
    "taqman:gate:gate_taqman_primer_tm" = "primary_literature",
    "taqman:gate:gate_taqman_amplicon" = "primary_literature",
    "castpcr:gate:gate_castpcr_as_tm" = "primary_literature",
    "castpcr:gate:gate_castpcr_probe_tm" = "primary_literature",
    "castpcr:gate:gate_castpcr_amplicon" = "primary_literature",
    family_probe_tm_delta = "primary_literature",
    "sanger:gate:gate_sanger_dist" = "primary_literature",
    family_sanger_read_distance = "primary_literature",
    "rpa:gate:gate_rpa_variant" = "primary_literature",
    inverted_repeat = "primary_literature",
    "caps:gate:gate_caps_resolvable" = "primary_literature",
    "arms:gate:gate_arms_bands" = "primary_literature",
    family_caps_band_gap = "primary_literature",
    family_size_delta = "primary_literature",
    "ssr:gate:gate_ssr_amplicon" = "primary_literature",
    "ssr:gate:gate_ssr_polymorphic" = "primary_literature",
    "indel_size:gate:gate_indel_size_amplicon" = "primary_literature",
    "indel_size:gate:gate_indel_size_gap" = "primary_literature",
    "hrm:gate:gate_hrm_tm" = "primary_literature",
    "indel_size:gate:gate_indel_size_tm" = "primary_literature",
    "ssr:gate:gate_ssr_tm" = "primary_literature",
    "caps:gate:gate_caps_tm" = "primary_literature",
    "sanger:gate:gate_sanger_tm" = "primary_literature",
    "duplex:gate:gate_duplex_tm" = "primary_literature",
    family_size_repeat_copies = "primary_literature",
    # a vendor default restated inside a paper
    "taqman:gate:gate_taqman_probe_tm" = "vendor_guideline",
    # ZW4.6 round 4 continued: the secondary-structure thresholds and the
    # primer-dimer / self-complementarity bars are all printed (see each
    # source's `access_detail` for the number and its offset).
    hairpin_dg = "primary_literature",
    cross_dimer_dg = "primary_literature",
    homodimer_dg = "primary_literature",
    "pair:gate:gate_structure" = "primary_literature",
    "indel_size:gate:gate_indel_size_structure" = "primary_literature",
    "dpcr:gate:gate_dpcr_probe_structure" = "primary_literature",
    "taqman:gate:gate_taqman_probe_structure" = "primary_literature",
    "snapshot:gate:gate_snapshot_structure" = "primary_literature",
    "rpa:gate:gate_rpa_structure" = "primary_literature",
    "duplex:gate:gate_duplex_complement" = "primary_literature",
    "duplex:gate:gate_duplex_selfcomp" = "primary_literature",
    "caps:gate:gate_caps_amplicon" = "primary_literature",
    "taqman:gate:gate_taqman_centrality" = "primary_literature",
    "rpa:gate:gate_rpa_amplicon" = "primary_literature",
    "rpa:gate:gate_rpa_probe" = "primary_literature",
    amplicon = "primary_literature",
    "snapshot:gate:gate_snapshot_comp" = "primary_literature",
    "crispr:gate:gate_crispr_amplicon" = "primary_literature",
    "snapshot:gate:gate_snapshot_placement" = "primary_literature",
    "snapshot:gate:gate_snapshot_tm" = "primary_literature",
    "indel_size:gate:gate_indel_size_unique" = "primary_literature",
    "taqman:gate:gate_taqman_probe_specific" = "primary_literature",
    "sanger:gate:gate_sanger_primer" = "primary_literature",
    "caps:gate:gate_caps_site" = "primary_literature",
    "ssr:gate:gate_ssr_offtarget" = "primary_literature",
    "ssr:gate:gate_ssr_flanks" = "primary_literature",
    family_lamp_fip_alt = "primary_literature",
    family_lamp_fip_variant = "primary_literature",
    family_caps_site_defined = "primary_literature",
    family_caps_site_separates = "primary_literature",
    family_caps_digest_spans = "primary_literature",
    family_caps_site_conserved = "primary_literature",
    family_lamp_inner_inside = "primary_literature",
    family_probe_spans = "primary_literature",
    # the review prints the 30-35 nt rule before noting PCR primers also work
    "rpa:gate:gate_rpa_primer_len" = "primary_literature",
    family_rpa_primer_length = "primary_literature",
    ssr = "primary_literature",
    offtarget_hits = "primary_literature",
    family_size_repeat_unit = "primary_literature",
    # ZW4.6 phases A/B: the target's rung.  Every tier above has a printed
    # reference point; the aggregation policy stays rule_of_thumb by
    # construction (CRITERION_TIER_POLICY_EVIDENCE), never promoted here.
    "as:tier:tier2" = "primary_literature",
    "pair:tier:tier1" = "primary_literature",
    "pair:tier:tier2" = "primary_literature",
    "pair:tier:tier3" = "primary_literature",
    "duplex:tier:tier1" = "primary_literature",
    "duplex:tier:tier2" = "primary_literature",
    "asq:tier:tier1" = "primary_literature",
    "caps:tier:tier1" = "primary_literature",
    "caps:tier:tier2" = "primary_literature",
    "caps:tier:tier3" = "primary_literature",
    "arms:tier:tier2" = "primary_literature",
    "arms:tier:tier3" = "primary_literature",
    "hrm:tier:tier1" = "primary_literature",
    "hrm:tier:tier2" = "primary_literature",
    "castpcr:tier:tier2" = "primary_literature",
    "castpcr:tier:tier3" = "primary_literature",
    "crispr:tier:guide" = "primary_literature",
    "dpcr:tier:tier1" = "primary_literature",
    "dpcr:tier:tier2" = "primary_literature",
    "dpcr:tier:tier3" = "primary_literature",
    "indel_size:tier:tier1" = "primary_literature",
    "indel_size:tier:tier3" = "primary_literature",
    "lamp:tier:tier2" = "primary_literature",
    "lamp:tier:tier3" = "primary_literature",
    "lamp:tier:tier4" = "primary_literature",
    "rpa:tier:tier2" = "primary_literature",
    "rpa:tier:tier3" = "primary_literature",
    "rpa:tier:tier4" = "primary_literature",
    "sanger:tier:tier2" = "primary_literature",
    "snapshot:tier:tier1" = "primary_literature",
    "snapshot:tier:tier2" = "primary_literature",
    "ssr:tier:tier1" = "primary_literature",
    "ssr:tier:tier2" = "primary_literature",
    "ssr:tier:tier3" = "primary_literature",
    "taqman:tier:tier2" = "vendor_guideline",
    "taqman:tier:tier3" = "primary_literature"
  )
}

# Where a metamorphic invariant deliberately does NOT hold, and why.  A
# criterion is absent when the rule holds for it; a row that carries a reason
# is one whose two swapped inputs are not interchangeable, so asserting
# symmetry over it would be asserting a falsehood (ZW0.1).
#
# Keys are criterion ids, values one sentence a reviewer can check against the
# design code - the reason names the structural fact that breaks the symmetry,
# never a bare "exception".
criterion_asymmetry <- function() {
  c(
    # KASP: `AS_ref` is built on the FAM tail and `AS_alt` on the HEX tail
    # (R/marker_assemble.R), while both sit at the same AS concentration
    # (R/marker_mix.R).  The gates grade the chosen allele-specific primer, not
    # the dye it lands on, so swapping REF/ALT changes which allele a fluorophore
    # reports even though neither primer changed.
    "as:gate:gate_discrimination" = paste0(
      "KASP fluorophore asymmetry: AS_ref carries the FAM tail and AS_alt the HEX ",
      "tail (marker_assemble.R), so swapping REF/ALT moves an allele to another ",
      "fluorophore though the AS recipe is identical."
    ),
    "as:gate:gate_wobble" = paste0(
      "KASP fluorophore asymmetry: the wobble gate grades the chosen AS primer ",
      "once, while its allele is pinned to FAM or HEX by the tail, so the pair is ",
      "not mirrorable across REF/ALT."
    ),
    "as:gate:gate_as_footprint" = paste0(
      "KASP fluorophore asymmetry: the footprint gate is judged on AS_ref/AS_alt ",
      "as placed on the FAM/HEX tails, so a REF/ALT swap re-assigns the graded ",
      "primer to the other dye slot."
    ),
    "as:gate:gate_as_window" = paste0(
      "KASP fluorophore asymmetry: the AS window gate reports the primer in its ",
      "fixed FAM or HEX slot, so the REF/ALT image of the probe set is a different ",
      "probe set, not the same one relabelled."
    ),
    # LAMP: FIP is built with the ALT allele at the 3' arm and its 3' end anchors
    # BIP (R/marker_lamp.R); `lamp_fuse()` also requires each inner primer's 5'
    # loop arm to be hotter than its 3' displacement arm.  FIP and BIP are
    # therefore neither interchangeable with each other nor mirrorable.
    "family_lamp_fip_alt" = paste0(
      "LAMP FIP/BIP asymmetry: only FIP reads the variant, with the ALT base at ",
      "the 3' terminus of its displacement arm, so BIP cannot host the same check."
    ),
    "family_lamp_inner_tm_balance" = paste0(
      "LAMP FIP/BIP asymmetry: each inner primer must have a 5' loop arm hotter ",
      "than its 3' displacement arm (lamp_fuse), so the FIP and BIP Tm patterns are ",
      "constrained separately and cannot be exchanged."
    ),
    "family_lamp_fip_variant" = paste0(
      "LAMP FIP/BIP asymmetry: the variant is a fixed feature of FIP's 3' arm, ",
      "and BIP is placed downstream of FIP (lamp_bip from fip$arm3), so the two ",
      "inner primers are positionally and functionally distinct."
    ),
    "lamp:gate:gate_lamp_variant" = paste0(
      "LAMP FIP/BIP asymmetry: discrimination is scored at the FIP 3' terminus ",
      "(or a loop primer), never at BIP, which lies downstream of the variant."
    ),
    "lamp:gate:gate_lamp_inner_tm" = paste0(
      "LAMP FIP/BIP asymmetry: FIP's 3' base is the ALT allele so it cannot always ",
      "be a G/C clamp, making the inner-primer Tm gate deliberately weaker for FIP ",
      "than for BIP."
    )
  )
}

# The numeric KASP-family bars, read from the limit object in force.
