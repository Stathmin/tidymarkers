# =============================================================================
#  marker_core.R  --  constants, geometry, cache, genome access
# =============================================================================
# =============================================================================
#  marker_functions.R  --  design + QC primitives for KASP assays
#  100% R; heavy thermodynamics is delegated to packages.
# =============================================================================
#
#  KASP = two allele-specific forward primers (FAM / HEX tails, identical except
#  for the 3'-terminal SNP base) plus one common reverse primer.  Taq cannot
#  extend a 3'-mismatched primer, so discrimination is competitive annealing.
#
#  Defaults follow LGC/KBioscience guidance and EasyKASP (Zhang, Yang & Wen,
#  BMC Bioinformatics 2025, 26:292) -- defaults, not dogma; validate in your lab.
#  The machine-checkable form is marker_design_checks()/KASP_LIMITS; the full
#  design manual is README section 3 and vignette("kasp-intro").
# =============================================================================


# --- constants ---------------------------------------------------------------

# Standard KASP 5' universal tails (added at synthesis; not used in Tm).
KASP_TAIL_FAM <- "GAAGGTGACCAAGTTCATGCT" # "FAM" channel
KASP_TAIL_HEX <- "GAAGGTCGGAGTCAACGGATT" # "HEX/VIC" channel

# Tm conditions (mM) for TmCalculator
PCR_COND <- list(Na = 50, Mg = 1.5, dNTPs = 0.6, DMSO = 0)
# TmCalculator renamed DNA_NN_Allawi_1998 to DNA_NN_Allawi_1997 in 1.1.3 and
# warns on the old spelling.  The parameters are unchanged; resolve the name
# the installed version knows so both 1.1.2 and >= 1.1.3 run without a warning.
NN_TABLE_ALLAWI <- local({
  choices <- eval(formals(TmCalculator::tm_nn)$nn_table)
  if ("DNA_NN_Allawi_1997" %in% choices) "DNA_NN_Allawi_1997" else "DNA_NN_Allawi_1998"
})
NN_TABLE <- "DNA_NN_SantaLucia_2004"
# The design Tm coordinate.  The KASP/AS/COM and dPCR windows are stated on this
# correction, so it stays the default until those windows are re-anchored
# together (see ZW3.4).  It is a historical GC-era monovalent correction
# applied to a nearest-neighbour model, not the NN-consistent one.  The
# NN-consistent reference model is "SantaLucia1998-2": the same SantaLucia 1998
# paper as the NN table, the pair Primer3/`oligotm` and Biopython apply.  It
# reproduces the printed soybean Primer3 panel to 0.08 C (ZW3.5) and is the
# model every external Tm comparison is stated in (ZW3.4).
TM_SALT_METHOD <- "Schildkraut2010"
# PCR conditions (molar) for DECIPHER::CalculateEfficiencyPCR
PCR_EFF <- list(temp = 60, P = 4e-7, ions = 0.2)

# Per-scheme reaction temperatures (C) for the dG engines, see scheme_conditions().
# PCR_COND holds only the salt / denaturant values; only `temperature` is
# scheme-specific here, and every value is one the code already assumes rather
# than a new buffer recipe:
#   pcr     the thermal-cycler chemistries anneal / read at the PCR
#           extension-anneal temperature pcr_efficiency() already passes to
#           DECIPHER (PCR_EFF$temp, 60 C)
#   lamp    the LAMP operating window recorded in R/marker_lamp.R:5 (60-65 C)
#   rpa     the RPA operating window recorded in R/marker_rpa.R (RPA_TEMP,
#           37-42 C), a crowded reaction this nearest-neighbour model does not
#           represent - its Tm is informational, never a gate
#   ambient a restriction digest / bead capture / capillary read with no heated
#           incubation; 25 C is an ASSUMPTION (room temperature), not a stated
#           protocol step
# The values are repeated here (not read from RPA_TEMP / the LAMP header, which
# collate later) and tests/testthat/test-oligo-dg.R pins them equal to those
# constants so the two cannot drift apart silently.
SCHEME_TEMPERATURE <- list(
  pcr     = c(PCR_EFF$temp, PCR_EFF$temp),
  lamp    = c(60, 65),
  rpa     = c(37, 42),
  ambient = c(25, 25)
)

# Scheme id -> family whose reaction conditions it uses.  Any id not listed is a
# PCR-style thermal-cycler chemistry.  `crispr` stays in the PCR family: only its
# readout is isothermal (collateral cleavage of an RPA-or-PCR amplicon), so 60 C
# is the model the rest of the family shares; it is separately named in
# TM_INFORMATIONAL below because its Tm must never gate.
SCHEME_TEMP_FAMILY <- c(
  snp = "pcr", indel = "pcr", multi = "pcr", duplex = "pcr", primer = "pcr",
  asq = "pcr", arms = "pcr", hrm = "pcr", castpcr = "pcr", crispr = "pcr",
  dpcr = "pcr", indel_size = "pcr", lax = "pcr", sanger = "pcr",
  taqman = "pcr",
  caps = "ambient",
  snapshot = "ambient",
  ssr = "ambient",
  lamp = "lamp",
  rpa = "rpa"
)

# Schemes whose reaction temperature this nearest-neighbour model does not
# represent, so their Tm is INFORMATIONAL and must never become a gate:
#   rpa    a crowded recombinase/SSB reaction (RPA_NO_TM in R/marker_rpa.R is an
#          open c(0, 200) on purpose - primer LENGTH is the gate, not Tm)
#   crispr an RPA pre-amplification feeds the Cas readout, so it inherits RPA's
#          open Tm window (R/marker_crispr.R screens its pre-pair with RPA_NO_TM)
TM_INFORMATIONAL <- c("rpa", "crispr")

# Buffer fields a scheme overrides, keyed by temperature family.  A field belongs
# here only when BOTH hold: an engine in this package models it (TmCalculator's
# tm_nn takes Na, K, Tris, Mg, dNTPs and DMSO; OligoArrayAux takes Na, Mg, dNTP
# and temperature) AND a cited source states the value.  The rows are empty on
# purpose: no scheme in the curated sources states a difference from PCR_COND,
# so the merge in scheme_conditions() is the seam, not an invented recipe.  A
# value added here moves every calibrated window (KASP 55-65 / 58-63 C, dPCR
# 44-50 C), so it is a ZW-scale re-anchor like changing PCR_COND itself.
SCHEME_BUFFER <- list(
  pcr     = list(),
  lamp    = list(),
  rpa     = list(),
  ambient = list()
)

# Additives a scheme actually uses that NO engine here can model, keyed by
# scheme id with the temperature families as the fallback.  Naming them on the
# returned object is the point: a Tm computed for that scheme is missing a term,
# never silently complete.
#   lamp    betaine is a common LAMP additive; it destabilises the duplex, so
#           the real Tm is lower than the computed one
#   rpa     the recombinase/SSB/crowding-agent buffer is not a PCR buffer
#   crispr  its pre-amplification is RPA, hence the same gap
# Neither TmCalculator nor OligoArrayAux takes a betaine or crowding term.
SCHEME_UNMODELLED <- list(
  pcr     = character(),
  lamp    = "betaine",
  rpa     = "crowding",
  ambient = character(),
  crispr  = "crowding"
)

# Ensembl Plants FASTA defline -> TAIR chromosome name; the only sequences kept.
ATH_CHROM_MAP <- c(
  "1" = "Chr1", "2" = "Chr2", "3" = "Chr3", "4" = "Chr4",
  "5" = "Chr5", "Mt" = "ChrM", "Pt" = "ChrC"
)

# Design targets / thresholds (named, so they are not magic numbers)
AS_TM_TARGET <- 58 # preferred AS-primer Tm (C); kept below the common primer
COM_TM_TARGET <- 60 # preferred common-primer Tm (C); KASP wants COM Tm >= AS Tm
HAIRPIN_DG_WARN <- -3 # flag hairpins with a reaction-temperature dG below this (kcal/mol)
SEED_K <- 15L # 3'-seed length for the fast specificity screen
SEED_MISMATCH <- 1L
# Largest pattern-hit set worth keeping in the per-session scan cache: a result
# bigger than this costs more memory than the rescan it would save.
HIT_CACHE_MAX <- 20000L
# Entry bound of that cache.  It must clear a whole design's working set: a
# 6-target design holds 357 distinct (pattern, chromosome, radius) answers, and
# a bound below that made the run evict and rescan its own answers (at 256, the
# last scans evicted the first and 26 of them were scanned twice).  The entries
# themselves are small - the design's 383 results held 42 hit ranges in total -
# so the bound is set for a multi-locus run rather than for memory.
HIT_CACHE_SLOTS <- 4096L

# Ensembl consequence wording -> design priority (lower = preferred)
# Sequence-Ontology term -> design priority (lower = preferred).  These are the
# canonical SO terms exactly as Ensembl/VEP spell them, so matching is by full
# token (not substring).  `SO_RANK_OTHER` sits one below the worst known term.
SO_TERM_RANK <- c(
  stop_gained              = 1L,
  stop_lost                = 1L,
  start_lost               = 1L,
  frameshift_variant       = 1L,
  splice_acceptor_variant  = 2L,
  splice_donor_variant     = 2L,
  missense_variant         = 3L,
  inframe_insertion        = 3L,
  inframe_deletion         = 3L,
  splice_region_variant    = 4L,
  synonymous_variant       = 5L,
  stop_retained_variant    = 5L,
  `5_prime_UTR_variant`    = 6L,
  `3_prime_UTR_variant`    = 6L,
  intron_variant           = 7L,
  upstream_gene_variant    = 8L,
  downstream_gene_variant  = 8L,
  intergenic_variant       = 9L
)
SO_RANK_OTHER <- max(SO_TERM_RANK) + 1L

# Ensembl VEP IMPACT ladder, used only when the variant table carries an `impact`
# column (the Ensembl `overlap/region` REST endpoint does not provide it).
IMPACT_RANK <- c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)

# KASP acceptance thresholds (consensus + EasyKASP 2025; the manual mirrors this)
KASP_LIMITS <- list(
  amplicon = c(50, 100), # bp  (EasyKASP: products < 100 bp)
  as_tm = c(55, 65), # C
  com_tm = c(58, 63), # C  (EasyKASP aims ~59.9-62.9)
  gc = c(0.40, 0.60), # fraction
  max_run = 4L, # longest homopolymer tolerated is max_run - 1
  hairpin_dg = -3, # hairpin dG at the reaction temperature (kcal/mol); must be >= this
  duplex_dg = -9, # cross-dimer dG at the reaction temperature; must be >= this
  homodimer_dg = -9, # self-dimer dG at the reaction temperature; must be >= this
  discrimination = 0.5, # Taq on- minus off-target efficiency
  discrimination_margin = 0.1, # near-miss band: [min_disc - margin, min_disc) reads as medium risk
  footprint_dropout = 0.10 # expected fraction of templates lost to footprint variants
)

# The AS-primer lengths the KASP engine designs (design_as_primers()'s default).
# The largest is the greatest reach a primer can have from the target - the
# half-width of the "primer window" haplotype_groups() keys on.
KASP_AS_LENGTHS <- 20:27L


# --- conditions / engines seam -----------------------------------------------

#' Reaction conditions for the thermodynamics engines
#'
#' One named object for the conditions a thermodynamic engine needs.  The units
#' are the ones the rest of the package already uses: `Na`, `Mg` and `dNTPs` in
#' mM, `DMSO` in percent, `temperature` in degrees C.  The defaults are exactly
#' the values used before this seam existed (`PCR_COND` plus the 37 C duplex
#' temperature of [oligo_dg()]), so a default call reproduces the package's
#' historical numbers.
#'
#' `temperature` is the reaction temperature used by the dG engines; the
#' nearest-neighbour Tm ([primer_tm()]) does not depend on it.
#'
#' @param Na,Mg,dNTPs monovalent salt, magnesium and dNTP concentrations (mM).
#' @param DMSO DMSO concentration (%).
#' @param temperature reaction temperature (C) for the dG engines.
#' @return a named list with elements `Na`, `Mg`, `dNTPs`, `DMSO`, `temperature`.
#' @examples
#' reaction_conditions()
#' reaction_conditions(Na = 100, temperature = 60)
#' @family engines
#' @export
reaction_conditions <- function(Na = PCR_COND$Na, Mg = PCR_COND$Mg, # nolint: object_name_linter.
                                dNTPs = PCR_COND$dNTPs, DMSO = PCR_COND$DMSO, # nolint: object_name_linter.
                                temperature = 37) {
  out <- list(Na = Na, Mg = Mg, dNTPs = dNTPs, DMSO = DMSO, temperature = temperature)
  for (nm in names(out)) {
    value <- out[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      abort_msg("`{nm}` must be a single finite number.")
    }
  }
  if (any(unlist(out[c("Na", "Mg", "dNTPs", "DMSO")]) < 0)) {
    abort_msg("Salt, dNTP and DMSO concentrations cannot be negative.")
  }
  out
}

#' Normalise a (possibly partial) reaction-conditions list
#'
#' Fills missing fields from [reaction_conditions()] and rejects unknown names,
#' so a typo cannot silently leave the old condition in place.
#' @keywords internal
#' @noRd
as_reaction_conditions <- function(conditions = reaction_conditions()) {
  if (is.null(conditions)) {
    return(reaction_conditions())
  }
  if (!is.list(conditions) || is.null(names(conditions)) || any(!nzchar(names(conditions)))) {
    abort_msg("`conditions` must be a named list, e.g. reaction_conditions().")
  }
  defaults <- reaction_conditions()
  unknown <- setdiff(names(conditions), names(defaults))
  if (length(unknown)) {
    abort_msg(
      "Unknown reaction condition(s): {.val {unknown}}.",
      i = "Allowed fields: {.val {names(defaults)}}."
    )
  }
  merged <- defaults
  for (nm in names(conditions)) merged[[nm]] <- conditions[[nm]]
  do.call(reaction_conditions, merged)
}

#' Reaction conditions for one marker scheme
#'
#' The single place a scheme id becomes a [reaction_conditions()] object, so the
#' dG engines can be evaluated at the temperature the chemistry actually runs at
#' instead of a global 37 C.  Two things are per-scheme: `temperature` (below)
#' and the buffer overrides in `SCHEME_BUFFER`.
#'
#' The buffer seam is deliberately thin.  A field is overridden only when an
#' engine here models it (TmCalculator's `tm_nn` takes Na, K, Tris, Mg, dNTPs and
#' DMSO; OligoArrayAux takes Na, Mg, dNTP and temperature) *and* a cited source
#' states the value.  No scheme in the curated sources states such a difference,
#' so every buffer field is `PCR_COND` today; the empty rows are the seam a
#' stated value plugs into.  Additives no engine can model are named on the
#' returned `unmodelled` attribute instead of being silently assumed absent.
#'
#' Encoded temperatures and their source:
#' * PCR-style chemistries (`snp`, `indel`, `multi`, `duplex`, `primer`, `asq`,
#'   `arms`, `hrm`, `castpcr`, `crispr`, `dpcr`, `indel_size`, `lax`, `sanger`,
#'   `taqman`) - `PCR_EFF$temp` (60 C), the anneal/extension temperature the
#'   package already uses for [pcr_efficiency()].
#' * `lamp` - 60-65 C, the LAMP operating window already recorded in
#'   `R/marker_lamp.R` (`LAMP_TEMP`).
#' * `rpa` - 37-42 C, the RPA operating window already recorded in
#'   `R/marker_rpa.R` (`RPA_TEMP`).  RPA buffers are crowded (recombinase, SSB,
#'   crowding agent), which this nearest-neighbour model does not represent: a
#'   PCR-style Tm is *informational* for RPA, never a gate (see [primer_tm()]; RPA
#'   primers are screened on length, with `RPA_NO_TM` deliberately open).
#'   `crispr` is mapped to the PCR family (its readout is what the family shares)
#'   but is likewise listed in `TM_INFORMATIONAL`: its isothermal pre-amplification
#'   is RPA and its pre-primers are already screened with `RPA_NO_TM`.
#' * `caps`, `snapshot`, `ssr` - no heated incubation is assumed by the code (a
#'   restriction digest, bead capture and capillary read are all at room
#'   temperature), so 25 C; this is an explicit *assumption*, not a stated step
#'   of any of those protocols.
#'
#' `temperature` is the top of each window (60 C, 65 C, 42 C, 25 C): a higher
#' temperature weakens every duplex, so it is the least punitive reading for a
#' dimer screen and the one that cannot hide a design under a gentler model.
#' The rest of the window is registered on the returned object's
#' `temperature_range` attribute.
#'
#' @param scheme a scheme id (see [marker_schemes()]) or a temperature family
#'   (`"pcr"`, `"lamp"`, `"rpa"`, `"ambient"`).  An unregistered id is treated as
#'   PCR-style, matching the package's historical behaviour.
#' @return a [reaction_conditions()] object with a `temperature_range` attribute
#'   (C, the method's operating window), a `tm_informational` attribute
#'   (`TRUE` for RPA and CRISPR, whose Tm must not gate) and an `unmodelled`
#'   character vector naming the scheme's additives that no engine here
#'   represents ("betaine" for LAMP, "crowding" for RPA and CRISPR, empty
#'   elsewhere).
#' @examples
#' scheme_conditions("lamp")
#' attr(scheme_conditions("rpa"), "temperature_range")
#' @seealso [reaction_conditions()], [primer_tm()], [dimer_thermo()]
#' @family engines
#' @export
scheme_conditions <- function(scheme = "pcr") {
  if (!is.character(scheme) || length(scheme) != 1L || is.na(scheme) || !nzchar(scheme)) {
    abort_msg("`scheme` must be a single scheme id.")
  }
  family <- if (scheme %in% names(SCHEME_TEMPERATURE)) {
    scheme
  } else if (scheme %in% names(SCHEME_TEMP_FAMILY)) {
    unname(SCHEME_TEMP_FAMILY[[scheme]])
  } else {
    "pcr"
  }
  window <- SCHEME_TEMPERATURE[[family]]
  cond <- do.call(
    reaction_conditions,
    c(list(temperature = window[2L]), SCHEME_BUFFER[[family]])
  )
  unmodelled <- if (scheme %in% names(SCHEME_UNMODELLED)) {
    SCHEME_UNMODELLED[[scheme]]
  } else {
    SCHEME_UNMODELLED[[family]]
  }
  structure(cond,
    temperature_range = window, scheme = scheme, family = family,
    tm_informational = scheme %in% TM_INFORMATIONAL, unmodelled = unmodelled
  )
}

#' Build one engine adapter
#' @keywords internal
#' @noRd
marker_engine <- function(name, capabilities, ...) {
  structure(list(name = name, capabilities = capabilities, ...), class = "marker_engine")
}

#' Marker engines and their capabilities
#'
#' The seam between the design code and the tools that compute the
#' thermodynamics and the genome-wide specificity.  Each adapter carries its
#' entry points (`tm`, `dimer`, `dg`; `pcr`, `seed`) and a capability table.
#' A capability the adapter does not declare is reported as `"unknown"` by
#' [capabilities()] - never `TRUE` - so an unverified feature cannot silently
#' pass a gate.
#'
#' @param thermo `"auto"` (the default, or `options(tidymarkers.thermo_engine)`)
#'   keeps OligoArrayAux as the `dg` engine when it is installed.  `"internal"`
#'   selects the pure-R adapters ([internal_duplex_dg], [internal_hairpin_dg]) so
#'   a machine without OligoArrayAux can still screen duplex dG and hairpin dG:
#'   `duplex_dg` and `hairpin_dg` are then `TRUE`, while `fold` stays `FALSE`.
#'   The hairpin adapter is a *screen*, not UNAFold: it scores a perfect stem plus
#'   one loop penalty, so it has no partition function and no mismatched, bulged
#'   or interior loops, and it predicts no structure.  On the panel in
#'   `tests/testthat/test-dg-internal.R` it is within mean 0.54 / max 1.48 kcal/mol
#'   of UNAFold at 37 C, and it under-reports that stability (its number is less
#'   negative), so the `hairpin_dg >= -3` gate reads more permissively than on
#'   UNAFold - never stricter.  `pcr_efficiency` needs UNAFold *through* DECIPHER,
#'   so it is `TRUE` only when OligoArrayAux is installed and DECIPHER loads - the
#'   internal adapter cannot supply it.
#'   Nothing falls back silently: the adapter is chosen by name or by option.
#' @return a named list (`thermo`, `specificity`) of `marker_engine` adapters,
#'   classed `marker_engines`.
#' @examples
#' engines <- marker_engines()
#' names(engines)
#' capabilities(engines$specificity)
#' capabilities(marker_engines(thermo = "internal")$thermo)
#' @family engines
#' @export
marker_engines <- function(thermo = getOption("tidymarkers.thermo_engine", "auto")) {
  thermo <- match.arg(thermo, c("auto", "internal"))
  internal <- identical(thermo, "internal")
  # probe once: the three capability flags below used to call
  # `has_oligoarrayaux()` separately and `oligoarrayaux_paths()` twice
  oaa <- !internal && has_oligoarrayaux()
  oaa_paths <- if (oaa) oligoarrayaux_paths() else NULL
  structure(
    list(
      thermo = marker_engine(
        "thermo",
        capabilities = list(
          tm = TRUE,
          salt = TRUE,
          temperature = TRUE,
          duplex_dg = internal || oaa,
          hairpin_dg = internal || oaa,
          fold = oaa && isTRUE(oaa_paths$ss),
          pcr_efficiency = oaa && requireNamespace("DECIPHER", quietly = TRUE)
        ),
        tm = primer_tm,
        dimer = dimer_thermo,
        dg = if (internal) internal_duplex_dg else oligo_dg
      ),
      specificity = marker_engine(
        "specificity",
        capabilities = list(
          in_silico_pcr = TRUE,
          seed_counts = TRUE,
          iupac = TRUE,
          mismatch = TRUE,
          bowtie = requireNamespace("Rbowtie", quietly = TRUE),
          blast = requireNamespace("rBLAST", quietly = TRUE) && has_blast(),
          fm = requireNamespace("Rsubread", quietly = TRUE) &&
            requireNamespace("Rsamtools", quietly = TRUE)
        ),
        pcr = in_silico_pcr,
        seed = specificity_seed
      )
    ),
    class = "marker_engines"
  )
}

#' Capabilities of a marker engine
#'
#' Queries the capability table of an adapter from [marker_engines()].  A
#' capability the adapter does not declare is `"unknown"`, not `TRUE`: the
#' absence of a declaration is not evidence that the engine can do it.
#'
#' @param engine a `marker_engine` adapter or the whole `marker_engines`
#'   registry.
#' @param capability a single capability name; `NULL` returns every declared
#'   capability.
#' @param ... reserved for methods.
#' @return with `capability = NULL`, a named list of declared capabilities (for
#'   the registry, one such list per engine).  With a capability, its value:
#'   `TRUE`, `FALSE`, or `"unknown"`.
#' @examples
#' engines <- marker_engines()
#' capabilities(engines$thermo, "tm")
#' capabilities(engines$thermo, "not_a_capability") # "unknown"
#' @family engines
#' @export
capabilities <- function(engine, capability = NULL, ...) {
  if (!is.null(capability) && (!is.character(capability) || length(capability) != 1L)) {
    abort_msg("`capability` must be a single capability name.")
  }
  UseMethod("capabilities")
}

#' @rdname capabilities
#' @export
capabilities.marker_engine <- function(engine, capability = NULL, ...) {
  if (is.null(capability)) {
    return(engine$capabilities)
  }
  value <- engine$capabilities[[capability]]
  if (is.null(value)) "unknown" else value
}

#' @keywords internal
#' @noRd
#' @exportS3Method
capabilities.default <- function(engine, capability = NULL, ...) {
  if (is.null(capability)) list() else "unknown"
}

#' @rdname capabilities
#' @export
capabilities.marker_engines <- function(engine, capability = NULL, ...) {
  if (is.null(capability)) {
    return(lapply(engine, capabilities))
  }
  values <- vapply(engine, function(one) {
    value <- capabilities(one, capability)
    if (isTRUE(value)) "true" else if (isFALSE(value)) "false" else "unknown"
  }, character(1))
  if (any(values == "true")) TRUE else if (any(values == "false")) FALSE else "unknown"
}

#' Does the run's thermo engine have a capability?
#'
#' The seam the design code asks through instead of probing the tool itself.
#' A [marker_run_context()] may carry an `engines` adapter ([marker_engines()]);
#' the call sites consult *that* table, so a caller who supplies a different
#' adapter gets that adapter's answer.  With no context (or no `engines` in it)
#' the declaration is taken from [marker_engines()], which is the live probe:
#' `duplex_dg` / `hairpin_dg` need OligoArrayAux unless
#' `options(tidymarkers.thermo_engine = "internal")` selects the pure-R adapter,
#' and `pcr_efficiency` needs it through DECIPHER.  Anything but a `TRUE`
#' in the table is a *no* - an undeclared or `"unknown"` capability is not
#' evidence that the engine can do it.
#'
#' @param capability a `thermo` capability name (`"tm"`, `"duplex_dg"`, ...).
#' @param engines optional [marker_run_context()] (its `engines` slot is read) or
#'   a `marker_engines` registry; `NULL` uses the live [marker_engines()].
#' @return a single logical.
#' @keywords internal
#' @noRd
engine_has <- function(capability, engines = NULL) {
  if (inherits(engines, "marker_run_context")) {
    engines <- engines$engines
  }
  if (is.null(engines)) {
    engines <- marker_engines()
  }
  isTRUE(capabilities(engines$thermo, capability))
}

#' OligoArrayAux duplex / hairpin free energy
#'
#' The DNA engine DECIPHER already uses: `hybrid-min` (duplex of two strands) or
#' `hybrid-ss-min` (intramolecular folding of one strand).  Returns dG in kcal/mol
#' (more negative = more stable); `NA` if OligoArrayAux is not on `PATH`, in which
#' case the missing-tool hint is printed once per session rather than once per
#' oligo.  A non-zero exit or an unparsable result is reported as a warning (not a
#' silent `NA`) on every call, and the run goes through `system2()` so failures
#' are explicit.
#'
#' @param seq1,seq2 oligos 5'->3'; with `seq2 = NULL` the hairpin (folding) dG of
#'   `seq1` alone is returned.
#' @param na sodium concentration (M).
#' @param temp temperature (C).
#' @return a single numeric dG in kcal/mol, or `NA_real_` when OligoArrayAux is
#'   unavailable or the run fails.  "No foldable structure" is `0` (no
#'   stabilisation relative to the unfolded state): UNAFold reports that case as
#'   a large positive sentinel (~1000 kcal/mol) from `-q`, which is mapped here
#'   rather than passed through.
#' @examples
#' \dontrun{
#' oligo_dg("GAAGGTGACCAAGTTCATGCT", "AGCATGAACTTGGTCACCTTC")
#' oligo_dg("GCGCGCGC") # hairpin (folding) dG
#' }
#' @family thermodynamics
#' @export
oligo_dg <- function(seq1, seq2 = NULL, na = 0.05, temp = 37) {
  bin <- if (is.null(seq2)) "hybrid-ss-min" else "hybrid-min"
  paths <- oligoarrayaux_paths()
  nm <- if (.Platform$OS.type == "windows") paste0(bin, ".exe") else bin
  exe <- if (is.null(paths)) Sys.which(bin) else fs::path(paths$bin, nm)
  if (!nzchar(exe) || !fs::file_exists(exe)) {
    # once per session: a cold machine asks this per oligo otherwise
    session_once(paste0("oaa_missing_", bin), {
      warn_msg("OligoArrayAux ('{bin}') was not found.")
    })
    return(NA_real_)
  }
  args <- c(
    "-n", "DNA", "-t", as.character(temp), "-T", as.character(temp),
    "-N", as.character(na), "-E", "-q"
  )
  args <- if (is.null(seq2)) c(args, seq1) else c(args, seq1, seq2)
  # -q makes UNAFold print a single number.  system2() is ~2x cheaper than a
  # processx object per oligo, and this is the hottest loop in design.
  # `env=` is not usable here: on Windows system2() appends "NAME=value" to the
  # command line instead of setting the variable (R's system2.Rd), so hybrid-min
  # reads the data path as a sequence and ignores the second oligo.  Set the
  # variable in the child's environment on every platform instead.
  dat <- if (is.null(paths)) NA_character_ else paths$dat
  # the resolver only knows its own prefix: an ambient UNAFOLDDAT still reaches
  # the child, so a missing share/ directory is not the same as missing tables
  if (!is.null(paths) && is.na(paths$dat) && !nzchar(Sys.getenv("UNAFOLDDAT"))) {
    # a run without the tables still returns a number, so name the gap once
    session_once("oaa_dat_missing", {
      warn_msg(
        "OligoArrayAux ('{bin}') has no UNAFold data tables:",
        "expected '{fs::path(fs::path_dir(paths$bin), 'share', 'oligoarrayaux')}'. dG may be meaningless."
      )
    })
  }
  out <- tryCatch(
    withr::with_envvar(
      oligoarrayaux_env(dat),
      suppressWarnings(system2(exe, args, stdout = TRUE, stderr = TRUE))
    ),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    warn_msg("OligoArrayAux ('{bin}') could not be run: {conditionMessage(out)}")
    return(NA_real_)
  }
  status <- attr(out, "status")
  if (!is.null(status) && as.integer(status) != 0L) {
    warn_msg("OligoArrayAux ('{bin}') exited with status {status}: {trimws(paste(out, collapse = ' '))}")
    return(NA_real_)
  }
  tok <- utils::tail(strsplit(trimws(paste(out, collapse = " ")), "[[:space:]]+")[[1L]], 1L)
  dg <- suppressWarnings(as.numeric(tok))
  if (length(dg) != 1L || !is.finite(dg)) {
    warn_msg("OligoArrayAux ('{bin}') returned no parsable dG (got '{tok}').")
    return(NA_real_)
  }
  # hybrid-ss-min / hybrid-min print a large POSITIVE value (998.701, 999.999,
  # 1002.39 were all observed with -q) when no fold or duplex can be found.
  # That is "no structure" (dG ~ 0 relative to the unfolded state), not a
  # binding free energy: passing it through put +999.999 kcal/mol into any mean
  # or sum over primers (found by the Primer3 cross-check, dev/primer3-crosscheck.R).
  if (dg > 100) {
    return(0)
  }
  dg
}


# --- small utilities ---------------------------------------------------------

#' Normalise chromosome names to "Chr1".."Chr5"/"ChrM"/"ChrC"
#' @keywords internal
#' @noRd
ath_chrom <- function(x) {
  key <- sub("^(chromosome|chr)", "", x, ignore.case = TRUE)
  out <- unname(ATH_CHROM_MAP[key])
  dplyr::if_else(is.na(out), paste0("Chr", key), out)
}

#' Read one genomic region as a DNAString
#' @keywords internal
#' @noRd
get_region <- function(genome, chrom, start, end) genome_seq(genome, chrom, start, end)


# --- primer geometry (single source of truth) --------------------------------
# All intervals are plus-strand, with `snp_pos` the SNP coordinate.  The SAME
# helpers are reused by design (in template coordinates: pass `snp_i`), by the
# assay assembly and by the plots, so the geometry cannot drift between them.

#' Plus-strand interval occupied by a primer
#'
#' `role = "AS"` = allele-specific primer (3'-terminal base at the SNP);
#' `role = "COM"` = the common primer, `gap` bases away on the other side.
#' @param snp_pos plus-strand SNP coordinate; the primer's 3'-terminal base is
#'   at `snp_pos` for `"AS"` and `gap + length` away for `"COM"`.
#' @param strand binding strand, `"+"` or `"-"`.
#' @param length primer length in nt.
#' @param gap distance from `snp_pos` to the near edge of a common primer.
#' @param role `"AS"` (allele-specific) or `"COM"` (common).
#' @return an [IRanges::IRanges] of width `length`, the primer's plus-strand
#'   interval.
#' @examples
#' primer_interval(100L, "+", 20L)
#' primer_interval(100L, "-", 20L, gap = 30L, role = "COM")
#' @family geometry
#' @export
primer_interval <- function(snp_pos, strand, length, gap = 0L, role = c("AS", "COM")) {
  role <- match.arg(role)
  if (identical(role, "AS")) {
    if (identical(strand, "+")) {
      IRanges::IRanges(snp_pos - length + 1L, snp_pos)
    } else {
      IRanges::IRanges(snp_pos, snp_pos + length - 1L)
    }
  } else {
    if (identical(strand, "+")) {
      IRanges::IRanges(snp_pos + gap, snp_pos + gap + length - 1L)
    } else {
      IRanges::IRanges(snp_pos - gap - length + 1L, snp_pos - gap)
    }
  }
}

#' Plus-strand interval of the AS primer's N 3'-terminal bases, SNP excluded
#'
#' `window` counts the 3'-terminal bases including the SNP itself, so the
#' returned interval has `window - 1` positions.
#' @param snp_pos coordinate of the variant's anchored base.
#' @param strand binding strand.
#' @param window number of 3'-terminal bases counted (includes the variant).
#' @return an [IRanges::IRanges] of the `window - 1` bases just inside the SNP.
#' @examples
#' as_3p_interval(100L, "+")
#' as_3p_interval(100L, "-")
#' @family geometry
#' @export
as_3p_interval <- function(snp_pos, strand, window = 5L) {
  if (identical(strand, "+")) {
    IRanges::IRanges(snp_pos - window + 1L, snp_pos - 1L)
  } else {
    IRanges::IRanges(snp_pos + 1L, snp_pos + window - 1L)
  }
}

#' Plus-strand coordinate of the AS primer's 3'-terminal base for a variant
#'
#' A variant is a substitution block `[v_start, v_start + ref_len - 1]` (a SNP is
#' the `ref_len = 1` case).  On the plus strand the AS primer's 3' end sits on
#' the block's last base; on the minus strand (where the primer is the reverse
#' complement) it sits on the block's first base.  All primer geometry is then
#' expressed through this single anchor, so SNP and MNP share one code path.
#'
#' `variant_offset` moves the block `0`, `1`, `2`, ... bases inside the 3' end:
#' `0` puts it on the terminal base (KASP), `1` on the penultimate base (ASQ /
#' PACE), `2` on the antepenultimate one.  The anchor is the primer's 3'-terminal
#' base, so only `variant_offset = 0` puts the anchor on the block: a larger
#' offset leaves the block behind and puts the anchor that many bases further
#' toward the primer's 3' terminus.  Measured on a 20-mer: `+` strand anchor
#' coordinate grows (`100, 101, 102`) while the block slides inwards
#' (`20, 19, 18`); `-` strand anchor shrinks (`100, 99, 98`) as the block moves
#' in (`1, 2, 3`).
#' @param v_start plus-strand start of the variant block.
#' @param ref_len length of the REF allele (block size; 1 for a SNP).
#' @param strand binding strand.
#' @param variant_offset bases between the variant block and the 3' terminus
#'   (`0` = terminal/KASP, `1` = penultimate/ASQ).
#' @return a single integer plus-strand coordinate.
#' @examples
#' as_anchor(100L, 1L, "+")
#' as_anchor(100L, 3L, "-")
#' as_anchor(100L, 1L, "+", variant_offset = 1L)
#' @family geometry
#' @export
as_anchor <- function(v_start, ref_len, strand, variant_offset = 0L) {
  if (identical(strand, "+")) {
    v_start + ref_len - 1L + variant_offset
  } else {
    v_start - variant_offset
  }
}

#' Index range of the variant block inside an AS-primer binding window
#'
#' For `strand = "+"` the block is at the 3' end of the window; for `"-"`
#' (reverse-complementary primer) it is at the 5' end.  `variant_offset` shifts
#' it that many bases inside the 3' terminus (see [as_anchor()]).
#' @keywords internal
#' @noRd
block_in_seg <- function(ref_len, primer_len, strand, variant_offset = 0L) {
  if (identical(strand, "+")) {
    IRanges::IRanges(primer_len - ref_len + 1L - variant_offset, primer_len - variant_offset)
  } else {
    IRanges::IRanges(1L + variant_offset, ref_len + variant_offset)
  }
}

#' Count variant positions inside an interval (IRanges)
#'
#' `exclude` drops positions (typically the target SNP) from the count.  One
#' *site* is one position however many ALT records it carries (a multi-allelic
#' VCF row is split into one row per allele upstream), so positions are
#' de-duplicated before counting - a site with two ALTs is not two variants.
#' @param interval an `IRanges` interval.
#' @param vpos integer vector of variant positions (or `NULL`).
#' @param exclude positions to drop (typically the target variant).
#' @return a single integer count (`NA_integer_` when `vpos` is `NULL`).
#' @examples
#' count_variants_in(IRanges::IRanges(10, 20), c(5L, 12L, 18L))
#' @family geometry
#' @export
count_variants_in <- function(interval, vpos, exclude = NA_integer_) {
  if (is.null(vpos)) {
    return(NA_integer_)
  }
  # `exclude` defaults to NA; `%in%` never yields NA for an NA pattern, so a
  # bare `!(vpos %in% exclude)` would keep every position and then build
  # zero-width ranges from a duplicated one.
  v <- unique(vpos[!is.na(vpos)])
  if (length(exclude) && !all(is.na(exclude))) {
    v <- v[!(v %in% exclude[!is.na(exclude)])]
  }
  if (!length(v)) {
    return(0L)
  }
  as.integer(length(IRanges::findOverlaps(IRanges::IRanges(v, v), interval)))
}

#' Weighted variant load of an interval (sum of per-variant weights)
#'
#' The count of [count_variants_in()] says *how many* known variants sit under a
#' primer; not all variants are equally bad, so when the variant table carries an
#' optional `pop_freq` (population allele frequency) we sum those instead of
#' counting.  Absent weights every variant counts 1, so the load degrades
#' gracefully to the plain count.  One *site* is one position: the ALT records
#' sharing a position (a split multi-allelic row) have their weights summed, so a
#' site's total non-reference frequency is counted once.
#' @keywords internal
#' @noRd
variant_load_in <- function(interval, vpos, vweight = NULL, exclude = NA_integer_) {
  if (is.null(vpos)) {
    return(NA_real_)
  }
  keep <- !is.na(vpos) & !(vpos %in% exclude)
  v <- vpos[keep]
  if (!length(v)) {
    return(0)
  }
  w <- if (is.null(vweight)) rep(1, length(vpos)) else as.numeric(vweight)
  w <- w[keep]
  # one site = one position: collapse the ALT rows sharing a position
  site <- tapply(w, v, sum)
  pos <- as.numeric(names(site))
  hits <- IRanges::findOverlaps(IRanges::IRanges(pos, pos), interval)
  sum(site[S4Vectors::queryHits(hits)])
}

#' Expected dropout under a primer's binding window (DECIPHER: hybridisation + Taq)
#'
#' For each known variant inside the window, build the variant-bearing template
#' and ask DECIPHER for the amplification efficiency; the expected dropout is the
#' frequency-weighted efficiency lost.  DECIPHER's `taqEfficiency` makes a
#' 3'-terminal mismatch far more costly than a mid-primer one, so position and
#' mismatch identity are handled by the library, and frequency by the weights.
#' Returns `0` when there is no (fit-in-window) variant.
#' @param primer primer sequence, 5'->3'.
#' @param plus_win the plus-strand sequence the primer covers (same length).
#' @param forward `TRUE` when the primer reads left-to-right on the plus strand.
#' @param var_offset 1-based offsets of the variants within `plus_win`.
#' @param var_alt the variant alleles (strings).
#' @param var_w per-variant weights (population frequency; `1` when unknown).
#' @keywords internal
#' @noRd
primer_dropout <- function(primer, plus_win, forward, var_offset, var_alt, var_w) {
  p <- dropout_pairs(primer, plus_win, forward, var_offset, var_alt, var_w)
  if (is.null(p)) {
    return(0)
  }
  eta <- pcr_efficiency(p$primer, p$target)
  min(1, sum(p$w * (1 - eta)))
}

#' Build the (primer, target, weight) units for a zone's expected dropout
#'
#' The primer is paired with each in-window variant allele in turn; the units are
#' kept so many zones' pairs can be pooled into one `pcr_efficiency()` call
#' (DECIPHER setup dominates when a zone holds only a few variants).  `NULL` when
#' no variant fits the window.
#' @keywords internal
#' @noRd
dropout_pairs <- function(primer, plus_win, forward, var_offset, var_alt, var_w) {
  if (!length(var_offset)) {
    return(NULL)
  }
  L <- nchar(plus_win)
  ok <- var_offset >= 1L & var_offset <= L & nchar(var_alt) >= 1L &
    (var_offset + nchar(var_alt) - 1L) <= L &
    !is.na(var_alt) & !is.na(var_w)
  if (!any(ok)) {
    return(NULL)
  }
  var_offset <- var_offset[ok]
  var_alt <- var_alt[ok]
  var_w <- var_w[ok]
  win <- Biostrings::DNAString(plus_win)
  alt_win <- Biostrings::DNAStringSet(vapply(
    seq_along(var_offset),
    function(i) {
      as.character(Biostrings::replaceAt(
        win, IRanges::IRanges(var_offset[i], width = nchar(var_alt[i])), var_alt[i]
      ))
    },
    character(1)
  ))
  target <- if (forward) Biostrings::reverseComplement(alt_win) else alt_win
  list(
    primer = rep(primer, length(target)),
    target = as.character(target),
    w = var_w
  )
}

#' Expected dropout of a window on the design template
#'
#' The one place that turns "a primer over the design template + the known
#' into `primer_dropout()`'s expected dropout, shared by the AS 3'
#' window and the COM footprint.  `0` when the window holds no (fit-in-window)
#' variant.
#' @keywords internal
#' @noRd
zone_dropout <- function(primer, start, end, forward, n_var, tmpl, off, vpos, vblock, valt, vweight) {
  if (is.na(n_var) || n_var == 0L) {
    return(0)
  }
  idx <- which(vpos >= start & vpos <= end & !(vpos %in% vblock))
  # one (position, ALT) unit each: a split multi-allelic row is not double-counted
  idx <- idx[!duplicated(paste(vpos[idx], valt[idx]))]
  if (!length(idx)) {
    return(0)
  }
  i1 <- max(start - off + 1L, 1L)
  i2 <- min(end - off + 1L, length(tmpl$seq))
  if (i2 < i1) {
    return(0)
  }
  w <- if (is.null(vweight)) rep(1, length(idx)) else vweight[idx]
  primer_dropout(
    primer, as.character(Biostrings::subseq(tmpl$seq, i1, i2)),
    forward, vpos[idx] - start + 1L, valt[idx], w
  )
}

#' Expected dropout of many zones, with every DECIPHER pair pooled into one call
#'
#' Numerically identical to calling `zone_dropout()` per zone, but the per-zone
#' (primer, target) units are concatenated and scored by a single
#' `pcr_efficiency()`.  With dozens of candidates each holding a few variants the
#' per-call DECIPHER setup dominates, so pooling is the difference between
#' seconds and minutes in a scan.  `forward` is a logical vector.
#' @keywords internal
#' @noRd
zone_dropouts <- function(primers, starts, ends, forwards, n_vars, tmpl, off, vpos,
                          vblock, valt, vweight, cores = 1L) {
  n <- length(primers)
  out <- numeric(n)
  if (n == 0L) {
    return(out)
  }
  weights <- vector("list", n)
  primer <- character(0)
  target <- character(0)
  span <- length(tmpl$seq)
  for (i in seq_len(n)) {
    if (is.na(n_vars[i]) || n_vars[i] == 0L) next
    idx <- which(vpos >= starts[i] & vpos <= ends[i] & !(vpos %in% vblock))
    idx <- idx[!duplicated(paste(vpos[idx], valt[idx]))]
    if (!length(idx)) next
    i1 <- max(starts[i] - off + 1L, 1L)
    i2 <- min(ends[i] - off + 1L, span)
    if (i2 < i1) next
    w <- if (is.null(vweight)) rep(1, length(idx)) else vweight[idx]
    p <- dropout_pairs(
      primers[i], as.character(Biostrings::subseq(tmpl$seq, i1, i2)),
      forwards[i], vpos[idx] - starts[i] + 1L, valt[idx], w
    )
    if (is.null(p)) next
    weights[[i]] <- p$w
    primer <- c(primer, p$primer)
    target <- c(target, p$target)
  }
  if (!length(primer)) {
    return(out)
  }
  eta <- pcr_efficiency(primer, target, cores = cores)
  k <- 1L
  for (i in seq_len(n)) {
    w <- weights[[i]]
    if (is.null(w)) next
    ei <- eta[k:(k + length(w) - 1L)]
    k <- k + length(w)
    out[i] <- min(1, sum(w * (1 - ei)))
  }
  out
}

#' Per-user cache directory for downloads and indices
#'
#' Everything heavy that the package fetches or builds (Ensembl FASTA/JSON, BLAST
#' and FM indices) lives here -- in the OS user cache, not in the project tree -- so
#' the working directory stays clean and the cache is shared across projects.
#' @return the cache directory path (created if needed).
#' @examples
#' marker_cache_dir()
#' @family input
#' @export
marker_cache_dir <- function() {
  d <- tools::R_user_dir("tidymarkers", which = "cache")
  fs::dir_create(d)
  d
}

#' Download (once) and load the TAIR10 genome as a DNAStringSet
#'
#' The FASTA comes from Ensembl Plants and is cached in the per-user cache
#' directory (`tools::R_user_dir("tidymarkers", "cache")`), not in the project.
#' Sequence names are normalised to Chr1..Chr5 / ChrM / ChrC (scaffolds
#' dropped), so they line up with the Ensembl gene records used elsewhere.
#' @param cache cache directory for the downloaded FASTA.
#' @param release Ensembl Plants release number.
#' @section Caching:
#' The download (~35 MB) happens once and is reused. The transfer is written to
#' a temporary file and moved into place only when it succeeds, so an
#' interrupted download cannot leave a truncated FASTA behind.
#' @return a `DNAStringSet` of the seven TAIR10 chromosomes (`Chr1`..`Chr5`,
#'   `ChrM`, `ChrC`).
#' @examples
#' \dontrun{
#' genome <- ath_genome()
#' length(genome)
#' }
#' @family input
#' @export
ath_genome <- function(cache = marker_cache_dir(), release = 63L) {
  fs::dir_create(cache)
  fa <- fs::path(cache, "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
  need <- c(paste0("Chr", 1:5), "ChrM", "ChrC")
  parse <- function() {
    g <- tryCatch(Biostrings::readDNAStringSet(fa), error = function(e) NULL)
    if (is.null(g)) {
      return(NULL)
    }
    nm <- sub("\\s.*$", "", names(g)) # keep only the first token of the defline
    dropped <- nm[!nm %in% names(ATH_CHROM_MAP)]
    if (length(dropped)) {
      cli::cli_alert_info(gettext(
        "Dropping {length(dropped)} non-chromosome sequence(s), e.g. {.val {utils::head(dropped, 3L)}}."
      ))
    }
    g <- g[nm %in% names(ATH_CHROM_MAP)]
    names(g) <- unname(ATH_CHROM_MAP[nm[nm %in% names(ATH_CHROM_MAP)]])
    g
  }
  if (!fs::file_exists(fa)) {
    tair10_download(fa, release)
  }
  g <- parse()
  if (is.null(g) || !all(need %in% names(g))) {
    # A truncated or renamed cache is otherwise permanent: file_exists() stays
    # TRUE forever, so refetch once before reporting the failure.
    if (fs::file_exists(fa)) {
      fs::file_delete(fa)
    }
    tair10_download(fa, release)
    g <- parse()
  }
  if (is.null(g)) {
    abort_msg(
      "ath_genome(): TAIR10 FASTA could not be read from {.path {fa}}.",
      i = "Delete the cached file and retry."
    )
  }
  missing <- setdiff(need, names(g))
  if (length(missing)) {
    abort_msg(
      "ath_genome(): TAIR10 FASTA is missing chromosome(s): {paste(missing, collapse = ', ')}.",
      i = "Ensembl deflines may have changed."
    )
  }
  g[need]
}

tair10_download <- function(fa, release) {
  url <- sprintf(paste0(
    "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-%d/fasta/",
    "arabidopsis_thaliana/dna/",
    "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
  ), release)
  tmp <- paste0(fa, ".part")
  on.exit(unlink(tmp), add = TRUE)
  cli::cli_alert_info(gettext("Downloading TAIR10 genome (~35 MB) from Ensembl Plants into {.path {dirname(fa)}} ..."))
  status <- utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
  if (!identical(status, 0L) || !fs::file_exists(tmp) || fs::file_size(tmp) == 0L) {
    abort_msg("ath_genome(): downloading TAIR10 from Ensembl Plants failed.")
  }
  fs::file_move(tmp, fa)
  invisible(fa)
}
