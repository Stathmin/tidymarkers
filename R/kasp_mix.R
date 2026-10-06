# =============================================================================
#  kasp_mix.R  --  reaction recipe
# =============================================================================
# =============================================================================
#  8. Reaction recipe (output / ordering)
# =============================================================================

#' Tidy KASP primer-mix + reaction recipe
#'
#' Standard KBioscience/LGC KASP chemistry: the two allele-specific primers are
#' used at half the concentration of the common primer, AS1:AS2:COM = 2:2:5 (the
#' usual working mix is 12/12/30 uM).  This returns the recipe as tidy data:
#'
#' * `stage = "primer mix"` - how to build the working primer mix from the
#'   orderable 100 uM oligo stocks;
#' * `stage = "reaction"` - the per-reaction setup (2x master mix, assay mix,
#'   DNA, water), scaled to `n_reactions` with a small `overfill` for pipetting.
#'
#' It is pure arithmetic, so it carries no chemistry assumptions beyond the
#' published LGC recipe; change `as_uM`/`com_uM` for other chemistries.
#'
#' @param n_reactions number of wells to prepare for.
#' @param mix_uL total volume of working primer mix to prepare.
#' @param stock_uM oligo stock concentration as ordered (IDT default 100 uM).
#' @param as_uM,com_uM working concentrations of each AS and the common primer.
#' @param reaction_uL final reaction volume per well.
#' @param master_mix_conc fold-concentrate of the supplied master mix.
#' @param dna_uL template DNA volume per reaction.
#' @param overfill fraction of extra volume to compensate pipetting loss.
#' @return tibble `stage`, `component`, `stock_uM`, `final_uM`, `per_rxn_uL`,
#'   `total_uL`.
#' @examples
#' kasp_mix_recipe(n_reactions = 96)
#' @family KASP output
#' @export
kasp_mix_recipe <- function(n_reactions = 96L, mix_uL = 100,
                            stock_uM = 100, as_uM = 12, com_uM = 30,
                            reaction_uL = 10, master_mix_conc = 2,
                            dna_uL = 1, overfill = 0.10) {
  ratios <- c(AS1_FAM = as_uM, AS2_HEX = as_uM, COM = com_uM)
  stock_share <- ratios / stock_uM # fraction of the mix
  water_share <- 1 - sum(stock_share)
  if (water_share < 0) {
    cli::cli_abort(c(
      glue::glue(gettext(
        "Requested working mix ({as_uM}/{as_uM}/{com_uM} uM) exceeds the {stock_uM} uM stock."
      )),
      i = gettext("Raise `stock_uM` or lower the targets.")
    ))
  }

  primer_mix <- tibble::tibble(
    component = c(names(ratios), "nuclease-free water"),
    stock_uM = c(rep(stock_uM, length(ratios)), NA_real_),
    final_uM = c(unname(ratios), NA_real_),
    per_rxn_uL = NA_real_,
    total_uL = c(unname(stock_share) * mix_uL, water_share * mix_uL),
    stage = "primer mix"
  )

  assay_per_rxn <- round(0.014 * reaction_uL, 3)
  master_per_rxn <- reaction_uL / master_mix_conc
  water_per_rxn <- reaction_uL - master_per_rxn - assay_per_rxn - dna_uL
  if (water_per_rxn < 0) {
    cli::cli_abort(glue::glue(gettext(
      "Reaction volumes exceed {reaction_uL} uL; check `dna_uL`/`reaction_uL`."
    )))
  }
  total_rxn <- ceiling(n_reactions * (1 + overfill))

  reaction <- tibble::tibble(
    component = c(
      glue::glue("{master_mix_conc}x KASP master mix"), "primer mix",
      "template DNA", "nuclease-free water"
    ),
    stock_uM = NA_real_,
    final_uM = NA_real_,
    per_rxn_uL = c(master_per_rxn, assay_per_rxn, dna_uL, water_per_rxn),
    total_uL = c(master_per_rxn, assay_per_rxn, dna_uL, water_per_rxn) * total_rxn,
    stage = "reaction"
  )

  dplyr::bind_rows(primer_mix, reaction)
}
