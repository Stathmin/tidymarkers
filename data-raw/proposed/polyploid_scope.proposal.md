# Polyploid scope: which of the package's polyploid claims are wheat-specific

Every sentence below cites the row it is re-checked against: `polyploid_*` keys
are rows in the four landed `inst/extdata/polyploid_*_*.csv` tables, `makhoul_*`
/ `grewal_*` / `conversion_rates` are the existing wheat tables. Nothing in `R/`
was changed by this pass; this is a list of claims to re-check, not a patch.

## 1. The subgenome axis is carried by the COMMON primer's 3' base

This rule comes from PolyMarker's source and Grewal's two-axis proof
(`polymarker_design_rules.csv` role `common-primer`;
`grewal_kasp_axes.csv`: the AS 3' base tracks the allele, never the chromosome).
It is a statement about *KASP construct design*.

- Wheat-specific in the sense that matters: no non-wheat source in this family
  reports a subgenome-specific KASP construct at all. `cotton2023` screened 11
  allotetraploid assays on cluster separation and records
  `homeolog_handled = "not stated"` for every marker; `gao2022` records the same
  for sugarcane.
- The one non-wheat row that does touch the axis weakens rather than supports
  the rule's necessity: `cotton2023` rows for `i60979Gt` and `i64854Gm` are
  `usable = TRUE` assays designed from *diploid* relatives, used unre-designed on
  an AADD panel. So the package should not state that an allotetraploid KASP
  without a subgenome step is expected to fail.
- `durum2019` row 2 (`homeolog_handled = "yes"`) shows the *SNP-selection* axis
  is separable from the *primer* axis: homoeolog-resolved SNP discovery (4,461
  hemi + 3,479 simple of 7,940 chip SNPs) with an ordinary array assay.

**Verdict:** keep the two-axis rule as a statement about wheat KASP constructs;
add that no non-wheat assay in the family demonstrates it, and that one
allotetraploid (cotton) succeeded without it.

## 2. `specificity_policy("homeolog")` — the paralog tolerance row

`R/marker_specificity.R` maps `homeolog` to the `paralog` policy
(max_mismatch = 4, seed_len = 8, seed_max = 1, three_prime_max = 1) and
documents it with wheat homoeologs. The rows that would re-calibrate it do not
exist in this family: `jaimes2024` is autopolyploid (`homeolog_handled = "not
applicable"`), and `cotton2023` / `gao2022` never separate homoeologs
(`"not stated"`). Nothing here can test the tolerance numbers.

**Verdict:** keep, but it must stay `calibrated = FALSE` with an explicit note
that the *only* supporting system is wheat (`makhoul_kasp_primers.csv`,
`grewal_kasp_axes.csv`), because applying a 4-mismatch homeolog tolerance to an
AADD cotton design is unvalidated.

## 3. "Homeolog co-amplification is the hard case for polyploids"

The package treats polyploidy as one problem. The family says there are two
axes, and the second one has a different remedy:

- Dosage: `jaimes2024` rows — 42 of 76 sugarcane SNPs got a KASP dosage; the 34
  failures are cluster-fit failures ("did not cluster within the genotypic
  assignment of a decaploid organism"), i.e. the failure is the *model*, not the
  primer. The `outcome_definition` column says in so many words that this is not
  an amplification outcome.
- Dosage, nested: `gao2022` rows 1–3 — 77 → 49 → 24, where the last step is a
  map-consistency filter, again not a primer-specificity step.

**Verdict:** the package needs a sentence that separates "no subgenome exists
therefore no homoeolog axis" (autopolyploid) from "subgenomes exist but were not
handled" (cotton, sugarcane hybrids). The `homeolog_handled` vocabulary already
carries `not applicable` vs `not stated`, and `test-polyploid-scope.R` pins the
distinction.

## 4. The naming: "homeolog" vs "homoeolog"

The package's public vocabulary is `homeolog` (`KASP_SPECIFICITY_MODES =
c("unique", "homeolog-aware", ...)`), while the non-wheat sources and BMC's
wheat sources spell it `homoeolog` (`durum2019` rows, Jaimes' text). This is a
search-and-documentation hazard for anyone re-checking these rows: a grep for
`homoeolog` misses the package's own rule text.

**Verdict:** state the alias once in the user-facing doc; do not rename.

## 5. Ploidy level changes the outcome metric, not just the rate

- `durum2019` rows 1–3: hexaploid-derived vs durum-derived assays on the same
  288-accession AABB panel — failed 6.5% vs 8.8%, polymorphic 30.2% vs 51.7% of
  non-failed. Adding the D genome does **not** raise no-calls; it lowers
  informativeness.
- `cotton2023`: 8 of 11 usable by cluster separation.
- `jaimes2024`: 42 of 76 by dosage assignment.

**Verdict:** the package must never present a single "polyploid conversion
rate". Three sources here measure three different things on three ploidies, and
`test-polyploid-scope.R` asserts no file carries a pooled-rate column; the
existing Z8 rule (no two studies pooled into one rate) extends to
`conversion_rates.csv`, where all four new keys must stay out until the
denominators are made commensurable.
