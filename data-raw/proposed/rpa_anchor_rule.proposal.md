# Proposal — a 3'-anchor mismatch weight scheme for RPA primer design

> **DECIDED 2026-10-21** — see `rpa_anchor_rule.decision.md`: implemented as `RPA_ANCHOR_NT` / `rpa_anchor_weights()` / `rpa_anchor_score()` in `R/marker_rpa.R`. The status
> line below is the state this file was written in; it is kept as the evidence
> record and is not edited after the decision.

Status: **proposal only.** No file under `R/` was touched. Everything below is
re-derivable with `Rscript data-raw/rpa_mismatch_class.R`, which re-reads
`.tmp/lit/rpa_mismatch_class/supplements/mmc1.xlsx` and `fulltext.xml` and
`stop()`s if any cell, count or arithmetic assertion no longer holds.

Source: Higgins M, et al. "Characterizing the Impact of Primer-Template
Mismatches on Recombinase Polymerase Amplification." J Mol Diagn 2022.
doi:10.1016/j.jmoldx.2022.08.005, PMC12179515 (CC-BY). Supplement S1 coef-
ficients are the calibration set; the paper's own numbers are asserted, not
assumed.

## 1. What the quantity is (and is not)

The landed number is `RS_logit`: the **Firth logistic-regression coefficient**
on a BINARY reaction-success call (fluorescence threshold reached or not),
relative to a fully complementary primer binding site, with log(copy number)
in the model. It is *not* a fold change, a Cq/Ct shift, or a probability — the
probability in `inst/extdata/rpa_mismatch_terminal.csv` is
`plogis(RS_INTERCEPT + RS_logit)` and is checked against the two probabilities
the paper prints verbatim (C–C 0.59, G–A 0.62).

The field papers (`rpa_asrpa_mdr_calls.csv`, `rpa_asrpa_sickle_calls.csv`)
report **call agreement against a reference method**, a different quantity. They
motivate the scheme and bound its usefulness; they do not calibrate a weight.

## 2. The anchor

The paper defines the anchor as "the pre-ante-penultimate to 3'-terminal
position" = **4 nt** (offset 0 = the primer's 3' terminal base). That width is
asserted three ways: the quoted sentence, the 9-character `????-????` pattern
width of all 315 S1 rows, and the offsets used by the 36 position classes.

## 3. Weight scheme

`w(position, mismatch)` is the penalty added to a primer's specificity score,
in the same log-odds units as `RS_logit`. Negative = worse.

| Position (offset) | Mismatch | Weight (RS_logit) | Table cell |
|---|---|---|---|
| 3 (pre-ante-penultimate) | any | −2.25 | S2 row `3n` (S1's `??C?-??C?` is −2.25, unstarred) |
| 2 (ante-penultimate) | any | −1.84 | S2 row `2n` (S1's `??T?-??C?` is −2.25, unstarred) |
| 1 (penultimate) | any | −1.71 | S2 row `1n` (the only non-terminal class significant on RS: **) |
| 0 (terminal) | A–C, T–G | 0 | S1 rows `???A-???C` / `???T-???G`, column H, unstarred |
| 0 | C–A, G–G, T–C | −2.32 | S1 rows `???C-???A` / `???G-???G` / `???T-???C`, column H |
| 0 | T–T | −2.54 | S1 row `???T-???T` |
| 0 | G–T | −3.19 | S1 row `???G-???T` |
| 0 | C–T | −3.40 | S1 row `???C-???T` |
| 0 | A–A | −3.55 | S1 row `???A-???A` |
| 0 | A–G | −3.84 | S1 row `???A-???G` |
| 0 | G–A | −4.31 | S1 row `???G-???A` |
| 0 | C–C | −4.43 | S1 row `???C-???C` (CI −5.51 to −3.44) |

Two design choices need justifying, because both are inferences from the table
rather than sentences in it:

1. **The four non-terminal positions are treated as one class per offset.**
   S1's 24 non-terminal single-mismatch rows are nucleotide-specific and mostly
   NOT significant (only 1n–G–A, 1n–G–G, 1n–C–T and 3n–T–T carry any star on
   RS), so the table does not support a per-base weight away from the terminus.
   The values in rows 1–3 above are the class coefficients from S2 (`3n`
   −2.27, `2n` −1.82, `1n` −1.71), read out of a different sheet than the
   per-base rows.
   **Open question:** S1 and S2 print slightly different values for the same
   class because they come from separate model fits; a future regression should
   pick one sheet and say which.
2. **Terminal weights are used as printed, including unstarred ones.** A–C and
   T–G are not distinguishable from zero (P ≥ 0.05 as printed), which is a
   finding — RPA is tolerant of exactly these two — not missing data. Encoding
   them as 0 rather than as "unknown" is what the paper's own discussion of
   PCR-like tolerance implies, but it is the paper's tolerance class, not a
   measured zero.

Beyond the 4-nt anchor, no mismatch is penalised: the paper's model has no
coefficients there.

## 4. Pair terms

Adjacent mismatches are much worse than either alone, and the table contains
those classes (`rpa_mismatch_position.csv`, 279 rows over 36 classes):

- terminal + penultimate (`??AA-??GG` family): −4.62 for T1n–A–G, but the
  worst cells are `??CA-??CC` / `??CC-??CA` / `??TC-??CC` at −6.77, which the
  paper calls complete inhibition
- terminal + ante-penultimate: −6.77 (`?G?A-?A?C`; the mildest 2n cell is
  `??AC-??CA` at −1.80, so this face of the table is the most spread out)
- terminal + pre-ante-penultimate: −6.60 (`C??A-A??A`; `??AA-??CA` is −3.55)

Proposed rule: when two adjacent mismatches fall inside the anchor, use the
**class coefficient** from `rpa_mismatch_position.csv` instead of summing two
single-mismatch weights. Summing terminal (−4.43) and penultimate (−1.71) for a
C–C/T1n pair gives −6.14 against a measured −6.77; summing is closer to the
truth than either term alone, but the class cell is the measured value.

## 5. What the scheme must NOT be used for

- **Predicting a call rate.** The AS-RPA field data show a designed assay
  missing 8 of 54 rifampicin-resistant isolates (14.8%) and 26 of 64
  isoniazid-resistant isolates. A −4.43 terminal penalty does not translate
  into an error rate; only a designed assay on real samples measures that.
- **Ranking alleles at a transition.** Sickle-cell AS-RPA reaches 100%
  sensitivity for both alleles with a design the paper attributes to the
  penultimate position, which this table scores at −1.71. One genotype cannot
  calibrate a weight, and the two papers are not measuring the same thing.
- **Anything outside the anchor**, or with more than two mismatches.

## 6. Validation this proposal still needs

1. Agree on S1 vs S2 for the non-terminal classes (see §3.1).
2. Check the weights against the 279 position-class cells: a fit of
   `observed ~ w(pos1) + w(pos2)` should reproduce the printed class
   coefficients within their CIs. The cells are already extracted; nothing new
   needs fetching.
3. Decide whether the unstarred terminal cells (A–C, T–G) are 0 or NA.
