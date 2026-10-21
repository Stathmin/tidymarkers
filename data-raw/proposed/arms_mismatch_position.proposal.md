# Proposal: the deliberate-mismatch position contradiction (−2 vs −3)

Status: **open, and still open on purpose.** No position is averaged, weighted or
preferred anywhere in this package. This file records the state of the evidence
after the batch-2 search (`feat/sci-arms2`), with a count of independent
measurements for each position.

Ledger: `inst/extdata/arms_mismatch_position.csv`
Builder: `data-raw/arms_mismatch_position.R` (asserts every quote against its file)

## Numbering, stated once

All positions below are **negative from the 3′ end**: last base = −1.

| paper's wording | our label |
|---|---|
| "second mismatch", "2nd base from the 3′ end", "penultimate" | −2 |
| "3rd base from the 3′ end" | −3 |
| Lefever's scan position 3 (0 = terminal) | −3 |
| "4th nt", "fourth from the 3′ end" | −4 |

Most of the historical confusion is that papers number from 1 upward while we
number from −1 downward, so "position two" and "the second base from the 3′ end"
are only the same thing by convention.

## Counts after this batch

| position | independent sources | of which measurements | of which primer_scan | of which rule |
|---|---|---|---|---|
| −2 | 2 | 0 | 1 (Zhu 2021, sorghum) | 1 (Ye 2001) |
| −3 | 2 | **2** (Liu 2012; Lefever 2019) | 0 | 0 |
| −4 | 1 | 0 | 1 (Zhu 2021, same paper, different primer) | 0 |

"Measurement" = the paper reports an outcome at each position it scanned.
"primer_scan" = the position is certain because it is read off the designed
primer sequences, but the per-position outcome is only in a figure, so the
ranking is the authors' assertion rather than a number we can re-check.

## What changed in this batch — and it is not what we expected

**Liu 2012 no longer supports −2. It supports −3.**

The package filed Liu 2012 as the largest "−2 equivalent" anchor. Its own
sentence ranks the positions:

> "polymorphism percents of the 2nd, 3rd and 4th base location closest to the
> SNP sites were 25.1% (139/554), 32.9% (189/575) and 29.1% (162/557)"

Recomputed: 139/554 = 25.1%, 189/575 = 32.9%, 162/557 = 29.1%, total 490/1686 =
29.1%, which is exactly the 1686 primers and 490 successes Liu reports. The
**3rd base from the 3′ end (−3) is the peak**, above both −2 and −4. The builder
asserts `which.max(pct) == 2L` so the filing cannot silently drift back.

Consequence: once Liu is counted where its data points, **two independent
controlled scans (Liu 1686 primers; Lefever 5 positions × 6 SNPs) both peak at
−3**, and neither controlled scan supports −2. −2's only quantitative-ish
support is Zhu 2021's sorghum primer set, where the −2 primer is the one called
specific but the per-primer Cq values exist only as Figure S2 (an image).

## What is deliberately NOT done

* **No averaging.** The counts above are the output; no preferred position is
  written to `z7_priors.csv` or anywhere else.
* **No silent resolution.** The builder `stop()`s if either side reaches zero
  sources, because that would be a scientific conclusion for a human to make
  and then record, not a side effect of a build.
* **No upgrade of Zhu 2021.** Its positions are solid (read from Table S1
  sequences), but its *ranking evidence* is figure-only, so it stays
  `primer_scan`, not `measurement`.

## One live conflict inside a single paper

Zhu 2021's prose says the mismatch "placed fourth from the 3′ end of wxt F4"
gave the best discrimination. Aligning the primers shows **wxt F4 carries two
mismatches, at −1 and −4**, so the supposed −4 test is really a −1 + −4
construct and cannot isolate −4. The paper's other primer, wxc F4, carries a
clean single mismatch at −2. Both are landed as separate rows. A paper that
reports "position −4 is best" on a primer that also carries −1 is not evidence
about −4 alone.

## Recommendation to the maintainer (a decision, not an edit)

1. Adopt **−3 as the default** in the deliberate-mismatch design rule, citing Liu
   2012 and Lefever 2019 as the two measurements that agree, and note Liu as a
   correction of the earlier reading.
2. Keep Ye 2001's −2 as a **historical rule of thumb**, not a competing
   measurement — it has no outcome behind it in this package.
3. Do not treat Zhu 2021 as support for −2 until its Figure S2 Cq values are
   extracted (they are plotted, not tabulated; the CSV marks them
   `outcome_available = FALSE`).

None of this has been applied. It requires editing `z7_priors.csv`, which is
outside this worker's allow-list.
