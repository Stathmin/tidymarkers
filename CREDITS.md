# Credits

Who moved the package, and what it cost them. This is the record the
`CONTRIBUTING.md` reviewer protocol writes into: a contribution is credited for
what it makes checkable, not for how many lines it touched.

Maintenance: add a row when a contribution lands; date it `YYYY-MM-DD`; link the
commit or pull request; keep the note to one line. A row never claims more than
the repository can show — if the evidence is a test, name the test.

## Contributors

| Person | Contribution | Lands as | Date | Note |
|---|---|---|---|---|
| Daniil "Stathmin" Ulyanov | Package author and maintainer: the design rules, the gate, the falsifiability test, the decision ladders | 311 commits | 2026-10-23 | Owns the whole history to date; the row exists so the gap is visible. |
| Alexey "alermol" Ermolaev | Reported two real defects from use: the fixed-trio path graded at a different Tm salt method than the `KASP_LIMITS` windows, and `marker_blast_index()` documented a database path it did not return while rebuilding on every call | bug report → `R/marker_assemble.R`, `R/marker_scan.R`, `tests/testthat/test-scan.R` | 2026-11-10 | The report is theirs; the fixes are not. The Tm gap measured 3.51 °C; the scan idempotence test (`<fasta>.nsq` exists → skip) is the check. |

## How a contribution is ranked here

The leaderboard below ranks by what a contribution unblocks, not by novelty.
That is the honest ordering for a package still filling its evidence base.

| Rank | Kind of work | Why it ranks there |
|---|---|---|
| 1 | A rule plus the measurement behind it | Changes what the package tells a user, and can be re-derived |
| 2 | A fixture that lets an existing criterion fail | Turns an untested claim into a test |
| 3 | A wiring fix — code that existed but nothing called | Restores value already paid for |
| 4 | A gate or tooling fix | Protects every future change |
| 5 | Prose that matches the code | Cheap, and the cheapest thing to get wrong |

## Not credited

- Bulk edits with no verifiable intent.
- A source quoted only by DOI or title: a citation is provenance, not proof.
- Anything whose acceptance criterion cannot be driven to `fail`.
