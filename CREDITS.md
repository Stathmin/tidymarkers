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
| fami | Proposed the acquisition and knowledge-graph change taken in `b1124a0b`: a re-keyed DOI stays one work, and the alias surface is scoped | proposal → `b1124a0b` (Merovingen) | 2026-10-26 | The proposal, not the commit, is theirs. The fix it prompted is the `_dois_disagree()` suffix rule in `src/acquisition/dedup.py`. |

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
