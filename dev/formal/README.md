# Formal model of the lexicographic candidate ordering

`lexicographic_order.py` restates the ordering that `R/marker_select.R` uses to
rank marker candidates, and states the algebra that ordering is supposed to
satisfy. A symbolic checker ([CrossHair](https://github.com/pschanely/CrossHair))
then searches for a counterexample instead of a human reading cases. The
checker proves the **model**; the mapping from the R code to the model is
argued below, not proven.

The artifact exists because `ROADMAP.md`'s `ZY13` asked for one algebraic
invariant proven against a real prover. It is referenced by `VALIDATION.md`.

## The R ordering, and the model of it

`order_by_margins()` (`R/marker_select.R:405`) is exactly

```r
dplyr::arrange(
  data,
  dplyr::across(dplyr::all_of(intersect(desc, names(data))), dplyr::desc),
  dplyr::across(dplyr::all_of(intersect(asc, names(data))))
)
```

so every `desc` column sorts high-to-low, then every `asc` column sorts
low-to-high, and a row's position is decided by the first column on which it
differs. `rank_candidates()` (`:419`) and `select_with_ladder()` (`:449`) both
go through it, including the no-winner fallback at `:484` and the gated
survivors at `:490` (the branch that asks `bncc_rank()` to rank them instead).

The model turns each sorted column into a **rank key** `d * v` with
`d = +1` for a `desc` column and `d = -1` for an `asc` one. Larger key is
always better, so the ordering becomes the lexicographic comparison over the
key tuple, with the first differing coordinate deciding. `lex_compare()`
(`lexicographic_order.py:50`) is that comparator; `lex_compare_by_refinement()`
(`:70`) states the same ordering the other way — first differing coordinate
decides, all-equal is a tie — so the checks can cross the two formulations
against each other.

## Verdicts

Six checks, run with `crosshair-tool 0.0.114` on Python 3.13.5. Each states
its precondition first (list lengths bounded, `direction` entries exactly `-1`
or `+1`) and its property after an assignment.

| check                                      | line | max iterations | verdict                  | wall    |
| ------------------------------------------ | ---- | -------------- | ------------------------ | ------- |
| `tie_is_reflexive`                         | 97   | 2000           | Confirmed over all paths | 1.2 s   |
| `ranking_is_asymmetric`                    | 106  | 2000           | Confirmed over all paths | 10.4 s  |
| `ranking_is_transitive`                    | 118  | 20000          | Confirmed over all paths | 233.4 s |
| `ranking_is_total`                         | 138  | 2000           | Confirmed over all paths | 10.0 s  |
| `comparator_matches_successive_refinement` | 151  | 2000           | Confirmed over all paths | 7.5 s   |
| `directions_only_rescale`                  | 163  | 2000           | Confirmed over all paths | 16.1 s  |

Reproduce one of them:

```sh
python3 -m venv .tmp/chvenv
.tmp/chvenv/bin/pip install crosshair-tool==0.0.114
.tmp/chvenv/bin/crosshair check --analysis_kind asserts --report_all \
  --per_path_timeout 60 --per_condition_timeout 900 \
  dev/formal/lexicographic_order.py:106
```

`--version` is not a CrossHair flag; get the version from `pip show
crosshair-tool`. Exit code 0 is not the signal — read the `info:` line and
accept only `Confirmed over all paths`; `Not confirmed` is a real failure.

The model also carries a dependency-free `_selftest()` (`:185`) with three
concrete pairs, so a reader without the prover can still run it:

```sh
python3 dev/formal/lexicographic_order.py    # -> selftest ok
```

`dev/check.R` runs exactly that, and skips when `python3` is absent. The six
CrossHair checks are not a gate step: they need a virtualenv CI does not have
and take about five minutes, so they are run by hand when the model changes.

## Bounded domain

The calltree search does not terminate on unbounded list lengths — it aborts
with `Aborted calltree search with UNKNOWN` (measured: 2138 iterations with
the cap at 2000). The preconditions therefore bound the lengths:

```python
MAX_KEY_LEN = 3
MAX_DIRECTION_LEN = 3
```

That makes the search finite: every key tuple of at most three integers and
every direction list of at most three entries of `-1` or `+1`. This is a
**bounded-domain proof**, not a proof for all lengths. Raising
`--max_uninteresting_iterations` was what moved `ranking_is_total` from
UNKNOWN to CONFIRMED at 10 s; transitivity needed 20000 and about four
minutes. A change that raises `MAX_KEY_LEN` must re-run the checks.

## What is argued, not proven

The checks prove the algebra of the Python comparator. Three steps take that
result to the R code, and none of them is a proof:

1. **`NA`/`NaN` never reach the comparator.** `dplyr::arrange()` sorts missing
   ranks to the end before comparing, in both directions. Measured on R 4.6.1:
   `tibble(k = c(2, NA, NaN, 1))` sorted with `desc(k)` gives `2 1 NA NaN`, and
   with `k` gives `1 2 NA NaN`. The comparator therefore only ever sees the
   ordinary values the model reasons over.
2. **The key type.** Coordinates are modelled as `int`, not `float`. The
   ordering uses only trichotomy of `>`, and CrossHair reasons exactly over
   unbounded integers, whereas its float support explores only special values
   (measured: on `assert x + 0.1 != x` it finds `x = inf` and nothing else). A
   `NaN` rank key has no trichotomy at all — it is outside this model's domain,
   and point 1 is why it cannot occur.
3. **Unequal lengths.** The comparator orders a strict prefix before the longer
   sequence, like Python tuples. That branch exists only because the model takes
   lists; at the R call site every key is a data-frame column, so all keys of a
   comparison have equal length and the branch is outside the package's domain.
   A column that is neither in `desc` nor in `asc` is not sorted on at all,
   which is why a short `direction` reads as `+1` for the remaining coordinates
   (`_key_at()`, `:36`) and why only the sorted prefix is modelled.

## The checks are not vacuous

A mutation shows the checker finds a real counterexample. Flip the tie test in
the comparator and check a property:

```sh
mkdir -p .tmp/chprobe
sed 's/if ka > kb:/if ka >= kb:/' dev/formal/lexicographic_order.py > .tmp/chprobe/mut.py
.tmp/chvenv/bin/crosshair check --analysis_kind asserts --report_all \
  --per_path_timeout 60 --per_condition_timeout 900 .tmp/chprobe/mut.py:106
```

Result: exit 1, `mut.py:114: error: AssertionError:  when calling
ranking_is_asymmetric([0], [0, 0], [])`. The mutant breaks asymmetry and the
checker says so.

One trap worth knowing: CrossHair's `asserts` analysis kind reads a function's
**leading** asserts as preconditions. A function whose whole body is asserts
therefore proves nothing — the run passes vacuously. Every check here makes an
assignment before the asserted property, and the precondition asserts are the
ones it is allowed to assume.

## Status

- Model: `lexicographic_order.py`, six checks confirmed over their bounded
  domain, mutation check re-run.
- Gate: `dev/check.R` runs the dependency-free `_selftest()` only.
- Not covered by the model: `bncc_rank()`'s Pareto step and the gate cascade.
  `ZY13` proved one invariant, as scoped; the rest of the ranking is covered by
  examples and the mutation tests in `dev/mutation-test.R`.
