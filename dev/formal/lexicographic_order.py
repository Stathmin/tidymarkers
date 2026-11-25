"""CrossHair model of the lexicographic candidate ordering.

`R/marker_select.R` orders candidates with `order_by_margins()`:
`dplyr::arrange(desc(tier_1), desc(tier_2), ..., asc(tie_1), ...)`.  This file
re-expresses that ordering as a comparator and states the algebra it is
supposed to satisfy, so a symbolic-execution checker can look for a
counterexample instead of a human looking at cases.

Two things are proven here and one is not.  Proven: the comparator's algebra
(reflexivity of ties, asymmetry, the two closure conditions of a total
preorder, totality, and agreement with the "first differing coordinate
decides" statement of the same ordering).  Not proven: that
`order_by_margins()` *is* this comparator -- that mapping is argued in
`dev/formal/README.md`.

Coordinates are modelled as `int`, not `float`.  The ordering only uses
trichotomy of `>` on the key type, and CrossHair reasons exactly over
unbounded integers, whereas its float support only explores a handful of
special values (measured: on `assert x + 0.1 != x` it found `x = inf` and
nothing else).  A float rank key that is NaN has no trichotomy at all, so it
is outside this model's domain; `dplyr::arrange()` sorts `NA`/`NaN` to the end
before any comparison happens, so the comparator never sees one.

CrossHair's `asserts` analysis kind reads the leading `assert` statements of a
function as its preconditions and checks the ones that come after the first
non-assert statement.  Every check below therefore leads with the real
precondition of the R call site -- `direction` entries are exactly `-1` (low
is better) or `+1` (high is better) -- and states the property after an
assignment, never as the first statement (a leading assert would be assumed,
not checked, and the run would pass vacuously).
"""

from typing import List


def _key_at(values: List[int], direction: List[int], i: int) -> int:
    """The rank key at coordinate `i`: larger is always better there.

    `direction[i]` is +1 for a "high is better" coordinate (a tier margin) and
    -1 for a "low is better" one (a tie-break).  A direction list shorter than
    the value list is read as +1 for the remaining coordinates, which is how
    `order_by_margins()` behaves: a column that is neither in `desc` nor in
    `asc` is not sorted on at all, and this model only covers the sorted
    prefix.
    """
    d = direction[i] if i < len(direction) else 1
    return d * values[i]


def lex_compare(a: List[int], b: List[int], direction: List[int]) -> int:
    """+1 if `a` ranks before `b`, -1 if after, 0 if the keys are equal.

    A pair of unequal lengths is ordered like Python tuples: the shared prefix
    decides first, and a strict prefix ranks before the longer sequence.
    """
    for i in range(max(len(a), len(b))):
        if i >= len(b):
            return 1
        if i >= len(a):
            return -1
        ka = _key_at(a, direction, i)
        kb = _key_at(b, direction, i)
        if ka > kb:
            return 1
        if ka < kb:
            return -1
    return 0


def lex_compare_by_refinement(a: List[int], b: List[int],
                              direction: List[int]) -> int:
    """The same ordering stated the other way: the first differing coordinate
    decides, and only an all-equal pair is a tie."""
    for i in range(max(len(a), len(b))):
        if i >= len(b):
            return 1
        if i >= len(a):
            return -1
        if a[i] != b[i]:
            return 1 if _key_at(a, direction, i) > _key_at(b, direction, i) else -1
    return 0


# --- the invariants ---------------------------------------------------------
#
# CrossHair explores list lengths without bound, so the calltree search never
# terminates on `range(max(len(a), len(b)))` and reports UNKNOWN.  The
# preconditions below bound the lengths, which makes the search finite: every
# key tuple of at most MAX_KEY_LEN entries and every direction list of at most
# MAX_DIRECTION_LEN entries, over all integers.  That is a bounded-domain
# proof, and `dev/formal/README.md` says so.

MAX_KEY_LEN = 3
MAX_DIRECTION_LEN = 3


def tie_is_reflexive(a: List[int], direction: List[int]) -> bool:
    assert len(a) <= MAX_KEY_LEN
    assert len(direction) <= MAX_DIRECTION_LEN
    assert all(d in (-1, 1) for d in direction)
    verdict = lex_compare(a, a, direction)
    assert verdict == 0
    return True


def ranking_is_asymmetric(a: List[int], b: List[int],
                          direction: List[int]) -> bool:
    assert len(a) <= MAX_KEY_LEN
    assert len(b) <= MAX_KEY_LEN
    assert len(direction) <= MAX_DIRECTION_LEN
    assert all(d in (-1, 1) for d in direction)
    ab = lex_compare(a, b, direction)
    ba = lex_compare(b, a, direction)
    assert not (ab == 1 and ba == 1)
    return True


def ranking_is_transitive(a: List[int], b: List[int], c: List[int],
                          direction: List[int]) -> bool:
    assert len(a) <= MAX_KEY_LEN
    assert len(b) <= MAX_KEY_LEN
    assert len(c) <= MAX_KEY_LEN
    assert len(direction) <= MAX_DIRECTION_LEN
    assert all(d in (-1, 1) for d in direction)
    ab = lex_compare(a, b, direction)
    bc = lex_compare(b, c, direction)
    ac = lex_compare(a, c, direction)
    # the three closure conditions of a total preorder
    if ab >= 0 and bc >= 0:
        assert ac >= 0
    if ab > 0 and bc > 0:
        assert ac > 0
    if ab == 0 and bc == 0:
        assert ac == 0
    return True


def ranking_is_total(a: List[int], b: List[int],
                     direction: List[int]) -> bool:
    assert len(a) <= MAX_KEY_LEN
    assert len(b) <= MAX_KEY_LEN
    assert len(direction) <= MAX_DIRECTION_LEN
    assert all(d in (-1, 1) for d in direction)
    ab = lex_compare(a, b, direction)
    ba = lex_compare(b, a, direction)
    assert ab in (-1, 0, 1)
    assert ba == -ab
    return True


def comparator_matches_successive_refinement(a: List[int], b: List[int],
                                             direction: List[int]) -> bool:
    assert len(a) <= MAX_KEY_LEN
    assert len(b) <= MAX_KEY_LEN
    assert len(direction) <= MAX_DIRECTION_LEN
    assert all(d in (-1, 1) for d in direction)
    left = lex_compare(a, b, direction)
    right = lex_compare_by_refinement(a, b, direction)
    assert left == right
    return True


def directions_only_rescale(a: List[int], b: List[int],
                            direction: List[int]) -> bool:
    """Mirroring each key and flipping its direction is the same ordering: the
    model depends on the keys only through their order."""
    assert len(a) <= MAX_KEY_LEN
    assert len(b) <= MAX_KEY_LEN
    assert len(direction) <= MAX_DIRECTION_LEN
    assert all(d in (-1, 1) for d in direction)
    n = max(len(a), len(b), len(direction))
    full = [direction[i] if i < len(direction) else 1 for i in range(n)]
    # negating the value and the direction together leaves every key unchanged
    mirrored_a = [-v for v in a]
    mirrored_b = [-v for v in b]
    mirrored_direction = [-d for d in full]
    left = lex_compare(a, b, direction)
    right = lex_compare(mirrored_a, mirrored_b, mirrored_direction)
    assert left == right
    return True


# --- how the R call sites are shaped ---------------------------------------

def _selftest() -> None:
    """A handful of concrete pairs, so the model is readable without a prover."""
    desc_then_asc = [1, -1]          # tier high-is-better, tie low-is-better
    assert lex_compare([1, 3], [1, 1], desc_then_asc) == -1
    assert lex_compare([2, 9], [1, 1], desc_then_asc) == 1
    assert lex_compare([1, 1], [1, 1], desc_then_asc) == 0


if __name__ == "__main__":
    _selftest()
    print("selftest ok")
