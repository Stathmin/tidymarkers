# Operator decisions

Dated decisions that close a blocker with a measurement or a policy, rather
than with code alone. They are not ADRs: an ADR records a module or dependency
choice, while these records say what the operator decided when the queue asked
for a machine or a gate that this host does not have.

| Record | Decision | Status |
| --- | --- | --- |
| [0001](0001-16gb-fixpur.md) | A 16 Gb fixpur is a disk job on per-chromosome shards, not a 21.5 GB RAM job; the actual run stays an operator action | accepted (measured proxy) |
| [0002](0002-ci-strict-gate.md) | CI fails when `styler`, `lintr`, `testthat` or the test files are absent; the local gate runs on the release leg | accepted |
