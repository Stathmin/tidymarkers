# =============================================================================
#  dev/check-criteria.R  --  every criterion must be able to fail (ZW1)
# =============================================================================
#  A criterion that cannot be driven to `fail` is a criterion nothing proves,
#  and a new one used to be caught only after the fact, by a test.  This is the
#  pull-request gate for that rule: it walks the registry (not a hand-kept list)
#  and requires every `gate` / `check` id to be
#
#    * exempt, with a non-empty reason (a genome / BLAST database / population
#      the in-memory suite does not carry), or
#    * backed by a fixture that names a REAL test in this repo.
#
#  The ledger (which id is exempt, which test makes it fail) is owned by
#  tests/testthat/helper-criterion-fixtures.R -- that file is the one source of
#  truth, so this gate and test-falsifiable.R cannot disagree about coverage.
#  Neither copies the table.
#
#  Run:  Rscript dev/check-criteria.R   (also run by dev/check.R)
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
setwd(root)

suppressMessages(pkgload::load_all(quiet = TRUE))
source(file.path("tests", "testthat", "helper-criterion-fixtures.R"))

ledger <- zw1_ledger()
problems <- zw1_ledger_violations(ledger)

if (length(problems)) {
  message("criterion fixtures (ZW1): ", length(problems), " problem(s)")
  for (p in problems) {
    message("  - ", p)
  }
  message(
    "\nEvery id in criterion_spec(KASP_LIMITS) with role 'gate' or 'check' must be ",
    "falsifiable: add a fail fixture to tests/testthat/test-falsifiable.R or an ",
    "exemption with a reason (tests/testthat/helper-criterion-fixtures.R)."
  )
  quit(status = 1L)
}

message(
  "criterion fixtures (ZW1): ",
  length(ledger$registered), " registered (gate + check) = ",
  nrow(ledger$fixtures), " fixture-backed + ",
  length(ledger$exempt), " exempt; every citation resolves to a real test"
)
