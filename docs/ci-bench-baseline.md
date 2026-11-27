# Capturing the CI benchmark baseline

`dev/bench.R --gate` compares a run against `dev/bench-baseline.csv`, keyed by
runner class (`TIDYMARKERS_BENCH_HOST`; see ADR 0014). The local class is
already in the file. The GitHub class `github-ubuntu-latest` cannot be measured
from a workstation, so its rows are captured on CI and committed by a human.

## Steps

1. On GitHub, run the **R-CMD-check** workflow with *Run workflow*
   (`workflow_dispatch`). The `bench-baseline` job runs on `ubuntu-latest`,
   sets `TIDYMARKERS_BENCH_HOST=github-ubuntu-latest`, and runs
   `Rscript dev/bench.R --update`.
2. Download the job artifact `bench-baseline-github-ubuntu-latest`. It contains
   the complete `dev/bench-baseline.csv` with one block per host: the job keeps
   existing rows and replaces only its own class.
3. Replace `dev/bench-baseline.csv` with the artifact and commit it.
4. Check locally:

   ```sh
   Rscript dev/bench.R --gate
   ```

   This host still compares against its own block; the CI job compares against
   `github-ubuntu-latest` on the next run.

## Why the job cannot commit

The workflow has read-only repository permission and the default branch is
protected. Uploading the file as an artifact keeps the measurement honest: a
human reviews the diff and commits it, instead of CI silently moving its own
baseline.

## Failure modes

- **No artifact**: the job failed before `--update`; read the job log, not this
  file.
- **Artifact has one host only**: it was produced before
  `dev/bench.R --update` learned to keep other hosts. Re-run on current `main`.
- **Gate says `not_run` locally**: the file has no block for this host; that is
  expected until someone runs `--update` here. It is not a pass.
