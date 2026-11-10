# Measurements — plant

Where plant's baselines live, what may re-bless them, how a build is validated, and how a speed claim is measured. [hazards.md](hazards.md) holds the invariants; [agents.md](../agents.md) the standing guide.

## Baselines

Each of these is compared bit for bit or at a stated tolerance by a test. None moves without a deliberate re-bless, in its own commit, with the before-and-after in the PR and, where the science changed, a `scientific_version` bump (agents.md §7).

| Baseline | What it pins | Read by | Re-blessed by |
|---|---|---|---|
| `tests/testthat/test_data/default_baseline.rds` | six runs that use no opt-in feature (FF16, K93, FF16 with a harvest, an invasion, a resume, TF24), recorded from the commit before #661 (`9f04309b`); compared only on macOS arm64, where it was recorded | `test-default-baseline.R`, through `default_baseline_cases()` in `helper-default-baseline.R` | `Rscript scripts/capture-default-baseline.R` against an installed build |
| `tests/testthat/reference/reference-gradient.tsv`, `reference-kinds.tsv` | the differenced reference the analytic leaf rows are refereed against, every step's reading, and which operating-point kinds each regime reached | `test-gradient-ladder-whole-run-difference.R` | `scripts/capture-reference-gradient.R`; read its header first, because the differenced implementation it used is gone |
| `tests/testthat/FF16_reference/` | FF16 against Daniel's original R implementation of the growth model | `test-strategy-ff16-reference-comparison.R` | not re-blessed: it is an independent implementation |
| `tests/testthat/_snaps/model-version.md` | each model's default `pars`, `control` and top-level science fields, against its `scientific_version` | `test-model-version.R` | `testthat::snapshot_accept("model-version")`, with the version bump |
| `tests/testthat/test_data/scenario_baseline.rds` | each TF24 hydraulic scenario's two classifications (`observed`, `persists`) and its offspring production at `SCENARIO_TOL` (1e-3) | `test-scenario-gateway.R`, only with `PLANT_RUN_SCENARIOS=1` | `make bless-scenarios`, with the before-and-after table in NEWS |

Three of these skip on CRAN (the default baseline, the drift guard, which is a snapshot, and the gateway), so a run with `NOT_CRAN` unset has no live gate on a version bump: check they ran before believing one.

## Validating a build

Run the suite against an installed plant, not the checkout: `R_LIBS=<lib> make test-installed`, and read the "loaded from" line it prints. agents.md §8 explains each setting it makes. Build against odelia and phylloptim at the shas pinned in `DESCRIPTION`'s `Remotes:`, and `rm -f src/*.o` first whenever a `LinkingTo` dependency changed (agents.md §6). Compare the pass count with a control run of the same invocation on a reference build, never with a remembered number.

The scenario gateway runs with `PLANT_RUN_SCENARIOS=1` and is needed for any change that touches TF24 or its environment.

## Measuring speed

Wall-clock on a Mac swings by up to 20% between rounds of the same binary, so a speed claim needs evidence that does not rest on it:

- **The compiled objects** of base and head, built the same way without `-g`: byte-identical objects (`cmp`) prove a comment-only or dead-code change; where bytes differ, hashes of each function's normalised disassembly say which functions moved. Always with a positive control: change one constant a translation unit reads, and confirm the comparison reports it.
- **Instructions retired** for a fixed script, from `/usr/bin/time -l`, which repeats to within about 0.3% where wall-clock drifts.
- **Interleaved medians**, base and head alternated in one session, each build in its own library and one `Rscript` per build, since two builds of plant cannot load into one R session.

The `profile-plant` skill has the profiling workflow. `run_plant_benchmarks()` times only FF16's run and schedule refinement; a change to TF24, to the reverse sweep or to anything phylloptim reaches needs those paths timed as well.

## Evaluation records kept in `notes/`

Each is a measured result, with the script that produced it. None runs on current code unless it says so.

| Record | What it decided |
|---|---|
| `notes/penman-monteith/step5-fick-vs-pm-evaluation.md` | the Penman-Monteith leaf energy balance changes leaf outputs enough, in the predicted places, to keep; re-run in October 2026 on plant `bd09645d` with `step5_fick_vs_pm.R` (its `Leaf(a = ...)` is now `phi_J`) |
| `notes/penman-monteith/sensitivity-results.md` | Penman-Monteith's answer depends on the boundary-layer resistance, net radiation and stand inputs nobody has calibrated: keep the core, default off; scripts `sensitivity_leaf.R`, `sensitivity_scm.R` |
| `notes/atls/gate-atls-vs-pm-evaluation.md` | the ATLS thermal-damage layer changes carbon gain materially and where predicted (July code); the method is `gate_atls_vs_pm.R` |
| `notes/atls/results-2026-07/` | the July viability runs, with `atls_demo_helpers-july.R` (July code) |
| `notes/scripts/soil_cascade/` | a standalone R replica of the soil cascade, needing no plant build, used for hazards 9 to 12 |
| `notes/scripts/continuous-patch-519/` | the continuous-patch pilot's fixed-point and long-patch runs |
