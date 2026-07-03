# TF24/TF24F Scenario Evaluation Framework

Status: implemented. Branch `feature/tf24-scenario-framework`, based on
`hydraulic_forest_expectations` (PR #555).

## Context

The TF24 / TF24F hydraulic strategy models (epic
[#424](https://github.com/traitecoevo/plant/issues/424)) fail under certain
drought/soil conditions:
[#549](https://github.com/traitecoevo/plant/issues/549) (non-finite `psi_soil`
at residual soil moisture) and #550 (growth-dependent mortality blow-up). PR
[#555](https://github.com/traitecoevo/plant/pull/555) adds
`model_scenarious_hydraulic.csv`, a table of hydraulic scenarios described
qualitatively (trait bins + environment descriptors) with an **expected**
outcome per row (model runs / model fails).

**Goal:** a rerunnable, branch-agnostic framework that, for each scenario,
translates the qualitative descriptors into concrete TF24 parameters +
environment config, runs the SCM, classifies the run as success/failure,
compares to the expected outcome, and emits a provenance-stamped **scorecard**.
Many scenarios are expected to fail on the current model — those are the
improvement targets (e.g. the NSC storage work in PR #554). The framework
doubles as long-term **gateway-check** infrastructure: re-runnable across
branches so scorecards are directly comparable as features land.

The model scored on this base is pre-NSC, so the first scorecard is the
baseline; re-running after #554 (NSC) / #548 (hyperpar) merge shows movement.
Provenance metadata (git sha/branch/dirty, package version) is recorded in every
scorecard.

## CSV review findings (logical flaws in the PR #555 CSV)

1. **Exact duplicate rows.** Row 1 == Row 4 and Row 5 == Row 8. 8 rows → 6
   unique. Copy-paste artifact.
2. **Fully confounded design.** All acquisitive/mesic traits co-vary as one
   syndrome and all conservative/xeric traits as another; environment (rainfall +
   drainage + seasonality) is confounded with the syndrome. A failure cannot be
   attributed to any single driver — paired archetypes, not a factorial.
3. **Only one expected success** (Row 2), differing from a failure (Row 6) *only*
   in root distribution in a wet, mildly-seasonal environment. Physically
   questionable that deep roots in wet soil cause failure — flag for confirmation.
4. **Scenario names mislead.** "High rainfall plants in arid environment" means
   high-rainfall-*adapted traits* in an arid environment (maladaptation test),
   not a high-rainfall environment.
5. **`Psi_crit,b,c` is one free trait.** In TF24 `c`, `b`, `psi_crit` are all
   derived from `p_50`; the column really means "choose `p_50`".
6. **No matched/adapted controls, no aseasonal baseline.**
7. **Binned thresholds leave gaps / units unstated** (0.5–1.0 m undefined;
   "amplitude = mean" ⇒ trough hits ~0, must be clamped ≥ 0).
8. **Binary expectation only** — no graduated target/metric.
9. **Cosmetic:** filename typo `scenarious`; UTF-8 BOM.

**Resolution:** correct the CSV in place (single source of truth): remove
duplicates (1), add matched/adapted controls and one-trait sweeps (3, 6),
clarify names (4), rename + drop BOM (9). Confounding (2) and binary expectation
(8) are documented as interpretation caveats in the scorecard rather than
mechanically altered. Pre-correction rows recoverable via git history.

## Decisions

1. **Scope:** full gateway infra (engine + editable mapping + corrected CSV +
   seasonal rainfall driver + scorecard report + opt-in baseline-diff gateway
   test + runner script + `make` target).
2. **`g1_TF24`:** expose as a settable `TF24_Pars` field first (C++ header + yml
   + regenerate + rebuild), so all trait columns can vary.
3. **CSV:** corrected in place, single source of truth.

## Architecture (layered)

```
inst/scenarios/model_scenarios_hydraulic.csv    <- corrected in place (git mv from repo root)
inst/scenarios/scenario_mapping.csv             <- editable qualitative->quantitative lookup
R/scenario_eval.R                               <- core engine + driver
inst/reports/scenario_scorecard.Rmd             <- human-facing scorecard report
scripts/run_scenario_gateway.R                  <- standalone rerunnable gateway
tests/testthat/test-scenario-eval.R             <- unit tests for engine functions
tests/testthat/test-scenario-gateway.R          <- opt-in baseline-diff test
tests/testthat/test_data/scenario_baseline.rds  <- recorded baseline scorecard
```

Four separable concerns: **semantics** (mapping data file) →
**translate/run/classify** (pure functions) → **iterate/score/provenance**
(driver) → **present/gate** (report + test + script).

## Mapping: qualitative → quantitative (`inst/scenarios/scenario_mapping.csv`)

Long/tidy lookup, one row per (csv_column × level):
`csv_column, level, target, param, value, note`. `target` routes to `trait` /
`env` (Environment field) / `driver` (extrinsic driver). First-pass magnitudes
use a single symmetric multiplier **k = 2** off the *Eucalyptus saligna*
defaults (High = default×2, Low = default÷2). Per-trait overrides are one cell.

| Column | Level | target | param | value |
|---|---|---|---|---|
| g1_TF24 | High / Low | trait | `g1_TF24` | 15 / 3.75 |
| Ks | High / Low | trait | `K_s` | 2 / 0.5 |
| Vcmax_25 | High / Low | trait | `vcmax_25` | 192 / 48 |
| LMA | High / Low | trait | `lma` | 0.3958 / 0.0989 |
| Wood density | High / Low | trait | `rho` | 1216 / 304 |
| Huber value | High / Low | trait | `theta` | 4.285e-4 / 1.071e-4 |
| Psi_crit,b,c | More/Less negative psi | trait | `p_50` | 3.70 / 0.925 → derive `c,b,psi_crit` |
| Root distribution | top / bottom heavy | trait | `root_depth_shape_eta` | 0.05 / 0.8 (verify sign) |
| Rainfall | <0.5m / >1m | driver | `rainfall_mean` | 0.4 / 1.5 (m/yr) |
| Drainage rate | High / Low | env | `K_sat` | 163.04 / 16.3 |
| Amplitude | =mean / half | driver | `rainfall_amp_frac` | 1.0 / 0.5 |

### Correctness details

- **Vulnerability curve.** `TF24_hyperpar` does **not** derive `c/b/psi_crit`
  from `p_50` today ([R/tf24.R](../R/tf24.R)); `prepare_strategy()` uses
  `pars.c/b/psi_crit` verbatim ([src/tf24_strategy.cpp](../src/tf24_strategy.cpp)).
  PR #548 (open) is the proper fix. **Interim:** `scenario_to_config` computes
  `c/b/psi_crit` in R from the chosen `p_50` and emits all four as traits. Marked
  for removal once #548 merges. Formulas
  ([tf24_strategy.h:62-64](../inst/include/plant/models/tf24_strategy.h)):
  `c = log(log(0.5)/log(0.12)) / (log(p_50) - log(5.16))`;
  `b = p_50 / (-log(0.5))^(1/c)`; `psi_crit = b * log(1/0.05)^(1/c)`.
- **`g1_TF24` exposure.** Currently in the "not exposed to R" block of
  `TF24_Strategy`, not `TF24_Pars`. Moved into `TF24_Pars` + added to
  `inst/RcppR6_classes.yml`, then `make RcppR6 && make full_compile`.

## Core engine — `R/scenario_eval.R`

```r
read_scenario_table(path = system.file("scenarios","model_scenarios_hydraulic.csv", package="plant"))
read_scenario_mapping(path = system.file("scenarios","scenario_mapping.csv", package="plant"))
scenario_to_config(row, mapping)   # -> list(traits, env, driver, expected)
build_scenario(config, max_patch_lifetime=100, ctrl=control(), birth_rate=list(1))  # -> list(p,env,ctrl)
classify_scm_run(p, env, ctrl)     # tryCatch + warn=1; -> list(status, offspring_production, finite, error_message, ...)
evaluate_scenario(row, mapping, ...)   # one-row scorecard tibble
run_scenarios(scenarios=read_scenario_table(), mapping=read_scenario_mapping(), ...)  # scorecard + attr metadata
scenario_run_metadata()            # git sha/branch/dirty, packageVersion, R, timestamp, platform
scenario_summary(scorecard)        # n, n_match, match_rate, expected-fail/success met
scenario_generate_report(scorecard, output_file=, ...)
```

Classification leans on the C++ layer failing fast on non-finite state; classify
on *any* error + finiteness, store the message as a diagnostic only. An
expected-failure that fails is a **match**.

`observed` is binary (`success` = finite, positive offspring; else `failure`),
but each row also carries a richer `outcome`: `persisted` (finite, offspring >
0), `extinct` (finite run, no offspring) or `crashed` (numerical failure — a
thrown error or non-finite state). This matters because the hydraulic/NSC work
(#549, #554) targets *crashes*, whereas `extinct` is often correct model
behaviour. The scorecard surfaces both so a "failure" is never mistaken for a
crash when it is really an extinction.

**Parallelism.** `run_scenarios(..., workers = N)` uses **fork-based**
parallelism (`parallel::mclapply`). This is deliberate: forked workers inherit
the current session's namespace and compiled `.so`, so it works when `plant` is
loaded for development via `pkgload::load_all` / `devtools::load_all`. A PSOCK /
`future::multisession` cluster spawns fresh R sessions that see only the
*installed* package — they would silently run the stale (or missing) build, so
they are not used. Forking is unavailable on Windows, where the run falls back
to sequential. Scenario runs are deterministic (no RNG), so parallel and
sequential scorecards are identical (verified).

## Gateway wiring

- `scripts/run_scenario_gateway.R` — `pkgload::load_all(".", recompile=FALSE)`,
  `run_scenarios()`, print summary, `saveRDS()` scorecard to an artifact path.
  `make scenarios` target.
- `tests/testthat/test-scenario-gateway.R` — opt-in
  (`PLANT_RUN_SCENARIOS=1`, `skip_on_cran()`). Baseline diff against
  `scenario_baseline.rds`, not "all pass": detects regressions AND improvements.
- `tests/testthat/test-scenario-eval.R` — fast unit tests (no full SCM).

## Verification

1. `make RcppR6 && make full_compile && make roxygen` — `g1_TF24` exposed;
   `TF24_Strategy()$pars$g1_TF24 == 7.5`.
2. `run_scenarios()` returns a scorecard with `expected/observed/match/
   offspring_production/error_message` and git provenance.
3. Spot-check an expected-failure row throws (cf. `test-strategy-tf24.R:333`) and
   an expected-success row returns finite positive `offspring_production`.
4. `Rscript scripts/run_scenario_gateway.R` prints summary + writes RDS;
   `scenario_generate_report()` renders HTML.
5. `PLANT_RUN_SCENARIOS=1 devtools::test(filter="scenario")` green.
6. Record baseline scorecard; note match_rate as the pre-NSC target to beat.

## Baseline results (this build: pre-NSC, `max_patch_lifetime = 100`)

Match rate **4/8 (50%)**. Recorded in
`tests/testthat/test_data/scenario_baseline.rds`.

| id | scenario | expected | observed | outcome | match |
|---|---|---|---|---|---|
| S01 | mesic / arid / shallow | failure | failure | extinct | ✓ |
| S02 | mesic / arid / deep | failure | failure | extinct | ✓ |
| S03 | xeric / wet-seasonal / shallow | success | failure | **crashed** | ✗ |
| S04 | xeric / wet-seasonal / deep | failure | success | persisted | ✗ |
| S05 | xeric / wet-extreme-seasonal / shallow | failure | failure | crashed | ✓ |
| S06 | xeric / wet-extreme-seasonal / deep | failure | success | persisted | ✗ |
| S07 | mesic / matched wet (control) | success | failure | extinct | ✗ |
| S08 | xeric / matched arid (control) | success | success | persisted | ✓ |

Targets surfaced by the framework:

- **S03 is the priority:** a *crash* (numerical failure) in a scenario expected
  to *succeed*. This is exactly the #549/#554 failure class and should be gone
  once NSC lands — re-run to confirm.
- **S07:** a matched mesic control goes extinct (finite run, no offspring).
  Likely a parameter/birth-rate tuning issue rather than a numerical bug; worth
  checking the mesic trait bundle and `birth_rate`.
- **S04, S06:** the model runs where the CSV expected failure — either it is
  more robust than the authors assumed, or the `k = 2` mapping is not extreme
  enough. Revisit the mapping magnitudes, or confirm the expectations.

## Open dependency

PR #548 (`tf24_hyperparameterisation`, open) is the proper home for the
`p_50 → c/b/psi_crit` derivation. When it merges, delete the interim R-side
block and pass `p_50` alone.
