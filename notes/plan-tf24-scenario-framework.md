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
scenario_summary(scorecard)        # two axes: viability (n_ran/n_crashed) + persistence (n_persists); CSV agreement reported but not headline
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

### Scoring: two axes, no single headline (#572, decided 2026-07-30)

`observed == "success"` tests `finite && total > 0`, which is far weaker than it
looks: with `birth_rate = 1` (what `build_scenario` sets) `offspring_production`
**is** R0, and at the blessed baseline five of the eight scenarios return R0
between 2e-15 and 6e-14 — numerically extinct — and were all recorded as
`persisted`. All eight "succeeded", none of the five expected failures failed,
and only one actually replaced itself. The gateway was returning no signal.

**Decision: options 3 + 1 of #572** — report both axes, and make persistence the
primary *ecological* metric, without rewriting the CSV or re-blessing the
baseline:

- **Numerical viability** (`n_ran` / `n_crashed` / `viability_rate`) keeps
  meaning what the CSV's "Model failure" means (#549, #550): did the model run.
  This is the gate the hydraulic/NSC work targets, and `status` / `outcome` are
  unchanged, so the blessed-baseline diff is still defined on the same columns.
- **Ecological persistence** (`n_persists` / `persistence_rate`) is judged at
  R0 >= 1 via `persists_at()`, added in #570.
- **CSV agreement** (`n_match` / `match_rate`) is still reported but is no longer
  the headline, in `scenario_summary()`, `scripts/run_scenario_gateway.R` and the
  scorecard report. Rationale: with the crashes fixed it mostly measures how well
  the CSV's crash predictions have aged, not model quality.

Deliberately **not** done: option 2 (rewriting the CSV's `Expectation` column to
be about persistence). That would discard the crash-prediction record and force a
re-blessing. If it is ever taken, the old expectations are recoverable from git
history and that should be said here.

**Threshold caveat.** R0 >= 1 is the textbook criterion but is evaluated here at a
*fixed* birth rate, so it carries the same density-dependence artefact that made
TF24t-vs-TF24 R0 gaps look enormous in unrelated work. All scenarios use
`birth_rate = 1`, so the comparison is internally consistent; if the gateway ever
varies birth rate, `persists` needs revisiting. `persists_at(total, finite,
threshold = 1)` takes the threshold as an argument for exactly that reason.

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

## PR #548 integration (merged into this branch)

PR #548 (`tf24_hyperparameterisation`) has been merged in. It moves the
hydraulic derivations into `TF24_hyperpar`, which changes how the CSV columns
map:

- The interim R-side `p_50 → c/b/psi_crit` derivation has been **removed** from
  `scenario_eval.R`; the hyperpar owns it.
- `p_50/c/b/psi_crit` are now derived from **`K_s`**, and **`g1_TF24`** from
  **`rho`** (`R/tf24.R` `extra` block). Both are therefore **removed as input
  traits** from `scenario_mapping.csv` — passing them would trip the hyperpar's
  overwrite guard. The CSV's pairings stay consistent (low `K_s` → more-negative
  `p_50`; high `rho` → lower `g1_TF24`), so the intended physiology is preserved
  via the `K_s` and `rho` mappings. Under #548's default the `g1_TF24` scaling
  exponent is 0, so `g1_TF24` is effectively constant unless recalibrated.
- Conflicts (all in the shared `g1_TF24`/`hk_s`/pars region) were resolved in
  favour of #548, which also exposes `g1_TF24` in `TF24_Pars` (superseding the
  standalone exposure this branch had added) and drops `hk_s` as a parameter.

## Caching: when to rerun

Scenario runs are deterministic, so a result is a pure function of
`(resolved config, max_patch_lifetime, model)`. `run_scenarios(cache = path)`
keys each row by a content hash of those and reruns only what changed:

- **Model change** (any C++ recompile or package R-source edit) moves
  `scenario_model_fingerprint()` → **all** rerun.
- **Mapping / scenario edit** changes only the affected rows' config hash → only
  those rerun.
- **Unrelated edits** (report, runner) → full reuse.

The key covers **all** inputs that affect a result: resolved config,
`max_patch_lifetime`, every `ctrl` field, and the model fingerprint. The
fingerprint is intentionally broad (package version + compiled `.so` md5 + md5
of all `R/*.R` + the scenario CSVs): the cache errs toward rerunning, since a
stale gateway result is worse than a redundant run.

## Baseline results (post-#548, `max_patch_lifetime = 100`): 5/8

| id | expected | observed | outcome | note |
|---|---|---|---|---|
| S01 | failure | failure | crashed | ✓ |
| S02 | failure | failure | crashed | ✓ |
| S03 | success | success | persisted | ✓ (**#548 fixed the earlier crash**) |
| S04 | failure | success | persisted | ✗ (runs where failure expected) |
| S05 | failure | failure | crashed | ✓ |
| S06 | failure | failure | crashed | ✓ |
| S07 | success | failure | crashed | ✗ (matched mesic control crashes — target) |
| S08 | success | failure | crashed | ✗ (matched xeric control crashes — target) |

#548 raised the match rate (4/8 → 5/8) and, notably, fixed S03. The two matched
controls (S07, S08) now crash — the clearest remaining targets, alongside S04.

### Superseded by the two-axis scoring (#572)

The 5/8 above is **stale**, and #557's body still quotes it. The crash fixes in
#546/#552/#554 landed and *every* scenario now runs. Re-measured on this branch
merged with `develop` at d16d1357 — i.e. after #574, #585, #590 and the
phylloptim port (#591) — at `max_patch_lifetime = 100`, with the head-of-#570
column kept alongside so the drift is visible:

| id | expected | observed | R0 (at #570) | R0 height | R0 birth date | persists |
|---|---|---|---|---|---|---|
| S01 | failure | success | 2.95e-01 | 2.59e-01 | **1.49e+00** | ✓ |
| S02 | failure | success | 4.98e-02 | 4.20e-02 | 9.96e-02 | ✗ |
| S03 | success | success | 6.49e-13 | 6.49e-13 | 5.50e-12 | ✗ |
| S04 | failure | success | 6.27e-14 | 6.27e-14 | 4.75e-13 | ✗ |
| S05 | failure | success | 2.33e-14 | 2.35e-14 | 5.89e-13 | ✗ |
| S06 | failure | success | 3.68e-14 | 3.68e-14 | 3.26e-13 | ✗ |
| S07 | success | success | **2.93e+01** | **3.02e+01** | **1.42e+03** | ✓ |
| S08 | success | success | 2.31e-15 | 2.17e-15 | 3.23e-14 | ✗ |

Two separate movements are folded into that table.

**#585's forward-model fixes** moved the three scenarios not already at machine
zero by a few percent each (the `at #570` → `height` columns). No verdict
changed.

**The density coordinate** (`height` → `birth date`) moved everything, by 2.4x
(S02) to 47x (S07), and moved S01 across R0 = 1. That is the change that
matters, and it is a correctness fix rather than a tuning choice — see below.

So CSV agreement is **3/8**, not 5/8 — and the drop is *good news misreported*:
the matched-control crashes that #557 called the clearest targets are gone. That
inversion (fixing the model lowers the score) is what showed the metric was
pointing the wrong way round and prompted #572.

Read on the two axes the gateway now reports: **numerical viability 8/8**, and
**persistence 2/8**. Both are informative, and neither is a "match rate".

### Density coordinate: the gateway runs in birth date (#590, decided 2026-08-05)

`scenario_control()` supplies `control(node_density_in_birth_date = TRUE)` and is
the default `ctrl` for `build_scenario`, `classify_scm_run`, `evaluate_scenario`
and `run_scenarios`. The package default (`FALSE`, height) is deliberately left
alone; only the gateway opts in.

Why the gateway and not the package: every scenario here is TF24, and TF24 is
exactly the model #590 identifies as the one the coordinates disagree on. The
compression term is the total derivative of growth along a cohort's own
trajectory, which equals `dg/dh` only when growth is a function of size. TF24's
reserve gate (#517) breaks that — the finite-difference probe moves height at
fixed *absolute* carbon, shifting the reserve fraction, whereas a real cohort
grows at roughly constant reserve fraction. So the height coordinate is not a
coarser approximation of the right derivative; it is an accurate derivative of
the wrong quantity. The diagnostic #590 gives is that refining the node schedule
does not close the gap for TF24, where for FF16 and K93 (size-only growth) the
two coordinates converge at ~2nd order.

Running the hydraulic gateway in the coordinate that is wrong for hydraulics
would score the model on an artefact, so the gateway opts in.

**What this cost, and what it caught.** Birth date is also ~2x faster here (36 s
against 70 s for the eight scenarios). And it exposed a hole in the baseline
guard: across a 47x swing in R0 and a persistence flip, `observed` did not change
on a *single* scenario, because it tests `finite && total > 0` and everything has
satisfied that since the crash fixes landed. `test-scenario-gateway.R` now diffs
`persists` alongside `observed`.

**The baseline is re-blessed** under birth date — reversing the "not re-blessed"
position taken earlier in this note, which was correct for a pure reporting
change and is not correct for a coordinate change that moves a verdict.
