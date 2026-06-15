# Issue #408 — Move cohort refining into C++ and simplify the SCM

Tracking doc for the refactor described in
[traitecoevo/plant#408](https://github.com/traitecoevo/plant/issues/408).

**Goal:** the SCM should be as simple as possible and no longer rely on R-side
calculations for any cohort refining. Move the cohort-splitting / schedule
refinement algorithm into C++, and collect reproduction/error information at the
node level (as it is produced) rather than reconstructing it after the run.

**Agreed approach:** phased PRs; the refinement loop lives as a method on
`SCM` (`SCM::refine_schedule()`); the final R API is a breaking cleanup that
consolidates onto a single `run_scm()` entry point.

---

## 1. What is being calculated today

Two intertwined machines are smeared across R and C++.

### 1a. The cohort-refinement loop (mostly R today)

[`build_schedule()`](../R/build_schedule.R) drives the loop. Each of up to
`schedule_nsteps` iterations:

1. Calls [`run_scm_error()`](../R/scm_support.R), which constructs a fresh `SCM`
   and **manually single-steps** it via `scm$run_next()`. After each step it
   records, *per species that got a node this step*, the vector
   `compute_competition_effect_error_by_node_for_species_i(idx)`.
2. Those per-step vectors are padded into a matrix and reduced with `max` per
   column → a **per-node competition (LAI) error**.
3. At the end it reads `scm$net_reproduction_ratio_errors` → a **per-node
   reproduction error**.
4. `total[node] = max(competition_error, reproduction_error)`.
5. `build_schedule` flags `split = total > schedule_eps`, then
   [`split_times()`](../R/build_schedule.R) bisects the interval *below* each
   flagged node (upwind scheme), re-sorts, and the loop repeats.

Both error signals use [`local_error_integration`](../src/util.cpp) (estimates
how much the middle point contributes to a trapezium integral):

- **Competition error** over node *heights*:
  `Species::r_compute_competition_effect_by_nodes_error`, scaled by total patch
  competition. Must be sampled *during* the run (running max across steps) —
  this is the only reason `run_scm_error` re-implements the step loop in R.
- **Reproduction error** over node *introduction times*:
  `SCM::r_net_reproduction_ratio_errors`, computed once at the end, scaled by
  total offspring.

### 1b. The reproduction accounting (C++, but with post-hoc lookups)

Per-node lifetime fitness is already integrated inside the ODE:
[`Node`](../inst/include/plant/node.h) accumulates
`offspring_produced_survival_weighted_dt = fecundity · survival_individual ·
pr_patch_survival / pr_patch_survival_at_birth`, surfaced as `node.fecundity()`.

[`SCM::net_reproduction_ratio_by_node_weighted`](../inst/include/plant/scm.h)
then weights post-hoc:

```cpp
... *= patch.survival_weighting->density(times[i]) * parameters.strategies[species_index].S_D;
```

It re-derives `times[i]` from `node_schedule.times(species)` and re-queries
`survival_weighting->density(t)` — **both knowable at node birth**. The node
already stores `pr_patch_survival_at_birth`, and for the Weibull regime
`density(t) = p0 · pr_survival(t)`, so `density(t_birth) =
p0 · pr_patch_survival_at_birth`. Only `p0` and `S_D` (constants) are missing,
and the introduction time isn't stored on the node at all.

`offspring_production` and `net_reproduction_ratios` are the same
trapezium-over-introduction-times integral, differing only in `scalars`
(birth-rate vs 1.0).

### 1c. Key insight

Every "look it up later" path (introduction time, patch-density weight) can be
replaced by **recording the value on the `Node` at introduction**. Once a node
owns `introduction_time` and `density_at_birth`, `Species` can produce both the
weighted-fitness vector and the integration x-axis itself, and `SCM` no longer
needs `node_schedule` or `survival_weighting` for any reproduction calc. The
refinement loop can then run entirely in C++, because the only thing forcing it
into R (per-step competition-error sampling) is just a running max that
`SCM::run()` can maintain directly.

---

## 2. Phased plan

### Phase 1 — Node-level bookkeeping (behaviour-preserving) — ✅ DONE

Nodes own their introduction time and disturbance weight; SCM stops reaching
back into `node_schedule` / `survival_weighting`.

Implemented:
- `node.h`: added `node_introduction_time` + `patch_density_at_birth` members,
  `set_introduction(time, patch_density)`, `introduction_time()`,
  `patch_density()`, and `weighted_fecundity(S_D)`.
- `patch.h` `introduce_new_nodes`: stamps each new node with `time()` and
  `survival_weighting->density(time())` via `Species::stamp_new_node`.
- `species.h`: added `stamp_new_node`, `net_reproduction_ratio_by_node_weighted`
  (applies `density_at_birth · S_D`), and `node_times`.
- `scm.h`: reproduction methods (`net_reproduction_ratio_by_node_weighted`,
  `net_reproduction_ratio_for_species`, `offspring_production`,
  `net_reproduction_ratios`, `r_net_reproduction_ratio_for_species`,
  `r_net_reproduction_ratio_errors`) now use `species.node_times()` and the
  species-weighted vector; all `node_schedule.times()` /
  `patch.survival_weighting->density()` lookups removed.
- No `RcppR6_classes.yml` change needed (no R-exposed signatures changed);
  `make compile` only.
- Verified green: test-scm, test-scm-support, test-schedule-build (141/186),
  test-strategy-ff16 + ff16-reference-comparison, k93, tf24, mutant,
  tidy-outputs, patch, species, all stochastic, environment-TF24.

- **`node.h`**: add `introduction_time` and `density_at_birth`; set them at
  introduction (e.g. `set_introduction(t, density)`). Add
  `weighted_fecundity(double S_D)` = `fecundity() · density_at_birth · S_D` and
  an `introduction_time` getter.
- **`patch.h`** `introduce_new_nodes`: pass `time()` and
  `survival_weighting->density(time())` into each new node.
- **`species.h`**: add `node_times()`; apply the `density_at_birth · S_D`
  weighting at node level (`S_D` reachable via `strategy`).
- **`scm.h`**: rewrite `net_reproduction_ratio_by_node_weighted`,
  `net_reproduction_ratio_for_species`, `offspring_production`,
  `net_reproduction_ratios` to use `species.node_times()` + the species-weighted
  vector. Drop the `patch.survival_weighting->density(...)` /
  `node_schedule.times(...)` dependencies.
- **Verify**: `make compile` + `devtools::test()`; the regression in
  `test-scm.R` and the FF16 reference baselines must stay green.

### Phase 2 — Error collection inside `SCM::run()` — ✅ DONE

A single run yields the refinement signal; delete the R-side step loop.

Implemented:
- `scm.h`: added `collect_errors` flag and `competition_error_by_node`
  (per-species running max, init `-Inf`, NA-skipping → mirrors
  `apply(., 2, max, na.rm=TRUE)`). `run()` folds in
  `collect_competition_errors(added)` after each `run_next()`.
  `combined_node_errors()` returns the per-node `max(competition, reproduction)`
  error — the exact signal `run_scm_error()$err$total` assembled in R.
- `RcppR6_classes.yml`: exposed `collect_errors` (field) and
  `combined_node_errors` (getter); `make RcppR6 && make full_compile`.
- Verified `combined_node_errors` == old `run_scm_error()$err$total` for FF16
  single / two-species / refined schedules (exact). Added durable regression
  `test-scm.R::"combined_node_errors collected in C++ matches per-step
  assembly"` (inlines the assembly so it survives Phase 4).
- `run_scm_error` is now redundant; removed in Phase 4.

- **`scm.h`**: add a `collect_errors` flag and per-species running state:
  `competition_error[species][node]` updated as an element-wise max after each
  `run_next()` (pad for newly added nodes). At `complete()`, compute the
  per-node reproduction error. Add `combined_node_errors()` → per-species
  `max(competition_error, reproduction_error)`, replacing what `run_scm_error`
  assembled in R.
- Makes `run_scm_error` redundant.
- **Verify**: temporary R test asserting `combined_node_errors()` matches the
  old `run_scm_error()$err$total` for an FF16 case.

### Phase 3 — Refinement loop + `split_times` in C++ — ✅ DONE

`SCM::refine_schedule()` owns the adaptive loop.

Implemented:
- `scm.h`: added `Control control` member (for `schedule_eps` /
  `schedule_nsteps`), a static `split_times(times, split)` (upwind bisection:
  insert `0.5*(t[j]+t[j-1])` for each flagged node), and `refine_schedule()`
  which loops `run()` (with `collect_errors`) → flag `combined_node_errors > eps`
  → bisect → `node_schedule.set_times`, up to `schedule_nsteps`, then records
  the refined `node_schedule_times` + `ode_times` back into `parameters`.
- `RcppR6_classes.yml`: exposed `refine_schedule`; `make RcppR6 && full_compile`.
- Verified C++ `refine_schedule` reproduces R `build_schedule` exactly (refined
  times, ode_times, offspring_production) for FF16 single / two-species /
  regression-case (186 nodes) and K93. Added durable test
  `test-schedule-build.R::"C++ refine_schedule matches R build_schedule"`.

- **`scm.h`**: add `refine_schedule()` looping up to `control.schedule_nsteps`:
  `reset()` → `run()` with `collect_errors` → break if all
  `combined_node_errors ≤ schedule_eps`, else bisect flagged intervals and
  `node_schedule.set_times(...)`. Port `split_times` as a C++ helper (upwind
  bisection: insert `t[i] - dt[i-1]/2`, keep sorted, never split the last
  interval). Reuse the single SCM instance rather than reconstructing.
- Keep `parameters.node_schedule_times` / `ode_times` in sync so the refined
  `Parameters` stays self-describing.
- Expose `refine_schedule` (+ accessors) in `RcppR6_classes.yml`;
  `make rebuild`.
- **Verify**: the `build_schedule` regression in `test-schedule-build.R`
  (141 / 186 length expectations) must reproduce before deleting the R loop.

### Phase 4 — Breaking R-API cleanup — ⏳ DONE (uncommitted, awaiting review)

One entry point; old names removed.

Implemented (left uncommitted for review):
- `R/scm_support.R`: `run_scm(p, env, ctrl, refine_schedule = FALSE,
  collect = FALSE, use_ode_times = FALSE)`. Returns the tidied results list when
  `collect = TRUE` (with refined `parameters` as `p`), otherwise the `SCM`
  object. `refine_schedule = TRUE` calls `scm$refine_schedule()`.
- Removed `run_scm_collect`, `run_scm_error`, and `R/build_schedule.R`
  (`build_schedule` + `split_times`). NAMESPACE now exports only `run_scm`;
  `man/build_schedule.Rd` deleted, `man/run_scm.Rd` updated.
- Updated callers/tests: `R/benchmark.R`, doc refs in `R/ff16.R` /
  `R/tidy_outputs.R`; tests `test-schedule-build.R` (rewritten for new API),
  `test-strategy-ff16.R`, `test-tidy-outputs.R`, `test-environment-TF24.R`,
  `test-scm-support.R`, `test-mutant.R`, `test-strategy-tf24.R` (comment).
- Full suite: 0 failures, 1832 pass, 1 pre-existing skip.

Vignettes:
- Mechanically updated (no `.orig`, edited directly): `example_analysis`,
  `plant`, `strategy_new`, `emergent`, `self_thinning`, `patch`,
  `models/strategy_K93`, `models/AWRA_soil_water_model`.
- `methods/node_spacing.Rmd.orig` rewritten to describe the new C++ approach
  (`run_scm(refine_schedule=TRUE)` / `SCM::refine_schedule()` /
  `combined_node_errors`), fixed the previously-broken chunks, and verified the
  runnable chunks execute. **`methods/node_spacing.Rmd` (generated) must be
  re-knit from the `.orig`** (do not hand-edit it).
- `methods/solving_dynamics.Rmd`: prose references updated.

Follow-ups for the user:
- Re-knit `methods/node_spacing.Rmd` from its `.orig`.
- Downstream **`plant.assembly`** uses `build_schedule()` and
  `run_scm_collect()` (`R/community_plant.R`, `scripts/example/ESA.Rmd`) — these
  break and need updating in that repo:
  `build_schedule(p, ctrl=ctrl)` -> `run_scm(p, ctrl=ctrl,
  refine_schedule=TRUE)$parameters`; `run_scm_collect(x)` ->
  `run_scm(x, collect=TRUE)`.

- Consolidate onto
  `run_scm(p, env, ctrl, collect = FALSE, refine_schedule = FALSE)` in
  `scm_support.R`: `refine_schedule` delegates to `scm$refine_schedule()`;
  `collect` returns tidied history + reproduction outputs.
- **Remove** `run_scm_collect`, `run_scm_error`, and R `build_schedule` /
  `split_times`.
- Update all callers/tests: `test-scm-support.R`, `test-schedule-build.R`,
  `test-tidy-outputs.R`, `test-mutant.R`, the strategy tests, `benchmark.R`,
  and doc refs in `tidy_outputs.R`. **Check downstream `plant.assembly`**, which
  likely calls these.
- `make rebuild` + `make roxygen` + full `devtools::test()`.

### Cross-cutting notes

- The stochastic path also exposes `offspring_production`
  (`R/stochastic.R`) — Phase 1 touches shared classes, so re-run stochastic
  tests.
- FF16 reference baselines in `tests/testthat/FF16_reference/` are the safety
  net; if any phase shifts numerics beyond tolerance, that's a bug, not an
  expected regen.
