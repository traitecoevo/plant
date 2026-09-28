# Plan: the prototype's gradient access points on the reverse-mode engine

Status: **open, order proposed, not yet agreed** (2026-09-28). The engine is #648's; the access points are the ones the #553 prototype shipped ([its list](https://github.com/traitecoevo/plant/pull/553#issuecomment-5186131586)). #553 is a prototype and will not merge, so anything a downstream user relied on from it has to be rebuilt here.

## Where #648 stands against the prototype

| prototype access point | what it gives | prototype models | in #648 |
|---|---|---|---|
| `stand_gradient(scm, metrics, traits, species, feedback)` | metrics × traits Jacobian of emergent stand outputs | FF16, TF24f all metrics; TF24 offspring only | **Partly.** TF24 only. Three census metrics at the final time: `leaf_area`, `mass_above_ground`, `area_stem`. No `species` or `feedback` argument. The one mode is the whole-run total, the analogue of `feedback = "resident"`; there is no frozen-canopy mode. FF16 fails with an `Rcpp::as` type error. |
| `offspring_production_gradient(scm, traits, species)` | d(offspring production)/d(trait): the invasion or selection gradient | FF16, TF24, TF24f | **Missing.** Offspring production is not a census metric. |
| `stand_state_jacobian(scm, traits, species)` | per-cohort state × trait Jacobian, an escape hatch for any downstream metric | FF16, TF24 | **Missing.** `census_trait_tangent` computes a directional derivative internally, but no per-cohort state Jacobian is exposed. `stand_census_state_adjoint` is a test-only export and is not this. |
| `grow_individual_to_size_gradient(individual, sizes, size_name, env, traits)` | d(time to size)/dθ and d(state at that time)/dθ for one plant in a fixed environment | FF16, TF24f | **Missing.** |
| `birth_rate_gradient(scm, metrics, species)` | d(census)/d(birth rate) and d(R0)/d(birth rate) | FF16 | **Missing.** Birth rate is not among the 48 columns. |

## What each needs on this engine, and a proposed order

Ordered by what a calibration or selection-gradient user needs first, and by what the engine already nearly does.

1. **Offspring production as a census metric (TF24).** Fecundity is an ODE state, so the stand's offspring production at the final time is a function of the final state. That is the shape a census metric already has: a row in `TF24_Strategy::census_metrics()`, then `stand_gradient(metrics = "offspring_production")` returns it and `offspring_production_gradient()` becomes a thin wrapper. Smallest step, and it is the fitness quantity.
2. **FF16 (then K93) as `Censusable`.** FF16 is already templated on the scalar but instantiates no active path. It needs `census_metrics()`, `ad_parameters()` and trait names, plus an active instantiation in the gradient translation units. That compile cost is the one §6.2 of the review measured; it should be re-measured.
3. **Birth rate as a column.** It enters through the introduction density, which is the insertion map the sweep already transposes, so it is a parameter of `apply_insertion` rather than a new mechanism. It is needed for d(R0)/d(birth rate).
4. **Frozen-canopy feedback, the rare-mutant invasion gradient.** This is the selection gradient adaptive dynamics needs. It sweeps a `run_mutant` replay with the resident light held passive. The replay path is already odelia's recording, so the sweep is the same; the difference is which System is active. ⚠️ `run_mutant` uses `step_to`, which subdivides on rejection, and odelia's sweep would transpose a subdivided row as one step (review §5.1 defect 1). A row-kind refusal in odelia has to land first.
5. **`stand_state_jacobian()`.** It exposes `census_trait_tangent`'s forward tangent per cohort and state, instead of contracting it with a census metric.
6. **`grow_individual_to_size_gradient()`.** It needs a single `Individual` in a fixed environment, the implicit function theorem on the time-to-size root, and a tangent of the state there. It is independent of the stand machinery.
7. **TF24f as `Censusable`.** The prototype's TF24f resident census gated at patch lifetime ≳ 5, so expect the same stiffness here.

## Decisions needed before starting

- Which of these are required before #648 merges, and which can follow on their own branches.
- Whether `stand_gradient()` keeps the prototype's `feedback = c("frozen", "resident")` argument, or frozen gets its own entry point.
- The `lma` question from the review (§4.2): columns are partials at fixed hyperparameters. If the chain rule through `TF24_hyperpar` lands, it lands before FF16 copies the same column convention.
