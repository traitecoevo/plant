## Plant (development version)

### Breaking changes

These change the R-facing interface and require updating downstream code. Each
entry gives the `old -> new` migration; the `plant-update-interface` skill
(`.claude/skills/plant-update-interface/`) reads this section to migrate
products using plant.

* **`TF24_Environment`'s `atm_kpa` driver now defaults to 101.3 kPa, not 100.5.**
  `scientific_version` for TF24 goes 7 -> 8, TF24f 7.1 -> 8.1. Migration: none
  required, but read this if you have TF24 results on disk.
  * The leaf model's ppm -> Pa conversion was the hard-coded constant `0.1013`, which
    is 101.3 kPa in disguise (`1e-6 * 101300 Pa`). The driver said 100.5, so the
    model's conductance side responded to 100.5 while Gamma*, Kc, Ko, Km and the ci
    root-find bounds silently assumed 101.3. Deriving the conversion from `atm_kpa`
    (phylloptim #15) made the model self-consistent and turned that disagreement into a
    **+2.4%** shift in output.
  * 100.5 came from `34d46ac2` (#446), an interface refactor that does not mention
    atmospheric pressure, with no recorded rationale, while every leaf-level test used
    101.3. Pinning the driver to the value the rest of the model already assumed
    reduces the **net** movement of this whole branch to **-0.20%**, and restores
    every pinned baseline to develop's own values -- including the seeded stochastic
    TF24 counts, which match exactly.
  * `atm_kpa` remains a driver: set it per site if you are modelling altitude.
    `env$extrinsic_drivers_set_constant("atm_kpa", <kPa>)`.

* **The collar bracket is now clamped to `root_psi_crit`** (phylloptim #24, #584).
  `scientific_version` for TF24 goes 6 -> 7, TF24f 6.1 -> 7.1. The clamp compared a
  magnitude against a signed potential, so it could never bind and the solver
  optimised over a root-collar potential the root system cannot supply. The window is
  **1.2 MPa wide at the defaults** (`psi_crit` 7.0855 vs `root_psi_crit` 5.8703), so
  water-limited runs change; no mesic run does, and the standard SCM scenario is
  bit-identical. No migration: no name or signature changes.

* **Water potential now has one representation everywhere: positive magnitudes in
  MPa** (phylloptim #25). `scientific_version` for TF24 goes 5 -> 6, TF24f 5.1 -> 6.1.
  Migration:
  * `leaf$root_collar_psi_` -> `leaf$opt_root_psi_`, **and its sign flips**: it is
    now the positive magnitude. Renamed rather than reused deliberately, so an old
    script gets an error instead of quietly reading the wrong sign.
  * the **`opt_root_psi` aux changes sign** for the same reason, and now agrees with
    TF24f's `opt_root_psi_state`, which always held the magnitude. Before, the aux
    reported the signed potential while the state was negated back from it — an
    inconsistency in plant's own outputs. Any stored output or plot that negated the
    aux must stop doing so.
  * `l$E_from_Soil_to_Root_Collar(collar, psi_soil)`, `l$find_root_psi(...)`,
    `l$find_psi_stem_from_psi_root(...)`, `l$dE_from_soil_dpsi_collar(...)` and
    `l$transpiration_to_psi_stem(E, psi_upstream)` keep their signatures but now
    take **positive magnitudes**, and the bracket ends of `find_root_psi` swap
    (wettest layer first, `psi_crit` second). Passing the old signed values raises
    an error rather than returning a wrong number — the leaf package validates it.
  * outputs move by about **-3.2e-4 relative** (one-species SCM offspring production
    83.9026 -> 83.8761). Not an equation change: the rewrite is exactly
    sign-symmetric, but boost's TOMS748 iterates depend on the bracket's
    orientation, and the SCM amplifies the resulting 1-3 ULP. Every pinned test
    value and the exact stochastic TF24 counts pass unchanged.

* The TF24 leaf gas-exchange and hydraulics model now comes from the standalone
  header-only [`phylloptim`](https://github.com/traitecoevo/phylloptim) package instead of
  a copy in this repo (phylloptim #9). `plant::Leaf` is an alias for `phylloptim::Leaf`, so
  C++ consumers that `#include <plant/leaf_model.h>` and read public `Leaf` members
  keep compiling unchanged. **`scientific_version` for TF24 goes 4 -> 5 and TF24f
  4.1 -> 5.1: TF24 output moves by about +2.4%** — see `tf24_strategy.h` for the
  measurement and its attribution. Migration:
  * `l$set_physiology(area_leaf =, mass_root_prop =, rho =, a_bio =, PPFD =,
    psi_soil =, soil_depth =, leaf_specific_conductance_max =, atm_vpd =, ca =,
    sapwood_volume_per_leaf_area =, leaf_temp =, atm_o2_kpa =, atm_kpa =)`
    -> `l$set_physiology(root_carbon_per_leaf_area =, PPFD =, psi_soil =,
    soil_depth =, leaf_specific_conductance_max =, atm_vpd =, ca =, leaf_temp =,
    atm_o2_kpa =, atm_kpa =)`. **This is a semantic change, not just a rename:**
    `area_leaf`, `rho`, `a_bio` and `sapwood_volume_per_leaf_area` were dead
    stores and are gone, and `root_carbon_per_leaf_area` is the old
    `mass_root_prop` **divided by `area_leaf`**. The leaf is purely intensive now;
    uptake is exactly homogeneous in that ratio, so passing the old absolute
    carbon silently gives a root system too weak by a factor of `1/area_leaf`.
  * `l$area_leaf_`, `l$rho_`, `l$a_bio_`, `l$sapwood_volume_per_leaf_area_`
    (removed) -> no equivalent; all four were assigned and never read.
  * `plant::umol_per_mol_to_Pa` (C++, removed) -> `Leaf::umol_per_mol_to_Pa_`, now
    derived per call from `atm_kpa` rather than hard-coded at 0.1013.
  * `l$initialize_integrator(rule, tol)` keeps its name but now sets the tolerance
    of the leaf package's header-only adaptive Simpson quadrature rather than
    configuring plant's compiled QAG; only `transpiration_full_integration`, a
    spline-fidelity diagnostic, was ever affected.
  * Requires **odelia at master** (commit `d8235d1` or later), not just the
    `>= 0.2.0` in DESCRIPTION: `Patch::ode_rates` is non-const since #585 and
    odelia 0.2.0's `r_ode_rates` takes the system by `const&`. plant's `develop`
    does not compile against 0.2.0 either — this is not new here, but it will bite
    anyone building from a released odelia.

* Penman-Monteith leaf energy balance added to TF24/TF24f behind an
  opt-in gate, **default off** (#523). With the gate off, runtime behaviour and
  outputs are bit-identical to before (verified against leaf-level and SCM
  baselines) and `scientific_version` is unchanged — but the change adds fields
  to the R-facing `TF24_Pars` and a new environment driver, so any snapshot /
  round-trip test that encodes the full TF24 parameter set or driver list will
  see new entries. Migration:
  * `TF24_Pars` gains two fields (default off / inert): `pars$use_energy_balance`
    (0 = off = today's `Tleaf = Tair`; non-zero = on) and `pars$d` (characteristic
    leaf dimension, m, for the aerodynamic resistance). No action needed unless
    you assert the exact `pars` field set.
  * `TF24_Environment` gains a `wind_speed` extrinsic driver (default `2.0`
    m s⁻¹), set like any other driver:
    `env$extrinsic_drivers_set_constant("wind_speed", U0)` /
    `env$extrinsic_drivers_set_variable("wind_speed", x =, y =)`.
  * The `Leaf` submodel exposes `use_energy_balance_`, `d_`, `wind_speed_`
    (settable) and `Tair_`/`Rn_`/`ra_` (readable) for leaf-level experimentation.
  * To enable: set `pars$use_energy_balance <- 1` (and optionally `pars$d`,
    the `wind_speed` driver) before running; leaf temperature is then solved from
    the energy balance and fed to the Farquhar temperature scaling.

* Strategy biological parameters are now stored in a nested `pars` sub-object
  rather than as flat fields on the strategy (#410). This applies to every
  strategy type (`FF16_Strategy`, `K93_Strategy`, `TF24_Strategy`). Migration —
  read/write any biological parameter through `$pars`:
  * `s$lma`    -> `s$pars$lma`
  * `s$rho`    -> `s$pars$rho`
  * `s$hmat`   -> `s$pars$hmat`
  * `s$eta`    -> `s$pars$eta`
  * `s$k_I`    -> `s$pars$k_I`
  * `s$S_D`    -> `s$pars$S_D`
  * …and likewise for every other strategy parameter (all FF16/TF24 allometry,
    production, mortality and hydraulic parameters; all K93 `b_*`/`c_*`/`d_*`
    and `height_0`). Construction by flat field also moves under `pars`:
    `FF16_Strategy(lma = 0.1)` -> `s <- FF16_Strategy(); s$pars$lma <- 0.1`
    (or build from traits via `add_strategies()`).
  * Top-level strategy fields are unchanged: `control`, `birth_rate_x`,
    `birth_rate_y`, `is_variable_birth_rate`, `collect_all_auxiliary`.
  * A `<Strategy>_Pars` constructor is exported per strategy (`FF16_Pars()`,
    `K93_Pars()`, `TF24_Pars()`).
* The trait -> strategy interface was renamed to be clearer and pipe-friendly,
  with the `Parameters` object now the first argument (#410). The old names
  remain as deprecated shims for one release. Migration:
  * `expand_parameters(traits, p, hyperpar, birth_rate_list =, keep_existing_strategies =)`
    -> `add_strategies(p, traits, hyperpar =, birth_rate =, keep_existing =)`
  * `mutant_parameters(traits, p, …)` -> `add_mutant(p, traits, …)`
  * `strategy_list(traits, p, hyperpar, birth_rate_list)`
    -> `generate_strategy(p, traits, hyperpar =, birth_rate =)`
  * argument `birth_rate_list` -> `birth_rate`; `keep_existing_strategies` ->
    `keep_existing`.

* SCM cohort-refinement now happens entirely in C++, and the R interface is
  consolidated onto a single `run_scm()`. `run_scm_collect()`, `run_scm_error()`
  and `build_schedule()` have been removed (#408, #459, #462). Migration:
  * `build_schedule(p, …)` -> `scm <- run_scm(p, …, refine_schedule = TRUE)`,
    then read the refined parameters as `scm$parameters`. The former
    `attr(p_new, "offspring_production")` side-channel is now
    `scm$offspring_production`.
  * `run_scm_collect(p, …)` -> `run_scm(p, …, collect = TRUE)`
  * `run_scm_error(p, …)` -> set `scm$collect_errors <- TRUE` then read
    `scm$combined_node_errors`
* The argument order of `run_scm()` changed; `use_ode_times` is now the last
  argument (after the new `refine_schedule` and `collect`).
* Numeric `Control` defaults are now set in the C++ `Control()` constructor
  (the pragmatic "fast" settings), and the two R preset helpers were removed
  (#463). Raw `Control()` now means **fast**, not accurate. Migration:
  * `fast_control()` -> `Control()` (or its lowercase alias `control()`)
  * `scm_base_control()` -> `Control()` / `control()`
  * for the old tight-tolerance behaviour -> `control_accurate()`
  * `scm_base_parameters()` is unaffected (it builds a `Parameters`, not a
    `Control`).
* Removed the TF24/TF24f parameter `hk_s` from the public interface and from
  the `Leaf` constructor; hydraulic cost now depends on `g1_TF24`, `beta2`, and
  the vulnerability curve terms only. Migration:
  * `TF24_Strategy()$pars$hk_s` -> removed (no equivalent)
  * `TF24f_Strategy()$pars$hk_s` -> removed (no equivalent)
  * `make_TF24_hyperpar(..., B_hks1 =, B_hks2 =)` -> adjust `g1_TF24` ~ `rho` scaling (new; no `hk_s` equivalent)
  * `make_TF24f_hyperpar(..., B_hks1 =, B_hks2 =)` -> adjust `g1_TF24` ~ `rho` scaling (new; no `hk_s` equivalent)
  * `Leaf(..., hk_s =)` -> remove `hk_s` argument
* The ODE solver and interpolator core were spun out into the standalone
  [odelia](https://github.com/traitecoevo/odelia) package, which plant now
  depends on (#456, #464). Plant's internal solver/interpolator headers
  (`inst/include/plant/ode_solver/*`, `interpolator.*`, `tk_spline.*`) were
  removed in favour of `odelia::ode::Solver` / `odelia/interpolator.hpp`.
  Migration: code that linked against plant's C++ numerics headers, or used the
  old `OdeRunner`/interpolator R6 types directly, must switch to odelia. Pure-R
  callers of `run_scm()` / `run_stochastic()` are unaffected.

The following earlier changes (since the 2.0.0 release) are also breaking and
were not previously recorded here:

* All fitness/equilibrium functionality was removed from plant; it now lives in
  the separate `regnans` package (#388). Removed `fitness_landscape()`,
  `solve_max_fitness()`, `viable_fitness()`, `fundamental_fitness()`,
  `assembly_parameters()`, `equilibrium_birth_rate()` and the `equilibrium_*`
  parameters/controls.
* The `Environment` object was removed from `Parameters`; pass it directly as
  the `env` argument to run functions, e.g.
  `run_scm(p, env = Environment("FF16"))` (#315). Likewise `Control` was removed
  from `Parameters` and is passed as the `ctrl` argument (#314).
* SCM & environment interface simplified (#446): environment construction is now
  the single `Environment("FF16")` constructor — removed `make_environment()`,
  `FF16_make_environment()`, `FF16_fixed_environment()`, `test_environment()`,
  and the helpers `scm_patch()`, `make_scm_integrate()`, `scm_state()`,
  `patch_to_internals()`, `scm_to_internals()`, `species_to_internals()`,
  `first()`, `last()`, `modify_list()`. State collection moved to C++ and
  `tidy_patch` runs by default.
* Node/Species/SCM competition methods renamed (#448):
  * `Species$competition_effects` -> `Species$compute_competition_effect_by_nodes`
  * `Species$competition_effects_error` -> `Species$compute_competition_effect_by_nodes_error`
  * `SCM$competition_effect_error` -> `SCM$compute_competition_effect_error_by_node_for_species_i`
  * the `Node$competition_effect` getter was removed.
* The exported-function surface was reduced (#420): the per-model R6 generators
  (`FF16_*`, `K93_*`, `TF24_*` `Node`/`Patch`/`SCM`/`Species`/`Stochastic*`) and
  many utilities (`Internals`, `OdeControl`, `bounds`, `check_bounds`,
  `clamp_domain`, `nlsolve`, `splinefun_log`/`splinefun_loglog`, `strategy`,
  `strategy_default`, `validate`, …) are no longer exported. Rename:
  `make_transparent()` -> `util_colour_set_opacity()`.
* The `Cohort` class/concept was renamed `Node` throughout — classes, methods,
  filenames, and R-side `coh` variables (#335).
* "Seed rain" terminology renamed for clarity (#297):
  * `seed_rain[_in]` -> `birth_rate`
  * `seed_rain[_out]` -> `net_reproduction_ratio` / `offspring_production`
  * `add_seeds()` -> `introduce_new_cohort()`
* Canopy/light renamed and generalised (#384): the `canopy` class became
  `Resource_spline` (`canopy.h` -> `resource_spline.h`); `compute_canopy` ->
  `light_availability$compute_environment`; `canopy_light_*` controls ->
  `light_availability_spline_*`; `get_canopy_at_height()` / `canopy_openness()`
  -> `get_value_at_height()`. The standalone `assimilation` class and its
  adaptive-integration options were removed.
* Extrinsic-driver API reworked (#340): per-driver `Environment` methods replaced
  by an `ExtrinsicDrivers` object at `env$extrinsic_drivers`, with
  `evaluate("rainfall", t)` / `evaluate_range("rainfall", c(...))`.
* ODE-stepper methods renamed (#413): `advance()` -> `advance_adaptive()`;
  `node_schedule_ode_times()` -> `ode_times()`. The corresponding `Parameters`
  field was also renamed: `p$node_schedule_ode_times` -> `p$ode_times`.
* Growth-optimisation routine generalised (#382):
  `FF16_solve_max_size_growth_rate_at_height()` ->
  `optimise_individual_rate_at_size_by_trait()` (new `size` / `size_name` args),
  with a height wrapper `optimise_individual_rate_at_height_by_trait()`.
* Removed the `FF16w` strategy, superseded by TF24 (`FF16w` -> `TF24`) (#438),
  and the `FF16r` strategy (#439); removed the experimental soil component from
  `FF16` (#441).
* Replaced the Leaf cost parameter `g1_TF24` with the Medlyn stomatal-conductance
  parameters `g0` (default 0.022) and `g1` (default 2.57) on the TF24 `Leaf`
  (#450, #451).
* `Disturbance` refactored into a `Disturbance_Regime` base with
  `No_Disturbance_Regime` / `Weibull_Disturbance_Regime` subclasses (#301);
  disturbance configuration moved from `Environment` to `Patch` (#290); the
  light-extinction coefficient `k_I` moved from `Parameters` onto the strategy
  (#293, #302).
* The stochastic runner's node-schedule interface was renamed to mirror `SCM`
  (#217, #506). Migration:
  * `StochasticPatchRunner$schedule`              -> `StochasticPatchRunner$node_schedule` (get and `<-` assignment)
  * `StochasticPatchRunner$set_schedule_times(x)` -> `StochasticPatchRunner$set_node_schedule_times(x)`
* `run_stochastic_collect()` no longer returns an `offspring_production` element
  in its result list — it is not defined for the finite-population model and was
  never populated (#498, #506). Migration:
  * `run_stochastic_collect(...)$offspring_production` -> removed (no equivalent)
* `Patch` computes its rates when they are read rather than when its state is
  set, so reading `$ode_rates` now evaluates and the separate compute entry point
  is gone (#585). Migration:
  * `patch$compute_rates()` (removed) -> `patch$ode_rates`, which computes at the
    state currently loaded and returns the result
  * Semantic change: `patch$ode_rates` was a cheap read of stored values and is
    now a full right-hand-side evaluation. `StochasticPatch$compute_rates()` is
    unchanged.

### New features

* **`Control$node_density_in_birth_date`** (default `FALSE`) carries the SCM's
  size distribution as a density in birth date instead of in height.

  The transport equation's compression term is the total derivative of the growth
  rate along a cohort's own trajectory, which equals `∂g/∂h` only when growth is
  a function of size. TF24's reserve gate breaks that: the finite-difference
  probe in `Node::growth_rate_gradient` moves height while holding *absolute*
  carbon fixed, so it shifts the reserve fraction `r = S/S_max`, whereas a cohort
  actually grows with `r` roughly constant. The probe is accurate about a
  quantity the plant never experiences, so this is a different derivative rather
  than a worse approximation of the right one.

  In birth-date coordinates the density rate is mortality alone (nothing moves an
  individual along the birth-date axis), the birth density is
  `birth_rate·pr_estab` with no division by the growth rate, and both resource
  integrals run over introduction times.

  With the flag off the solver is bit-identical to before. With it on:

  - FF16 and K93 have size-only growth, so both coordinates must converge to the
    same answer, and do. Relative difference in offspring production over
    successive halvings of the node spacing: FF16 `1.0e-2`, `4.2e-3`, `1.2e-3`,
    `2.9e-4`; K93 `2.6e-3`, `5.9e-4`, `1.4e-4`, `3.5e-5`. That is the ~2nd order
    of the trapezium rule, and confirms the coordinate change is a change of
    coordinates.
  - The gap is almost entirely the *height* coordinate's error. Across those same
    four schedules FF16's birth-date answer moves 18.6956 → 18.6996 while its
    height answer climbs 18.5071 → 18.6941 toward it. At the default schedule the
    birth-date answer is already within `2e-4` of converged and the height answer
    is `1.0e-2` away — roughly 50x more accurate for the same number of nodes.
  - TF24, whose growth is not a function of size alone, is the case the two
    coordinates genuinely disagree on, and the diagnostic is that refining the
    schedule does *not* close the gap the way it does for FF16 and K93. At
    `lma = 0.0825`, `hmat = 5`, `birth_rate = 20`,
    `max_patch_lifetime = 30`, over the same three schedules: height
    474 → 587 → 696, birth date 3495 → 3729 → 3840, a ratio of 7.4 → 6.3 → 5.5.
    Both are still moving, so neither figure is a converged value; the point is
    that they are not converging *to each other*, which is what a wrong
    compression term looks like as against a coarse quadrature.

  Multi-species runs behave the same way. The established two-species FF16 and
  three-species K93 cases converge per species at the same ~2nd order (FF16
  `4.4e-3`/`7.1e-3` → `2.8e-4`/`6.0e-4`; K93 all three species `3.1e-3`–`6.2e-3`
  → `1.9e-4`–`3.7e-4`), and K93's competitively marginal first species — ~90x
  below the dominant one — converges no worse in relative terms than the
  dominant ones. So the shared competition profile does not disadvantage a
  marginal species under the coordinate change.

  **TF24 two-species changes the ecological outcome, not just the numbers.** On
  the two-species case of `test-strategy-tf24.R` (`lma` 0.0825 / 0.10,
  `max_patch_lifetime = 30`), the height coordinate excludes the slower species
  (offspring production `1.4e-4` against `503` for the faster) while the
  birth-date coordinate has them coexisting at comparable abundance (`1004`
  against `2707`), and refining the schedule does not move either toward the
  other. The conclusion recorded in that test — that reserve-gated growth
  (#517) largely excludes the slower species — is therefore coordinate
  dependent, and needs re-deriving before it is relied on.

  `Node::growth_rate_gradient` is not called, which also removes one leaf solve
  per cohort per Runge-Kutta stage.

  **What R receives is unchanged in meaning.** The coordinate is an internal
  choice, so `log_density` is converted back to a density in *height* before it
  leaves C++ — `Species$log_densities`, `Patch$state`, and therefore
  `run_scm(collect = TRUE)`, `tidy_outputs.R`'s `density = exp(log_density)`,
  `interpolate_to_heights()` and the plots all keep their existing meaning. The
  conversion is `N = ν / |dh/dτ|`, with the Jacobian formed by central
  differences of adjacent node heights against their birth dates (the boundary
  node supplies the extra point at the young end). Measured against an otherwise
  identical height-coordinate run, the median node agrees to `1.7e-3` in log
  space — 0.17% in density — improving to `4.4e-4` and `1.1e-4` over two
  schedule refinements. See Known issues for where this is weakest.

  The boundary node needs no differencing at all: `dh/dτ = −g(H₀)` at birth, so
  `compute_initial_conditions()` now records the birth growth rate on every node
  (`Node$growth_rate_at_birth`) and the Jacobian uses it there directly. It is
  exact *and* current for that node, because the boundary node is re-evaluated
  every step; for an introduced node the recorded rate is frozen at its own
  birth and so cannot serve as its present-day Jacobian. This removes a 31%
  error on the boundary node — it sits a whole introduction interval from its
  neighbour, which is the worst case for a one-sided difference. Recording it
  also lets the two coordinates' boundary conditions be checked against each
  other: `exp(log_density) · g(H₀)` on the height path must equal
  `birth_rate · pr_estab`, which the birth-date path carries directly, and that
  is now a test.

  The quantity actually integrated is reported alongside rather than lost:
  `Species$log_densities_state` (and a `log_density_state` row in `Patch$state`,
  present only on the birth-date path), plus `Species$height_jacobian` for the
  conversion itself. `Node$log_density` stays unconverted, since forming
  `|dh/dτ|` needs the neighbouring nodes and a `Node` does not have them.
  `export_patch_state()` resumes from `patch$ode_state`, which is the raw state,
  so resume is unaffected by any of this — pinned by a test.

  Two invariants the birth-date axis needs, which the height axis did not:

  - The boundary node's birth date is the current time, so
    `Patch::compute_environment()` refreshes it before building the profile.
    `compute_rates()` (which stamps it) runs *after* the `set_ode_state()` that
    rebuilds the environment, so reading the stamp there would use the previous
    derivs call's time. The measured effect on FF16 offspring production is below
    `1e-6` — the boundary node carries almost no leaf area — but it made the
    spatial quadrature a function of the ODE step size.
  - Introduction times are the quadrature grid, so repeated ones span zero width
    and drop out of the integrals. A scheduled run cannot produce them; a patch
    seeded or resumed *without* per-node times gives every node the boundary
    node's birth date, which would silently zero the competition profile.
    `Patch::check_birth_dates_distinct()` now rejects that with an actionable
    message. A faithful `export_patch_state()` / `set_initial_state()`
    round-trip carries the times and is unaffected.

  `Species::compute_competition()` keeps the sorted-grid fallback of #574 for the
  height coordinate only: that path integrates in height, so sending the
  birth-date coordinate down it would swap coordinates mid-run. The birth-date
  abscissa cannot invert (introduction times are fixed at birth), but the
  *height* early exit is skipped when the height ordering has broken, since a
  node below the query height can then be followed by a taller one.
* **A dry TF24f patch no longer aborts the whole run on the ci root-find.**
  `Leaf::dprofit_droot_collar_psi` — TF24f's exact AD/IFT gradient — called
  `psi_stem_to_ci()` before testing for hydraulic shut-down. In shut-down,
  `psi_upstream >= psi_stem` (both positive magnitudes), so
  `gc = const · transpiration` goes negative, the residual stops crossing zero
  over (Γ*, ca], and TOMS748 throws *"a and b do not bracket the root"* rather
  than returning non-finite — which is what the existing `isfinite` guard on the
  next line was written to catch, and cannot see. One dry patch therefore killed
  the run. The condition is now tested *before* the solve, matching what
  `set_leaf_states_rates_from_psi_stem()` has always done, and returns a zero
  gradient. Reproduced at 5 soil layers, θ = 0.005–0.03 with 1 m yr⁻¹ rainfall
  (ψ_stem = 1.23 against ψ_upstream = 5.92 MPa). Behaviour changes only in states
  that previously threw, so no scientific version bump.
* **Adaptive interpolation now names a non-finite target instead of blaming
  resolution.** `AdaptiveInterpolator::check_err()` compares against NaN, and
  every NaN comparison is false, so a single NaN or Inf from the target made its
  interval permanently unacceptable: refinement halved the spacing to
  `max_depth` and then reported *"Interpolated function as refined as currently
  possible"*. That message cost real debugging time on a TF24 patch. A
  non-finite value now fails immediately, naming the point and the value; the
  resolution-limit message additionally reports the spacing reached, the limit,
  and the tolerances missed. A new `test_adaptive_interpolator()` test hook
  drives the refiner from R (the production caller passes a C++ lambda).
* **The scenario scorecard now reports `persists`.** Whether a strategy replaces
  itself, R0 >= 1, which is a different question from whether the run completed.
  It matters because `status`/`outcome` test `total > 0`: at the current baseline
  five of the eight hydraulic scenarios return R0 between 2e-15 and 6e-14 —
  numerically extinct — and were all recorded as `"persisted"`. Only **1 of 8**
  scenarios persists. That is the main reason the gateway discriminates so little
  (3/8, with no expected failure failing). Reported alongside, not folded into,
  the existing classification: the scenario CSV's "Model failure" means the model
  *breaks numerically* (#549/#550), not that the strategy dies out, and the
  blessed baseline diff is defined on the existing columns.
  `scenario_summary()` gains `n_persists` and tolerates scorecards recorded
  before the column existed.
* **A hydraulically shut-down plant no longer draws water from the soil
  (`TF24@v4`, `TF24f@v4.1`).** Both shut-down exits in
  `Leaf::find_root_collar_psi` set `profit_` directly and bypass
  `profit_psi_stem_TF`, and the first returns before any
  `E_from_Soil_to_Root_Collar` call in that solve. Because `Leaf` is a value
  member reused across every `compute_rates` call for an individual, every leaf
  output they did not assign kept the **previous step's** value. That was not
  just a reporting problem: `soil_consumption_` feeds
  `TF24_Strategy::evapotranspiration_dt` and hence the patch water balance, so a
  plant whose stomata had closed carried on extracting its last wet-step uptake.
  Measured: an individual moved from θ = 0.30 to θ = 0.02 (past ψ_crit) reported
  `E_up_` and `transpiration` *identical* to its wet step. Both exits now zero
  the transport chain (`transpiration_`, `stom_cond_CO2_`, `E_up_`,
  `soil_consumption_`), and recovery on rewetting is tested. Note the water
  budget still *closed* in the buggy state — what was recorded as depleted was
  what was removed — so the conservation tests could not catch this; only the
  physics was wrong. Because water-limited runs change (scenario gateway:
  offspring production moves by up to 5e-3 relative on 5 of 8 scenarios, with
  every success/failure classification unchanged), the TF24 scientific version
  is bumped **3 → 4** and `TF24f` tracks to **4.1**. This invalidates `logpile`
  caches for water-limited TF24 runs, which is the intended, safe direction.
* **Root hydraulic parameters are now settable from R.** `root_b`, `root_c`,
  `root_psi_crit` and `rooting_depth_max` move into `TF24_Pars` (they were fixed
  members of `TF24_Strategy` and a file-static constant in
  `src/tf24_strategy.cpp`, so unreachable from R). Root shutoff was pinned at
  ψ ≈ 5.87 MPa, too conservative for taxa that operate below it (e.g. *Acacia
  aneura*), and rooting depth at 1.5 m. **Defaults are unchanged**, so this is
  an interface change only and the TF24 scientific version does *not* move.
  Two cautions: `root_psi_crit` is derived from `root_b` and `root_c` exactly as
  `psi_crit` is from `b` and `c`, so set it whenever you set those; and
  `rooting_depth_max` beyond the soil column depth (`TF24_Environment$depth`,
  default 1.5 m) gains nothing, because the layers do not exist.
* **`check_driver_interpolation()`** reports how a driver's control points
  survive interpolation — evaluated range, fraction of the series that goes
  negative, undershoot area, and the interpolated vs supplied integral — and
  warns when a non-negative series interpolates negative. Extrinsic drivers use
  a cubic spline, which is a poor fit for intermittent forcing: a realistic
  daily rainfall series with a ~10% wet-day fraction evaluates negative at ~45%
  of points, reaching −5.7 m yr⁻¹. Worth running on any site-forcing series
  before committing to a long run.
* **`assimilation` is now reported for TF24/TF24f.** The auxiliary variable was
  listed in `aux_names()` but never written — no index was resolved and no
  `set_aux` call referenced it — so the slot reported whatever `Internals`
  happened to hold, and carbon uptake was unavailable as a model output. It now
  carries net CO₂ assimilation at the optimal operating point, per unit leaf
  area (µmol CO₂ m⁻² s⁻¹). **Net, not gross:** `Leaf::assim_colimited()`
  subtracts dark respiration, so gross = `assimilation` + R_d with
  R_d = 0.015·vcmax at the acclimated vcmax. It is integrated over the crown
  under `deep-crown` shading alongside the other leaf outputs, and the two
  hydraulic shut-down exits now set it explicitly (to −R_d) instead of leaving a
  stale probe value, keeping `profit == assimilation − hydraulic cost` true in
  every branch.

* **NSC storage pool for TF24 (`TF24@v3`, `TF24f@v3.1`).** TF24 now carries a
  non-structural-carbohydrate storage state so growth and mortality respond to
  *buffered* carbon rather than instantaneous net production (#517, #554).
  Growth/reproduction are reserve-gated (a smooth logistic on relative reserves
  `r = S/S_max`), and mortality is now `a_dG1·exp(-a_dG2·r)` — bounded in
  `[a_dG1·e⁻ᵃ_dG2, a_dG1]` — which is the root-cause fix for the SCM
  cohort-density blow-up (#550). New parameters `a_st1`/`a_st2`/`a_st3` (storage
  capacity per unit sapwood, growth half-on reserve fraction, birth fill).
  Because this changes the simulation output for identical inputs, the TF24
  scientific version is bumped **2 → 3**; `TF24f` inherits the state and its
  compound version auto-tracks to **3.1**. This invalidates the `logpile` cache
  for both models. See the staging guide in `overstorey-staging/guides/`.
* **Per-model scientific versioning.** Each model now carries a scientific
  version — an integer that is independent of the package `Version` and is
  bumped only when the model's equations or default parameters change the
  simulation output for identical inputs (not for refactors, performance, or
  interface changes). Read it with `model_version("FF16")` (an integer) or
  `model_id("FF16")` (`"FF16@v1"`). The number is authored as the
  `scientific_version` constant on each strategy class in
  `inst/include/plant/models/*_strategy.h`, next to the equations it versions.
  Downstream tools (e.g. `logpile`) use it to decide when archived simulations
  must be re-run: reruns follow scientific changes, not every software release.
  A drift-guard test (`tests/testthat/test-model-version.R`) fails when a
  model's default parameters change without a version bump. Starting versions:
  `FF16@v1`, `K93@v1`, `TF24@v2` (a published result used the pre-versioning
  "v1" science). `TF24f`, being a fast *approximation* of TF24, carries a
  compound version `"<TF24 version>.<approximation revision>"` (`TF24f@v2.1`):
  the major component auto-tracks TF24 so a TF24 change also invalidates TF24f,
  and the minor tracks changes to the approximation itself.
* `run_scm()` can start a patch from **pre-existing nodes** rather than always
  growing from empty — to resume an exported run or to seed an arbitrary initial
  size distribution at patch age 0 (#499, revives #304). New
  `export_patch_state(scm, step)` captures a patch's full state; the initial
  condition rides on the `Parameters` object (`initial_state`,
  `n_initial_cohorts`, …) so `run_scm(p)` just works and stays serialisable.
* Selectable **canopy shading models** for FF16/TF24 via `control$shading_model`
  (#490, #417): `""` (each strategy's own default), `deep-crown`, `mean-light`,
  `crown-centre`, `flat-top-soft-box`, `flat-top-box`, `ppa` (with
  `ppa_layer_optical_depth` / `ppa_layer_smoothing`). Dispatch is bound once in
  `prepare_strategy()`; deep-crown is bit-for-bit unchanged.
* Alternative **fixed-step forward-Euler** ODE integration alongside the
  adaptive Cash-Karp solver, selected with `control$fixed_time_step` (years;
  `0` = adaptive RKCK, the default) (#489). Combining `fixed_time_step > 0`
  with a mutant run / `save_RK45_cache` / `use_ode_times` errors clearly.
* Multi-layer **root water uptake** for the TF24 strategy (soil→root→stem→leaf
  hydraulic pathway following Potkay et al. 2021), replacing the previous
  single-value soil-water treatment (#488).
* Added the **TF24 strategy** — a leaf-level water-use/hydraulics model with a
  `TF24_Environment` whose soil-layer count and parameters can be set after
  construction (#445), including a Medlyn stomatal-conductance model on the leaf
  (`solve_medlyn_ci_*`, `medlyn_model_gs`) (#450).
* New **mutant-fitness method**: caches resident environments at each ODE step so
  mutant fitness reuses the residents' resource shadow. Adds a `save_history` /
  `save_RK45_cache` control, `environment_history` / `patch_step_history` on the
  SCM, and a `mutant_parameters()` method (#362, #379).
* Per-species **time-varying birth rates** via the extrinsic-driver mechanism:
  `set_constant_birth_rate(p, i, k)`, `set_interpolated_birth_rate(p, i, x, y)`,
  and `Parameters` fields `birth_rate_x` / `birth_rate_y` /
  `is_constant_birth_rate` (#334).
* `run_scm()` collected output now exports the extrinsic environment drivers per
  timestep (`out$env[[t]]$canopy`, `…$rainfall`); new drivers export
  automatically (#347).
* `expand_state()` generalised to apply to any strategy (#443).
* New strategy parameter `recruitment_decay` — exponential decline of
  establishment probability with patch age (#330).
* Per-individual `consumption_rates` (water extraction across cohorts/soil
  layers) feeding the water-enabled environment's soil balance (#329); FF16
  rainfall getters/setters + spline exposed to R (#324); environment auxiliary
  variables (#323) and environment state variables (#305) exposed/added.
* Added an HTML report (plots + analyses) for the FF16 strategy (#350).

### Known issues

* On the birth-date path, the *worst-case* node of the reported height density
  does not improve with schedule refinement, even though the typical node does
  (see the reporting note under New features). `|dh/dτ|` is a ratio of two
  differences that both shrink as the schedule is refined, so cancellation error
  sets a floor, and refinement adds nodes in the near-empty tail where that floor
  is worst. Aggregate and plotting use is sound; individual node densities far
  out in the tail are not. Note the exact `−g(H₀)` Jacobian does not help here:
  it is exact only at a node's own birth, and the worst cases are interior nodes
  in a compressed region of the size distribution (log density ~ −11), not the
  boundary. Removing this properly means evolving `∂h/∂τ` along the
  characteristic, whose rate is `∂g/∂h · ∂h/∂τ` — which is the compression term
  this change exists to avoid, so it would have to be an opt-in extra state
  wanted only for reporting.

* `stochastic_schedule()` passes `patch_area` into
  `stochastic_arrival_times()`'s third positional argument, which is `delta_t`.
  Arrival rates therefore never scale with patch area, and the binning interval
  is set to the area instead. Passing it by name is the fix, but
  `test-stochastic-patch-runner.R` runs at `patch_area = 50` with seed-pinned
  expectations, so correcting it makes that file ~50x heavier and needs its
  parameters and baselines revisited. Probes that need a correct schedule build
  their own (see `plant-dev` `probes/12-oracle.R`).

### Minor changes & bug fixes

* `SpeciesBase::control()` called `strategy->get_control()`, which does not
  exist on any strategy. The member had never been instantiated, so the error
  had never been compiled; it now reads `strategy->control`.
* **The TF24 rainfall driver is floored at zero.** Because drivers are
  interpolated with a cubic spline, an intermittent series undershoots below
  every supplied value, and negative rainfall gave negative infiltration and an
  unphysical drying rate. That failed two different ways depending on soil
  wetness: above residual moisture the water really was removed, while at or
  below residual the guard in `compute_rates` clamped the rate, so the removal
  was recorded in `sum_rainfall` but never applied and the water budget stopped
  closing. Drylands sit at residual for much of the year, so the second case is
  the common one. Note the floor bounds the *sign* only and is not a correction
  to the interpolation: the spline conserves the integral exactly (undershoot is
  compensated by overshoot), so discarding the negative lobes raises total
  rainfall by the undershoot area — ~7% for a realistic daily series. The remedy
  is not to spline an intermittent series; see `check_driver_interpolation()`.
  Runs with constant or smooth seasonal rainfall are unaffected (the scenario
  gateway's sinusoidal driver never goes negative), so the TF24 scientific
  version does not move.
* Corrected the comment on `Leaf::assim_colimited()`, which claimed "no dark
  respiration included at the moment" while the code subtracts `R_d_`. The
  function returns a *net* rate; the comment now says so, along with how to
  recover the gross rate.
* Water-budget test coverage extended from a single layer to 1, 5 and 15 layers,
  now including root uptake in the balance (the previous check ran for 0.01 yr,
  where uptake was negligible), across the saturated-to-dry range and at and
  below the residual-moisture floor. Closure holds to round-off (relative
  residual < 1e-12) in every case. A companion test documents *why* the residual
  clamp never leaks in practice: with `n_psi = 6.57` the conductivity exponent
  is 2·n_psi+3 ≈ 16, so K(θ) collapses and ψ(θ) diverges far above θ_r —
  transport has already stopped, making the clamp a safety net rather than an
  active mass sink.

* `run_scm()` now fails with an actionable message when the SCM size-density
  (characteristic) equations run away under extreme forcing (e.g. severe
  seasonal drought in TF24): a cohort density overflowing to `+Inf`, or the
  density-weighted resource uptake driving a soil-water state non-finite. The
  guard fires each ODE step, before the non-finite value propagates into the
  competition integral or physiology, replacing the previous opaque
  `Detected non-finite contribution` / `non-finite psi_soil` errors (#550).
* `Node` now records its introduction time and patch-age density at the moment
  it is introduced, so reproduction and integration-error calculations no longer
  re-derive these from the node schedule / disturbance regime after the run.
* `SCM::run()` accumulates the per-node refinement error when `collect_errors`
  is set, exposed as `combined_node_errors`; `SCM::refine_schedule()` runs the
  adaptive node-introduction loop in C++.
* The resource spline is now floored at 0, fixing spurious negative light
  values (notably the K93 light spline at high `k_I`) (#253, #497).
* `interpolate_to_heights()` no longer silently drops individuals in the largest
  size class on a coarse grid; it appends the actual largest individual per time
  step instead (#352, #497).
* Fixed `tidyselect` deprecation warnings from `integrate_over_size_distribution()`
  by using string column names in `across()`/`rename()` selection contexts
  (#375, #501).
* Account for variable patch size in light/competition calculations (#422).
* Fixed the argument order in `net_mass_production_dt()` (#389).
* Fixed integration of density (#345).
* Fixed `scm_support` to use the correct error name `net_reproduction_ratio_errors`
  and guard the zero-offspring case (#447).
* The ODE solver now continues when the minimum step size is reached (#413);
  refined the `interpolate_to_heights()` error check (#437).
* `run_stochastic_collect()` returned empty output: it read a `state` accessor
  that `StochasticPatchRunner` no longer exposed, so every step silently
  collected `NULL` (#498, #506). It now reads `StochasticPatch$state` (a new
  accessor mirroring `Patch$state`), and each species' per-step state matrix
  carries an `is_alive` attribute, so the height/survival trajectory is
  populated again.

### Internals & performance

* The scenario gateway's seasonal rainfall driver places its spline knots **per
  year** (48 yr⁻¹) rather than per run, so the realised seasonality no longer
  depends on `max_patch_lifetime`. The previous `max(200, mpl × 6)` gave 6 knots
  per annual cycle at the default `mpl = 100`. Measured, that was more accurate
  than it looked — cubic interpolation of a sine at 6 points/cycle is accurate to
  5.4e-3 on a peak of 3.0 (0.18%) and conserves the annual total to 1e-5 % — so
  this is robustness, not a bug fix; the error is 4th-order in the spacing and
  would degrade at larger `mpl`. Gateway effect: offspring production moves by up
  to 6.2% relative (demography amplifies the 0.18% driver change), with every
  success/failure classification unchanged, so the baseline needs no re-blessing.
* Each strategy now keeps its biological parameters in a value-member struct
  (`FF16_Pars`/`K93_Pars`/`TF24_Pars`) exposed to R as a nested `pars` list;
  derived/computed members (eta_c, height_0, the TF24 Leaf model, solver
  tolerances) stay plain. Hot-path access stays inlined and FF16 output is
  bit-identical (#410).
* `FF16_expand_state`/`TF24_expand_state` no longer duplicate the C++ allometry
  formulas; they call the strategy's own functions via
  `FF16/TF24_strategy_expand_allometry` (`src/strategy_expand.cpp`). Output
  columns are unchanged; `K93_expand_state` remains a no-op. Fixed a latent bug
  in the (previously dead) `mass_above_ground` (it summed root and omitted
  heartwood); it now returns leaf+bark+sapwood+heartwood, matching the rate
  version and the expand_state output (#254).
* `Individual` no longer knows any strategy's state/aux layout: it passes the
  whole `Internals` to `strategy->compute_competition(z, vars)` /
  `net_mass_production_dt(env, vars)`, and each strategy reads its own slots
  (#266, #254). A new strategy with a different layout can reuse the generic
  `Individual` unchanged.
* Hot-path optimisation of FF16/TF24 (~2.7–3.5× FF16, ~1.3× TF24) (#471), K93
  SCM perf (#493), and shared `pow(area_leaf, a_l2)` in FF16 `compute_rates`
  (#361, #494). See the `profile-plant` skill for the methodology.
* `Species` and `StochasticSpecies` now share storage, the ODE state plumbing,
  and per-element serialisation through a compile-time `SpeciesBase`; the
  stochastic species stores a `StochasticNode` element (Individual + `alive`) in
  place of a parallel `is_alive` vector, and stochastic deaths are applied to the
  solver's live system in place (#217, #506). Internal only — FF16 is
  bit-identical and shows no timing change.
* Earlier internals (since 2.0.0): default compilation set to `-O2` (#365);
  assimilation decoupled into a per-strategy `Assimilator` (#313); extrinsic
  drivers refactored into an `ExtrinsicDrivers` class (#334, #340); test-only
  reference plant relocated into the tests (#418); unused logging removed (#436);
  Makefile dev-loop speedups (#319); general `R CMD check` cleanup (#440).

### Documentation & tooling

* Narrative docs (former `vignettes/`, theory, dated posts) migrated to the
  [Overstorey](https://traitecoevo.github.io/overstorey/) site; the pkgdown
  site is now the function reference only (#496).
* The strategy-scaffolder workflow and a profiling workflow are now captured
  as the `plant-new-strategy` and `profile-plant` skills (#495, #492).
* Relicensed from GPL-2 to **AGPL-3** (#457).
* Upgraded the minimum C++ standard from C++14 — currently C++20 (#442).
* Expanded the K93 (#421) and self-thinning (#369) vignettes; added a draft
  `extrinsic_drivers` vignette (#340).

## Plant 2.0.0 release notes

v2.0.0 was released on 25/02/2021

### Major changes

* Improved templating of strategies and environments to allow for inheritance and re-use.
  See: `FF16r_strategy` for example of method overloading and `K93_environment` for environment inheritance.
* Added two new models: Kohyama 1993 (K93) and a soil water strategy
* Recovered the FF16r strategy
* Decoupled patch, environments and strategies by moving several routines to strategies, e.g. assimilation
* Hyperparameterisation now handled in R only.
* Added a strategy implementation vignette

### Minor changes

* Moved strategy defaults to header
* Moved strategy and environment specific files to `inst/include/plant/models`
* Renamed several functions, including:
  * `germination` -> `establishment`
  * `Plant` -> `Individual`
  * `area_leaf` -> `competition_effect`
  * `area_leaf_above` -> `compute_competition`
  * `vars_phys` -> `rates`
  * `scm_vars` -> `compute rates`
* Updated vignettes, documentation and tests
* Switch Ubuntu versions in Travis pipeline
* Added `r_init_interpolators` methods for environment initialisation
* A few fixes to scaffolder and tests
* Increased lenience on integration test
* Fixed compiler warnings
* Removed PlantPlus

A full account of changes from the previous version is available on GitHub: [v1.2.1...v2.0.0](https://github.com/traitecoevo/plant/compare/v1.2.1...v2.0.0)

## Plant 1.2.1  Release Notes

v1.2.1 was released on 20/09/2019

### Major Changes

### Minor Changes

- update Makefile to use `pkgbuild` instead of `devtools` for building dll ( because of upstream changes )
- switching to using `remotes` instead of `devtools` for installing from github ( because of upstream changes )

A full account of changes from the previous version is available on Github: [v1.2.0...v1.2.1](https://github.com/traitecoevo/plant/compare/v1.2.0...v1.2.1)

## Plant 1.2.0  Release Notes

v1.2.0 was released on 20/03/2018

### Major Changes

- New strategy scaffolder ensures higher level operations work across all strategies and for new strategies. This includes `lcp_whole_plant` `XXX_PlantPlus`, `grow_plant_to_size`, `grow_plant_to_height`, `grow_plant_to_time`
- website now builds via `pkgdown` package (used to use staticdocs but this is deprecated)
- simplified workflow for building website
- converted supporting materials from tex to Rmd

### Minor Changes

- start documenting notes for developers in `inst/docs/developer_notes.Rmd`
- add `CITATION` file
- Address many issues in documentation and package setup causing rcmdcheck to fail

A full account of changes from the previous version is available on Github: [v1.1.0...v1.2.0](https://github.com/traitecoevo/plant/compare/v1.1.0...v1.2.0)

## Plant 1.1.0 Release Notes

v1.1.0 was released on 2/02/2018

### Major Changes

- Now compiles and runs on Windows machines (requires R 3.3.0 or newer)
- Further details on installation
- Enable assembly_parameters to accept more named arguments

### Minor Changes

- Add Appveyor for build tests on Windows machines
- Update tests to use latest version of testthat
- Remove package traitecoevo/callr, previously used to make system calls  
- Makefile: Add roxygen and RcppR6 to compile target
- roxygen & Rcpp updates
- Added a `NEWS.md` file to track changes to the package.

A full account of changes from the previous version is available on Github: [v1.0.0...v1.1.0](https://github.com/traitecoevo/plant/compare/v1.0.0...v1.1.0)

## Plant 1.0.0 Release Notes

v1.0.0 was released on 23/02/2016

This version corresponds to the paper describing the package:

Falster, DS, RG FitzJohn, Å Brännström, U Dieckmann, M Westoby (2016) plant: A package for modelling forest trait ecology and evolution. Methods in Ecology and Evolution 7: 136-146, doi: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)

A full account of changes from the previous version is available on Github: [v0.2.2...v1.0.0](https://github.com/traitecoevo/plant/compare/v0.2.2...v1.0.0)


## Plant 0.2.2 Release Notes

v0.2.2 was released on 1/06/2015

Draft paper about package submitted to Methods in Ecology & Evolution.

### Major changes

First stable release of advanced implementation of plant
