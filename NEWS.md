## Plant (development version)

### Breaking changes

These change the R-facing interface and require updating downstream code. Each
entry gives the `old -> new` migration; the `plant-update-interface` skill
(`.claude/skills/plant-update-interface/`) reads this section to migrate
products using plant.

* **phylloptim 0.7.0 -> 0.8.0.** Bumped deliberately, which is what the `==` pin
  below exists to force. 0.8.0 (traitecoevo/phylloptim#133) reports the seated
  curve's `lambda_emergent` through `leaf_solve()` and fixes a stale
  `lambda_emergent_` after the two terminal exits on a reused `Leaf`.

  **Inert for plant, and by construction rather than by luck:** the header diff is
  24 lines, being one element appended to `operating_point_values()` and two
  `lambda_emergent_ = NA` writes. plant reads neither that vector nor that member.
  Verified anyway — suite unchanged at 2985 pass, and the scenario scorecard's
  `offspring_production` bit-identical across all 8 scenarios against the same
  plant source built on 0.7.0.

  ⚠️ Newly available and **not** yet used here: `lambda_emergent` is the seated
  curve's own `(dC/dpsi)/(dE/dpsi)`, where `marginal_cost_water()` is TF24's price
  whatever curve ran. On `TF24_floor` the two differ by the price floor
  (32996 vs 82996 at `lambda_o = 5e4`), so anything reporting a marginal cost of
  water from this model wants the emergent one.

* **TF24's leaf parameters carry phylloptim's names (#634).** `TF24_Pars` renamed:
  `p_50 -> stem_P50`, `c -> stem_c`, `b -> stem_b`, `beta2 -> TF24_beta2`,
  `g1_TF24 -> TF24_cost_scale`. `Leaf()`'s arguments move with them, and
  `root_p_50 -> root_P50`. Same values, same equations: the model-version snapshot
  changes names only, with **no value moved** on either TF24 or TF24f.

  ⚠️ **`g1_TF24` was never a stomatal slope.** It is the TF24 hydraulic *cost
  scale*, handed to `Leaf()` positionally, while the leaf separately has a real
  Medlyn `g1`. One name per quantity, spelled the same on both sides, is what
  makes the hand-over in `prepare_strategy()` checkable by eye.

  `TF24_hyperpar` emits the new names, so anything reading `ret[, "g1_TF24"]` or
  `ret[, "p_50"]` needs updating. `stem_b` and `psi_crit` are TF24_hyperpar's
  *reporting* copies -- phylloptim derives its own from `(stem_P50, stem_c)` -- so
  setting either does not change the curve. `root_b` keeps its name: it is the
  Weibull *scale*, a different quantity from `root_P50 = root_b*(ln 2)^(1/root_c)`,
  and plant still converts at hand-over.

* **`TF24_Pars` gains `TF24_floor_lambda_o` (#634), defaulting to 0.** TF24 and
  TF24f are now seated on phylloptim's `TF24_floor` cost curve, and this is its
  one parameter: the price of water as transpiration goes to zero, in
  `umol C (kg H2O)^-1`.

  Every conductance-loss cost -- TF24's included -- prices water at zero as `E`
  goes to zero, since no conductivity is lost when nothing flows, so water is free
  precisely when it is abundant. `TF24_floor` splits the cost into a part depending
  on potential alone and a part linear in the flux, `Theta(E) = Theta~(psi) +
  lambda_o*E`, with `Theta~` being TF24's own cost at TF24's own traits. **TF24 is
  `TF24_floor` at `lambda_o = 0` at identical parameter values**, so "is TF24
  missing a price of water?" is a one-restriction question.

  **Nothing moves at the default.** phylloptim asserts the reduction bit-for-bit
  rather than to a tolerance, and plant's whole suite passes with every pinned
  baseline unchanged. `scientific_version` is deliberately **not** bumped: no
  equation and no default output changes.

  A new aux, **`shadow_cost`**, is reported beside `profit` -- so `aux_size` goes
  11 -> 12 (12 -> 13 with `collect_all_auxiliary`). It is `lambda_o * E` on this
  curve and exactly 0.0 on the other seven, and zero on another curve means "this
  curve does not separate the two", not "this curve's cost is all realised carbon".

  ⚠️ **The price is a shadow price, not a cost, and growth is billed on
  `profit + shadow_cost`.** `profit` is the *objective* and deducts both terms;
  `Theta~` is carbon actually forgone, while `lambda_o` is the value of water in
  its best alternative use, and paying it loses no carbon. Growing on the
  objective would tax the plant by carbon it never spent -- measured on a 5 m
  plant at PPFD 1800 and theta 0.25, the objective understates the carbon kept by
  2.3% at `lambda_o = 1e4`, 9.4% at 5e4, 15.4% at 1e5 and 22.6% at 2e5.

  Its scale is set by the leaf, not chosen freely: phylloptim's own marginal cost
  of water runs ~9e4 to 3e5 in these units at its defaults, so a value far outside
  that band pins the optimum against a bracket bound.

* **`phylloptim` and `odelia` are pinned with `==`, not `>=` (#634).**
  `LinkingTo: odelia (== 0.3.1), phylloptim (== 0.8.0)`. plant's regression
  baselines are bit-exact and both packages are compiled *into* plant, so under
  `>=` a later upstream release silently changes plant's arithmetic and the first
  sign is a red baseline with no local change to explain it -- a bisect across two
  repositories. Under `==` the same upgrade fails at install time, naming the
  version, before anything is built. Upgrading is then a deliberate act: bump the
  pin, rebuild, read the baseline diff, record it.

* **Migrated to phylloptim 0.6.0 and odelia 0.3.1 (#622).** `Leaf()` now takes
  `p_50` and `root_p_50` instead of `b`/`psi_crit` and `root_b`/`root_psi_crit`:
  phylloptim parameterises both vulnerability curves on P50 and derives the rest.

  ⚠️ **`b` is not `p_50`.** They are different quantities, both in MPa, so swapping
  them compiles, runs, and silently describes a curve 1.1465x too wide at TF24's
  defaults. Migration: pass `p_50` directly where plant already has it; convert a
  literal `b` with `p_50 = b * (log(2))^(1 / c)`. Verified -- handing over
  `(p_50, c)` reproduces plant's `b`, `psi_crit`, `root_b` and `root_psi_crit` bit
  for bit.

  Also: `optimise_psi_stem_TF()` / `optimise_psi_stem_Sperry()` -> configure then
  solve, `l$set_model("TF24", "stem")` followed by `l$optimise()`;
  `hydraulic_cost_Sperry()` removed; `profit_psi_stem_Sperry()` is ProfitMax
  underneath; `lambda_analytical_` removed. `psi_crit` and `stem_b` are newly
  readable, since they are now derived rather than supplied.

  **TF24 output moves by +0.36%** (net reproduction 3.671011e-25 -> 3.684168e-25
  on the reference run; 1055 -> 1065 ODE steps). Bisected to one upstream constant:
  `kg_to_mol_h2o` changed from a hard-coded `55.4939` to `1/molar_mass_h2o =
  55.509298`, +2.77e-4 relative. It enters the E -> g_sc conversion, so stomatal
  conductance carries exactly that shift and ci, assimilation and profit follow;
  the demographic integration amplifies it. The new value is the more accurate one
  (water's molar mass is 0.018015 kg/mol). odelia 0.2.1 -> 0.3.1 and every
  pre-#126 phylloptim commit were measured inert.

  FF16 and K93 are untouched. `scientific_version` was deliberately **not** bumped,
  so archived TF24 results from before this change carry the same `TF24@v8` tag
  despite differing by ~0.36%.

* **`Leaf$set_physiology()` takes `root_network`, not `root_carbon_per_leaf_area`,
  and `Leaf()` no longer takes `beta_R_H` or `beta_R_V`.** Requires
  phylloptim >= 0.2.0 (phylloptim #33). Migration:
  `l$set_physiology(root_carbon_per_leaf_area = x, ..., soil_depth = d, ...)` ->
  `l$set_physiology(root_network = phylloptim::root_network_from_carbon(x, soil_depth = d, beta_R_H = ..., beta_R_V = ...), ..., soil_depth = d, ...)`;
  drop `beta_R_H`/`beta_R_V` from `Leaf()` calls. `RootNetwork()` builds a network
  from resistances directly, for a caller who has those rather than a carbon
  profile.

  The leaf's supply solve reads two vectors -- `r_R_H_min` and `r_R_V_sum` -- and
  nothing in it touches root carbon, the 1/3 : 2/3 root split, the layer thickness
  or either `beta_R_*` constant. Those are a root-ARCHITECTURE model, which is now
  `TF24_Strategy`'s: `beta_R_H` and `beta_R_V` are still TF24 parameters at the
  same values, and `net_mass_production_dt` calls
  `phylloptim::root_network_from_carbon` into a strategy member before each
  `set_physiology`. Same move `leaf_specific_conductance_max` already made -- plant
  computes `kmax` from height and passes a scalar, because which
  conductance-versus-height model is in force is not the leaf's business.

  **No TF24 results change.** Verified bit-identical over 18 operating points
  (6 heights x 3 soil-moisture profiles) x 26 states, rates and auxiliaries,
  against `develop`; phylloptim's own 288-point golden file is bit-identical too.
  `scientific_version` is therefore unchanged.

* **`run_stochastic_collect()`'s environment field is `env`, not `light_env`.**
  Migration: `out$light_env -> out$env`. The old name was never produced by
  anything — `StochasticPatch::r_get_state()` had its environment leg commented
  out and the collector asked for a name that did not exist — so every element of
  that field was `NULL` and nothing can have depended on its contents. Code that
  tests for the name, or indexes the collected list positionally, needs updating.
  `env` is what the patch reports and what `run_scm()`'s collected output already
  calls it. See the corresponding entry under Minor changes.

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
    TF24 counts, which match exactly. (Those counts were subsequently re-derived when
    the stochastic arrival schedule started scaling with patch area — see Minor
    changes — so the file no longer shows the matching values. The equivalence that
    was established here still holds; it is just no longer the thing the test pins.)
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

* **The NSC storage pool is bounded by the shape of its own flow (`TF24@v9`,
  `TF24f@v9.1`).** `dS/dt` was `net_flux > 0 ? net_flux : floor_gate * net_flux`
  -- one signed flux with a gate applied in one direction. It is now a charge and
  a drain, each non-negative without a test and each limited by the room the
  other has (#609):

  ```
  charge = Ppos * (1 - G)          drain = Ppos - P
  dS/dt  = charge * (1 - r)  -  drain * r   =   charge - (charge + drain) * r
  ```

  The split is exact -- `charge - drain = P - Ppos*G` is the old net flux -- so
  what changes is where the limits sit. Both bounds now follow from the form,
  with no scale left to choose: at `r = 0` the rate is `charge >= 0`, at `r = 1`
  it is `-drain <= 0`, and past either boundary the corresponding limiter goes
  negative and pushes back. The pool is a first-order filter on production, so
  `d(dS/dt)/dS` is the constant `-(charge + drain)/S_max`.

  **Storage previously had no upper bound at all**; `r = min(S/S_max, 1)` clipped
  the *read* while the state ran past capacity. Measured on a one-species stand
  at the full 105.32-year lifetime: the state reached **1.035** of capacity and
  **52.7 per cent** of cohorts sat at or above it, so for half the stand `dr/dS`
  was exactly zero and the mortality that reads `r` could not respond. It is now
  0.000 per cent at every lifetime from 2 to 105.32, with the reserve fraction
  peaking at 0.800 -- the value a newborn is seeded at -- and a median of 0.62
  where it was exactly 1.

  **The derivative was discontinuous where the net flux changed sign**, which is
  what motivated the change. `d(dS/dt)/dP` stepped by `1/floor_gate =
  (r + 1e-3)/r` across that point -- measured 1.972 against a predicted 2.000 at
  a reserve fraction of 1e-3, and unbounded as the pool empties. The new rate is
  smooth there.

  **A negative state is now refused rather than floored**, which is what lets
  both read-clamps go. The exact flow cannot leave `[0, S_max]` but a finite
  Runge-Kutta step can, and it did: about 3 per cent of records at the default
  height at maturity, deepest at 4.5 per cent of capacity. The rate now calls
  `odelia::util::stop_domain` there and `Patch::ode_state_valid` declares the
  bound to the stepper, so the step is rejected, shrunk and retried instead of
  committed -- requires **odelia >= 0.3.1**. Crossings go to none, on the
  birth-date coordinate and the height one alike. `max(S, 0)` and
  `min(S/S_max, 1)` are gone with the branch, which is #609's own closing
  requirement that both consumers read the same quantity.

  **The drain is limited in proportion to what the pool holds, and that is a
  change to the drawdown as well as to the bound.** A plant at half reserves
  draws at half the deficit rate, where the old gate drew at essentially the
  full rate until the pool was within 0.1 per cent of empty. What it does *not*
  change is how much carbon goes unpaid -- what a plant pays from reserves is
  exactly the pool's depletion in either form -- so the integral over a drawdown
  is the same and only its schedule differs, which is the buffering #517 asked
  for.

  **A narrower drain limiter is what makes this shape necessary, and it is
  measured rather than argued.** Limiting the drain by `r/(r + D)` puts an
  *attracting fixed point* at `r ~ D` whose relaxation time is `S_max D / drain`
  -- 0.64 hours at the stiffest record on the reference stand, a 0.41 m seedling
  with a storage capacity of 1e-6 kg, against a solver stepping in days. The
  shipped rate has the same knee but falls through it into the region where
  `max(S, 0)` makes the rate identically zero, so it never resolves it; a rate
  that holds the bound has to. Accepted ODE steps over one 10-year run:

  | rate | accepted steps | wall |
  |---|---|---|
  | shipped `max(S,0)` and `net>0 ?:` | 940 | 8 s |
  | drain limited by `r/(r + 1e-3)` | 12 865 | 114 s |
  | drain limited by `r` | **488** | **3 s** |

  So the form that holds the bound is also **1.9 times faster than the one that
  shipped**, and the step count stops growing with height at maturity. Note this
  is not a case for a change of state basis: an eigenvalue at a fixed point is a
  property of the vector field and not of the chart, so `S = w^2` or `log S`
  would leave it exactly where it is.

  **Offspring production moves on every fixture.** The height-coordinate answers
  fall by a factor of about 2.7 (82.09 -> 30.22, 67.54 -> 23.20) and the
  birth-date ones by 19 and 27 per cent (287.16 -> 233.06, 59.53 -> 43.64). Most
  of that is the fill limiter rather than the drain, since the median cohort of a
  mature stand previously sat clipped at a reserve fraction of 1. The height
  coordinate amplifies any change to the pool because its compression term
  differences growth against height at fixed absolute carbon, so it reads the
  reserve fraction moving where the birth-date density rate never asks.

  **Carbon is still not conserved across the block.** The pool is capped by
  withholding the surplus rather than by spending it, so at capacity the charge
  the gate withheld -- `Ppos * (1 - G)`, about 1.2e-4 of production -- leaves the
  budget instead of charging the pool; at the empty boundary the unmet part of
  the deficit does the same. Both are properties of limiting a flux rather than
  spending it, and routing either into growth or into tissue loss is a larger
  modelling change than #609 proposes.

  **Not addressed here: `storage_prod_eps` is an absolute rate against fluxes
  spanning six orders.** It is 1e-4 kg/yr, and a plant's own maintenance flux on
  the reference stand runs from 4.60e-05 to 2.88e+01 kg/yr -- so the smoothing
  constant is between 3.5e-06 and **2.17** of a plant's whole turnover, and
  exceeds it on **39.8 per cent** of records. A plant at zero net production is
  therefore credited with `eps/2 = 5e-5` kg/yr it does not have, which is 0.5 per
  cent of an adult's maintenance and 29 per cent of a seedling's, and `Ppos`
  reaches **21.3 times `|P|` with the opposite sign** where the pool is
  stiffest. #609 predicted this and asked for it to be checked against a real
  germination state before acting; it now has been. Making it a fraction of a
  flux the plant has moves establishment and every plant near its compensation
  point, so it wants its own change and its own re-blessing rather than being
  folded in here where one census movement would have two causes.

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
  dependent, and needs re-deriving before it is relied on. *(Re-derived and
  retracted; see the entry under Minor changes & bug fixes below.)*

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

* **A dense TF24 stochastic run throws at the default ODE step cap** (#599). The
  soil water balance is stiff — the conductivity curve's exponent is
  `2*n_psi+3 ≈ 16`, so `K_sat/dz` is ~543/yr at the defaults — and at
  `Control()`'s `ode_step_size_max` of 5 yr the explicit RKCK stepper diverges
  rather than losing accuracy. A failing run reaches
  `theta = [-51.4, 52.0, 0.146, nan, nan]`, and once the rates are non-finite the
  step is not even rejected, because `adjust_step_size` derives its error ratio
  from a NaN. Downstream, the leaf's collar root-find is handed a soil potential
  no retention curve can produce and throws
  (`find_root_psi ... do not bracket the root`). Measured over seeds 1–40 at
  `patch_area = 1` (~105 individuals m⁻²): 17 of 40 runs throw on the code before
  this release's stochastic-solver work, 5 of 40 after it. Setting
  `ctrl$ode_step_size_max <- 0.05` gives 0 of 60; the relationship is not monotone
  (0.2 fails *more* often than 5 does), so that value is measured rather than
  derived. It is left as a user setting rather than clamped in the library,
  because the right fix is a per-environment bound or odelia's stiff RODAS stepper
  for environments carrying stiff state, not a global default that every model
  pays for. `test-stochastic-patch-runner.R` sets it for its own TF24 runs so the
  suite is deterministic; that is a test setting and changes nothing for callers.

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

### Minor changes & bug fixes

* **An out-of-domain hydraulic lookup now names the spline, the point, the domain and the caller** (#576), via odelia >= 0.2.2 and phylloptim >= 0.2.1, both already required. The message used to be odelia's bare "Extrapolation disabled and evaluation point outside of interpolated domain.", which named none of the four — and since the leaf reads the same transport spline from four places and holds a second spline that is its inverse, that sentence did not distinguish the cases that matter. Localising #576 meant instrumenting the call sites by hand to find out which spline was read and at what value; the answer (the *lower* end of the inverse, where the demanded flux is negative) is what showed the spline had been right to refuse. Error text only; no model behaviour changes. `test-leaf.r` now asserts the spline, which end was missed and the caller, rather than the whole sentence — which is phylloptim's to word.
* **TF24f's AD gradient branch mirrors the finite-difference branch — 21% faster, same numbers** (#576). `TF24f_Strategy::solve_leaf`'s two gradient methods were asymmetric. The finite-difference branch has, since #530, run `prepare_collar_solve()` once itself and shared it across its profit evaluations, and taken a zero gradient when that call reports the operating point *forced* by feasibility handling rather than chosen from an interval. The AD branch — the default, hence TF24f-only — went through `evaluate_root_collar_psi()`, which hides that return value, so it re-derived the soil-side caches on each of its three leaf evaluations and asked for a gradient even at a collar potential the feasibility analysis had already rejected. The AD branch now does what the FD branch does.

  Measured on a 10-year dry-start TF24f SCM run: **bit-identical `offspring_production` at all 26 points** of a 13-rainfall × two-gradient-branch sweep, and **2.61 s → 2.10 s** at rainfall 1 m/yr and **1.74 s → 1.41 s** at 0.04 (minimum of five; the FD branch is the control and does not move, 2.20 → 2.23 s and 1.19 → 1.21 s). The forced-operating-point exit is an ordinary state rather than an edge case, and since the storage-pool change (#619) a common one: it fires on 3.9% of leaf solves at rainfall 1 m/yr and **50.0%** at 0.04.

  This began as the fix for #576, where the AD branch died in a dry rainfall window (0.01 and 0.03–0.06 m/yr at θ = 0.005, five layers) with "Extrapolation disabled and evaluation point outside of interpolated domain" while TF24 ran — asking, in shutdown, for the stem potential that would carry a *negative* sap flux, i.e. wetter than saturation. **That throw no longer reproduces**, on this commit or on develop: phylloptim 0.6.0 closed it upstream. The regression test added here is therefore a guard against its return, not a reproducer, and the reason to make this change is now the redundant work and the asymmetry rather than the crash.

* **The scenario gateway scores numerical viability and persistence as separate
  axes, and no longer leads with a match rate** (#572). `scenario_summary()`
  classified a run as a success on `finite && total > 0`; with `birth_rate = 1`
  (what `build_scenario()` sets) `offspring_production` *is* R0, so five of the
  eight hydraulic scenarios returned R0 between 2e-15 and 6e-14 — numerically
  extinct — and were all recorded as `persisted`. All eight "succeeded", none of
  the five expected failures failed, and only one replaced itself: the gateway
  returned no signal. It now reports **numerical viability** (`n_ran`,
  `n_crashed`, `viability_rate` — did the model run, which is what the CSV's
  "Model failure" means) and **ecological persistence** (`n_persists`,
  `persistence_rate`, at R0 >= 1) as separate, labelled axes, in the summary, in
  `scripts/run_scenario_gateway.R` and in the scorecard report. `n_match` /
  `match_rate` are still reported but demoted to *agreement with the CSV's
  crash-era expectations*: with the crashes fixed (#546/#552/#554) the match rate
  mostly measures how well those predictions have aged, not model quality — which
  is why fixing the model *lowered* it from 5/8 to 3/8. `status` and `outcome`
  are unchanged and the CSV is not rewritten.

* **The scenario gateway integrates in birth date, not height** (#590), via a
  new `scenario_control()` that supplies the `Control` every entry point in the
  framework now defaults to. The package default is untouched.

  Every scenario in the gateway is TF24, and TF24 is the model the two density
  coordinates genuinely disagree on: the compression term is the total
  derivative of growth along a cohort's trajectory, which equals `dg/dh` only
  when growth is a function of size, and TF24's reserve gate (#517) breaks that.
  So the gateway had been scoring the hydraulic model through a compression term
  that is wrong for it.

  Measured on the eight scenarios at `max_patch_lifetime = 100`, the coordinate
  change raises R0 on every one — S02 2.4x, S06 8.9x, S08 15x, S05 25x, S07 47x
  — and moves **S01 across R0 = 1** (2.59e-01 to 1.49e+00). Persistence goes
  **1/8 to 2/8**. Numerical viability stays 8/8 and CSV agreement stays 3/8.
  Birth date is also about twice as fast here (36 s against 70 s).

  Note what that implies about the old guard: `observed` did not change on a
  *single* scenario across a 47x swing in R0. It tests `finite && total > 0`,
  which every scenario has satisfied since the crash fixes landed, so it was
  pinned at 8/8 and could not move. The baseline diff in
  `test-scenario-gateway.R` therefore now diffs `persists` as well as
  `observed`, and **the baseline is re-blessed** under the new coordinate.

* **Retracted: "reserve-gated growth largely excludes the slower species".**
  `test-strategy-tf24.R` recorded that as a finding about #517's NSC reserve
  gate. It is a property of the density coordinate, and #590 flagged it as
  needing re-deriving. Re-derived here, at the test's own configuration
  (`lma` 0.0825 / 0.10, `max_patch_lifetime = 5`):

  | coordinate | fast | slow | ratio |
  |---|---|---|---|
  | height | 67.32 | 2.775e-04 | 2.4e5 |
  | birth date | 287.2 | 59.53 | 4.8 |

  So the two species **coexist at comparable abundance**; they are not
  separated by five orders of magnitude. At `max_patch_lifetime = 30` the ratio
  is 2.7 (2707 against 1004, matching #590) — i.e. it *narrows* with patch
  lifetime, where progressive exclusion would widen it.

  What survives the retraction is the weaker claim: the faster species still
  leads, by roughly 3-5x. Reserve-gated growth disadvantages the slower species
  without excluding it. How much of even that gap is attributable to #517 as
  against the trait difference itself is not measured here, and should not be
  read into these numbers.

  The evidence that this is a wrong derivative rather than an under-resolved
  one is refinement. Over two halvings of the node spacing (88 → 175 → 349
  nodes) the birth-date answers are already converged — fast
  287.2/287.1/287.2, slow 59.53/59.66/59.69 — while the height answers are
  still climbing (fast 67.32/73.18/74.98) and the exclusion ratio does not
  shrink toward the birth-date one, it *grows*: 2.43e5/2.49e5/2.51e5.
  Quadrature error closes under refinement; a different derivative does not.

  The height-coordinate assertions are kept, since they pin what the package
  default still does, and a birth-date case is added alongside so the
  coordinate that is correct for TF24 is actually covered by a test rather than
  only described in a comment.

* **`run_stochastic_collect()` now reports the environment.**
  `StochasticPatch::r_get_state()` had its environment leg commented out, so
  unlike `Patch::r_get_state()` it returned only `time` and `species`. The R
  collector then read `light_env`, a name nothing had ever produced, so every
  element of that field came back `NULL`. The field is now `env`, matching what
  the patch reports and what `run_scm()`'s collected output calls it, and it is
  populated. This was inert while the stochastic solver held its environment at
  the initial state; now that the environment is integrated, the soil trajectory
  is real and worth reporting. **Breaking:** the returned list field is renamed
  from `light_env` to `env`. Nothing could have depended on its contents, since
  it was always `NULL`, but code testing for the name will need updating.

* **The stochastic arrival schedule now scales with patch area.**
  `stochastic_schedule()` passed `patch_area` into `stochastic_arrival_times()`'s
  third positional argument, which is `delta_t`, leaving `patch_area` at its
  default of 1. Two things followed. The arrival rate never scaled with area, so
  the expected number of arrivals came out as `max_patch_lifetime × birth_rate`
  whatever the patch size, and a 50 m² patch was seeded like a 1 m² one. The
  binning interval was also silently set to the area, which for a large patch
  left only a handful of intervals over which a variable birth rate was
  averaged. Arrivals now scale linearly with `patch_area` as intended, and the
  interval keeps its 0.1 yr default. `run_stochastic_collect()` is the only
  caller; its runs change accordingly, and even at the default `patch_area = 1`
  they are unchanged in expectation but not bit-identical, because the finer
  binning draws from the RNG differently. The seeded baseline in
  `test-stochastic-patch-runner.R` was re-derived and its `patch_area` reduced
  from 50 to the default 1, which keeps the ~105-individual stand that test has
  always actually run rather than the ~5300 that 50 m² now implies.
* **An empty stochastic patch no longer discards the environment's integrated
  state.** `StochasticPatch::compute_environment()` calls
  `Environment::clear_environment()` when the patch holds no individuals, which
  is right for the competition profile — nothing is casting shade — but TF24's
  override also restored the soil states and cumulative-flux accumulators to the
  values the run began with. Now that the environment is part of the ODE system
  that discarded what the solver had integrated, on every derivatives evaluation
  while the patch was empty. `clear_environment()` is now the competition
  profile alone; restoring state moved to a new `clear_state()` that only
  `Environment::clear()` calls, so resetting for a new run is unchanged. This
  matters most for a schedule whose first arrival is some way into the run,
  which is `run_stochastic_collect()`'s default.
* **The stochastic solver now integrates the environment.** `StochasticPatch`'s
  ODE system was the species alone: `ode_size()`, `set_ode_state()`,
  `ode_state()` and `ode_rates()` did not chain through the environment as
  `Patch`'s always have, and `compute_rates()` never accumulated resource
  consumption, so `Environment::compute_rates()` was never called. An
  environment carrying ODE state was therefore held at its initial value for a
  whole run. TF24 carries nine such states — five soil-moisture layers and four
  cumulative fluxes — while FF16 and K93 carry none, which is why this went
  unnoticed. Stochastic TF24 runs change: soil water recharges from 0.214 to the
  drainage equilibrium 0.3106 and is then drawn down as leaf area grows,
  tracking the SCM's own trajectory on the same drivers to ~1e-4, and total leaf
  area at patch ages 1–3 moves from 11–16% below the SCM to within 1.4% of it.
  FF16 and K93 stochastic runs are bit-identical, given the guard described in
  the next entry. The consumption vector is sized by the environment's
  `n_resources()` rather than its ODE width, so TF24's four diagnostic flux slots
  are not mistaken for resources.
* **`StochasticPatchRunner::reset()` integrates to the first arrival with error
  control.** `reset()` advanced from time zero to the first scheduled arrival in a
  single `advance_fixed` step, which for `run_stochastic_collect()`'s default
  first arrival — uniform on (0, 50) — is one step of up to fifty years. That was
  harmless while nothing was integrated over it. Once the environment's own states
  are, it drives soil water out of range and nothing establishes: fixing the state
  handling alone produces nine failures and a runner that introduces no
  individuals. The leg now uses `advance_adaptive`, as every other advance in the
  runner already does — but only when there is state to integrate. An empty patch
  whose environment carries none (FF16, K93) keeps the fixed step, and that guard
  is what preserves bit-identity for those two. `advance_adaptive` walks the
  step-size controller up to `ode_step_size_max` even on a zero-width system, and
  `step_size_last` survives into the first real step, because only
  `SolverInternal::step()` writes it and `step_to()` — which `advance_fixed`
  drives — does not. Ungated, the first step after the first arrival therefore
  begins from a rejected five-year attempt rather than
  `ode_step_size_initial`, and FF16's collected trajectory moves by up to 1e-5
  relative: inside `ode_tol_rel = 1e-4`, so a valid realisation either way, but
  not identical.
* **A stochastic recruit's strategy-specific initial states are seeded from its
  birth environment.** `StochasticSpecies::introduce_new_node(environment)`
  computed the new individual's rates without first calling
  `Individual::set_initial_states()`, which the deterministic path calls from
  `Node::compute_initial_conditions()`. A TF24 seedling was born with an empty
  carbohydrate store rather than `a_st3` of its storage capacity. FF16 and K93
  do not override `set_initial_states()`, so they are unaffected.

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

* **CI is faster and the test suite is leaner.** The suite ran 147 s serially,
  with one file (`test-tf24-arid-corner.R`) at 24% of it, and tests were ~65% of
  each `R CMD check` job. Now 124 s (-16%), and no single file is over 19%.
  * `test-tf24-arid-corner.R` 35.5 s -> 12.7 s. It computed the same dry-start
    TF24 run four times; the completed run is now memoised and shared (the
    `short_run()` pattern from `test-density-coordinate.R`). #571's nine-point
    soil-moisture sweep is cut to the three values bracketing the
    residual-moisture floor, the full set being recorded in #571.
  * Dropped `"a failed light spline reports the patch state that caused it"`. Since
    #574 fixed the cause its `skip_if()` fired every run, so it paid for a full
    SCM run and asserted nothing. See the comment in the file for how to bring it
    back as an unconditional test.
  * Dropped `"TF24 water budget closes independently of layer count"`: it re-ran
    the same three configurations as the test above it to re-assert a bound that
    test already checks per layer.
  * Benchmarks no longer run on every push and PR — `run_plant_benchmarks()` was
    costing 3-5 min per trigger on a macOS runner while nothing compared its
    output against a baseline. Now weekly plus `workflow_dispatch`, on Linux.
  * Both workflows cancel superseded in-flight runs on the same ref, and
    `R-CMD-check` gained a `ubuntu-latest` job — commented out since CI was first
    set up (#295), so plant had never had a green Linux run.
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
