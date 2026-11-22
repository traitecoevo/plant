# Hazards — plant

The invariants in this code that a later edit can break without any test failing at once, each with the reason it holds. [agents.md](../agents.md) §11 is the short checklist and points here; [measurements.md](measurements.md) says where the baselines live and how a change is measured.

**The numbers are stable.** Comments may cite "hazard N", so a hazard keeps its number for life. A resolved hazard keeps its slot and says what replaced it; a new one is appended at the end.

## Solver and events

### 1. What "bit-identical" means for a run

A run is bit-identical when its set of stop times and its actions are unchanged. Node introductions are already stop-and-jump events, so a change that adds or moves no stop time and no action must leave every number unchanged, and is checked with `identical()` on raw arrays. Adding an event adds a stop time, which changes the adaptive step sequence and so legitimately moves the trajectory at the solver's tolerance. Say which of the two a change is, so a reviewer does not read an expected movement as a regression.

### 2. The step size is inherited across a leg boundary

odelia's `step_size_last` survives `set_state_from_system()` and every `advance_*` call, so a leg that starts after a discrete change (a pulse, a harvest, an introduction) begins at whatever step the previous quiet leg grew to, up to `ode_step_size_max` (5 years by default). It is what triggers TF24's failures under rainfall pulses (PR #608), and pulses create far more leg boundaries than introductions do. Setting `ode_step_size_max = 0.05` on a pulsed run is a measured guard, not a fix. The same channel means swapping one `advance_*` for another moves every later step, even on an empty ODE system.

### 3. A guard a trial stage can reach must reject the step

`phylloptim::util::infeasible_error` and `odelia::util::DomainError` are siblings, both deriving from `std::runtime_error`, so odelia's step handler cannot see phylloptim's throw: it escapes the step and kills the run. `TF24_Model::solve_leaf()` translates it. Any new call into phylloptim from a rate evaluation has to go through the same translation.

The same holds for plant's own guards on the rate path. A Runge–Kutta trial stage evaluates rates at states the run will never keep, and a stage can carry a newborn to a log density of hundreds or below its stem's terminal segment. A guard such a stage can reach (`Patch::check_finite_ode_state`, `Species`' competition contribution, `TF24_Model::update_leaf_specific_conductance_max`) throws `odelia::util::stop_domain`, which rejects the step so the stepper retries it smaller; written as `util::stop` it kills the run at a state no accepted step would reach, which was the robustness sweep's most common failure. The rule is about which guards, not where they fire: a guard a trial stage can reach is a `DomainError` wherever it fires, and a guard only preparation or an introduction can reach stays a `util::stop`. Where the step cannot shrink the message survives: at the minimum step odelia stops with "Cannot leave an invalid state at t = ...: <message>", and outside a step (at preparation, at an introduction, from R) a `DomainError` reaches R as an error of that class with the message. On a pinned run a row with a recorded step size is taken at that size (`step_by`), so a guard there reaches R as itself; a grid-only row goes through `step_to`, which splits the interval on a `DomainError`, so that run now finishes and a later sweep or replay refuses the split row.

### 4. The event queue is built by a linear scan

`Schedule::entry_at` walks the schedule from the start for every event it places, so building a queue is quadratic in its length. It does not show at the hundreds of entries of a normal run; at the 10^4 to 10^5 events of a long pulsed run it will.

## The reverse sweep

### 5. `Patch::assign_from` must equal a rebind

`Patch::assign_from(from)` must leave the patch equal to `from.rebind<value_type>()` in every member a recording reads, for a `from` built from the same parameters. It does not rebuild the parameters, node schedule or disturbance regime. odelia's sweep holds one active patch and re-seats it at every width through this hook (`AssignsFrom` in odelia's `ode_interface.hpp`); `rebind()` constructs a fresh patch and calls `assign_from`, so a member added to one is added to both. A member that a rate evaluation writes and a rebind would reset, but that `assign_from` leaves alone, gives finite wrong gradients and raises nothing. To check a change, build once with the hook renamed, so odelia falls back to a rebind per width, and compare `stand_gradient()` bit for bit.

Three parts of this are easy to undo by tidying:

- Each species' strategy is re-assigned in `assign_from`, not kept. A TF24 strategy holds a `phylloptim::Leaf` reused across evaluations; a rebind copies the passive strategy's, and so does `assign_from`.
- The environment is rebound afresh, so no cache a previous width's evaluations left (TF24's soil-potential and driver memos, the clamp counters) reaches the next.
- `Patch::rebind()` copies every `Parameters` field that `validate()` reads, the fitness-window fields included, or the rebound patch runs a different regime with nothing raised.

### 6. The light spline's coefficients are double at every scalar

odelia's `hermite_spline` keeps its span coefficients in double, so a read at an active scalar is the passive read with the knot data's partials attached. Hand the spline its knot vectors by `std::move`: a copy is a recorded statement per knot.

### 7. plant compiles odelia's tape layout

`XAD_REDUCED_MEMORY` is set in odelia's vendored `Config.hpp`, and plant compiles odelia's tape layout into its own objects. plant must be built against the odelia it links, which the `==` pin in `LinkingTo` enforces.

### 8. Four things tried on the sweep's cost and decided against

Each looks like a win and is not; try one again only with a new reason:

- Moving stage states into the System (about 1% in a spike). XAD's move leaves the value and drops the slot, so a System that reads an iterator entry twice gets a value with no derivative and nothing reports it.
- Field-reduction rows (2% at one metric, 6% at five). They are a second copy of `Species::field_splits_nodes` and `CanopyShape::Q_and_q` that must agree with the first bit for bit.
- Building the excluded field only below the boundary node (1%). It threads a parameter through `Patch`, three environments and `ResourceSpline`.
- A Gram–Schmidt basis of output adjoints. The only dependence seen is one row proportional to another, which `adjoint_multiples` catches in a pairwise test.

## Soil water

### 9. The inter-layer cascade is free drainage, and has no internode distance

The flux from layer `i` into layer `i+1` is the donor's free drainage `K(theta_i)`, Darcy's `q = -K(dpsi/dz + 1)` with the matric gradient dropped, less what a nearly full receiver refuses (hazard 12). Neither part carries an internode distance, so unequal layer widths need nothing rescaled, and mass conservation `sum(dz[i] * dtheta_i/dt) = infiltration - drainage - uptake` is exact at any width distribution (`test-tf24-water-budget.R`, on four graded profiles). `water_flux[i]` is the flux out of the base of layer `i` as accepted by the layer below, not `K(theta_i)`. If a capillary or diffusive flux is added, as soil evaporation will want, its distance is `z_mid[i+1] - z_mid[i]`, not `dz[i]`, and the interfacial conductivity is a distance-weighted harmonic or geometric mean of the two cells. The two distances agree only for equal layers, so a uniform-profile budget test passes on the mistake and the graded one does not.

### 10. Two routes to the soil geometry must stay as they are

- `TF24_Model` passes `phylloptim::layer_thicknesses(z)` to the root network, at both call sites (the rate evaluation reads the copy `TF24_Environment` computes once when its geometry is set, `z_thickness`; the gradient path takes it from the leaf's own profile), and must not switch to `environment.dz`: on the default 1.5 m over 5 layers, layers 3 to 5 differ from `depth / n` in the last bit, so the switch moves every TF24 number. Make it at the next rebaseline.
- `set_soil_number_of_depths` builds `z[i] = (i+1) * delta_z` and must not be unified with `set_soil_layer_widths`' running sum. `z` feeds `Q()` and the gravitational head, so unifying them re-baselines every TF24 number for no gain. The geometry test in `test-environment-TF24.R` compares the two.

### 11. A thin wet layer has a stiffness floor

Layer widths are validated for being finite, positive and not thinner with depth; nothing bounds how thin the top layer may be. The cascade's diagonal `K'(theta)/dz` grows as the layer thins and as it wets: at 0.99 θ_sat it is 1.8e4 yr⁻¹ at the default 0.30 m (stable step 1.1e-4 yr) and 2.6e5 yr⁻¹ at 0.02 m (7.6e-6 yr, about 4 minutes); at 0.5 θ_sat and 0.02 m it is 8.5. So the hazard is the post-rain transient. `Control()$ode_step_size_min` is 1e-6 yr, which the stable step reaches at a layer of about 2.6 mm; below that the controller cannot be stable at all, which is a hard floor on layer thickness in this model. `soil_widths_graded()`'s documentation (now in `fieldtrial`) carries the same warning. A thin top layer also shifts runoff: infiltration keys off layer 0 alone, through `(theta_0/theta_sat)^b_infil`, and `a_infil` and `b_infil` were calibrated against a 30 cm surface layer, so a thinner layer 0 saturates sooner in a pulse and sheds more.

### 12. Saturation is a barrier the balance keeps, through the receiver's refusal

A layer takes of its supply from above `supply - w (supply - smooth_min(supply, out))`, with `w = (theta / theta_sat)^b_infil` and `out` what the layer passes on below (`TF24_Environment::accepted_inflow`). At saturation its inflow is at most its outflow, so its rate is not positive while its roots draw on it or leave it alone; above saturation `w > 1` pushes it back. What a layer refuses stays above it, and at the surface is runoff (`sum_rainfall - sum_infiltration - sum_pulse_runoff`), so the budget closes. Three things hold this up:

- **The fluxes are found from the bottom up.** `out` is the flux the layer below actually accepts, not the layer's own `K(theta)`. Comparing the supply with the receiver's own drainage instead overfills a layer above a full one: tried, a middle layer reached 4 against 0.428.
- **The smooth minimum never exceeds the smaller argument**, which is what makes the inflow at saturation at most the outflow. A smooth minimum that can overshoot (a hyperbola with an additive width) breaks the barrier, and one that can go negative for widely separated arguments pumps water upward. Conductivities span ~25 orders of magnitude, so it is evaluated from the larger argument.
- **The weight's slope at saturation sets the stiffness of a perched layer**, one held full above a slower one: about `b_infil * (supply - out) / (theta_sat * dz)`, of the same order as the soil's own `K'(theta_sat) / dz` at `b_infil = 8`. Measured on a column whose K_sat falls to 0.5 m/yr at the base, under `Control()`: 8261 steps over 10 years at 8, 12434 at 16, 22335 at 32. `control_tf24()`, which substeps the soil, took 377 over 3 years. A steeper weight leaves a column that can pass on its supply closer to free drainage (at 8 the default TF24 stand's offspring is 3% below free drainage alone; at 16, 0.05%), at that cost.

The one inflow not limited is water the roots release into a layer (a negative uptake): released into a full layer it overfills it, and `Patch::ode_state_valid` stops the run, naming the layer. Error control can also carry a step past saturation; the same check rejects it and the step is retried smaller.

Do not change the soil physics to fix a pulse failure. What kills a pulsed run is an inherited step (hazard 2) or an exception thrown where the stepper cannot reject it (hazard 3); PR #608 has the measurements. Clamps in the water balance were measured and do not help a pulsed run (5 of 40 failing seeds became 4 of 40), and a θ ≤ θ_sat clamp in `psi_from_soil_moist` would hide an overfilled layer rather than prevent it. What the pulse feature did need, a cap on the pulse jump at free capacity with the excess to runoff, is in `TF24_Environment`.

## The continuous patch

### 13. A cohort leaves the integrated set only by passing the absorbing height

Each species has one top cohort pinned at the height where its strategy allocates 1 − ε of its maximum to reproduction (`Control$absorb_epsilon`). At every introduction instant a live cohort that has reached that height moves its share of every trapezium into the top cohort and stays in the node list as a frozen zero-density grid point. There is deliberately no survival floor, for any cohort, scored ones included: retiring a cohort discards the biomass old trees still hold. The cost is accepted: on the pulsed TF24 stand at b*, 65 of 101 window cohorts are still integrated at year 350, which makes the window run 2.5 times what it cost under a floor.

### 14. The top layer's limits

- The layer is birth-date coordinate only, and is refused on the height coordinate.
- The reverse sweep differentiates the layer because absorption runs inside the insertion map it transposes (`Patch::apply_insertion`), at the sweep's scalar. So everything an absorbed cohort carries into the layer, and the height the layer is pinned at (`asymptotic_height()`, from `hmat` and `a_f2`), must stay at the strategy's scalar. A `to_passive` there drops a column's dependence and leaves every number finite. `test-absorb.R` referees `hmat`, `a_f2` and `S_D` for this reason.
- With the layer on, the sweep refuses a harvest or climate extreme that reached any cohort (`require_transposable_events`). Unlike a node, which carries log density, the layer carries its density and weights as they are, so an event scales them, and the sweep crosses an event as the identity. An event at an introduction instant also falls between the recorded state and the absorption the sweep replays from it. Measured across a harvest on the layer: `S_D`'s fitness column 17% off, every number finite. To lift the refusal, the event's scaling has to be transposed, and absorption has to be replayed from the state after the event.
- Refining the schedule with the layer on is refused, because `refine_schedule` pairs per-node errors with schedule times.
- The absorbing height is a selection on state. The sweep holds it as the run chose it. A trait perturbation that changes which cohorts merge makes a pinned-step finite difference disagree with the sweep at that kink, so a referee must check that both of its runs absorbed the same cohorts.
- The layer also makes the canopy top a `max()`. `height_max()` is the larger of the layer's height and the tallest live cohort's, and a cohort growing through the absorbing height between two introductions switches it, so the light field's knots follow a different height from that stage on. The sweep differentiates the branch the run took. A difference whose step moves the switch past a stage straddles a kink that its absorption record does not show: on the resumed FF16 stand in `test-absorb.R`, a relative step of 1e-5 in `hmat` put the net reproduction ratio 6e-5 off the sweep, and 1e-6 agreed to 1e-7. Shrink the step until the difference stops moving before reading a disagreement as a gradient error.
- A younger TF24 cohort can overtake an older one (the reserve gate). Absorption reads each cohort's own height, so the younger goes first and the older when it arrives.
- A spent layer is emptied, not left subnormal: at each introduction instant a density or weight below the smallest normal double is put to zero, keeping the seeds it accrued. Left alone, a decay step at the bottom of the subnormal range rounds back to where it started, and the layer pinned `height_max` and the light field's knots at the absorbing height. Anything that changes how the layer decays has to keep that.
- A cohort born with no density is frozen at birth, so a species at birth rate zero would read zero per-capita output; `Patch::push_nodes` refuses it wherever the regime weights the cohort.
- R-side integrals (`integrate_over_size_distribution()`, `tidy_species()`, `tidy_patch()`) see the nodes only; `tidy_species()` warns when the layer carries density. Read the layer from `species$top_layer` or the snapshots' `top_density`.

### 15. Long runs and resumes

- A resumed run is not bit-identical to a straight one at the seam; say so beside any number that crosses one.
- A resume with no introductions left stops where it starts, because the deterministic runner's `$complete` counts schedule entries. It is refused when the run begins; keep at least one introduction ahead of a resume.
- `make_initial_state()` on the birth-date coordinate stamps every node at time 0 and is refused for more than one node per species; it is not the resume entry point.
- `history` under `collect` is a Patch copy per leg; never turn it on for a long run.
- An invader must stop wherever the residents' run stopped. A pulse leaves the rates stale until the runner re-reads the system at its stop, so an invader that stepped through one began its next step from rates taken before it, and a copy of a pulsed TF24 resident scored 36% low. `run_mutant()` therefore gives invaders every action of the residents' run; a filter on event types there brings that back.
- A species' first cohort takes establishment only from its birth on (`Species`' open interval), while a later cohort's cell reaches back to the cohort before. An invader's first cohort therefore differs from the residents' cohort of the same date unless the residents' is also their first. `run_mutant()` refuses a window that weights it; a resume is one way to get there, an invader schedule that starts at the window is another.
- The introduction grid can make a fixed point a periodic orbit: between introductions the field swings (1.8% in held R for FF16 at b* with a cohort a year), and a field held at an introduction is the top of the swing. Check a held field inside an interval, by ending the resident's run off the grid, and refine the grid until the swing is below what the comparison needs (#665).

## TF24 variants

### 16. TF24t

- `TF24t_Model::rebind` is deleted and must stay deleted until TF24t has a derivative of its own. The inherited one returns a `TF24_Model<U>`, so it would rebind a TF24t patch as a TF24 patch: no thermal states, no damage, every number finite.
- The thermal parameters are top-level strategy fields, not `TF24_Pars`, as TF24f's are. A trait matrix reaches them as numeric top-level fields (`generate_strategy` assigns them beside `pars`), and because they are outside `pars` the drift guard snapshots them only through the explicit `top_level_science` list in `test-model-version.R`: a new thermal parameter must be added there.
- TF24t refuses the DeepCrown shading model: the leaf's damage rates are read at its leaf temperature after the crown loop, which under DeepCrown would be the last node's.
- At non-damaging settings a rate evaluation is bit-identical to TF24 with `pars$use_energy_balance = 1` (alpha = 0, k_i = 0, topt_offset = 0, all costs 0).
- The leaf's damage rates are per day; TF24t multiplies them by `DAYS_PER_YEAR = 365.25` for the yearly ODE clock, except leaf turnover and growth dilution, which are already per year.
- `Individual$set_initial_states` must run before a TF24t plant is grown with `run(individual, ...)`. The growers do not seed, and an unseeded storage pool starts empty, where a plant losing carbon has every step refused.

### 17. Flexible allometry (TF24's leaf and sapwood states)

- Height's monotonicity holds by construction: `dh/dt` is a product of non-negative factors. It has to stay structural, because the solver cannot check it.
- Every expression in the carbon-to-growth and carbon-to-mortality maps must be smooth: no `min`, `max` or `?:`. Shape the flow so bounds follow from its form, as the storage rate in `TF24_Model::compute_rates` does.
- `competition_effect` is load-bearing. `compute_rates` never calls `area_leaf`; it reads the aux written only by `update_dependent_aux`, so an allometry change not routed through that function passes every test that calls `area_leaf` directly.
- `Individual::set_ode_state` refreshes aux per slot, in slot order, so an aux that depends on several states sees the old value of any later slot when an earlier one fires. TF24's `update_dependent_aux` therefore recomputes the carried leaf area, in full, on both the height slot and the leaf-departure slot; whichever fires last leaves it consistent. Keep any such aux a full recomputation on every slot it depends on.
- The allocation derivatives are differentiated by hand and move together: `dmass_sapwood_darea_leaf` carries an `(a_l2 + 1)` factor that is the product rule against `H(A)`, `dmass_bark_darea_leaf` inherits it, and `darea_leaf_dmass_live` is their reciprocal sum.
- Canopy shape stays `Q(z/height)`, so a smaller leaf area at unchanged height is a sparser crown of the same depth, not a shallower one.
- `set_initial_states` is not virtual; TF24f must keep chaining to the base, or the leaf-area seed is silently dropped.
- Height crossings are routine on the height coordinate once leaf area is free; `scale_node_densities` selects by height, and the competition reduction runs over the abscissa set rather than the node list so a crossing cannot give it a negative width (`test-canopy-methods.R`'s permutation test).
- `test-model-version.R` snapshots only `pars` and `control` defaults, so adding a state without a parameter passes it. Bumping `scientific_version` for that is review discipline.
- The six flexible-allometry parameters carry no gradient column (they are `undifferentiable`): at `a_pl0 = 0` four of them multiply zero, and a zero column reads as an answer. A gradient taken with `a_sw > 0` is refused, because the sapwood sensitivity is a finite difference through two further leaf solves that carry no rows.
- The replacement gate's logistic has two branches by the sign of its exponent. A healthy plant sits tens of gate widths above the centre, where `1/(1 + e^x)` overflows: the value survives, its derivative is `inf/inf`, and a NaN partial times `a_pl0 = 0` is still NaN.
- The tracked return is seeded at birth only when `a_pl0 > 0`; it reaches no rate otherwise, and seeding costs a leaf solve per birth.
- `set_initial_states` seats a newborn (it rewrites the seed's height), so it cannot seed reserves at an arbitrary size. Tests use `tf24_seat()` in `tests/testthat/helper-tf24-seat.R`.

### 17b. TF24w's cached column potentials

- The strategy object is shared by every cohort, and it caches the soil potentials of one column (`psi_soil_own_`) for the leaf solve to read through `soil_water_potential_for_leaf`. Every call that reaches a leaf solve must fill them first from the column it means: `compute_rates` from the cohort's water states, and `establishment_probability` and `set_initial_states` from the newborn's column. A new path to a leaf solve that skips this reads whichever column was filled last, with no error unless the cache is still empty.

## Adding a state

### 18. A state held at zero error loosens every run's tolerance

odelia's step rule reads the error as an RMS over every state. A state that is inert at the defaults (zero error every step) dilutes that mean, so adding one loosens the effective tolerance of every run that carries it, and default runs move inside the ODE tolerance though nothing they compute changed. For TF24's three allometry states the movement is 2e-7 on the default-run baseline and 1.1e-5 on a substepped-soil stand. Expect it, measure it, and say so in the PR; it is not a bug in the new state.
