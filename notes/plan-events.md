# A general discrete-events mechanism (#628)

**Status: in progress.** This is the working plan for [#628](https://github.com/traitecoevo/plant/issues/628), the event-structure subtask of [#522](https://github.com/traitecoevo/plant/issues/522). The constraints in it bind anyone editing the schedule, the SCM run loop, or the TF24 soil state while the work is open.

Two other issues consume what is built here: [#627](https://github.com/traitecoevo/plant/issues/627) wants removal that can target size classes, which is the `harvest` action below; [#601](https://github.com/traitecoevo/plant/issues/601) wants the same queue for finite births and deaths when the stochastic and deterministic solvers are unified.

## Context

`plant` has exactly one kind of discrete thing that happens during a run: a node (cohort) introduction. `SCM::run_next_impl` ([scm.h:213-291](../inst/include/plant/scm.h#L213-L291)) pops the introductions due at the current time, mutates the patch, re-reads the ODE state into the solver, and integrates to the next introduction. That *is* the standard "stop the integrator, apply a jump, resume" pattern — it is just hard-wired to one action.

Several strands of science want the same machinery for other jumps:

- **Rainfall pulses** (#522, the motivating case). TF24 takes rainfall as a smooth cubic spline, which smears out exactly the variability that matters ecologically — large pulses.
- **Births with counts, and deaths** (#601). Unifying the stochastic and deterministic solvers needs both.
- **Management and partial disturbance** — harvest, thinning, damage that does not clear the patch (#627).
- **Climate extremes** — heat, in the case raised by @itowers1 on #522.

The deliverable is **a generic interface plus a deliberately crude worked example of each event type**, so the team refines concrete cases rather than arguing an abstraction in the void. The science in each action is placeholder; the interface is the product.

## Decisions

1. **Every type gets a real but crude implementation**, each a small replaceable unit.
2. **Events are instantaneous to the outer solver** — see below, which is the load-bearing idea.
3. **Node introduction becomes just another event type.** One queue, one representation.
4. **The queue lives in plant, owned by the runner**, as `NodeSchedule` is today. No odelia change.
5. **Events are a `run_scm()` argument**, not a field on `Parameters`.

## What "instantaneous" means, and why it already covers heatwaves

An event is instantaneous **from the outer solver's point of view**: patch time does not advance across it, and odelia refuses to move the clock during one anyway ([ode_solver_internal.hpp:438-452](../../odelia/inst/include/odelia/ode_solver_internal.hpp#L438-L452) — `set_time` will not move backwards or sideways).

But an action is free to compute its jump **however it likes**. It may apply a closed-form rule, or it may internally sub-integrate its own fast sub-model over the event's nominal real-world duration — a fortnight at half-hourly steps — with demography frozen throughout. Both are the same thing to the solver: state in, state out, clock unmoved.

This is what makes the heatwave case fit **with no widening of the solver contract at all**. A `ClimateExtreme` event carries a nominal duration, runs its damage sub-model internally, and hands back a jump. The event interface does not need a duration concept; the *action* does, and that is a private matter for that action.

The accepted approximation this buys, stated plainly: patch age does not advance during an event's nominal duration. Ten fourteen-day heatwaves in a 105.32 yr run leave 0.38 yr — **0.36%** — of real time unaccounted.

## The bit-identity rule

Node introductions are *already* stop-and-jump events, so the honest statement is not "runs with no events" — it is:

> **A run is bit-identical when its set of stop times and actions is unchanged.** The migration below changes neither, so it must be bit-identical. Adding any new event adds a stop time, which changes the adaptive step sequence and therefore legitimately changes the trajectory at solver-tolerance level.

Both halves matter. The first is the guardrail on the migration commit. The second must be said in the PR so a reviewer does not read an expected change as a regression.

## Design

### 1. The event type: tagged POD, templated dispatch

`Patch<T,E>` is templated over four `(Strategy, Environment)` pairs. A polymorphic `Action` with a virtual `apply(Patch<T,E>&)` would force the action classes — and therefore the queue — to be templated too, turning today's single R-facing `NodeSchedule` binding ([RcppR6_classes.yml:268-309](../inst/RcppR6_classes.yml#L268-L309)) into four.

So: **keep the queue plain data; put the dispatch in the templated layer.**

```cpp
enum class EventType   { ResourcePulse, ClimateExtreme, Harvest, NodeIntroduction };
enum class EventTarget { Patch, Environment, Species };

class ScheduleEvent {
public:
  EventType type;
  EventTarget target;            // what it acts on
  size_t target_index;           // which species, or which resource
  std::vector<double> params;    // per-type payload, incl. nominal duration where relevant
  std::vector<double> times;     // [t_intro, ...extra ode times..., t_end] -- semantics UNCHANGED
  double time_introduction() const { return times.front(); }
  double time_end()          const { return times.back(); }
};
```

Each event says when it happens, what kind of thing it is, what it acts on, and the values it needs — the four things #628 asks for. Type and target are separate because the same action serves different scopes: harvesting a whole patch and harvesting one species run the same code over a different set of nodes.

**The names are deliberately taxa- and model-agnostic.** This layer is shared by every strategy and environment, so nothing in it should be specific to plants, to water or to temperature: a resource pulse is water in TF24 and could be anything countable in a size-structured animal model; a climate extreme is heat in one model and could be cold, salinity or hypoxia in another. Model-specific vocabulary belongs in the model, where it is accurate — `TF24_Environment::add_water_pulse()` and R's `rainfall_pulse()` are the same action under the name that reads correctly there. The pulse rides plant's *existing* generic abstraction, `Environment::n_resources()`, so `target_index` names the resource just as it names the species elsewhere.

**There is deliberately no per-cohort target.** A cohort has no stable address across a run — nodes are appended and never removed, and `refine_schedule()` changes how many exist — so "cohort 7" in a schedule written before the run is not well defined. Selecting particular cohorts is expressed as a predicate on their state (a size band) in the action's parameters, which is both well defined and what #627 actually asks for.

`target_index` stays a named field rather than folding into `params`, so the run loop's `ret.push_back(...)` and `collect_competition_errors(added)` survive untouched.

Dispatch is `Patch<T,E>::apply_event(const ScheduleEvent&)` switching on `type`. Closed extension is the right trade for five in-tree types, and it keeps every object copyable — which matters because `next_event()` returns by value and `odelia::ode::Solver` owns the system by value.

### 2. Events are a `run_scm()` argument

The tempting justification for hanging events off `Parameters` is that it keeps a run self-describing. **That justification does not survive contact with the code**: `Parameters` does not carry the `Environment`, and the environment is where all the continuous drivers live ([R/scenario_eval.R:146-194](../R/scenario_eval.R#L146-L194) is the canonical example — rainfall goes on `env`, not `p`). A resumed or re-run simulation already requires `env` to be supplied separately, so adding events to `Parameters` would buy a self-description that is already incomplete.

So events become a first-class argument:

```r
run_scm(p, env = NULL, ctrl = control(), events = NULL,
        refine_schedule = FALSE, collect = FALSE, use_ode_times = FALSE)
```

and the C++ constructor gains it: `SCM(parameters, environment, events, control)`.

**Node introductions live in the same object.** `events = NULL` means exactly `events(node_introductions(p))` — the defaults `Parameters::validate()` already generates ([parameters.h:118-129](../inst/include/plant/parameters.h#L118-L129)). `p$node_schedule_times` is retained as a **deprecated view** onto the introduction subset: existing scripts keep running unchanged and silently, documented in `NEWS.md` as scheduled for removal. It must not warn — [tests/testthat/setup.R:4](../tests/testthat/setup.R#L4) makes warnings errors across the whole suite.

`refine_schedule()` ([scm.h:334-372](../inst/include/plant/scm.h#L334-L372)) rewrites only the introduction subset and leaves other events untouched. **This is a concrete API requirement**: `Schedule::set_times()` today clears and rebuilds everything ([node_schedule.cpp:59-68](../src/node_schedule.cpp#L59-L68)) and must become merge-preserving.

### 3. Ordering at a shared time

Two events at one instant need a deterministic, documented order. Today `add_time` ([node_schedule.cpp:273-283](../src/node_schedule.cpp#L273-L283)) inserts equal-time events before existing ones, which is an accident, not a decision. Sort by `(time, type_rank)`:

1. **Environment events** (a resource pulse) — set the external conditions.
2. **Demographic removals** (harvest, partial disturbance).
3. **Node introductions** last, so a newborn's initial conditions ([node.h:183-213](../inst/include/plant/node.h#L183-L213)) are computed against the post-event environment.

With introductions only, this ordering is inert — so it costs nothing in bit-identity.

`Patch::check_birth_dates_distinct()` ([patch.h:377-396](../inst/include/plant/patch.h#L377-L396)) is per-species, so two *different* species at one time is already fine and no new constraint arises.

### 4. The rewritten run loop

```
const double t0 = time();
Event e = schedule.next_event();

if (e.time_introduction() > t0) { ...resume-gap branch, UNCHANGED... }

std::vector<size_t> introduced;
std::vector<Event>  actions;
while (true) {
  if (!util::identical(t0, e.time_introduction())) util::stop("Start time not what was expected");
  if (e.type == EventType::NodeIntroduction) introduced.push_back(e.species_index);
  else                                       actions.push_back(e);
  schedule.pop();
  if (e.time_end() > t0 || complete()) break;
  e = schedule.next_event();
}

for (auto& a : actions) sys.apply_event(a);          // in type order
if (!introduced.empty()) sys.introduce_new_nodes(introduced);
solver.set_state_from_system();

...three integration modes on e.time_end(), UNCHANGED...
return introduced;
```

The break condition still keys off `e.time_end() > t0`, which `reset()` ([node_schedule.cpp:87-104](../src/node_schedule.cpp#L87-L104)) derives from the *next* event's introduction time. That derivation is type-blind, so a mixed queue automatically shortens the preceding leg — the behaviour we want, for free. `distribute_ode_times()` ([node_schedule.cpp:119-140](../src/node_schedule.cpp#L119-L140)) is likewise unchanged.

`complete()` is purely `remaining() == 0` ([scm.h:388-390](../inst/include/plant/scm.h#L388-L390)), so events *must* live on this queue or the loop never sees them. That is the reason for one queue rather than two.

### 5. The five actions

Harvest and climate extremes share one primitive, differing only in how a per-node survival factor is chosen:

```cpp
// applies factor phi in (0,1] to one node, keeping the two accountings in step
void scale_node_density(node_type& n, double phi) {
  n.set_log_density(n.get_log_density() + std::log(phi));
  n.individual.set_state("mortality", n.individual.state(MORTALITY_INDEX) - std::log(phi));
}
```

**Both states must move.** On the birth-date coordinate the code maintains an exact identity: `log_density_dt = -rate(MORTALITY)` ([node.h:156](../inst/include/plant/node.h#L156)) with initial conditions `mortality = -log(pr_estab)` ([node.h:194](../inst/include/plant/node.h#L194)) and `log_density = log(birth_rate*pr_estab)` ([node.h:209](../inst/include/plant/node.h#L209)) — integrating gives **`log_density(t) = log(birth_rate) - mortality(t)`**. Meanwhile `survival_individual = exp(-state(MORTALITY))` weights fecundity independently ([node.h:162,171-173](../inst/include/plant/node.h#L162-L173)). Move only `log_density` and fecundity is not discounted; move only `mortality` and the standing density is not reduced. That identity is directly assertable and becomes a test.

| # | Event | Params | Crude rule | The decision to refine |
|---|---|---|---|---|
| 1 | **NodeIntroduction** | — | `Patch::introduce_new_nodes()`, unchanged ([patch.h:780-791](../inst/include/plant/patch.h#L780-L791)) | nothing; must not move |
| 2 | **ResourcePulse** | `amount`, resource named by `target_index` | jump in that pool, capped at free capacity, excess shed | should it also pass the saturation-excess term? |
| 3 | **Harvest** | `fraction`, `size_min`, `size_max` | `phi = 1-f` for nodes in the size band | removal vs partial biomass loss |
| 4 | **ClimateExtreme** | `intensity`, `duration`, `threshold`, `sensitivity` | internal sub-integration over `duration`, accumulating dose to `phi` | the whole thing — see below |

**Harvest, thinning and partial disturbance are one type.** They differ only in selectivity: leave the size band at its defaults and it is an across-the-board knock-down; set `size_min` and it takes everything above a size; set both and it thins one size class, which is what #627 asks for. One action, three ways of asking for it — so there is one implementation, not three.

## What the events actually did

Each applied event is recorded — time, type, target, what was requested, and what was achieved — and read back as `scm$event_log` ([#628](https://github.com/traitecoevo/plant/issues/628)).

This is not bookkeeping for its own sake. **The two are routinely different**: a pulse is capped at what the pool can hold, so the amount that lands is often less than the amount asked for; harvesting a size band removes whatever was in that band, which is not knowable in advance. A pulse reports `{accepted, shed}`; harvest and climate extremes report `{fraction_applied, nodes_affected}`. Without the log the shortfall is only inferable from an accumulator, which is no way to answer "what did this run do".

Node introductions are not logged: they are the schedule, not an intervention, and a hundred and forty of them would bury the three that matter.

The log lives on the runner rather than the patch, because the patch is copied into `history` once per step and a log that grew with the run would be copied with it every time.

**Rainfall pulse.** A new `TF24_Environment::add_water_pulse(double depth)`:

```
capacity = max(0, (theta_sat_0 - theta_0) * dz[0]);
accepted = min(depth, capacity);   excess = depth - accepted;
vars.set_state(0, theta_0 + accepted / dz[0]);
vars.state(n+0) += depth;      // sum_rainfall
vars.state(n+1) += accepted;   // sum_infiltration
vars.state(n+4) += excess;     // sum_pulse_runoff  (NEW)
psi_soil_cache_valid_ = false;
```

The cap is required, not optional: the jump is applied *outside* the integrator, so no error estimate and no step rejection protect it. PR #608's measured table shows a realistic ~13 mm event already exceeds layer 0's free capacity from a moderately wet start. A hard `min()` is fine precisely because it is outside the integrator — the "no kinks in the rates" concern applies only to continuous derivatives.

Two traps:

- **Do not use `set_soil_water_state()`** ([tf24_environment.h:506-519](../inst/include/plant/models/tf24_environment.h#L506-L519)) — it zeroes the accumulators *and* resets `initial_states`.
- The accumulators are **rate-integrated** (`vars.set_rate`, [tf24_environment.h:381-384](../inst/include/plant/models/tf24_environment.h#L381-L384)), so a pulse's contribution must be a direct **state** increment, not a rate.

On whether the pulse should also pass the empirical infiltration term at [tf24_environment.h:336-338](../inst/include/plant/models/tf24_environment.h#L336-L338): **no, for v1.** That term is a *rate*-based partition tuned for continuous forcing; applying it to an instantaneous depth double-counts with the capacity cap. Flagged as the first refinement.

**The 5th accumulator must be pulse-only, and this is not cosmetic.** odelia's controller takes `rmax` over *all* state components, so a new component with a non-zero rate can change accepted step sizes and move TF24 trajectories. Give `sum_pulse_runoff` **rate exactly 0** and state changed only by pulse events: with no pulses it stays exactly 0, contributes exactly 0 to the error norm, and TF24 stays bit-identical — no `scientific_version` bump, no re-bless. Accumulating *continuous* runoff too is a deliberate separate science change. **Verify this with `identical()`; do not assert it.**

**Temperature extreme is the weakest of the five, deliberately.** plant runs its TF24 leaf at a constant 25 °C ([tf24_environment.h:82-105](../inst/include/plant/models/tf24_environment.h#L82-L105)), so no physiological response exists to drive. [notes/penman-monteith/implementation-plan.md:428-434](penman-monteith/implementation-plan.md#L428-L434) already sets out the timescale mismatch and its two candidate resolutions. The crude version is the sub-integration described above: step a scalar damage variable over `duration` at a sub-daily step under a simple diurnal temperature profile, convert accumulated damage to `phi`, apply. It is a **worked example of the sub-integration pattern** rather than defensible thermal biology — the real response belongs in the ATLS damage state (#566), which plugs into exactly this hook.

**Harvest ordering constraint.** `Species::compute_competition` assumes decreasing heights, enforced at [species.h:684-695](../inst/include/plant/species.h#L684-L695). Any *monotone* height rule (scale by a fraction; cap at a threshold) preserves that; only a non-monotone rule breaks it, and the `compute_competition_unordered` fallback ([species.h:431-465](../inst/include/plant/species.h#L431-L465)) is height-coordinate-only. Crude harvest sidesteps the question by removing individuals rather than shortening them — biomass reduction is the refinement, and it must additionally call `invalidate_height_scan()` ([species.h:192](../inst/include/plant/species.h#L192)).

### 6. R interface

Per-type vectorised constructors combined by `events()`:

```r
ev <- events(
  node_introductions(p),                                        # the existing default schedule
  rainfall_pulse(time = c(1.5, 3.2, 7.8), depth = c(0.013, 0.005, 0.050)),
  harvest(time = 20, fraction = 0.5, size_min = 10),
  harvest(time = 30, fraction = 0.2, size_min = 2, size_max = 5),
  climate_extreme(time = 60, intensity = 45, duration = 14/365)
)
res <- run_scm(p, env = env, ctrl = ctrl, events = ev, collect = TRUE)
```

Plus pulse-series generators, since the decision on #522 is that the user bakes the whole series in up front — regular, drawn, or from an observed record. [R/stochastic.R:7-62](../R/stochastic.R#L7-L62) is the existing precedent for R-side event generation, and `check_driver_interpolation()` ([R/drivers.R](../R/drivers.R)) is the diagnostic that motivates pulses in the first place.

[R/tidy_outputs.R:51-57](../R/tidy_outputs.R#L51-L57) gains `"sum_pulse_runoff"` as a 5th positional name — it must stay in lockstep with the C++ state order. [tests/testthat/test-environment-TF24.R:175](../tests/testthat/test-environment-TF24.R#L175) currently *derives* runoff as `sum_rainfall - sum_infiltration`; that stays, since the new accumulator is pulse-only and measures a different thing.

## Commit sequence

1. This plan.
2. **Generalise the queue; migrate node introduction onto it.** `EventType` tag, merge-preserving `set_times()`, rewritten `run_next_impl`. **Zero behaviour change.** The risky commit, isolated.
3. **Plumb events through.** `Events` wire format, `SCM` constructor argument, R `events` argument, ordering rule.
4. **Rainfall pulse.** `aux_num` 4 to 5 (pulse-only, rate 0), `add_water_pulse`, the action, `tidy_outputs`, water-balance test.
5. **Targets, the event log, and the removal actions.** `EventTarget`, `scm$event_log`, the shared density primitive, harvest and climate extremes.
6. **Docs.** Series generators, `NEWS.md`, a vignette section.
7. **#505 rename** `NodeSchedule` to `Schedule`, pure mechanical, kept out of the risky diffs.

## Things found along the way

- **`r_set_max_time()` read `events.back()` on an empty list**, which is undefined behaviour and the *normal* path — both `make_node_schedule()` and `node_schedule_default()` set `max_time` before adding any times. It had always read harmless garbage; changing the event's member layout turned it into a segfault. Fixed, with a test.
- **A pulse suppresses the continuous infiltration that follows it.** Wetting layer 0 raises the saturation-excess term at [tf24_environment.h:336-338](../inst/include/plant/models/tf24_environment.h#L336-L338), so a run gains strictly less infiltration than the pulse delivered. Real behaviour, and worth knowing before reading a water budget.
- **Pulsed columns end *drier* than unpulsed ones.** `K(theta)` rises as `theta^16.14`, so the extra water drains fast and the wetter interval costs more drainage than it stores. Final storage is therefore the wrong thing to test a pulse with; the cumulative fluxes are the right thing.
- **Positional flux slots are fragile, as predicted.** Adding the fifth accumulator silently changed what `tail(patch$ode_rates, 1)` means — a stochastic test was reading root uptake that way and started reading a pinned zero.

## Verification

**Bit-identity** (commit 2, and 3-5 for runs whose stop times are unchanged):

```r
for (type in c("FF16", "K93", "TF24")) {
  res <- run_scm(p_for(type), collect = TRUE)
  stopifnot(identical(res, readRDS(baseline_path(type))))   # raw arrays, not by eye
}
```

Baselines are captured on `develop` **before** commit 2. `git stash` cannot revert your own commits, so a later "before" measurement needs a separate checkout.

**Regression exposure to re-check by hand:** [test-patch.R:52-57](../tests/testthat/test-patch.R#L52-L57) and [:88-91](../tests/testthat/test-patch.R#L88-L91) (literal FF16 ODE vectors; TF24 soil inits `rep(0,4)` to `rep(0,5)`), [test-scm.R:81-85](../tests/testthat/test-scm.R#L81-L85) and [:126-152](../tests/testthat/test-scm.R#L126-L152) (`expect_identical` determinism), [test-schedule-build.R:29-30](../tests/testthat/test-schedule-build.R#L29-L30) (literal lengths 141/148), [test-environment-TF24.R:64,235,242](../tests/testthat/test-environment-TF24.R#L64).

**New tests:** queue merge and ordering units; per-action units; the `log_density(t) = log(birth_rate) - mortality(t)` identity across a partial disturbance; pulse water balance (`dtheta*dz + runoff == depth`) and saturation (a pulse into a nearly-full layer 0 routes almost everything to runoff); a sub-integration event that leaves the outer clock unmoved; an end-to-end run using all five types. [test-initial-state.R](../tests/testthat/test-initial-state.R) is the best structural template.

**Scientific versioning** ([agents.md](../agents.md)): no bump expected — this adds a capability inert when unused. If the TF24 `identical()` check fails, that is a real output change needing a TF24 v2 to v3 bump plus `make bless-scenarios`. Do not bless to make a test pass.

## Non-goals, recorded deliberately

- **Runtime insertion.** Designed for (a single-event `insert()` that keeps `add_time`'s currently-discarded iterator hint), left unexposed. #601 needs it for stochastic deaths.
- **State-triggered events** (integrate until a hazard crosses a threshold). Needs rootfinding, which odelia does not have. #601's death design decides whether it is ever needed.
- **The stochastic tower.** `StochasticPatchRunner` pops only one event per time ([stochastic_patch_runner.h:84-87](../inst/include/plant/stochastic_patch_runner.h#L84-L87)); migrating it is #601.
- **Disturbance as a patch-clearing event.** `Disturbance_Regime` is a pure survival *weighting* ([patch.h:745](../inst/include/plant/patch.h#L745), [:784-785](../inst/include/plant/patch.h#L784-L785)) that never touches state — genuinely a different mechanism from a partial-disturbance event, not a duplicate.

## Two hazards

1. **Step-size inheritance.** `step_size_last` survives `set_state_from_system()`, so a leg starting after a discrete change begins at whatever step the previous quiet leg grew to — up to `h_max = 5`. PR #608 section 3 measured this as the trigger for TF24 blow-ups, and **pulses create far more leg boundaries than introductions do.** The single thing most likely to bite. Decide whether an event resets the step size.
2. **odelia's domain hooks are not adopted.** odelia 0.3.1 has `ode_state_valid()` and `util::DomainError`, which reject a step that leaves a declared domain (odelia #55/#56). plant pins **0.2.1** ([DESCRIPTION](../DESCRIPTION)) and implements neither. A capped pulse still parks layer 0 at exactly theta_sat, where `K'(theta)` is about 1.5e4 yr^-1. Worth bumping and adopting before the pulse commit lands.

Both are items 2-4 of PR #608's recommendation, whose design doc is now in the repo at [`plan-tf24-soil-redistribution.md`](plan-tf24-soil-redistribution.md) — read it for the measurements behind the capacity cap and the step-size hazard.

**Hazard 2 is now addressed**, in part: plant is on odelia 0.3.1 (#633) and declares `ode_state_valid()`, and `solve_leaf()` translates phylloptim's `infeasible_error` into a rejectable `DomainError`. Hazard 1, the inherited step size, is untouched.
