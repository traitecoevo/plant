# Flexible allometry for TF24: leaf area and sapwood as ODE states (#516)

Design note for the epic #516, covering tasks #513 (sketch), #514 (draft implementation) and #407 (break the fixed relationship). Written before the code so the decisions are reviewable and so a second session on this branch can see what is required rather than optional.

Issues: #516 (epic) · #513 · #514 · #407 · #512 (empirical review, still open) · #517 (NSC, the precedent) · #609/#619 (storage bounds) · #620 (`storage_prod_eps`, open) · #590 (birth-date coordinate)

## Problem

TF24 inherits FF16's fixed allometry, `area_leaf(h) = (h/a_l1)^(1/a_l2)`, refreshed on every height write. This was a deliberate and productive simplification: one size dimension makes the state a scalar, makes the size-density PDE one-dimensional, and is much of why the SCM is fast. It is also what the model has been calibrated against.

What it cannot represent is the first thing a struggling plant does — thin its canopy. Three consequences hold in the current code:

- Leaf turnover is charged as a cost and the tissue is then restored by the allometry, so replacement is unconditional. A plant cannot economise by carrying less leaf.
- `dh/dt >= 0` is structural (`growth_flux >= 0` times two positive derivatives), so shrinking is not expressible at all.
- At the empty-reserve boundary the unmet part of the carbon deficit leaves the budget. NEWS records this and names "routing it into tissue loss" as deferred work — and that sink is a leaf-area state, which makes this change a completion of the NSC work rather than a new direction.

## Decisions

| Decision | Choice |
|---|---|
| Recovery | Gradual relaxation toward the preferred trajectory, not priority rebuild |
| Scope | Leaf area **and** sapwood area. Root:leaf is a separate issue |
| Gate signal | NSC relative reserves `r = S/S_max`, reusing the existing logistic |
| Density coordinate | Birth-date for TF24, via a tri-state `Control` field |
| Sapwood target | `theta * area_leaf`, written as a function. Optimality-derived targets are a separate epic |
| Delivery | One PR, numbered commits, each recording its own number movements |

## Required, not optional

**The birth-date coordinate is a prerequisite.** The height-coordinate density equation subtracts a compression term `∂g/∂h` obtained by a finite difference that perturbs height while holding every other state at its absolute value. That equals the true compression only when growth is a function of size alone. NEWS already records the reserve gate breaking this. A leaf-area state makes it qualitatively worse, because height's rate is defined *through* leaf area — so perturbing height at fixed leaf area changes growth only via self-shading, which is near zero in a uniform light environment, and the probe returns approximately the wrong thing by construction. Under the birth-date coordinate the density ODE is `-mortality`, which is valid for a state vector of any dimension.

**Sapwood must be freed with leaf area, not after it.** Mortality reads only `r = S/S_max`, with `S_max = a_st1 * mass_sapwood` and `area_sapwood = theta * area_leaf`, so `S_max` is proportional to leaf area. Free leaf area alone and shedding shrinks the denominator of `r`, raising apparent reserves and dropping mortality to its floor at zero carbon cost, while pushing `r` above 1 — the hazard flagged in #609, which the storage rate's form no longer bounds. Physically a leaf's sapwood dissolving with the leaf is also wrong: sapwood leaves the pool only by becoming heartwood, at `k_s`.

Freeing sapwood fixes all of it emergently. The canopy thins at up to `k_l = 0.457`/yr while sapwood decays only at `k_s = 0.2`/yr, so per-leaf-area maintenance cost rises as the canopy thins — the penalty that makes shedding self-limiting — and per-leaf-area hydraulic supply improves, which is real drought acclimation reaching phylloptim through `leaf_specific_conductance_max`.

**Keep `S_max` on sapwood.** Bark does not help: `area_bark = a_b1 * theta * area_leaf`, the same proportionality, so capacity on bark has exactly the shrinking problem capacity on sapwood has. Bark helps only if bark is repinned to sapwood rather than to leaf area, and then it is equivalent. Heartwood would make capacity monotone but NSC sits in living parenchyma, so sapwood is the right tissue. Once sapwood is a state, capacity declines only at `k_s`, which is enough.

## Mechanism

Two new states, `state_size()` 6 -> 8. **Each is integrated as a log departure from its preferred value, not as a level:**

- `log_area_leaf_departure` = `ln(A / A*(h))`, where `A*(h) = (h/a_l1)^(1/a_l2)` is the preferred canopy for the current height;
- `log_area_sapwood_departure` = `ln(A_s / (theta * A))`, added in the sapwood commit.

Leaf area itself is then derived, `A = A*(h) * exp(phi)`, and cached in the `competition_effect` aux the model already reads.

**Why the departure rather than the level.** Integrating `A` directly makes the pair `(h, A)` a *redundant* system: `A = A*(h)` is an invariant of the exact flow but Runge–Kutta preserves it only to truncation error, so height's rate would read a drifted `A` and exactness would be unreachable *by construction* rather than through any coding mistake. On the departure coordinate the allometric motion stays analytic and only the departure is integrated, so with plasticity off `dphi/dt` is exactly `0`, `phi` stays exactly `0`, and `exp(0) == 1` makes `A` bit-identical to today's `area_leaf(height)`.

That exactness survives the solver because odelia's error control is a **max** norm over the states (`ode_control.hpp`, `rmax = std::max(r, rmax)`): a state with zero error contributes zero, so step sizes are untouched. An RMS norm would have diluted the mean and changed the steps, and the exactness claim would have been false — worth re-checking if odelia's controller ever changes.

Three consequences worth having:

- **Zero is the right default.** `Internals` zero-initialises states and `phi = 0` *is* "on the preferred trajectory", so a freshly constructed `Individual`, a `make_initial_state()` seeding, and an export/resume are all correct with no seeding code. Integrating levels would have started every plant at `A = 0`.
- **Neither state belongs in `non_negative_states()`** — a log ratio is signed, so there is no bound to violate and nothing to declare.
- **Scale-free by construction**, which is what keeps the plasticity parameters from repeating `storage_prod_eps`'s mistake (#620) of being an absolute constant against fluxes spanning six orders.

`height` stays a state and stays monotone: its rate remains a product of non-negative factors, so height never decreases even while the canopy is shedding.

**Conditional replacement of leaf turnover.** A smooth replacement fraction `f(r)`, logistic, centred below the growth threshold `a_st2`:

- the leaf turnover charged to the budget becomes `f(r) * k_l * m_leaf`, other pools unchanged;
- leaf area declines by the remainder, `dA/dt` gaining `-(1 - f(r)) * k_l * A`.

These are two halves of one decision and conserve carbon exactly: carbon not spent is `(1-f) * k_l * lma * A`, and leaf mass lost is `lma * (1-f) * k_l * A`. **Assert this as an identity, not a tolerance.** It also bounds thinning naturally at `k_l`, because the mechanism is declining to replace what died rather than active shedding.

**Priority ladder.** Growth is already gated at `a_st2` by `G(r)`. Gating replacement at a lower centre gives a smooth ladder out of two logistics with no branching: reproduction and growth are cut first, leaf replacement next, then reserves empty and mortality rises.

**Gradual rebuild.** `growth_flux = Ppos * G` splits smoothly between extension along the preferred trajectory and closing the leaf-area gap, with the split a smooth function of the gap that vanishes when it closes, so the gap closes asymptotically with no breakpoint.

**Sapwood.** `dA_s/dt = build - k_s * A_s`, preferred value `theta * A`. The loss term stays unconditional, matching the existing behaviour that sapwood-to-heartwood conversion proceeds regardless of carbon status: sapwood is not sheddable, so its plasticity enters only through the allocation of new growth, as a slow ratchet.

## Invariants a future edit could break

- **Height's monotonicity is by construction, and nothing can catch it if that breaks.** `dh/dt` is a product of non-negative factors, which is why height never decreases. That has to stay true structurally, because the solver cannot check it: a step is rejected either by a `DomainError` thrown from inside `compute_rates` (as the negative-storage check does) or by `ode_state_valid()` refusing the state the step landed on — and the latter is a predicate over the flattened state vector alone, so it can express "height >= 0" but cannot see the previous step and so cannot express "height is non-decreasing". `non_negative_states()` will not help here. If the product structure ever has to go, the replacement is an explicit `stop_domain` on a negative computed rate, not a state bound.
- **A new call into phylloptim must translate its exceptions.** `phylloptim::util::infeasible_error` and `odelia::util::DomainError` are siblings, both deriving from `std::runtime_error`, so odelia's handler cannot see phylloptim's throw: it escapes the step and kills the solve having taken none (#608 measured this). `src/tf24_strategy.cpp` and `src/tf24f_strategy.cpp` catch and re-throw as `stop_domain` for that reason. Any new leaf-solving path needs the same translation.
- **Every new expression must be smooth.** No `min`, `max` or `?:` in the carbon-to-growth or carbon-to-mortality map. The house pattern is at `src/tf24_strategy.cpp:201-297`: shape the flow so bounds follow from the form rather than clamping the read. `test-strategy-tf24.R:427-520` asserts a derivative-difference *convergence rate* rather than a value at one step width, precisely because a clamp passes bound assertions while leaving the rate kinked; `:522-560` asserts `|lambda|/own_rate < 10`. New gates need both.
- **The exactness anchor.** With plasticity off the model should be bit-identical to today. This is a diagnostic goal subordinate to the science — if a better formulation breaks it, take the formulation and re-bless — but it must be *checked early*, while the diff is small enough to attribute a discrepancy. Do not discover a lost exactness late and relax a tolerance to hide it.
- **`competition_effect` is load-bearing, not diagnostic.** `compute_rates` never calls `area_leaf`; it reads the cached aux written only by `update_dependent_aux`. An allometry change not routed through that function passes every test that calls `area_leaf` directly while the solver keeps the old behaviour.
- **`set_ode_state` refreshes aux per slot, in slot order.** A height-triggered write at slot 0 would read a stale leaf area still at slot 6. Do not simply add a second branch — write the aux after all states are in, or derive it from whichever state is authoritative.
- **The derivative family is hand-differentiated and must move in lockstep.** `dmass_sapwood_darea_leaf` carries an `(a_l2 + 1)` factor that is the product rule against `H(A)`; `dmass_bark_darea_leaf` inherits it; `darea_leaf_dmass_live` is their reciprocal sum. `a_l1`/`a_l2` appear on exactly four lines of TF24 code.
- **Canopy shape stays `Q(z/height)`**, so a smaller leaf area at unchanged height is a *sparser* crown of the same depth, not a shallower one. That is the intended reading of thinning, but it only becomes visible once the two decouple.
- **`set_initial_states` is non-virtual.** TF24f must keep chaining to the base or the leaf-area seed is silently dropped.
- **Height crossings become routine.** The birth-date integral is safe, but `describe_nodes_near`'s "node heights are NOT decreasing" message becomes a false accusation, `r_set_heights()` still stops, and `scale_node_densities` selects by height.
- **The version drift-guard will not catch a missing bump.** `test-model-version.R` snapshots only `pars`/`control` defaults, so adding a state without a parameter passes silently. Bumping `scientific_version` is reviewer discipline.

## Not in this work

- **Optimality-derived allometric targets** — own epic under #516. The criterion is stationarity of the plant's own *growth rate*, not leaf profit: against profit sapwood is always worth building, because nothing charges the plant for the carbon the sapwood consumed, whereas putting growth in the objective internalises the trade-off. The signal must be long-run, and sapwood's own inertia may supply that for free, since `dA_s/dt = k(target - A_s)` makes `A_s` an exponentially-weighted moving average of its target. Needs `∂profit/∂kmax` from phylloptim, whose gradient outputs are currently disjoint from the `profit_` plant bills. The same criterion is how preferred *leaf* area would be derived, which is the larger prize.
- **Root:leaf plasticity** — roots are sheddable like leaves, so it reuses this machinery with no new mechanism.
- **`storage_prod_eps` rescaling** (#620) — inherited, and this work is sensitive to it near the compensation point. Its own commit and its own number check; do not silently absorb it.

### Deciding the objective by invasion analysis

*Which* growth rate the optimality criterion maximises — `dh/dt`, `dA/dt`, or the relative `(dA/dt)/A` — should not be argued from first principles. All three are **proxies for fitness**, and a proxy is worth only what its correlation with lifetime offspring production is in the environment the plant actually experiences. plant already computes that currency, so this is an empirical question inside the model rather than a modelling preference.

The sharper form is not "which proxy wins" but **which proxy attains the ESS**. Two separable questions, in this order:

1. **What sapwood:leaf ratio is evolutionarily stable?** Treat `theta` as a trait and run the standard resident–mutant invasion analysis. This needs **no proximate rule at all** — evolution supplies the set-point — and it needs no new code: `theta` is already a `TF24_Pars` field and `add_strategies()` accepts it as a trait (verified, two distinct values across two strategies), and `run_mutant()` works for TF24 as of #643.

   Run this **first**, because it prices the rest of the epic. If the fitness landscape in `theta` is close to flat, an optimality-derived target is not worth building whatever it costs. If it is sharply peaked, the ESS is the number every candidate proxy then has to reproduce.

2. **Does a plastic rule beat the best fixed ratio, and which proxy does it best?** Only here is a rule needed. The test is whether a plastic strategy invades a resident sitting at the ESS ratio from (1). If it cannot, the plasticity is not earning its keep and the fixed ratio stands.

**The trap that would silently drain this of power.** Under constant conditions all three proxies are monotone in the same carbon surplus and will very likely agree, possibly exactly. The differences live in the *fluctuating* case, which is the whole motivation for plasticity — so the comparison has to run under drought or deep shade (the scenario gateway's regimes), not a benign constant environment. This is the same shape as the trap already recorded for shading models: a stand-alone individual sees uniform light, so every shading model agrees exactly and the test has no power. **A null result in a constant environment means the experiment was posed wrongly, not that the choice does not matter.**

**Scope of what it establishes.** Invasion analysis answers "given TF24's physiology, which rule is evolutionarily stable". That is well posed, but the answer is a property of the model, not evidence about real plants — using the model to select its own mechanism is legitimate only if the conclusion is stated at that scope. The reality check is #512: observed Huber values and how they shift with light and water. Both are needed and they answer different questions.

**Cost is the binding constraint.** TF24 is expensive per step and an ESS search over one trait is many resident-plus-mutant runs under a fluctuating driver. Expect this to want HPC rather than a laptop.

**Evaluate mutants one at a time, not in a batch.** Measured on a wet one-species TF24 stand (`mpl = 20`, default traits, `theta` as the trait), comparing the resident's own R0 against the resident slot the mutant run recomputes:

| mutants in the call | resident slot | relative shift |
|---|---|---|
| 1 (a twin of the resident) | 5.7636882085232299 | 5.5e-07 |
| 9 (spanning 0.25x to 4x) | 5.6314597794715402 | **2.3e-02** |

In both cases the twin mutant agrees with the resident slot *exactly*, so this is not mutants leaking into the competitive environment — that part works. The likely cause is that mutants are integrated in the **same ODE system** as the resident, so their error enters the shared max-norm step-size control and changes the accepted steps, and with them the resident's own trajectory. That is a hypothesis: the batch here varied both the count and the extremeness of the mutants (some sit at R0 ~ 1e-18), and those were not separated.

Either way the practical rule is the same, and it is not a small correction: a batched landscape carries a per-cent error that depends on **which other mutants happened to be in the batch**, which is exactly the kind of contamination an ESS search cannot tolerate, since it compares nearby strategies. One mutant per call costs one run each and buys ~1e-7.

### What stage 1 already found

Run on a wet (1.5 m/yr) and an arid (0.4 m/yr) constant stand, `mpl = 20`, default traits, `theta` swept over 0.25x to 4x of its default 1/4669.

**The landscape is nowhere near flat, so the epic is worth building.** On the wet stand R0 spans about **1e20** across the sweep. Fitness is steeply asymmetric: *raising* theta above the default collapses R0 by 13 orders of magnitude and more, while lowering it toward ~0.7x is strongly favoured — the peak mutant reaches R0 ~ 200 against the resident's 5.76. Sapwood is expensive (its mass scales with theta *and* with height, so respiration and turnover scale with it) while the hydraulic benefit saturates, and the model prices that very sharply.

**Two things it does not establish.** The peak is one step of an iteration, not an ESS: for that the resident has to be moved to the peak and the scan repeated until the peak sits at the resident, and the resident should be at demographic equilibrium (R0 = 1) rather than the 5.76 it sits at here. And the arid stand is **not an ESS analysis at all** — every strategy on it has R0 between 1e-11 and 2e-9, so nothing persists and the ranking is among strategies that all go extinct. Its apparent preference for *less* sapwood under drought is therefore not evidence of a drought response; that comparison needs a stand where something survives.

**The extremely seasonal case does not complete**, failing on the negative-storage domain boundary at t ≈ 2.94 with the sub-step already at its minimum. That is the pre-existing #550/#609 failure mode under an extreme driver, not something this work introduced.

**Worth its own look:** TF24's default theta is far from its own model's fitness optimum, by ~1.4x in the trait and ~35x in fitness on a wet stand. Since theta is a *calibrated* value, the more likely reading is that the model over-prices sapwood rather than that the parameter is wrong — but either way it is a finding about TF24 rather than about this epic, and it deserves separating from it.

## Where this stands after the first behavioural validation

The mechanism is in and exact when off (`a_pl0 = 0`), and it thins a canopy in
the intended direction. Three things measured on it say the default should
**not** be turned on yet, and each is a modelling decision rather than a bug.

**Measured, single plant at h = 5 m with reserves held empty, `a_pl0 = 1`:**

| departure `phi` | `dphi/dt` | `dpsi/dt` | `dh/dt` | production per leaf area |
|---|---|---|---|---|
| 0.00 | −0.422 | +0.237 | 0.632 | 1.765 |
| −0.50 | −0.290 | +0.105 | 0.383 | 1.894 |
| −1.00 | −0.272 | +0.087 | 0.242 | 1.989 |

Thinning throughout, the Huber value rising, height still climbing but more
slowly, and production per leaf area rising 13 per cent as the canopy thins --
the benefit the mechanism exists to deliver.

**1. Rebuilding has to be gated on the same signal as replacement.** The first
version made the rebuild share a function of the gap alone, so a starving plant
withheld replacement and spent its growth flux rebuilding at the same time.
Measured, it rebuilt at ~1.5/yr against shedding's 0.42/yr: the canopy thinned
5 per cent and stalled, and the Huber value *fell*. Fixed by multiplying the
share by `replacement`, so a plant rebuilds only to the degree it is willing to
maintain. This one is settled and done, and is why the factor is load-bearing.

**2. The response is far too small at the default parameters.** On a strongly
seasonal stand (mean rain 0.7 and 1.1 m/yr, amplitude 0.9 of mean) the canopy
thins by only **1.6 to 3.0 per cent** (`phi` reaching −0.030 and −0.016). That
is not "thinning out the canopy" in any meaningful sense. The cause is that the
NSC gate keeps reserves well above `a_pl1 = 0.05` except in permanent deficit,
so the replacement gate almost never opens. Either `a_pl1` has to sit much
closer to the growth gate `a_st2 = 0.10`, or the gate should read something
other than the reserve fraction. Undecided.

**3. The sapwood departure has no restoring force, and drifts the wrong way.**
Its only terms are `+withheld*(k_l - k_s)` and `−rebuild_rel`. Nothing relaxes
it toward zero, and because rebuilding buys leaf area *without* sapwood, and
rebuilding is fast, the second term dominates: on both seasonal stands `psi`
ended in [−0.10, 0], i.e. plants finished hydraulically **under**-built rather
than over-built. Rebuilding a canopy ought to buy the sapwood that supplies it.

Fixing it means splitting the growth flux three ways -- extension, leaf
rebuild, sapwood rebuild -- with the shares summing to at most one. That is the
three-degrees-of-freedom allocation problem against a single budget constraint,
and it needs a decision rather than a default.

**4. A stiffness hazard, identified and not fixed.** `rebuild_rel` is an
absolute carbon rate divided by leaf area, which is correct -- it is what turns
a carbon flux into a relative rate, and the absolute `dA/dt` stays bounded by
the carbon. But in the departure coordinate it scales as `exp(-phi)`, so it
grows without bound as the canopy empties. This is structurally the same
fast-attracting-boundary shape that #609's analysis identified for `log S`, and
the same remedy applies if it bites: a basis whose rate vanishes at the
boundary. It has not caused a failure in any run here, and a plant deep enough
in deficit to reach it is dying, but a cohort in that state is carried by the
SCM and can slow or reject steps.

**On fitness.** Turning plasticity on lowered R0 in every stand tested: −0.07
per cent on a wet stand where it barely engages, and factors of 0.76 and 0.83
on the seasonal ones. All the seasonal stands have R0 far below 1, so nothing
persists on them and the ranking is among strategies that all die out -- the
same caveat as the theta scan, and the reason these numbers do not settle
whether the mechanism helps. A stand where something survives, and where the
deficit is transient rather than permanent, is what would.

## Open

- Which growth rate the optimality criterion should maximise is **not** a modelling choice to be argued — see "Deciding the objective by invasion analysis" above. It is settled by experiment, and the cheap half of that experiment can run before any of this code exists.

## Settled

- **The new gate centres and the rebuild rate are `TF24_Pars` fields, and therefore settable as traits.** Nothing further is needed for that: plant's trait mechanism reaches nested `pars`, so a parameter that lives there can be varied through `trait_matrix()`/`add_strategies()` immediately — verified on `theta`, which takes two distinct values across two strategies with no hyperpar involvement. A **hyperpar** entry is a different thing and is *not* wanted yet: that machinery is for a parameter *derived* from another trait (as `k_l` is from `lma`), and there is no trade-off to encode until #512 says what one would look like. Defaults reproduce the current model.
