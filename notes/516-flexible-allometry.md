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

Two new states, `state_size()` 6 -> 8: `area_leaf` and `area_sapwood`, both added to `non_negative_states()`. Positivity is a property of the form rather than a guard — every loss term is proportional to the state itself, so each decays exponentially and cannot cross zero. `height` stays a state and stays monotone: its rate remains a product of non-negative factors, so height never decreases even while the canopy is shedding.

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

## Open

- Which growth rate the optimality criterion should maximise is **not** a modelling choice to be argued — see "Deciding the objective by invasion analysis" above. It is settled by experiment, and the cheap half of that experiment can run before any of this code exists.

## Settled

- **The new gate centres and the rebuild rate are `TF24_Pars` fields, and therefore settable as traits.** Nothing further is needed for that: plant's trait mechanism reaches nested `pars`, so a parameter that lives there can be varied through `trait_matrix()`/`add_strategies()` immediately — verified on `theta`, which takes two distinct values across two strategies with no hyperpar involvement. A **hyperpar** entry is a different thing and is *not* wanted yet: that machinery is for a parameter *derived* from another trait (as `k_l` is from `lma`), and there is no trade-off to encode until #512 says what one would look like. Defaults reproduce the current model.
