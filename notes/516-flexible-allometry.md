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

**Sapwood — and a correction to how this note first described it.** It said sapwood is *not sheddable*, so its plasticity could only enter through the allocation of new growth. That is wrong, and the implementation is better for it: sapwood is continuously lost to heartwood at `k_s`, and that loss is currently replaced *implicitly* by the allometry, exactly as leaf turnover is. So the same "decline to replace" mechanism applies, and one uniform gate covers both pools rather than two special cases.

The asymmetry is then a consequence rather than an assumption. Because the two departures differ only by their rate constants,

    dpsi/dt = withheld * (k_l - k_s) - rebuild

so withholding both raises sapwood per leaf area at up to `k_l - k_s` = 0.457 − 0.2 = 0.257/yr. The Huber value rises under stress because leaves turn over faster than sapwood does, not because anything says sapwood is special. Heartwood accumulation stays ungated either way: the tissue becomes heartwood whether or not it is replaced.

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

**3. The sapwood departure had no restoring force. FIXED.** Its only terms were
`+withheld*(k_l - k_s)` and `−rebuild_rel`, neither of which returns it to zero,
and because rebuilding bought leaf area *without* sapwood the second dominated:
on both seasonal stands `psi` ended in [−0.10, 0], leaving plants hydraulically
**under**-built. Two changes, and neither needed the three-way flux split that
first looked unavoidable:

- **Rebuilding buys the whole package** -- leaf, fine root, and the sapwood and
  bark that supply it -- at the ratio the plant *prefers* rather than the one it
  has, so closing a canopy gap no longer dilutes the Huber value.
- **Sapwood renews a scaled share**, `replacement * exp(-psi / a_pl2)`, so a
  stem carrying more conducting area than its canopy needs declines to renew the
  excess even when carbon is ample. That is the restoring force.

**`psi >= 0` is now invariant, by the form rather than by a guard.** At `psi = 0`
the two replacement fractions coincide, so `dpsi/dt = (1 - replacement) *
(k_l - k_s) >= 0` because `k_l > k_s`: the boundary flow points inward, the "too
little sapwood" branch is unreachable, and the scaling factor stays in (0, 1] so
no over-replacement can arise to cost unbounded carbon. ⚠️ **If `k_s` were ever
raised above `k_l` that argument reverses.**

Measured after the fix, on the same seasonal stands: `psi` in **[0, +0.040]** and
**[0, +0.043]** -- the right sign -- and `dpsi/dt` holds near +0.21 across the
gap range instead of decaying to +0.087. Leaf departure reaches −0.052 against
sapwood's +0.043, so leaf is the larger of the two, with the ratio set by
`k_l / (k_l - k_s)` ~ 1.8 under pure withholding: the speeds are governed by the
turnover constants, as they should be.

**What the fix does NOT do: bound either departure under a permanent deficit.**
With reserves empty, `withheld * k_l` (0.42/yr) outruns `withheld_sapwood * k_s`
(0.2/yr) whatever `psi` is, so `psi` keeps rising -- measured at +0.22/yr even at
`psi = 6`. This is arithmetic rather than a defect: a plant thinning its canopy
at 46 per cent a year while losing conducting area at 20 per cent really does
raise its Huber value without limit. But it means a plant in *permanent* deficit
has both departures diverging linearly, `phi` toward −inf and `psi` toward +inf.
What removes such a plant is mortality, not the allometry. Whether that is
acceptable, or whether the departures need a bound, is the live question -- and
it is the same question as (4) below.

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

## Second exploration: why thinning does not buy survival (2026-09-08)

Five measurements, and they move the problem out of the allometry entirely.

### 1. The sheddable costs are large. RETRACTION.

An earlier note here and in the demo said sapwood respiration dominates a drought
deficit, so thinning could not reach it. **That is wrong.** Leaf, fine root and
bark all scale with leaf area — only sapwood is sticky. At 10 m:

| leaf | bark | sapwood | root | reachable by thinning |
|---|---|---|---|---|
| 49.6% | 5.4% | 23.5% | 21.5% | **76.5%** |

(Roots follow leaf area *exactly*, `m_r = a_r1 * A`, so they shrink instantly with
the canopy rather than decaying at `k_r`. Same simplification as bark. Neither has
a state.)

### 2. But thinning cuts income just as fast

Assimilation scales with leaf area too. Net production per unit leaf area is flat
as the canopy falls from 100% to 1.8%, and slightly *worsens* as the stem it still
carries is spread over fewer leaves:

| canopy fullness | 1.00 | 0.50 | 0.20 | 0.10 | 0.05 | 0.018 |
|---|---|---|---|---|---|---|
| net production per m² leaf | −3.81 | −3.69 | −3.77 | −3.92 | −4.15 | −4.68 |
| stem per leaf | 1.0 | 1.5 | 2.5 | 3.6 | 5.4 | 9.5 |

A plant whose leaves cost more than they earn cannot fix that by having fewer
leaves. Thinning slows the burn; it never reaches balance.

### 3. There is no drought in which it changes the outcome

| soil water | 0.20 | 0.17 | 0.15 | 0.13 | 0.11 |
|---|---|---|---|---|---|
| survivorship, 4-yr drought | 0.93 | 0.93 | 0.93 | 4.9e-10 | 3.4e-10 |
| gate engaged? | no | no | no | yes | yes |

A cliff, not a gradient. Above 0.15 reserves stay full and the gate never opens;
below it the plant dies whatever it does. Thinning *early* — raising `a_pl1` from
0.05 to 0.30, which was the obvious fix — was measured across this sweep and buys
a factor of **1.0 to 1.13**, while costing height.

### 4. The root cause is the reserve pool, not the gate

`dS/dt = charge*(1-r) - drain*r` has **no resting point**: positive production
fills the pool to capacity, negative empties it, and nothing sits between. So `r`
is close to a binary readout of the *sign* of production. Mortality then reads it
through `a_dG1 * exp(-a_dG2 * r)` with `a_dG2 = 20`, which is sharp enough to turn
a near-binary input into a near-binary output spanning eight orders.

**A mechanism that makes a plant modestly cheaper to run has nothing to act on.**
The buffer does not buffer. That is an NSC-design property, not an allometry one.

### 5. The departures are unbounded, and it bites in ORDINARY conditions

Withholding gives a constant *relative* thinning rate `u * k_l`, so the canopy
decays exponentially with no floor. `k_l` is a trait, so this arrives sooner for
fast-leaved species:

| lma | k_l /yr | max thinning /yr | canopy after 2 yr | after 4 yr |
|---|---|---|---|---|
| 0.05 | 4.80 | 4.44 | 1.4e-04 | 2.0e-08 |
| 0.0825 | 2.04 | 1.88 | 0.023 | 5.4e-04 |
| 0.1979 | 0.46 | 0.42 | 0.43 | 0.185 |
| 0.40 | 0.14 | 0.13 | 0.78 | 0.60 |

⚠️ **And it is not confined to drought.** In a plain WET stand (rainfall 1.5 m/yr,
`lma = 0.0825`, self-thinning over 25 years) suppressed cohorts show exactly the
intended gradual decline at first — a 0.45 m seedling at 0.63 of its preferred
canopy — but the deeply shaded ones run away: at t = 20 a 10 m suppressed cohort
sits at canopy **0.000** with stem-per-leaf **1.3e9**. Mortality has already
removed them (survivorship 1e-34) so `R0` is unaffected (ratio 0.9992), but the
states are meaningless and numerically hazardous.

**This is a NUMERICAL failure, not a biological one — the distinction matters for
the fix.** The ordering across LMA is right and is a success: a fast-leaved species
*should* shed faster and die sooner, and `u * k_l` is the correct rate law. What
fails is that the state runs far past any meaningful floor and into arithmetic that
stops meaning anything — leaf area at 1e-8 m², production per leaf area at −9.2e7,
`stem_per_leaf` at 1.3e9 — all on the same vanishing denominator, and
`rebuild_rel` scales as `1/A` too. Those divisions sit in a solver step shared with
live cohorts.

But underneath it is a real biological gap: **nothing in the model kills a plant
for having lost its canopy.** Mortality reads reserves alone, so TF24 will carry a
plant to 2e-8 of its canopy provided its reserve fraction says otherwise. The fix
is therefore not a clamp on the departure but a mortality term that responds to
canopy state — which would also give mortality something *continuous* to read in
place of the near-binary reserve fraction.

### What the literature says, and it contradicts the current design

A global synthesis of 121 studies and 177 species finds **minimum NSC around 46
per cent of the seasonal maximum**, with depletion of total NSC "rare"
(Martínez-Vilalta et al. 2016, Ecol. Monogr., doi:10.1002/ecm.1231). Soluble
sugars in particular are "kept above some critical threshold". Beech under three
years of drought *increased* NSC for two years before starch fell, and it was the
**dead** trees that were "virtually empty" — the authors conclude that
"maintaining an active C storage function at the expense of growth was certainly
key to survival" (Chuste et al. 2019, Trees, doi:10.1007/s00468-019-01923-5).

So real plants **defend** the pool, sacrificing growth to do it, and are found
empty only when dead. TF24 does the reverse: growth is gated but storage is the
residual, so the pool drains to zero and the plant is scored dead by a fraction
that had no floor to begin with.

Also relevant: time-to-starch-depletion predicted mortality better than any
instantaneous level (Trueba et al. 2024, Ann. For. Sci.,
doi:10.1186/s13595-024-01246-7) — which is an argument for what mortality should
read.

### Where this leaves the design

Confirmed with Daniel: **the canopy thinning IS the tissue consumption that pays
the deficit.** That reframes the gate. Withholding should be sized to the
*shortfall* rather than to reserves — withhold what is needed to close the deficit,
up to what turnover makes available — so the pool acquires a resting point by
construction and `r` stops being binary. One change, addressing items 4 and 5
together, since a demand-driven withholding also stops once the deficit closes.

Open, and needing a decision:

- **What bounds the departures** when the deficit exceeds what withholding can
  close. Mortality removes the plant, but not before the states are absurd.
- **What mortality should read.** A reserve fraction cannot see a plant that has
  made itself cheaper to run. Time-to-depletion at the current burn rate can, and
  has empirical support.
- **Whether storage should outrank growth**, as the literature says it does.

## The empirical case against a reserve-gated trigger (2026-09-08)

Two sources read in full: Robinson's mulga chapter (573 *Acacia aneura*, 49 dune
sites, drone imagery 2017-2024 across a severe drought) and Manzoni et al. 2015
(*Adv. Water Resour.* 84:37-51), an optimality model of leaf phenology.

**They say the trigger should be the instantaneous carbon balance, not the size of
the reserve pool.** Four independent lines:

1. **The optimal trigger is a flux, derived not assumed.** Manzoni maximises
   cumulative season carbon gain over leaf area and gets `dG/dt = 0` -> `A_net = 0`
   (Eq. 7, p.40): shed when the marginal leaf stops paying for its own
   respiration. In water-potential units that is `psi_s ~ -2 MPa`, about day 20 of
   a dry-down — *"optimal leaf loss is a later and more drastic measure of water
   loss regulation than stomatal closure"* (p.48).
2. **It fires while the plant is still carbon-rich.** Manzoni's Fig. 4D has
   cumulative gain at its *seasonal maximum* when shedding begins. **Waiting for
   reserves to run down IS the evergreen strategy** — and that strategy runs
   `A_net` to -10 umol m^-2 s^-1 and its season gain from +150 to **-50 gC m^-2**.
   It is the loser at long drought. This is the sharpest contradiction of what we
   built.
3. **Shedding is what healthy plants do.** 85% of Robinson's trees shed, including
   the healthiest; the died and survived canopy-density distributions are
   *indistinguishable* in 2018 and only separate from 2019. Shedding is not a
   marker of imminent death.
4. **The size result runs backwards for a reserve rule.** Under reserve exhaustion
   the biggest stores should shed last and die when they shed. Robinson finds
   large trees shed *less* (*"a large canopy area apparently buffering against the
   most severe drought-driven leaf loss"*) **and** survive best — canopy area is
   the strongest mortality predictor (log-odds -0.03, p = 8.4e-5). The explanation
   offered is water access, not carbon stores.

### What the data require

- **Leaf area must reach ~zero and come back, in a plant that stays alive.**
  Manzoni's optimum is literally `L = 0` for all three shedding strategies;
  Robinson's methods treat complete de-greening followed by re-leafing as
  documented mulga behaviour, and only score a tree dead if it fails to re-green
  in a later wet year. ⚠️ **This kills the idea of canopy-dependent mortality**
  floated above — losing the canopy must not itself be lethal.
- **Decline is slow, recovery fast.** Robinson: 2-3 years down, one season back.
- **Recovery's cost is a discrete re-flush, not a rate.** Manzoni prices it at
  `gamma = LMA * f_C / Y ~ 48 gC m^-2` per unit LAI, so a full canopy is ~95
  against a 600-800 gC m^-2 season — 12-16 per cent. We charge continuous turnover
  instead and so miss the flushing cost that drives the whole result.
- **Duration matters more than intensity.** The switch is at `T_d/T ~ 0.1`;
  evergreen gain crosses zero at `T_d/T ~ 0.45`.
- **Gradual and instantaneous shedding give identical carbon gain** in Manzoni, so
  gradual is free in the optimality currency and is what the field data show.

### Consequence for this work

**The gate reads the wrong variable.** Raising `a_pl1` failed not because the
threshold was in the wrong place but because reserves are the wrong signal — they
are near-binary (see above), and they turn over *after* the decision should have
been made. Switching the gate to the instantaneous carbon balance would:

- fire early, while the plant still has reserves, which is what the empirical
  record shows and what Daniel proposed;
- leave the pool defended rather than drained, so `r` stops being binary and gains
  a resting point without any change to `dS/dt`;
- make the response continuous in a quantity that is itself continuous, so there is
  something for it to act on.

⚠️ **One structural tension it creates.** If leaf area must be able to reach zero,
the log-departure coordinate cannot represent the endpoint — `phi -> -inf` — and
`rebuild_rel ~ 1/A` diverges there. The coordinate was chosen to make the resting
model exact, and it does; but it assumes the canopy stays a fixed fraction away
from nothing. A fully-shedding species breaks that assumption, and choosing between
exactness at rest and representability at zero is now an open design question.

## When does shedding actually pay? An analytical answer (2026-09-08)

Daniel's framing: shedding is only useful if it improves the plant's net carbon
budget, and that should be designable without stand simulations. It is.

### The decomposition

At fixed height, net production splits by whether a cost scales with leaf area.
Leaf, fine root and bark all do; sapwood does not:

$$P = c\,\bar a\!\left(\tfrac{A_s}{A}\right) A \;-\; g(h)\,A \;-\; s(h, A_s)$$

**Verified bit-exactly** against `net_mass_production_dt` (relative error 1e-16 at
four size/soil combinations). `g(h)` is a constant per unit leaf area; `s` is the
sapwood burden and is independent of `A`.

### The criterion

Differentiating at **fixed sapwood** — which is what real thinning does, since
sapwood only decays at `k_s`, and which is the step my earlier probe got wrong by
letting `A_s` follow `A` down the pipe model:

$$\frac{\partial P}{\partial A}\bigg|_{A_s} = c\,\bar a\,(1-\eta) - g(h),
\qquad \eta \equiv \frac{k_{\max}}{\bar a}\frac{\partial \bar a}{\partial k_{\max}}$$

because `kmax` is proportional to `1/A` at fixed sapwood. With
`kappa = g/(c*abar)` the leaf-side cost as a fraction of leaf-side gain,

> **shedding pays iff `eta > 1 - kappa`.**

`eta` is the elasticity of per-leaf assimilation to hydraulic supply. Reading it:

- `eta = 0` (no hydraulic feedback) recovers `kappa > 1`, which is exactly
  Manzoni's `A_net = 0`. So the criterion currently implemented is the special
  case where leaf area does not affect supply.
- `eta = 1` (assimilation proportional to supply) means shedding **always** pays:
  halving leaf area leaves total assimilation unchanged and halves the cost.
- In between, shedding pays *earlier* than Manzoni's rule, by exactly `eta`.

**This closes the gap** where a tall plant died without shedding. `kmax ~ 1/h`, so
`eta` rises with height (measured 0.20 at 5 m, 0.43 at 10 m, 0.63 at 15 m, 0.84 at
20 m) — tall plants have most to gain from thinning, and the current criterion,
which sets `eta = 0`, cannot see it.

**And it rules out gating on whole-plant solvency.** At h = 15, theta = 0.20 the
plant is solvent (`P = +3.4`) and shedding still improves the budget. Solvency and
marginal benefit are different questions, and only the second says whether
shedding pays.

### Two cautions on the numbers

⚠️ **The `profit` / `shadow_cost` auxes do not faithfully report the assimilation
net production used.** The decomposition reconstructs `P` exactly at a resting
state but is 2 per cent out at a departed one (5.71 against 6.09 at h = 10,
theta = 0.20, phi = −0.05). So `abar` recovered from those auxes is unreliable off
the trajectory, and every `eta` above is an estimate rather than a measurement.
The *direct* finite difference of `P` is unambiguous and is what the sign
conclusions rest on. Worth its own look — it may be the same class of issue as
the leaf's carried state.

⚠️ **The apparent size of the benefit is an artefact of a badly-placed Huber
value.** Shedding 10 per cent of the canopy at fixed sapwood raised `P` from 3.42
to 8.33 for a healthy 15 m plant — but `P` is a small residual between large
numbers precisely because the sapwood burden is large, so a modest absolute gain
looks dramatic. The plant is sitting far from its optimal sapwood-to-leaf ratio,
which is what gives thinning such leverage.

### What this implies for the order of work

The leverage above is the **optimality-derived sapwood target** asking to be built
(the epic split out earlier). If sapwood tracked its optimum, a plant would not sit
at a Huber value this wrong, and the thinning question would be posed cleanly
instead of against a misconfigured stem.

It also depends on the height-hydraulics relation, which **#617 (epic #615)** is
replacing: it derives the height exponent from conduit widening and the Huber
profile rather than assuming resistance linear in height. TF24 currently drops
per-leaf assimilation 46-fold between 5 m and 15 m, and that number is doing most
of the work in everything above. Any `eta` measured now is measured against a
height-resistance relation that is about to change.

**So: sapwood optimality first, then re-measure `eta`, then implement the
criterion.** Implementing `eta` against the current stem model would be fitting to
a relation that #617 replaces.

## Sapwood tracking its optimum: the design (2026-09-08)

Agreed direction. Four decisions, and one of them removes the cost objection.

**1. The objective is GROWTH RATE, not net production.** Corrected by Daniel, and
it matters here in a way it does not for shedding. Growth is
`F * f_g * darea_leaf_dmass_live(A*)`, and that conversion factor is evaluated at
`A*(h)` -- the leaf area the plant's HEIGHT implies -- so at fixed height it is a
constant. Therefore `d(growth)/dA` and `dP/dA` have the same sign and **the
shedding criterion is unaffected**. For sapwood it is decisive, because extra
sapwood must be paid for in forgone leaf-area growth.

**2. First fix an inconsistency that biases the optimum.**
`dmass_sapwood_darea_leaf` uses `pars.theta`, not the plant's actual Huber value.
So a plant carrying 3.3x the sapwood per leaf currently pays the *pipe-model*
price for new growth and gets the hydraulic benefit free. That is much of why P
peaked as far out as 3.3x. Fixing it moves the growth optimum well below the P
optimum. Its own commit, its own number check -- it changes behaviour on its own.

**3. `psi` integrating the marginal return IS the slow, time-averaged response.**
`dpsi/dt = a_sw * R_s`, with `R_s` the dimensionless marginal growth return of
sapwood (zero at the optimum), is an integral controller: it converges exactly
onto `R_s = 0`, and a small `a_sw` makes it slow. **No extra tracked state is
needed, because integration is averaging** -- which is what the tracked marginal
for the shedding gate had to be built for, and is not needed twice.

**4. ⚠️ RETRACTED — see "What the implementation found" below. The derivative is NOT nearly free; the envelope theorem does not apply to the collar potential, which is root-found rather than maximised. The reasoning as originally written:**

**The derivative is nearly free, by the envelope theorem.** `R_s` needs
`d(profit)/d(kmax)`, which looked like a second leaf solve per rate evaluation --
a doubling of the hot path. It is not. `profit_` is already *maximised* over the
collar potential, so the derivative with respect to a parameter is the partial
derivative of the objective at the optimum already found:

    d(profit*)/d(kmax) = (partial profit / partial kmax) at psi*

the indirect term vanishing because `partial profit / partial psi = 0` there. It
holds at the KKT corner too: the active constraint is `psi <= psi_crit`, and
`psi_crit` is built from `stem_b`/`stem_c`, not from `kmax`, so the multiplier
term drops out. **One evaluation of a closed form at a known point, not a
re-optimisation.**

This is what phylloptim's IFT gradient machinery (#4) exists for; the known gap is
that its gradient outputs are `A`, `gc`, `psi_stem` and `collar` rather than the
`profit_` plant bills -- plant #614, and the same sensitivity the optimality epic
already needed.

### Order

1. Construction cost reads the actual Huber value, not `theta`. Behaviour change,
   own commit.
2. Re-measure the optimum against growth rate. Expect it well below 3.3x.
3. Add `d(profit)/d(kmax)` to phylloptim by the envelope theorem.
4. The integral controller, with `a_sw = 0` by default so exactness holds.

⚠️ And #617 (epic #615) replaces the height-resistance relation this is all
measured against, so step 2's number should be re-taken after it lands.

### What the implementation found (2026-09-08)

Steps 1, 2 and 4 are done. Step 3 is not needed as an upstream change: `d(profit)/d(kmax)` is measured inside plant. plant #614 remains worth having for the analytic version, and is now worth MORE than it looked, for the reason below.

**⚠️ RETRACTION: the envelope theorem does not apply, and point 4 above is wrong.** The shortcut was to re-evaluate the leaf at the collar potential already found, on the argument that the indirect term vanishes at an optimum. That holds for `opt_psi_stem_`, which is genuinely chosen to maximise profit. It does NOT hold for `opt_root_psi_`, which is found by a **root find** — the collar potential at which the soil-root network's supply matches the leaf's demand. A constraint is not an optimum, so `d(profit)/d(collar)` is not zero and the indirect term is real.

Measured: holding the collar fixed under-reports `d(profit)/d(kmax)` by **15-19 per cent**, which put the controller's zero at 1.22x the pipe-model ratio where growth actually peaks at 1.28x. A full re-solve matches a finite difference of the model's own assimilation to **0.09 per cent** and lands the zero on the peak in all twelve height/soil cells tried.

**Why it survived so long:** in wet soil the growth peak is extremely flat (at 16 m, `dh/dt` differs by 6e-4 relative between 1.13x and 1.22x), so a 0.025 scan grid could not resolve a 0.03 offset in psi and the zero appeared to land exactly on the peak. It only showed up in dry soil, on a 0.01 grid. **A coarse grid does not merely lose precision here; it manufactures an exact-looking agreement.**

**Consequence:** `a_sw > 0` costs a second full leaf solve per rate evaluation (three in total: the original, the perturbed, and one to restore the leaf's members for the auxes). The default `a_sw = 0` skips the block entirely, so the hot path is untouched unless acclimation is on. An analytic `d(profit)/d(kmax)` from phylloptim would remove that cost.

**The controller.** `dpsi/dt` gains `a_sw * R_s`, with

    R_s = availability * (dP/dpsi - P * dmass_sapwood_darea_leaf * darea_leaf_dmass_live) / maintenance

The second term in the bracket is the price of the extra stem, paid in forgone leaf area; dropping it optimises production instead, whose optimum sits far out near 3.3x. `a_sw = 0` is the default and the model is then **bit-exact** — verified by building at HEAD and re-running the scenario gateway, which returned the same eight numbers to every digit.

**Four things the implementation got wrong first, each worth stating because each looked right:**

1. **⚠️ Sapwood turnover is GATED, so its `psi`-derivative is not `k_s * m_s`.** What the budget is charged is `replacement_sapwood * k_s * m_s` with `replacement_sapwood ∝ exp(-psi/a_pl2)` against `m_s ∝ exp(psi)`, so the charge falls as `exp((1 - 1/a_pl2) psi)` — steeply, `a_pl2` being well below 1. Differentiating the ungated cost drops a term worth **+6.8 kg/yr of 13.1** and moved the controller's zero from the true 1.35x down to 1.05x. This is the single largest error in the derivation and it is invisible to any sign or monotonicity check.

2. **The denominator must be strictly positive, and `profit` is not.** `profit + shadow_cost` is net of leaf respiration and goes negative in exactly the drought where the controller most needs a sign. Normalising by it silently switched the controller off at soil 0.13. It now divides by the maintenance bill (`a_bio * a_y * respiration + turnover`), which is positive whenever the plant has tissue.

3. **Re-proportioning must be paid for.** Nothing debits the budget for the `a_sw` drift, so without a gate a plant in carbon deficit keeps thickening its stem on carbon it does not have. **The reserve gate `G` is NOT the right switch**: at soil 0.13 reserves are still full so `G ~ 1`, and what has gone to zero is `Ppos`. The factor used is `growth_flux / (growth_flux + maintenance)` — a share of throughput, so smooth, in `[0, 1)`, vanishing with the growth flux, and strictly positive so **it cannot move the zero**, only the speed.

4. **Deep-crown is refused, not skipped.** The sensitivity is measured at one radiation; deep-crown integrates profit over crown positions, so that difference is not its derivative. Skipping it leaves the controller with only its negative cost term, and the stem shrinks without bound while the run looks plausible the whole way down. `prepare_strategy()` now throws on `a_sw > 0` with deep-crown.

**On the FD step.** The first diagnosis was that `dk = 1e-6 * kmax` sat in the solver's noise floor. It did not — widening it 1000x moved `R_s` by 0.3%, which also rules out a constant offset between a point value and an integral. The step is now `1e-3` relative on the reasoning that it should clear the leaf's internal tolerances, but the recorded finding is that **the step was never the problem**; the missing turnover term was.

**Measured, at h = 10 m, soil 0.20, against the pipe-model ratio:**

| stem/leaf | dh/dt | R_s |
|---|---|---|
| 1.00 | 1.191 | +0.264 |
| 1.22 | 1.393 | +0.089 |
| **1.35** | **1.411** | **−0.002** |
| 1.49 | 1.396 | −0.085 |
| 1.82 | 1.303 | −0.221 |

The zero-crossing lands on the argmax of `dh/dt` — which is the defining correctness property, since a controller built on production instead is still monotone, still converges, and is still wrong. `test-strategy-tf24.R` asserts it on the grid rather than by interpolation so that the production version fails rather than passing on a tolerance.

At soil 0.15 growth is still rising at 1.82x and `R_s` stays positive throughout; at soil 0.13 the plant has no carbon and `R_s` is ~1e-11. In an SCM run on a wet stand, `a_sw = 0.5` raises R0 from 97.3 to 154.0 (+58%), which is the size of the prize in leaving the pipe model.

## Shedding's benefit is a steep function of SIZE — re-measured after #617 (2026-09-08)

**The prediction below was made before #617 landed and it held.** It said the benefit band sits at the height-resistance relation's viability ceiling, and that the path integral should push the ceiling up and take the band with it. Measured after the merge:

| | before #617 | after #617 |
|---|---|---|
| viability ceiling (P > 0 in wet soil) | between 17 and 18 m | between 30 and 36 m |
| shedding benefit first appears | ~16 m | ~20 m |
| largest benefit found | 6.1e9x (at 18 m, and both arms dead) | 84x (at 30 m, 4 yr @ soil 0.14) |
| benefit at soil 0.16 for 4 yr, any height | up to 4.4e9x | **exactly 1.00 at every height to 24 m** |

The ceiling roughly doubled, which is what conduit widening is for. Three things follow.

**1. TF24 is now far more drought-tolerant, and the earlier "cliff" is gone.** A 4-year drought at soil 0.16 no longer kills anything up to 24 m (survivorship ~0.92 at every height) and the canopy never thins, so shedding has nothing to do. Before #617 that same drought took a 20 m plant to 2.9e-19.

**2. The size dependence survives, shifted.** At soil 0.14 the benefit is 1.00 at 10 m, 1.04 at 16 m, 2.4 at 20 m, 21 at 24 m and 84 at 30 m. Below ~13 m the canopy still does not thin at all: the marginal leaf keeps paying for itself, which is the gate working as designed.

**3. Shedding still never converts death into survival.** Zero cells of the post-#617 grid (heights 16-30, soils 0.14-0.12, 4 yr) have the fixed plant dead and the shedding plant alive; the large ratios still sit between two very small numbers. And the benefit is largest at *moderate* drought (0.14), not the harshest — at soil 0.12 everything dies whatever it does, so the ratio returns to ~1.

⚠️ **The numbers in the original section below are PRE-#617 and superseded.** They are kept because the reasoning they support — that the band tracks the viability ceiling — is what the merge then confirmed, and because the size of the shift is the clearest available statement of how much the height-resistance relation was driving this.

## Shedding's benefit is a steep function of SIZE — the pre-#617 measurement (2026-09-08)

Re-asking "does shedding buy survival?" against a well-configured stem, as the earlier section said had to be done. The answer is more interesting than either yes or no, and **the earlier null result was measured at one size**.

Measured on a single plant: wet acclimation, mortality integral zeroed at the start of the drought so only the drought's own hazard is compared, then drought, then recovery.

| height | 4 | 8 | 10 | 13 | 16 | 18 | 20 |
|---|---|---|---|---|---|---|---|
| survival ratio, shed / fixed (4 yr @ 0.16) | 1.00 | 1.00 | 1.00 | 1.03 | 9.5 | — | 4.4e9 |
| canopy retained at the end | 1.00 | 1.00 | 1.00 | 1.00 | 0.74 | — | 0.38 |

**Below about 13 m shedding does literally nothing** — the ratio is 1.000 and the canopy never leaves 1.000. Not a weak effect: no effect. The gate reads the *marginal* leaf's balance, and a short plant's marginal leaf still pays for itself even in drought. This is why the earlier sweep, run at 10 m and below, found 1.001-1.06x and concluded thinning does not buy survival. That conclusion was right about the size it was taken at and wrong as a general statement.

**⚠️ But the large ratios are not what they look like.** At 20 m and above the plant has `P < 0 in WET soil` (-17.9 kg/yr at 20 m, -84.2 at 24 m). Those cells are not a tree surviving a drought; they are a tree that this height-resistance relation cannot sustain at all, which shedding then rescues. A ratio of 4e9 between 3e-19 and 1e-9 is two kinds of dead.

**The viability ceiling and the benefit band are the same thing.** `P` in wet soil crosses zero between 17 and 18 m:

| height | 14 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|
| P at soil 0.30 (kg/yr) | +9.6 | +6.6 | +3.2 | **-1.8** | -8.7 | -17.9 |
| shed/fixed survival, 1 yr @ 0.15 | 2.6 | 4.8 | 10.6 | 6.1e9 | 4.9e9 | 2.1e9 |

So the whole of shedding's leverage lives in a narrow band straddling the height at which the model stops being able to keep a tree alive. Of a 90-cell grid (heights 8-18, soil 0.20-0.14, droughts 1-4 yr), **zero cells** had shedding take a plant from dead to comfortably alive by a strict bar; the best real rescue is 18 m, 1 yr at soil 0.15, where survivorship goes from 2.2e-11 to **0.137**.

### What this means for the work

1. **The demo's "thinning does not buy survival" section needs qualifying, not deleting.** It is true below 13 m and false above 16 m, and the reason is the one the analytical criterion already gives: eta grows with height because `k_max ∝ 1/h`, so tall plants are exactly the ones for which thinning pays. The criterion predicted this; the behavioural test was taken at a size where the criterion says the answer is no.

2. **⚠️ Every number here is provisional against #617**, and more so than anything else measured on this branch. The band's location IS the height-resistance relation. #617 replaces resistance-linear-in-height with a path integral over conduit widening, and basipetal widening is precisely what stops resistance growing linearly with height in real trees — so it should move the ceiling up and take the band with it. Re-measure this table after #617, before drawing any ecological conclusion from it.

3. **Sapwood acclimation helps more than shedding does, in the range where plants are viable.** At 10 m, `a_sw = 0.5` raised drought survivorship about 10x (7.7e-3 to 7.9e-2 at soil 0.16) where shedding gave 1.002x. Both mechanisms matter, but they matter at different sizes.

## Open

- Which growth rate the optimality criterion should maximise is **not** a modelling choice to be argued — see "Deciding the objective by invasion analysis" above. It is settled by experiment, and the cheap half of that experiment can run before any of this code exists.

## Settled

- **The new gate centres and the rebuild rate are `TF24_Pars` fields, and therefore settable as traits.** Nothing further is needed for that: plant's trait mechanism reaches nested `pars`, so a parameter that lives there can be varied through `trait_matrix()`/`add_strategies()` immediately — verified on `theta`, which takes two distinct values across two strategies with no hyperpar involvement. A **hyperpar** entry is a different thing and is *not* wanted yet: that machinery is for a parameter *derived* from another trait (as `k_l` is from `lma`), and there is no trade-off to encode until #512 says what one would look like. Defaults reproduce the current model.
