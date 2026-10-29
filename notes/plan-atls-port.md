# Porting the ATLS thermal-damage layer to develop (#566)

**Status: in progress (started 2026-10-21); both halves implemented and redesigned to three states (TF24t thermal revision 4, phylloptim#162 at 60aef98); item 4 blocked.** The ATLS layer was built in July on `feature/tf24-atls-thermal-damage` (tip `bee0a02a`, TF24t@v3.3) against plant's own leaf model. #591 then moved the leaf into phylloptim and deleted 23 of the files that branch changes, so it cannot merge. This note is the plan for re-implementing it on current code.

## Decisions

1. **Faithful port, not a redesign.** Port the v3.3 equations, states, parameters and defaults as committed. Design changes (damage entering the optimiser's profit, faster acclimation, separating T_crit from T_opt, condition-dependent repair) come afterwards as separate PRs, each measured against this port as its baseline. Decided by Daniel, 2026-10-21.
2. **Two halves, two repos.** Leaf-level machinery goes into phylloptim's `Leaf` behind `use_thermal_damage_` (default off); the TF24t strategy, its ODE states and whole-plant costs stay in plant.
3. **Default off is bit-identical.** With `use_thermal_damage_` false, every phylloptim golden file and every existing plant result is unchanged bit-for-bit. Assert with `expect_identical`, and pair each inertness test with a non-zero check.
4. **The plant half is stacked on #678** (odelia v0.8.0 + phylloptim v0.10.3), because phylloptim master requires odelia >= 0.8.0. Rebase onto develop after #678 squash-merges.
5. **No driver rename.** The July branch renamed the `leaf_temp` driver to `air_temp`. develop still calls it `leaf_temp`, and with the energy balance on phylloptim already reads it as air temperature. TF24t reads `get_leaf_temp()` as air temperature; the rename is a separate TF24-wide change.
6. **Respiration follows current phylloptim.** On the July branch `R_d = 0.015 × Vcmax` after the damage discount, so damaged leaves respired less. phylloptim#86 made `R_d_25` a trait with a Tjoelker Q10, independent of Vcmax, for every model. The port keeps that and adds the ATLS maintenance costs on top, so damage no longer lowers R_d. This is the one intended divergence from v3.3.
7. **The `operating_leaf_temp` aux is dropped.** develop's TF24 already reports `Tleaf` (#645); TF24t's running-mean state relaxes toward the leaf's `Tleaf_`.
8. **TF24t refuses the DeepCrown shading model.** The leaf's damage rates are evaluated at its `Tleaf_` after the crown loop, which under DeepCrown would give the last node's values. MeanLight (default) and CrownCentre evaluate one leaf and are supported.
9. **No gradients through TF24t.** In plant the thermal parameters are declared `no_gradient` and TF24t has no stand gradient. (Superseded in phylloptim: its trait gradient now carries the layer; see "Design changes after the port".)
10. **Versions.** phylloptim takes the next minor after #120 (#120 is 0.11.0, so this is 0.12.0; renumber if they merge in the other order). TF24t's compound version stays `<TF24 version>.<thermal_revision>`; the port was `thermal_revision = 3`, the three-state redesign is 4.

## Interface between the halves

phylloptim `Leaf` gains (all R-bound):

- parameters set by the strategy: `use_thermal_damage_`, `topt_offset_`, `impair_H_`, `impair_Tmid_`, `k_i_`, `k_i_Ea_`, `k_i_Tref_`, `k_rec_`, `m_rep_`, `t_rep_cut_`, `k_mat_`, `acclim_alpha_`, `acclim_beta_`, `acclim_T_`, `acclim_s_`, `dTopt_max_`, `dTimpair_max_`, `K_A_`, and the four `c_*` maintenance costs (`k_i_baseline_` keeps the leaf default);
- states handed over before each solve: `A_`, `I_`, `D_`;
- outputs read back after the solve: `k_imp_eff_`, `k_rec_eff_` and `functional_`, at the operating leaf temperature;
- `damage_rates(T, A, I, D)`, the per-day right-hand side, which TF24t multiplies by 365.25 and adds dilution of D to; and the static `softplus_`, used for the birth equilibrium of A.

## Invariants a later edit could break

- The temperature cache (`photo_temp_key`) must not serve a cached value when `use_thermal_damage_` is on, because `A_`, `I_` and `D_` change every step.
- The functional fraction `1 − I − D` is floored at `1e-6` where it scales capacity, so the leaf solve stays well-conditioned as I + D approaches 1.
- Every rate constant is per day; TF24t multiplies by `DAYS_PER_YEAR = 365.25` for the yearly ODE clock, except leaf turnover and growth dilution, which are already per year.
- Impairment enters the benefit (`Vcmax` and `Jmax` both scaled by F), not the hydraulic cost.

## Order of work

- [x] phylloptim: leaf-level layer, tests, NEWS, PR — [phylloptim#162](https://github.com/traitecoevo/phylloptim/pull/162), open; plant pins its head `60aef98`.
- [x] plant: TF24t on the templated TF24 strategy, registration, R helpers, tests, model-version snapshot; pin phylloptim to the PR head while both are open — #694 (`feature/tf24t-atls-port`), a draft stacked on #678. Re-pin to phylloptim's squash sha after #162 merges, and rebase onto develop after #678 merges.
- [x] plant: carry over the July notes (`notes/atls/`) and the July viability results, marked as computed on the July code — in `964638d2`. The July staging docs (`overstorey_staging/TF24t_leaf_ATLS_model.qmd`, `TF24t_leaf_ATLS_demo.qmd`, `atls_demo_helpers.R`) are ported to the current API on the same branch; the demo's two SCM sections are not re-rendered (see below).
- [ ] Re-run the materiality gate and the low-density viability boundary on the ported code, on the birth-date coordinate. **Blocked on the TF24 baseline below**: as set up in July (lma 0.0825, birth rate 1) the stand is inviable whenever the energy balance is on, so the comparison crosses no boundary until it is re-set.

## Invariants the plant half adds

- **`TF24t_Strategy::rebind_from` is deleted, and must stay deleted** until TF24t has a derivative of its own. The inherited one returns a `TF24_Strategy<U>`, so it would rebind a TF24t patch as a TF24 patch — no thermal states, no layer, every number finite. Deleted, Patch's `at_scalar` and odelia's `Rebindable` read false and every AD route (stand gradient, sweep, an implicit stepper's Jacobian) refuses at compile time. `census_backend()` names TF24t with the reason.
- **The thermal parameters are top-level strategy fields, not `TF24_Pars`**, as TF24f's are. A trait matrix cannot reach them (`add_strategies` refuses them by name); set them on the strategy. The drift guard snapshots them through `top_level_science` in `test-model-version.R`, so a new TF24t parameter goes there too, and a changed default needs a `thermal_revision` bump.
- **At non-damaging settings a rate evaluation is bit-identical to TF24 with `pars$use_energy_balance = 1`** (alpha = 0, k_i = 0, topt_offset = 0, all costs 0). The one operation that is not an identity by construction is the ΔS round trip at zero shift, which is exact at phylloptim's default temperature parameters. A free SCM run is not identical (the four extra states sit in the step-size controller's error norm: 260 against 266 steps); TF24 pinned to TF24t's own step program is. `test-strategy-tf24t.R` asserts both.
- **`Individual$set_initial_states` must run before an `IndividualRunner`** is built from a TF24t plant. The runner does not seed, and an unseeded storage pool starts empty, where a plant losing carbon has every step that dips below it refused.

## Known issues carried over, not fixed here

- k_imp drives impairment but never enters the step's profit, so the optimiser cannot avoid impairment.
- At the default `k_mat = 0.02 /day` the damaged fraction ratchets to inviability in every climate; the July runs used 0.00125.
- At the default switch (midpoint 34.5 °C, impair_H 200 kJ/mol) about 12% of leaf is impaired at equilibrium at a 25 °C leaf, which costs 12% of capacity. The port's recoverable pool did the same; the switch needs recalibrating.
- Acclimation (about 100 days to equilibrate, 6 °C cap) is slower and smaller than Piet Arnold's evidence (hours to days, about 12 °C).
- Repair depends on temperature only; no within-day dose or overnight recovery.

## Design changes after the port

Each lands as its own PR, measured against this port as its baseline (decision 1).

**Decided (Daniel, 2026-10-29): the damage onset gets its own parameters, independent of Jmax's.** Implemented on phylloptim#162 (folded in from #165). On the port, φ_d reuses the Jmax curve's deactivation terms (`jmax_H_d_`, `jmax_d_S_`, in `deactivated_fraction` and `d_S_damage_shifted`), so recalibrating the photosynthetic temperature response silently moves the damage onset. The impairment equilibrium gets its own enthalpy `impair_H_` and midpoint `impair_Tmid_`, written as a centred Arrhenius term, K = exp[(impair_H_/R)(1/Tmid − 1/T)], and K/(1 + K) a switch on the impairment rate.

**Decided (2026-10-29): process-agnostic names.** The three states are functional ⇌ impaired → damaged. phylloptim#162 reports the impairment rate `k_imp_eff_` (the port reported `phi_d_`); "deactivation" is kept for the Vcmax/Jmax curves' own instantaneous term.

The interpretation this fixes: the peaked-Arrhenius Vcmax and Jmax curves are the **instantaneous, fully reversible** temperature response, and nothing about them persists past the time step. Impairment is a **separate process in addition to it**, and is the only route into lasting loss (I, D), which then scales capacity through F. So a hot leaf loses capacity twice over, by different mechanisms on different timescales, and that is intended rather than a double count. Required with it:

- Default `impair_H_` = 200 kJ/mol and `impair_Tmid_` = 200000/650 K (34.54 °C), so the defaults reproduce the port's φ_d to rounding (measured ≤ 1.3e-15).
- φ_d stops being "the deactivated fraction of the curve the capacities follow"; the comments and test in phylloptim that claimed it are rewritten in #162.
- Acclimation shifts the switch midpoint by its **own** saturating amount, `dTimpair_max * A/(A + K_A)`, with `dTimpair_max` defaulting to `dTopt_max`'s 6 °C so the default reproduces the port. Decided 2026-10-29.

**Every new parameter defaults to the value that reproduces the port**, so each change below is checked against the port at defaults before any default moves. Decided 2026-10-29.

**Decided (2026-10-29): the irreversible step is temperature-dependent.** `k_i` becomes a centred Arrhenius rate, `k_i_` · exp[(`k_i_Ea_`/R)(1/T_ref − 1/T)] at `k_i_Tref_`, multiplying the switch inside `k_imp_eff_` (functional ⇌ impaired → damaged is unchanged). Default `k_i_Ea_` = 0 gives exactly `k_i_`. Implemented on phylloptim#162; the protection baseline `k_i_ref_` is renamed `k_i_baseline_` there so it cannot be read as the Arrhenius reference. The switch's ceiling of 1 keeps the impairment rate bounded only while E_i is moderate: an uncapped Arrhenius damage rate reaches ~2.4e3 /yr at 45 °C and ~2.4e4 /yr at 55 °C on the per-year clock, against ~180 /yr capped, so watch stiffness as E_i is fitted.

**Decided (2026-10-29): `k_mat` stays constant in temperature**, unless recovery time courses at several temperatures say otherwise.

**Decided (2026-10-29): the default `k_mat` changes.** 0.02 /day ratchets permanent damage to inviability in every climate, which makes it unusable as a default. The July runs' 0.00125 /day is the candidate, pending a fit. This is a default change, so it is its own PR and moves TF24t results.

**Decided (2026-10-29): acclimation kinetics are open to fitting, and the defaults are too slow.** What sets them today:

- Speed is β alone: A relaxes on a 1/β = 100-day timescale. Faster acclimation (hours to days in the project's data) means raising β, and α with it if the equilibrium A = (α/β)·softplus(T_air − t_accl) is to stay put.
- Size is `dTopt_max`, the shift as A → ∞, through δ = topt_offset + dTopt_max·A/(A + K_A). It is a parameter, not a hard cap, and 6 °C is a placeholder against the data's ~12 °C. The only hard bound is the [250, 350] K clamp on the shifted optimum, a numerical guard.
- α, β, t_accl, softplus_s, dTopt_max and K_A all need to be in whatever fitting interface is built; plant currently declares them `no_gradient`, so fitting them means finite differences or a TF24t gradient.
- Open: whether A and K_A are separately identifiable from shift data alone. Only the product of the equilibrium A and 1/K_A sets the shift, so one may need fixing.

**Decided (2026-10-29): the leaf-side changes above are part of phylloptim#162 itself**, not a separate PR: drafted as #165 and folded into #162's branch (one 0.12.0 NEWS entry). The `k_mat` default and the acclimation defaults are plant-side and go separately.

**Done in plant (#694, re-pinned to 60aef98):** TF24t exposes and copies `impair_H`, `impair_Tmid`, `dTimpair_max`, `k_i_Ea` and `k_i_Tref`, which are in `top_level_science` in `test-model-version.R`, and `thermal_revision` is 4.

**Decided (Daniel, 2026-10-29): the leaf's damage dynamics move into phylloptim, so leaf experiments can be fitted with phylloptim alone.** This revises decision 2 ("two halves, two repos"). Implemented in phylloptim#162 (`Leaf::damage_rates`, `leaf_damage_trajectory()`) and called by TF24t in #694.

- **Moves to phylloptim:** the rate equations for acclimation, impairment and recoverable damage, maturation into permanent damage, and their parameters (`k_mat`, `alpha`, `beta`, `t_accl`, `softplus_s`); a pure rates function of leaf temperature and the states, returning per-day rates; and a trajectory integrator over a temperature time series, bound to R.
- **Stays in plant:** clearing permanent damage by leaf turnover and canopy growth, the per-day to per-year clock, the running-mean leaf temperature, the whole-plant construction and induction costs, and crown aggregation. TF24t calls the leaf's rates function and adds these.
- **Comparison to data lives outside phylloptim**, in an analysis repo: the mapping from model states to measurements (e.g. Fv/Fm) is not phylloptim's.

**Decided (Daniel, 2026-10-29): phylloptim's gradients include the damage states and parameters where present.** Implemented in #162. With the gate on, the IFT trait gradient takes A, I, D and the layer's parameters as inputs (34 in all); with it off, they are inert. `leaf_damage_gradient()` gives exact forward-mode sensitivities of the trajectory to the dynamics parameters and initial states, for fitting. plant's TF24t still has no gradient (decision 9).

**Decided (Daniel, 2026-10-29): three states, as fractions of leaf area.** functional F ⇌ impaired I → damaged D, with F = 1 − I − D. Implemented in #162 and #694. A four-state version (a fast impaired state plus recoverable damage) was considered the same day and dropped: the fast state and recoverable damage are merged into one impaired state, so the only instantaneous, reversible loss is the Vcmax/Jmax peaked Arrhenius, and impairment is rate-and-state.

```
dA/dt = alpha * softplus(T_leaf - t_accl; s) - beta * A
dI/dt = k_imp(T) * F - (k_rec_eff(T) + k_mat) * I        (I -> F by repair)
dD/dt = k_mat * I                                       (plant subtracts dilution)
k_imp(T) = k_i * exp[(E_i/R)(1/T_ref - 1/T)] * K/(1 + K),  K = exp[(impair_H/R)(1/Tmid - 1/T)]
Vcmax, Jmax scaled by F (floored at 1e-6); R_d not scaled.
```

- **K/(1+K) is a temperature switch on the impairment rate, not a state or a fraction of the leaf.** Its midpoint moves with acclimation by its own `dTimpair_max`.
- **This is the port's structure renamed:** I is the old I_r and D the old I_p. At the same forcing temperature, the default rates reproduce the port.
- **Acclimation is forced by LEAF temperature**, not the air-temperature driver the port read. Its birth equilibrium still uses the air, because the leaf is not solved until the first rate evaluation.
- **Impaired and damaged area still respires** (R_d unscaled, as decision 6).
- **No quasi-steady-state option is needed**: the impaired state relaxes at the repair rate, ~0.25 /day, not in minutes.
- **The induction cost is charged per unit leaf area** (`c_accl_induct` in kg m⁻² per unit positive dA/dt), not per plant as in July. Found while implementing leaf forcing, 2026-10-29: A is seeded at the birth AIR temperature's equilibrium but now forced by the warmer leaf, so a newborn at 28 °C air has dA/dt = 16 /yr. Per plant, 0.01 × 16 = 0.16 kg/yr against a net production of 9e-5 starved every seedling and slowed a one-year stand from 0.9 s to over ten minutes. Per m² it is about a fifth of a seedling's net production. Awaiting Daniel's confirmation; the alternative is to keep per-plant units and seed A at the leaf temperature, which removes only the birth transient.

**Under consideration, not decided:**

- Damage in the stomatal decision. Awaiting a separate analysis; potentially yes.
- Within-day dose and overnight recovery. Needs time-varying temperature drivers first; future work.
- Replace the repair gate `logistic(−m_rep (T − t_rep_cut))` with a peaked Arrhenius (rising limb plus deactivation). The current gate equals a Boltzmann deactivation term with H = m_rep R T_cut² = 337 kJ/mol to within 9e-4 day⁻¹ of 0.25. Decide on whether recovery data show repair speeding up with warmth below the cut-off.

## Found while porting the plant half

- **TF24 with the energy balance on is inviable in the July viability setup.** Single species, lma 0.0825, birth rate 1, full patch lifetime: low-density R0 is 1e-8 at 15 °C, 5e-15 at 20 °C, 7e-16 at 24 °C and 1e-18 at 28 °C, against 56 at 25 °C with the balance off (July: 57 at 24 °C, crossing 1 near 28 °C). TF24t inherits it, since it forces the balance on. The PM re-evaluation (#693) measured the cause on develop: it is the larger leaf-to-air vapour deficit (1.7–2.4 kPa at the leaf against 1 kPa in the air), not leaf temperature as such. Fick run at PM's leaf temperature grows normally, and adding PM's leaf deficit reproduces PM's one-year growth to within 3% at 10–20 °C. Boundary-layer resistance and net radiation, both uncalibrated, swing stand offspring by orders of magnitude. This blocks the last item above until the energy balance is calibrated.
- **The operating leaf runs hotter than in July.** An undamaged lone plant in full light at 25 °C air sits 4–6 °C above it (July: about 4 °C above 27 °C air), and 9 °C once damage has closed its stomata. At default TF24t settings a lone plant at 25 °C air is dying within 0.6 years (I_p 0.74, storage empty, cumulative mortality 2.2; measured on the port). At 15 °C air it grows.
- **The stochastic path reads stale thermal state.** `establishment_probability(env)` (no internals) solves the species' shared leaf with whatever A, I and D the last rate evaluation of any individual left there. The SCM uses the newborn's own rates; the stochastic runner does not. Same as July.
