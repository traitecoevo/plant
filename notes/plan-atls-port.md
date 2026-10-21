# Porting the ATLS thermal-damage layer to develop (#566)

**Status: in progress (started 2026-10-21); both halves implemented, item 4 blocked.** The ATLS layer was built in July on `feature/tf24-atls-thermal-damage` (tip `bee0a02a`, TF24t@v3.3) against plant's own leaf model. #591 then moved the leaf into phylloptim and deleted 23 of the files that branch changes, so it cannot merge. This note is the plan for re-implementing it on current code.

## Decisions

1. **Faithful port, not a redesign.** Port the v3.3 equations, states, parameters and defaults as committed. Design changes (damage entering the optimiser's profit, faster acclimation, separating T_crit from T_opt, condition-dependent repair) come afterwards as separate PRs, each measured against this port as its baseline. Decided by Daniel, 2026-10-21.
2. **Two halves, two repos.** Leaf-level machinery goes into phylloptim's `Leaf` behind `use_thermal_damage_` (default off); the TF24t strategy, its ODE states and whole-plant costs stay in plant.
3. **Default off is bit-identical.** With `use_thermal_damage_` false, every phylloptim golden file and every existing plant result is unchanged bit-for-bit. Assert with `expect_identical`, and pair each inertness test with a non-zero check.
4. **The plant half is stacked on #678** (odelia v0.8.0 + phylloptim v0.10.3), because phylloptim master requires odelia >= 0.8.0. Rebase onto develop after #678 squash-merges.
5. **No driver rename.** The July branch renamed the `leaf_temp` driver to `air_temp`. develop still calls it `leaf_temp`, and with the energy balance on phylloptim already reads it as air temperature. TF24t reads `get_leaf_temp()` as air temperature; the rename is a separate TF24-wide change.
6. **Respiration follows current phylloptim.** On the July branch `R_d = 0.015 × Vcmax` after the damage discount, so damaged leaves respired less. phylloptim#86 made `R_d_25` a trait with a Tjoelker Q10, independent of Vcmax, for every model. The port keeps that and adds the ATLS maintenance costs on top, so damage no longer lowers R_d. This is the one intended divergence from v3.3.
7. **The `operating_leaf_temp` aux is dropped.** develop's TF24 already reports `Tleaf` (#645); TF24t's running-mean state relaxes toward the leaf's `Tleaf_`.
8. **TF24t refuses the DeepCrown shading model.** The leaf's `phi_d_`, `k_rec_eff_` and `Tleaf_` are read after the crown loop, which under DeepCrown would give the last node's values. MeanLight (default) and CrownCentre evaluate one leaf and are supported.
9. **No gradients through the thermal layer.** phylloptim's gradient, derivative-row and closed-form paths refuse `use_thermal_damage_` with an error rather than return a derivative that ignores it. In plant, the thermal parameters are declared `no_gradient`.
10. **Versions.** phylloptim takes the next minor after #120 (#120 is 0.11.0, so this is 0.12.0; renumber if they merge in the other order). TF24t's compound version stays `<TF24 version>.<thermal_revision>` with `thermal_revision = 3`.

## Interface between the halves

phylloptim `Leaf` gains (all R-bound):

- inputs set by the strategy: `use_thermal_damage_`, `topt_offset_`, `m_rep_`, `t_rep_cut_`, `k_i_`, `k_i_ref_`, `k_rec_`, `dTopt_max_`, `K_A_`, `A_`, `I_r_`, `I_p_`, `c_acclim_maint_`, `c_repair_maint_`, `c_protect_maint_`, `c_repair_flux_`;
- outputs read back after the solve: `phi_d_` (deactivated fraction at the operating leaf temperature) and `k_rec_eff_` (repair rate gated off above `t_rep_cut_`);
- static helpers `logistic_` and `softplus_`, which TF24t's acclimation forcing also uses.

`k_mat` is read only by the strategy's ODE, so it lives in plant alone.

## Invariants a later edit could break

- The temperature cache (`photo_temp_key`) must not serve a cached value when `use_thermal_damage_` is on, because `A_`, `I_r_` and `I_p_` change every step.
- The survival fraction `1 − I_r − I_p` is floored at `1e-6` so the leaf solve stays well-conditioned as damage approaches 1.
- Every rate constant is per day; TF24t multiplies by `DAYS_PER_YEAR = 365.25` for the yearly ODE clock, except leaf turnover and growth dilution, which are already per year.
- Damage enters the benefit (`Vcmax` and `Jmax` both scaled by survival), not the hydraulic cost.

## Order of work

- [x] phylloptim: leaf-level layer, tests, NEWS, PR — [phylloptim#162](https://github.com/traitecoevo/phylloptim/pull/162), open; plant pins its head `c57c6d33`.
- [x] plant: TF24t on the templated TF24 strategy, registration, R helpers, tests, model-version snapshot; pin phylloptim to the PR head while both are open — on `feature/tf24t-atls-port`, stacked on #678. Re-pin to phylloptim's squash sha after #162 merges, and rebase onto develop after #678 merges.
- [x] plant: carry over the July notes (`notes/atls/`) and the July viability results, marked as computed on the July code — in `964638d2`. The July staging docs (`overstorey_staging/TF24t_leaf_ATLS_model.qmd`, `TF24t_leaf_ATLS_demo.qmd`, `atls_demo_helpers.R`) are ported to the current API on the same branch; the demo's two SCM sections are not re-rendered (see below).
- [ ] Re-run the materiality gate and the low-density viability boundary on the ported code, on the birth-date coordinate. **Blocked on the TF24 baseline below**: as set up in July (lma 0.0825, birth rate 1) the stand is inviable whenever the energy balance is on, so the comparison crosses no boundary until it is re-set.

## Invariants the plant half adds

- **`TF24t_Strategy::rebind_from` is deleted, and must stay deleted** until TF24t has a derivative of its own. The inherited one returns a `TF24_Strategy<U>`, so it would rebind a TF24t patch as a TF24 patch — no thermal states, no layer, every number finite. Deleted, Patch's `at_scalar` and odelia's `Rebindable` read false and every AD route (stand gradient, sweep, an implicit stepper's Jacobian) refuses at compile time. `census_backend()` names TF24t with the reason.
- **The thermal parameters are top-level strategy fields, not `TF24_Pars`**, as TF24f's are. A trait matrix cannot reach them (`add_strategies` refuses them by name); set them on the strategy. The drift guard snapshots them through `top_level_science` in `test-model-version.R`, so a new TF24t parameter goes there too, and a changed default needs a `thermal_revision` bump.
- **At non-damaging settings a rate evaluation is bit-identical to TF24 with `pars$use_energy_balance = 1`** (alpha = 0, k_i = 0, topt_offset = 0, all costs 0). The one operation that is not an identity by construction is the ΔS round trip at zero shift, which is exact at phylloptim's default temperature parameters. A free SCM run is not identical (the four extra states sit in the step-size controller's error norm: 260 against 266 steps); TF24 pinned to TF24t's own step program is. `test-strategy-tf24t.R` asserts both.
- **`Individual$set_initial_states` must run before an `IndividualRunner`** is built from a TF24t plant. The runner does not seed, and an unseeded storage pool starts empty, where a plant losing carbon has every step that dips below it refused.

## Known issues carried over, not fixed here

- φ_d drives damage accumulation but never enters the step's profit, so the optimiser cannot avoid damage.
- At the default `k_mat = 0.02 /day` permanent damage ratchets to inviability in every climate; the July runs used 0.00125.
- Acclimation (about 100 days to equilibrate, 6 °C cap) is slower and smaller than Piet Arnold's evidence (hours to days, about 12 °C).
- Repair depends on temperature only; no within-day dose or overnight recovery.

## Found while porting the plant half

- **TF24 with the energy balance on is inviable in the July viability setup.** Single species, lma 0.0825, birth rate 1, full patch lifetime: low-density R0 is 1e-8 at 15 °C, 5e-15 at 20 °C, 7e-16 at 24 °C and 1e-18 at 28 °C, against 56 at 25 °C with the balance off (July: 57 at 24 °C, crossing 1 near 28 °C). TF24t inherits it, since it forces the balance on. Cause not investigated; candidates are the hotter operating leaf and phylloptim's Q10 respiration rising with it. This blocks the last item above.
- **The operating leaf runs hotter than in July.** A lone plant in full light at 25 °C air sits 6–9 °C above it (July: about 4 °C above 27 °C air). At default TF24t settings a lone plant at 25 °C air is dead within 0.6 years (I_p 0.74, storage empty), through damage, closing stomata and a hotter leaf. At 15 °C air it grows.
- **The stochastic path reads stale thermal state.** `establishment_probability(env)` (no internals) solves the species' shared leaf with whatever A, I_r and I_p the last rate evaluation of any individual left there. The SCM uses the newborn's own rates; the stochastic runner does not. Same as July.
