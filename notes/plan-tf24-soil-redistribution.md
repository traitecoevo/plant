# Stage A design pass: soil influx & redistribution for rainfall pulses (#522, #599)

**Status: complete. Conclusion — do not change the soil physics.** The redistribution change
this pass was convened to design is not needed and would cost more than it buys. The measured
faults lie in the integrator's error handling and in where exceptions are thrown, not in the
water balance. The one physics-adjacent change the pulse feature genuinely needs is a capacity
cap on the pulse jump itself.

Backing experiments: `notes/scripts/soil_cascade/` (standalone R replica of
`tf24_environment.h`'s `compute_rates` and `soil_K_from_soil_theta`, plus odelia's exact
`OdeControl` with plant's defaults). No plant build required.

Defaults used throughout: θ_sat = 0.428, θ_res = 1e-2, K_sat = 163.0411 m yr⁻¹,
n_psi = 6.57 (so the conductivity exponent is 2·n_psi+3 = 16.14), a_infil = 1, b_infil = 8,
depth = 1.5 m over 5 layers (dz = 0.3 m), and odelia control
`tol_rel = tol_abs = 1e-4, a_y = 1, a_dydt = 0, h_min = h_init = 1e-6, h_max = 5`.

## 1. The premise was wrong: θ_sat is already a barrier

The concern was that `water_input[i] = water_flux[i-1] = K(θ_{i-1})` ignores layer `i`'s
capacity to accept it, so a saturated layer over a dry one drives θ past saturation. The
inflow is indeed capacity-blind — but the *outflow* closes the loop, because both fluxes
saturate at the same `K_sat`:

| θ₁ (θ₁/θ_sat) | inflow K(θ₀=θ_sat) | outflow K(θ₁) | rate₁ (yr⁻¹) |
|---|---|---|---|
| 0.0500 (0.117) | 163.041 | 0.000 | **+543.47** |
| 0.3500 (0.818) | 163.041 | 6.339 | +522.34 |
| 0.3852 (0.900) | 163.041 | 29.769 | +444.24 |
| 0.4237 (0.990) | 163.041 | 138.627 | +81.38 |
| 0.4280 (1.000) | 163.041 | 163.041 | **0.00** |

At θ = θ_sat inflow and outflow are both `K_sat` and the rate is exactly zero. Since
`soil_K_from_soil_theta` clamps its argument to `[0, θ_sat]`
([tf24_environment.h:403](../inst/include/plant/models/tf24_environment.h#L403)), the rate cannot
be positive above saturation either. **In exact arithmetic θ cannot exceed θ_sat.** Every
observed excursion is explicit-integrator overshoot:

| step h (yr) | θ₁ after one Euler step from θ₁ = 0.05 |
|---|---|
| 1e-5 | 0.055 |
| 1e-4 | 0.104 |
| 1e-3 | 0.593 — past saturation |
| 1e-2 | 5.485 |
| 0.08 | **43.53** |

That last row reproduces the `θ = [-51.4, 52.0, …]` recorded on
`origin/fix/soil-saturation-guard` (commit `32240725`) to within the arithmetic — confirming
the mechanism is step size, not the water balance. The boundary layer is narrow because of the
exponent 16.14: the rate is still +444 yr⁻¹ at 90% of saturation and only collapses inside the
final ~1%.

## 2. Error control already rejects those steps

Replicating odelia's controller exactly, on the same worst case:

| h | θ₁ after | max &#124;yerr&#124; | verdict | h_next |
|---|---|---|---|---|
| 1e-4 | 0.079 | 5.4e-5 | accept | 1e-4 |
| 1e-3 | 0.159 | 8.7e-3 | **reject** | 3.8e-4 |
| 1e-2 | −0.562 | 0.395 | **reject** | 2.0e-3 |
| 0.08 | −4.844 | 3.16 | **reject** | 0.016 |

Integrating the worst case for 0.5 yr from a **fresh** solver completes cleanly in 15 steps
(21 at `h_max = 0.05`) with `max θ = 0.428` exactly. `h_max` barely matters. So neither the
physics nor the controller is at fault when the controller is allowed to do its job.

## 3. What actually kills the run: an inherited step, then one of two faults

`step_size_last` persists across `set_state_from_system()` — documented at length at
[stochastic_patch_runner.h:146-159](../inst/include/plant/stochastic_patch_runner.h#L146-L159).
So a leg that begins after a discrete change (an introduction, and in future a pulse) starts at
whatever step size the *previous, quiet* leg had grown to, which may be `h_max = 5`. Varying
only that starting step, with a leaf model that either throws or returns NaN on an
out-of-domain ψ probe:

| h at leg start | leaf throws | leaf returns NaN |
|---|---|---|
| 1e-6 (fresh) | completes, 16 steps | completes, 16 steps |
| 1e-3 | **dies at step 0** | **dies, NaN state accepted** |
| 0.08 | **dies at step 0** | **dies, NaN state accepted** |
| 1.0 / 5.0 | **dies at step 0** | **dies, NaN state accepted** |

Two independent faults, and the trigger is the inherited step:

**Fault 1 — odelia's error reduction does not propagate non-finiteness, so a step carrying a
NaN is accepted.** In `ode_control.hpp:75-84` the reduction is `rmax = std::max(r, rmax)`.
`std::max(a, b)` is `(a < b) ? b : a`, and NaN comparisons are unordered, so it returns NaN for
`a = NaN` — but on the *next* element a finite `a` returns that finite value and **wipes the
NaN**. The NaN component is therefore never accounted for, and there are two ways this ends in
an accepted step:

- **(a) the NaN survives to the end of the loop** (last element, or all of them): `rmax` stays
  NaN, both `rmax > 1.1` and `rmax < 0.5` are false, and control falls through to the branch
  that reports no shrink. The caller
  ([ode_solver_internal.hpp:302-327](../../odelia/inst/include/odelia/ode_solver_internal.hpp#L302-L327))
  branches solely on `step_size_shrank()`, so the step is committed.
- **(b) the NaN is wiped by finite elements whose own ratios are small**: `rmax` is finite and
  passes, so the step is accepted *carrying* a NaN.

Measured, both accept: `yerr = [1e-3, 1e-3, 1e-3, 1e-3, NaN]` and
`yerr = [NaN, 1e-12, 1e-12, 1e-12, 1e-12]`. **(b) is the mode the coupled runs hit.** `Inf`
rejects correctly in either position, via `> 1.1`.

⚠️ Note this makes the fault *positional*, which is why it survived so long — and why it bites
this caller in particular. `Patch` chains its ODE state **species first, environment last**, so a
stiff soil block sits in the trailing indices, i.e. exactly mode (a)'s window; and a single
poisoned soil layer among finite neighbours is mode (b). Fixed in **odelia PR #54** by breaking
on the first non-finite ratio, which removes the positional dependence. With that, the same runs
complete in 14–15 steps with 9–12 rejections and `max θ = 0.428`.

⚠️ **R's `max()` propagates NaN; C++'s `std::max` does not.** The replica scripts model
`std::max` explicitly (`cxx_max` in `exp2_ctl.R`). An earlier version of this note used R's
semantics and consequently overstated the fault as "any NaN is accepted" — it is narrower and
positional, per the above.

**Fault 2 — a throw during a stage evaluation preempts rejection.** `Patch::set_ode_state`
calls `check_finite_ode_state()`
([patch.h:399-416](../inst/include/plant/patch.h#L399-L416), invoked at *every* RK stage), and
phylloptim's collar root-find throws `find_root_psi(find_root_crit=1) failed`
(`leaf_model.hpp:1356-1359`). Either fires before an error estimate exists, so the step is
never rejected — the run dies having taken **zero** steps. Fixing fault 1 does not help this at
all; in the table above the throw column is unchanged by it. Note `check_finite_ode_state`
exists deliberately (#550) to turn an opaque downstream error into an actionable one, so the
answer is not to delete it but to evaluate it at accepted-step boundaries rather than inside
stage evaluations.

## 4. Why the intuitive fix is the wrong one

The criterion used here — **preserve the Jacobian's lower-triangular structure, so the spectrum
stays diagonal-readable** — is not new to this note. It is the same principle
[`427-multipatch-water.md`](427-multipatch-water.md) (PR #565) applies *between* patches, where
donor control makes the cross-patch Jacobian triangular in topographic order so a downslope
sweep replaces a stiff coupled solve, under the stated house philosophy of *"design out
stiffness rather than reach for implicit solvers"*. This pass applies it *within* the column and
supplies measured backing for two claims that doc asserts qualitatively:

- **"Triangularising the coupling doesn't remove fast local drainage: near saturation K(θ) is
  large → stiff diagonal."** Confirmed and quantified: the diagonal is ≈1.5e4 yr⁻¹ at
  saturation (2.0e4 by the analytic `K′/dz`), and the structure is *exactly* triangular —
  max |upper-triangular entry| is 0 to machine precision.
- **"Closed budget for small/medium events; spill/deep-drainage loss term for large events."**
  Confirmed as necessary, with the threshold: see §5 — a ~13 mm event already exceeds layer 0's
  free capacity from a moderately wet start.

The convergence matters because it means the within-column and between-patch water work can share
one criterion and one justification, rather than each arguing stiffness from scratch.

Two candidate schemes were implemented and measured against that criterion.

| scheme | max &#124;upper-tri&#124; | lower-tri? | max &#124;λ&#124; |
|---|---|---|---|
| baseline (donor-only, as shipped) | 0 | yes | 1.546e4 |
| receiver rejects excess, water exits | 0 | yes | 1.546e4 |
| donor throttles when receiver full | 4.67e3 | **no** | 1.518e4 |

The criterion works as intended, and it vindicates the original donor-only choice: making the
donor's outflow depend on the receiver's wetness — the literal reading of "a layer can't
overfill the one below" — introduces super-diagonal coupling. Receiver-side rejection avoids
that. But receiver-side rejection is not worth having either:

- **It makes the spectrum stiffer**, by 1.20× at 90% of saturation and 1.45× at 99%
  (5045 vs 4188; 2.57e4 vs 1.77e4). The limiter adds its own diagonal term on top of the
  outflow's.
- **It changes behaviour substantially.** At θ₁/θ_sat = 0.99 the layer's rate goes from
  **+81 yr⁻¹ to −420 yr⁻¹**: the limiter rejects almost all inflow while the outflow is still
  138 m yr⁻¹, so the layer drains hard and the rejected water leaves the column. That is a real
  hydrological change (bypass/lateral loss) requiring a `scientific_version` bump and a
  scenario re-bless.
- **It buys nothing for the stated goal**, because §1 shows θ was never unbounded in the field
  and §2 shows the controller already rejects the overshooting steps.

Paying stiffness plus a re-bless to fix a bound that already holds is the wrong trade.

## 5. What the pulse feature does need

A pulse jump is applied **between** solver legs, outside the integrator, so no error estimate
and no rejection can protect it. It must be capped explicitly. Free capacity of layer 0 is
`(θ_sat − θ₀)·dz₀`:

| θ₀ | free capacity | 5 mm | 13.3 mm | 50 mm |
|---|---|---|---|---|
| 0.100 | 98.4 mm | 0.117 | 0.144 | 0.267 |
| 0.214 (default init) | 64.2 mm | 0.231 | 0.258 | 0.381 |
| 0.350 | 23.4 mm | 0.367 | 0.394 | **0.517** |
| 0.400 | 8.4 mm | 0.417 | **0.444** | **0.567** |

So a realistic dryland event (≈13 mm) already exceeds saturation from a moderately wet layer 0,
and a large event exceeds it from 0.35. **Cap the accepted depth at free capacity and route the
excess to a runoff accumulator.** Because the jump is outside the integrator a hard `min()` is
fine here — the "no kinks in the rates" concern applies only to continuous derivatives.

This is local to the pulse action, needs no change to `compute_rates`, and leaves every
existing TF24 trajectory bit-identical.

## 6. Recommendation

**Drop the planned physics change (Stage B).** Replace it with:

1. **Cap the pulse jump at free capacity** with an excess-to-runoff accumulator. Local, cheap,
   required, and the only item on this list that blocks the pulse feature. Note it adds a
   fifth accumulator, so `R/tidy_outputs.R:50-57` (positional flux names) and
   `test-patch.R:54-55` (pins the 5+4 ODE layout) both move.
2. **odelia: treat a non-finite error estimate as a rejection.** A controller bug fix,
   unrelated to events, and guarded — non-finite `yerr` cannot occur on a run that currently
   passes, so nothing existing moves.
3. **Make an out-of-domain probe rejectable rather than fatal.** Move
   `check_finite_ode_state()` off the per-stage path, and have the collar root-find signal an
   unphysical probe in a way that becomes a rejected step. Without this, item 2 alone leaves
   the throw path dying at step 0.
4. **Reconsider inheriting `step_size_last` across a discrete change.** It is the trigger for
   everything in §3, and a pulse creates far more leg boundaries than introductions do.

Items 2–4 span odelia, plant and phylloptim, so they need a decision before implementation.
They are also the real fix for **#599**, which should be re-scoped accordingly: the cause is
not a missing saturation bound in the water balance, and `32240725`'s clamps were measured
ineffective (5/40 → 4/40) because they treat a symptom of the step size rather than the
error-handling faults above.

**Not needed, on this evidence:** receiver-capacity limiting in `compute_rates`, a θ ≤ θ_sat
clamp in `psi_from_soil_moist`, and any change to the inter-layer cascade. Also not needed for
this: an implicit/IMEX stepper — worth revisiting for *speed* (the diagonal is ≈1.5e4 yr⁻¹
regardless), but it is not what stands between us and a working pulse.

### Bearing on the multi-patch water design (#427 / PR #565)

That doc's process 1 (event-time surface redistribution) and process 3 (donor-controlled
subsurface flux) are unaffected by anything here — donor control is exactly the right choice and
this pass reinforces it. Two things carry across:

- Its process 2 (within-patch vertical bucket, N = 15) inherits the faults in §3, not a
  water-balance defect. **More layers make the trigger worse**, because the diagonal scales as
  `1/dz`. Measured at 0.99 θ_sat: **1.77e4 yr⁻¹ at N = 5 (dz = 0.3 m) → 5.30e4 at N = 15
  (dz = 0.1 m)**, a factor of 3 as expected (analytic `K′/dz` at saturation: 2.05e4 → 6.15e4).
  So the step a leg *inherits* matters more at N = 15, not less.
- Its per-patch stiffness self-diagnosis ("wet ones stiff, dry Mulga patches not") is sound, but
  should key off the diagonal `K′(θ)/dz`, which is measurable directly, rather than off observed
  failures — because from a cold start even the worst column integrates cleanly, so absence of
  failure is not evidence of absence of stiffness.

Finally, fix the stale comment at
[tf24_environment.h:363](../inst/include/plant/models/tf24_environment.h#L363), which describes
"the explicit fixed-step solver" — there is none; the correct comment 30 lines later
contradicts it, and this comment is the origin of #522's fixed-step framing.
