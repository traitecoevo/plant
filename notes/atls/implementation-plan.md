# Plan: Modular thermal-response strategy layer for TF24 (ATLS / Lumry–Eyring, #566)

## Context

TF24 currently runs at a constant `leaf_temp = 25 °C` and has no representation of
heat damage or thermal acclimation. We want a leaf thermal-response layer in which
**different thermal strategies are points in a continuous trait space** — differing
along four axes (tolerance, acclimation, repair, avoidance-via-cooling) — all
plugging into the *same* underlying physics/biochemistry, so we can ask which
strategy wins under a given climate (esp. Australian heatwaves).

This formalises epic #566 and builds directly on #523 (Penman–Monteith leaf energy
balance). PM is a hard prerequisite: the damage kinetics are driven by `Tleaf`,
which only becomes physically meaningful once the energy balance is in place.

**Decisions locked with the user:**
1. Acclimation acts on **both `T_opt` and `T_crit`, via two separate ODE states**
   (`A_opt`, `A_crit`) with independent kinetics.
2. Damage integration = **single midday evaluation**: evaluate the leaf once at the
   midday (peak-temperature) operating point as the representative daily condition,
   giving a quasi-steady damage factor `N(Tleaf_midday, A_crit)`. No fast ODE, no
   diurnal integration.
3. Build **all four axes at once**.
4. **Costs from the start** — each axis carries a trade-off so competition is
   meaningful.

## The central distinction: `T_opt` vs `T_crit`

- **`T_opt`** — optimum of photosynthetic *capacity*. Lives in the reversible
  **peaked-Arrhenius** (Medlyn 2002) response of Vcmax/Jmax (`peak_arrh_curve`,
  [leaf_model.cpp:246-247](../../src/leaf_model.cpp#L246-L247)). The peak is set by
  the two deactivation terms: `H_d` (deactivation energy, high-T collapse; ≈200
  kJ/mol) and `d_S` (entropy of deactivation; ≈650 J/mol/K). `T_opt` rises as `d_S`
  falls (roughly `T_opt ≈ H_d/d_S`; defaults give ≈30 °C). So `d_S` is the knob for
  shifting `T_opt`. No damage; fully reversible.
- **`T_crit`** — *damage* threshold. Lives in the **Lumry–Eyring** kinetics: the
  temperature where irreversible loss runs away, entering through the damage switch
  `S_damage = 1/(1+exp(−m·(Tleaf − T_crit(A_crit))))`. **The switch dominates, not
  Arrhenius:** with the ATLS default rate activation energies near zero
  (`E_d ≈ 400–500 J/mol` vs `RT ≈ 2500`), the Arrhenius factor is nearly flat across
  the biological range, so damage is effectively `k_0·S_damage(Tleaf; T_crit)` — a
  sharp turn-on near `T_crit`. This matches the demo behaviour; `T_crit` and slope
  `m` are the real controls, the `E_d` terms near-inert (keep for ATLS fidelity).

Two distinct axes; a strategy can move one without the other. In the ATLS demo,
acclimation moved only `T_crit`; here we let acclimation move both, independently.

## Shared core equations (every strategy plugs into these)

1. **Leaf temperature** — PM energy balance from #523,
   `Tleaf = Tair + (Rn − λE)·ra/(ρcp)`
   ([leaf_model.cpp:1093](../../src/leaf_model.cpp#L1093), branch
   `feature/penman-monteith-leaf-energy-balance`), with `E` pinned by the existing
   Sperry hydraulic supply. Gate `use_energy_balance`.
2. **Reversible capacity** — existing peaked-Arrhenius Vcmax/Jmax and Arrhenius
   Ko/Kc/Γ* at `Tleaf` ([leaf_model.cpp:246-252](../../src/leaf_model.cpp#L246-L252)).
   Sets `T_opt`.
3. **Lumry–Eyring damage** — protein pool `N ⇌ D → I`, evaluated **once at the
   midday operating point** at quasi-steady to yield a functional fraction
   `N ∈ (0,1]` that stands in for the day. `N` downscales electron transport via
   `jmax → jmax·N` at [leaf_model.cpp:247](../../src/leaf_model.cpp#L247). Sets
   `T_crit`.
4. **Acclimation** — two slow relaxation ODE states `A_opt`, `A_crit`, driven by the
   thermal regime, decaying over days. Slow enough to live at the demographic
   timescale (unlike the fast damage kinetics, which the midday evaluation collapses
   to a quasi-steady factor).

Reference equations/params: external ATLS demo at
`/Users/z2209343/GitHub/projects/ATLS_leaf_demo/` (`theory.qmd`,
`src/leaf_thermal_system.hpp`, `src/leaf_thermal_pars.hpp`; defaults `T_crit_0=38`,
`ΔTmax=6`, `K_A=1`, `k_d1_0=k_d2_0=k_r1_0=864 d⁻¹`, `E_d1=400`, `E_d2=500`,
`E_r1=500`, `m_switch=1`, `m_rep=0.4`, `T_rep_cut=45`, `T_accl=30`, `α=0.02`,
`β=0.01`).

## The four strategy axes → concrete knobs + costs

Cost principle (user): the **intrinsic optimum** (offset = 0) is the temperature
suiting the default protein structure. Moving it constitutively is a **construction
cost** (building more complex/stable proteins). **Instantaneous costs live with the
dynamic responses** (acclimation, repair), not with constitutive tolerance.

| Axis | Moves | Mechanism (knob) | Cost (from start) |
|---|---|---|---|
| **Tolerance** | `T_opt`, `T_crit` (constitutive) | trait offsets from the intrinsic optimum: shift Vcmax/Jmax `d_S` (→ `T_opt`); raise baseline `T_crit_0` (→ `T_crit`) | **construction cost only**: higher `T_opt`/`T_crit` ⇒ more complex/stable proteins ⇒ higher leaf construction cost. No instantaneous cost |
| **Acclimation** | `T_opt`, `T_crit` (inducible) | `A_opt` shifts `d_S`; `A_crit` shifts `T_crit(A_crit)=T_crit_0+ΔTmax·A/(A+K_A)`; both rise with heat, decay | **maintenance** respiration ∝ (`A_opt+A_crit`) **+ induction/construction** cost ∝ positive `dA/dt` (building the proteins) |
| **Repair** | `T_crit` (recovery) | trait scaling refold rate `k_r1_0` (D→N); keeps `D` low, starves D→I | **standing** maintenance ∝ repair *capacity* `k_r1_0` **+ activity** cost ∝ realized repair *flux* `k_r1·D` |
| **Avoidance** | `Tleaf` | higher transpiration/`g` → lower `Tleaf` via PM | **endogenous**: water cost already priced by Sperry hydraulics; no new term |

Mechanistic notes:
- Acclimation/repair act on the **repair/aggregation** side (`k_r1↑`, damage-switch
  shift), *not* the unfolding rate; constitutive `T_crit` tolerance is the one lever
  on intrinsic thermostability (`T_crit_0`).
- Cost split honours the principle above: **construction** (one-off, tied to
  building the tissue/protein pool) for constitutive tolerance and the *inducible*
  part of acclimation; **maintenance** (standing, ∝ pool held) for acclimation and
  repair capacity; **activity** (∝ realized flux) for actual repair.
- The kcat–thermostability trade-off (thermostable enzymes are slower, so tolerance
  could instead reduce `vcmax_25`/`jmax_25`) is a noted *alternative* to the
  construction-cost framing — set aside per the user, revisit if construction cost
  proves too weak to constrain the trait space.
- Avoidance is self-selecting once `N` feeds the assimilation the psi_stem optimiser
  sees: cooler midday leaf → higher `N` → higher carbon → optimiser values cooling.

## Architecture: what's an ODE state vs algebraic

- **New ODE states (slow, on the strategy):** `A_opt`, `A_crit`. Relaxation
  dynamics `dA/dt = α·softplus(regime − T_accl) − β·A` (softplus, not `max(0,·)`, so
  the forcing is C¹ — see numerical-stability section).
- **Algebraic (fast, on the Leaf):** the functional fraction `N`, computed each
  `set_physiology` from the midday `Tleaf` operating point + current `A_crit`. No
  `D`/`I` ODE states are integrated.

## Code changes

**Model shape — new `TF24t_Strategy : public TF24_Strategy`** (mirror TF24f exactly;
[tf24f_strategy.h](../../inst/include/plant/models/tf24f_strategy.h)). Keeps base
TF24 + #523 PM untouched and backward-compatible. Nested hierarchy: TF24 ⊂ TF24+PM
(`use_energy_balance`) ⊂ TF24t (+PM+ATLS, `use_thermal_damage`).
- `state_size()` = `TF24_Strategy::state_size() + 2`; `state_names()` append
  `"acclim_topt"`, `"acclim_tcrit"` (indices 0..N-1 unchanged).
- `refresh_indices()` → cache `state_idx_acclim_topt/_tcrit` (call base first).
- `compute_rates()` → delegate to base, then set the two `dA/dt` rates from the
  thermal regime.
- `set_initial_states()` → base first, then seed `A_opt=A_crit=0` (or acclimated
  equilibrium for the birth environment).

**Leaf-level (gated, backward-compatible), in
[leaf_model.h](../../inst/include/plant/leaf_model.h) /
[leaf_model.cpp](../../src/leaf_model.cpp):**
- Add `use_thermal_damage_` member + `use_thermal_damage` `TF24_Pars` field,
  mirroring the existing `use_energy_balance_` gate
  ([leaf_model.h:240](../../inst/include/plant/leaf_model.h#L240),
  [tf24_strategy.h:101](../../inst/include/plant/models/tf24_strategy.h#L101)).
  Default off → identical to current TF24 (bit-for-bit; unlike TF24f).
- New `Leaf::thermal_damage_factor(A_crit, thermal regime)` → time-averaged `N`,
  and apply the `A_opt`-driven `d_S` shift, inside the temp block at
  [leaf_model.cpp:246-252](../../src/leaf_model.cpp#L246-L252). Apply `N` to
  `jmax_` at [leaf_model.cpp:247](../../src/leaf_model.cpp#L247).
- Apply the tolerance **construction cost** at leaf build (raise leaf C cost with
  `T_opt`/`T_crit` offsets, in the carbon economy, not as a `vcmax_25` penalty), and
  add the acclimation/repair **maintenance + activity** terms to `R_d_`
  ([leaf_model.cpp:251](../../src/leaf_model.cpp#L251)): `R_d_ += r_accl_maint·(A_opt+A_crit) + r_repair_maint·k_r1_0 + c_repair_flux·(k_r1·D)`,
  plus an induction cost ∝ positive `dA/dt` charged on the acclimation states.
- Extend `set_physiology(...)` signature to receive `A_opt`, `A_crit` from the
  strategy ([leaf_model.cpp:204](../../src/leaf_model.cpp#L204); call site
  [tf24_strategy.cpp:440](../../src/tf24_strategy.cpp#L440)).
- **Cache:** the temp block memoises on `(leaf_temp_, atm_o2_kpa_)`
  ([leaf_model.cpp:243-256](../../src/leaf_model.cpp#L243-L256)); extend the cache
  key to include `A_opt`, `A_crit` (and the midday `Tleaf`), or the damage factor
  and `d_S` shift will be stale.

**Parameters — `TF24_Pars`** (C++ struct
[tf24_strategy.h:20-97](../../inst/include/plant/models/tf24_strategy.h#L20-L97) +
YAML [RcppR6_classes.yml](../../inst/RcppR6_classes.yml) `TF24_Pars` list, both must
match): add tolerance offsets (`topt_offset`, `tcrit_0`), LE rate constants
(`k_d1_0,k_d2_0,k_r1_0,E_*`, switch slopes), acclimation kinetics
(`alpha_opt,beta_opt,dTopt_max,alpha_crit,beta_crit,dTcrit_max,K_A,T_accl`), repair
scaling, and cost coefficients (construction `c_build_topt`,`c_build_tcrit`;
acclimation `r_accl_maint`,`c_accl_induct`; repair `r_repair_maint`,`c_repair_flux`).
Expose
`TF24t`-specific scalars in the strategy YAML `list:` (mirror TF24f's
`k_acclim`/`psi_fd_step`, [RcppR6_classes.yml] `TF24f_Strategy`).

**Environment — `TF24_Environment`** ([tf24_environment.h](../../inst/include/plant/models/tf24_environment.h)):
supply the **midday** air temperature for the damage evaluation (the existing
radiation/PPFD forcing already represents peak conditions). Reuse the existing
`leaf_temp`/`Tair` driver at the midday operating point — no diurnal-amplitude
driver needed. Mirror the constant-driver pattern
([tf24_environment.h:86](../../inst/include/plant/models/tf24_environment.h#L86)).

**Register `TF24t` everywhere** the model list appears in
[RcppR6_classes.yml](../../inst/RcppR6_classes.yml) (SCM/Species/Patch/Node/
Parameters/IndividualRunner blocks), then regenerate: `make RcppR6` then
`make rebuild`. The scaffolder
[scripts/new_strategy_scaffolder.R](../../scripts/new_strategy_scaffolder.R) (skill
`plant-new-strategy`) can generate the boilerplate.

## Modeling decisions to settle during implementation (recommended defaults)

1. **Irreversible `I`.** True steady state has `I→1` (monotonic). *Default:* base
   `N` on the reversible `N⇌D` equilibrium (repairable) at the midday operating
   point; log cumulative `I` as a diagnostic hook for future leaf-turnover/mortality
   coupling. Keeps `N` bounded and physiological.
2. **Avoidance coupling into `N`.** `N` is evaluated at the midday operating point,
   so a plant's cooling investment lowers midday `Tleaf` and raises `N` directly —
   avoidance pays off inside the existing psi_stem optimiser with no diurnal
   machinery.
3. **`T_opt` ↔ `d_S` mapping.** Expose a target `T_opt` offset and invert to a `d_S`
   shift for Vcmax and Jmax (peaked-Arrhenius peak location), rather than exposing
   `d_S` directly.

## Numerical stability & robustness (requirement)

The layer feeds the ODE integrator and the `run_mutant` replay cache; a single
NaN/Inf poisons the whole SCM run. So it must be smooth in every driver and trait
and degrade gracefully on outlandish values.

Design rules:
- **Smooth everywhere (C¹).** Every temperature/threshold response is a smooth
  function of drivers and traits. Use **softplus** `(1/s)·log1p(exp(s·x))` for the
  acclimation forcing (not `max(0,·)`); keep the logistic `S_damage`/`S_repair`
  switches. No `if`-branches on driver/trait *values* inside the response.
- **Bounded / saturating forms.** `T_crit(A)=T_crit_0+ΔTmax·A/(A+K_A)` saturates;
  switches ∈ [0,1]; `N=1−D−I ∈ (0,1]`; acclimation states self-limit via `−β·A`
  decay. Each new expression must have a finite limit as any input → ±∞.
- **Guard the transcendentals.** Clamp logistic/Arrhenius/softplus arguments to a
  safe range (e.g. |arg| ≲ 30) before `exp` so extreme `Tleaf`, slope `m`, or trait
  offsets can't overflow; use `log1p`/`expm1` near 0.
- **Protect denominators & parameter domains.** `K_A>0` (guard `A+K_A`),
  `k_r1_0≥0`; clamp the `d_S` shift so the peaked-Arrhenius stays peaked and finite
  (`H_d−Ha>0`, `T_opt` finite) — reject/clamp values that would invert the curve;
  construction-cost-adjusted leaf C cost / `vcmax_25` stay strictly positive.
- **Finite-output contract.** `thermal_damage_factor` and the acclimation rates
  return finite, bounded values for *any* finite input; debug-assert (release-mode
  clamp) that no NaN/Inf reaches `ode_rates`.
- **Stiffness.** Large switch slopes `m` / acclimation gains approach a step and
  stiffen the integrator; keep defaults moderate and document `m` as the
  smoothness↔sharpness knob. The midday-only damage evaluation (no fast `D/I` ODE)
  already removes the stiffest timescale.

## Build order (all four axes, phased for testability)

1. Leaf-level LE damage `N` (midday quasi-steady) + `use_thermal_damage` gate + `jmax·N`;
   constitutive tolerance (`T_opt` via `d_S`, `T_crit_0`) + repair (`k_r1_0`).
   Validate `N` against the ATLS demo.
2. `TF24t_Strategy` with `A_opt`, `A_crit` ODE states + acclimation kinetics.
3. Cost structure (tolerance construction cost; acclimation/repair maintenance +
   activity terms on `R_d_`).
4. Wire the midday evaluation (midday `Tair` / operating point → `N`).
5. Registration + RcppR6 regen + R interface.

## Verification

- **Backward compat:** `use_thermal_damage=0` ⇒ TF24t reproduces TF24(+PM)
  bit-for-bit (contrast TF24f, which does not). Use the `origin/chore/refresh-pre-pm-baselines`
  discipline (#568).
- **Unit tests** (`tests/testthat/`, mirror the demo's
  `ATLS_leaf_demo/tests/testthat/test-leaf_thermal.R`): `thermal_damage_factor` `N(Tleaf, A_crit)`
  matches ATLS steady-state values; acclimation ODEs relax correctly (rise under
  heat, decay after); cost terms move `vcmax_25`/`R_d_` as specified.
- **Single-plant integration:** run TF24t under a hot day (high midday `Tair`) —
  confirm `N<1`, downscaled assimilation, and `A_opt`/`A_crit` raising
  `T_opt`/`T_crit` over days; confirm a higher-transpiration strategy runs cooler at
  midday with higher `N`; confirm construction/maintenance/activity costs debit the
  carbon balance as specified.
- **Evaluation gate (#566 / mirror #523 step-5 factorial,
  `notes/penman-monteith/step5_fick_vs_pm.R`):** does the damage feedback materially
  change annual carbon gain / competitive outcome vs PM-only under representative
  Australian heatwaves? Keep only if material.
- **Property tests (smoothness & bounds):** sweep each driver (midday
  `Tair`/`Tleaf`, radiation) and each trait (`T_opt`/`T_crit` offsets, `k_r1_0`,
  slope `m`, acclimation gains) across a wide grid including extremes; assert
  outputs finite, `N∈(0,1]`, monotone where expected, and C¹-smooth (finite-
  difference derivative shows no spike at thresholds).
- **Edge cases:** `Tleaf≪T_crit` (`N≈1`); `Tleaf≫T_crit` (`N` small but >0); `A=0`;
  `A` huge; `k_r1_0=0` (no repair); zero transpiration (avoidance off); `K_A→0⁺`;
  `m` very large (near-step); negative and very large trait offsets.
- **Outlandish inputs:** `Tair=±100 °C`, offsets of ±50 °C, gains orders of
  magnitude off — model returns sensible bounded output (saturates, no NaN/Inf) and
  a single-plant run still completes.
- Build/test loop: `make rebuild` then `devtools::test()`.

## Constraints

- **Branching:** new feature branch + PR; never commit/push to `master`.
- **Depends on #523** — build on `feature/penman-monteith-leaf-energy-balance`
  (base it on that branch or land after PM merges).
- **odelia:** build plant against odelia **master** (mutant run cache-hook names),
  not a feature branch — odelia is currently on `pr-41`. The new states flow through
  the standard `set_ode_state`/`ode_rates` path, so `run_mutant` replay works
  automatically (the `Replayable` hooks are compile-time-optional).

---

## Implementation progress & handoff

Branch: **`feature/tf24-atls-thermal-damage`** (based on
`feature/penman-monteith-leaf-energy-balance`). Built against the
**currently-installed odelia (`pr-41`)** for now (user decision; the odelia-master
note above still applies before landing).

**Phase 1 — leaf damage core — DONE, committed `b2cbd5f7`.**
- `leaf_model.{h,cpp}`: `use_thermal_damage_` gate (default off ⇒ bit-identical);
  switch-dominated quasi-steady `N = k_r/(k_r+k_d)` (`thermal_damage_factor`),
  `t_crit(A_crit)`, `d_S_shifted` (T_opt shift via Medlyn d_S inversion), guarded
  `logistic_`/`softplus_`; gated temp block applies `jmax → jmax·N`.
- New Leaf fields exposed via RcppR6. Tests: `test-leaf-thermal.R` (14 checks).

**Phase 2 — TF24t strategy + acclimation ODE states — DONE, committed `7aa008ff`.**
- `TF24t_Strategy : TF24_Strategy` (`tf24t_strategy.{h,cpp}`): appends
  `acclim_topt`, `acclim_tcrit`; `dA/dt = alpha·softplus(T−T_accl) − beta·A`,
  seeded at environmental equilibrium; `prepare_strategy` turns on the leaf layer
  and copies thermal traits; `compute_rates` feeds live states into the leaf.
- Reuses parent `TF24_Pars`; thermal/acclimation knobs are TF24t strategy scalars.
  `strategy_version.cpp` has a TF24t compound-version case. Leaf temp-cache guard
  extended so varying `A_opt`/`A_crit` isn't masked. Scaffolded build wiring
  (reuse-environment mode); fixed the scaffolder's `g1_TF24` rename; TF24t excluded
  from generic cross-strategy helper lists. Tests: `test-strategy-tf24t.R` (15).
- Regression: 437 shared checks pass (leaf, leaf-thermal, tf24, environment,
  model-version, tf24f).

**Phase 3 — cost structure — IN PROGRESS. Locked design + as-built map:**

Two homes for the costs, mirroring the physics:
- **Leaf-level respiration (`R_d_`)** in `update_temperature_dependent_params`,
  gated by `use_thermal_damage_`, coefficients default 0 (a bare Leaf pays
  nothing — Phase 1 leaf tests untouched), set nonzero by TF24t:
  - acclimation **maintenance** ∝ held load `(A_opt+A_crit)` → `c_acclim_maint_`
  - repair **standing maintenance** ∝ repair *capacity* `k_r1_0` (paid even when
    cold) → `c_repair_maint_`
  - repair **activity** ∝ realized refold flux `k_r1·(1−N)` (D≈1−N) →
    `c_repair_flux_`
- **Whole-plant carbon balance** via a virtual override
  `TF24t_Strategy::net_mass_production_dt` (base `net_mass_production_dt` is
  already `virtual`, dispatches correctly even when the base `compute_rates`
  calls it — so no base-class edit and base TF24 stays bit-identical):
  - tolerance **construction cost** — a one-off build cost per leaf amortised
    through leaf turnover: `(c_build_topt·max(0,topt_offset) +
    c_build_tcrit·max(0,tcrit_0−38))·turnover_leaf(mass_leaf)`. Applies at
    establishment too (establishment_probability shares the virtual). Uses the
    ATLS baseline `T_crit_0=38` as the intrinsic reference.
  - acclimation **induction cost** ∝ smoothed positive build rate of each
    acclimation state: `c_accl_induct·(pos(dA_opt/dt)+pos(dA_crit/dt))`, with
    `pos(x)=½(x+√(x²+ε²))` (the same C∞ positive-part used for net production).
    Read from `leaf.A_opt_`/`leaf.A_crit_` + env temp — which `compute_rates`
    sets before base `compute_rates` calls the override, so they are the current
    states. (At establishment the leaf states read 0, giving a negligible
    build-from-zero snapshot cost at default temps — documented, not fixed.)

- **Cost coefficients ship modest nonzero, documented, R-settable, flagged as
  CALIBRATION TARGETS** (real trade-off out of the box). Defaults chosen so a
  *default* TF24t (no extra tolerance offset) pays no construction cost — you
  only pay for tolerance you buy — while acclimation/repair costs bite whenever
  those axes are engaged.

- **Single-plant check (user request):** confirm an individual TF24t plant runs
  through the Individual/OdeRunner path (build `Individual("TF24t", env)`, step
  the `OdeRunner`), grows, stays finite, and that raising the cost coefficients
  measurably reduces growth (costs debit carbon). Added as a test in
  `test-strategy-tf24t.R`.

Units note (flagged, not a Phase-3 fix): acclimation kinetics are labelled day⁻¹
but integrated on the SCM's yearly clock, so `dA/dt` magnitude — and hence the
induction-cost scale — is uncertain; the calibration-target coefficients absorb
this until the timescale is reconciled (Phase 2 follow-up).

**Phase 4 — midday evaluation wiring:** make `N` evaluate at the midday operating
point (add a midday `Tair` driver / use the midday operating point), rather than
whatever `leaf_temp`/PM `Tleaf` is currently supplied.

**Phase 5 — leftover:** add `TF24t` to the `test-model-version.R` `models` vector
+ accept the new `_snaps/model-version.md` entry (registration + R interface are
already done via the scaffolder).

**Verification gate (#566):** does the damage feedback materially change annual
carbon gain / competitive outcome vs PM-only under representative Australian
heatwaves (mirror the #523 step-5 factorial).

**Build/test loop:** `make rebuild` then, per-file,
`Rscript -e 'pkgload::load_all(".", compile=FALSE); testthat::test_file("tests/testthat/<f>")'`.
Use `run_scm(..., refine_schedule = FALSE)` for fast smoke tests.
